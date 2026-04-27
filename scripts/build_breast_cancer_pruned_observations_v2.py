#!/usr/bin/env python3
"""Create a clone-aware breast-cancer observations file with variant detail rows.

The script reads breast_cancer_clone_proportions.csv and observations.csv, drops
genes whose clone VAF is below a threshold at a given sequencing timepoint, and
rewrites the remaining genomic sequencing content into standardized observation
groups per variant.

Non-genomic observations are kept unchanged.

Genomic sequencing rows are normalized into one group per retained variant using:
- the selected MAF row for that patient-gene
- the clone VAF for that patient-gene at that sequencing timepoint

Each group reuses the original DATE, PATIENT, and ENCOUNTER values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

from build_breast_cancer_clones_v2 import (
    BASELINE_POSITIVE_SUFFIX,
    GENETIC_ASSESSMENT_DESCRIPTION,
    HER2_IHC_DESCRIPTION,
    HER2_POSITIVE_VALUE,
    newest_observations_csv,
    POST_TREATMENT_PREFIX,
)


DEFAULT_PROPORTIONS_NAME = "breast_cancer_clone_proportions.csv"
DEFAULT_OUTPUT_NAME = "observations_pruned_by_clone_vaf.csv"
DEFAULT_DRIVER_VARIANTS_NAME = "breast_cancer_driver_variants_from_maf.csv"
DEFAULT_NON_DISRUPTIVE_VARIANTS_NAME = "breast_cancer_non_disruptive_variants_from_maf.csv"
HER2_FISH_DESCRIPTION = "HER2 [Presence] in Breast cancer specimen by FISH"
GENE_ALIASES = {"HER2": "ERBB2"}

GENE_CODE = "48018-6"
GENE_DESCRIPTION = "Gene studied [ID]"
GENOMIC_CHANGE_CODE = "81290-9"
GENOMIC_CHANGE_DESCRIPTION = "Genomic DNA change (gHGVS)"
LOCAL_VARIANT_CLASS_CODE = "83005-9"
LOCAL_VARIANT_CLASS_DESCRIPTION = "Variant classification"
VARIANT_TYPE_CODE = "48019-4"
VARIANT_TYPE_DESCRIPTION = "DNA change type"
VAF_CODE = "81258-6"
VAF_DESCRIPTION = "Sample variant allelic frequency [NFr]"
CLONAL_ROLE_CODE = "CLONAL_ROLE"
CLONAL_ROLE_DESCRIPTION_SUFFIX = " - clonal role"
FOUNDING_ROLE_VALUE = "founding"
SUBCLONAL_ROLE_VALUE = "subclonal"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a copy of observations.csv where breast-cancer sequencing rows "
            "are replaced by standardized variant observation groups."
        )
    )
    parser.add_argument(
        "--clone-proportions",
        type=Path,
        help=(
            "Path to breast_cancer_clone_proportions.csv. Defaults to the newest "
            "output_runs/*/csv/breast_cancer_clone_proportions.csv."
        ),
    )
    parser.add_argument(
        "--observations",
        type=Path,
        help=(
            "Path to observations.csv. Defaults to observations.csv in the same "
            "directory as the clone proportions file."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Output CSV path. Defaults to observations_pruned_by_clone_vaf.csv in "
            "the same directory as the input observations.csv."
        ),
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=5.0,
        help="Drop genes at any timepoint with clone VAF below this percentage. Default: 5.0",
    )
    parser.add_argument(
        "--driver-variants",
        type=Path,
        help=(
            "Path to breast_cancer_driver_variants_from_maf.csv. Defaults to "
            "scripts/breast_cancer_driver_variants_from_maf.csv."
        ),
    )
    parser.add_argument(
        "--non-disruptive-variants",
        type=Path,
        help=(
            "Path to breast_cancer_non_disruptive_variants_from_maf.csv. Defaults to "
            "scripts/breast_cancer_non_disruptive_variants_from_maf.csv."
        ),
    )
    return parser.parse_args()


def newest_clone_proportions_csv(repo_root: Path) -> Path:
    candidates = sorted(
        repo_root.glob(f"output_runs/*/csv/{DEFAULT_PROPORTIONS_NAME}"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(f"No {DEFAULT_PROPORTIONS_NAME} found under output_runs/*/csv/")
    return candidates[0]


def extract_baseline_gene(value: str) -> str | None:
    if value.endswith(BASELINE_POSITIVE_SUFFIX):
        return value[: -len(BASELINE_POSITIVE_SUFFIX)]
    return None


def canonical_gene_name(gene: str) -> str:
    return GENE_ALIASES.get(gene, gene)


def is_genomic_source_row(row: dict[str, str]) -> bool:
    description = row["DESCRIPTION"].strip()
    return (
        description == GENETIC_ASSESSMENT_DESCRIPTION
        or description == HER2_IHC_DESCRIPTION
        or description == HER2_FISH_DESCRIPTION
        or description.startswith(POST_TREATMENT_PREFIX)
    )


def source_row_to_gene(row: dict[str, str]) -> str | None:
    description = row["DESCRIPTION"].strip()
    value = row["VALUE"].strip()

    if description == GENETIC_ASSESSMENT_DESCRIPTION:
        return extract_baseline_gene(value)

    if description == HER2_IHC_DESCRIPTION and value == HER2_POSITIVE_VALUE:
        return "HER2"

    if description.startswith(POST_TREATMENT_PREFIX):
        return description[len(POST_TREATMENT_PREFIX) :].strip()

    return None


def load_low_vaf_targets(
    clone_proportions_path: Path, threshold: float
) -> tuple[dict[str, set[str]], dict[str, dict[str, set[str]]], Counter]:
    baseline_genes_to_remove: dict[str, set[str]] = defaultdict(set)
    dated_genes_to_remove: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    stats = Counter()

    with clone_proportions_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            patient_id = row["patient_id"].strip()
            genes = [gene.strip() for gene in row["genes"].split(";") if gene.strip()]
            if not genes:
                continue

            for label in ("t0", "t1", "t2"):
                date = row[f"{label}_date"].strip()
                value = row[f"{label}_vaf_pct"].strip()
                if not date or not value:
                    continue

                try:
                    vaf = float(value)
                except ValueError:
                    continue

                if vaf >= threshold:
                    continue

                stats[f"{label}_clone_rows_below_threshold"] += 1
                stats[f"{label}_genes_marked_for_removal"] += len(genes)
                if label == "t0":
                    baseline_genes_to_remove[patient_id].update(genes)
                else:
                    dated_genes_to_remove[patient_id][date].update(genes)

    return baseline_genes_to_remove, dated_genes_to_remove, stats


def load_clone_vaf_map(
    clone_proportions_path: Path,
) -> dict[str, dict[str, dict[str, float]]]:
    gene_vafs: dict[str, dict[str, dict[str, float]]] = defaultdict(lambda: defaultdict(dict))

    with clone_proportions_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            patient_id = row["patient_id"].strip()
            genes = [gene.strip() for gene in row["genes"].split(";") if gene.strip()]
            if not genes:
                continue

            for label in ("t0", "t1", "t2"):
                date = row[f"{label}_date"].strip()
                value = row[f"{label}_vaf_pct"].strip()
                if not date or not value:
                    continue
                try:
                    vaf = float(value)
                except ValueError:
                    continue

                for gene in genes:
                    gene_vafs[patient_id][date][gene] = vaf
                    canonical_gene = canonical_gene_name(gene)
                    if canonical_gene != gene:
                        gene_vafs[patient_id][date][canonical_gene] = vaf

    return gene_vafs


def load_founding_gene_map(clone_proportions_path: Path) -> dict[str, set[str]]:
    founding_genes: dict[str, set[str]] = defaultdict(set)

    with clone_proportions_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["clone_type"].strip() != "founding":
                continue
            patient_id = row["patient_id"].strip()
            genes = [gene.strip() for gene in row["genes"].split(";") if gene.strip()]
            for gene in genes:
                founding_genes[patient_id].add(gene)
                canonical_gene = canonical_gene_name(gene)
                founding_genes[patient_id].add(canonical_gene)

    return founding_genes


def lookup_gene_vaf(
    row: dict[str, str],
    gene: str,
    gene_vafs: dict[str, dict[str, dict[str, float]]],
) -> float | None:
    patient_id = row["PATIENT"].strip()
    date = row["DATE"].strip()
    patient_vafs = gene_vafs.get(patient_id, {})
    exact = patient_vafs.get(date, {}).get(gene)
    if exact is not None:
        return exact

    description = row["DESCRIPTION"].strip()
    if description == GENETIC_ASSESSMENT_DESCRIPTION or description == HER2_IHC_DESCRIPTION:
        for candidate_date in sorted(patient_vafs):
            vaf = patient_vafs[candidate_date].get(gene)
            if vaf is not None:
                return vaf
    return None


def row_should_be_removed(
    row: dict[str, str],
    baseline_genes_to_remove: dict[str, set[str]],
    dated_genes_to_remove: dict[str, dict[str, set[str]]],
) -> bool:
    patient_id = row["PATIENT"].strip()
    date = row["DATE"].strip()
    gene = source_row_to_gene(row)
    if not gene:
        return False

    if row["DESCRIPTION"].strip() == GENETIC_ASSESSMENT_DESCRIPTION or row["DESCRIPTION"].strip() == HER2_IHC_DESCRIPTION:
        return gene in baseline_genes_to_remove.get(patient_id, set())

    return gene in dated_genes_to_remove.get(patient_id, {}).get(date, set())


def load_variant_rows(variant_path: Path) -> dict[str, list[dict[str, str]]]:
    rows_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    with variant_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            gene = row["helper_gene"].strip()
            if gene:
                rows_by_gene[gene].append(row)
    return rows_by_gene


def choose_patient_gene_variant(
    patient_id: str,
    gene: str,
    driver_variants: dict[str, list[dict[str, str]]],
    non_disruptive_variants: dict[str, list[dict[str, str]]],
) -> tuple[dict[str, str], str]:
    canonical_gene = canonical_gene_name(gene)
    pool = driver_variants.get(canonical_gene)
    source = "driver"
    if not pool:
        pool = non_disruptive_variants.get(canonical_gene)
        source = "non_disruptive"
    if not pool:
        return {}, "missing"

    seed_input = f"{patient_id}|{canonical_gene}".encode("utf-8")
    seed = int(hashlib.sha256(seed_input).hexdigest(), 16)
    rng = random.Random(seed)
    return dict(rng.choice(pool)), source


def build_observation_row(
    base_row: dict[str, str],
    code: str,
    description: str,
    value: str,
    value_type: str,
    units: str = "",
) -> dict[str, str]:
    return {
        "DATE": base_row["DATE"],
        "PATIENT": base_row["PATIENT"],
        "ENCOUNTER": base_row["ENCOUNTER"],
        "CATEGORY": "laboratory",
        "CODE": code,
        "DESCRIPTION": description,
        "VALUE": value,
        "UNITS": units,
        "TYPE": value_type,
    }


def format_vaf_value(vaf_pct: float | None) -> str:
    if vaf_pct is None:
        return ""
    return f"{vaf_pct / 100.0:.4f}"


def build_variant_observation_group(
    base_row: dict[str, str],
    gene: str,
    variant_row: dict[str, str],
    vaf_pct: float | None,
    clonal_role: str,
) -> list[dict[str, str]]:
    output_gene = variant_row.get("Hugo_Symbol", "").strip() or canonical_gene_name(gene)
    genome_change = variant_row.get("Genome_Change", "").strip()
    variant_classification = variant_row.get("Variant_Classification", "").strip()
    block_label = f"{output_gene} mutation found"

    rows = [
        build_observation_row(base_row, GENE_CODE, block_label, output_gene, "text"),
        build_observation_row(
            base_row,
            GENOMIC_CHANGE_CODE,
            f"{block_label} - genomic DNA change (gHGVS)",
            genome_change,
            "text",
        ),
        build_observation_row(
            base_row,
            LOCAL_VARIANT_CLASS_CODE,
            f"{block_label} - variant classification",
            variant_classification,
            "text",
        ),
        build_observation_row(
            base_row,
            VAF_CODE,
            f"{block_label} - sample variant allelic frequency",
            format_vaf_value(vaf_pct),
            "numeric",
        ),
        build_observation_row(
            base_row,
            CLONAL_ROLE_CODE,
            f"{block_label}{CLONAL_ROLE_DESCRIPTION_SUFFIX}",
            clonal_role,
            "text",
        ),
    ]

    return rows


def is_shadowed_alias_gene(gene: str, genes_for_date: dict[str, float]) -> bool:
    for alias, canonical in GENE_ALIASES.items():
        if gene == canonical and alias in genes_for_date:
            return True
    return False


def transform_observations(
    observations_path: Path,
    baseline_genes_to_remove: dict[str, set[str]],
    dated_genes_to_remove: dict[str, dict[str, set[str]]],
    gene_vafs: dict[str, dict[str, dict[str, float]]],
    founding_genes: dict[str, set[str]],
    driver_variants: dict[str, list[dict[str, str]]],
    non_disruptive_variants: dict[str, list[dict[str, str]]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[str], Counter]:
    kept_rows: list[dict[str, str]] = []
    removed_rows: list[dict[str, str]] = []
    stats = Counter()
    patient_variant_cache: dict[tuple[str, str], tuple[dict[str, str], str]] = {}
    post_treatment_templates: dict[tuple[str, str, str], dict[str, str]] = {}
    post_treatment_seen_genes: dict[tuple[str, str, str], set[str]] = defaultdict(set)

    with observations_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        for row in reader:
            if not is_genomic_source_row(row):
                kept_rows.append(row)
                continue

            description = row["DESCRIPTION"].strip()
            if description in {HER2_IHC_DESCRIPTION, HER2_FISH_DESCRIPTION}:
                kept_rows.append(row)
                stats["her2_baseline_rows_kept"] += 1

                if description == HER2_FISH_DESCRIPTION or row["VALUE"].strip() != HER2_POSITIVE_VALUE:
                    continue

                gene = "HER2"
                if row_should_be_removed(row, baseline_genes_to_remove, dated_genes_to_remove):
                    continue

                patient_id = row["PATIENT"].strip()
                vaf_pct = lookup_gene_vaf(row, gene, gene_vafs)
                if vaf_pct is None:
                    stats["her2_variant_groups_skipped_without_vaf"] += 1
                    continue
                cache_key = (patient_id, gene)
                if cache_key not in patient_variant_cache:
                    patient_variant_cache[cache_key] = choose_patient_gene_variant(
                        patient_id,
                        gene,
                        driver_variants,
                        non_disruptive_variants,
                    )
                variant_row, variant_source = patient_variant_cache[cache_key]
                clonal_role = (
                    FOUNDING_ROLE_VALUE
                    if gene in founding_genes.get(patient_id, set())
                    else SUBCLONAL_ROLE_VALUE
                )
                group_rows = build_variant_observation_group(
                    row, gene, variant_row, vaf_pct, clonal_role
                )
                kept_rows.extend(group_rows)
                stats["variant_groups_emitted"] += 1
                stats["standardized_variant_rows_added"] += len(group_rows)
                stats[f"variant_groups_emitted_{variant_source}"] += 1
                if vaf_pct is None:
                    stats["variant_groups_missing_vaf"] += 1
                continue

            gene = source_row_to_gene(row)

            if gene is None:
                removed_rows.append(row)
                stats["genomic_source_rows_suppressed"] += 1
                continue

            if description.startswith(POST_TREATMENT_PREFIX):
                encounter_key = (
                    row["PATIENT"].strip(),
                    row["DATE"].strip(),
                    row["ENCOUNTER"].strip(),
                )
                post_treatment_templates.setdefault(encounter_key, dict(row))
                post_treatment_seen_genes[encounter_key].add(gene)

            if row_should_be_removed(row, baseline_genes_to_remove, dated_genes_to_remove):
                removed_rows.append(row)
                if description == GENETIC_ASSESSMENT_DESCRIPTION:
                    stats["baseline_rows_removed"] += 1
                elif description.startswith(POST_TREATMENT_PREFIX):
                    stats["post_treatment_rows_removed"] += 1
                elif description == HER2_IHC_DESCRIPTION:
                    stats["her2_rows_removed"] += 1
                continue

            patient_id = row["PATIENT"].strip()
            vaf_pct = lookup_gene_vaf(row, gene, gene_vafs)
            if vaf_pct is None:
                removed_rows.append(row)
                stats["genomic_rows_suppressed_without_clone_vaf"] += 1
                continue

            cache_key = (patient_id, gene)
            if cache_key not in patient_variant_cache:
                patient_variant_cache[cache_key] = choose_patient_gene_variant(
                    patient_id,
                    gene,
                    driver_variants,
                    non_disruptive_variants,
                )
            variant_row, variant_source = patient_variant_cache[cache_key]

            clonal_role = (
                FOUNDING_ROLE_VALUE
                if gene in founding_genes.get(patient_id, set())
                else SUBCLONAL_ROLE_VALUE
            )
            group_rows = build_variant_observation_group(
                row, gene, variant_row, vaf_pct, clonal_role
            )
            kept_rows.extend(group_rows)
            stats["variant_groups_emitted"] += 1
            stats["standardized_variant_rows_added"] += len(group_rows)
            stats[f"variant_groups_emitted_{variant_source}"] += 1
            if vaf_pct is None:
                stats["variant_groups_missing_vaf"] += 1

        for encounter_key, template_row in sorted(post_treatment_templates.items()):
            patient_id, date, _ = encounter_key
            genes_for_date = gene_vafs.get(patient_id, {}).get(date, {})
            if not genes_for_date:
                continue

            seen_genes = post_treatment_seen_genes.get(encounter_key, set())
            for gene in sorted(genes_for_date):
                if gene in seen_genes:
                    continue
                if is_shadowed_alias_gene(gene, genes_for_date):
                    continue
                if gene in dated_genes_to_remove.get(patient_id, {}).get(date, set()):
                    continue

                cache_key = (patient_id, gene)
                if cache_key not in patient_variant_cache:
                    patient_variant_cache[cache_key] = choose_patient_gene_variant(
                        patient_id,
                        gene,
                        driver_variants,
                        non_disruptive_variants,
                    )
                variant_row, variant_source = patient_variant_cache[cache_key]
                vaf_pct = genes_for_date.get(gene)
                if vaf_pct is None:
                    continue

                synthetic_row = dict(template_row)
                synthetic_row["DESCRIPTION"] = f"{POST_TREATMENT_PREFIX}{gene}"
                clonal_role = (
                    FOUNDING_ROLE_VALUE
                    if gene in founding_genes.get(patient_id, set())
                    else SUBCLONAL_ROLE_VALUE
                )
                group_rows = build_variant_observation_group(
                    synthetic_row, gene, variant_row, vaf_pct, clonal_role
                )
                kept_rows.extend(group_rows)
                stats["variant_groups_emitted"] += 1
                stats["standardized_variant_rows_added"] += len(group_rows)
                stats[f"variant_groups_emitted_{variant_source}"] += 1
                stats["synthetic_post_treatment_variant_groups_added"] += 1

        stats["total_rows_removed"] = len(removed_rows)
        stats["total_rows_kept"] = len(kept_rows)
        return kept_rows, removed_rows, fieldnames, stats


def write_rows(output_path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]

    clone_proportions_path = args.clone_proportions or newest_clone_proportions_csv(repo_root)
    if not clone_proportions_path.exists():
        raise FileNotFoundError(f"Clone proportions file not found: {clone_proportions_path}")

    observations_path = args.observations or clone_proportions_path.with_name("observations.csv")
    if not observations_path.exists():
        if args.observations is None:
            observations_path = newest_observations_csv(repo_root)
        if not observations_path.exists():
            raise FileNotFoundError(f"Observations file not found: {observations_path}")

    output_path = args.output or observations_path.with_name(DEFAULT_OUTPUT_NAME)
    driver_variants_path = args.driver_variants or (
        repo_root / "scripts" / DEFAULT_DRIVER_VARIANTS_NAME
    )
    non_disruptive_variants_path = args.non_disruptive_variants or (
        repo_root / "scripts" / DEFAULT_NON_DISRUPTIVE_VARIANTS_NAME
    )
    if not driver_variants_path.exists():
        raise FileNotFoundError(f"Driver variant file not found: {driver_variants_path}")
    if not non_disruptive_variants_path.exists():
        raise FileNotFoundError(
            f"Non-disruptive variant file not found: {non_disruptive_variants_path}"
        )

    baseline_genes_to_remove, dated_genes_to_remove, target_stats = load_low_vaf_targets(
        clone_proportions_path, args.threshold
    )
    gene_vafs = load_clone_vaf_map(clone_proportions_path)
    founding_genes = load_founding_gene_map(clone_proportions_path)
    driver_variants = load_variant_rows(driver_variants_path)
    non_disruptive_variants = load_variant_rows(non_disruptive_variants_path)
    kept_rows, removed_rows, fieldnames, transform_stats = transform_observations(
        observations_path,
        baseline_genes_to_remove,
        dated_genes_to_remove,
        gene_vafs,
        founding_genes,
        driver_variants,
        non_disruptive_variants,
    )
    write_rows(output_path, fieldnames, kept_rows)

    print(f"Clone proportions: {clone_proportions_path}")
    print(f"Input observations: {observations_path}")
    print(f"Driver variants: {driver_variants_path}")
    print(f"Non-disruptive variants: {non_disruptive_variants_path}")
    print(f"Wrote standardized observations: {output_path}")
    print(f"Threshold: {args.threshold:.2f}%")
    print(f"Patients with baseline removals: {len(baseline_genes_to_remove)}")
    print(f"Patients with dated removals: {len(dated_genes_to_remove)}")
    print(
        "Clone rows below threshold: "
        f"{sum(target_stats[key] for key in target_stats if key.endswith('_clone_rows_below_threshold'))}"
    )
    print(f"Rows removed from source observations: {transform_stats['total_rows_removed']}")
    print(f"Baseline mutation rows removed: {transform_stats['baseline_rows_removed']}")
    print(f"Post-treatment trend rows removed: {transform_stats['post_treatment_rows_removed']}")
    print(f"HER2 rows removed: {transform_stats['her2_rows_removed']}")
    print(f"HER2 baseline rows kept unchanged: {transform_stats['her2_baseline_rows_kept']}")
    print(f"Non-positive or FISH genomic rows suppressed: {transform_stats['genomic_source_rows_suppressed']}")
    print(
        "Genomic rows suppressed without clone VAF: "
        f"{transform_stats['genomic_rows_suppressed_without_clone_vaf']}"
    )
    print(
        "HER2 variant groups skipped without clone VAF: "
        f"{transform_stats['her2_variant_groups_skipped_without_vaf']}"
    )
    print(f"Variant groups emitted: {transform_stats['variant_groups_emitted']}")
    print(
        "Variant groups emitted from driver file: "
        f"{transform_stats['variant_groups_emitted_driver']}"
    )
    print(
        "Variant groups emitted from non-disruptive file: "
        f"{transform_stats['variant_groups_emitted_non_disruptive']}"
    )
    print(
        "Variant groups emitted without a matching variant row: "
        f"{transform_stats['variant_groups_emitted_missing']}"
    )
    print(f"Variant groups missing VAF: {transform_stats['variant_groups_missing_vaf']}")
    print(f"Standardized variant rows added: {transform_stats['standardized_variant_rows_added']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
