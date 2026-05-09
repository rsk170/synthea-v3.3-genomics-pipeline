#!/usr/bin/env python3
"""Create per-patient, per-clone complete MAF files.

Each output MAF combines:
- driver variants for one reconstructed clone, using clone membership from
  breast_cancer_clone_groups.csv and exact variant annotations recovered from
  observations_pruned_by_clone_vaf.csv
- a non-overlapping subset of passenger variants assigned to that patient from
  breast_cancer_assigned_passenger_mutations.tsv

Driver rows are recovered from the exact ICGC-backed source libraries using
gene + genomic change + variant classification keys. Passenger rows are split
randomly, but reproducibly, across the patient's clones so that each passenger
mutation appears in at most one clone-level MAF file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import random
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_PRUNED_NAME = "observations_pruned_by_clone_vaf.csv"
DEFAULT_ASSIGNED_PASSENGERS_NAME = "breast_cancer_assigned_passenger_mutations.tsv"
DEFAULT_CLONE_GROUPS_NAME = "breast_cancer_clone_groups.csv"
DEFAULT_CLONE_PROPORTIONS_NAME = "breast_cancer_clone_proportions.csv"
DEFAULT_CIVIC_DRIVER_VARIANTS_NAME = "civic_breast_cancer_driver_variants_from_maf.csv"
DEFAULT_DRIVER_VARIANTS_NAME = "breast_cancer_driver_variants_from_maf.csv"
DEFAULT_NON_DISRUPTIVE_VARIANTS_NAME = "breast_cancer_non_disruptive_variants_from_maf.csv"
DEFAULT_MAF_NAME = "final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf"

GENE_CODE = "48018-6"
GENOMIC_CHANGE_CODE = "81290-9"
VARIANT_CLASS_CODE = "83005-9"
GENE_FOUND_SUFFIX = " mutation found"
GENOMIC_CHANGE_SUFFIX = " mutation found - genomic DNA change (gHGVS)"
VARIANT_CLASS_SUFFIX = " mutation found - variant classification"
GENE_ALIASES = {"HER2": "ERBB2"}

ASSIGNED_PATIENT_ID = "assigned_patient_id"
MUTATION_ORIGIN = "mutation_origin"
SOURCE_TUMOR_SAMPLE_BARCODE = "source_tumor_sample_barcode"
SOURCE_DONOR_ID = "source_donor_id"

HELPER_FIELDS = [
    MUTATION_ORIGIN,
    SOURCE_TUMOR_SAMPLE_BARCODE,
    SOURCE_DONOR_ID,
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create per-patient complete MAF files for each reconstructed clone by "
            "combining clone-specific driver variants from "
            "observations_pruned_by_clone_vaf.csv with a non-overlapping split of "
            "assigned passenger mutations."
        )
    )
    parser.add_argument(
        "--pruned-observations",
        type=Path,
        help=(
            "Path to observations_pruned_by_clone_vaf.csv. Defaults to the newest "
            "output_runs/*/csv/observations_pruned_by_clone_vaf.csv."
        ),
    )
    parser.add_argument(
        "--assigned-passengers",
        type=Path,
        help=(
            "Path to breast_cancer_assigned_passenger_mutations.tsv. Defaults to "
            "the same directory as the pruned observations file."
        ),
    )
    parser.add_argument(
        "--clone-groups",
        type=Path,
        help=(
            "Path to breast_cancer_clone_groups.csv. Defaults to the same "
            "directory as the pruned observations file."
        ),
    )
    parser.add_argument(
        "--clone-proportions",
        type=Path,
        help=(
            "Path to breast_cancer_clone_proportions.csv. Defaults to the same "
            "directory as the pruned observations file."
        ),
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
        "--civic-driver-variants",
        type=Path,
        help=(
            "Path to civic_breast_cancer_driver_variants_from_maf.csv. Defaults "
            "to scripts/civic_breast_cancer_driver_variants_from_maf.csv."
        ),
    )
    parser.add_argument(
        "--non-disruptive-variants",
        type=Path,
        help=(
            "Path to breast_cancer_non_disruptive_variants_from_maf.csv. Defaults "
            "to scripts/breast_cancer_non_disruptive_variants_from_maf.csv."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Directory where maf_files/ should be written. Defaults to a maf_files "
            "subdirectory under the run root (sibling of csv/)."
        ),
    )
    parser.add_argument(
        "--maf",
        type=Path,
        help=(
            "Path to the full ICGC MAF for exact fallback recovery when a driver "
            "variant is not present in the curated variant tables. Defaults to "
            "scripts/final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf."
        ),
    )
    return parser.parse_args()


def newest_pruned_observations_csv(repo_root: Path) -> Path:
    candidates = sorted(
        repo_root.glob(f"output_runs/*/csv/{DEFAULT_PRUNED_NAME}"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(f"No {DEFAULT_PRUNED_NAME} found under output_runs/*/csv/")
    return candidates[0]


def canonical_gene_name(gene: str) -> str:
    return GENE_ALIASES.get(gene, gene)


def extract_gene_from_description(description: str, suffix: str) -> str | None:
    if description.endswith(suffix):
        return description[: -len(suffix)].strip()
    return None


def derive_maf_fieldnames(fieldnames: list[str]) -> list[str]:
    if "Hugo_Symbol" not in fieldnames or "Donor_ID" not in fieldnames:
        raise ValueError("Expected Hugo_Symbol and Donor_ID in MAF-like source file header")
    start = fieldnames.index("Hugo_Symbol")
    end = fieldnames.index("Donor_ID")
    return fieldnames[start : end + 1]


def variant_sort_key(row: dict[str, str]) -> tuple:
    project = row.get("Project_Code", "") or row.get("helper_project_code", "")
    sample = row.get("Tumor_Sample_Barcode", "") or row.get("sample_id", "")
    donor = row.get("Donor_ID", "")
    helper_key = row.get("helper_variant_key", "")
    return (
        0 if "breast" in project.lower() else 1,
        project,
        sample,
        donor,
        helper_key,
    )


def stable_seed(*parts: str) -> int:
    digest = hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()
    return int(digest, 16)


def clone_sort_key(clone_id: str) -> tuple[int, str]:
    prefix, _, suffix = clone_id.partition("_")
    if prefix == "clone" and suffix.isdigit():
        return (int(suffix), clone_id)
    return (10**9, clone_id)


def load_variant_index(
    path: Path, source_name: str
) -> tuple[dict[tuple[str, str, str], list[dict[str, str]]], list[str]]:
    index: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        maf_fieldnames = derive_maf_fieldnames(fieldnames)
        for row in reader:
            gene = canonical_gene_name(
                row.get("helper_gene", "").strip() or row.get("Hugo_Symbol", "").strip()
            )
            genomic_change = row.get("Genome_Change", "").strip()
            variant_class = row.get("Variant_Classification", "").strip()
            if not gene or not genomic_change or not variant_class:
                continue
            normalized = dict(row)
            normalized["_source_variant_library"] = source_name
            index[(gene, genomic_change, variant_class)].append(normalized)
    return index, maf_fieldnames


def load_assigned_passengers(
    path: Path,
) -> tuple[dict[str, list[dict[str, str]]], list[str]]:
    passengers_by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        maf_fieldnames = derive_maf_fieldnames(fieldnames)
        for row in reader:
            patient_id = row.get(ASSIGNED_PATIENT_ID, "").strip()
            if not patient_id:
                continue
            passengers_by_patient[patient_id].append(dict(row))
    return passengers_by_patient, maf_fieldnames


def load_pruned_driver_variants(
    path: Path,
) -> dict[str, dict[str, dict[str, str]]]:
    variant_state: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            code = row["CODE"].strip()
            description = row["DESCRIPTION"].strip()
            patient_id = row["PATIENT"].strip()

            if code not in {GENE_CODE, GENOMIC_CHANGE_CODE, VARIANT_CLASS_CODE}:
                continue

            if code == GENE_CODE:
                gene = row["VALUE"].strip() or extract_gene_from_description(description, GENE_FOUND_SUFFIX)
            elif code == GENOMIC_CHANGE_CODE:
                gene = extract_gene_from_description(description, GENOMIC_CHANGE_SUFFIX)
            else:
                gene = extract_gene_from_description(description, VARIANT_CLASS_SUFFIX)

            if not gene:
                continue

            canonical_gene = canonical_gene_name(gene)
            record = variant_state[patient_id].setdefault(
                canonical_gene,
                {
                    "event_gene_name": gene,
                    "genomic_changes": Counter(),
                    "variant_classifications": Counter(),
                },
            )
            if code == GENE_CODE:
                if not record["event_gene_name"]:
                    record["event_gene_name"] = gene
            elif code == GENOMIC_CHANGE_CODE and row["VALUE"].strip():
                record["genomic_changes"][row["VALUE"].strip()] += 1
            elif code == VARIANT_CLASS_CODE and row["VALUE"].strip():
                record["variant_classifications"][row["VALUE"].strip()] += 1

    variants_by_patient: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for patient_id, patient_records in variant_state.items():
        for gene, record in patient_records.items():
            genomic_changes: Counter = record["genomic_changes"]  # type: ignore[assignment]
            variant_classes: Counter = record["variant_classifications"]  # type: ignore[assignment]
            variants_by_patient[patient_id][gene] = {
                "gene": gene,
                "event_gene_name": str(record["event_gene_name"] or gene),
                "genomic_change": genomic_changes.most_common(1)[0][0] if genomic_changes else "",
                "variant_classification": variant_classes.most_common(1)[0][0]
                if variant_classes
                else "",
            }
    return variants_by_patient


def load_clone_groups(path: Path) -> dict[str, list[dict[str, object]]]:
    clones_by_patient: dict[str, list[dict[str, object]]] = defaultdict(list)
    genes_seen_by_patient: dict[str, dict[str, str]] = defaultdict(dict)

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            patient_id = row["patient_id"].strip()
            clone_id = row["clone_id"].strip()
            genes = [
                canonical_gene_name(gene.strip())
                for gene in row.get("genes", "").split(";")
                if gene.strip()
            ]
            for gene in genes:
                existing_clone = genes_seen_by_patient[patient_id].get(gene)
                if existing_clone and existing_clone != clone_id:
                    raise ValueError(
                        f"Gene {gene} appears in multiple clones for patient {patient_id}: "
                        f"{existing_clone} and {clone_id}"
                    )
                genes_seen_by_patient[patient_id][gene] = clone_id

            clones_by_patient[patient_id].append(
                {
                    "clone_id": clone_id,
                    "clone_type": row.get("clone_type", "").strip(),
                    "parent_clone_id": row.get("parent_clone_id", "").strip(),
                    "signature_key": row.get("signature_key", "").strip(),
                    "timepoint_dates": row.get("timepoint_dates", "").strip(),
                    "genes": genes,
                }
            )

    for patient_id in clones_by_patient:
        clones_by_patient[patient_id].sort(
            key=lambda clone: clone_sort_key(str(clone["clone_id"]))
        )
    return clones_by_patient


def load_clone_proportion_metadata(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    metadata: dict[tuple[str, str], dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            patient_id = row["patient_id"].strip()
            clone_id = row["clone_id"].strip()
            dates: list[str] = []
            vaf_parts: list[str] = []
            for index in range(3):
                date = row.get(f"t{index}_date", "").strip()
                vaf = row.get(f"t{index}_vaf_pct", "").strip()
                if date:
                    dates.append(date)
                if date or vaf:
                    label = f"t{index}"
                    if date and vaf:
                        vaf_parts.append(f"{label}:{date}={vaf}")
                    elif vaf:
                        vaf_parts.append(f"{label}:{vaf}")
                    else:
                        vaf_parts.append(f"{label}:{date}")
            metadata[(patient_id, clone_id)] = {
                "clone_timepoint_dates": ";".join(dates),
                "clone_timepoint_vaf_pct": ";".join(vaf_parts),
            }
    return metadata


def collect_missing_keys_from_clones(
    clones_by_patient: dict[str, list[dict[str, object]]],
    pruned_variants_by_patient: dict[str, dict[str, dict[str, str]]],
    indexes: list[dict[tuple[str, str, str], list[dict[str, str]]]],
) -> set[tuple[str, str, str]]:
    missing: set[tuple[str, str, str]] = set()
    for patient_id, clones in clones_by_patient.items():
        patient_variants = pruned_variants_by_patient.get(patient_id, {})
        for clone in clones:
            for gene in clone["genes"]:  # type: ignore[index]
                record = patient_variants.get(str(gene))
                if not record:
                    continue
                genomic_change = record["genomic_change"].strip()
                variant_classification = record["variant_classification"].strip()
                if not genomic_change or not variant_classification:
                    continue
                key = (canonical_gene_name(str(gene)), genomic_change, variant_classification)
                if any(key in index for index in indexes):
                    continue
                missing.add(key)
    return missing


def load_exact_maf_fallback(
    maf_path: Path,
    needed_keys: set[tuple[str, str, str]],
) -> tuple[dict[tuple[str, str, str], list[dict[str, str]]], list[str]]:
    fallback_index: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    if not needed_keys:
        return fallback_index, []

    with maf_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        maf_fieldnames = derive_maf_fieldnames(fieldnames)
        for row in reader:
            gene = canonical_gene_name(row.get("Hugo_Symbol", "").strip())
            genomic_change = row.get("Genome_Change", "").strip()
            variant_class = row.get("Variant_Classification", "").strip()
            key = (gene, genomic_change, variant_class)
            if key not in needed_keys:
                continue
            normalized = dict(row)
            normalized["_source_variant_library"] = "maf_exact_fallback"
            fallback_index[key].append(normalized)
    return fallback_index, maf_fieldnames


def choose_variant_row(
    gene: str,
    genomic_change: str,
    variant_classification: str,
    indexes: list[dict[tuple[str, str, str], list[dict[str, str]]]],
) -> dict[str, str] | None:
    key = (canonical_gene_name(gene), genomic_change, variant_classification)
    for index in indexes:
        matches = index.get(key)
        if matches:
            return sorted(matches, key=variant_sort_key)[0]
    return None


def split_passengers_across_clones(
    patient_id: str,
    clones: list[dict[str, object]],
    passenger_rows: list[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    passengers_by_clone = {str(clone["clone_id"]): [] for clone in clones}
    if not clones or not passenger_rows:
        return passengers_by_clone

    shuffled_rows = list(passenger_rows)
    rng = random.Random(stable_seed(patient_id, "clone_passenger_split"))
    rng.shuffle(shuffled_rows)

    ordered_clone_ids = [str(clone["clone_id"]) for clone in clones]
    for index, row in enumerate(shuffled_rows):
        clone_id = ordered_clone_ids[index % len(ordered_clone_ids)]
        passengers_by_clone[clone_id].append(row)

    return passengers_by_clone


def synthetic_maf_filename(clone_id: str, clone_type: str) -> str:
    safe_clone_type = clone_type or "clone"
    return f"{clone_id}_{safe_clone_type}.maf"


def build_driver_output_row(
    variant_row: dict[str, str],
    maf_fieldnames: list[str],
) -> dict[str, str]:
    output = {
        MUTATION_ORIGIN: "driver",
        SOURCE_TUMOR_SAMPLE_BARCODE: variant_row.get("Tumor_Sample_Barcode", ""),
        SOURCE_DONOR_ID: variant_row.get("Donor_ID", ""),
    }
    for field in maf_fieldnames:
        output[field] = variant_row.get(field, "")
    return output


def build_passenger_output_rows(
    passenger_rows: list[dict[str, str]],
    maf_fieldnames: list[str],
) -> list[dict[str, str]]:
    output_rows: list[dict[str, str]] = []
    for row in passenger_rows:
        output = {
            MUTATION_ORIGIN: "passenger",
            SOURCE_TUMOR_SAMPLE_BARCODE: row.get(
                SOURCE_TUMOR_SAMPLE_BARCODE, row.get("Tumor_Sample_Barcode", "")
            ),
            SOURCE_DONOR_ID: row.get(SOURCE_DONOR_ID, row.get("Donor_ID", "")),
        }
        for field in maf_fieldnames:
            output[field] = row.get(field, "")
        output_rows.append(output)
    return output_rows


def write_complete_maf(
    output_path: Path,
    fieldnames: list[str],
    driver_rows: list[dict[str, str]],
    passenger_rows: list[dict[str, str]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in driver_rows:
            writer.writerow(row)
        for row in passenger_rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]

    pruned_path = args.pruned_observations or newest_pruned_observations_csv(repo_root)
    assigned_passengers_path = args.assigned_passengers or pruned_path.with_name(
        DEFAULT_ASSIGNED_PASSENGERS_NAME
    )
    clone_groups_path = args.clone_groups or pruned_path.with_name(DEFAULT_CLONE_GROUPS_NAME)
    clone_proportions_path = args.clone_proportions or pruned_path.with_name(
        DEFAULT_CLONE_PROPORTIONS_NAME
    )
    driver_variants_path = args.driver_variants or (
        repo_root / "scripts" / DEFAULT_DRIVER_VARIANTS_NAME
    )
    civic_variants_path = args.civic_driver_variants or (
        repo_root / "scripts" / DEFAULT_CIVIC_DRIVER_VARIANTS_NAME
    )
    non_disruptive_variants_path = args.non_disruptive_variants or (
        repo_root / "scripts" / DEFAULT_NON_DISRUPTIVE_VARIANTS_NAME
    )
    maf_path = args.maf or (repo_root / "scripts" / DEFAULT_MAF_NAME)
    output_dir = args.output_dir or pruned_path.parent.parent / "maf_files"

    for required_path in [
        pruned_path,
        assigned_passengers_path,
        clone_groups_path,
        clone_proportions_path,
        driver_variants_path,
        civic_variants_path,
        non_disruptive_variants_path,
        maf_path,
    ]:
        if not required_path.exists():
            raise FileNotFoundError(f"Required file not found: {required_path}")

    driver_index, maf_fieldnames = load_variant_index(driver_variants_path, "driver")
    civic_index, civic_maf_fieldnames = load_variant_index(civic_variants_path, "civic")
    non_disruptive_index, non_disruptive_maf_fieldnames = load_variant_index(
        non_disruptive_variants_path, "non_disruptive"
    )
    passengers_by_patient, passenger_maf_fieldnames = load_assigned_passengers(
        assigned_passengers_path
    )
    pruned_variants_by_patient = load_pruned_driver_variants(pruned_path)
    clones_by_patient = load_clone_groups(clone_groups_path)
    clone_proportion_metadata = load_clone_proportion_metadata(clone_proportions_path)
    needed_maf_fallback_keys = collect_missing_keys_from_clones(
        clones_by_patient,
        pruned_variants_by_patient,
        [driver_index, civic_index, non_disruptive_index],
    )
    maf_fallback_index, fallback_maf_fieldnames = load_exact_maf_fallback(
        maf_path, needed_maf_fallback_keys
    )

    if maf_fieldnames != civic_maf_fieldnames or maf_fieldnames != non_disruptive_maf_fieldnames:
        raise ValueError("Driver source libraries do not share the same MAF column layout")
    if maf_fieldnames != passenger_maf_fieldnames:
        raise ValueError("Passenger file MAF columns do not match driver source columns")
    if fallback_maf_fieldnames and maf_fieldnames != fallback_maf_fieldnames:
        raise ValueError("Raw MAF columns do not match driver source columns")

    fieldnames = [*HELPER_FIELDS, *maf_fieldnames]
    stats = Counter()

    for patient_id in sorted(clones_by_patient):
        clones = clones_by_patient[patient_id]
        patient_variants = pruned_variants_by_patient.get(patient_id, {})
        passenger_rows = passengers_by_patient.get(patient_id, [])
        passenger_rows_by_clone = split_passengers_across_clones(
            patient_id, clones, passenger_rows
        )

        if not passenger_rows:
            stats["patients_with_no_passenger_rows"] += 1

        patient_dir = output_dir / patient_id

        for clone in clones:
            clone_id = str(clone["clone_id"])
            clone["clone_timepoint_vaf_pct"] = clone_proportion_metadata.get(
                (patient_id, clone_id), {}
            ).get("clone_timepoint_vaf_pct", "")
            if not clone.get("timepoint_dates"):
                clone["timepoint_dates"] = clone_proportion_metadata.get(
                    (patient_id, clone_id), {}
                ).get("clone_timepoint_dates", "")

            driver_rows: list[dict[str, str]] = []
            for gene in clone["genes"]:  # type: ignore[index]
                record = patient_variants.get(str(gene))
                if record is None:
                    stats["clone_driver_genes_missing_from_pruned_observations"] += 1
                    continue

                genomic_change = record["genomic_change"].strip()
                variant_classification = record["variant_classification"].strip()
                if not genomic_change or not variant_classification:
                    stats["driver_records_skipped_missing_fields"] += 1
                    continue

                variant_row = choose_variant_row(
                    str(gene),
                    genomic_change,
                    variant_classification,
                    [driver_index, civic_index, non_disruptive_index, maf_fallback_index],
                )
                if variant_row is None:
                    stats["driver_records_missing_exact_match"] += 1
                    continue

                driver_rows.append(
                    build_driver_output_row(variant_row, maf_fieldnames)
                )
                stats["driver_rows_written"] += 1

            passenger_output_rows = build_passenger_output_rows(
                passenger_rows_by_clone.get(clone_id, []),
                maf_fieldnames,
            )

            output_path = patient_dir / synthetic_maf_filename(
                clone_id, str(clone["clone_type"])
            )
            write_complete_maf(output_path, fieldnames, driver_rows, passenger_output_rows)
            stats["maf_files_written"] += 1
            stats["passenger_rows_written"] += len(passenger_output_rows)

        if clones:
            stats["patients_written"] += 1
            stats["clones_written"] += len(clones)

    print(f"Pruned observations: {pruned_path}")
    print(f"Clone groups: {clone_groups_path}")
    print(f"Clone proportions: {clone_proportions_path}")
    print(f"Assigned passengers: {assigned_passengers_path}")
    print(f"Raw MAF exact fallback source: {maf_path}")
    print(f"Output directory: {output_dir}")
    print(f"Patients with MAF output: {stats['patients_written']}")
    print(f"Patients with no passenger rows: {stats['patients_with_no_passenger_rows']}")
    print(f"Clone MAF files written: {stats['maf_files_written']}")
    print(f"Clones written: {stats['clones_written']}")
    print(f"Driver rows written: {stats['driver_rows_written']}")
    print(f"Passenger rows written: {stats['passenger_rows_written']}")
    print(
        "Clone driver genes missing from pruned observations: "
        f"{stats['clone_driver_genes_missing_from_pruned_observations']}"
    )
    print(
        "Driver rows skipped (missing genomic change or variant classification): "
        f"{stats['driver_records_skipped_missing_fields']}"
    )
    print(
        "Driver rows missing exact source match: "
        f"{stats['driver_records_missing_exact_match']}"
    )
    print(f"Exact-fallback keys requested from raw MAF: {len(needed_maf_fallback_keys)}")


if __name__ == "__main__":
    main()
