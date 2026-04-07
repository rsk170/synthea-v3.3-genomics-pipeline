#!/usr/bin/env python3
"""Build minimal breast cancer clone groupings from Synthea observations.csv.

The script reconstructs per-patient clone groupings from:
- baseline genomic profile genes from "Genetic variant assessment"
- baseline HER2 positivity from HER2 IHC observations
- later post-treatment clone-trend observations

Clone construction rules:
- all genes start at baseline state "present"
- likely non-treatment random drivers seed early founding/branch clones
- later informative genes are grouped into late clones by exact full signature
- output rows represent new genes added at each clone plus the cumulative genotype
- patients with no post-treatment clone-trend rows are skipped
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path


GENETIC_ASSESSMENT_DESCRIPTION = "Genetic variant assessment"
HER2_IHC_DESCRIPTION = "HER2 [Presence] in Breast cancer specimen by Immune stain"
HER2_POSITIVE_VALUE = "Positive (qualifier value)"
ER_IHC_DESCRIPTION = "Estrogen receptor Ag [Presence] in Breast cancer specimen by Immune stain"
PR_IHC_DESCRIPTION = "Progesterone receptor Ag [Presence] in Breast cancer specimen by Immune stain"
STAGE_DESCRIPTION = "Stage group.clinical Cancer"
RESPONSE_DESCRIPTION = "Response to cancer treatment"
IMPROVING_VALUE = "Improving (qualifier value)"
WORSENING_VALUE = "Worsening (qualifier value)"
POST_TREATMENT_PREFIX = "Post-treatment clone trend "
BASELINE_POSITIVE_SUFFIX = " mutation positive"
DEFAULT_OUTPUT_NAME = "breast_cancer_clone_groups.csv"
ADVANCED_STAGES = {
    "Stage 3 (qualifier value)",
    "Stage 3A (qualifier value)",
    "Stage 3B (qualifier value)",
    "Stage 3C (qualifier value)",
    "Stage 4 (qualifier value)",
}
TREATMENT_BUDGET_BY_TARGET = {
    1: 1,
    2: 1,
    3: 2,
    4: 3,
    5: 4,
    6: 4,
    7: 5,
    8: 6,
    9: 6,
    10: 7,
    11: 8,
    12: 8,
    13: 9,
    14: 10,
    15: 11,
    16: 11,
}
GROUP_ORDER = ["chemo", "endocrine", "her2", "cdk46"]
SENSITIVE_POOLS = {
    "chemo": [
        ["CYP19A1", "BIRC5", "ABCB1", "CYP1B1", "SLCO1B1"],
        ["CYP19A1", "BIRC5", "ABCB1", "CYP1B1", "SLCO1B1"],
    ],
    "endocrine": [["PGR"], []],
    "her2": [
        ["BRAF", "PIK3CA", "ESR1", "PGR", "ESR2", "ERBB3"],
        ["BRAF", "PIK3CA", "ESR1", "PGR", "ESR2", "ERBB3"],
    ],
    "cdk46": [
        ["ESR1", "CDKN2A", "CCND1", "PGR", "ESR2", "NRAS"],
        ["ESR1", "CDKN2A", "CCND1", "PGR", "ESR2", "NRAS"],
    ],
}
RESISTANCE_POOLS = {
    "chemo": ["NCOR2"],
    "endocrine": ["NF1"],
    "her2": ["MET"],
    "cdk46": ["RB1"],
}
NONINFORMATIVE_STATES = {"unknown", "missing"}
FOUNDER_TIER_GENES = ["TP53", "PIK3CA", "GATA3", "CDH1", "PTEN"]
CONTEXT_FOUNDER_GENES = ["MAP3K1", "CCND1", "RB1", "KMT2C", "CBFB", "CDKN2A", "NF1", "BRCA2"]
FALLBACK_PRIORITY_GENES = [
    "BRCA2",
    "MAP2K4",
    "NCOR1",
    "NOTCH2",
    "PARK2",
    "PRKCI",
    "MECOM",
    "FOXQ1",
    "YWHAZ",
    "MYC",
    "ZNF703",
    "MCL1",
    "LINC00290",
    "PGR",
]
AMPLIFICATION_LIKE_GENES = {"YWHAZ", "MYC", "ZNF703", "MCL1", "CCND1"}
RECEPTOR_MARKER_GENES = {"PGR"}
TREATMENT_RELATED_GENES = set(
    gene
    for pools in list(SENSITIVE_POOLS.values()) + [[pool] for pool in RESISTANCE_POOLS.values()]
    for pool in pools
    for gene in pool
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create minimal per-patient breast cancer clone groupings from "
            "Synthea observations.csv."
        )
    )
    parser.add_argument(
        "--observations",
        type=Path,
        help="Path to observations.csv. Defaults to the newest output_runs/*/csv/observations.csv.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Output CSV path. Defaults to breast_cancer_clone_groups.csv in the "
            "same directory as the input observations.csv."
        ),
    )
    return parser.parse_args()


def newest_observations_csv(repo_root: Path) -> Path:
    candidates = sorted(
        repo_root.glob("output_runs/*/csv/observations.csv"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError("No observations.csv found under output_runs/*/csv/")
    return candidates[0]


def extract_baseline_gene(value: str) -> str | None:
    if value.endswith(BASELINE_POSITIVE_SUFFIX):
        return value[: -len(BASELINE_POSITIVE_SUFFIX)]
    return None


def signature_key(signature: tuple[str, ...]) -> str:
    return "|".join(f"t{i}:{state}" for i, state in enumerate(signature))


def mixed_signature_key(signatures: list[tuple[str, ...]]) -> str:
    unique = {sig for sig in signatures}
    if len(unique) == 1:
        return signature_key(signatures[0])
    return "mixed:multiple-signatures"


def trend_priority(later_states: list[str]) -> int:
    if any(state in {"increasing", "decreasing"} for state in later_states):
        return 0
    if any(state == "stable" for state in later_states):
        return 1
    return 2


def treatment_budget(total_driver_genes: int) -> int:
    return TREATMENT_BUDGET_BY_TARGET.get(total_driver_genes, 0)


def infer_survival(patient_data: dict) -> str:
    if WORSENING_VALUE in patient_data["responses"]:
        return "no"
    if IMPROVING_VALUE in patient_data["responses"]:
        return "yes"
    return "yes"


def eligible_groups(patient_data: dict, survival: str) -> set[str]:
    groups: set[str] = set()
    er_positive = patient_data["er_status"] == HER2_POSITIVE_VALUE or patient_data["er_status"] == "Positive (qualifier value)"
    pr_positive = patient_data["pr_status"] == "Positive (qualifier value)"
    her2_positive = patient_data["her2_status"] == HER2_POSITIVE_VALUE
    triple_negative = (
        patient_data["her2_status"] != HER2_POSITIVE_VALUE
        and patient_data["er_status"] != "Positive (qualifier value)"
        and patient_data["pr_status"] != "Positive (qualifier value)"
    )
    advanced_stage = any(stage in ADVANCED_STAGES for stage in patient_data["stage_values"])

    if triple_negative or advanced_stage:
        groups.add("chemo")
    if er_positive or pr_positive:
        groups.add("endocrine")
    if her2_positive:
        groups.add("her2")
    if er_positive and pr_positive:
        groups.add("cdk46")
    if survival == "no":
        return groups
    return groups


def choose_gene_from_pool(
    pool: list[str],
    available_genes: set[str],
    gene_trends: dict[str, dict[str, str]],
    trend_dates: list[str],
) -> str | None:
    candidates = [gene for gene in pool if gene in available_genes]
    if not candidates:
        return None
    candidates.sort(
        key=lambda gene: (
            trend_priority([gene_trends.get(gene, {}).get(date, "missing") for date in trend_dates]),
            pool.index(gene),
            gene,
        )
    )
    return candidates[0]


def infer_random_driver_genes(patient_data: dict, trend_dates: list[str]) -> set[str]:
    driver_genes = set(patient_data["baseline_driver_genes"])
    if not driver_genes:
        return set()

    budget = treatment_budget(len(driver_genes))
    if budget <= 0:
        return driver_genes

    available_genes = set(driver_genes)
    selected_treatment_genes: list[str] = []
    survival = infer_survival(patient_data)
    groups = eligible_groups(patient_data, survival)

    if survival == "yes":
        for pass_index in (0, 1):
            for group in GROUP_ORDER:
                if len(selected_treatment_genes) >= budget:
                    break
                if group not in groups:
                    continue
                pool = SENSITIVE_POOLS[group][pass_index]
                if not pool:
                    continue
                gene = choose_gene_from_pool(
                    pool, available_genes, patient_data["gene_trends"], trend_dates
                )
                if gene is None:
                    continue
                selected_treatment_genes.append(gene)
                available_genes.remove(gene)
            if len(selected_treatment_genes) >= budget:
                break
    else:
        for group in GROUP_ORDER:
            if len(selected_treatment_genes) >= budget:
                break
            if group not in groups:
                continue
            pool = RESISTANCE_POOLS[group]
            gene = choose_gene_from_pool(pool, available_genes, patient_data["gene_trends"], trend_dates)
            if gene is None:
                continue
            selected_treatment_genes.append(gene)
            available_genes.remove(gene)

    return available_genes


def subtype_flags(patient_data: dict) -> dict[str, bool]:
    er_positive = patient_data["er_status"] == "Positive (qualifier value)"
    pr_positive = patient_data["pr_status"] == "Positive (qualifier value)"
    her2_positive = patient_data["her2_status"] == HER2_POSITIVE_VALUE
    aggressive = (
        (not er_positive and not pr_positive and not her2_positive)
        or any(stage in ADVANCED_STAGES for stage in patient_data["stage_values"])
    )
    luminal_like = er_positive or pr_positive
    lobular_like = "CDH1" in patient_data["baseline_driver_genes"]
    return {
        "er_positive": er_positive,
        "pr_positive": pr_positive,
        "her2_positive": her2_positive,
        "aggressive": aggressive,
        "luminal_like": luminal_like,
        "lobular_like": lobular_like,
    }


def founder_category(gene: str) -> int:
    if gene in FOUNDER_TIER_GENES:
        return 0
    if gene in CONTEXT_FOUNDER_GENES:
        return 1
    if gene in FALLBACK_PRIORITY_GENES:
        return 2
    return 3


def founder_rank(gene: str, flags: dict[str, bool], candidate_genes: set[str]) -> tuple:
    category = founder_category(gene)
    subtype_score = 0

    if flags["lobular_like"]:
        if gene == "CDH1":
            subtype_score -= 30
        elif gene == "PIK3CA" and "CDH1" in candidate_genes:
            subtype_score -= 20

    if flags["aggressive"] and gene in {"TP53", "PTEN"}:
        subtype_score -= 25

    if flags["luminal_like"] and gene in {"PIK3CA", "GATA3", "MAP3K1"}:
        subtype_score -= 15

    if category >= 2 and gene in TREATMENT_RELATED_GENES:
        subtype_score += 25
    if category >= 2 and gene in RECEPTOR_MARKER_GENES:
        subtype_score += 25
    if category >= 2 and gene in AMPLIFICATION_LIKE_GENES:
        subtype_score += 15

    fallback_index = (
        FALLBACK_PRIORITY_GENES.index(gene) if gene in FALLBACK_PRIORITY_GENES else len(FALLBACK_PRIORITY_GENES)
    )
    tier_index = (
        FOUNDER_TIER_GENES.index(gene)
        if gene in FOUNDER_TIER_GENES
        else CONTEXT_FOUNDER_GENES.index(gene)
        if gene in CONTEXT_FOUNDER_GENES
        else fallback_index
    )
    return (category, subtype_score, tier_index, gene)


def select_origin_genes(
    patient_id: str,
    patient_data: dict,
    candidate_genes: set[str],
    gene_signatures: dict[str, tuple[str, ...]],
    fallback_genes: set[str] | None = None,
) -> tuple[list[str], bool]:
    del patient_id  # Origin selection is deterministic from gene content and subtype context.

    usable_genes = set(candidate_genes)
    if not usable_genes and fallback_genes:
        usable_genes = set(fallback_genes)
    if not usable_genes:
        return [], True

    non_treatment_genes = {gene for gene in usable_genes if gene not in TREATMENT_RELATED_GENES}
    if non_treatment_genes:
        usable_genes = non_treatment_genes

    flags = subtype_flags(patient_data)
    ranked = sorted(usable_genes, key=lambda gene: founder_rank(gene, flags, usable_genes))
    best_category = founder_category(ranked[0])
    origin_signature = gene_signatures[ranked[0]]

    if best_category > 1:
        return [ranked[0]], False

    selected = [ranked[0]]
    for gene in ranked[1:]:
        if len(selected) >= 2:
            break
        if founder_category(gene) > 1:
            break
        if gene_signatures[gene] != origin_signature:
            continue
        selected.append(gene)

    return sorted(selected), False


def deterministic_shuffle(patient_id: str, genes: list[str], salt: str = "") -> list[str]:
    ordered = sorted(genes)
    rng = random.Random(f"{patient_id}:{salt}")
    rng.shuffle(ordered)
    return ordered


def allocate_noninformative_buckets(
    patient_id: str, genes: list[str], max_clones: int
) -> tuple[list[list[str]], list[str]]:
    if not genes:
        return [], []
    if max_clones <= 0:
        return [], sorted(genes)

    ordered = deterministic_shuffle(patient_id, genes, "origin-branch")
    bucket_count = min(max_clones, max(1, (len(ordered) + 1) // 2))
    buckets = [[] for _ in range(bucket_count)]

    index = 0
    remaining = len(ordered)
    for bucket_index in range(bucket_count):
        if remaining <= 0:
            break

        if bucket_index == 0:
            target = 2 if remaining >= 2 else 1
        else:
            minimum_needed_for_rest = bucket_count - bucket_index - 1
            target = 2 if remaining - 2 >= minimum_needed_for_rest else 1

        take = min(target, remaining)
        buckets[bucket_index].extend(ordered[index : index + take])
        index += take
        remaining -= take

    overflow_genes = ordered[index:]
    if overflow_genes:
        if bucket_count > 1:
            overflow_order = list(range(1, bucket_count))
            rng = random.Random(f"{patient_id}:overflow")
            rng.shuffle(overflow_order)
            for offset, gene in enumerate(overflow_genes):
                buckets[overflow_order[offset % len(overflow_order)]].append(gene)
            overflow_genes = []
        else:
            origin_capacity = max(0, 3 - len(buckets[0]))
            if origin_capacity:
                buckets[0].extend(overflow_genes[:origin_capacity])
                overflow_genes = overflow_genes[origin_capacity:]

    return [sorted(bucket) for bucket in buckets if bucket], overflow_genes


def load_patient_data(observations_path: Path) -> dict[str, dict]:
    patients: dict[str, dict] = defaultdict(
        lambda: {
            "genomic_dates": set(),
            "her2_dates": set(),
            "baseline_driver_genes": set(),
            "special_baseline_genes": set(),
            "trend_dates": set(),
            "gene_trends": defaultdict(dict),
            "her2_status": None,
            "er_status": None,
            "pr_status": None,
            "stage_values": set(),
            "responses": set(),
        }
    )

    with observations_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            patient_id = row["PATIENT"].strip()
            description = row["DESCRIPTION"].strip()
            value = row["VALUE"].strip()
            date = row["DATE"].strip()
            patient = patients[patient_id]

            if description == GENETIC_ASSESSMENT_DESCRIPTION:
                patient["genomic_dates"].add(date)
                gene = extract_baseline_gene(value)
                if gene:
                    patient["baseline_driver_genes"].add(gene)
                continue

            if description == HER2_IHC_DESCRIPTION:
                patient["her2_status"] = value
                if value == HER2_POSITIVE_VALUE:
                    patient["her2_dates"].add(date)
                    patient["special_baseline_genes"].add("HER2")
                continue

            if description == ER_IHC_DESCRIPTION:
                patient["er_status"] = value
                continue

            if description == PR_IHC_DESCRIPTION:
                patient["pr_status"] = value
                continue

            if description == STAGE_DESCRIPTION:
                patient["stage_values"].add(value)
                continue

            if description == RESPONSE_DESCRIPTION:
                patient["responses"].add(value)
                continue

            if description.startswith(POST_TREATMENT_PREFIX):
                gene = description[len(POST_TREATMENT_PREFIX) :].strip()
                if gene:
                    patient["trend_dates"].add(date)
                    patient["gene_trends"][gene][date] = value

    return patients


def baseline_date(patient_data: dict) -> str | None:
    if patient_data["genomic_dates"]:
        return min(patient_data["genomic_dates"])
    if patient_data["her2_dates"]:
        return min(patient_data["her2_dates"])
    if patient_data["trend_dates"]:
        return min(patient_data["trend_dates"])
    return None


def build_clone_rows(patients: dict[str, dict]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for patient_id in patients:
        patient = patients[patient_id]
        trend_dates = sorted(patient["trend_dates"])
        t0 = baseline_date(patient)
        genes = sorted(
            set(patient["baseline_driver_genes"])
            | set(patient["special_baseline_genes"])
            | set(patient["gene_trends"].keys())
        )

        if not genes:
            continue

        if not trend_dates:
            continue

        gene_signatures = {
            gene: tuple(
                ["present", *[patient["gene_trends"].get(gene, {}).get(date, "missing") for date in trend_dates]]
            )
            for gene in genes
        }

        if len(genes) <= 2 and len({gene_signatures[gene] for gene in genes}) == 1:
            rows.append(
                {
                    "patient_id": patient_id,
                    "clone_id": "clone_1",
                    "clone_type": "founding",
                    "parent_clone_id": "",
                    "timepoint_dates": ";".join([d for d in [t0, *trend_dates] if d]),
                    "signature_key": signature_key(gene_signatures[genes[0]]),
                    "gene_count": str(len(genes)),
                    "genes": ";".join(genes),
                    "cumulative_gene_count": str(len(genes)),
                }
            )
            continue

        random_driver_genes = infer_random_driver_genes(patient, trend_dates)
        origin_genes, implicit_origin = select_origin_genes(
            patient_id, patient, random_driver_genes, gene_signatures, set(genes)
        )
        remaining_random_genes = set(random_driver_genes) - set(origin_genes)
        noninformative_random_genes = sorted(
            gene
            for gene in remaining_random_genes
            if all(state in NONINFORMATIVE_STATES for state in gene_signatures[gene][1:])
        )
        branch_signature_groups: dict[tuple[str, ...], list[str]] = defaultdict(list)
        for gene in noninformative_random_genes:
            branch_signature_groups[gene_signatures[gene]].append(gene)

        if origin_genes and not implicit_origin:
            origin_signature = gene_signatures[origin_genes[0]]
            origin_capacity = max(0, 3 - len(origin_genes))
            matching_branch_genes = branch_signature_groups.get(origin_signature, [])
            if origin_capacity and matching_branch_genes:
                absorbed = deterministic_shuffle(
                    patient_id, matching_branch_genes, "origin-absorb"
                )[:origin_capacity]
                origin_genes = sorted(set(origin_genes) | set(absorbed))
                remaining_matching = sorted(
                    gene for gene in matching_branch_genes if gene not in set(absorbed)
                )
                if remaining_matching:
                    branch_signature_groups[origin_signature] = remaining_matching
                else:
                    branch_signature_groups.pop(origin_signature, None)

        late_signature_groups: dict[tuple[str, ...], list[str]] = defaultdict(list)

        for gene in genes:
            if gene in noninformative_random_genes or gene in origin_genes:
                continue
            late_signature_groups[gene_signatures[gene]].append(gene)

        for signature in list(branch_signature_groups):
            if signature in late_signature_groups:
                late_signature_groups[signature].extend(branch_signature_groups.pop(signature))

        clone_rows: list[dict[str, str]] = []
        patient_timepoints = ";".join([d for d in [t0, *trend_dates] if d])
        sorted_late_signatures = sorted(late_signature_groups, key=signature_key)
        branch_groups = [
            (signature, sorted(branch_signature_groups[signature]))
            for signature in sorted(branch_signature_groups, key=signature_key)
        ]

        clone_counter = 1
        cumulative_ancestry: list[str] = []

        clone_rows.append(
            {
                "patient_id": patient_id,
                "clone_id": "clone_1",
                "clone_type": "founding",
                "parent_clone_id": "",
                "timepoint_dates": patient_timepoints,
                "signature_key": "implicit:founder"
                if implicit_origin
                else signature_key(gene_signatures[origin_genes[0]]),
                "gene_count": str(len(origin_genes)),
                "genes": ";".join(origin_genes),
                "cumulative_gene_count": str(len(origin_genes)),
            }
        )
        cumulative_ancestry = sorted(origin_genes)
        last_early_clone_id = "clone_1"
        clone_counter += 1

        for signature, bucket_genes in branch_groups:
            clone_id = f"clone_{clone_counter}"
            cumulative_ancestry = sorted(set(cumulative_ancestry) | set(bucket_genes))
            clone_rows.append(
                {
                    "patient_id": patient_id,
                    "clone_id": clone_id,
                    "clone_type": "branch",
                    "parent_clone_id": last_early_clone_id,
                    "timepoint_dates": patient_timepoints,
                    "signature_key": signature_key(signature),
                    "gene_count": str(len(bucket_genes)),
                    "genes": ";".join(bucket_genes),
                    "cumulative_gene_count": str(len(cumulative_ancestry)),
                }
            )
            last_early_clone_id = clone_id
            clone_counter += 1

        base_ancestry = list(cumulative_ancestry)
        for signature in sorted_late_signatures:
            grouped_genes = sorted(late_signature_groups[signature])
            cumulative_genes = sorted(set(base_ancestry) | set(grouped_genes))
            clone_rows.append(
                {
                    "patient_id": patient_id,
                    "clone_id": f"clone_{clone_counter}",
                    "clone_type": "late",
                    "parent_clone_id": last_early_clone_id,
                    "timepoint_dates": patient_timepoints,
                    "signature_key": signature_key(signature),
                    "gene_count": str(len(grouped_genes)),
                    "genes": ";".join(grouped_genes),
                    "cumulative_gene_count": str(len(cumulative_genes)),
                }
            )
            clone_counter += 1

        rows.extend(clone_rows)

    return rows


def write_rows(output_path: Path, rows: list[dict[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "patient_id",
                "clone_id",
                "clone_type",
                "parent_clone_id",
                "timepoint_dates",
                "signature_key",
                "gene_count",
                "genes",
                "cumulative_gene_count",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]

    observations_path = args.observations or newest_observations_csv(repo_root)
    if not observations_path.exists():
        raise FileNotFoundError(f"Observations file not found: {observations_path}")

    output_path = args.output or observations_path.with_name(DEFAULT_OUTPUT_NAME)

    patients = load_patient_data(observations_path)
    rows = build_clone_rows(patients)
    write_rows(output_path, rows)

    patient_count = len({row["patient_id"] for row in rows})
    print(f"Input observations: {observations_path}")
    print(f"Wrote clone groups: {output_path}")
    print(f"Patients with clone output: {patient_count}")
    print(f"Total clone rows: {len(rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
