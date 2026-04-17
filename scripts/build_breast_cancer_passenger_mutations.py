#!/usr/bin/env python3
"""Build breast-cancer passenger mutation files from the ICGC MAF.

Passenger mutations are defined here as breast-cancer MAF rows whose Hugo_Symbol
does not appear in the driver TSV gene column.

The script writes two derived files without modifying the original MAF:
1. a breast passenger-only MAF pool with the exact original MAF columns
2. an assigned passenger mutation file where each synthetic patient is mapped to
   one source tumor sample and receives all passenger rows from that sample
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import random
from collections import Counter, defaultdict
from pathlib import Path

from build_breast_cancer_clones import newest_observations_csv, patient_order_from_observations


DEFAULT_DRIVER_TSV_NAME = "TableS3_panorama_driver_mutations_ICGC_samples.controlled.tsv"
DEFAULT_MAF_NAME = "final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf"
DEFAULT_PASSENGER_POOL_NAME = "breast_cancer_passenger_only_from_maf.maf"
DEFAULT_ASSIGNED_NAME = "breast_cancer_assigned_passenger_mutations.tsv"
DEFAULT_PROJECT_PREFIX = "Breast-"

ASSIGNED_PATIENT_ID = "assigned_patient_id"
SOURCE_TUMOR_SAMPLE_BARCODE = "source_tumor_sample_barcode"
SOURCE_DONOR_ID = "source_donor_id"
SOURCE_PROJECT_CODE = "source_project_code"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a breast passenger-only MAF by removing driver genes from the "
            "ICGC MAF, then assign one source sample's passenger mutation set to "
            "each synthetic patient."
        )
    )
    parser.add_argument(
        "--maf",
        type=Path,
        help=(
            "Path to the tab-delimited ICGC MAF. Defaults to "
            "scripts/final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf."
        ),
    )
    parser.add_argument(
        "--driver-tsv",
        type=Path,
        help=(
            "Path to the driver TSV whose gene column defines driver genes. "
            "Defaults to scripts/TableS3_panorama_driver_mutations_ICGC_samples.controlled.tsv."
        ),
    )
    parser.add_argument(
        "--patients",
        type=Path,
        help=(
            "Path to patients.csv for the synthetic cohort. Defaults to the newest "
            "output_runs/*/csv/patients.csv."
        ),
    )
    parser.add_argument(
        "--observations",
        type=Path,
        help=(
            "Path to observations.csv used to define patient ordering and to limit "
            "assignment to patients with genomic sequencing. Defaults to the newest "
            "output_runs/*/csv/observations.csv."
        ),
    )
    parser.add_argument(
        "--passenger-maf",
        type=Path,
        help=(
            "Output path for the passenger-only MAF. Defaults to "
            "scripts/breast_cancer_passenger_only_from_maf.maf."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Output path for assigned passenger mutations. Defaults to "
            "breast_cancer_assigned_passenger_mutations.tsv next to patients.csv."
        ),
    )
    parser.add_argument(
        "--project-prefix",
        default=DEFAULT_PROJECT_PREFIX,
        help=(
            "Only retain MAF rows whose Project_Code starts with this prefix. "
            f"Default: {DEFAULT_PROJECT_PREFIX}"
        ),
    )
    return parser.parse_args()


def newest_patients_csv(repo_root: Path) -> Path:
    candidates = sorted(
        repo_root.glob("output_runs/*/csv/patients.csv"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError("No patients.csv found under output_runs/*/csv/")
    return candidates[0]


def load_driver_genes(driver_tsv_path: Path) -> set[str]:
    driver_genes: set[str] = set()
    with driver_tsv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row.get("gene", "").strip()
            if gene and gene.lower() != "x":
                driver_genes.add(gene)
    if not driver_genes:
        raise ValueError(f"No driver genes were loaded from {driver_tsv_path}")
    return driver_genes


def build_passenger_maf(
    maf_path: Path,
    driver_genes: set[str],
    output_path: Path,
    project_prefix: str,
) -> tuple[list[str], Counter, dict[str, str], dict[str, str]]:
    sample_counts: Counter = Counter()
    sample_to_donor: dict[str, str] = {}
    sample_to_project: dict[str, str] = {}
    stats: Counter = Counter()

    with maf_path.open(newline="", encoding="utf-8") as source, output_path.open(
        "w", newline="", encoding="utf-8"
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in {maf_path}")
        fieldnames = list(reader.fieldnames)
        writer = csv.DictWriter(target, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for row in reader:
            stats["maf_rows_total"] += 1
            project_code = row.get("Project_Code", "").strip()
            if not project_code.startswith(project_prefix):
                continue

            stats["project_rows_seen"] += 1
            gene = row.get("Hugo_Symbol", "").strip()
            if gene in driver_genes:
                stats["driver_rows_removed"] += 1
                continue

            writer.writerow(row)
            stats["passenger_rows_written"] += 1

            sample = row.get("Tumor_Sample_Barcode", "").strip()
            if sample:
                sample_counts[sample] += 1
                sample_to_donor[sample] = row.get("Donor_ID", "").strip()
                sample_to_project[sample] = project_code

    if not sample_counts:
        raise ValueError(
            "No passenger rows remained after driver-gene filtering. "
            "Check the MAF path, project prefix, and driver TSV."
        )

    return fieldnames, stats, sample_to_donor, sample_to_project, sample_counts


def load_patient_ids(patients_path: Path) -> list[str]:
    patient_ids: list[str] = []
    with patients_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            patient_id = row.get("Id", "").strip()
            if patient_id:
                patient_ids.append(patient_id)
    if not patient_ids:
        raise ValueError(f"No patient ids were loaded from {patients_path}")
    return patient_ids


def load_ordered_genomic_patient_ids(
    patients_path: Path, observations_path: Path
) -> list[str]:
    patient_ids = set(load_patient_ids(patients_path))
    ordered_patient_ids = [
        patient_id
        for patient_id in patient_order_from_observations(observations_path)
        if patient_id in patient_ids
    ]
    if not ordered_patient_ids:
        raise ValueError(
            "No patients with genomic sequencing were found in observations.csv "
            f"for cohort {patients_path.parent}"
        )
    return ordered_patient_ids


def choose_source_sample(patient_id: str, source_samples: list[str]) -> str:
    digest = hashlib.sha256(patient_id.encode("utf-8")).hexdigest()
    seed = int(digest[:16], 16)
    rng = random.Random(seed)
    return source_samples[rng.randrange(len(source_samples))]


def assign_source_samples(
    patient_ids: list[str], available_samples: Counter
) -> dict[str, str]:
    source_samples = sorted(available_samples)
    if not source_samples:
        raise ValueError("No source samples are available for passenger assignment.")

    assignments: dict[str, str] = {}
    for patient_id in patient_ids:
        assignments[patient_id] = choose_source_sample(patient_id, source_samples)
    return assignments


def write_assigned_passenger_rows(
    passenger_maf_path: Path,
    maf_fieldnames: list[str],
    patient_ids_in_order: list[str],
    assignments: dict[str, str],
    sample_to_donor: dict[str, str],
    sample_to_project: dict[str, str],
    output_path: Path,
) -> Counter:
    reverse_assignments: dict[str, list[str]] = defaultdict(list)
    for patient_id, sample in assignments.items():
        reverse_assignments[sample].append(patient_id)
    patient_order = {patient_id: index for index, patient_id in enumerate(patient_ids_in_order)}

    stats: Counter = Counter()
    helper_fieldnames = [
        ASSIGNED_PATIENT_ID,
        SOURCE_TUMOR_SAMPLE_BARCODE,
        SOURCE_DONOR_ID,
        SOURCE_PROJECT_CODE,
    ]
    ordered_rows: list[tuple[int, int, dict[str, str]]] = []
    sequence_number = 0

    with passenger_maf_path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source, delimiter="\t")

        for row in reader:
            sample = row.get("Tumor_Sample_Barcode", "").strip()
            assigned_patients = reverse_assignments.get(sample)
            if not assigned_patients:
                continue

            for patient_id in assigned_patients:
                output_row = {
                    ASSIGNED_PATIENT_ID: patient_id,
                    SOURCE_TUMOR_SAMPLE_BARCODE: sample,
                    SOURCE_DONOR_ID: sample_to_donor.get(sample, row.get("Donor_ID", "").strip()),
                    SOURCE_PROJECT_CODE: sample_to_project.get(sample, row.get("Project_Code", "").strip()),
                }
                output_row.update(row)
                ordered_rows.append((patient_order[patient_id], sequence_number, output_row))
                sequence_number += 1
                stats["assigned_rows_written"] += 1

    ordered_rows.sort(key=lambda item: (item[0], item[1]))

    with output_path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(
            target,
            fieldnames=[*helper_fieldnames, *maf_fieldnames],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for _, _, row in ordered_rows:
            writer.writerow(row)

    stats["patients_assigned"] = len(assignments)
    stats["source_samples_used"] = len(set(assignments.values()))
    return stats


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]

    maf_path = args.maf or repo_root / "scripts" / DEFAULT_MAF_NAME
    driver_tsv_path = args.driver_tsv or repo_root / "scripts" / DEFAULT_DRIVER_TSV_NAME
    patients_path = args.patients or newest_patients_csv(repo_root)
    observations_path = args.observations or newest_observations_csv(repo_root)
    passenger_maf_path = args.passenger_maf or repo_root / "scripts" / DEFAULT_PASSENGER_POOL_NAME
    assigned_output_path = args.output or patients_path.with_name(DEFAULT_ASSIGNED_NAME)

    driver_genes = load_driver_genes(driver_tsv_path)
    (
        maf_fieldnames,
        passenger_stats,
        sample_to_donor,
        sample_to_project,
        sample_counts,
    ) = build_passenger_maf(
        maf_path=maf_path,
        driver_genes=driver_genes,
        output_path=passenger_maf_path,
        project_prefix=args.project_prefix,
    )

    patient_ids = load_ordered_genomic_patient_ids(patients_path, observations_path)
    assignments = assign_source_samples(patient_ids, sample_counts)
    assigned_stats = write_assigned_passenger_rows(
        passenger_maf_path=passenger_maf_path,
        maf_fieldnames=maf_fieldnames,
        patient_ids_in_order=patient_ids,
        assignments=assignments,
        sample_to_donor=sample_to_donor,
        sample_to_project=sample_to_project,
        output_path=assigned_output_path,
    )

    print(f"Driver genes loaded: {len(driver_genes)}")
    print(f"Passenger-only MAF written: {passenger_maf_path}")
    print(f"Assigned passenger mutations written: {assigned_output_path}")
    print(f"Observations order source: {observations_path}")
    for key in (
        "maf_rows_total",
        "project_rows_seen",
        "driver_rows_removed",
        "passenger_rows_written",
    ):
        print(f"{key}: {passenger_stats[key]}")
    print(f"available_source_samples: {len(sample_counts)}")
    print(f"patients_assigned: {assigned_stats['patients_assigned']}")
    print(f"source_samples_used: {assigned_stats['source_samples_used']}")
    print(f"assigned_rows_written: {assigned_stats['assigned_rows_written']}")


if __name__ == "__main__":
    main()
