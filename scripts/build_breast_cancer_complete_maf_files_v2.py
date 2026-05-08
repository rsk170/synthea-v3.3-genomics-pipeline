#!/usr/bin/env python3
"""Create per-patient, per-sequencing-event complete MAF files.

Each output MAF combines:
- driver variants present in observations_pruned_by_clone_vaf.csv for that
  specific sequencing event
- passenger variants assigned to the patient from
  breast_cancer_assigned_passenger_mutations.tsv

Driver rows are recovered from the exact ICGC-backed source libraries using
gene + genomic change + variant classification keys. Passenger rows are copied
from the already-assigned passenger mutation file and included unchanged across
all sequencing events for the same patient.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_PRUNED_NAME = "observations_pruned_by_clone_vaf.csv"
DEFAULT_ASSIGNED_PASSENGERS_NAME = "breast_cancer_assigned_passenger_mutations.tsv"
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
SYNTHETIC_TUMOR_SAMPLE_BARCODE = "synthetic_tumor_sample_barcode"
SEQUENCING_EVENT_DATE = "sequencing_event_date"
MUTATION_ORIGIN = "mutation_origin"
SOURCE_VARIANT_LIBRARY = "source_variant_library"
SOURCE_TUMOR_SAMPLE_BARCODE = "source_tumor_sample_barcode"
SOURCE_DONOR_ID = "source_donor_id"
SOURCE_PROJECT_CODE = "source_project_code"

HELPER_FIELDS = [
    ASSIGNED_PATIENT_ID,
    SYNTHETIC_TUMOR_SAMPLE_BARCODE,
    SEQUENCING_EVENT_DATE,
    MUTATION_ORIGIN,
    SOURCE_VARIANT_LIBRARY,
    SOURCE_TUMOR_SAMPLE_BARCODE,
    SOURCE_DONOR_ID,
    SOURCE_PROJECT_CODE,
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create per-patient complete MAF files for each sequencing event by "
            "combining driver variants from observations_pruned_by_clone_vaf.csv "
            "with assigned passenger mutations."
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


def load_pruned_driver_events(
    path: Path,
) -> dict[str, dict[str, dict[str, object]]]:
    events_by_patient: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            code = row["CODE"].strip()
            description = row["DESCRIPTION"].strip()
            patient_id = row["PATIENT"].strip()
            encounter_id = row["ENCOUNTER"].strip()
            date = row["DATE"].strip()

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
            event = events_by_patient[patient_id].setdefault(
                encounter_id,
                {"date": date, "drivers": {}},
            )
            drivers = event["drivers"]  # type: ignore[assignment]
            if canonical_gene not in drivers:
                drivers[canonical_gene] = {
                    "gene": canonical_gene,
                    "event_gene_name": gene,
                    "genomic_change": "",
                    "variant_classification": "",
                }
            driver = drivers[canonical_gene]
            if code == GENE_CODE:
                driver["event_gene_name"] = gene
            elif code == GENOMIC_CHANGE_CODE:
                driver["genomic_change"] = row["VALUE"].strip()
            elif code == VARIANT_CLASS_CODE:
                driver["variant_classification"] = row["VALUE"].strip()

    return events_by_patient


def collect_missing_keys_from_events(
    events_by_patient: dict[str, dict[str, dict[str, object]]],
    indexes: list[dict[tuple[str, str, str], list[dict[str, str]]]],
) -> set[tuple[str, str, str]]:
    missing: set[tuple[str, str, str]] = set()
    for patient_events in events_by_patient.values():
        for event in patient_events.values():
            drivers = event["drivers"]  # type: ignore[assignment]
            for gene, record in drivers.items():
                genomic_change = str(record["genomic_change"]).strip()
                variant_classification = str(record["variant_classification"]).strip()
                if not genomic_change or not variant_classification:
                    continue
                key = (canonical_gene_name(gene), genomic_change, variant_classification)
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


def synthetic_maf_filename(sequence_index: int, date: str, encounter_id: str) -> str:
    safe_date = date.replace("-", "").replace(":", "")
    return f"sequencing_{sequence_index}_{safe_date}_{encounter_id}.maf"


def build_driver_output_row(
    patient_id: str,
    encounter_id: str,
    event_date: str,
    variant_row: dict[str, str],
    maf_fieldnames: list[str],
) -> dict[str, str]:
    output = {
        ASSIGNED_PATIENT_ID: patient_id,
        SYNTHETIC_TUMOR_SAMPLE_BARCODE: encounter_id,
        SEQUENCING_EVENT_DATE: event_date,
        MUTATION_ORIGIN: "driver",
        SOURCE_VARIANT_LIBRARY: variant_row.get("_source_variant_library", "driver"),
        SOURCE_TUMOR_SAMPLE_BARCODE: variant_row.get("Tumor_Sample_Barcode", ""),
        SOURCE_DONOR_ID: variant_row.get("Donor_ID", ""),
        SOURCE_PROJECT_CODE: variant_row.get("Project_Code", ""),
    }
    for field in maf_fieldnames:
        output[field] = variant_row.get(field, "")
    return output


def build_passenger_output_rows(
    patient_id: str,
    encounter_id: str,
    event_date: str,
    passenger_rows: list[dict[str, str]],
    maf_fieldnames: list[str],
) -> list[dict[str, str]]:
    output_rows: list[dict[str, str]] = []
    for row in passenger_rows:
        output = {
            ASSIGNED_PATIENT_ID: patient_id,
            SYNTHETIC_TUMOR_SAMPLE_BARCODE: encounter_id,
            SEQUENCING_EVENT_DATE: event_date,
            MUTATION_ORIGIN: "passenger",
            SOURCE_VARIANT_LIBRARY: "passenger_assigned",
            SOURCE_TUMOR_SAMPLE_BARCODE: row.get(SOURCE_TUMOR_SAMPLE_BARCODE, row.get("Tumor_Sample_Barcode", "")),
            SOURCE_DONOR_ID: row.get(SOURCE_DONOR_ID, row.get("Donor_ID", "")),
            SOURCE_PROJECT_CODE: row.get(SOURCE_PROJECT_CODE, row.get("Project_Code", "")),
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
    events_by_patient = load_pruned_driver_events(pruned_path)
    needed_maf_fallback_keys = collect_missing_keys_from_events(
        events_by_patient, [driver_index, civic_index, non_disruptive_index]
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

    for patient_id in sorted(events_by_patient):
        passenger_rows = passengers_by_patient.get(patient_id)
        if not passenger_rows:
            stats["patients_skipped_no_passengers"] += 1
            continue

        patient_events = sorted(
            events_by_patient[patient_id].items(),
            key=lambda item: (str(item[1]["date"]), item[0]),
        )
        patient_dir = output_dir / patient_id
        patient_sequence_count = 0

        for sequence_index, (encounter_id, event) in enumerate(patient_events, start=1):
            event_date = str(event["date"])
            drivers = event["drivers"]  # type: ignore[assignment]
            driver_rows: list[dict[str, str]] = []

            for gene in sorted(drivers):
                record = drivers[gene]
                genomic_change = str(record["genomic_change"]).strip()
                variant_classification = str(record["variant_classification"]).strip()
                if not genomic_change or not variant_classification:
                    stats["driver_records_skipped_missing_fields"] += 1
                    continue

                variant_row = choose_variant_row(
                    gene,
                    genomic_change,
                    variant_classification,
                    [driver_index, civic_index, non_disruptive_index, maf_fallback_index],
                )
                if variant_row is None:
                    stats["driver_records_missing_exact_match"] += 1
                    continue

                driver_rows.append(
                    build_driver_output_row(
                        patient_id, encounter_id, event_date, variant_row, maf_fieldnames
                    )
                )
                stats["driver_rows_written"] += 1

            passenger_output_rows = build_passenger_output_rows(
                patient_id, encounter_id, event_date, passenger_rows, maf_fieldnames
            )

            output_path = patient_dir / synthetic_maf_filename(
                sequence_index, event_date, encounter_id
            )
            write_complete_maf(output_path, fieldnames, driver_rows, passenger_output_rows)
            patient_sequence_count += 1
            stats["maf_files_written"] += 1
            stats["passenger_rows_written"] += len(passenger_output_rows)

        if patient_sequence_count:
            stats["patients_written"] += 1

    print(f"Pruned observations: {pruned_path}")
    print(f"Assigned passengers: {assigned_passengers_path}")
    print(f"Raw MAF exact fallback source: {maf_path}")
    print(f"Output directory: {output_dir}")
    print(f"Patients with MAF output: {stats['patients_written']}")
    print(f"Patients skipped (no passenger file rows): {stats['patients_skipped_no_passengers']}")
    print(f"MAF files written: {stats['maf_files_written']}")
    print(f"Driver rows written: {stats['driver_rows_written']}")
    print(f"Passenger rows written: {stats['passenger_rows_written']}")
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
    try:
        main()
    except Exception as exc:  # pragma: no cover - CLI error path
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
