#!/usr/bin/env python3
"""Proof-of-concept VAF-informed passenger-to-clone assignment.

This script does not modify the existing maf_files/ directory. It reads the
assigned passenger mutation TSV and breast_cancer_clone_proportions.csv, then
assigns each passenger mutation to a clone using either nearest-VAF matching or
a ranked split where low-VAF passengers are reserved for low-abundance clones.

The output is intended for inspection, not as a replacement for the production
MAF-generation pipeline.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


ASSIGNED_PATIENT_ID = "assigned_patient_id"
SOURCE_TUMOR_SAMPLE_BARCODE = "source_tumor_sample_barcode"
SOURCE_DONOR_ID = "source_donor_id"
MUTATION_ORIGIN = "mutation_origin"
PASSENGER_VAF_FIELD = "i_VAF"
T_ALT_COUNT_FIELD = "t_alt_count"
T_REF_COUNT_FIELD = "t_ref_count"

HELPER_FIELDS = [
    MUTATION_ORIGIN,
    SOURCE_TUMOR_SAMPLE_BARCODE,
    SOURCE_DONOR_ID,
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Proof-of-concept split of assigned breast cancer passenger mutations "
            "into clone-specific passenger MAFs using clone VAFs."
        )
    )
    parser.add_argument(
        "--assigned-passengers",
        type=Path,
        required=True,
        help="Path to breast_cancer_assigned_passenger_mutations.tsv.",
    )
    parser.add_argument(
        "--clone-proportions",
        type=Path,
        required=True,
        help="Path to breast_cancer_clone_proportions.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Output directory for proof-of-concept files. Defaults to "
            "passenger_vaf_split_poc/ under the run root."
        ),
    )
    parser.add_argument(
        "--timepoint",
        choices=["t0", "t1", "t2", "max", "mean"],
        default="t0",
        help=(
            "Clone VAF value to match against passenger VAF. Use t0/t1/t2 for a "
            "specific sequencing timepoint, max for each clone's maximum observed "
            "VAF, or mean for its mean observed VAF. Default: t0."
        ),
    )
    parser.add_argument(
        "--method",
        choices=["nearest", "ranked"],
        default="nearest",
        help=(
            "Assignment method. nearest assigns each passenger to the clone with "
            "the closest selected clone VAF. ranked sorts passengers by VAF and "
            "allocates the lowest VAF passengers to the lowest-abundance clones, "
            "with clone bin sizes proportional to selected clone VAF. Default: nearest."
        ),
    )
    parser.add_argument(
        "--write-mafs",
        action="store_true",
        help="Write passenger-only clone MAF files in addition to assignment TSVs.",
    )
    return parser.parse_args()


def default_output_dir(assigned_passengers: Path) -> Path:
    csv_dir = assigned_passengers.parent
    if csv_dir.name == "csv":
        return csv_dir.parent / "passenger_vaf_split_poc"
    return csv_dir / "passenger_vaf_split_poc"


def parse_float(value: str) -> float | None:
    cleaned = (value or "").strip()
    if not cleaned or cleaned.upper() in {"NA", "N/A", "NONE", "NULL"}:
        return None
    try:
        return float(cleaned)
    except ValueError:
        return None


def derive_maf_fieldnames(fieldnames: list[str]) -> list[str]:
    if "Hugo_Symbol" not in fieldnames or "Donor_ID" not in fieldnames:
        raise ValueError("Expected Hugo_Symbol and Donor_ID in assigned passenger TSV header")
    start = fieldnames.index("Hugo_Symbol")
    end = fieldnames.index("Donor_ID")
    return fieldnames[start : end + 1]


def passenger_vaf_pct(row: dict[str, str]) -> tuple[float | None, str]:
    direct_vaf = parse_float(row.get(PASSENGER_VAF_FIELD, ""))
    if direct_vaf is not None:
        if 0 <= direct_vaf <= 1:
            return direct_vaf * 100.0, PASSENGER_VAF_FIELD
        return direct_vaf, PASSENGER_VAF_FIELD

    alt_count = parse_float(row.get(T_ALT_COUNT_FIELD, ""))
    ref_count = parse_float(row.get(T_REF_COUNT_FIELD, ""))
    if alt_count is not None and ref_count is not None and (alt_count + ref_count) > 0:
        return (alt_count / (alt_count + ref_count)) * 100.0, "alt_ref_counts"

    return None, "missing"


def clone_vaf_for_mode(row: dict[str, str], timepoint: str) -> tuple[float | None, str]:
    if timepoint in {"t0", "t1", "t2"}:
        return parse_float(row.get(f"{timepoint}_vaf_pct", "")), row.get(f"{timepoint}_date", "")

    values = [
        parse_float(row.get(f"{label}_vaf_pct", ""))
        for label in ("t0", "t1", "t2")
    ]
    present = [value for value in values if value is not None]
    if not present:
        return None, ""
    if timepoint == "max":
        return max(present), "max"
    if timepoint == "mean":
        return statistics.mean(present), "mean"
    raise ValueError(f"Unsupported timepoint mode: {timepoint}")


def load_clone_targets(path: Path, timepoint: str) -> dict[str, list[dict[str, object]]]:
    clones_by_patient: dict[str, list[dict[str, object]]] = defaultdict(list)

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            patient_id = row["patient_id"].strip()
            clone_vaf_pct, clone_date = clone_vaf_for_mode(row, timepoint)
            if clone_vaf_pct is None:
                continue
            clones_by_patient[patient_id].append(
                {
                    "patient_id": patient_id,
                    "clone_id": row["clone_id"].strip(),
                    "clone_type": row["clone_type"].strip(),
                    "clone_vaf_pct": clone_vaf_pct,
                    "clone_date": clone_date,
                }
            )

    for patient_id in list(clones_by_patient):
        clones_by_patient[patient_id].sort(
            key=lambda clone: (str(clone["clone_id"]), str(clone["clone_type"]))
        )

    return clones_by_patient


def choose_clone(
    clones: list[dict[str, object]], passenger_vaf: float | None
) -> tuple[dict[str, object] | None, float | None, str]:
    if not clones:
        return None, None, "missing_clone_proportions"

    if passenger_vaf is None:
        selected = max(
            clones,
            key=lambda clone: (float(clone["clone_vaf_pct"]), str(clone["clone_id"])),
        )
        return selected, None, "missing_passenger_vaf_largest_clone"

    selected = min(
        clones,
        key=lambda clone: (
            abs(float(clone["clone_vaf_pct"]) - passenger_vaf),
            str(clone["clone_id"]),
        ),
    )
    distance = abs(float(selected["clone_vaf_pct"]) - passenger_vaf)
    return selected, distance, "nearest_vaf"


def clone_vaf_sort_key(clone: dict[str, object]) -> tuple[float, str, str]:
    return (
        float(clone["clone_vaf_pct"]),
        str(clone["clone_id"]),
        str(clone["clone_type"]),
    )


def passenger_sort_key(item: tuple[dict[str, str], float, str, int]) -> tuple[float, str, str, str, int]:
    row, passenger_vaf, _passenger_vaf_source, index = item
    return (
        passenger_vaf,
        row.get("Chromosome", ""),
        row.get("Start_position", ""),
        row.get("Genome_Change", ""),
        index,
    )


def proportional_clone_counts(
    total_rows: int, clones: list[dict[str, object]]
) -> list[int]:
    if total_rows == 0:
        return [0 for _clone in clones]

    weights = [max(float(clone["clone_vaf_pct"]), 0.0) for clone in clones]
    weight_total = sum(weights)
    if weight_total <= 0:
        base = total_rows // len(clones)
        counts = [base for _clone in clones]
        for index in range(total_rows - sum(counts)):
            counts[index] += 1
        return counts

    raw_counts = [(total_rows * weight) / weight_total for weight in weights]
    counts = [math.floor(raw_count) for raw_count in raw_counts]
    remaining = total_rows - sum(counts)
    remainder_order = sorted(
        range(len(clones)),
        key=lambda index: (
            raw_counts[index] - counts[index],
            float(clones[index]["clone_vaf_pct"]),
            str(clones[index]["clone_id"]),
        ),
        reverse=True,
    )
    for index in remainder_order[:remaining]:
        counts[index] += 1
    return counts


def build_nearest_assignments(
    passenger_rows: list[dict[str, str]],
    clones_by_patient: dict[str, list[dict[str, object]]],
) -> list[dict[str, object]]:
    assignments: list[dict[str, object]] = []
    for row in passenger_rows:
        patient_id = row.get(ASSIGNED_PATIENT_ID, "").strip()
        passenger_vaf, passenger_vaf_source = passenger_vaf_pct(row)
        clone, distance, reason = choose_clone(clones_by_patient.get(patient_id, []), passenger_vaf)
        assignments.append(
            {
                "row": row,
                "clone": clone,
                "distance": distance,
                "reason": reason,
                "passenger_vaf": passenger_vaf,
                "passenger_vaf_source": passenger_vaf_source,
            }
        )
    return assignments


def build_ranked_assignments(
    passenger_rows: list[dict[str, str]],
    clones_by_patient: dict[str, list[dict[str, object]]],
) -> list[dict[str, object]]:
    rows_by_patient: dict[str, list[tuple[int, dict[str, str]]]] = defaultdict(list)
    for index, row in enumerate(passenger_rows):
        rows_by_patient[row.get(ASSIGNED_PATIENT_ID, "").strip()].append((index, row))

    assignments: list[dict[str, object]] = []
    for patient_id, indexed_rows in rows_by_patient.items():
        clones = clones_by_patient.get(patient_id, [])
        if not clones:
            for _index, row in indexed_rows:
                passenger_vaf, passenger_vaf_source = passenger_vaf_pct(row)
                assignments.append(
                    {
                        "row": row,
                        "clone": None,
                        "distance": None,
                        "reason": "missing_clone_proportions",
                        "passenger_vaf": passenger_vaf,
                        "passenger_vaf_source": passenger_vaf_source,
                    }
                )
            continue

        clones_sorted = sorted(clones, key=clone_vaf_sort_key)
        largest_clone = max(
            clones,
            key=lambda clone: (float(clone["clone_vaf_pct"]), str(clone["clone_id"])),
        )
        numeric_rows: list[tuple[dict[str, str], float, str, int]] = []

        for index, row in indexed_rows:
            passenger_vaf, passenger_vaf_source = passenger_vaf_pct(row)
            if passenger_vaf is None:
                assignments.append(
                    {
                        "row": row,
                        "clone": largest_clone,
                        "distance": None,
                        "reason": "missing_passenger_vaf_largest_clone",
                        "passenger_vaf": passenger_vaf,
                        "passenger_vaf_source": passenger_vaf_source,
                    }
                )
                continue
            numeric_rows.append((row, passenger_vaf, passenger_vaf_source, index))

        numeric_rows.sort(key=passenger_sort_key)
        counts = proportional_clone_counts(len(numeric_rows), clones_sorted)
        start = 0
        for clone, count in zip(clones_sorted, counts):
            for row, passenger_vaf, passenger_vaf_source, _index in numeric_rows[start : start + count]:
                distance = abs(float(clone["clone_vaf_pct"]) - passenger_vaf)
                assignments.append(
                    {
                        "row": row,
                        "clone": clone,
                        "distance": distance,
                        "reason": "ranked_vaf_bin",
                        "passenger_vaf": passenger_vaf,
                        "passenger_vaf_source": passenger_vaf_source,
                    }
                )
            start += count

    return assignments


def safe_maf_filename(clone_id: str, clone_type: str) -> str:
    safe_type = clone_type or "clone"
    return f"{clone_id}_{safe_type}_passengers.maf"


def passenger_output_row(row: dict[str, str], maf_fieldnames: list[str]) -> dict[str, str]:
    output = {
        MUTATION_ORIGIN: "passenger",
        SOURCE_TUMOR_SAMPLE_BARCODE: row.get(SOURCE_TUMOR_SAMPLE_BARCODE, ""),
        SOURCE_DONOR_ID: row.get(SOURCE_DONOR_ID, row.get("Donor_ID", "")),
    }
    for field in maf_fieldnames:
        output[field] = row.get(field, "")
    return output


def write_passenger_mafs(
    output_dir: Path,
    rows_by_patient_clone: dict[tuple[str, str, str], list[dict[str, str]]],
    maf_fieldnames: list[str],
) -> int:
    fieldnames = [*HELPER_FIELDS, *maf_fieldnames]
    files_written = 0

    for (patient_id, clone_id, clone_type), rows in sorted(rows_by_patient_clone.items()):
        patient_dir = output_dir / "passenger_maf_files" / patient_id
        patient_dir.mkdir(parents=True, exist_ok=True)
        output_path = patient_dir / safe_maf_filename(clone_id, clone_type)
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow(passenger_output_row(row, maf_fieldnames))
        files_written += 1

    return files_written


def write_summary(
    output_path: Path,
    clone_rows: list[dict[str, str]],
    rows_by_patient_clone: dict[tuple[str, str, str], list[dict[str, str]]],
) -> None:
    fieldnames = [
        "patient_id",
        "clone_id",
        "clone_type",
        "passenger_count",
        "mean_passenger_vaf_pct",
        "min_passenger_vaf_pct",
        "max_passenger_vaf_pct",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in clone_rows:
            key = (row["patient_id"], row["clone_id"], row["clone_type"])
            passenger_rows = rows_by_patient_clone.get(key, [])
            vafs = [passenger_vaf_pct(passenger)[0] for passenger in passenger_rows]
            present_vafs = [vaf for vaf in vafs if vaf is not None]
            writer.writerow(
                {
                    "patient_id": row["patient_id"],
                    "clone_id": row["clone_id"],
                    "clone_type": row["clone_type"],
                    "passenger_count": len(passenger_rows),
                    "mean_passenger_vaf_pct": (
                        f"{statistics.mean(present_vafs):.4f}" if present_vafs else ""
                    ),
                    "min_passenger_vaf_pct": (
                        f"{min(present_vafs):.4f}" if present_vafs else ""
                    ),
                    "max_passenger_vaf_pct": (
                        f"{max(present_vafs):.4f}" if present_vafs else ""
                    ),
                }
            )


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir or default_output_dir(args.assigned_passengers)

    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory already exists and is not empty: {output_dir}. "
            "Choose a new --output-dir so existing proof-of-concept results are not overwritten."
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    clones_by_patient = load_clone_targets(args.clone_proportions, args.timepoint)

    assignment_path = output_dir / "passenger_vaf_clone_assignments.tsv"
    summary_path = output_dir / "passenger_vaf_clone_summary.tsv"
    rows_by_patient_clone: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    summary_clone_rows: dict[tuple[str, str, str], dict[str, str]] = {}

    total_passengers = 0
    assigned_passengers = 0
    missing_clone_patients: set[str] = set()
    maf_fieldnames: list[str] | None = None

    assignment_fields = [
        "assigned_patient_id",
        "clone_id",
        "clone_type",
        "assignment_method",
        "timepoint_mode",
        "clone_date_or_mode",
        "clone_vaf_pct",
        "passenger_vaf_pct",
        "passenger_vaf_source",
        "abs_distance_pct",
        "assignment_reason",
        "source_tumor_sample_barcode",
        "source_donor_id",
        "Hugo_Symbol",
        "Chromosome",
        "Start_position",
        "End_position",
        "Variant_Classification",
        "Variant_Type",
        "Genome_Change",
    ]

    with args.assigned_passengers.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {args.assigned_passengers}")
        maf_fieldnames = derive_maf_fieldnames(list(reader.fieldnames))
        passenger_rows = [dict(row) for row in reader]

    if args.method == "nearest":
        assignments = build_nearest_assignments(passenger_rows, clones_by_patient)
    elif args.method == "ranked":
        assignments = build_ranked_assignments(passenger_rows, clones_by_patient)
    else:
        raise ValueError(f"Unsupported assignment method: {args.method}")

    with assignment_path.open("w", newline="", encoding="utf-8") as assignment_handle:
        assignment_writer = csv.DictWriter(
            assignment_handle,
            fieldnames=assignment_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        assignment_writer.writeheader()

        for assignment in assignments:
            row = assignment["row"]
            total_passengers += 1
            patient_id = row.get(ASSIGNED_PATIENT_ID, "").strip()
            passenger_vaf = assignment["passenger_vaf"]
            passenger_vaf_source = str(assignment["passenger_vaf_source"])
            clone = assignment["clone"]
            distance = assignment["distance"]
            reason = str(assignment["reason"])

            if clone is None:
                missing_clone_patients.add(patient_id)
                assignment_writer.writerow(
                    {
                        "assigned_patient_id": patient_id,
                        "assignment_method": args.method,
                        "timepoint_mode": args.timepoint,
                        "passenger_vaf_pct": f"{passenger_vaf:.4f}" if passenger_vaf is not None else "",
                        "passenger_vaf_source": passenger_vaf_source,
                        "assignment_reason": reason,
                        "source_tumor_sample_barcode": row.get(SOURCE_TUMOR_SAMPLE_BARCODE, ""),
                        "source_donor_id": row.get(SOURCE_DONOR_ID, ""),
                        "Hugo_Symbol": row.get("Hugo_Symbol", ""),
                        "Chromosome": row.get("Chromosome", ""),
                        "Start_position": row.get("Start_position", ""),
                        "End_position": row.get("End_position", ""),
                        "Variant_Classification": row.get("Variant_Classification", ""),
                        "Variant_Type": row.get("Variant_Type", ""),
                        "Genome_Change": row.get("Genome_Change", ""),
                    }
                )
                continue

            assigned_passengers += 1
            clone_id = str(clone["clone_id"])
            clone_type = str(clone["clone_type"])
            clone_key = (patient_id, clone_id, clone_type)
            rows_by_patient_clone[clone_key].append(dict(row))
            summary_clone_rows.setdefault(
                clone_key,
                {
                    "patient_id": patient_id,
                    "clone_id": clone_id,
                    "clone_type": clone_type,
                },
            )

            assignment_writer.writerow(
                {
                    "assigned_patient_id": patient_id,
                    "clone_id": clone_id,
                    "clone_type": clone_type,
                    "assignment_method": args.method,
                    "timepoint_mode": args.timepoint,
                    "clone_date_or_mode": str(clone["clone_date"]),
                    "clone_vaf_pct": f"{float(clone['clone_vaf_pct']):.4f}",
                    "passenger_vaf_pct": f"{passenger_vaf:.4f}" if passenger_vaf is not None else "",
                    "passenger_vaf_source": passenger_vaf_source,
                    "abs_distance_pct": f"{distance:.4f}" if distance is not None else "",
                    "assignment_reason": reason,
                    "source_tumor_sample_barcode": row.get(SOURCE_TUMOR_SAMPLE_BARCODE, ""),
                    "source_donor_id": row.get(SOURCE_DONOR_ID, ""),
                    "Hugo_Symbol": row.get("Hugo_Symbol", ""),
                    "Chromosome": row.get("Chromosome", ""),
                    "Start_position": row.get("Start_position", ""),
                    "End_position": row.get("End_position", ""),
                    "Variant_Classification": row.get("Variant_Classification", ""),
                    "Variant_Type": row.get("Variant_Type", ""),
                    "Genome_Change": row.get("Genome_Change", ""),
                }
            )

    if maf_fieldnames is None:
        raise ValueError(f"No passenger rows found in {args.assigned_passengers}")

    write_summary(summary_path, list(summary_clone_rows.values()), rows_by_patient_clone)
    maf_files_written = 0
    if args.write_mafs:
        maf_files_written = write_passenger_mafs(output_dir, rows_by_patient_clone, maf_fieldnames)

    print(f"Assigned passengers: {args.assigned_passengers}")
    print(f"Clone proportions: {args.clone_proportions}")
    print(f"Timepoint mode: {args.timepoint}")
    print(f"Assignment method: {args.method}")
    print(f"Output directory: {output_dir}")
    print(f"Assignment table: {assignment_path}")
    print(f"Summary table: {summary_path}")
    print(f"Total passenger rows: {total_passengers}")
    print(f"Assigned passenger rows: {assigned_passengers}")
    print(f"Patients missing clone proportions: {len(missing_clone_patients)}")
    if args.write_mafs:
        print(f"Passenger-only MAF files written: {maf_files_written}")
    else:
        print("Passenger-only MAF writing skipped; pass --write-mafs to enable it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
