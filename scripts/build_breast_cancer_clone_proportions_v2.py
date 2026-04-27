#!/usr/bin/env python3
"""Assign per-timepoint clonal VAFs to breast cancer clone groups.

The script reads breast_cancer_clone_groups.csv and emits one row per clone with
wide VAF columns for t0, t1, and t2 when present. For each patient and
timepoint, clone VAFs sum to 100.00.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from pathlib import Path

from build_breast_cancer_clones_v2 import (
    DEFAULT_OUTPUT_NAME as CLONE_GROUPS_DEFAULT_OUTPUT_NAME,
    build_clone_rows,
    load_patient_data,
    newest_observations_csv,
    patient_order_from_observations,
    write_rows as write_clone_group_rows,
)


DEFAULT_INPUT_NAME = "breast_cancer_clone_groups.csv"
DEFAULT_OUTPUT_NAME = "breast_cancer_clone_proportions.csv"
LOW_CLONE_MAX = 4.2
LOW_CLONE_MIN = 0.8
ROUNDING_UNITS = 10000  # basis points of 1%


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create per-timepoint breast cancer clonal VAFs from clone groups, "
            "or build clone groups from observations.csv first."
        )
    )
    parser.add_argument(
        "--observations",
        type=Path,
        help=(
            "Path to observations.csv. If provided, clone groups are rebuilt first. "
            "If neither --observations nor --clone-groups is given, the newest "
            "output_runs/*/csv/observations.csv is used."
        ),
    )
    parser.add_argument(
        "--clone-groups",
        type=Path,
        help=(
            "Path to breast_cancer_clone_groups.csv. If omitted, clone groups are "
            "built from observations.csv first."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Output CSV path. Defaults to breast_cancer_clone_proportions.csv in "
            "the same directory as the input clone groups file."
        ),
    )
    return parser.parse_args()


def newest_clone_groups_csv(repo_root: Path) -> Path:
    candidates = sorted(
        repo_root.glob(f"output_runs/*/csv/{DEFAULT_INPUT_NAME}"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(f"No {DEFAULT_INPUT_NAME} found under output_runs/*/csv/")
    return candidates[0]


def parse_signature_key(signature_key: str, timepoint_count: int) -> tuple[str, ...]:
    if signature_key == "implicit:founder":
        if timepoint_count <= 0:
            return tuple()
        return tuple(["present", *["unknown"] * (timepoint_count - 1)])

    states_by_index: dict[int, str] = {}
    for part in signature_key.split("|"):
        if ":" not in part:
            continue
        label, state = part.split(":", 1)
        if label.startswith("t") and label[1:].isdigit():
            states_by_index[int(label[1:])] = state

    if not states_by_index:
        return tuple(["present", *["unknown"] * max(0, timepoint_count - 1)])

    max_index = max(states_by_index)
    states = []
    for index in range(max(max_index + 1, timepoint_count)):
        if index == 0:
            states.append(states_by_index.get(index, "present"))
        else:
            states.append(states_by_index.get(index, "unknown"))
    return tuple(states)


def load_clone_groups(path: Path) -> dict[str, list[dict]]:
    patients: dict[str, list[dict]] = {}

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            patient_id = row["patient_id"].strip()
            timepoint_dates = [part.strip() for part in row["timepoint_dates"].split(";") if part.strip()]
            timepoint_count = len(timepoint_dates)
            parsed_row = {
                "patient_id": patient_id,
                "clone_id": row["clone_id"].strip(),
                "clone_type": row["clone_type"].strip(),
                "parent_clone_id": row["parent_clone_id"].strip(),
                "timepoint_dates": timepoint_dates,
                "signature_key": row["signature_key"].strip(),
                "gene_count": int(row["gene_count"]),
                "genes": row["genes"].strip(),
                "cumulative_gene_count": int(row["cumulative_gene_count"]),
            }
            parsed_row["states"] = parse_signature_key(parsed_row["signature_key"], timepoint_count)
            patients.setdefault(patient_id, []).append(parsed_row)

    return patients


def round_percentages(values: list[float]) -> list[float]:
    if not values:
        return []

    total = sum(values)
    if total <= 0:
        equal_units = ROUNDING_UNITS // len(values)
        remainder = ROUNDING_UNITS - (equal_units * len(values))
        units = [equal_units] * len(values)
        for index in range(remainder):
            units[index] += 1
        return [unit / 100 for unit in units]

    normalized = [(value / total) * 100 for value in values]
    raw_units = [value * 100 for value in normalized]
    base_units = [math.floor(unit) for unit in raw_units]
    remainder = ROUNDING_UNITS - sum(base_units)

    ranked_indices = sorted(
        range(len(values)),
        key=lambda index: (raw_units[index] - base_units[index], normalized[index], -index),
        reverse=True,
    )
    for index in ranked_indices[:remainder]:
        base_units[index] += 1

    return [unit / 100 for unit in base_units]


def midpoint_state(clone: dict) -> str:
    return clone["states"][1] if len(clone["states"]) > 2 else "unknown"


def stable_band(p0: float) -> tuple[float, float]:
    delta = max(2.0, min(3.5, p0 * 0.08))
    return max(0.0, p0 - delta), min(100.0, p0 + delta)


def donor_floor(clone: dict, start: dict[str, float], final: dict[str, float]) -> float:
    clone_id = clone["clone_id"]
    p0 = start[clone_id]
    p2 = final[clone_id]
    state = midpoint_state(clone)

    if state == "unknown":
        return 0.5
    if state == "increasing":
        return max(0.5, p0 + max(1.0, 0.40 * (p2 - p0)))
    if state == "decreasing":
        return max(0.5, p2 + max(1.0, 0.15 * (p0 - p2)))
    if state == "stable":
        return stable_band(p0)[0]
    return 0.5


def adjust_midpoint_stable_clones(
    clones: list[dict],
    proportions: dict[str, float],
    start: dict[str, float],
    final: dict[str, float],
) -> dict[str, float]:
    adjusted = dict(proportions)
    donors_by_priority = sorted(
        clones,
        key=lambda clone: (
            {"unknown": 0, "increasing": 1, "decreasing": 2, "stable": 3}.get(midpoint_state(clone), 4),
            -start[clone["clone_id"]],
            clone["clone_id"],
        ),
    )

    for clone in clones:
        if midpoint_state(clone) != "stable":
            continue

        clone_id = clone["clone_id"]
        low, high = stable_band(start[clone_id])

        if adjusted[clone_id] < low:
            need = round(low - adjusted[clone_id], 2)
            for donor in donors_by_priority:
                donor_id = donor["clone_id"]
                if donor_id == clone_id or need <= 0:
                    continue

                floor = donor_floor(donor, start, final)
                available = round(adjusted[donor_id] - floor, 2)
                if available <= 0:
                    continue

                transfer = min(need, available)
                transfer = math.floor(transfer * 100) / 100
                if transfer <= 0:
                    continue

                adjusted[donor_id] = round(adjusted[donor_id] - transfer, 2)
                adjusted[clone_id] = round(adjusted[clone_id] + transfer, 2)
                need = round(need - transfer, 2)

        elif adjusted[clone_id] > high:
            excess = round(adjusted[clone_id] - high, 2)
            recipients = [donor for donor in donors_by_priority if donor["clone_id"] != clone_id]
            for recipient in recipients:
                if excess <= 0:
                    break
                recipient_id = recipient["clone_id"]
                adjusted[recipient_id] = round(adjusted[recipient_id] + 0.01, 2)
                adjusted[clone_id] = round(adjusted[clone_id] - 0.01, 2)
                excess = round(excess - 0.01, 2)

    return adjusted


def midpoint_buffer_clone(clones: list[dict], excluded_ids: set[str]) -> dict | None:
    founding = next(
        (
            clone
            for clone in clones
            if clone["clone_type"] == "founding" and clone["clone_id"] not in excluded_ids
        ),
        None,
    )
    if founding is not None:
        return founding

    fallback_order = {"unknown": 0, "stable": 1, "increasing": 2, "decreasing": 3}
    candidates = [clone for clone in clones if clone["clone_id"] not in excluded_ids]
    if not candidates:
        return None

    return min(
        candidates,
        key=lambda clone: (
            fallback_order.get(midpoint_state(clone), 4),
            -clone["gene_count"],
            clone["clone_id"],
        ),
    )


def adjust_midpoint_decreasing_clones(
    clones: list[dict],
    proportions: dict[str, float],
    start: dict[str, float],
) -> dict[str, float]:
    adjusted = dict(proportions)
    for clone in clones:
        if midpoint_state(clone) != "decreasing":
            continue

        clone_id = clone["clone_id"]
        ceiling = round(start[clone_id], 2)
        if adjusted[clone_id] <= ceiling:
            continue

        buffer_clone = midpoint_buffer_clone(clones, {clone_id})
        if buffer_clone is None:
            continue

        buffer_id = buffer_clone["clone_id"]
        excess = round(adjusted[clone_id] - ceiling, 2)
        adjusted[clone_id] = ceiling
        adjusted[buffer_id] = round(adjusted[buffer_id] + excess, 2)

    return adjusted


def adjust_final_stable_clones(
    clones: list[dict],
    final: dict[str, float],
    start: dict[str, float],
) -> dict[str, float]:
    stable_ids = {
        clone["clone_id"]
        for clone in clones
        if clone["states"] and clone["states"][-1] == "stable"
    }
    if not stable_ids:
        return dict(final)

    adjusted = dict(final)
    clamped_stable: dict[str, float] = {}
    for clone in clones:
        clone_id = clone["clone_id"]
        if clone_id not in stable_ids:
            continue
        low, high = stable_band(start[clone_id])
        clamped_stable[clone_id] = min(max(adjusted[clone_id], low), high)

    fixed_low_ids = {
        clone["clone_id"]
        for clone in clones
        if clone["clone_id"] not in stable_ids
        and clone["states"]
        and clone["states"][-1] == "decreasing"
    }
    fixed_low = {clone_id: adjusted[clone_id] for clone_id in fixed_low_ids}
    flexible_clones = [
        clone
        for clone in clones
        if clone["clone_id"] not in stable_ids and clone["clone_id"] not in fixed_low_ids
    ]
    remaining = max(0.0, 100.0 - sum(clamped_stable.values()) - sum(fixed_low.values()))

    if not flexible_clones:
        combined = dict(clamped_stable)
        combined.update(fixed_low)
        ordered = [combined[clone["clone_id"]] for clone in clones]
        rounded = round_percentages(ordered)
        return {clone["clone_id"]: value for clone, value in zip(clones, rounded)}

    flexible_total = sum(adjusted[clone["clone_id"]] for clone in flexible_clones)
    if flexible_total <= 0:
        even_share = remaining / len(flexible_clones)
        for clone in flexible_clones:
            adjusted[clone["clone_id"]] = even_share
    else:
        scale = remaining / flexible_total
        for clone in flexible_clones:
            clone_id = clone["clone_id"]
            adjusted[clone_id] = adjusted[clone_id] * scale

    for clone_id, value in clamped_stable.items():
        adjusted[clone_id] = value
    for clone_id, value in fixed_low.items():
        adjusted[clone_id] = value

    ordered = [adjusted[clone["clone_id"]] for clone in clones]
    rounded = round_percentages(ordered)
    return {clone["clone_id"]: value for clone, value in zip(clones, rounded)}


def type_factor(clone_type: str, endpoint: str) -> float:
    if endpoint == "start":
        return {"founding": 1.7, "branch": 1.2, "late": 1.0}.get(clone_type, 1.0)
    return {"founding": 1.25, "branch": 1.05, "late": 1.1}.get(clone_type, 1.0)


def has_future_increase(states: tuple[str, ...]) -> bool:
    return any(state == "increasing" for state in states[1:])


def has_future_decrease(states: tuple[str, ...]) -> bool:
    return any(state == "decreasing" for state in states[1:])


def start_weight(clone: dict) -> float:
    final_state = clone["states"][-1] if clone["states"] else "unknown"
    state_factor = {
        "decreasing": 2.1,
        "stable": 1.5,
        "unknown": 1.35,
        "increasing": 0.4,
    }.get(final_state, 1.2)
    return max(0.1, clone["gene_count"] * type_factor(clone["clone_type"], "start") * state_factor)


def final_weight(clone: dict) -> float:
    final_state = clone["states"][-1] if clone["states"] else "unknown"
    state_factor = {
        "increasing": 2.25,
        "stable": 1.5,
        "unknown": 1.2,
        "decreasing": 0.4,
    }.get(final_state, 1.0)
    return max(0.1, clone["gene_count"] * type_factor(clone["clone_type"], "final") * state_factor)


def low_share(rng: random.Random, clone: dict) -> float:
    gene_bonus = min(0.45, 0.12 * max(0, clone["gene_count"] - 1))
    return min(LOW_CLONE_MAX, LOW_CLONE_MIN + gene_bonus + rng.random() * 2.6)


def endpoint_anchor_clone(clones: list[dict], endpoint: str) -> dict | None:
    if endpoint == "start":
        candidates = [clone for clone in clones if has_future_decrease(clone["states"])]
        weight_fn = start_weight
    else:
        candidates = [clone for clone in clones if has_future_increase(clone["states"])]
        weight_fn = final_weight

    if not candidates:
        return None

    return max(candidates, key=lambda clone: (weight_fn(clone), clone["gene_count"], clone["clone_id"]))


def all_clones_strictly_decreasing(clones: list[dict]) -> bool:
    if not clones:
        return False
    return all(
        clone["states"] and clone["states"][-1] == "decreasing"
        for clone in clones
    )


def allocate_endpoint_percentages(
    patient_id: str, clones: list[dict], endpoint: str
) -> dict[str, float]:
    if endpoint == "start":
        low_ids = {clone["clone_id"] for clone in clones if has_future_increase(clone["states"])}
        weight_fn = start_weight
    else:
        if all_clones_strictly_decreasing(clones):
            founding_clone = next(
                (clone for clone in clones if clone["clone_type"] == "founding"),
                None,
            )
            if founding_clone is not None:
                low_ids = {
                    clone["clone_id"]
                    for clone in clones
                    if clone["clone_id"] != founding_clone["clone_id"]
                }
            else:
                low_ids = {
                    clone["clone_id"]
                    for clone in clones
                    if clone["states"] and clone["states"][-1] == "decreasing"
                }
        else:
            low_ids = {
                clone["clone_id"]
                for clone in clones
                if clone["states"] and clone["states"][-1] == "decreasing"
            }
        weight_fn = final_weight

    low_clones = [clone for clone in clones if clone["clone_id"] in low_ids]
    nonlow_clones = [clone for clone in clones if clone["clone_id"] not in low_ids]

    if not nonlow_clones and low_clones:
        anchor = max(low_clones, key=weight_fn)
        low_ids.remove(anchor["clone_id"])
        low_clones = [clone for clone in clones if clone["clone_id"] in low_ids]
        nonlow_clones = [anchor]
    elif not nonlow_clones:
        low_ids = set()
        low_clones = []
        nonlow_clones = list(clones)

    rng = random.Random(f"{patient_id}:{endpoint}")
    proportions: dict[str, float] = {}
    low_total = 0.0
    for clone in low_clones:
        share = low_share(rng, clone)
        proportions[clone["clone_id"]] = share
        low_total += share

    remaining = max(0.0, 100.0 - low_total)
    weights = [weight_fn(clone) for clone in nonlow_clones]
    anchor_clone = endpoint_anchor_clone(nonlow_clones, endpoint)
    if anchor_clone is not None and len(nonlow_clones) > 1:
        anchor_index = next(
            index for index, clone in enumerate(nonlow_clones) if clone["clone_id"] == anchor_clone["clone_id"]
        )
        other_weights = [weight for index, weight in enumerate(weights) if index != anchor_index]
        if other_weights:
            max_other = max(other_weights)
            weights[anchor_index] = max(weights[anchor_index], max_other + 0.35)

    weight_total = sum(weights)
    if weight_total <= 0:
        even_share = remaining / max(1, len(nonlow_clones))
        for clone in nonlow_clones:
            proportions[clone["clone_id"]] = even_share
    else:
        for clone, weight in zip(nonlow_clones, weights):
            proportions[clone["clone_id"]] = remaining * (weight / weight_total)

    ordered = [proportions[clone["clone_id"]] for clone in clones]
    rounded = round_percentages(ordered)
    return {clone["clone_id"]: value for clone, value in zip(clones, rounded)}


def interpolate_midpoint_percentages(
    patient_id: str, clones: list[dict], start: dict[str, float], final: dict[str, float]
) -> dict[str, float]:
    rng = random.Random(f"{patient_id}:mid")
    raw_values = []
    for clone in clones:
        p0 = start[clone["clone_id"]]
        p2 = final[clone["clone_id"]]
        state = clone["states"][1] if len(clone["states"]) > 2 else "unknown"

        if state == "increasing":
            raw = (0.15 * p0) + (0.85 * p2)
        elif state == "decreasing":
            raw = (0.30 * p0) + (0.70 * p2)
        elif state == "stable":
            raw = p0
        elif state == "unknown":
            raw = (0.70 * p0) + (0.30 * p2)
        else:
            raw = (p0 + p2) / 2

        if state == "stable":
            raw *= 0.995 + (rng.random() * 0.01)
        else:
            raw *= 0.98 + (rng.random() * 0.04)
        raw_values.append(raw)

    rounded = round_percentages(raw_values)
    midpoint = {clone["clone_id"]: value for clone, value in zip(clones, rounded)}
    midpoint = adjust_midpoint_stable_clones(clones, midpoint, start, final)
    return adjust_midpoint_decreasing_clones(clones, midpoint, start)


def proportions_for_patient(patient_id: str, clones: list[dict]) -> dict[str, list[float]]:
    timepoint_count = len(clones[0]["timepoint_dates"]) if clones else 0
    if timepoint_count < 2:
        return {}

    start = allocate_endpoint_percentages(patient_id, clones, "start")
    final = adjust_final_stable_clones(
        clones,
        allocate_endpoint_percentages(patient_id, clones, "final"),
        start,
    )

    if timepoint_count == 2:
        return {
            clone["clone_id"]: [start[clone["clone_id"]], final[clone["clone_id"]]]
            for clone in clones
        }

    if timepoint_count == 3:
        mid = interpolate_midpoint_percentages(patient_id, clones, start, final)
        return {
            clone["clone_id"]: [
                start[clone["clone_id"]],
                mid[clone["clone_id"]],
                final[clone["clone_id"]],
            ]
            for clone in clones
        }

    raise ValueError(
        f"Unsupported number of timepoints for patient {patient_id}: {timepoint_count}. "
        "Expected 2 or 3."
    )


def build_output_rows(
    patients: dict[str, list[dict]], patient_order: list[str] | None = None
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    if patient_order is None:
        ordered_patient_ids = list(patients)
    else:
        ordered_patient_ids = [patient_id for patient_id in patient_order if patient_id in patients]
        ordered_patient_ids.extend(
            patient_id for patient_id in patients if patient_id not in set(ordered_patient_ids)
        )

    for patient_id in ordered_patient_ids:
        clones = patients[patient_id]
        clone_proportions = proportions_for_patient(patient_id, clones)
        if not clone_proportions:
            continue

        timepoint_dates = clones[0]["timepoint_dates"]
        for clone in clones:
            proportions = clone_proportions[clone["clone_id"]]
            row = {
                "patient_id": clone["patient_id"],
                "clone_id": clone["clone_id"],
                "clone_type": clone["clone_type"],
                "parent_clone_id": clone["parent_clone_id"],
                "timepoint_dates": ";".join(clone["timepoint_dates"]),
                "signature_key": clone["signature_key"],
                "gene_count": str(clone["gene_count"]),
                "genes": clone["genes"],
                "t0_date": timepoint_dates[0] if len(timepoint_dates) > 0 else "",
                "t0_vaf_pct": f"{proportions[0]:.2f}" if len(proportions) > 0 else "",
                "t1_date": timepoint_dates[1] if len(timepoint_dates) > 1 else "",
                "t1_vaf_pct": f"{proportions[1]:.2f}" if len(proportions) > 1 else "",
                "t2_date": timepoint_dates[2] if len(timepoint_dates) > 2 else "",
                "t2_vaf_pct": f"{proportions[2]:.2f}" if len(proportions) > 2 else "",
            }
            rows.append(row)

    return rows


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
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
                "t0_date",
                "t0_vaf_pct",
                "t1_date",
                "t1_vaf_pct",
                "t2_date",
                "t2_vaf_pct",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def validate_rows(rows: list[dict[str, str]]) -> None:
    grouped: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        patient_id = row["patient_id"]
        for label in ("t0", "t1", "t2"):
            date = row[f"{label}_date"]
            value = row[f"{label}_vaf_pct"]
            if not date or not value:
                continue
            grouped.setdefault((patient_id, label), []).append(float(value))

    for (patient_id, label), values in grouped.items():
        total = round(sum(values), 2)
        if total != 100.00:
            raise ValueError(
                f"VAFs for patient {patient_id} at {label} sum to {total:.2f}, not 100.00"
            )


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]

    built_from_observations = args.observations is not None or args.clone_groups is None
    clone_groups_path: Path

    if built_from_observations:
        observations_path = args.observations or newest_observations_csv(repo_root)
        if not observations_path.exists():
            raise FileNotFoundError(f"Observations file not found: {observations_path}")

        patients = load_patient_data(observations_path)
        patient_order = patient_order_from_observations(observations_path)
        clone_group_rows = build_clone_rows(patients, patient_order=patient_order)
        clone_groups_path = observations_path.with_name(CLONE_GROUPS_DEFAULT_OUTPUT_NAME)
        write_clone_group_rows(clone_groups_path, clone_group_rows)
    else:
        clone_groups_path = args.clone_groups or newest_clone_groups_csv(repo_root)
        if not clone_groups_path.exists():
            raise FileNotFoundError(f"Clone groups file not found: {clone_groups_path}")

    output_path = args.output or clone_groups_path.with_name(DEFAULT_OUTPUT_NAME)

    patients = load_clone_groups(clone_groups_path)
    observations_for_order = clone_groups_path.with_name("observations.csv")
    patient_order = None
    if observations_for_order.exists():
        patient_order = patient_order_from_observations(observations_for_order)

    rows = build_output_rows(patients, patient_order=patient_order)
    validate_rows(rows)
    write_rows(output_path, rows)

    patient_count = len({row["patient_id"] for row in rows})
    if built_from_observations:
        print(f"Input observations: {observations_path}")
        print(f"Rebuilt clone groups: {clone_groups_path}")
    else:
        print(f"Input clone groups: {clone_groups_path}")
    print(f"Wrote clone proportions: {output_path}")
    print(f"Patients with VAF output: {patient_count}")
    print(f"Total clone rows: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
