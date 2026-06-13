#!/usr/bin/env python3
"""Summarize breast-cancer Synthea benchmark logs into a CSV table."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path


RUN_RE = re.compile(r"^bc_p(?P<size>\d+)_seed(?P<seed>\d+)$")
RECORDS_RE = re.compile(r"Records: total=(\d+), alive=(\d+), dead=(\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a CSV summary from benchmark_bc_generation log folders."
    )
    parser.add_argument(
        "--bench-root",
        type=Path,
        default=Path("output_runs/benchmark_bc_generation_single"),
        help="Benchmark root containing logs/ and run output directories.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional CSV output path. If omitted, writes to stdout.",
    )
    return parser.parse_args()


def parse_elapsed_seconds(value: str) -> float | None:
    value = value.strip()
    parts = value.split(":")
    try:
        if len(parts) == 2:
            minutes, seconds = parts
            return int(minutes) * 60 + float(seconds)
        if len(parts) == 3:
            hours, minutes, seconds = parts
            return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    except ValueError:
        return None
    return None


def parse_status_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def parse_time_log(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        cleaned = line.strip()
        if cleaned.startswith("User time (seconds):"):
            values["user_seconds"] = cleaned.split(":", 1)[1].strip()
        elif cleaned.startswith("System time (seconds):"):
            values["system_seconds"] = cleaned.split(":", 1)[1].strip()
        elif cleaned.startswith("Elapsed (wall clock) time"):
            elapsed = cleaned.split("):", 1)[1].strip()
            values["elapsed_raw"] = elapsed
            elapsed_seconds = parse_elapsed_seconds(elapsed)
            if elapsed_seconds is not None:
                values["elapsed_seconds"] = f"{elapsed_seconds:.2f}"
        elif cleaned.startswith("Maximum resident set size (kbytes):"):
            kb = cleaned.split(":", 1)[1].strip()
            values["max_rss_kb"] = kb
            try:
                values["max_rss_mb"] = f"{int(kb) / 1024:.2f}"
            except ValueError:
                pass
        elif cleaned.startswith("Exit status:"):
            values["time_exit_status"] = cleaned.split(":", 1)[1].strip()
    return values


def parse_stdout(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("Population:"):
            values["population_requested"] = line.split(":", 1)[1].strip()
        elif line.startswith("Seed:"):
            values["synthea_seed"] = line.split(":", 1)[1].strip()
        elif line.startswith("Provider Seed:"):
            values["provider_seed"] = line.split(":", 1)[1].strip()
        elif line.startswith("Reference Time:"):
            values["reference_time"] = line.split(":", 1)[1].strip()
        elif line.startswith("RNG="):
            values["rng_count"] = line.split("=", 1)[1].strip()
        elif line.startswith("Clinician RNG="):
            values["clinician_rng_count"] = line.split("=", 1)[1].strip()
        else:
            records = RECORDS_RE.search(line)
            if records:
                values["records_total"] = records.group(1)
                values["records_alive"] = records.group(2)
                values["records_dead"] = records.group(3)
    return values


def dir_size(path: Path) -> int | None:
    if not path.exists():
        return None
    total = 0
    for root, _dirs, files in os.walk(path):
        for filename in files:
            try:
                total += (Path(root) / filename).stat().st_size
            except OSError:
                pass
    return total


def csv_patient_rows(path: Path) -> int | None:
    patients = path / "csv" / "patients.csv"
    if not patients.exists():
        return None
    with patients.open(encoding="utf-8", errors="replace") as handle:
        return max(sum(1 for _line in handle) - 1, 0)


def summarize_run(bench_root: Path, log_dir: Path) -> dict[str, str]:
    run = log_dir.name
    match = RUN_RE.match(run)
    row: dict[str, str] = {
        "run": run,
        "cohort_size": match.group("size") if match else "",
        "seed": match.group("seed") if match else "",
    }
    row.update(parse_status_file(log_dir / "run_status.txt"))
    row.update(parse_stdout(log_dir / "synthea.stdout.log"))
    row.update(parse_time_log(log_dir / "synthea.time.log"))

    out_dir = bench_root / run
    total_bytes = dir_size(out_dir)
    csv_bytes = dir_size(out_dir / "csv")
    patients = csv_patient_rows(out_dir)

    if total_bytes is not None:
        row["storage_total_bytes"] = str(total_bytes)
        row["storage_total_mb"] = f"{total_bytes / (1024 * 1024):.2f}"
    if csv_bytes is not None:
        row["storage_csv_bytes"] = str(csv_bytes)
        row["storage_csv_mb"] = f"{csv_bytes / (1024 * 1024):.2f}"
    if patients is not None:
        row["patients_csv_rows"] = str(patients)

    if "exit_status" not in row and "time_exit_status" in row:
        row["exit_status"] = row["time_exit_status"]
    if "timed_out" not in row:
        row["timed_out"] = "true" if row.get("exit_status") in {"124", "137"} else "false"
    return row


def main() -> int:
    args = parse_args()
    logs_root = args.bench_root / "logs"
    if not logs_root.exists():
        raise FileNotFoundError(f"No logs directory found: {logs_root}")

    rows = [
        summarize_run(args.bench_root, log_dir)
        for log_dir in sorted(logs_root.iterdir())
        if log_dir.is_dir()
    ]

    fieldnames = [
        "run",
        "cohort_size",
        "seed",
        "exit_status",
        "timed_out",
        "elapsed_seconds",
        "elapsed_raw",
        "max_rss_mb",
        "max_rss_kb",
        "storage_total_mb",
        "storage_csv_mb",
        "patients_csv_rows",
        "records_total",
        "records_alive",
        "records_dead",
        "rng_count",
        "clinician_rng_count",
        "age_range",
        "keep_file",
        "max_attempts",
        "reference_date",
        "end_date",
        "time_limit",
        "command",
    ]

    extra_fields = sorted({key for row in rows for key in row} - set(fieldnames))
    all_fields = [*fieldnames, *extra_fields]

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        handle = args.output.open("w", newline="", encoding="utf-8")
    else:
        handle = None

    output_handle = handle if handle is not None else os.sys.stdout
    writer = csv.DictWriter(output_handle, fieldnames=all_fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

    if handle is not None:
        handle.close()
        print(f"Wrote {len(rows)} rows to {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
