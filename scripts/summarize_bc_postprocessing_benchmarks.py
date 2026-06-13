#!/usr/bin/env python3
"""Summarize breast-cancer post-processing benchmark logs."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from pathlib import Path


STEP_LABELS = {
    "01_clone_groups": "Clone reconstruction",
    "02_clone_proportions": "Clone proportions",
    "03_pruned_observations": "Pruned observations",
    "04_passenger_mutations": "Passenger assignment",
    "05_complete_maf_files": "Complete MAF files",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create post-processing benchmark tables, plot, and draft text."
    )
    parser.add_argument(
        "--bench-root",
        type=Path,
        default=Path("output_runs/benchmark_bc_postprocessing"),
        help="Benchmark root containing logs/<run_label>/ folders.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory. Defaults to <bench-root>/analytics.",
    )
    return parser.parse_args()


def parse_elapsed_seconds(value: str) -> float | None:
    parts = value.strip().split(":")
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


def parse_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def parse_step_statuses(path: Path) -> dict[str, str]:
    statuses: dict[str, str] = {}
    if not path.exists():
        return statuses
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("step_exit_status="):
            _, value = line.split("=", 1)
            fields = value.split()
            if len(fields) == 2:
                statuses[fields[0]] = fields[1]
    return statuses


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


def read_du_bytes(path: Path) -> str:
    if not path.exists():
        return ""
    text = path.read_text(encoding="utf-8", errors="replace").strip()
    if not text:
        return ""
    return text.split()[0]


def bytes_to_mb(value: str) -> str:
    if not value:
        return ""
    try:
        return f"{int(value) / (1024 * 1024):.2f}"
    except ValueError:
        return ""


def float_or_zero(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def summarize_steps(bench_root: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    logs_root = bench_root / "logs"
    if not logs_root.exists():
        raise FileNotFoundError(f"No logs directory found: {logs_root}")

    step_rows: list[dict[str, str]] = []
    run_rows: list[dict[str, str]] = []

    for log_dir in sorted(path for path in logs_root.iterdir() if path.is_dir()):
        run_status = parse_key_value_file(log_dir / "run_status.txt")
        statuses = parse_step_statuses(log_dir / "run_status.txt")
        run_label = run_status.get("run_label", log_dir.name)

        for step in STEP_LABELS:
            row = {
                "run_label": run_label,
                "step": step,
                "step_label": STEP_LABELS[step],
                "exit_status": statuses.get(step, ""),
            }
            row.update(parse_time_log(log_dir / f"{step}.time.log"))
            storage_total_bytes = read_du_bytes(log_dir / f"{step}.storage_total_bytes.txt")
            storage_csv_bytes = read_du_bytes(log_dir / f"{step}.storage_csv_bytes.txt")
            storage_maf_bytes = read_du_bytes(log_dir / f"{step}.storage_maf_files_bytes.txt")
            row["storage_total_bytes_after_step"] = storage_total_bytes
            row["storage_total_mb_after_step"] = bytes_to_mb(storage_total_bytes)
            row["storage_csv_bytes_after_step"] = storage_csv_bytes
            row["storage_csv_mb_after_step"] = bytes_to_mb(storage_csv_bytes)
            row["storage_maf_files_bytes_after_step"] = storage_maf_bytes
            row["storage_maf_files_mb_after_step"] = bytes_to_mb(storage_maf_bytes)
            step_rows.append(row)

        completed_step_rows = [
            row for row in step_rows if row["run_label"] == run_label and row.get("elapsed_seconds")
        ]
        total_elapsed = sum(float_or_zero(row.get("elapsed_seconds", "")) for row in completed_step_rows)
        peak_rss = max((float_or_zero(row.get("max_rss_mb", "")) for row in completed_step_rows), default=0.0)
        final_storage_bytes = read_du_bytes(log_dir / "05_complete_maf_files.storage_total_bytes.txt")
        final_csv_bytes = read_du_bytes(log_dir / "05_complete_maf_files.storage_csv_bytes.txt")
        final_maf_bytes = read_du_bytes(log_dir / "05_complete_maf_files.storage_maf_files_bytes.txt")

        run_rows.append(
            {
                "run_label": run_label,
                "source_run": run_status.get("source_run", ""),
                "patients_csv_rows": run_status.get("patients_csv_rows", ""),
                "clone_groups_rows": run_status.get("clone_groups_rows", ""),
                "clone_proportions_rows": run_status.get("clone_proportions_rows", ""),
                "pruned_observations_rows": run_status.get("pruned_observations_rows", ""),
                "assigned_passenger_rows": run_status.get("assigned_passenger_rows", ""),
                "maf_files": run_status.get("maf_files", ""),
                "passenger_pool_mb": bytes_to_mb(run_status.get("passenger_pool_bytes", "")),
                "total_elapsed_seconds": f"{total_elapsed:.2f}",
                "total_elapsed_minutes": f"{total_elapsed / 60:.2f}",
                "peak_step_rss_mb": f"{peak_rss:.2f}",
                "peak_step_rss_gb": f"{peak_rss / 1024:.2f}",
                "final_storage_total_mb": bytes_to_mb(final_storage_bytes),
                "final_storage_csv_mb": bytes_to_mb(final_csv_bytes),
                "final_storage_maf_files_mb": bytes_to_mb(final_maf_bytes),
            }
        )

    return step_rows, run_rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = sorted({key for row in rows for key in row})
    preferred = [
        "run_label",
        "step",
        "step_label",
        "exit_status",
        "elapsed_seconds",
        "elapsed_raw",
        "max_rss_mb",
        "max_rss_kb",
        "storage_total_mb_after_step",
        "storage_csv_mb_after_step",
        "storage_maf_files_mb_after_step",
    ]
    fieldnames = [field for field in preferred if field in fieldnames] + [
        field for field in fieldnames if field not in preferred
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(run_rows: list[dict[str, str]], step_rows: list[dict[str, str]], path: Path) -> None:
    lines = [
        "## Post-processing run summary",
        "",
        "| Run | Patients | Total time (min) | Peak RAM (GB) | Assigned passenger rows | MAF files | Final storage (MB) |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in run_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["run_label"],
                    row["patients_csv_rows"],
                    row["total_elapsed_minutes"],
                    row["peak_step_rss_gb"],
                    row["assigned_passenger_rows"],
                    row["maf_files"],
                    row["final_storage_total_mb"],
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Step summary",
            "",
            "| Step | Time (min) | Peak RAM (GB) | Storage after step (MB) |",
            "| --- | ---: | ---: | ---: |",
        ]
    )
    for row in step_rows:
        elapsed = float_or_zero(row.get("elapsed_seconds", "")) / 60
        rss = float_or_zero(row.get("max_rss_mb", "")) / 1024
        lines.append(
            "| "
            + " | ".join(
                [
                    row["step_label"],
                    f"{elapsed:.2f}",
                    f"{rss:.2f}",
                    row.get("storage_total_mb_after_step", ""),
                ]
            )
            + " |"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(step_rows: list[dict[str, str]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib.pyplot as plt

    rows = [row for row in step_rows if row.get("elapsed_seconds")]
    labels = [row["step_label"] for row in rows]
    minutes = [float(row["elapsed_seconds"]) / 60 for row in rows]

    fig_height = max(3.8, 0.45 * len(labels) + 1.2)
    fig, ax = plt.subplots(figsize=(7.2, fig_height), dpi=180)
    ax.barh(labels, minutes, color="#1b9e77")
    ax.set_xlabel("Execution time (minutes)")
    ax.set_title("Breast cancer post-processing benchmark")
    ax.grid(True, axis="x", linewidth=0.6, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()
    for index, value in enumerate(minutes):
        ax.text(value, index, f" {value:.2f}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def write_text(run_rows: list[dict[str, str]], step_rows: list[dict[str, str]], path: Path) -> None:
    if not run_rows:
        path.write_text("", encoding="utf-8")
        return

    run = run_rows[0]
    slowest = max(step_rows, key=lambda row: float_or_zero(row.get("elapsed_seconds", "")))
    text = (
        "Post-processing benchmark draft text\n"
        "====================================\n\n"
        "Methods: Downstream breast-cancer genomics post-processing was benchmarked "
        "after Synthea generation by copying the generated CSV files into a separate "
        "benchmark workspace and running the clone reconstruction, clone-proportion "
        "assignment, pruned-observation generation, passenger-mutation assignment, "
        "and complete clone-level MAF generation steps sequentially. Each step was "
        "measured with GNU /usr/bin/time -v, and storage was measured after each "
        "step from the benchmark output directory. The original Synthea output was "
        "not modified.\n\n"
        f"Results: For benchmark run {run['run_label']} with {run['patients_csv_rows']} "
        f"patient rows, post-processing took {run['total_elapsed_minutes']} minutes in "
        f"total. The maximum peak resident memory across steps was {run['peak_step_rss_gb']} "
        f"GB. The pipeline produced {run['assigned_passenger_rows']} assigned passenger "
        f"mutation rows and {run['maf_files']} clone-level MAF files. Final benchmark "
        f"output storage was {run['final_storage_total_mb']} MB. The slowest step was "
        f"{slowest['step_label']} ({float_or_zero(slowest.get('elapsed_seconds', '')) / 60:.2f} minutes).\n"
    )
    path.write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir or args.bench_root / "analytics"
    output_dir.mkdir(parents=True, exist_ok=True)

    step_rows, run_rows = summarize_steps(args.bench_root)
    write_csv(output_dir / "postprocessing_step_summary.csv", step_rows)
    write_csv(output_dir / "postprocessing_run_summary.csv", run_rows)
    write_markdown(run_rows, step_rows, output_dir / "postprocessing_benchmark_table.md")
    write_plot(step_rows, output_dir / "postprocessing_step_times.png")
    write_text(run_rows, step_rows, output_dir / "postprocessing_methods_results_draft.txt")

    print(f"Wrote post-processing analytics to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
