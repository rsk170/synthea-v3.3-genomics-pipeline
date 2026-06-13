#!/usr/bin/env python3
"""Create manuscript-ready benchmark tables and plots."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create supplementary table, figure, and draft text from benchmark_summary.csv."
    )
    parser.add_argument(
        "--bench-root",
        type=Path,
        default=Path("output_runs/benchmark_bc_generation_final"),
        help="Benchmark root containing benchmark_summary.csv and machine metadata.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory. Defaults to <bench-root>/analytics.",
    )
    return parser.parse_args()


def read_rows(summary_path: Path) -> list[dict[str, str]]:
    with summary_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return sorted(rows, key=lambda row: int(row["cohort_size"]))


def float_or_none(value: str) -> float | None:
    value = (value or "").strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def fmt_float(value: float | None, digits: int = 2) -> str:
    if value is None or math.isnan(value):
        return ""
    return f"{value:.{digits}f}"


def extract_machine_summary(bench_root: Path) -> dict[str, str]:
    summary: dict[str, str] = {}
    lscpu = bench_root / "machine_lscpu.txt"
    memory = bench_root / "machine_memory.txt"
    java = bench_root / "java_version.txt"

    if lscpu.exists():
        for line in lscpu.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("Model name:"):
                summary["cpu_model"] = line.split(":", 1)[1].strip()
            elif line.startswith("CPU(s):"):
                summary["cpu_threads"] = line.split(":", 1)[1].strip()
            elif line.startswith("Core(s) per socket:"):
                summary["cpu_cores_per_socket"] = line.split(":", 1)[1].strip()

    if memory.exists():
        for line in memory.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("Mem:"):
                fields = line.split()
                if len(fields) >= 2:
                    summary["memory_total"] = fields[1]

    if java.exists():
        first_line = java.read_text(encoding="utf-8", errors="replace").splitlines()
        if first_line:
            summary["java_version"] = first_line[0].strip()

    return summary


def write_table(rows: list[dict[str, str]], output_path: Path) -> None:
    fieldnames = [
        "requested_cohort_size",
        "generated_patient_rows",
        "records_total",
        "records_alive",
        "records_dead",
        "elapsed_minutes",
        "elapsed_hms",
        "peak_ram_gb",
        "peak_ram_mb",
        "csv_storage_mb",
        "total_storage_mb",
        "seed",
        "max_attempts",
        "timed_out",
        "exit_status",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            elapsed_seconds = float_or_none(row.get("elapsed_seconds", ""))
            max_rss_mb = float_or_none(row.get("max_rss_mb", ""))
            writer.writerow(
                {
                    "requested_cohort_size": row.get("cohort_size", ""),
                    "generated_patient_rows": row.get("patients_csv_rows", ""),
                    "records_total": row.get("records_total", ""),
                    "records_alive": row.get("records_alive", ""),
                    "records_dead": row.get("records_dead", ""),
                    "elapsed_minutes": fmt_float(elapsed_seconds / 60 if elapsed_seconds else None),
                    "elapsed_hms": row.get("elapsed_raw", ""),
                    "peak_ram_gb": fmt_float(max_rss_mb / 1024 if max_rss_mb else None),
                    "peak_ram_mb": fmt_float(max_rss_mb),
                    "csv_storage_mb": fmt_float(float_or_none(row.get("storage_csv_mb", ""))),
                    "total_storage_mb": fmt_float(float_or_none(row.get("storage_total_mb", ""))),
                    "seed": row.get("seed", ""),
                    "max_attempts": row.get("max_attempts", ""),
                    "timed_out": row.get("timed_out", ""),
                    "exit_status": row.get("exit_status", ""),
                }
            )


def write_markdown_table(rows: list[dict[str, str]], output_path: Path) -> None:
    headers = [
        "Requested cohort",
        "Generated rows",
        "Wall time (min)",
        "Peak RAM (GB)",
        "CSV storage (MB)",
        "Total storage (MB)",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        elapsed_seconds = float_or_none(row.get("elapsed_seconds", ""))
        max_rss_mb = float_or_none(row.get("max_rss_mb", ""))
        values = [
            row.get("cohort_size", ""),
            row.get("patients_csv_rows", ""),
            fmt_float(elapsed_seconds / 60 if elapsed_seconds else None),
            fmt_float(max_rss_mb / 1024 if max_rss_mb else None),
            fmt_float(float_or_none(row.get("storage_csv_mb", ""))),
            fmt_float(float_or_none(row.get("storage_total_mb", ""))),
        ]
        lines.append("| " + " | ".join(values) + " |")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(rows: list[dict[str, str]], output_path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib.pyplot as plt

    x_values = [int(row["cohort_size"]) for row in rows]
    y_values = [float(row["elapsed_seconds"]) / 60 for row in rows]
    generated = [row.get("patients_csv_rows", "") for row in rows]

    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=180)
    ax.plot(x_values, y_values, marker="o", linewidth=2.0, color="#2166ac")
    ax.scatter(x_values, y_values, s=42, color="#2166ac", zorder=3)

    for x, y, generated_rows in zip(x_values, y_values, generated):
        ax.annotate(
            f"n={generated_rows}",
            (x, y),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
            fontsize=8,
            color="#333333",
        )

    ax.set_xlabel("Requested breast cancer cohort size")
    ax.set_ylabel("Execution time (minutes)")
    ax.set_title("Synthea breast cancer cohort generation benchmark")
    ax.grid(True, axis="y", linewidth=0.6, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks(x_values)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def write_text(rows: list[dict[str, str]], machine: dict[str, str], output_path: Path) -> None:
    sizes = ", ".join(row["cohort_size"] for row in rows)
    min_elapsed = min(float(row["elapsed_seconds"]) / 60 for row in rows)
    max_elapsed = max(float(row["elapsed_seconds"]) / 60 for row in rows)
    max_ram = max(float(row["max_rss_mb"]) / 1024 for row in rows if row.get("max_rss_mb"))
    max_csv = max(float(row["storage_csv_mb"]) for row in rows if row.get("storage_csv_mb"))
    max_total = max(float(row["storage_total_mb"]) for row in rows if row.get("storage_total_mb"))
    max_attempts = sorted({row.get("max_attempts", "") for row in rows if row.get("max_attempts", "")})

    hardware = []
    if machine.get("cpu_model"):
        hardware.append(machine["cpu_model"])
    if machine.get("cpu_threads"):
        hardware.append(f"{machine['cpu_threads']} logical CPUs")
    if machine.get("memory_total"):
        hardware.append(f"{machine['memory_total']} RAM")

    text = (
        "Benchmark draft text\n"
        "====================\n\n"
        "Methods: Breast cancer cohort-generation benchmarks were run with Synthea using "
        "the breast cancer keep filter, an age range of 45-90 years, fixed Synthea and "
        "clinician seeds, and fixed reference/end dates. Runtime and peak resident memory "
        "were measured with GNU /usr/bin/time -v, and output storage was measured from the "
        "generated output directory. "
    )
    if max_attempts:
        text += f"The maximum number of candidate-generation attempts per retained patient was {', '.join(max_attempts)}. "
    if hardware:
        text += f"Benchmarks were performed on {'; '.join(hardware)}. "
    if machine.get("java_version"):
        text += f"The Java runtime was {machine['java_version']}. "

    text += (
        "\n\n"
        f"Results: Requested cohort sizes of {sizes} were evaluated. Wall-clock execution "
        f"time ranged from {min_elapsed:.2f} to {max_elapsed:.2f} minutes. Peak RAM use "
        f"was at most {max_ram:.2f} GB. CSV output storage was at most {max_csv:.2f} MB "
        f"and total output storage was at most {max_total:.2f} MB. The number of retained "
        "patient rows can differ from the requested population size when the keep-filter "
        "and maximum-attempt limit are applied; therefore, both requested size and generated "
        "patient rows are reported in the supplementary table.\n"
    )

    output_path.write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    summary_path = args.bench_root / "benchmark_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(
            f"Missing {summary_path}. Run scripts/summarize_bc_generation_benchmarks.py first."
        )

    output_dir = args.output_dir or args.bench_root / "analytics"
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = [
        row
        for row in read_rows(summary_path)
        if row.get("exit_status") == "0" and row.get("elapsed_seconds")
    ]
    if not rows:
        raise ValueError(f"No completed benchmark rows found in {summary_path}")

    machine = extract_machine_summary(args.bench_root)
    write_table(rows, output_dir / "bc_generation_benchmark_table.csv")
    write_markdown_table(rows, output_dir / "bc_generation_benchmark_table.md")
    write_plot(rows, output_dir / "bc_generation_execution_time.png")
    write_text(rows, machine, output_dir / "benchmark_methods_results_draft.txt")

    print(f"Wrote analytics to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
