#!/usr/bin/env python3
"""Internal coherence checks for synthetic clinical-genomic outputs.

This script is intended for infrastructure testing of synthetic data generated
by the breast-cancer Synthea/genomics workflow. It checks whether clinical CSVs,
clone outputs, MAF files, and optional FASTQ/metadata files are internally
linked in plausible ways. These are not validations of clinical predictive
accuracy, biological realism, or downstream discovery utility.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


STATUSES = ("PASS", "WARNING", "FAIL", "NA")
TIMEPOINT_LABELS = ("t0", "t1", "t2", "baseline", "pre", "post")
GENOMIC_KEYWORDS = (
    "gene",
    "genomic",
    "genom",
    "variant",
    "mutation",
    "clone",
    "clonal",
    "driver",
    "sequenc",
    "tumor",
    "tumour",
)
THERAPY_DRIVER_GENES = {
    "HER2",
    "ERBB2",
    "ESR1",
    "PGR",
    "PIK3CA",
    "TP53",
    "BRCA1",
    "BRCA2",
    "PTEN",
    "CDK4",
    "CDK6",
    "CCND1",
    "CHEK2",
}
THERAPY_TERMS = (
    "trastuzumab",
    "pertuzumab",
    "tamoxifen",
    "anastrozole",
    "letrozole",
    "exemestane",
    "palbociclib",
    "ribociclib",
    "abemaciclib",
    "olaparib",
    "her2",
    "erbb2",
    "estrogen",
    "progesterone",
)
GENE_ALIASES = {
    "HER2": {"HER2", "ERBB2"},
    "ERBB2": {"HER2", "ERBB2"},
}


@dataclass
class Metric:
    category: str
    metric: str
    status: str
    observed_value: str
    expected_or_threshold: str
    details: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate internal clinical-genomic coherence of one generated synthetic "
            "output directory. Reports metrics only; does not modify input files."
        )
    )
    parser.add_argument("--input-dir", type=Path, required=True, help="Generated output directory.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for validation reports.")
    parser.add_argument(
        "--count-fastq-reads",
        action="store_true",
        help="Optionally count FASTQ reads using gzip-aware line counting. This can be slow.",
    )
    return parser.parse_args()


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def status_from_bool(condition: bool, fail_status: str = "FAIL") -> str:
    return "PASS" if condition else fail_status


def add_metric(
    metrics: list[Metric],
    category: str,
    metric: str,
    status: str,
    observed_value: object,
    expected_or_threshold: object,
    details: object,
) -> None:
    status = status if status in STATUSES else "WARNING"
    metrics.append(
        Metric(
            category=category,
            metric=metric,
            status=status,
            observed_value=str(observed_value),
            expected_or_threshold=str(expected_or_threshold),
            details=str(details),
        )
    )


def recursive_files(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*") if path.is_file())


def choose_first(paths: Iterable[Path], patterns: Iterable[str]) -> Path | None:
    normalized_patterns = [pattern.lower() for pattern in patterns]
    for path in paths:
        name = path.name.lower()
        if any(pattern in name for pattern in normalized_patterns):
            return path
    return None


def discover_files(input_dir: Path) -> dict[str, object]:
    files = recursive_files(input_dir)
    csv_like = [p for p in files if p.suffix.lower() in {".csv", ".tsv", ".txt"}]
    mafs = [p for p in files if p.suffix.lower() == ".maf"]
    fastqs = [
        p
        for p in files
        if p.name.lower().endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))
    ]
    manifests = [
        p
        for p in files
        if p.suffix.lower() in {".csv", ".tsv", ".json", ".txt"}
        and any(token in p.name.lower() for token in ("manifest", "metadata", "sequencing"))
    ]
    canonical_clone_mafs = [
        p
        for p in mafs
        if "maf_files" in {part.lower() for part in p.parts}
    ]
    fallback_clone_mafs = [p for p in mafs if p.name.lower().startswith("clone_")]
    clone_mafs = canonical_clone_mafs if canonical_clone_mafs else fallback_clone_mafs
    report_mafs = clone_mafs if clone_mafs else mafs

    return {
        "all_files": files,
        "patients": choose_first(csv_like, ("patients.csv", "persons.csv")),
        "conditions": choose_first(csv_like, ("conditions.csv",)),
        "procedures": choose_first(csv_like, ("procedures.csv",)),
        "medications": choose_first(csv_like, ("medications.csv",)),
        "observations": choose_first(
            [p for p in csv_like if p.name.lower() == "observations.csv"], ("observations.csv",)
        ),
        "variant_observations": choose_first(
            csv_like,
            (
                "observations_pruned_by_clone_vaf",
                "variant_observations",
                "enhanced_observations",
                "pruned_observations",
            ),
        ),
        "clone_groups": choose_first(csv_like, ("breast_cancer_clone_groups", "clone_groups")),
        "clone_proportions": choose_first(csv_like, ("breast_cancer_clone_proportions", "clone_proportions")),
        "assigned_passengers": choose_first(csv_like, ("breast_cancer_assigned_passenger", "assigned_passenger")),
        "maf_files": report_mafs,
        "clone_maf_files": clone_mafs,
        "fastq_files": fastqs,
        "manifest_files": manifests,
    }


def infer_sep(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".maf"}:
        return "\t"
    return ","


def safe_read_table(path: Path | None, nrows: int | None = None) -> tuple[pd.DataFrame | None, str | None]:
    if path is None:
        return None, "missing file"
    try:
        if path.stat().st_size == 0:
            return pd.DataFrame(), None
    except OSError as exc:
        return None, str(exc)
    try:
        return (
            pd.read_csv(path, sep=infer_sep(path), dtype=str, keep_default_na=False, nrows=nrows),
            None,
        )
    except Exception as exc:  # noqa: BLE001 - report file-level validation failures, do not crash.
        return None, str(exc)


def normalized_columns(df: pd.DataFrame | None) -> dict[str, str]:
    if df is None:
        return {}
    return {normalize_name(col): col for col in df.columns}


def find_col(df: pd.DataFrame | None, candidates: Iterable[str]) -> str | None:
    cmap = normalized_columns(df)
    for candidate in candidates:
        normalized = normalize_name(candidate)
        if normalized in cmap:
            return cmap[normalized]
    return None


def find_cols_containing(df: pd.DataFrame | None, tokens: Iterable[str]) -> list[str]:
    if df is None:
        return []
    normalized_tokens = [normalize_name(token) for token in tokens]
    result = []
    for col in df.columns:
        ncol = normalize_name(col)
        if any(token in ncol for token in normalized_tokens):
            result.append(col)
    return result


def nonempty_set(series: pd.Series | None) -> set[str]:
    if series is None:
        return set()
    return {str(value).strip() for value in series.dropna().astype(str) if str(value).strip()}


def patient_column(df: pd.DataFrame | None) -> str | None:
    return find_col(
        df,
        (
            "PATIENT",
            "patient",
            "patient_id",
            "person",
            "person_id",
            "Id",
            "ID",
            "assigned_patient_id",
        ),
    )


def clone_column(df: pd.DataFrame | None) -> str | None:
    return find_col(df, ("clone_id", "clone", "Clone_ID", "clone"))


def gene_column(df: pd.DataFrame | None) -> str | None:
    return find_col(df, ("gene", "genes", "Hugo_Symbol", "hugo_symbol", "symbol"))


def split_genes(value: object) -> set[str]:
    if value is None:
        return set()
    text = str(value).strip()
    if not text:
        return set()
    return {gene.strip().upper() for gene in re.split(r"[;,|]\s*", text) if gene.strip()}


def gene_equivalents(gene: str) -> set[str]:
    normalized = str(gene).strip().upper()
    return GENE_ALIASES.get(normalized, {normalized})


def gene_present_in_set(gene: str, genes: set[str]) -> bool:
    normalized_genes = {str(value).strip().upper() for value in genes}
    return bool(gene_equivalents(gene) & normalized_genes)


def extract_genes_from_dataframe(df: pd.DataFrame | None, known_genes: set[str]) -> set[str]:
    if df is None or df.empty:
        return set()
    found: set[str] = set()
    gcol = gene_column(df)
    if gcol:
        for value in df[gcol]:
            found.update(split_genes(value))
    if not known_genes:
        return found
    text_cols = [
        col
        for col in df.columns
        if normalize_name(col) in {"description", "value", "code", "category", "type"}
        or any(token in normalize_name(col) for token in ("gene", "variant", "genomic", "mutation"))
    ]
    pattern = re.compile(r"\b(" + "|".join(re.escape(gene) for gene in sorted(known_genes)) + r")\b", re.I)
    for col in text_cols:
        sample = df[col].astype(str)
        for value in sample:
            found.update(match.upper() for match in pattern.findall(value))
    return found


def parse_clone_from_maf_path(path: Path) -> tuple[str | None, str | None]:
    patient_id = path.parent.name if path.parent.name.lower() != "maf_files" else None
    match = re.search(r"(clone_\d+)", path.stem, re.I)
    clone_id = match.group(1) if match else None
    return patient_id, clone_id


def load_maf_summaries(maf_paths: list[Path]) -> tuple[list[dict[str, object]], dict[tuple[str, str], set[str]], set[str]]:
    summaries: list[dict[str, object]] = []
    genes_by_patient_clone: dict[tuple[str, str], set[str]] = defaultdict(set)
    all_maf_genes: set[str] = set()
    required = {
        "hugo_symbol": "Hugo_Symbol",
        "chromosome": "Chromosome",
        "start_position": "Start_position",
        "reference_allele": "Reference_Allele",
        "tumor_seq_allele2": "Tumor_Seq_Allele2",
    }

    for path in maf_paths:
        df, error = safe_read_table(path)
        patient_id, clone_id = parse_clone_from_maf_path(path)
        summary: dict[str, object] = {
            "path": str(path),
            "patient_id": patient_id or "",
            "clone_id": clone_id or "",
            "rows": 0,
            "error": error or "",
            "missing_required": [],
            "missing_core_rows": 0,
            "duplicate_exact_variants": 0,
            "passenger_rows": "",
            "driver_rows": "",
        }
        if df is None:
            summaries.append(summary)
            continue
        summary["rows"] = len(df)
        cmap = normalized_columns(df)
        missing = [label for normalized, label in required.items() if normalized not in cmap]
        summary["missing_required"] = missing

        gene_col = find_col(df, ("Hugo_Symbol", "hugo_symbol", "gene"))
        chrom_col = find_col(df, ("Chromosome", "chromosome"))
        start_col = find_col(df, ("Start_position", "Start_Position", "start"))
        ref_col = find_col(df, ("Reference_Allele", "reference_allele", "ref"))
        alt_col = find_col(df, ("Tumor_Seq_Allele2", "alternate_allele", "alt"))
        if all([gene_col, chrom_col, start_col, ref_col, alt_col]):
            core_cols = [gene_col, chrom_col, start_col, ref_col, alt_col]
            missing_core = df[core_cols].eq("").any(axis=1).sum()
            summary["missing_core_rows"] = int(missing_core)
            duplicate_count = int(df.duplicated(subset=core_cols).sum())
            summary["duplicate_exact_variants"] = duplicate_count
        if gene_col:
            genes = {gene.upper() for gene in nonempty_set(df[gene_col])}
            all_maf_genes.update(genes)
            if patient_id and clone_id:
                genes_by_patient_clone[(patient_id, clone_id)].update(genes)
        origin_col = find_col(df, ("mutation_origin", "origin", "variant_origin"))
        if origin_col:
            origins = df[origin_col].astype(str).str.lower()
            summary["passenger_rows"] = int((origins == "passenger").sum())
            summary["driver_rows"] = int((origins == "driver").sum())
        summaries.append(summary)
    return summaries, genes_by_patient_clone, all_maf_genes


def rows_with_keywords(df: pd.DataFrame | None) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    text_cols = [
        col
        for col in df.columns
        if normalize_name(col) in {"description", "value", "category", "code", "type"}
    ]
    if not text_cols:
        return pd.DataFrame()
    combined = df[text_cols].astype(str).agg(" ".join, axis=1).str.lower()
    mask = combined.apply(lambda text: any(keyword in text for keyword in GENOMIC_KEYWORDS))
    return df.loc[mask]


def add_file_discovery_metrics(metrics: list[Metric], files: dict[str, object]) -> None:
    expected = {
        "patients": "patients/persons CSV",
        "conditions": "conditions CSV",
        "procedures": "procedures CSV",
        "medications": "medications CSV",
        "observations": "observations CSV",
        "variant_observations": "variant-aware/pruned observations CSV",
        "clone_groups": "clone group output",
        "clone_proportions": "clone proportion output",
        "assigned_passengers": "assigned passenger mutation TSV",
    }
    for key, label in expected.items():
        path = files.get(key)
        add_metric(
            metrics,
            "file_discovery",
            label,
            "PASS" if path else "NA",
            path or "not found",
            "file present when relevant",
            "Optional project outputs are reported as NA if absent.",
        )
    add_metric(
        metrics,
        "file_discovery",
        "clone-level MAF files",
        "PASS" if files["clone_maf_files"] else "NA",
        len(files["clone_maf_files"]),  # type: ignore[arg-type]
        ">=1 when MAF output was generated",
        "Detected recursively under maf_files/ or clone_*.maf names.",
    )
    add_metric(
        metrics,
        "file_discovery",
        "FASTQ files",
        "PASS" if files["fastq_files"] else "NA",
        len(files["fastq_files"]),  # type: ignore[arg-type]
        ">=1 when sequencing reads were generated",
        "FASTQ checks are optional and are skipped when no FASTQ files are present.",
    )


def clinical_genomic_checks(
    metrics: list[Metric],
    patients_df: pd.DataFrame | None,
    observations_df: pd.DataFrame | None,
    variant_df: pd.DataFrame | None,
    clone_df: pd.DataFrame | None,
    prop_df: pd.DataFrame | None,
    patient_ids: set[str],
    clone_genes: set[str],
    all_maf_genes: set[str],
) -> dict[str, object]:
    variant_patient_ids: set[str] = set()
    variant_gene_rows = pd.DataFrame()

    if variant_df is None:
        add_metric(metrics, "clinical_genomic_linkage", "variant-aware observation patient IDs", "NA", "no variant-aware observations file", "valid patient IDs", "No pruned/enhanced observations file was detected.")
    else:
        pcol = patient_column(variant_df)
        if not pcol:
            add_metric(metrics, "clinical_genomic_linkage", "variant-aware observation patient IDs", "WARNING", "missing patient column", "patient/person identifier column", "Could not detect a patient ID column.")
        else:
            variant_patient_ids = nonempty_set(variant_df[pcol])
            if patient_ids:
                invalid = sorted(variant_patient_ids - patient_ids)
                status = status_from_bool(not invalid)
                add_metric(metrics, "clinical_genomic_linkage", "variant-aware observation patient IDs", status, f"{len(invalid)} invalid of {len(variant_patient_ids)} unique IDs", "0 invalid IDs", f"Invalid examples: {invalid[:5]}")
            else:
                add_metric(metrics, "clinical_genomic_linkage", "variant-aware observation patient IDs", "NA", f"{len(variant_patient_ids)} IDs observed", "patients.csv available", "Patient ID universe unavailable.")

        variant_gene_rows = rows_with_keywords(variant_df)
        if pcol and patient_ids and not variant_gene_rows.empty:
            driver_ids = nonempty_set(variant_gene_rows[pcol])
            invalid = sorted(driver_ids - patient_ids)
            add_metric(metrics, "clinical_genomic_linkage", "driver/genomic observation patient linkage", status_from_bool(not invalid), f"{len(invalid)} invalid of {len(driver_ids)} patients", "0 invalid IDs", f"Rows with genomic keywords: {len(variant_gene_rows)}")
        elif variant_df is not None:
            add_metric(metrics, "clinical_genomic_linkage", "driver/genomic observation patient linkage", "NA", "not enough structured rows/columns", "genomic rows linked to patient IDs", "The observations file did not expose structured driver rows beyond standard Synthea columns.")

    observed_genes = extract_genes_from_dataframe(variant_df, clone_genes | all_maf_genes)
    if observed_genes:
        missing_from_outputs = sorted(observed_genes - clone_genes - all_maf_genes)
        status = status_from_bool(not missing_from_outputs, fail_status="WARNING")
        add_metric(metrics, "clinical_genomic_linkage", "genes in observations represented in clone/MAF outputs", status, f"{len(observed_genes)} genes detected; {len(missing_from_outputs)} not in outputs", "all detected genes appear in clone or MAF output", f"Missing examples: {missing_from_outputs[:10]}")
    else:
        add_metric(metrics, "clinical_genomic_linkage", "genes in observations represented in clone/MAF outputs", "NA", "no structured/recognized genes detected", "detectable gene mentions", "Could not extract gene symbols from observation rows.")

    exact_field_candidates = {
        "gene": ("gene", "Hugo_Symbol", "hugo_symbol"),
        "chromosome": ("chromosome", "Chromosome"),
        "position": ("position", "start_position", "Start_position", "Start_Position"),
        "reference_allele": ("reference_allele", "Reference_Allele", "ref"),
        "alternate_allele": ("alternate_allele", "Tumor_Seq_Allele2", "alt"),
        "clone_id": ("clone_id", "clone"),
    }
    if variant_df is None:
        add_metric(metrics, "clinical_genomic_linkage", "exact variant fields in observations", "NA", "no variant-aware observations file", "gene/chromosome/position/ref/alt/clone columns", "Optional check skipped.")
    else:
        present_cols = {
            name: find_col(variant_df, candidates)
            for name, candidates in exact_field_candidates.items()
        }
        missing = [name for name, col in present_cols.items() if col is None]
        if missing:
            add_metric(metrics, "clinical_genomic_linkage", "exact variant fields in observations", "NA", f"missing fields: {', '.join(missing)}", "explicit gene/chromosome/position/ref/alt/clone columns", "The detected observations file uses standard Synthea observation columns rather than structured variant columns.")
        else:
            columns = [col for col in present_cols.values() if col]
            missing_rows = int(variant_df[columns].eq("").any(axis=1).sum())
            add_metric(metrics, "clinical_genomic_linkage", "exact variant fields in observations", status_from_bool(missing_rows == 0), missing_rows, "0 rows with missing exact variant fields", f"Checked columns: {columns}")

    if observations_df is None or prop_df is None:
        add_metric(metrics, "clinical_genomic_linkage", "sequencing observation dates match clone proportion dates", "NA", "missing observations or clone proportions", "date overlap where structured timepoints exist", "Optional check skipped.")
    else:
        date_col = find_col(observations_df, ("DATE", "date", "time"))
        prop_date_cols = find_cols_containing(prop_df, ("date",))
        if not date_col or not prop_date_cols:
            add_metric(metrics, "clinical_genomic_linkage", "sequencing observation dates match clone proportion dates", "NA", "missing date columns", "date columns in observations and clone proportions", "Optional check skipped.")
        else:
            genomic_rows = rows_with_keywords(observations_df)
            obs_dates = set(str(value).strip() for value in genomic_rows[date_col].astype(str) if str(value).strip())
            prop_dates = set()
            for col in prop_date_cols:
                prop_dates.update(str(value).strip() for value in prop_df[col].astype(str) if str(value).strip())
            overlap = obs_dates & prop_dates
            if prop_dates and obs_dates:
                add_metric(metrics, "clinical_genomic_linkage", "sequencing observation dates match clone proportion dates", "PASS" if overlap else "WARNING", f"{len(overlap)} overlapping dates", ">=1 overlapping structured sequencing/proportion date", f"Observation genomic dates: {len(obs_dates)}; clone proportion dates: {len(prop_dates)}")
            else:
                add_metric(metrics, "clinical_genomic_linkage", "sequencing observation dates match clone proportion dates", "NA", "no comparable dates detected", "sequencing observation dates and clone proportion dates", "Standard observation rows may not expose sequencing timepoints directly.")

    return {
        "variant_patient_ids": variant_patient_ids,
        "observed_genes": observed_genes,
    }


def clone_checks(
    metrics: list[Metric],
    clone_df: pd.DataFrame | None,
    prop_df: pd.DataFrame | None,
    clone_maf_paths: list[Path],
) -> tuple[dict[tuple[str, str], set[str]], set[str], int]:
    clone_genes_by_patient_clone: dict[tuple[str, str], set[str]] = defaultdict(set)
    all_clone_genes: set[str] = set()
    clone_count = 0

    if clone_df is None:
        add_metric(metrics, "clone_consistency", "clone definition file", "NA", "not found", "clone group output present", "Clone checks requiring clone definitions were skipped.")
        return clone_genes_by_patient_clone, all_clone_genes, clone_count

    pcol = patient_column(clone_df)
    ccol = clone_column(clone_df)
    type_col = find_col(clone_df, ("clone_type", "type"))
    parent_col = find_col(clone_df, ("parent_clone_id", "parent", "parent_id"))
    genes_col = find_col(clone_df, ("genes", "gene"))

    if not pcol or not ccol:
        add_metric(metrics, "clone_consistency", "clone definition required columns", "FAIL", "missing patient_id or clone_id", "patient and clone columns present", f"Columns: {list(clone_df.columns)}")
        return clone_genes_by_patient_clone, all_clone_genes, clone_count

    clone_count = len(clone_df)
    for _, row in clone_df.iterrows():
        patient_id = str(row[pcol]).strip()
        clone_id = str(row[ccol]).strip()
        genes = split_genes(row[genes_col]) if genes_col else set()
        clone_genes_by_patient_clone[(patient_id, clone_id)].update(genes)
        all_clone_genes.update(genes)

    if type_col:
        missing_founding = []
        for patient_id, group in clone_df.groupby(pcol):
            clone_types = {str(value).lower().strip() for value in group[type_col]}
            if "founding" not in clone_types:
                missing_founding.append(patient_id)
        add_metric(metrics, "clone_consistency", "patients with founding clone", status_from_bool(not missing_founding), f"{len(missing_founding)} patients missing founding clone", "0", f"Examples: {missing_founding[:5]}")
    else:
        add_metric(metrics, "clone_consistency", "patients with founding clone", "NA", "clone_type column not found", "clone_type/founding information", "Optional check skipped.")

    if type_col and parent_col:
        invalid_parent = []
        clones_by_patient = {
            patient_id: set(group[ccol].astype(str).str.strip())
            for patient_id, group in clone_df.groupby(pcol)
        }
        for _, row in clone_df.iterrows():
            clone_type = str(row[type_col]).lower().strip()
            if clone_type in {"late", "branch", "subclone", "subclonal"}:
                patient_id = str(row[pcol]).strip()
                parent = str(row[parent_col]).strip()
                if not parent or parent not in clones_by_patient.get(patient_id, set()):
                    invalid_parent.append((patient_id, str(row[ccol]).strip(), parent))
        add_metric(metrics, "clone_consistency", "late/subclonal parent clone IDs", status_from_bool(not invalid_parent), f"{len(invalid_parent)} invalid parent links", "0", f"Examples: {invalid_parent[:5]}")
    else:
        add_metric(metrics, "clone_consistency", "late/subclonal parent clone IDs", "NA", "parent/type columns not found", "parent clone IDs where present", "Optional check skipped.")

    clone_keys = set(clone_genes_by_patient_clone)
    if prop_df is not None:
        ppcol = patient_column(prop_df)
        pccol = clone_column(prop_df)
        if ppcol and pccol:
            prop_keys = {(str(row[ppcol]).strip(), str(row[pccol]).strip()) for _, row in prop_df.iterrows()}
            missing = sorted(prop_keys - clone_keys)
            add_metric(metrics, "clone_consistency", "clone IDs in proportions present in clone definitions", status_from_bool(not missing), f"{len(missing)} missing of {len(prop_keys)} proportion clone IDs", "0 missing clone IDs", f"Examples: {missing[:5]}")
        else:
            add_metric(metrics, "clone_consistency", "clone IDs in proportions present in clone definitions", "NA", "missing patient/clone columns", "patient and clone columns in proportions", "Optional check skipped.")
    else:
        add_metric(metrics, "clone_consistency", "clone IDs in proportions present in clone definitions", "NA", "clone proportions not found", "clone proportions available", "Optional check skipped.")

    maf_keys = set()
    for path in clone_maf_paths:
        patient_id, clone_id = parse_clone_from_maf_path(path)
        if patient_id and clone_id:
            maf_keys.add((patient_id, clone_id))
    if maf_keys:
        missing = sorted(maf_keys - clone_keys)
        add_metric(metrics, "clone_consistency", "clone IDs in MAF filenames present in clone definitions", status_from_bool(not missing, fail_status="WARNING"), f"{len(missing)} missing of {len(maf_keys)} MAF clone IDs", "0 missing clone IDs", f"Examples: {missing[:5]}")
    else:
        add_metric(metrics, "clone_consistency", "clone IDs in MAF filenames present in clone definitions", "NA", "no clone MAF filenames detected", "clone-level MAF files", "Optional check skipped.")

    if prop_df is None:
        add_metric(metrics, "clone_consistency", "clone proportions are non-negative", "NA", "clone proportions not found", "non-negative VAF/proportion columns", "Optional check skipped.")
        add_metric(metrics, "clone_consistency", "clone proportions sum per patient/timepoint", "NA", "clone proportions not found", "approximately 100% or 1.0", "Optional check skipped.")
        return clone_genes_by_patient_clone, all_clone_genes, clone_count

    value_cols = [
        col for col in prop_df.columns
        if normalize_name(col).endswith(("vaf_pct", "proportion", "fraction"))
        or "vaf" in normalize_name(col)
    ]
    numeric_values = {}
    negative_cells = 0
    zero_active_cells = 0
    for col in value_cols:
        values = pd.to_numeric(prop_df[col], errors="coerce")
        numeric_values[col] = values
        negative_cells += int((values.dropna() < 0).sum())
        date_col = find_col(prop_df, (col.replace("vaf_pct", "date"), col.replace("proportion", "date")))
        if date_col:
            active = prop_df[date_col].astype(str).str.strip() != ""
            zero_active_cells += int(((values.fillna(0) <= 0) & active).sum())
        else:
            zero_active_cells += int((values.dropna() <= 0).sum())
    add_metric(metrics, "clone_consistency", "clone proportions are non-negative", status_from_bool(negative_cells == 0), negative_cells, "0 negative values", f"Checked columns: {value_cols}")
    add_metric(metrics, "clone_consistency", "active clone proportions are non-zero", status_from_bool(zero_active_cells == 0, fail_status="WARNING"), zero_active_cells, "0 zero/negative active clone proportions", "Active is inferred from non-empty timepoint date columns when available.")

    ppcol = patient_column(prop_df)
    if ppcol and value_cols:
        bad_sums = []
        sums = []
        for patient_id, group in prop_df.groupby(ppcol):
            for col, values in numeric_values.items():
                group_values = values.loc[group.index].dropna()
                if group_values.empty:
                    continue
                total = float(group_values.sum())
                sums.append(total)
                expected = 100.0 if total > 1.5 else 1.0
                tolerance = 1.0 if expected == 100.0 else 0.01
                if abs(total - expected) > tolerance:
                    bad_sums.append((patient_id, col, round(total, 4)))
        status = status_from_bool(not bad_sums, fail_status="WARNING")
        observed = f"{len(bad_sums)} bad sums across {len(sums)} patient/timepoint sums"
        expected = "100 +/- 1 percentage point or 1.0 +/- 0.01"
        add_metric(metrics, "clone_consistency", "clone proportions sum per patient/timepoint", status, observed, expected, f"Examples: {bad_sums[:5]}")
    else:
        add_metric(metrics, "clone_consistency", "clone proportions sum per patient/timepoint", "NA", "missing patient/proportion columns", "patient and numeric clone proportion columns", "Optional check skipped.")

    return clone_genes_by_patient_clone, all_clone_genes, clone_count


def maf_checks(
    metrics: list[Metric],
    maf_summaries: list[dict[str, object]],
    clone_genes_by_patient_clone: dict[tuple[str, str], set[str]],
    genes_by_patient_clone: dict[tuple[str, str], set[str]],
) -> None:
    if not maf_summaries:
        add_metric(metrics, "variant_maf_consistency", "MAF files discovered", "NA", "0", ">=1 when MAF output generated", "MAF-dependent checks skipped.")
        return

    missing_required_files = [s for s in maf_summaries if s.get("missing_required")]
    add_metric(metrics, "variant_maf_consistency", "required MAF columns present", status_from_bool(not missing_required_files), f"{len(missing_required_files)} files missing required columns of {len(maf_summaries)}", "0 files missing required columns", f"Examples: {[(Path(str(s['path'])).name, s['missing_required']) for s in missing_required_files[:5]]}")

    missing_core_rows = sum(int(s.get("missing_core_rows") or 0) for s in maf_summaries)
    total_rows = sum(int(s.get("rows") or 0) for s in maf_summaries)
    add_metric(metrics, "variant_maf_consistency", "MAF rows have gene symbols and coordinates", status_from_bool(missing_core_rows == 0, fail_status="WARNING"), f"{missing_core_rows} rows with missing core fields of {total_rows}", "0 missing core rows", "Core fields: gene, chromosome, start position, reference allele, alternate allele.")

    duplicate_rows = sum(int(s.get("duplicate_exact_variants") or 0) for s in maf_summaries)
    add_metric(metrics, "variant_maf_consistency", "duplicated exact variants within clone MAFs", status_from_bool(duplicate_rows == 0, fail_status="WARNING"), duplicate_rows, "0 duplicated exact variants", "Duplicates are checked within each clone MAF by gene/chrom/start/ref/alt.")

    missing_clone_genes = []
    for key, clone_genes in clone_genes_by_patient_clone.items():
        maf_genes = genes_by_patient_clone.get(key, set())
        for gene in clone_genes:
            if gene and not gene_present_in_set(gene, maf_genes):
                missing_clone_genes.append((key[0], key[1], gene))
    add_metric(metrics, "variant_maf_consistency", "clone definition genes appear in corresponding clone MAF", status_from_bool(not missing_clone_genes, fail_status="WARNING"), f"{len(missing_clone_genes)} missing clone genes", "0 missing genes", f"Examples: {missing_clone_genes[:10]}")

    passenger_known = [s for s in maf_summaries if s.get("passenger_rows") != ""]
    if passenger_known:
        passenger_rows = sum(int(s.get("passenger_rows") or 0) for s in passenger_known)
        add_metric(metrics, "variant_maf_consistency", "passenger mutations present in MAF output", "PASS" if passenger_rows > 0 else "WARNING", passenger_rows, ">0 passenger rows if passengers were assigned", "Uses mutation_origin column where present.")
    else:
        add_metric(metrics, "variant_maf_consistency", "passenger mutations present in MAF output", "NA", "mutation_origin column not found", "mutation_origin/passenger labels", "Cannot distinguish driver and passenger rows.")

    per_clone_counts = {
        (s.get("patient_id", ""), s.get("clone_id", "")): int(s.get("rows") or 0)
        for s in maf_summaries
    }
    per_patient_counts: Counter[str] = Counter()
    for (patient_id, _clone_id), count in per_clone_counts.items():
        if patient_id:
            per_patient_counts[str(patient_id)] += count
    add_metric(metrics, "variant_maf_consistency", "variants per clone", "PASS", f"{len(per_clone_counts)} clone MAFs; min={min(per_clone_counts.values()) if per_clone_counts else 0}; max={max(per_clone_counts.values()) if per_clone_counts else 0}", "reported distribution", "Counts include driver and passenger rows.")
    add_metric(metrics, "variant_maf_consistency", "variants per patient", "PASS", f"{len(per_patient_counts)} patients with MAFs; min={min(per_patient_counts.values()) if per_patient_counts else 0}; max={max(per_patient_counts.values()) if per_patient_counts else 0}", "reported distribution", "Counts include all clone MAF rows per patient.")


def strip_read_marker(path: Path) -> tuple[str, str | None]:
    name = path.name
    lowered = name.lower()
    read = None
    if re.search(r"(^|[_\-.])r?1([_\-.]|$)", lowered):
        read = "R1"
    elif re.search(r"(^|[_\-.])r?2([_\-.]|$)", lowered):
        read = "R2"
    prefix = re.sub(r"([_\-.])r?[12]([_\-.])", r"\1READ\2", lowered)
    prefix = re.sub(r"(_r?[12])(?=\.f(ast)?q)", "_READ", prefix)
    prefix = re.sub(r"(\.fastq|\.fq)(\.gz)?$", "", prefix)
    return prefix, read


def fastq_read_count(path: Path) -> int:
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    line_count = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_count, _line in enumerate(handle, start=1):
            pass
    return line_count // 4


def fastq_checks(metrics: list[Metric], fastq_paths: list[Path], manifest_paths: list[Path], count_reads: bool) -> None:
    if not fastq_paths and not manifest_paths:
        add_metric(metrics, "sequencing_fastq_consistency", "sequencing/FASTQ outputs discovered", "NA", "0 FASTQ files and 0 manifest files", "sequencing outputs present when generated", "FASTQ checks skipped.")
        return

    add_metric(metrics, "sequencing_fastq_consistency", "FASTQ files discovered", "PASS" if fastq_paths else "NA", len(fastq_paths), ">=1 if FASTQs generated", f"Manifest/metadata files detected: {len(manifest_paths)}")

    paired_by_prefix: dict[str, set[str]] = defaultdict(set)
    for path in fastq_paths:
        prefix, read = strip_read_marker(path)
        if read:
            paired_by_prefix[prefix].add(read)
    complete_pairs = [prefix for prefix, reads in paired_by_prefix.items() if {"R1", "R2"} <= reads]
    incomplete_pairs = [prefix for prefix, reads in paired_by_prefix.items() if reads and {"R1", "R2"} - reads]
    add_metric(metrics, "sequencing_fastq_consistency", "R1/R2 files paired by filename convention", status_from_bool(not incomplete_pairs, fail_status="WARNING"), f"{len(complete_pairs)} complete pairs; {len(incomplete_pairs)} incomplete prefixes", "0 incomplete prefixes", f"Incomplete examples: {incomplete_pairs[:5]}")

    normal_paths = [p for p in fastq_paths if "normal" in p.name.lower() or "germline" in p.name.lower()]
    tumor_paths = [p for p in fastq_paths if any(token in p.name.lower() for token in ("tumor", "tumour", "t0", "t1", "t2"))]
    normal_pairs = {
        prefix for prefix, reads in paired_by_prefix.items()
        if {"R1", "R2"} <= reads and any(token in prefix for token in ("normal", "germline"))
    }
    tumor_pairs = {
        prefix for prefix, reads in paired_by_prefix.items()
        if {"R1", "R2"} <= reads and any(token in prefix for token in ("tumor", "tumour", "t0", "t1", "t2"))
    }
    add_metric(metrics, "sequencing_fastq_consistency", "matched normal FASTQ R1/R2 pair exists", "PASS" if normal_pairs else ("WARNING" if normal_paths else "NA"), len(normal_pairs), ">=1 normal R1/R2 pair when FASTQs include normal sample", f"Normal-like FASTQ files: {len(normal_paths)}")
    add_metric(metrics, "sequencing_fastq_consistency", "tumour timepoint FASTQ R1/R2 pairs exist", "PASS" if tumor_pairs else ("WARNING" if tumor_paths else "NA"), len(tumor_pairs), ">=1 tumour/timepoint R1/R2 pair when tumour FASTQs generated", f"Tumour/timepoint-like FASTQ files: {len(tumor_paths)}")

    depth_checked = False
    for manifest in manifest_paths:
        if manifest.suffix.lower() not in {".csv", ".tsv", ".txt"}:
            continue
        df, _error = safe_read_table(manifest)
        if df is None or df.empty:
            continue
        cols = normalized_columns(df)
        tumor_depth_col = next((cols[c] for c in cols if "tumor" in c and "depth" in c), None)
        clone_fraction_col = next((cols[c] for c in cols if ("clone" in c and ("fraction" in c or "proportion" in c or "vaf" in c))), None)
        clone_depth_col = next((cols[c] for c in cols if "clone" in c and "depth" in c), None)
        if tumor_depth_col and clone_fraction_col and clone_depth_col:
            tumor_depth = pd.to_numeric(df[tumor_depth_col], errors="coerce")
            fraction = pd.to_numeric(df[clone_fraction_col], errors="coerce")
            clone_depth = pd.to_numeric(df[clone_depth_col], errors="coerce")
            if fraction.max(skipna=True) and fraction.max(skipna=True) > 1.5:
                fraction = fraction / 100.0
            expected = tumor_depth * fraction
            valid = pd.DataFrame({"expected": expected, "observed": clone_depth}).dropna()
            if not valid.empty:
                rel_error = ((valid["observed"] - valid["expected"]).abs() / valid["expected"].replace(0, pd.NA)).dropna()
                bad = int((rel_error > 0.10).sum())
                add_metric(metrics, "sequencing_fastq_consistency", "target clone depths match tumour depth x clone fraction", status_from_bool(bad == 0, fail_status="WARNING"), f"{bad} rows >10% relative error of {len(valid)} checked", "<=10% relative error", f"Metadata source: {manifest}")
                depth_checked = True
                break
    if not depth_checked:
        add_metric(metrics, "sequencing_fastq_consistency", "target clone depths match tumour depth x clone fraction", "NA", "no suitable depth metadata columns", "tumour depth, clone fraction, clone depth columns", "Optional depth check skipped.")

    if count_reads and fastq_paths:
        counts = {}
        for path in fastq_paths:
            try:
                counts[str(path)] = fastq_read_count(path)
            except Exception as exc:  # noqa: BLE001
                counts[str(path)] = f"ERROR: {exc}"
        numeric_counts = [value for value in counts.values() if isinstance(value, int)]
        add_metric(metrics, "sequencing_fastq_consistency", "FASTQ read counts", "PASS" if numeric_counts else "WARNING", f"{len(numeric_counts)} files counted", "read counting requested by user", json.dumps(dict(list(counts.items())[:10]))[:1000])
    else:
        add_metric(metrics, "sequencing_fastq_consistency", "FASTQ read counts", "NA", "not requested", "--count-fastq-reads", "Read counting is intentionally opt-in because FASTQ files may be large.")


def proxy_checks(
    metrics: list[Metric],
    clone_genes_by_patient_clone: dict[tuple[str, str], set[str]],
    genes_by_patient_clone: dict[tuple[str, str], set[str]],
    variant_patient_ids: set[str],
    medications_df: pd.DataFrame | None,
    observations_df: pd.DataFrame | None,
) -> None:
    patients_with_proxy = set()
    for (patient_id, clone_id), clone_genes in clone_genes_by_patient_clone.items():
        treatment_genes = clone_genes & THERAPY_DRIVER_GENES
        if not treatment_genes:
            continue
        maf_genes = genes_by_patient_clone.get((patient_id, clone_id), set())
        if treatment_genes & maf_genes and (not variant_patient_ids or patient_id in variant_patient_ids):
            patients_with_proxy.add(patient_id)
    if clone_genes_by_patient_clone and genes_by_patient_clone:
        add_metric(metrics, "proxy_integrated_use_case", "patients with therapy-relevant driver, clone assignment, observation, and MAF entry", "PASS" if patients_with_proxy else "WARNING", len(patients_with_proxy), ">=1 patient supports integrated query", "This is a proxy infrastructure query, not a clinical-validity claim.")
    else:
        add_metric(metrics, "proxy_integrated_use_case", "patients with therapy-relevant driver, clone assignment, observation, and MAF entry", "NA", "missing clone or MAF mappings", "clone and MAF outputs available", "Proxy query skipped.")

    therapy_mentions = 0
    if medications_df is not None and not medications_df.empty:
        text = medications_df.astype(str).agg(" ".join, axis=1).str.lower()
        therapy_mentions += int(text.apply(lambda value: any(term in value for term in THERAPY_TERMS)).sum())
    receptor_mentions = 0
    if observations_df is not None and not observations_df.empty:
        text = observations_df.astype(str).agg(" ".join, axis=1).str.lower()
        receptor_mentions += int(text.apply(lambda value: any(term in value for term in ("her2", "erbb2", "estrogen", "progesterone", " er ", " pr "))).sum())
    if therapy_mentions or receptor_mentions:
        gene_support = sum(
            1 for genes in clone_genes_by_patient_clone.values()
            if genes & {"HER2", "ERBB2", "ESR1", "PGR"}
        )
        add_metric(metrics, "proxy_integrated_use_case", "HER2/ER/PR clinical terms have genomic context", "PASS" if gene_support else "WARNING", f"{therapy_mentions} therapy rows; {receptor_mentions} receptor-like rows; {gene_support} clone gene contexts", "therapy/receptor mentions plus HER2/ER/PR-related genes where present", "Proxy check only; terminology is detected by simple text search.")
    else:
        add_metric(metrics, "proxy_integrated_use_case", "HER2/ER/PR clinical terms have genomic context", "NA", "no therapy/receptor terms detected", "therapy/receptor fields or medication names", "Required terms were not found in medications/observations.")


def write_outputs(metrics: list[Metric], output_dir: Path, summary: dict[str, object]) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "tsv": output_dir / "coherence_metrics.tsv",
        "md": output_dir / "coherence_metrics.md",
        "json": output_dir / "coherence_metrics.json",
    }
    rows = [asdict(metric) for metric in metrics]

    with paths["tsv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["category", "metric", "status", "observed_value", "expected_or_threshold", "details"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    status_counts = Counter(metric.status for metric in metrics)
    md_lines = [
        "# Synthetic Output Coherence Metrics",
        "",
        "These are internal coherence checks for infrastructure-testing synthetic data. They do not validate clinical predictive accuracy or biological discovery utility.",
        "",
        "## Summary",
        "",
    ]
    for key, value in summary.items():
        md_lines.append(f"- **{key}**: {value}")
    md_lines.append("- **metric_status_counts**: " + ", ".join(f"{status}={status_counts.get(status, 0)}" for status in STATUSES))
    md_lines.extend(
        [
            "",
            "## Metrics",
            "",
            "| Category | Metric | Status | Observed | Expected/Threshold | Details |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for metric in metrics:
        md_lines.append(
            "| "
            + " | ".join(
                str(value).replace("|", "\\|").replace("\n", " ")
                for value in (
                    metric.category,
                    metric.metric,
                    metric.status,
                    metric.observed_value,
                    metric.expected_or_threshold,
                    metric.details,
                )
            )
            + " |"
        )
    paths["md"].write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    with paths["json"].open("w", encoding="utf-8") as handle:
        json.dump({"summary": summary, "status_counts": dict(status_counts), "metrics": rows}, handle, indent=2)

    return paths


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    metrics: list[Metric] = []
    files = discover_files(input_dir)
    add_file_discovery_metrics(metrics, files)

    patients_df, patients_error = safe_read_table(files["patients"])  # type: ignore[arg-type]
    observations_df, observations_error = safe_read_table(files["observations"])  # type: ignore[arg-type]
    variant_df, variant_error = safe_read_table(files["variant_observations"])  # type: ignore[arg-type]
    clone_df, clone_error = safe_read_table(files["clone_groups"])  # type: ignore[arg-type]
    prop_df, prop_error = safe_read_table(files["clone_proportions"])  # type: ignore[arg-type]
    medications_df, medications_error = safe_read_table(files["medications"])  # type: ignore[arg-type]

    for label, error in (
        ("patients", patients_error),
        ("observations", observations_error),
        ("variant observations", variant_error),
        ("clone groups", clone_error),
        ("clone proportions", prop_error),
        ("medications", medications_error),
    ):
        if error and error != "missing file":
            add_metric(metrics, "file_reading", f"{label} readable", "WARNING", error, "file can be parsed", "Parsing failed; dependent checks may be NA.")

    patient_ids = set()
    if patients_df is not None:
        pcol = patient_column(patients_df)
        if pcol:
            patient_ids = nonempty_set(patients_df[pcol])
            add_metric(metrics, "clinical_genomic_linkage", "patients analysed", "PASS", len(patient_ids), ">=1", f"Patient column: {pcol}")
        else:
            add_metric(metrics, "clinical_genomic_linkage", "patients analysed", "WARNING", "patient ID column not detected", "patient ID column present", f"Columns: {list(patients_df.columns)}")
    else:
        add_metric(metrics, "clinical_genomic_linkage", "patients analysed", "NA", "patients/persons CSV missing", "patients.csv or persons.csv", "Patient-universe checks are limited.")

    clone_genes_by_patient_clone, clone_genes, clone_count = clone_checks(
        metrics,
        clone_df,
        prop_df,
        files["clone_maf_files"],  # type: ignore[arg-type]
    )
    maf_summaries, maf_genes_by_patient_clone, all_maf_genes = load_maf_summaries(
        files["clone_maf_files"]  # type: ignore[arg-type]
    )
    clinical_context = clinical_genomic_checks(
        metrics,
        patients_df,
        observations_df,
        variant_df,
        clone_df,
        prop_df,
        patient_ids,
        clone_genes,
        all_maf_genes,
    )
    maf_checks(metrics, maf_summaries, clone_genes_by_patient_clone, maf_genes_by_patient_clone)
    fastq_checks(
        metrics,
        files["fastq_files"],  # type: ignore[arg-type]
        files["manifest_files"],  # type: ignore[arg-type]
        args.count_fastq_reads,
    )
    proxy_checks(
        metrics,
        clone_genes_by_patient_clone,
        maf_genes_by_patient_clone,
        clinical_context.get("variant_patient_ids", set()),  # type: ignore[arg-type]
        medications_df,
        observations_df,
    )

    timepoint_count = 0
    if prop_df is not None:
        timepoint_count = len([col for col in prop_df.columns if normalize_name(col).endswith("date")])
    summary = {
        "input_dir": str(input_dir),
        "patients_analyzed": len(patient_ids),
        "timepoint_columns_detected": timepoint_count,
        "clones_detected": clone_count,
        "clone_maf_files_detected": len(files["clone_maf_files"]),  # type: ignore[arg-type]
        "fastq_files_detected": len(files["fastq_files"]),  # type: ignore[arg-type]
    }
    paths = write_outputs(metrics, output_dir, summary)
    counts = Counter(metric.status for metric in metrics)

    print("Synthetic output coherence validation complete")
    print(f"Patients analysed: {summary['patients_analyzed']}")
    print(f"Timepoint columns detected: {summary['timepoint_columns_detected']}")
    print(f"Clones detected: {summary['clones_detected']}")
    print(f"Clone MAF files detected: {summary['clone_maf_files_detected']}")
    print(f"FASTQ files detected: {summary['fastq_files_detected']}")
    print("Metric statuses: " + ", ".join(f"{status}={counts.get(status, 0)}" for status in STATUSES))
    print(f"TSV: {paths['tsv']}")
    print(f"Markdown: {paths['md']}")
    print(f"JSON: {paths['json']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
