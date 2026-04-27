#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_COLUMNS = [
    "drug_name",
    "type",
    "civic_molecular_profile_id",
    "civic_molecular_profile_url",
    "resolved_profile_name",
    "match_level",
    "match_context",
    "resolution_note",
    "gene",
    "variant_name",
    "chromosome",
    "start_position",
    "end_position",
    "reference_allele",
    "tumor_seq_allele2",
    "variant_classification",
    "protein_change",
    "cdna_change",
    "sample_id",
    "tumor_type",
    "project_code",
]


BATCH_DUPLICATE_SKIPS = 0


REQUESTS = [
    {
        "id": "104",
        "url": "https://civicdb.org/molecular-profiles/104",
        "resolved_profile_name": "PIK3CA E545K",
        "gene": "PIK3CA",
        "drug_name": "Lapatinib",
        "type": "resistance",
        "profile_kind": "simple",
        "variants": [
            {
                "variant_name": "E545K",
                "protein_change": "E545K",
                "chromosome": "3",
                "start_position": "178936091",
                "end_position": "178936091",
                "reference_allele": "G",
                "tumor_seq_allele2": "A",
                "cdna_change": "NM_006218.3:c.1633G>A",
            }
        ],
        "accepted_protein_changes": ["E545K"],
        "resolution_note": "Resolved directly from CIViC GraphQL API as PIK3CA E545K before MAF search.",
    },
    {
        "id": "37",
        "url": "https://civicdb.org/molecular-profiles/37",
        "resolved_profile_name": "ERBB2 L755_T759del",
        "gene": "ERBB2",
        "drug_name": "Lapatinib",
        "type": "resistance",
        "profile_kind": "simple",
        "variants": [],
        "accepted_protein_changes": ["L755_T759del"],
        "resolution_note": (
            "Resolved directly from CIViC GraphQL API as ERBB2 L755_T759del. "
            "CIViC exposes an exact deletion HGVS (NC_000017.10:g.37880219_37880233del), "
            "but no exact matching local MAF row was found."
        ),
        "not_found_reason": "no exact genomic MAF match",
    },
    {
        "id": "1208",
        "url": "https://civicdb.org/molecular-profiles/1208",
        "resolved_profile_name": "PIK3CA K111N",
        "gene": "PIK3CA",
        "drug_name": "Lapatinib",
        "type": "resistance",
        "profile_kind": "simple",
        "variants": [
            {
                "variant_name": "K111N",
                "protein_change": "K111N",
                "chromosome": "3",
                "start_position": "178916946",
                "end_position": "178916946",
                "reference_allele": "G",
                "tumor_seq_allele2": "C",
                "cdna_change": "NM_006218.3:c.333G>C",
            }
        ],
        "accepted_protein_changes": ["K111N"],
        "resolution_note": "Resolved directly from CIViC GraphQL API as PIK3CA K111N before MAF search.",
    },
    {
        "id": "2204",
        "url": "https://civicdb.org/molecular-profiles/2204",
        "resolved_profile_name": "ERBB2 T798I",
        "gene": "ERBB2",
        "drug_name": "Neratinib",
        "type": "resistance",
        "profile_kind": "simple",
        "variants": [
            {
                "variant_name": "T798I",
                "protein_change": "T798I",
                "chromosome": "17",
                "start_position": "37881064",
                "end_position": "37881064",
                "reference_allele": "C",
                "tumor_seq_allele2": "T",
                "cdna_change": "NM_001289937.1:c.2393C>T",
            }
        ],
        "accepted_protein_changes": ["T798I"],
        "resolution_note": "Resolved directly from CIViC GraphQL API as ERBB2 T798I before MAF search.",
    },
    {
        "id": "634",
        "url": "https://civicdb.org/molecular-profiles/634",
        "resolved_profile_name": "RB1 M695FS*26",
        "gene": "RB1",
        "drug_name": "Palbociclib",
        "type": "resistance",
        "profile_kind": "simple",
        "variants": [],
        "accepted_protein_changes": ["M695FS*26"],
        "resolution_note": (
            "Resolved directly from CIViC GraphQL API as RB1 M695FS*26. "
            "CIViC exposes an exact insertion HGVS (NC_000013.10:g.49033946_49033947insA), "
            "but no exact matching local MAF row was found."
        ),
        "not_found_reason": "no exact genomic MAF match",
    },
    {
        "id": "777",
        "url": "https://civicdb.org/molecular-profiles/777",
        "resolved_profile_name": "RB1 Mutation",
        "gene": "RB1",
        "drug_name": "Ribociclib",
        "type": "resistance",
        "profile_kind": "non_snv_biomarker",
        "variants": [],
        "accepted_protein_changes": ["Mutation"],
        "resolution_note": (
            "Resolved directly from CIViC GraphQL API as the generic biomarker RB1 Mutation. "
            "CIViC exposes only a broad genomic span and no exact reference/alternate allele."
        ),
        "not_found_reason": "not compatible with somatic SNV/indel MAF matching",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Append exact ICGC MAF rows matching the curated CIViC Lapatinib, "
            "Neratinib, Palbociclib, and Ribociclib resistance batch to "
            "civic_breast_cancer_driver_variants_from_maf.csv."
        )
    )
    repo_root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--maf",
        type=Path,
        default=repo_root / "scripts" / "final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf",
        help="Path to the tab-delimited ICGC public MAF.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=repo_root / "scripts" / "civic_breast_cancer_driver_variants_from_maf.csv",
        help="Path to the output CSV to update in place.",
    )
    parser.add_argument(
        "--not-found-log",
        type=Path,
        default=repo_root / "scripts" / "civic_breast_cancer_driver_variants_from_maf.not_found.txt",
        help="Path to the text note listing profiles without exact MAF matches.",
    )
    return parser.parse_args()


def variant_key_from_variant(variant: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (
        variant["chromosome"],
        variant["start_position"],
        variant["end_position"],
        variant["reference_allele"],
        variant["tumor_seq_allele2"],
    )


def variant_key_from_maf(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (
        row["Chromosome"],
        row["Start_position"],
        row["End_position"],
        row["Reference_Allele"],
        row["Tumor_Seq_Allele2"],
    )


def helper_variant_key(row: dict[str, str]) -> str:
    return "|".join(
        [
            row.get("Hugo_Symbol", ""),
            row.get("Chromosome", ""),
            row.get("Start_position", ""),
            row.get("End_position", ""),
            row.get("Reference_Allele", ""),
            row.get("Tumor_Seq_Allele2", ""),
            row.get("Variant_Classification", ""),
        ]
    )


def is_breast_project(project_code: str) -> bool:
    return "breast" in project_code.strip().lower()


def load_existing_rows(output_path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not output_path.exists():
        return [], []
    with output_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader], list(reader.fieldnames or [])


def backfill_row(row: dict[str, str]) -> dict[str, str]:
    output = dict(row)
    for field in REQUIRED_COLUMNS:
        output.setdefault(field, "")

    output["gene"] = output["gene"] or output.get("civic_gene", "") or output.get("Hugo_Symbol", "")
    output["chromosome"] = output["chromosome"] or output.get("Chromosome", "")
    output["start_position"] = output["start_position"] or output.get("Start_position", "")
    output["end_position"] = output["end_position"] or output.get("End_position", "")
    output["reference_allele"] = output["reference_allele"] or output.get("Reference_Allele", "")
    output["tumor_seq_allele2"] = output["tumor_seq_allele2"] or output.get("Tumor_Seq_Allele2", "")
    output["variant_classification"] = output["variant_classification"] or output.get("Variant_Classification", "")
    output["sample_id"] = output["sample_id"] or output.get("Tumor_Sample_Barcode", "")
    output["project_code"] = output["project_code"] or output.get("Project_Code", "")
    output["tumor_type"] = output["tumor_type"] or output["project_code"]
    if not output.get("match_source"):
        if output.get("match_context") == "breast_cancer":
            output["match_source"] = "breast_cancer_maf_exact"
        elif output.get("match_context") == "pan_cancer_fallback":
            output["match_source"] = "pan_cancer_maf_exact"
    if not output.get("match_confidence") and output.get("match_level"):
        output["match_confidence"] = "high"
    output.setdefault("match_level", "")
    output.setdefault("match_context", "")
    output.setdefault("resolution_note", "")
    return output


def load_maf_hits(maf_path: Path) -> tuple[list[str], dict[tuple[str, str, str, str, str], list[dict[str, str]]]]:
    keys = {
        variant_key_from_variant(variant)
        for request in REQUESTS
        for variant in request["variants"]
    }
    hits: dict[tuple[str, str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    with maf_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            key = variant_key_from_maf(row)
            if key in keys:
                hits[key].append(dict(row))
    return fieldnames, hits


def dedupe_key(row: dict[str, str]) -> tuple[str, str, str, str, str, str, str, str]:
    return (
        row.get("chromosome", ""),
        row.get("start_position", ""),
        row.get("end_position", ""),
        row.get("reference_allele", ""),
        row.get("tumor_seq_allele2", ""),
        row.get("sample_id", "") or row.get("Tumor_Sample_Barcode", ""),
        row.get("drug_name", ""),
        row.get("type", ""),
    )


def select_best_exact_match(
    variants: list[dict[str, str]],
    hits: dict[tuple[str, str, str, str, str], list[dict[str, str]]],
) -> tuple[dict[str, str] | None, dict[str, str] | None, str | None]:
    first_pan_match: tuple[dict[str, str], dict[str, str], str] | None = None
    for variant in variants:
        exact_rows = hits.get(variant_key_from_variant(variant), [])
        if not exact_rows:
            continue
        breast_rows = [row for row in exact_rows if is_breast_project(row.get("Project_Code", ""))]
        if breast_rows:
            return variant, breast_rows[0], "breast_cancer"
        if first_pan_match is None:
            first_pan_match = (variant, exact_rows[0], "pan_cancer_fallback")
    if first_pan_match is not None:
        return first_pan_match
    return None, None, None


def build_resolution_note(request: dict[str, object], matched_variant: dict[str, str]) -> str:
    accepted = request.get("accepted_protein_changes")
    base_note = str(request["resolution_note"])
    if accepted:
        accepted_text = ", ".join(accepted)
        return (
            f"{base_note} Accepted protein changes from CIViC: {accepted_text}. "
            f"Saved exact genomic match for {matched_variant['protein_change']}."
        )
    return base_note


def build_output_row(
    request: dict[str, object],
    matched_variant: dict[str, str],
    maf_row: dict[str, str],
    match_context: str,
) -> dict[str, str]:
    legacy_match_source = (
        "breast_cancer_maf_exact" if match_context == "breast_cancer" else "pan_cancer_maf_exact"
    )
    row = dict(maf_row)
    row.update(
        {
            "civic_molecular_profile_id": str(request["id"]),
            "civic_molecular_profile_url": str(request["url"]),
            "resolved_profile_name": str(request["resolved_profile_name"]),
            "drug_name": str(request["drug_name"]),
            "type": str(request["type"]),
            "match_level": "exact_allele",
            "match_context": match_context,
            "resolution_note": build_resolution_note(request, matched_variant),
            "match_source": legacy_match_source,
            "match_confidence": "high",
            "gene": str(request["gene"]),
            "variant_name": matched_variant["variant_name"],
            "chromosome": maf_row.get("Chromosome", ""),
            "start_position": maf_row.get("Start_position", ""),
            "end_position": maf_row.get("End_position", ""),
            "reference_allele": maf_row.get("Reference_Allele", ""),
            "tumor_seq_allele2": maf_row.get("Tumor_Seq_Allele2", ""),
            "variant_classification": maf_row.get("Variant_Classification", ""),
            "protein_change": matched_variant["protein_change"],
            "cdna_change": matched_variant.get("cdna_change", ""),
            "sample_id": maf_row.get("Tumor_Sample_Barcode", ""),
            "tumor_type": maf_row.get("Project_Code", ""),
            "project_code": maf_row.get("Project_Code", ""),
            "civic_submitted_variant": str(request["resolved_profile_name"]),
            "civic_gene": str(request["gene"]),
            "civic_match_rule": "exact_allele",
            "civic_match_status": "matched",
            "helper_gene": maf_row.get("Hugo_Symbol", ""),
            "helper_project_code": maf_row.get("Project_Code", ""),
            "helper_variant_key": helper_variant_key(maf_row),
            "helper_is_coding_or_disruptive": "yes",
        }
    )
    for field in REQUIRED_COLUMNS:
        row.setdefault(field, "")
    return row


def write_log(log_path: Path, lines: list[str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()

    existing_rows, existing_fields = load_existing_rows(args.output)
    rows = [backfill_row(row) for row in existing_rows]
    before_count = len(rows)
    existing_keys = {dedupe_key(row) for row in rows}

    maf_fields, hits = load_maf_hits(args.maf)

    fieldnames = list(existing_fields)
    for field in REQUIRED_COLUMNS:
        if field not in fieldnames:
            fieldnames.append(field)
    for field in maf_fields:
        if field not in fieldnames:
            fieldnames.append(field)
    for field in (
        "civic_submitted_variant",
        "civic_gene",
        "civic_match_rule",
        "civic_match_status",
        "helper_gene",
        "helper_project_code",
        "helper_variant_key",
        "helper_is_coding_or_disruptive",
    ):
        if field not in fieldnames:
            fieldnames.append(field)

    stats = Counter()
    matched_breast = 0
    matched_pan = 0
    skipped_existing = 0
    not_found_profiles: list[str] = []
    accepted_change_lines: list[str] = []

    for request in REQUESTS:
        if request.get("accepted_protein_changes"):
            accepted_change_lines.append(
                f"{request['id']} accepted protein changes: "
                + ", ".join(request["accepted_protein_changes"])
            )

        matched_variant, maf_row, match_context = select_best_exact_match(request["variants"], hits)
        if maf_row is None or matched_variant is None or match_context is None:
            stats["not_found"] += 1
            reason = str(request.get("not_found_reason", "no exact genomic MAF match"))
            not_found_profiles.append(
                f"{request['id']} | {request['url']} | {reason} | {request['resolved_profile_name']}"
            )
            continue

        output_row = build_output_row(request, matched_variant, maf_row, match_context)
        key = dedupe_key(output_row)
        if key not in existing_keys:
            rows.append(output_row)
            existing_keys.add(key)
        else:
            skipped_existing += 1

        if match_context == "breast_cancer":
            matched_breast += 1
        else:
            matched_pan += 1

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    summary = [
        f"total unique profiles searched: {len(REQUESTS)}",
        f"matched in breast cancer: {matched_breast}",
        f"matched in pan-cancer fallback: {matched_pan}",
        f"skipped as duplicate within batch: {BATCH_DUPLICATE_SKIPS}",
        f"skipped as already present: {skipped_existing}",
        f"not found: {stats['not_found']}",
        f"output rows before: {before_count}",
        f"output rows after: {len(rows)}",
        "",
        "MAF matching note: the local ICGC MAF does not contain protein-change or cDNA-change columns, so matches were saved at exact_allele level after CIViC profile resolution.",
        "Therapy-specific duplicate rows are allowed when drug_name or type differs.",
        "",
        "accepted protein changes after API resolution:",
    ]
    summary.extend(accepted_change_lines or ["none"])
    summary.append("")
    summary.append("not found details:")
    summary.extend(not_found_profiles or ["none"])

    write_log(args.not_found_log, summary)

    print(f"Total unique profiles searched: {len(REQUESTS)}")
    print(f"Matched in breast cancer: {matched_breast}")
    print(f"Matched in pan-cancer fallback: {matched_pan}")
    print(f"Skipped as duplicate within batch: {BATCH_DUPLICATE_SKIPS}")
    print(f"Skipped as already present: {skipped_existing}")
    print(f"Not found: {stats['not_found']}")
    print(f"Output rows before: {before_count}")
    print(f"Output rows after: {len(rows)}")
    print(f"Not found log: {args.not_found_log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
