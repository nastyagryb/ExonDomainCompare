#!/usr/bin/env python3
"""
prepare_interpro_clean_fasta_v2.py

Prepare selected FGFR2 protein FASTA records for InterProScan/Web submission.

Rationale
---------
InterProScan operates on submitted amino-acid sequences. For a comparative,
role-aware FGFR2 workflow, protein FASTA records need short stable identifiers,
non-redundant submission sets, and a mapping table that can expand InterPro
results back to all selected transcript roles.

This script therefore creates:
  1. a clean all-record FASTA with short IDs,
  2. a non-redundant clean FASTA deduplicated by exact amino-acid sequence,
  3. optional split FASTA files for batch submission,
  4. mapping tables from clean/unique IDs back to original selected proteins,
  5. warnings, metadata and small HTML/Markdown reports.

The script intentionally does not infer biology from InterPro results. It only
prepares a reproducible, traceable input layer for downstream domain annotation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

SCRIPT_NAME = "prepare_interpro_clean_fasta_v2.py"
SCRIPT_VERSION = "2.0"

# IUPAC one-letter amino-acid codes plus common ambiguity/unknown symbols.
# U (selenocysteine), O (pyrrolysine), B/Z/J/X are kept but reported.
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYXBZJUO")
AMBIGUOUS_AA = set("XBZJUO")
STOP_SYMBOLS = set("*")


@dataclass
class FastaRecord:
    original_index: int
    original_header: str
    original_fasta_id: str
    sequence: str
    sequence_hash: str
    sequence_length: int
    clean_id: str = ""
    unique_id: str = ""
    is_unique_representative: str = "0"
    duplicate_group_size: int = 1
    species_canonical: str = ""
    source_db: str = ""
    selection_role: str = ""
    transcript_id: str = ""
    protein_id: str = ""
    isoform: str = ""


@dataclass
class WarningRow:
    severity: str
    warning_type: str
    record_id: str
    original_header: str
    message: str


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def write_tsv(path: Path, rows: List[Dict[str, object]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def parse_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    header: Optional[str] = None
    seq_parts: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def wrap_sequence(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def safe_token(value: str, fallback: str = "NA") -> str:
    value = str(value or "").strip()
    if not value:
        return fallback
    value = re.sub(r"[^A-Za-z0-9_.:=,+-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or fallback


def canonicalize_sequence(seq: str) -> str:
    return re.sub(r"\s+", "", seq).upper()


def sequence_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii", errors="ignore")).hexdigest()


def parse_header(header: str) -> Dict[str, str]:
    """Parse both legacy and v2 protein-export headers.

    Supported examples:
      fgfr2prot_000001|species=homo_sapiens|source=Ensembl|role=reference|transcript=ENST...|protein=ENSP...|isoform=IIIb
      homo_sapiens|Ensembl|reference|transcript=ENST...|protein=ENSP...|isoform=IIIb
    """
    parts = header.split("|")
    out = {
        "original_header": header,
        "original_fasta_id": parts[0].strip() if parts else "",
        "species_canonical": "",
        "source_db": "",
        "selection_role": "",
        "transcript_id": "",
        "protein_id": "",
        "isoform": "",
    }

    key_values: Dict[str, str] = {}
    positional: List[str] = []
    for part in parts[1:]:
        if "=" in part:
            key, value = part.split("=", 1)
            key_values[key.strip().lower()] = value.strip()
        else:
            positional.append(part.strip())

    if key_values:
        out["species_canonical"] = key_values.get("species", key_values.get("species_canonical", ""))
        out["source_db"] = key_values.get("source", key_values.get("source_db", ""))
        out["selection_role"] = key_values.get("role", key_values.get("selection_role", ""))
        out["transcript_id"] = key_values.get("transcript", key_values.get("transcript_id", ""))
        out["protein_id"] = key_values.get("protein", key_values.get("protein_id", ""))
        out["isoform"] = key_values.get("isoform", key_values.get("iii_isoform", ""))
    else:
        # Legacy header format: species|source|role|transcript=... can still have key-value fields later.
        out["species_canonical"] = parts[0].strip() if len(parts) > 0 else ""
        out["source_db"] = parts[1].strip() if len(parts) > 1 else ""
        out["selection_role"] = parts[2].strip() if len(parts) > 2 else ""

    # Legacy positional first three after ID if no key-values for species/source/role.
    if not out["species_canonical"] and len(positional) >= 1:
        out["species_canonical"] = positional[0]
    if not out["source_db"] and len(positional) >= 2:
        out["source_db"] = positional[1]
    if not out["selection_role"] and len(positional) >= 3:
        out["selection_role"] = positional[2]

    # Scan all parts for transcript/protein/isoform keys even in mixed legacy headers.
    for part in parts[1:]:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key = key.strip().lower()
        value = value.strip()
        if key == "transcript" and not out["transcript_id"]:
            out["transcript_id"] = value
        elif key == "protein" and not out["protein_id"]:
            out["protein_id"] = value
        elif key == "isoform" and not out["isoform"]:
            out["isoform"] = value
    return out


def load_export_report(path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    """Load optional protein_export_report.tsv by FASTA/output ID when available."""
    if not path:
        return {}
    rows = read_tsv(path)
    index: Dict[str, Dict[str, str]] = {}
    for row in rows:
        for key in ["fasta_id", "clean_id", "output_id", "protein_export_id"]:
            value = row.get(key, "")
            if value:
                index[value] = row
        original_header = row.get("fasta_header", "") or row.get("output_header", "")
        if original_header:
            index[original_header.split("|", 1)[0]] = row
    return index


def validate_sequence(record_id: str, header: str, seq: str, min_len: int, max_len: int) -> Tuple[List[WarningRow], bool]:
    warnings: List[WarningRow] = []
    ok = True
    if not seq:
        warnings.append(WarningRow("error", "empty_sequence", record_id, header, "Sequence is empty."))
        return warnings, False

    invalid = sorted(set(seq) - VALID_AA - STOP_SYMBOLS)
    if invalid:
        warnings.append(
            WarningRow("error", "invalid_amino_acid_symbols", record_id, header, f"Invalid symbols found: {''.join(invalid)}")
        )
        ok = False

    stops = seq.count("*")
    if stops:
        warnings.append(WarningRow("error", "stop_symbol_in_protein_sequence", record_id, header, f"Found {stops} '*' stop symbol(s)."))
        ok = False

    ambiguous = sorted(set(seq) & AMBIGUOUS_AA)
    if ambiguous:
        warnings.append(
            WarningRow("warning", "ambiguous_amino_acid_symbols", record_id, header, f"Ambiguous or rare symbols present: {''.join(ambiguous)}")
        )

    if len(seq) < min_len:
        warnings.append(WarningRow("warning", "protein_shorter_than_expected", record_id, header, f"Length {len(seq)} < {min_len}."))
    if len(seq) > max_len:
        warnings.append(WarningRow("warning", "protein_longer_than_expected", record_id, header, f"Length {len(seq)} > {max_len}."))
    return warnings, ok


def build_records(input_fasta: Path, prefix: str, min_len: int, max_len: int) -> Tuple[List[FastaRecord], List[WarningRow]]:
    records: List[FastaRecord] = []
    warnings: List[WarningRow] = []
    seen_headers: Counter[str] = Counter()

    for idx, (header, raw_seq) in enumerate(parse_fasta(input_fasta), start=1):
        seq = canonicalize_sequence(raw_seq)
        parsed = parse_header(header)
        original_id = parsed.get("original_fasta_id", f"record_{idx}")
        clean_id = f"{safe_token(prefix)}_{idx:04d}"

        seen_headers[header] += 1
        if seen_headers[header] > 1:
            warnings.append(WarningRow("warning", "duplicate_original_header", original_id, header, "Original FASTA header occurs more than once."))

        seq_warnings, is_valid = validate_sequence(original_id, header, seq, min_len, max_len)
        warnings.extend(seq_warnings)
        if not is_valid:
            # Keep invalid records in mapping? No. InterPro input should not contain invalid proteins.
            # The strict mode decides later whether invalid records abort.
            continue

        record = FastaRecord(
            original_index=idx,
            original_header=header,
            original_fasta_id=original_id,
            sequence=seq,
            sequence_hash=sequence_hash(seq),
            sequence_length=len(seq),
            clean_id=clean_id,
            species_canonical=parsed.get("species_canonical", ""),
            source_db=parsed.get("source_db", ""),
            selection_role=parsed.get("selection_role", ""),
            transcript_id=parsed.get("transcript_id", ""),
            protein_id=parsed.get("protein_id", ""),
            isoform=parsed.get("isoform", ""),
        )
        records.append(record)

    return records, warnings


def assign_unique_ids(records: List[FastaRecord], prefix: str) -> List[FastaRecord]:
    hash_to_unique: Dict[str, str] = {}
    hash_to_group: Dict[str, List[FastaRecord]] = defaultdict(list)
    unique_records: List[FastaRecord] = []

    for record in records:
        if record.sequence_hash not in hash_to_unique:
            uid = f"{safe_token(prefix)}_U{len(hash_to_unique) + 1:04d}"
            hash_to_unique[record.sequence_hash] = uid
            record.is_unique_representative = "1"
            unique_records.append(record)
        record.unique_id = hash_to_unique[record.sequence_hash]
        hash_to_group[record.sequence_hash].append(record)

    for group in hash_to_group.values():
        size = len(group)
        for record in group:
            record.duplicate_group_size = size
    return unique_records


def write_fasta(path: Path, rows: List[Tuple[str, str]], width: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in rows:
            handle.write(f">{header}\n{wrap_sequence(seq, width)}\n")


def write_reports(outdir: Path, prefix: str, records: List[FastaRecord], unique_records: List[FastaRecord], warnings: List[WarningRow],
                  args: argparse.Namespace, output_paths: Dict[str, str]) -> None:
    report_md = outdir / f"{prefix.lower()}_interpro_prepare_report.md"
    report_html = outdir / f"{prefix.lower()}_interpro_prepare_report.html"

    status = "PASS" if not any(w.severity == "error" for w in warnings) else "REVIEW"
    duplicate_count = len(records) - len(unique_records)
    warn_counter = Counter(w.warning_type for w in warnings)

    md_lines = [
        f"# InterPro FASTA preparation report",
        "",
        f"**Status:** {status}",
        f"**Input records accepted:** {len(records)}",
        f"**Unique protein sequences:** {len(unique_records)}",
        f"**Duplicate records collapsed for InterPro:** {duplicate_count}",
        f"**Split size:** {args.split_size}",
        "",
        "## Output files",
    ]
    for key, value in output_paths.items():
        md_lines.append(f"- **{key}:** `{Path(value).name}`")
    md_lines.extend(["", "## Warning counts"])
    if warn_counter:
        for warning_type, count in sorted(warn_counter.items()):
            md_lines.append(f"- {warning_type}: {count}")
    else:
        md_lines.append("- none")
    md_lines.extend([
        "",
        "## Interpretation",
        "The unique FASTA should be submitted to InterProScan. The mapping table is required to expand InterPro results back to all selected FGFR2 transcript roles.",
    ])
    report_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    body = "\n".join(f"<p>{html.escape(line)}</p>" if line else "" for line in md_lines)
    report_html.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>InterPro FASTA preparation report</title>"
        "<style>body{font-family:Arial,sans-serif;max-width:900px;margin:2rem auto;line-height:1.45}code{background:#f2f2f2;padding:2px 4px}</style>"
        f"</head><body>{body}</body></html>",
        encoding="utf-8",
    )


def write_task9_outputs(outdir: Path, prefix: str, records: List[FastaRecord], unique_records: List[FastaRecord],
                        warnings: List[WarningRow], unique_fasta: Path, split_paths: List[str],
                        args: argparse.Namespace) -> Dict[str, str]:
    """Write the pre-InterPro summary, run instructions and input manifest.

    These outputs make the FASTA preparation a self-contained, reproducible input
    layer for a *later* InterProScan run. No InterPro results are assumed here.
    """
    low = prefix.lower()
    summary_tsv = outdir / f"{low}_interpro_prepare_summary.tsv"
    instructions_md = outdir / "interproscan_run_instructions.md"
    manifest_tsv = outdir / "interproscan_input_manifest.tsv"

    invalid = sum(1 for w in warnings if w.severity == "error" and w.warning_type in
                  ("empty_sequence", "invalid_amino_acid_symbols", "stop_symbol_in_protein_sequence"))
    species = sorted({r.species_canonical for r in records if r.species_canonical})
    isoforms = sorted({r.isoform for r in records if r.isoform})
    roles = sorted({r.selection_role for r in records if r.selection_role})
    iso_counts = Counter(r.isoform or "unknown" for r in records)
    role_counts = Counter(r.selection_role or "unknown" for r in records)
    species_iso = {sp: sorted({r.isoform for r in records if r.species_canonical == sp and r.isoform}) for sp in species}
    both_iso_species = sum(1 for sp in species if set(species_iso[sp]) >= {"IIIb", "IIIc"})
    lengths = [r.sequence_length for r in unique_records] or [0]

    summary_rows: List[Dict[str, object]] = [
        {"metric": "total_selected_proteins", "value": len(records)},
        {"metric": "unique_sequences", "value": len(unique_records)},
        {"metric": "duplicates_collapsed", "value": len(records) - len(unique_records)},
        {"metric": "invalid_sequences_rejected", "value": invalid},
        {"metric": "warning_rows", "value": len(warnings)},
        {"metric": "error_rows", "value": sum(1 for w in warnings if w.severity == "error")},
        {"metric": "species_covered", "value": len(species)},
        {"metric": "species_with_both_isoforms", "value": both_iso_species},
        {"metric": "isoforms_covered", "value": ";".join(isoforms) if isoforms else "none"},
        {"metric": "roles_covered", "value": ";".join(roles) if roles else "none"},
        {"metric": "records_isoform_IIIb", "value": iso_counts.get("IIIb", 0)},
        {"metric": "records_isoform_IIIc", "value": iso_counts.get("IIIc", 0)},
        {"metric": "records_role_reference", "value": role_counts.get("reference", 0)},
        {"metric": "unique_seq_len_min", "value": min(lengths)},
        {"metric": "unique_seq_len_max", "value": max(lengths)},
        {"metric": "unique_seq_len_mean", "value": round(sum(lengths) / len(lengths), 1)},
        {"metric": "interpro_input_fasta", "value": unique_fasta.name},
        {"metric": "interpro_status", "value": "interpro_ready_input_prepared" if (len(unique_records) > 0 and invalid == 0)
            else ("interpro_input_missing" if len(unique_records) == 0 else "interpro_pending")},
    ]
    write_tsv(summary_tsv, summary_rows, ["metric", "value"])

    # Per-file input manifest with checksums so a later run is reproducible.
    def sha256_file(p: Path) -> str:
        h = hashlib.sha256()
        h.update(p.read_bytes())
        return h.hexdigest()

    manifest_rows: List[Dict[str, object]] = []
    manifest_rows.append({
        "file": unique_fasta.name, "role": "interproscan_primary_input",
        "n_sequences": len(unique_records), "sha256": sha256_file(unique_fasta) if unique_fasta.exists() else "",
        "description": "Non-redundant amino-acid FASTA to submit to InterProScan (one record per unique sequence).",
    })
    for sp in split_paths:
        p = Path(sp)
        n = sum(1 for _ in parse_fasta(p)) if p.exists() else 0
        manifest_rows.append({
            "file": p.name, "role": "interproscan_split_input", "n_sequences": n,
            "sha256": sha256_file(p) if p.exists() else "",
            "description": f"Split batch ({args.split_size} sequences max) of the unique FASTA for size-limited submission.",
        })
    write_tsv(manifest_tsv, manifest_rows, ["file", "role", "n_sequences", "sha256", "description"])

    instructions = f"""# Running InterProScan on the prepared FGFR2 input

This directory contains an **InterProScan-ready** input set. InterProScan has **not**
been executed yet; no domain annotations exist in this thesis output so far.

## Input

- Primary FASTA: `{unique_fasta.name}` ({len(unique_records)} unique protein sequences)
- ID mapping (reconstruct results back to species/isoform/role/transcript/protein): `{low}_interpro_id_mapping.tsv`
- Unique-sequence mapping (duplicate groups): `{low}_interpro_unique_mapping.tsv`
- Preparation summary: `{low}_interpro_prepare_summary.tsv`
- Warnings: `{low}_interpro_prepare_warnings.tsv`

The sequence IDs are short and stable (e.g. `{prefix}_U0001`) and contain no spaces
or unstable characters, so they are safe for InterProScan output parsing.

## Run command (local InterProScan)

```bash
interproscan.sh \\
  -i {unique_fasta.name} \\
  -f TSV,GFF3,JSON \\
  -appl Pfam,SMART,PROSITE,PRINTS,CDD \\
  -goterms -pa -iprlookup \\
  -cpu 4 \\
  -o fgfr2_interproscan_results
```

For the InterProScan web service, submit `{unique_fasta.name}` directly.

## Reconstructing results

After InterProScan finishes, join its result IDs (`{prefix}_U####`) back to the
biological metadata using `{low}_interpro_id_mapping.tsv` (column `unique_id`) and,
for full duplicate expansion, `{low}_interpro_unique_mapping.tsv`. Each unique ID maps
to one or more (species, isoform, role, transcript_id, protein_id, original_header,
sequence_hash) tuples.

## Status

`interpro_status = interpro_ready_input_prepared` — domain annotation is **pending**.
"""
    instructions_md.write_text(instructions, encoding="utf-8")
    return {
        "prepare_summary_tsv": str(summary_tsv),
        "interproscan_run_instructions_md": str(instructions_md),
        "interproscan_input_manifest_tsv": str(manifest_tsv),
    }


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Prepare clean, traceable FGFR2 protein FASTA files for InterProScan.")
    ap.add_argument("--input", required=True, help="Protein FASTA exported from selected FGFR2 proteins.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--prefix", default="FGFR2", help="Prefix for clean IDs.")
    ap.add_argument("--split_size", type=int, default=50, help="Number of unique sequences per split FASTA file.")
    ap.add_argument("--min_protein_len", type=int, default=500, help="Warning threshold for short FGFR2 proteins.")
    ap.add_argument("--max_protein_len", type=int, default=1200, help="Warning threshold for long FGFR2 proteins.")
    ap.add_argument("--wrap", type=int, default=60, help="FASTA line width.")
    ap.add_argument("--protein_export_report", default="", help="Optional protein_export_report.tsv for additional metadata checks.")
    ap.add_argument("--strict", action="store_true", help="Abort on invalid sequences or missing accepted records.")
    args = ap.parse_args(argv)

    input_fasta = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = safe_token(args.prefix, "FGFR2")

    if args.split_size < 1:
        raise ValueError("--split_size must be >= 1")
    if args.min_protein_len < 1 or args.max_protein_len < args.min_protein_len:
        raise ValueError("Invalid protein length thresholds")

    records, warnings = build_records(input_fasta, prefix, args.min_protein_len, args.max_protein_len)
    if not records:
        warnings.append(WarningRow("error", "no_valid_records", "", "", "No valid FASTA records were accepted."))

    unique_records = assign_unique_ids(records, prefix)

    # Optional export report is loaded to verify that original FASTA IDs are traceable.
    export_report_index = load_export_report(Path(args.protein_export_report) if args.protein_export_report else None)
    if args.protein_export_report:
        for record in records:
            if record.original_fasta_id not in export_report_index:
                warnings.append(
                    WarningRow("warning", "original_id_not_found_in_export_report", record.original_fasta_id, record.original_header,
                               "Original FASTA ID was not found in the provided protein export report.")
                )

    has_errors = any(w.severity == "error" for w in warnings)
    if args.strict and has_errors:
        # Still write warnings and metadata before aborting.
        write_tsv(outdir / f"{prefix.lower()}_interpro_prepare_warnings.tsv", [asdict(w) for w in warnings],
                  ["severity", "warning_type", "record_id", "original_header", "message"])
        raise SystemExit("Strict mode: invalid FASTA records or fatal preparation errors detected.")

    full_fasta = outdir / f"{prefix.lower()}_interpro_clean_all.fasta"
    unique_fasta = outdir / f"{prefix.lower()}_interpro_clean_unique.fasta"
    mapping_tsv = outdir / f"{prefix.lower()}_interpro_id_mapping.tsv"
    unique_mapping_tsv = outdir / f"{prefix.lower()}_interpro_unique_mapping.tsv"
    warnings_tsv = outdir / f"{prefix.lower()}_interpro_prepare_warnings.tsv"
    metadata_json = outdir / "run_metadata.json"

    write_fasta(full_fasta, [(r.clean_id, r.sequence) for r in records], args.wrap)
    write_fasta(unique_fasta, [(r.unique_id, r.sequence) for r in unique_records], args.wrap)

    split_paths: List[str] = []
    for start in range(0, len(unique_records), args.split_size):
        chunk = unique_records[start : start + args.split_size]
        part_no = start // args.split_size + 1
        part_path = outdir / f"{prefix.lower()}_interpro_unique_part{part_no:02d}.fasta"
        write_fasta(part_path, [(r.unique_id, r.sequence) for r in chunk], args.wrap)
        split_paths.append(str(part_path))

    mapping_fields = [
        "clean_id", "unique_id", "is_unique_representative", "duplicate_group_size", "sequence_hash", "sequence_length",
        "original_index", "original_fasta_id", "species_canonical", "source_db", "selection_role", "transcript_id",
        "protein_id", "isoform", "original_header",
    ]
    write_tsv(mapping_tsv, [asdict(r) for r in records], mapping_fields)

    unique_rows: List[Dict[str, object]] = []
    by_uid: Dict[str, List[FastaRecord]] = defaultdict(list)
    for r in records:
        by_uid[r.unique_id].append(r)
    for uid, group in by_uid.items():
        rep = group[0]
        unique_rows.append({
            "unique_id": uid,
            "sequence_hash": rep.sequence_hash,
            "sequence_length": rep.sequence_length,
            "duplicate_group_size": len(group),
            "representative_clean_id": rep.clean_id,
            "all_clean_ids": ";".join(r.clean_id for r in group),
            "all_original_fasta_ids": ";".join(r.original_fasta_id for r in group),
            "all_species": ";".join(sorted(set(r.species_canonical for r in group if r.species_canonical))),
            "all_roles": ";".join(sorted(set(r.selection_role for r in group if r.selection_role))),
        })
    write_tsv(unique_mapping_tsv, unique_rows, [
        "unique_id", "sequence_hash", "sequence_length", "duplicate_group_size", "representative_clean_id",
        "all_clean_ids", "all_original_fasta_ids", "all_species", "all_roles",
    ])

    write_tsv(warnings_tsv, [asdict(w) for w in warnings], ["severity", "warning_type", "record_id", "original_header", "message"])

    output_paths = {
        "clean_all_fasta": str(full_fasta),
        "clean_unique_fasta": str(unique_fasta),
        "mapping_tsv": str(mapping_tsv),
        "unique_mapping_tsv": str(unique_mapping_tsv),
        "warnings_tsv": str(warnings_tsv),
    }
    task9_paths = write_task9_outputs(outdir, prefix, records, unique_records, warnings, unique_fasta, split_paths, args)
    output_paths.update(task9_paths)
    write_reports(outdir, prefix, records, unique_records, warnings, args, output_paths)

    metadata = {
        "script_name": SCRIPT_NAME,
        "script_version": SCRIPT_VERSION,
        "run_datetime_utc": datetime.now(timezone.utc).isoformat(),
        "input_fasta": str(input_fasta),
        "outdir": str(outdir),
        "prefix": prefix,
        "split_size": args.split_size,
        "min_protein_len": args.min_protein_len,
        "max_protein_len": args.max_protein_len,
        "wrap": args.wrap,
        "strict": bool(args.strict),
        "protein_export_report": args.protein_export_report,
        "input_records_accepted": len(records),
        "unique_sequences": len(unique_records),
        "duplicates_collapsed": len(records) - len(unique_records),
        "warning_rows": len(warnings),
        "error_rows": sum(1 for w in warnings if w.severity == "error"),
        "split_files": split_paths,
        "outputs": output_paths,
    }
    metadata_json.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"Input records accepted: {len(records)}")
    print(f"Unique protein sequences: {len(unique_records)}")
    print(f"Duplicates collapsed: {len(records) - len(unique_records)}")
    print(f"Warnings: {len(warnings)}")
    print(f"Clean unique FASTA for InterProScan: {unique_fasta}")
    print(f"Mapping table: {mapping_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
