#!/usr/bin/env python3
"""Build the generic, EXPLORATORY event-region evidence layer (never blocking).

This converts the raw pairwise isoform-difference table
(``event_candidate_regions.tsv``) into a stable, user-facing evidence schema
(``event_region_evidence.tsv``). The evidence layer is intentionally
gene/event-AGNOSTIC and never asserts a validated event region:

  * every row has ``evidence_status = exploratory`` (unless a later, curated
    collector such as UniProt appends ``curated_annotation`` rows),
  * confidence stays conservative (copied from the raw scan),
  * NO candidate is ever labelled a validated event and NO event views are
    activated by this file.

Isoform A and Isoform B are simply two protein isoforms of the SAME gene that
were compared; a "candidate region" is a sequence difference between them.

Input (from a core run's results/core_gene_analysis/):
  * event_candidate_regions.tsv  (optional; produced by scan_isoform_event_candidates.py)

Output (same directory):
  * event_region_evidence.tsv    (always written; header-only if no candidates)

Usage:
  python -m exondomaincompare.framework.build_event_region_evidence --run-id <run_id>
  python -m exondomaincompare.framework.build_event_region_evidence --core-dir <dir> [--analysis-id ID] [--gene-symbol SYM]
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from exondomaincompare.config import discover_repository_root, load_config
from exondomaincompare.runs.registry import RegistryError, resolve_run_record

PROJECT_ROOT = discover_repository_root(__file__)
RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)

# Stable user-facing evidence schema.
EVIDENCE_COLUMNS = [
    "analysis_id", "gene_symbol", "species_id", "event_candidate_id",
    "evidence_source", "evidence_status",
    "transcript_a", "transcript_b", "protein_a", "protein_b",
    "region_start_aa", "region_end_aa", "region_length_aa",
    "event_type_candidate", "exon_aligned",
    "confidence", "confidence_reason", "notes",
]

RAW_CANDIDATES_FILE = "event_candidate_regions.tsv"
EVIDENCE_FILE = "event_region_evidence.tsv"


def read_tsv(p: Path) -> List[Dict[str, str]]:
    if not Path(p).is_file():
        return []
    with open(p, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, columns: List[str], rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in columns})


def _int_or_blank(v: Any) -> Any:
    try:
        if v is None or v == "":
            return ""
        return int(float(v))
    except (TypeError, ValueError):
        return ""


def _is_exon_aligned(candidate_type: str, evidence: str) -> bool:
    return candidate_type.startswith("exon_aligned") or "exon_map" in (evidence or "")


def _base_event_type(candidate_type: str) -> str:
    """Strip the exon_aligned_ prefix so the biological difference type is explicit."""
    ct = (candidate_type or "").strip()
    return ct[len("exon_aligned_"):] if ct.startswith("exon_aligned_") else ct


def convert_candidates(rows: List[Dict[str, str]], analysis_id: str,
                       gene_symbol: str) -> List[Dict[str, Any]]:
    """Convert raw pairwise candidate rows into the exploratory evidence schema."""
    out: List[Dict[str, Any]] = []
    for i, r in enumerate(rows, start=1):
        candidate_type = (r.get("candidate_type", "") or "").strip()
        evidence = (r.get("evidence", "") or "").strip()
        exon_aligned = _is_exon_aligned(candidate_type, evidence)
        species_id = r.get("species_id", "")
        gene = r.get("gene_symbol", "") or gene_symbol
        start = _int_or_blank(r.get("candidate_start_aa"))
        end = _int_or_blank(r.get("candidate_end_aa"))
        length = _int_or_blank(r.get("candidate_length_aa"))
        # A stable-ish id (region + pair index) so downstream can reference a row.
        cand_id = f"{gene}:{species_id}:{start}-{end}:{i:03d}"
        # Distinguish protein-comparison and exon-aligned evidence.
        source = ("exon_aligned_isoform_difference" if exon_aligned
                  else "protein_isoform_difference")
        confidence = (r.get("confidence", "") or "low").strip()
        confidence_reason = (
            "Exon-aligned isoform difference (region matches an exon block); "
            "still exploratory, not a validated event."
            if exon_aligned else
            "Protein-level isoform sequence difference only; exploratory, "
            "not a validated event.")
        out.append({
            "analysis_id": analysis_id,
            "gene_symbol": gene,
            "species_id": species_id,
            "event_candidate_id": cand_id,
            "evidence_source": source,
            "evidence_status": "exploratory",
            "transcript_a": r.get("transcript_a", ""),
            "transcript_b": r.get("transcript_b", ""),
            "protein_a": r.get("protein_a", ""),
            "protein_b": r.get("protein_b", ""),
            "region_start_aa": start,
            "region_end_aa": end,
            "region_length_aa": length,
            "event_type_candidate": _base_event_type(candidate_type),
            "exon_aligned": "yes" if exon_aligned else "no",
            "confidence": confidence,
            "confidence_reason": confidence_reason,
            "notes": r.get("notes", "") or "Exploratory candidate; not a validated event region.",
        })
    return out


def build_evidence(core_dir: Path, analysis_id: str, gene_symbol: str) -> Dict[str, Any]:
    """Write event_region_evidence.tsv from event_candidate_regions.tsv.

    Always writes the file (header-only when there are no candidates), so a
    valid-but-empty evidence file is produced for genes with no differences.
    """
    raw = read_tsv(core_dir / RAW_CANDIDATES_FILE)
    rows = convert_candidates(raw, analysis_id, gene_symbol)
    out = core_dir / EVIDENCE_FILE
    write_tsv(out, EVIDENCE_COLUMNS, rows)
    return {"out": out, "n_evidence": len(rows), "n_raw": len(raw)}


def _resolve_run_meta(run_dir: Path) -> Dict[str, str]:
    rc = {}
    p = run_dir / "run_config.json"
    if p.is_file():
        try:
            rc = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            rc = {}
    return {
        "analysis_id": rc.get("analysis_id", "") or rc.get("case_study", ""),
        "gene_symbol": rc.get("gene_symbol", ""),
    }


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Convert raw isoform candidates into the exploratory event-region evidence layer.")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run-id", help="Run id under runs/.")
    g.add_argument("--core-dir", help="Explicit results/core_gene_analysis directory.")
    ap.add_argument("--analysis-id", default="", help="Analysis id (defaults from run_config).")
    ap.add_argument("--gene-symbol", default="", help="Gene symbol (defaults from run_config).")
    args = ap.parse_args(argv)

    if args.run_id:
        try:
            record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
        except RegistryError as exc:
            raise SystemExit(str(exc)) from exc
        if record is None:
            raise SystemExit(f"Run not found: {args.run_id}")
        if record.read_only:
            raise SystemExit(
                "Run is registered read-only; copy it before writing derived evidence.")
        run_dir = record.path
        core_dir = run_dir / "results" / "core_gene_analysis"
        meta = _resolve_run_meta(run_dir)
        analysis_id = args.analysis_id or meta["analysis_id"]
        gene_symbol = args.gene_symbol or meta["gene_symbol"]
    else:
        core_dir = Path(args.core_dir)
        if not core_dir.is_absolute():
            core_dir = PROJECT_ROOT / core_dir
        analysis_id = args.analysis_id
        gene_symbol = args.gene_symbol

    core_dir.mkdir(parents=True, exist_ok=True)
    res = build_evidence(core_dir, analysis_id, gene_symbol)
    print(f"OK  event evidence: {res['n_evidence']} exploratory row(s) "
          f"from {res['n_raw']} raw candidate(s).")
    print(f"    wrote: {res['out']}")
    print("    All rows are exploratory; no validated event region is asserted.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
