#!/usr/bin/env python3
"""
Build the FGFR2 orthology evidence table.

Build a dedicated FGFR2 orthology evidence table that integrates:
  * gene-symbol evidence (Step 2 genes.tsv),
  * paralog-panel evidence (multi-vertebrate FGFR1/2/3/4 screen),
  * protein QC (Step 5b validation, used as QC only),
  * sequence-calibrated IIIb/IIIc direction status (Step 4),
  * joint detection of both isoforms (Step 10 pair QC).

This table does NOT change IIIb/IIIc labels. It only summarises orthology
confidence and a recommended use per transcript/protein.

Outputs:
  fgfr2_orthology_evidence.tsv          (one row per species/isoform/transcript/protein)
  fgfr2_orthology_species_summary.tsv   (one row per species)
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional

SCRIPT_VERSION = "1.0"

EVID_COLS = [
    "species", "isoform", "transcript_id", "protein_id", "fgfr2_gene_symbol_evidence",
    "reciprocal_or_best_hit_evidence", "paralog_panel_status", "protein_qc_status",
    "sequence_calibrated_direction_status", "both_isoforms_detected", "orthology_confidence",
    "orthology_status", "orthology_warning", "recommended_use",
]

_PROTEIN_OK = {"protein_validated_expected_isoform",
               "protein_supports_expected_isoform_below_strict_threshold"}
_PROTEIN_TRANSCRIPT_ONLY = {"transcript_level_evidence_only_protein_unavailable"}
_PROTEIN_REVIEW = {"protein_conflicts_expected_isoform", "protein_ambiguous_inconclusive",
                   "unresolved_segment_extraction", "systematic_isoform_inversion_review"}


def read_tsv(path: Optional[Path]) -> List[Dict[str, str]]:
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return []
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def norm(s: str) -> str:
    return str(s or "").strip().lower()


def strip_ver(acc: str) -> str:
    return str(acc or "").split(".")[0]


def main() -> int:
    ap = argparse.ArgumentParser(description="Build the FGFR2 orthology evidence table.")
    ap.add_argument("--genes", type=Path, required=True, help="Step 2 genes.tsv")
    ap.add_argument("--paralog_detailed", type=Path, required=True, help="Multi-vertebrate paralog screen detailed")
    ap.add_argument("--protein_validation_summary", type=Path, required=True, help="Step 5b final selected validation summary")
    ap.add_argument("--isoform_evidence", type=Path, required=True, help="Step 4 isoform evidence")
    ap.add_argument("--pair_qc", type=Path, required=True, help="Step 10 pair-level QC summary")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    genes = {norm(r.get("species_canonical")): r for r in read_tsv(args.genes)}

    # Protein QC by (species, isoform, protein) with fallbacks.
    prot_qc: Dict[tuple, str] = {}
    prot_qc_sp_iso: Dict[tuple, List[str]] = defaultdict(list)
    for r in read_tsv(args.protein_validation_summary):
        sp = norm(r.get("species")); iso = r.get("expected_isoform_final") or r.get("isoform") or ""
        pid = strip_ver(r.get("protein_id"))
        vs = str(r.get("validation_status", "")).strip()
        prot_qc[(sp, iso, pid)] = vs
        prot_qc_sp_iso[(sp, iso)].append(vs)

    # Direction status per species (Step 4).
    dir_by_sp: Dict[str, set] = defaultdict(set)
    for r in read_tsv(args.isoform_evidence):
        sp = norm(r.get("species_canonical") or r.get("species"))
        st = str(r.get("direction_validation_status", "")).strip()
        if sp and st:
            dir_by_sp[sp].add(st)

    # Both isoforms per species (Step 10).
    both_by_sp: Dict[str, bool] = {}
    for r in read_tsv(args.pair_qc):
        sp = norm(r.get("species_canonical") or r.get("species"))
        both_by_sp[sp] = norm(r.get("has_both_isoforms")) in ("true", "1", "yes")

    rows: List[Dict[str, object]] = []
    for r in read_tsv(args.paralog_detailed):
        sp = norm(r.get("species")); iso = r.get("isoform", ""); pid = strip_ver(r.get("protein_id"))
        g = genes.get(sp, {})
        sym_found = g.get("gene_symbol_found", "")
        sym_evidence = ("gene_symbol_confirmed_FGFR2" if sym_found.upper() == "FGFR2"
                        else (f"gene_symbol_mismatch_{sym_found}" if sym_found else "gene_symbol_unavailable"))

        best_gene = r.get("best_paralog_gene", "")
        best_hit_ev = ("best_panel_hit_is_FGFR2" if best_gene == "FGFR2"
                       else (f"best_panel_hit_is_{best_gene}" if best_gene else "no_panel_hit"))
        paralog_status = r.get("paralog_status", "paralog_evidence_unavailable")

        pqc = prot_qc.get((sp, iso, pid))
        if pqc is None:
            cand = prot_qc_sp_iso.get((sp, iso), [])
            pqc = cand[0] if cand else "protein_qc_unavailable"

        dir_set = sorted(dir_by_sp.get(sp, set()))
        dir_status = dir_set[0] if len(dir_set) == 1 else (";".join(dir_set) if dir_set else "direction_unresolved_no_sequence")
        both = both_by_sp.get(sp, False)

        # Derive orthology status.
        warnings: List[str] = []
        sym_ok = sym_evidence == "gene_symbol_confirmed_FGFR2"
        if not sym_ok:
            warnings.append(sym_evidence)
        paralog_ok = paralog_status in ("fgfr2_high_confidence_multi_vertebrate", "fgfr2_supported_low_margin",
                                        "fgfr2_supported_human_only")
        paralog_ambig = paralog_status in ("ambiguous_fgfr_paralog_review", "non_fgfr2_best_hit_review")
        protein_ok = pqc in _PROTEIN_OK
        protein_transcript_only = pqc in _PROTEIN_TRANSCRIPT_ONLY
        protein_review = pqc in _PROTEIN_REVIEW
        if paralog_status not in ("fgfr2_high_confidence_multi_vertebrate",):
            warnings.append(f"paralog:{paralog_status}")
        if protein_review:
            warnings.append(f"protein_qc:{pqc}")
        if "unresolved" in dir_status or "ambiguous" in dir_status:
            warnings.append(f"direction:{dir_status}")
        if not both:
            warnings.append("both_isoforms_not_jointly_detected")

        if paralog_ambig:
            status = "fgfr2_ortholog_ambiguous_paralog_review"
            conf = "low"
        elif not paralog_ok and not protein_ok and not protein_transcript_only:
            status = "fgfr2_ortholog_unresolved"
            conf = "low"
        elif sym_ok and paralog_status == "fgfr2_high_confidence_multi_vertebrate" and protein_ok and both:
            status = "fgfr2_ortholog_high_confidence"
            conf = "high"
        elif protein_transcript_only and (sym_ok or paralog_ok):
            status = "fgfr2_ortholog_transcript_level_only"
            conf = "medium"
        elif (sym_ok or paralog_ok) and not protein_review:
            status = "fgfr2_ortholog_supported_with_warnings"
            conf = "medium"
        elif sym_ok or paralog_ok:
            status = "fgfr2_ortholog_supported_with_warnings"
            conf = "medium"
        else:
            status = "fgfr2_ortholog_unresolved"
            conf = "low"

        recommended = {
            "fgfr2_ortholog_high_confidence": "main_text_primary_claim",
            "fgfr2_ortholog_supported_with_warnings": "main_text_with_footnote",
            "fgfr2_ortholog_transcript_level_only": "supplementary_transcript_level",
            "fgfr2_ortholog_ambiguous_paralog_review": "exclude_from_primary_claims",
            "fgfr2_ortholog_unresolved": "exclude_from_primary_claims",
        }[status]

        rows.append({
            "species": sp, "isoform": iso, "transcript_id": r.get("transcript_id", ""), "protein_id": r.get("protein_id", ""),
            "fgfr2_gene_symbol_evidence": sym_evidence,
            "reciprocal_or_best_hit_evidence": best_hit_ev,
            "paralog_panel_status": paralog_status,
            "protein_qc_status": pqc,
            "sequence_calibrated_direction_status": dir_status,
            "both_isoforms_detected": "true" if both else "false",
            "orthology_confidence": conf,
            "orthology_status": status,
            "orthology_warning": ";".join(warnings) if warnings else "",
            "recommended_use": recommended,
        })

    rows.sort(key=lambda r: (str(r["species"]), str(r["isoform"])))
    args.outdir.mkdir(parents=True, exist_ok=True)
    with open(args.outdir / "fgfr2_orthology_evidence.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=EVID_COLS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in EVID_COLS})

    # Species summary (worst-wins).
    _RANK = {"fgfr2_ortholog_high_confidence": 5, "fgfr2_ortholog_supported_with_warnings": 4,
             "fgfr2_ortholog_transcript_level_only": 3, "fgfr2_ortholog_ambiguous_paralog_review": 2,
             "fgfr2_ortholog_unresolved": 1}
    by_sp: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for r in rows:
        by_sp[str(r["species"])].append(r)
    _REC = {"fgfr2_ortholog_high_confidence": "main_text_primary_claim",
            "fgfr2_ortholog_supported_with_warnings": "main_text_with_footnote",
            "fgfr2_ortholog_transcript_level_only": "supplementary_transcript_level",
            "fgfr2_ortholog_ambiguous_paralog_review": "exclude_from_primary_claims",
            "fgfr2_ortholog_unresolved": "exclude_from_primary_claims"}
    summ: List[Dict[str, object]] = []
    for sp, rs in sorted(by_sp.items()):
        statuses = [str(r["orthology_status"]) for r in rs]
        worst = min(statuses, key=lambda s: _RANK.get(s, 0))
        conf = next((str(r["orthology_confidence"]) for r in rs if r["orthology_status"] == worst), "")
        summ.append({
            "species": sp, "n_records": len(rs),
            "orthology_status_species": worst,
            "orthology_confidence_species": conf,
            "gene_symbol_evidence": rs[0]["fgfr2_gene_symbol_evidence"],
            "paralog_panel_status_set": ";".join(sorted({str(r["paralog_panel_status"]) for r in rs})),
            "both_isoforms_detected": rs[0]["both_isoforms_detected"],
            "recommended_use": _REC[worst],
        })
    with open(args.outdir / "fgfr2_orthology_species_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        cols = ["species", "n_records", "orthology_status_species", "orthology_confidence_species",
                "gene_symbol_evidence", "paralog_panel_status_set", "both_isoforms_detected", "recommended_use"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        w.writerows(summ)

    meta = {
        "script_version": SCRIPT_VERSION,
        "n_records": len(rows),
        "orthology_status_counts": dict(Counter(str(r["orthology_status"]) for r in rows)),
        "species_orthology_status_counts": dict(Counter(str(s["orthology_status_species"]) for s in summ)),
    }
    (args.outdir / "fgfr2_orthology_evidence_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"[OK] orthology evidence: {len(rows)} records; status={meta['orthology_status_counts']}")
    print(f"     species-level: {meta['species_orthology_status_counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
