#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


AUDIT_COLS = [
    "species", "isoform", "transcript_id", "protein_id", "strand",
    "cds_feature_id", "cds_genomic_start", "cds_genomic_end", "cds_order_in_transcript",
    "cds_phase_raw", "cds_phase_interpreted", "cds_left_phase", "cds_right_phase",
    "cds_left_boundary_precision", "cds_right_boundary_precision",
    "legacy_cds_boundary_precision", "cds_boundary_precision_refined",
    "cds_boundary_confidence", "cds_phase_source", "cds_phase_warning",
    "reason_if_unknown", "reason_if_split",
    "transcript_cds_reconstruction_status", "protein_translation_check_status",
]

# Explainability columns propagated into coordinate outputs.
PROPAGATE_COLS = [
    "reason_if_unknown", "reason_if_split",
    "transcript_cds_reconstruction_status", "protein_translation_check_status",
]


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, rows: List[Dict[str, object]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def _int(v, default=None):
    try:
        return int(float(str(v).strip()))
    except Exception:
        return default


def norm_tx(tx: str) -> str:
    t = (tx or "").strip()
    for pref in ("rna-", "transcript:", "transcript-"):
        if t.lower().startswith(pref):
            t = t[len(pref):]
    return t


def parse_protein_lengths(fasta: Optional[Path]) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    if not fasta or not fasta.exists():
        return lengths
    pid, seqlen = None, 0
    with open(fasta, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if pid:
                    lengths[pid] = max(lengths.get(pid, 0), seqlen)
                pid, seqlen = None, 0
                for tok in line[1:].strip().split("|"):
                    if tok.startswith("protein="):
                        pid = tok.split("=", 1)[1].strip()
            else:
                seqlen += len(line.strip())
        if pid:
            lengths[pid] = max(lengths.get(pid, 0), seqlen)
    return lengths


def reconstruct_transcript(cds_list: List[Dict[str, str]], protein_len: Optional[int]) -> Tuple[str, str]:
    if not cds_list:
        return "no_cds_features_for_transcript", "protein_sequence_unavailable"
    total_bp = sum(_int(c.get("cds_length_bp"), 0) or 0 for c in cds_list)
    if total_bp <= 0:
        return "cds_features_incomplete", "protein_sequence_unavailable"
    if total_bp % 3 != 0:
        recon = "cds_length_not_multiple_of_3"
    else:
        recon = "reconstructed_from_cds_coordinates"
    if protein_len is None or protein_len <= 0:
        return recon, "protein_sequence_unavailable"
    expected = total_bp // 3  # may include stop codon
    if abs(expected - protein_len) <= 1:
        trans = "cds_protein_length_consistent"
    else:
        trans = "cds_protein_length_mismatch_review"
    return recon, trans


def reason_unknown(refined: str, source_db: str, phase_raw: str, cds_id: str,
                   recon_status: str) -> str:
    if "unknown" not in (refined or "").lower():
        return "not_unknown"
    if not (cds_id or "").strip():
        return "cds_feature_unmatched"
    if recon_status == "no_cds_features_for_transcript":
        return "nucleotide_sequence_unavailable"
    if (source_db or "").strip().lower().startswith("ensembl"):
        return "phase_not_propagated_from_source"
    if not (phase_raw or "").strip():
        return "missing_gff3_phase"
    return "missing_gff3_phase"


def reason_split(left_prec: str, right_prec: str) -> str:
    sl = "split" in (left_prec or "").lower()
    sr = "split" in (right_prec or "").lower()
    if sl and sr:
        return "both_boundaries_split_codons"
    if sl:
        return "left_boundary_splits_codon"
    if sr:
        return "right_boundary_splits_codon"
    return "not_split"


def species_explainability(rows: List[Dict[str, object]]) -> str:
    unknown = [r for r in rows if r["reason_if_unknown"] != "not_unknown"]
    split = [r for r in rows if r["reason_if_split"] != "not_split"]
    if not unknown and not split:
        return "all_boundaries_known"
    parts = []
    if unknown:
        reasons = sorted({str(r["reason_if_unknown"]) for r in unknown})
        parts.append(f"unknown_x{len(unknown)}:{'/'.join(reasons)}")
    if split:
        parts.append(f"split_x{len(split)}")
    return "; ".join(parts)


def main() -> int:
    ap = argparse.ArgumentParser(description="Explainable CDS phase/boundary audit (Part A).")
    ap.add_argument("--coordinate_audit", type=Path, required=True)
    ap.add_argument("--cds_features", type=Path, required=True)
    ap.add_argument("--proteins", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--update_coordinate_audit", action="store_true")
    ap.add_argument("--update_exon_cds_mapping", type=Path, default=None)
    ap.add_argument("--update_pair_qc", type=Path, default=None)
    args = ap.parse_args()

    coord = read_tsv(args.coordinate_audit)
    prot_len = parse_protein_lengths(args.proteins)

    cds_by_tx: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for c in read_tsv(args.cds_features):
        cds_by_tx[norm_tx(c.get("transcript_id_source"))].append(c)
    for tx in cds_by_tx:
        cds_by_tx[tx].sort(key=lambda c: _int(c.get("cds_rank"), 0) or 0)

    audit: List[Dict[str, object]] = []
    for r in coord:
        tx = (r.get("transcript_id_source") or "").strip()
        pid = (r.get("protein_id") or "").strip()
        src = (r.get("resolver_source_db") or r.get("source_db") or "").strip()
        refined = r.get("cds_boundary_precision_refined", "")
        phase_raw = r.get("resolver_phase", "")
        cds_id = r.get("resolver_cds_id_source") or r.get("resolver_internal_cds_id") or ""
        left_p = r.get("cds_left_boundary_precision", "")
        right_p = r.get("cds_right_boundary_precision", "")

        recon, trans = reconstruct_transcript(cds_by_tx.get(norm_tx(tx), []), prot_len.get(pid))
        r_unknown = reason_unknown(refined, src, phase_raw, cds_id, recon)
        r_split = reason_split(left_p, right_p)
        phase_interp = (r.get("cds_phase_value") or "").strip() or "uninterpretable"

        audit.append({
            "species": r.get("species_canonical", ""),
            "isoform": r.get("inferred_isoform", ""),
            "transcript_id": tx, "protein_id": pid,
            "strand": r.get("resolver_strand", ""),
            "cds_feature_id": cds_id,
            "cds_genomic_start": r.get("resolver_start", ""),
            "cds_genomic_end": r.get("resolver_end", ""),
            "cds_order_in_transcript": r.get("resolver_cds_rank", ""),
            "cds_phase_raw": phase_raw,
            "cds_phase_interpreted": phase_interp,
            "cds_left_phase": r.get("cds_left_phase", ""),
            "cds_right_phase": r.get("cds_right_phase", ""),
            "cds_left_boundary_precision": left_p,
            "cds_right_boundary_precision": right_p,
            "legacy_cds_boundary_precision": r.get("legacy_cds_boundary_precision", ""),
            "cds_boundary_precision_refined": refined,
            "cds_boundary_confidence": r.get("cds_boundary_confidence", ""),
            "cds_phase_source": r.get("cds_phase_source", ""),
            "cds_phase_warning": r.get("cds_phase_warning", ""),
            "reason_if_unknown": r_unknown,
            "reason_if_split": r_split,
            "transcript_cds_reconstruction_status": recon,
            "protein_translation_check_status": trans,
        })

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir / "cds_phase_boundary_audit.tsv", audit, AUDIT_COLS)

    # legacy vs refined
    legacy = Counter(str(r.get("legacy_cds_boundary_precision") or "unknown") for r in audit)
    refined_c = Counter(str(r.get("cds_boundary_precision_refined") or "unknown") for r in audit)
    cats = sorted(set(legacy) | set(refined_c))
    crows = [{"precision_category": c, "legacy_count": legacy.get(c, 0),
              "refined_count": refined_c.get(c, 0),
              "delta_refined_minus_legacy": refined_c.get(c, 0) - legacy.get(c, 0)} for c in cats]
    crows.append({"precision_category": "TOTAL_rows", "legacy_count": sum(legacy.values()),
                  "refined_count": sum(refined_c.values()), "delta_refined_minus_legacy": 0})
    write_tsv(args.outdir / "cds_phase_boundary_legacy_vs_refined_counts.tsv", crows,
              ["precision_category", "legacy_count", "refined_count", "delta_refined_minus_legacy"])

    exp_rows: List[Dict[str, object]] = []

    def add_block(dim: str, counter: Counter):
        for k, v in sorted(counter.items(), key=lambda kv: (-kv[1], kv[0])):
            exp_rows.append({"dimension": dim, "category": k, "count": v,
                             "fraction": round(v / max(1, len(audit)), 3)})

    add_block("reason_if_unknown", Counter(str(r["reason_if_unknown"]) for r in audit))
    add_block("reason_if_split", Counter(str(r["reason_if_split"]) for r in audit))
    add_block("cds_boundary_precision_refined", refined_c)
    add_block("cds_boundary_confidence", Counter(str(r["cds_boundary_confidence"]) for r in audit))
    add_block("transcript_cds_reconstruction_status",
              Counter(str(r["transcript_cds_reconstruction_status"]) for r in audit))
    add_block("protein_translation_check_status",
              Counter(str(r["protein_translation_check_status"]) for r in audit))
    # source-db breakdown for unknown cases
    unknown_src = Counter(
        (next((c.get("source_db") for c in cds_by_tx.get(norm_tx(r["transcript_id"]), [])), "unknown") or "unknown")
        for r in audit if r["reason_if_unknown"] != "not_unknown")
    add_block("unknown_by_source_db", unknown_src)
    write_tsv(args.outdir / "cds_phase_boundary_explainability_summary.tsv", exp_rows,
              ["dimension", "category", "count", "fraction"])

    key = {(str(r["species"]).lower(), str(r["isoform"])): r for r in audit}

    def propagate_rowwise(path: Path):
        rows = read_tsv(path)
        sp_col = "species_canonical" if rows and "species_canonical" in rows[0] else "species"
        iso_col = "inferred_isoform" if rows and "inferred_isoform" in rows[0] else "isoform"
        for r in rows:
            a = key.get((str(r.get(sp_col, "")).lower(), str(r.get(iso_col, ""))))
            if a:
                for c in PROPAGATE_COLS:
                    r[c] = a[c]
        fields = list(rows[0].keys()) if rows else []
        for c in PROPAGATE_COLS:
            if c not in fields:
                fields.append(c)
        write_tsv(path, rows, fields)

    if args.update_coordinate_audit:
        propagate_rowwise(args.coordinate_audit)
    if args.update_exon_cds_mapping:
        propagate_rowwise(args.update_exon_cds_mapping)

    if args.update_pair_qc and args.update_pair_qc.exists():
        per_sp: Dict[str, List[Dict[str, object]]] = defaultdict(list)
        for r in audit:
            per_sp[str(r["species"]).lower()].append(r)
        pq = read_tsv(args.update_pair_qc)
        for r in pq:
            sp = str(r.get("species_canonical") or r.get("species") or "").lower()
            r["cds_boundary_explainability_summary"] = species_explainability(per_sp.get(sp, []))
        fields = list(pq[0].keys()) if pq else []
        if "cds_boundary_explainability_summary" not in fields:
            fields.append("cds_boundary_explainability_summary")
        write_tsv(args.update_pair_qc, pq, fields)

    print(f"[OK] cds_phase_boundary_audit.tsv rows={len(audit)}")
    print(f"     reason_if_unknown={dict(Counter(r['reason_if_unknown'] for r in audit))}")
    print(f"     reason_if_split={dict(Counter(r['reason_if_split'] for r in audit))}")
    print(f"     reconstruction={dict(Counter(r['transcript_cds_reconstruction_status'] for r in audit))}")
    print(f"     translation_check={dict(Counter(r['protein_translation_check_status'] for r in audit))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
