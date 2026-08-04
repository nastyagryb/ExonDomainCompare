#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List


COLS = [
    "species", "isoform", "transcript_id", "protein_id",
    "original_boundary_precision", "original_reason_if_unknown", "rescue_attempted",
    "phase_found_in_original_cds_features", "phase_found_in_source_gff3",
    "phase_inferred_from_cds_reconstruction", "phase_inferred_from_translation",
    "reconstructed_cds_nt_length", "reconstructed_protein_length", "selected_protein_length",
    "translation_check_status", "rescued_left_boundary_precision",
    "rescued_right_boundary_precision", "rescued_boundary_precision",
    "rescue_status", "rescue_warning",
]

CONSISTENT_RECON = {
    "cds_reconstruction_matches_protein",
    "cds_reconstruction_matches_with_terminal_stop_offset",
}


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, rows, fields) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def _int(v, d=None):
    try:
        return int(float(str(v).strip()))
    except Exception:
        return d


def norm_tx(tx: str) -> str:
    t = (tx or "").strip()
    for pref in ("rna-", "transcript:", "transcript-"):
        if t.lower().startswith(pref):
            t = t[len(pref):]
    return t


def combine(left_split: bool, right_split: bool) -> str:
    if left_split and right_split:
        return "codon_split_both_sides"
    if left_split or right_split:
        return "codon_split_one_side"
    return "codon_boundary_exact"


def main() -> int:
    ap = argparse.ArgumentParser(description="Rescue unknown-phase IIIb/IIIc cassettes.")
    ap.add_argument("--cds_audit", type=Path, required=True)
    ap.add_argument("--cds_features", type=Path, required=True)
    ap.add_argument("--cassette_map", type=Path, required=True)
    ap.add_argument("--reconstruction_audit", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    audit = read_tsv(args.cds_audit)
    cmap = {(r["species"], r["isoform"]): r for r in read_tsv(args.cassette_map)}
    recon = {(r["species"], r["isoform"]): r for r in read_tsv(args.reconstruction_audit)}

    cds_by_tx: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for c in read_tsv(args.cds_features):
        cds_by_tx[norm_tx(c.get("transcript_id_source"))].append(c)
    for tx in cds_by_tx:
        cds_by_tx[tx].sort(key=lambda b: _int(b.get("cds_rank"), 0) or 0)

    out = []
    for a in audit:
        sp, iso = a["species"], a["isoform"]
        tx = a["transcript_id"]
        pid = a["protein_id"]
        orig = a.get("cds_boundary_precision_refined", "")
        reason = a.get("reason_if_unknown", "not_unknown")
        m = cmap.get((sp, iso), {})
        rc = recon.get((sp, iso), {})
        rank = _int(m.get("matched_cds_rank"))
        blocks = cds_by_tx.get(norm_tx(tx), [])

        # does the cassette block carry an explicit phase in the parsed CDS model?
        cass_blk = next((b for b in blocks if _int(b.get("cds_rank")) == rank), None)
        phase_in_features = bool((cass_blk or {}).get("phase", "").strip())

        needs_rescue = ("unknown" in orig.lower()
                        or reason in ("phase_not_propagated_from_source",
                                      "missing_gff3_phase", "nucleotide_sequence_unavailable"))

        rec_status = rc.get("reconstruction_status", "")
        consistent = rec_status in CONSISTENT_RECON

        row = {
            "species": sp, "isoform": iso, "transcript_id": tx, "protein_id": pid,
            "original_boundary_precision": orig, "original_reason_if_unknown": reason,
            "rescue_attempted": "true" if needs_rescue else "false",
            "phase_found_in_original_cds_features": str(phase_in_features).lower(),
            "phase_found_in_source_gff3": "false",  # not re-fetched here (Part C handles NCBI fetch)
            "phase_inferred_from_cds_reconstruction": "false",
            "phase_inferred_from_translation": "false",
            "reconstructed_cds_nt_length": rc.get("reconstructed_cds_nt_length", ""),
            "reconstructed_protein_length": rc.get("reconstructed_protein_length", ""),
            "selected_protein_length": rc.get("selected_protein_length", ""),
            "translation_check_status": rec_status,
            "rescued_left_boundary_precision": a.get("cds_left_boundary_precision", ""),
            "rescued_right_boundary_precision": a.get("cds_right_boundary_precision", ""),
            "rescued_boundary_precision": orig,
            "rescue_status": "", "rescue_warning": "",
        }

        if not needs_rescue:
            row["rescue_status"] = "not_needed_already_resolved"
            out.append(row)
            continue

        # 1) explicit source phase present -> already resolvable from original phase
        if phase_in_features and "unknown" not in orig.lower():
            row["rescue_status"] = "rescued_from_original_phase"
            out.append(row)
            continue

        # 2) infer from cumulative CDS reconstruction (only if length-consistent)
        if rank is not None and blocks and consistent:
            before = sum(_int(b.get("cds_length_bp"), 0) or 0
                         for b in blocks if (_int(b.get("cds_rank"), 0) or 0) < rank)
            blk_len = _int((cass_blk or {}).get("cds_length_bp"), 0) or 0
            through = before + blk_len
            left_split = (before % 3) != 0
            right_split = (through % 3) != 0
            rl = "codon_split_codon" if left_split else "codon_boundary_exact"
            rr = "codon_split_codon" if right_split else "codon_boundary_exact"
            row["phase_inferred_from_cds_reconstruction"] = "true"
            row["rescued_left_boundary_precision"] = rl
            row["rescued_right_boundary_precision"] = rr
            row["rescued_boundary_precision"] = combine(left_split, right_split)
            row["rescue_status"] = "rescued_from_cds_reconstruction"
            out.append(row)
            continue

        # 3) not rescuable -> keep uncertain, label clearly (coordinate still resolved)
        if not blocks:
            row["rescued_boundary_precision"] = "nucleotide_sequence_unavailable"
            row["rescue_status"] = "not_rescuable_sequence_unavailable"
            row["rescue_warning"] = "transcript absent from local CDS model"
        elif not consistent:
            row["rescued_boundary_precision"] = "phase_not_available_but_coordinate_resolved"
            row["rescue_status"] = "not_rescuable_translation_mismatch"
            row["rescue_warning"] = f"reconstruction not length-consistent ({rec_status})"
        else:
            row["rescued_boundary_precision"] = "phase_not_available_but_coordinate_resolved"
            row["rescue_status"] = "not_rescuable_phase_absent"
        out.append(row)

    write_tsv(args.outdir / "cds_phase_rescue_audit.tsv", out, COLS)
    from collections import Counter
    print(f"[OK] cds_phase_rescue_audit.tsv rows={len(out)}")
    print(f"     rescue_status={dict(Counter(r['rescue_status'] for r in out))}")
    print(f"     rescued_boundary_precision={dict(Counter(r['rescued_boundary_precision'] for r in out))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
