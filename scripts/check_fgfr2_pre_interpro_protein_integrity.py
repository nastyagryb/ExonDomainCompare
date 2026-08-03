#!/usr/bin/env python3
"""
Check pre-InterPro protein integrity.

Pre-InterPro protein integrity QC: validate that selected FGFR2 proteins are biologically
plausible InterProScan inputs (valid alphabet, no internal stop, plausible length, cassette
within bounds). Feeds the boundary robustness score, master and readiness reports.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402


COLS = ["species", "isoform", "protein_id", "transcript_id", "sequence_length",
        "valid_amino_acid_alphabet", "contains_internal_stop", "terminal_stop_present",
        "sequence_hash", "duplicate_group_id", "has_III_region_coordinates",
        "cassette_within_protein_bounds", "expected_length_range_status",
        "length_outlier_status", "interpro_ready", "protein_integrity_status",
        "protein_integrity_warning"]
SUMMARY_COLS = ["metric", "value"]

# FGFR2 full-length plausible range (vertebrate orthologs; partial models flagged)
LEN_MIN, LEN_MAX = 600, 900


def build_protein_lookup(faa: Path):
    by_key: Dict[Tuple[str, str, str], str] = {}
    by_pid: Dict[str, str] = {}
    for hid, seq in M.read_fasta(faa):
        meta = {t.split("=", 1)[0]: t.split("=", 1)[1] for t in hid.split("|") if "=" in t}
        sp = (meta.get("species") or "").lower()
        s = M.clean_alignment_seq(seq)
        pid = meta.get("protein") or ""
        if pid:
            by_pid.setdefault(pid, s)
            by_key.setdefault((sp, meta.get("isoform") or "", pid), s)
    return by_key, by_pid


def main() -> int:
    ap = argparse.ArgumentParser(description="Check pre-InterPro protein integrity.")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    out = dirs["protein_integrity"]

    coord = M.read_tsv(M.require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv"))
    cmap = {(r["species"].lower(), r["isoform"]): r
            for r in M.read_tsv(M.require(base, "fgfr2_cassette_cds_block_map.tsv"))}
    by_key, by_pid = build_protein_lookup(M.require(base, "selected_fgfr2_proteins.faa"))
    # The upstream exon-structure classifier's ``inferred_isoform`` is provisional: where
    # the marker-validated sequence contradicts it, the reconciliation stage corrects it.
    # Upstream FASTA headers and the cassette map are keyed by that provisional label, so
    # it stays the *join* key — but the label this QC table reports must be the same final
    # biological label the MSA manifests, truth table, boundary model and Gallery carry,
    # otherwise a downstream join by (species, isoform) silently swaps IIIb and IIIc.
    recon = M.load_label_reconciliation(base)

    # duplicate groups by sequence hash
    rows: List[Dict[str, object]] = []
    hash_group: Dict[str, str] = {}
    for c in coord:
        sp = (c.get("species_canonical") or "").lower()
        up_iso = c.get("inferred_isoform") or ""
        iso = M.final_label(recon, sp, up_iso)
        pid = c.get("protein_id") or ""
        tx = c.get("transcript_id_source") or ""
        seq = by_key.get((sp, up_iso, pid)) or by_pid.get(pid) or ""
        if not seq:
            rows.append({"species": sp, "isoform": iso, "protein_id": pid, "transcript_id": tx,
                         "sequence_length": 0, "valid_amino_acid_alphabet": "false",
                         "contains_internal_stop": "", "terminal_stop_present": "",
                         "sequence_hash": "", "duplicate_group_id": "", "has_III_region_coordinates": "",
                         "cassette_within_protein_bounds": "", "expected_length_range_status": "",
                         "length_outlier_status": "", "interpro_ready": "false",
                         "protein_integrity_status": "missing_sequence_fail",
                         "protein_integrity_warning": "no selected protein sequence"})
            continue
        body = seq[:-1] if seq.endswith("*") else seq
        terminal_stop = seq.endswith("*")
        internal_stop = "*" in body
        bad = M.invalid_residues(body)
        valid_alpha = (not bad)
        slen = len(body)
        h = M.sha256_text(body)
        dg = hash_group.setdefault(h, f"dgrp_{len(hash_group)+1:03d}")
        iii_s = M.to_int(c.get("III_region_start_aa"))
        iii_e = M.to_int(c.get("III_region_end_aa"))
        has_iii = "true" if (iii_s is not None and iii_e is not None) else "false"
        m = cmap.get((sp, up_iso), {})
        ce = M.to_int(m.get("matched_protein_end_aa")) or M.to_int(c.get("native_protein_end_aa"))
        within = "true" if (ce is not None and ce <= slen) else ("false" if ce is not None else "")
        if slen < LEN_MIN:
            lstat, lrange = "short_outlier_review", "below_expected_range"
        elif slen > LEN_MAX:
            lstat, lrange = "long_outlier_review", "above_expected_range"
        else:
            lstat, lrange = "within_expected_range", "within_expected_range"
        warns = []
        if internal_stop:
            warns.append("internal_stop_codon")
        if bad:
            warns.append(f"invalid_residues:{bad}")
        if within == "false":
            warns.append("cassette_beyond_protein_bounds")
        if lstat != "within_expected_range":
            warns.append(lstat)
        interpro_ready = "true" if (valid_alpha and not internal_stop) else "false"
        if not valid_alpha or internal_stop:
            status = "invalid_sequence_review"
        elif lstat != "within_expected_range" or within == "false":
            status = "protein_length_outlier_review"
        elif warns:
            status = "protein_integrity_pass_with_minor_warning"
        else:
            status = "protein_integrity_pass"
        rows.append({"species": sp, "isoform": iso, "protein_id": pid, "transcript_id": tx,
                     "sequence_length": slen, "valid_amino_acid_alphabet": str(valid_alpha).lower(),
                     "contains_internal_stop": str(internal_stop).lower(),
                     "terminal_stop_present": str(terminal_stop).lower(), "sequence_hash": h,
                     "duplicate_group_id": dg, "has_III_region_coordinates": has_iii,
                     "cassette_within_protein_bounds": within,
                     "expected_length_range_status": lrange, "length_outlier_status": lstat,
                     "interpro_ready": interpro_ready, "protein_integrity_status": status,
                     "protein_integrity_warning": ";".join(warns)})

    M.write_tsv(out / "fgfr2_pre_interpro_protein_integrity_qc.tsv", rows, COLS)
    status_counts = Counter(r["protein_integrity_status"] for r in rows)
    summary = [{"metric": "n_proteins", "value": len(rows)},
               {"metric": "n_interpro_ready", "value": sum(1 for r in rows if r["interpro_ready"] == "true")},
               {"metric": "n_duplicate_groups", "value": len(hash_group)}]
    for k, v in sorted(status_counts.items()):
        summary.append({"metric": f"status::{k}", "value": v})
    M.write_tsv(out / "fgfr2_pre_interpro_protein_integrity_summary.tsv", summary, SUMMARY_COLS)
    print(f"[OK] protein integrity: {len(rows)} proteins; status={dict(status_counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
