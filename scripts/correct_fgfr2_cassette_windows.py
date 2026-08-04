#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M

try:
    from Bio.Align import PairwiseAligner, substitution_matrices
    _ALIGNER = PairwiseAligner()
    _ALIGNER.mode = "local"
    _ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
    _ALIGNER.open_gap_score = -11
    _ALIGNER.extend_gap_score = -1
except Exception:
    _ALIGNER = None


AUDIT_COLS = ["species", "isoform", "protein_id", "transcript_id", "original_start_aa",
              "original_end_aa", "original_cassette_seq", "corrected_start_aa", "corrected_end_aa",
              "corrected_cassette_seq", "anchor_identity", "anchor_coverage", "length_delta",
              "correction_status", "reference_window_corrected", "correction_warning"]
MANIFEST_COLS = ["msa_input_id", "species", "display_species_name", "isoform", "upstream_label",
                 "legacy_label", "final_isoform_label", "validated_exon_type",
                 "label_consistency_status", "label_reconciliation_action", "protein_id",
                 "transcript_id", "sequence_type", "source_file", "native_start_aa",
                 "native_end_aa", "sequence_length", "sequence_hash", "recommended_use",
                 "final_display_class", "extraction_status", "extraction_warning",
                 "reference_window_corrected", "correction_status"]
CANON = {"IIIb": 51, "IIIc": 42}
CANON_TOL = 6


def build_protein_lookup(faa: Path):
    by_key, by_pid = {}, {}
    for hid, seq in M.read_fasta(faa):
        meta = {t.split("=", 1)[0]: t.split("=", 1)[1] for t in hid.split("|") if "=" in t}
        s = M.ungapped(M.clean_alignment_seq(seq))
        pid = meta.get("protein") or ""
        if pid:
            by_pid.setdefault(pid, s)
            by_key.setdefault(((meta.get("species") or "").lower(), meta.get("isoform") or "", pid), s)
    return by_key, by_pid


def consensus_anchor(cassettes: List[str], canon_len: int) -> str:
    clean = [c for c in cassettes if abs(len(c) - canon_len) <= CANON_TOL]
    if not clean:
        clean = cassettes
    # medoid = clean sequence closest to canonical length, then highest summed identity
    medoid = min(clean, key=lambda c: abs(len(c) - canon_len))
    if _ALIGNER is None:
        return medoid
    # build consensus by aligning each clean seq to medoid and voting per medoid position
    votes: List[Counter] = [Counter() for _ in medoid]
    for c in clean:
        try:
            aln = _ALIGNER.align(medoid, c)[0]
            idx = aln.indices  # 2 x L array of indices into (medoid, c); -1 for gap
            for mi, ci in zip(idx[0], idx[1]):
                if mi >= 0 and ci >= 0:
                    votes[mi][c[ci]] += 1
        except Exception:
            continue
    cons = "".join((v.most_common(1)[0][0] if v else medoid[i]) for i, v in enumerate(votes))
    return cons


def realign_extract(protein: str, anchor: str) -> Optional[Tuple[int, int, float, float]]:
    if _ALIGNER is None or not protein or not anchor:
        return None
    try:
        aln = _ALIGNER.align(anchor, protein)[0]
    except Exception:
        return None
    idx = aln.indices  # row0=anchor, row1=protein
    prot_positions = [p for p in idx[1] if p >= 0]
    if not prot_positions:
        return None
    start, end = min(prot_positions) + 1, max(prot_positions) + 1
    ident = sum(1 for ai, pi in zip(idx[0], idx[1])
                if ai >= 0 and pi >= 0 and anchor[ai] == protein[pi])
    aligned_cols = sum(1 for ai, pi in zip(idx[0], idx[1]) if ai >= 0 and pi >= 0)
    identity = ident / aligned_cols if aligned_cols else 0.0
    coverage = aligned_cols / len(anchor) if anchor else 0.0
    return start, end, round(identity, 4), round(coverage, 4)


def main() -> int:
    ap = argparse.ArgumentParser(description="Anchor-correct cassette windows (Part A pre-step).")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--min_identity", type=float, default=0.5)
    ap.add_argument("--min_coverage", type=float, default=0.6)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    inp, maps = dirs["inputs"], dirs["maps"]

    if _ALIGNER is None:
        print("[WARN] Biopython PairwiseAligner unavailable; cassette window correction skipped.",
              file=sys.stderr)
        return 0

    coord = M.read_tsv(M.require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv"))
    master = {r["species"].lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}
    by_key, by_pid = build_protein_lookup(M.require(base, "selected_fgfr2_proteins.faa"))
    recon = M.load_label_reconciliation(base)

    # original extracted cassettes from current manifests
    orig_man = {iso: {r["msa_input_id"]: r for r in
                      M.read_tsv(inp / f"fgfr2_{iso}_cassette_msa_input_manifest.tsv")}
                for iso in ("IIIb", "IIIc")}
    orig_cass = {iso: dict(M.read_fasta(inp / f"fgfr2_{iso}_cassette_msa_input.faa"))
                 for iso in ("IIIb", "IIIc")}
    # Anchor to curated UniProt-anchored references (consistent window definition for BOTH
    # isoforms); fall back to within-isoform consensus only if curated refs are unavailable.
    curated = dict(M.read_fasta(inp / "curated_human_FGFR2_IIIb_IIIc_cassette_reference.faa"))
    anchors = {}
    for iso in ("IIIb", "IIIc"):
        cref = next((s for cid, s in curated.items() if f"|{iso}|" in cid), "")
        anchors[iso] = cref or consensus_anchor(list(orig_cass[iso].values()), CANON[iso])

    # backup originals once
    backup = inp / "original_uncorrected"
    backup.mkdir(exist_ok=True)
    for iso in ("IIIb", "IIIc", "IIIb_IIIc_combined"):
        for suff in ("msa_input.faa", "msa_input_manifest.tsv"):
            src = inp / (f"fgfr2_{iso}_cassette_{suff}" if iso != "IIIb_IIIc_combined"
                         else f"fgfr2_IIIb_IIIc_combined_cassette_{suff}")
            if src.exists() and not (backup / src.name).exists():
                shutil.copyfile(src, backup / src.name)

    audit: List[Dict[str, object]] = []
    new_seq: Dict[Tuple[str, str], str] = {}     # (sid, iso) -> corrected cassette
    new_coords: Dict[str, Tuple[int, int, str, str]] = {}  # sid -> (start,end,status,flag)

    for c in coord:
        sp = (c.get("species_canonical") or "").lower()
        up_iso = c.get("inferred_isoform") or ""
        iso = M.final_label(recon, sp, up_iso)  # FINAL biological label
        pid = c.get("protein_id") or ""
        tx = c.get("transcript_id_source") or ""
        protein = by_key.get((sp, up_iso, pid)) or by_pid.get(pid) or ""
        man = orig_man.get(iso, {})
        sid = f"{sp}|{iso}|{pid}|{M.recommended_use_token(master.get(sp, {}).get('recommended_use',''))}"
        orow = man.get(sid, {})
        o_start = M.to_int(orow.get("native_start_aa"))
        o_end = M.to_int(orow.get("native_end_aa"))
        o_seq = orig_cass.get(iso, {}).get(sid, "")
        res = realign_extract(protein, anchors[iso]) if protein else None
        if res is None:
            audit.append({"species": sp, "isoform": iso, "protein_id": pid, "transcript_id": tx,
                          "original_start_aa": o_start or "", "original_end_aa": o_end or "",
                          "original_cassette_seq": o_seq, "corrected_start_aa": "",
                          "corrected_end_aa": "", "corrected_cassette_seq": "", "anchor_identity": "",
                          "anchor_coverage": "", "length_delta": "",
                          "correction_status": "correction_failed_no_protein",
                          "reference_window_corrected": "false",
                          "correction_warning": "protein sequence unavailable for realignment"})
            continue
        cs, ce, ident, cov = res
        corr_seq = protein[cs - 1:ce]
        accept = (ident >= args.min_identity and cov >= args.min_coverage
                  and abs(len(corr_seq) - CANON[iso]) <= 3 * CANON_TOL)
        if not accept:
            status, flag, warn = ("reference_window_unresolved_divergent_exon", "false",
                                  f"selected protein does not match {iso} consensus cassette "
                                  f"(identity={ident}, coverage={cov}); possible upstream "
                                  f"isoform-selection/label inconsistency; kept original window")
            use_start, use_end, use_seq = (o_start or cs), (o_end or ce), (o_seq or corr_seq)
        elif o_start is not None and abs(cs - o_start) <= 2 and abs(len(corr_seq) - len(o_seq or "")) <= 2:
            status, flag, warn = "window_confirmed", "false", ""
            use_start, use_end, use_seq = cs, ce, corr_seq
        else:
            status, flag, warn = ("reference_window_corrected", "true",
                                  f"window re-centered to {iso} consensus cassette "
                                  f"(was {o_start}-{o_end}, now {cs}-{ce})")
            use_start, use_end, use_seq = cs, ce, corr_seq
        new_seq[(sid, iso)] = use_seq
        new_coords[sid] = (use_start, use_end, status, flag)
        audit.append({"species": sp, "isoform": iso, "protein_id": pid, "transcript_id": tx,
                      "original_start_aa": o_start or "", "original_end_aa": o_end or "",
                      "original_cassette_seq": o_seq, "corrected_start_aa": cs, "corrected_end_aa": ce,
                      "corrected_cassette_seq": corr_seq, "anchor_identity": ident,
                      "anchor_coverage": cov, "length_delta": len(corr_seq) - (len(o_seq) if o_seq else 0),
                      "correction_status": status, "reference_window_corrected": flag,
                      "correction_warning": warn})

    M.write_tsv(maps / "fgfr2_cassette_window_correction_audit.tsv", audit, AUDIT_COLS)

    # rewrite cassette inputs + manifests with corrected windows
    combined_items: List[Tuple[str, str]] = []
    combined_man: List[Dict[str, object]] = []
    for iso in ("IIIb", "IIIc"):
        items, man_rows = [], []
        for sid, orow in orig_man[iso].items():
            seq = new_seq.get((sid, iso), orig_cass[iso].get(sid, ""))
            coords = new_coords.get(sid)
            row = dict(orow)
            if coords and seq:
                row["native_start_aa"], row["native_end_aa"] = coords[0], coords[1]
                row["sequence_length"] = len(seq)
                row["sequence_hash"] = M.sha256_text(seq)
                row["correction_status"] = coords[2]
                row["reference_window_corrected"] = coords[3]
            else:
                row["correction_status"] = "no_correction_applied"
                row["reference_window_corrected"] = "false"
            man_rows.append(row)
            if seq:
                items.append((sid, seq))
                combined_items.append((sid, seq))
                cm = dict(row)
                cm["sequence_type"] = "combined_cassette"
                combined_man.append(cm)
        M.write_fasta(inp / f"fgfr2_{iso}_cassette_msa_input.faa", items)
        M.write_tsv(inp / f"fgfr2_{iso}_cassette_msa_input_manifest.tsv", man_rows, MANIFEST_COLS)
    M.write_fasta(inp / "fgfr2_IIIb_IIIc_combined_cassette_msa_input.faa", combined_items)
    M.write_tsv(inp / "fgfr2_IIIb_IIIc_combined_cassette_msa_input_manifest.tsv",
                combined_man, MANIFEST_COLS)

    from collections import Counter as C
    sc = C(a["correction_status"] for a in audit)
    print(f"[OK] cassette window correction: {dict(sc)}")
    print(f"     anchors: IIIb={anchors['IIIb'][:30]}... IIIc={anchors['IIIc'][:30]}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
