#!/usr/bin/env python3
"""
analyze_fgfr2_isoform_discriminating_residues.py  (MSA boundary-robustness sprint, Part 8)

Identify cassette alignment positions that are conserved within IIIb, conserved within
IIIc, and different between IIIb and IIIc, using the COMBINED IIIb+IIIc cassette MSA as a
shared coordinate frame. This SUPPORTS the sequence-calibrated IIIb/IIIc distinction; it
NEVER relabels any isoform. Both main-only and all-species versions are produced.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fgfr2_msa_common as M  # noqa: E402


COLS = ["alignment_col", "human_IIIb_aa_if_available", "human_IIIc_aa_if_available",
        "IIIb_major_aa", "IIIc_major_aa", "IIIb_major_aa_fraction", "IIIc_major_aa_fraction",
        "IIIb_gap_fraction", "IIIc_gap_fraction", "between_isoform_difference",
        "discriminating_score", "position_class", "position_warning"]
SUMMARY_COLS = ["analysis_set", "position_class", "count"]


def col_stats(chars: List[str]):
    n = len(chars)
    residues = [c for c in chars if c not in M.GAP_CHARS]
    gap_frac = round((n - len(residues)) / n, 4) if n else 1.0
    if not residues:
        return "", 0.0, gap_frac
    cnt = Counter(residues)
    maj, mn = cnt.most_common(1)[0]
    return maj, round(mn / len(residues), 4), gap_frac


def classify(b_maj, b_frac, b_gap, c_maj, c_frac, c_gap) -> Tuple[str, str]:
    if b_gap >= 0.5 and c_gap >= 0.5:
        return "gap_rich_review", "both_isoforms_gap_rich"
    if b_gap >= 0.5 and c_frac >= 0.7:
        return "IIIc_specific_conserved", "IIIb_absent_or_gap_rich"
    if c_gap >= 0.5 and b_frac >= 0.7:
        return "IIIb_specific_conserved", "IIIc_absent_or_gap_rich"
    if b_frac >= 0.7 and c_frac >= 0.7:
        if b_maj != c_maj:
            return "isoform_discriminating_conserved", ""
        return "shared_conserved", ""
    return "variable", ""


def analyze(items: List[Tuple[str, str]], href_b: str, href_c: str):
    if not items:
        return []
    L = max(len(s) for _, s in items)
    b = [s for sid, s in items if sid.split("|")[1] == "IIIb"]
    c = [s for sid, s in items if sid.split("|")[1] == "IIIc"]
    hb = dict(items).get(href_b, "") if href_b else ""
    hc = dict(items).get(href_c, "") if href_c else ""
    rows = []
    for col in range(L):
        bc = [s[col] if col < len(s) else "-" for s in b]
        cc = [s[col] if col < len(s) else "-" for s in c]
        b_maj, b_frac, b_gap = col_stats(bc)
        c_maj, c_frac, c_gap = col_stats(cc)
        diff = 1 if (b_maj and c_maj and b_maj != c_maj) else 0
        disc = round(min(b_frac, c_frac) * diff * (1 - max(b_gap, c_gap)), 4)
        pclass, warn = classify(b_maj, b_frac, b_gap, c_maj, c_frac, c_gap)
        hba = hb[col] if hb and col < len(hb) and hb[col] not in M.GAP_CHARS else ""
        hca = hc[col] if hc and col < len(hc) and hc[col] not in M.GAP_CHARS else ""
        rows.append({"alignment_col": col + 1, "human_IIIb_aa_if_available": hba,
                     "human_IIIc_aa_if_available": hca, "IIIb_major_aa": b_maj,
                     "IIIc_major_aa": c_maj, "IIIb_major_aa_fraction": b_frac,
                     "IIIc_major_aa_fraction": c_frac, "IIIb_gap_fraction": b_gap,
                     "IIIc_gap_fraction": c_gap, "between_isoform_difference": diff,
                     "discriminating_score": disc, "position_class": pclass,
                     "position_warning": warn})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description="Isoform-discriminating residue analysis (Part 8).")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    cons_dir = dirs["conservation"]

    aln = dirs["alignments"] / "fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa"
    items_all = [(i, M.clean_alignment_seq(s)) for i, s in M.read_fasta(aln)]
    man = {r["msa_input_id"]: r for r in
           M.read_tsv(dirs["inputs"] / "fgfr2_IIIb_IIIc_combined_cassette_msa_input_manifest.tsv")}
    main_ids = {sid for sid in dict(items_all)
                if M.is_main_use(man.get(sid, {}).get("recommended_use", ""))}
    items_main = [(i, s) for i, s in items_all if i in main_ids]

    href_b = next((i for i, _ in items_all if i.startswith("homo_sapiens|IIIb")), "")
    href_c = next((i for i, _ in items_all if i.startswith("homo_sapiens|IIIc")), "")

    rows_main = analyze(items_main, href_b, href_c)
    rows_all = analyze(items_all, href_b, href_c)

    M.write_tsv(cons_dir / "fgfr2_IIIb_IIIc_discriminating_positions_main_only.tsv", rows_main, COLS)
    M.write_tsv(cons_dir / "fgfr2_IIIb_IIIc_discriminating_positions_all_species.tsv", rows_all, COLS)
    # base file == main-only (primary, less distorted by review/supplement sequences)
    M.write_tsv(cons_dir / "fgfr2_IIIb_IIIc_discriminating_positions.tsv", rows_main, COLS)

    summary: List[Dict[str, object]] = []
    for label, rows in (("main_only", rows_main), ("all_species", rows_all)):
        cc = Counter(r["position_class"] for r in rows)
        for k, v in sorted(cc.items()):
            summary.append({"analysis_set": label, "position_class": k, "count": v})
    M.write_tsv(cons_dir / "fgfr2_IIIb_IIIc_discriminating_positions_summary.tsv",
                summary, SUMMARY_COLS)

    n_disc = sum(1 for r in rows_main if r["position_class"] == "isoform_discriminating_conserved")
    print(f"[OK] discriminating residues: main_only cols={len(rows_main)} "
          f"(isoform_discriminating_conserved={n_disc}); all_species cols={len(rows_all)}")
    print(f"     analysis set reported: main_only (primary) + all_species (supplement)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
