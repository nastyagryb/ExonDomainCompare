#!/usr/bin/env python3
"""
analyze_fgfr2_msa_conservation.py  (MSA boundary-robustness sprint, Part 6)

Per-column conservation and gap scores for each MSA, plus region-level summaries
(cassette / left & right boundary windows / full sequence) per species/isoform.

conservation_score = 1 - normalized Shannon entropy over non-gap residues. The gap
fraction is reported SEPARATELY and never hidden inside the conservation score.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fgfr2_msa_common as M  # noqa: E402


ALN = [
    ("full_length", "fgfr2_full_length_protein_msa.aln.faa",
     "fgfr2_full_length_msa_column_conservation.tsv"),
    ("IIIb_cassette", "fgfr2_IIIb_cassette_msa.aln.faa",
     "fgfr2_IIIb_msa_column_conservation.tsv"),
    ("IIIc_cassette", "fgfr2_IIIc_cassette_msa.aln.faa",
     "fgfr2_IIIc_msa_column_conservation.tsv"),
    ("IIIb_IIIc_combined_cassette", "fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa",
     "fgfr2_IIIb_IIIc_combined_msa_column_conservation.tsv"),
]
COL_COLS = ["msa_name", "alignment_col", "n_sequences", "n_non_gap", "gap_fraction",
            "major_aa", "major_aa_fraction", "shannon_entropy", "normalized_entropy",
            "conservation_score", "conservation_class", "human_reference_aa_if_available",
            "in_boundary_window", "in_cassette_region"]
REGION_COLS = ["species", "isoform", "msa_name", "region_type", "mean_conservation_score",
               "mean_gap_fraction", "n_alignment_columns", "conservation_status",
               "conservation_warning"]
LOG2_20 = math.log2(20)


def human_id(items: List[Tuple[str, str]], isoform: Optional[str]) -> Optional[str]:
    cand = [sid for sid, _ in items if sid.startswith("homo_sapiens|")]
    if isoform:
        pref = [c for c in cand if c.split("|")[1] == isoform]
        if pref:
            return pref[0]
    return cand[0] if cand else None


def column_stats(col_chars: List[str]):
    residues = [c for c in col_chars if c not in M.GAP_CHARS]
    n = len(col_chars)
    n_ng = len(residues)
    gap_fraction = round((n - n_ng) / n, 4) if n else 1.0
    if not residues:
        return n, 0, gap_fraction, "", 0.0, 0.0, 1.0, 0.0
    cnt = Counter(residues)
    major_aa, major_n = cnt.most_common(1)[0]
    major_frac = round(major_n / n_ng, 4)
    ent = 0.0
    for _, c in cnt.items():
        p = c / n_ng
        ent -= p * math.log2(p)
    norm_ent = round(ent / LOG2_20, 4)
    cons = round(1.0 - norm_ent, 4)
    return n, n_ng, gap_fraction, major_aa, major_frac, round(ent, 4), norm_ent, cons


def cons_class(cons: float, gap_fraction: float) -> str:
    if gap_fraction >= 0.5:
        return "gap_rich_review"
    if cons >= 0.8:
        return "highly_conserved"
    if cons >= 0.5:
        return "moderately_conserved"
    return "variable"


def main() -> int:
    ap = argparse.ArgumentParser(description="MSA per-column conservation (Part 6).")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    cons_dir, maps = dirs["conservation"], dirs["maps"]

    proj = M.read_tsv(maps / "fgfr2_exon_boundary_msa_projection.tsv")
    proj_by = {(p["species"].lower(), p["isoform"]): p for p in proj}

    # human reference cassette region (full-length) for in_cassette_region flag
    human_iiib = proj_by.get(("homo_sapiens", "IIIb"), {})
    fl_cass_start = M.to_int(human_iiib.get("full_length_msa_start_col"))
    fl_cass_end = M.to_int(human_iiib.get("full_length_msa_end_col"))

    col_cons_cache: Dict[str, Dict[int, Dict[str, object]]] = {}

    for msa_name, aln_file, out_name in ALN:
        items = [(i, M.clean_alignment_seq(s)) for i, s in M.read_fasta(dirs["alignments"] / aln_file)]
        if not items:
            M.write_tsv(cons_dir / out_name, [], COL_COLS)
            continue
        L = max(len(s) for _, s in items)
        n_seq = len(items)
        iso_for_human = ("IIIb" if msa_name == "IIIb_cassette" else
                         "IIIc" if msa_name == "IIIc_cassette" else "IIIb")
        href = human_id(items, iso_for_human)
        href_seq = dict(items).get(href, "") if href else ""
        is_cassette_msa = "cassette" in msa_name
        rows: List[Dict[str, object]] = []
        col_map: Dict[int, Dict[str, object]] = {}
        for c in range(L):
            chars = [s[c] if c < len(s) else "-" for _, s in items]
            _n, n_ng, gf, maj, majf, ent, nent, cons = column_stats(chars)
            klass = cons_class(cons, gf)
            href_aa = (href_seq[c] if href_seq and c < len(href_seq) else "")
            # boundary/cassette flags
            if is_cassette_msa:
                in_cass = "true"
                in_bw = "true" if (c < 5 or c >= L - 5) else "false"
            else:
                in_cass = ("true" if (fl_cass_start and fl_cass_end
                                      and fl_cass_start - 1 <= c <= fl_cass_end - 1) else "false")
                in_bw = "false"
                for bcol in (fl_cass_start, fl_cass_end):
                    if bcol and abs((c + 1) - bcol) <= 5:
                        in_bw = "true"
            row = {"msa_name": msa_name, "alignment_col": c + 1, "n_sequences": n_seq,
                   "n_non_gap": n_ng, "gap_fraction": gf, "major_aa": maj,
                   "major_aa_fraction": majf, "shannon_entropy": ent, "normalized_entropy": nent,
                   "conservation_score": cons, "conservation_class": klass,
                   "human_reference_aa_if_available": href_aa if href_aa not in M.GAP_CHARS else "",
                   "in_boundary_window": in_bw, "in_cassette_region": in_cass}
            rows.append(row)
            col_map[c + 1] = {"cons": cons, "gap": gf}
        M.write_tsv(cons_dir / out_name, rows, COL_COLS)
        col_cons_cache[msa_name] = col_map

    # ---- region-level summary per species/isoform ----
    full_cols = col_cons_cache.get("full_length", {})

    def region_stat(cols: List[int]):
        vals = [(full_cols[c]["cons"], full_cols[c]["gap"]) for c in cols if c in full_cols]
        if not vals:
            return "", "", 0
        mc = round(sum(v[0] for v in vals) / len(vals), 4)
        mg = round(sum(v[1] for v in vals) / len(vals), 4)
        return mc, mg, len(vals)

    region_rows: List[Dict[str, object]] = []
    for (sp, iso), p in proj_by.items():
        fs = M.to_int(p.get("full_length_msa_start_col"))
        fe = M.to_int(p.get("full_length_msa_end_col"))
        regions = {}
        if fs and fe and fe >= fs:
            regions["cassette"] = list(range(fs, fe + 1))
            regions["left_boundary_window"] = list(range(max(1, fs - 5), fs + 6))
            regions["right_boundary_window"] = list(range(max(1, fe - 5), fe + 6))
        # full_sequence: span of this species' residues in full-length MSA
        regions["full_sequence"] = list(range(1, len(full_cols) + 1)) if full_cols else []
        for rtype, cols in regions.items():
            mc, mg, ncol = region_stat(cols)
            warn = ""
            if mg != "" and mg >= 0.5:
                status = f"{rtype}_gap_rich_review"
                warn = f"mean_gap_fraction={mg}"
            elif mc != "" and mc >= 0.8:
                status = f"{rtype}_highly_conserved"
            elif mc != "" and mc >= 0.5:
                status = f"{rtype}_moderately_conserved"
            elif mc == "":
                status = f"{rtype}_unavailable"
            else:
                status = f"{rtype}_variable"
            region_rows.append({"species": sp, "isoform": iso, "msa_name": "full_length",
                                "region_type": rtype, "mean_conservation_score": mc,
                                "mean_gap_fraction": mg, "n_alignment_columns": ncol,
                                "conservation_status": status, "conservation_warning": warn})
    M.write_tsv(cons_dir / "fgfr2_msa_region_conservation_summary.tsv", region_rows, REGION_COLS)

    print(f"[OK] column conservation for {len(ALN)} MSAs; "
          f"region summary rows={len(region_rows)}")
    for msa_name, _, out_name in ALN:
        cm = col_cons_cache.get(msa_name, {})
        if cm:
            mc = round(sum(v["cons"] for v in cm.values()) / len(cm), 3)
            print(f"     {msa_name}: cols={len(cm)} mean_conservation={mc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
