#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402
from fgfr2 import reference_projection as RP  # noqa: E402

try:
    from Bio.Align import substitution_matrices
    _BLOSUM = substitution_matrices.load("BLOSUM62")
except Exception:  # noqa: BLE001
    _BLOSUM = None


MAP_COLS = ["msa_name", "alignment_col", "human_reference_isoform", "human_reference_residue_index",
            "human_reference_aa", "species", "isoform", "residue_aa", "is_gap",
            "native_cassette_residue_index", "native_protein_aa_position", "mapping_status"]
AGREE_COLS = ["species", "isoform", "final_isoform_label", "upstream_label",
              "label_consistency_status", "human_reference_isoform", "human_reference_residue_index",
              "alignment_col", "human_reference_aa", "species_aa", "agreement_class",
              "substitution_class", "blosum62_score_if_available", "is_gap", "is_review_species",
              "recommended_use", "human_reference_source", "human_reference_mapping_method",
              "position_warning"]
SP_SUM_COLS = ["species", "isoform", "final_isoform_label", "upstream_label",
               "label_consistency_status", "n_reference_positions", "n_identical", "n_conservative",
               "n_nonconservative", "n_gap_or_missing", "percent_identical",
               "percent_identical_or_conservative", "agreement_status", "agreement_warning"]
POS_SUM_COLS = ["human_reference_isoform", "human_reference_residue_index", "alignment_col",
                "human_reference_aa", "n_species", "n_identical", "n_conservative",
                "n_nonconservative", "n_gap", "percent_identical",
                "percent_identical_or_conservative", "position_agreement_class"]
DISC_COLS = ["human_reference_residue_index", "combined_alignment_col", "alignment_col",
             "human_IIIb_aa", "human_IIIc_aa", "human_IIIb_aa_if_available",
             "human_IIIc_aa_if_available", "human_IIIb_reference_index",
             "human_IIIc_reference_index", "human_reference_source",
             "IIIb_major_aa", "IIIc_major_aa",
             "IIIb_major_aa_fraction", "IIIc_major_aa_fraction", "IIIb_gap_fraction",
             "IIIc_gap_fraction", "between_isoform_difference", "discriminating_score",
             "position_class", "informative_column", "gap_rich_excluded_from_main_plot",
             "position_warning"]
SEG_COLS = ["species", "isoform", "segment_type", "n_positions", "percent_identical",
            "percent_identical_or_conservative", "percent_nonconservative", "gap_fraction",
            "segment_conservation_status", "segment_warning"]


def blosum(a: str, b: str) -> Optional[float]:
    if _BLOSUM is None or not a or not b:
        return None
    for k in ((a, b), (b, a)):
        try:
            return float(_BLOSUM[k])
        except Exception:
            continue
    return None


def classify(human_aa: str, sp_aa: str, is_gap: bool) -> Tuple[str, str, object]:
    if is_gap or sp_aa in ("", "-"):
        return "gap_or_missing", "gap", ""
    if not human_aa:
        return "unmapped_review", "unknown", ""
    sc = blosum(human_aa, sp_aa)
    if sp_aa == human_aa:
        return "identical_to_human", "identical", sc if sc is not None else ""
    if sc is None:
        return "nonconservative_substitution", "unknown", ""
    if sc > 0:
        return "conservative_substitution", "conservative", sc
    if sc == 0:
        return "nonconservative_substitution", "semi_conservative", sc
    return "nonconservative_substitution", "nonconservative", sc


def col_to_ungapped(seq: str) -> Dict[int, int]:
    out = {}
    idx = 0
    for col, ch in enumerate(seq):
        if ch not in M.GAP_CHARS:
            idx += 1
            out[col] = idx
    return out




def load_manifest(inp: Path, name: str) -> Dict[str, Dict[str, str]]:
    return {r["msa_input_id"]: r for r in M.read_tsv(inp / name)}


def main() -> int:
    ap = argparse.ArgumentParser(description="Reference-guided cassette agreement (Parts B-E).")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--window", type=int, default=5, help="boundary window size (aa)")
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    inp, maps, cons = dirs["inputs"], dirs["maps"], dirs["conservation"]
    W = args.window

    master = {r["species"].lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}
    man_iiib = load_manifest(inp, "fgfr2_IIIb_cassette_msa_input_manifest.tsv")
    man_iiic = load_manifest(inp, "fgfr2_IIIc_cassette_msa_input_manifest.tsv")
    man_comb = load_manifest(inp, "fgfr2_IIIb_IIIc_combined_cassette_msa_input_manifest.tsv")

    def aln(name):
        items = [(i, M.clean_alignment_seq(s)) for i, s in
                 M.read_fasta(dirs["alignments"] / name)]
        return items

    iiib_aln = aln("fgfr2_IIIb_cassette_msa_linsi.aln.faa") or aln("fgfr2_IIIb_cassette_msa.aln.faa")
    iiic_aln = aln("fgfr2_IIIc_cassette_msa_linsi.aln.faa") or aln("fgfr2_IIIc_cassette_msa.aln.faa")
    comb_aln = (aln("fgfr2_IIIb_IIIc_combined_cassette_msa_linsi.aln.faa")
                or aln("fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa"))

    all_agree_rows: Dict[str, List[Dict[str, object]]] = {"IIIb": [], "IIIc": [], "combined": []}
    all_map_rows: Dict[str, List[Dict[str, object]]] = {"IIIb": [], "IIIc": [], "combined": []}

    def is_review(ruse: str) -> str:
        return "false" if M.is_main_use(ruse) else "true"

    def process_single(items, manifest, isoform, msa_name, ref_label):
        projection = RP.resolve(isoform, items)
        href = projection.by_column
        seq_idx = {sid: col_to_ungapped(s) for sid, s in items}
        for sid, s in items:
            mrow = manifest.get(sid, {})
            sp = mrow.get("species", sid.split("|")[0])
            ruse = mrow.get("recommended_use", "")
            nstart = M.to_int(mrow.get("native_start_aa"), 0) or 0
            for col in range(len(s)):
                ch = s[col]
                is_gap = ch in M.GAP_CHARS
                h_idx, h_aa = href.get(col, (None, ""))
                sp_idx = seq_idx[sid].get(col)
                native_pos = (nstart + sp_idx - 1) if (sp_idx and nstart) else ""
                if h_idx is not None and not is_gap:
                    mstat = "mapped_to_human_reference_position"
                elif h_idx is not None and is_gap:
                    mstat = "deletion_relative_to_human"
                elif h_idx is None and not is_gap:
                    mstat = "insertion_relative_to_human"
                else:
                    mstat = "human_reference_gap"
                all_map_rows[ref_label].append({
                    "msa_name": msa_name, "alignment_col": col + 1,
                    "human_reference_isoform": isoform,
                    "human_reference_residue_index": h_idx if h_idx is not None else "",
                    "human_reference_aa": h_aa, "species": sp, "isoform": isoform,
                    "residue_aa": "" if is_gap else ch, "is_gap": str(is_gap).lower(),
                    "native_cassette_residue_index": sp_idx if sp_idx else "",
                    "native_protein_aa_position": native_pos, "mapping_status": mstat})
                # agreement only at human-reference positions
                if h_idx is not None:
                    aclass, sclass, sc = classify(h_aa, "" if is_gap else ch, is_gap)
                    all_agree_rows[ref_label].append({
                        "species": sp, "isoform": isoform, "final_isoform_label": isoform,
                        "upstream_label": mrow.get("upstream_label", ""),
                        "label_consistency_status": mrow.get("label_consistency_status", ""),
                        "human_reference_isoform": isoform,
                        "human_reference_residue_index": h_idx, "alignment_col": col + 1,
                        "human_reference_aa": h_aa, "species_aa": "" if is_gap else ch,
                        "agreement_class": aclass, "substitution_class": sclass,
                        "blosum62_score_if_available": sc, "is_gap": str(is_gap).lower(),
                        "is_review_species": is_review(ruse), "recommended_use": ruse,
                        "human_reference_source": projection.source,
                        "human_reference_mapping_method": projection.method,
                        "position_warning": ("" if projection.available
                                             else "human_reference_unavailable")})

    process_single(iiib_aln, man_iiib, "IIIb", "IIIb_cassette_linsi", "IIIb")
    process_single(iiic_aln, man_iiic, "IIIc", "IIIc_cassette_linsi", "IIIc")

    comb_proj = {iso: RP.resolve(iso, [(i, s) for i, s in comb_aln
                                       if (man_comb.get(i, {}).get("isoform")
                                           or (i.split("|")[1] if "|" in i else "")) == iso])
                 for iso in ("IIIb", "IIIc")}
    comb_href = {iso: comb_proj[iso].by_column for iso in ("IIIb", "IIIc")}
    comb_idx = {sid: col_to_ungapped(s) for sid, s in comb_aln}
    for sid, s in comb_aln:
        mrow = man_comb.get(sid, {})
        sp = mrow.get("species", sid.split("|")[0])
        iso = mrow.get("isoform", sid.split("|")[1] if "|" in sid else "")
        ruse = mrow.get("recommended_use", "")
        nstart = M.to_int(mrow.get("native_start_aa"), 0) or 0
        href = comb_href.get(iso, {})
        for col in range(len(s)):
            ch = s[col]
            is_gap = ch in M.GAP_CHARS
            h_idx, h_aa = href.get(col, (None, ""))
            sp_idx = comb_idx[sid].get(col)
            native_pos = (nstart + sp_idx - 1) if (sp_idx and nstart) else ""
            if h_idx is not None and not is_gap:
                mstat = "mapped_to_human_reference_position"
            elif h_idx is not None and is_gap:
                mstat = "deletion_relative_to_human"
            elif h_idx is None and not is_gap:
                mstat = "insertion_relative_to_human"
            else:
                mstat = "human_reference_gap"
            all_map_rows["combined"].append({
                "msa_name": "combined_cassette_linsi", "alignment_col": col + 1,
                "human_reference_isoform": iso,
                "human_reference_residue_index": h_idx if h_idx is not None else "",
                "human_reference_aa": h_aa, "species": sp, "isoform": iso,
                "residue_aa": "" if is_gap else ch, "is_gap": str(is_gap).lower(),
                "native_cassette_residue_index": sp_idx if sp_idx else "",
                "native_protein_aa_position": native_pos, "mapping_status": mstat})
            if h_idx is not None:
                aclass, sclass, sc = classify(h_aa, "" if is_gap else ch, is_gap)
                all_agree_rows["combined"].append({
                    "species": sp, "isoform": iso, "final_isoform_label": iso,
                    "upstream_label": mrow.get("upstream_label", ""),
                    "label_consistency_status": mrow.get("label_consistency_status", ""),
                    "human_reference_isoform": iso,
                    "human_reference_residue_index": h_idx, "alignment_col": col + 1,
                    "human_reference_aa": h_aa, "species_aa": "" if is_gap else ch,
                    "agreement_class": aclass, "substitution_class": sclass,
                    "blosum62_score_if_available": sc, "is_gap": str(is_gap).lower(),
                    "is_review_species": is_review(ruse), "recommended_use": ruse,
                    "human_reference_source": comb_proj[iso].source if iso in comb_proj else "",
                    "human_reference_mapping_method": (comb_proj[iso].method
                                                       if iso in comb_proj else ""),
                    "position_warning": ""})

    M.write_tsv(maps / "fgfr2_IIIb_human_reference_msa_coordinate_map.tsv", all_map_rows["IIIb"], MAP_COLS)
    M.write_tsv(maps / "fgfr2_IIIc_human_reference_msa_coordinate_map.tsv", all_map_rows["IIIc"], MAP_COLS)
    M.write_tsv(maps / "fgfr2_combined_human_reference_msa_coordinate_map.tsv", all_map_rows["combined"], MAP_COLS)
    M.write_tsv(cons / "fgfr2_IIIb_human_reference_residue_agreement.tsv", all_agree_rows["IIIb"], AGREE_COLS)
    M.write_tsv(cons / "fgfr2_IIIc_human_reference_residue_agreement.tsv", all_agree_rows["IIIc"], AGREE_COLS)
    M.write_tsv(cons / "fgfr2_combined_human_reference_residue_agreement.tsv", all_agree_rows["combined"], AGREE_COLS)

    def summary_by_species(agree_rows) -> List[Dict[str, object]]:
        by = defaultdict(list)
        for r in agree_rows:
            by[(r["species"], r["isoform"])].append(r)
        out = []
        for (sp, iso), rs in by.items():
            n = len(rs)
            n_id = sum(1 for r in rs if r["agreement_class"] == "identical_to_human")
            n_co = sum(1 for r in rs if r["agreement_class"] == "conservative_substitution")
            n_nc = sum(1 for r in rs if r["agreement_class"] == "nonconservative_substitution")
            n_gap = sum(1 for r in rs if r["agreement_class"] == "gap_or_missing")
            pid = round(n_id / n, 4) if n else 0.0
            pic = round((n_id + n_co) / n, 4) if n else 0.0
            gapf = (n_gap / n) if n else 1.0
            if gapf >= 0.5:
                status, warn = "gap_rich_review", f"gap_fraction={round(gapf,3)}"
            elif pic >= 0.9:
                status, warn = "high_reference_agreement", ""
            elif pic >= 0.7:
                status, warn = "moderate_reference_agreement", ""
            else:
                status, warn = "low_reference_agreement_review", f"percent_id_or_cons={pic}"
            out.append({"species": sp, "isoform": iso,
                        "final_isoform_label": iso, "upstream_label": rs[0].get("upstream_label", ""),
                        "label_consistency_status": rs[0].get("label_consistency_status", ""),
                        "n_reference_positions": n,
                        "n_identical": n_id, "n_conservative": n_co, "n_nonconservative": n_nc,
                        "n_gap_or_missing": n_gap, "percent_identical": pid,
                        "percent_identical_or_conservative": pic, "agreement_status": status,
                        "agreement_warning": warn})
        return sorted(out, key=lambda d: (M.to_int(master.get(d["species"].lower(), {}).get("phylo_order"), 999), d["isoform"]))

    M.write_tsv(cons / "fgfr2_IIIb_reference_agreement_summary_by_species.tsv",
                summary_by_species(all_agree_rows["IIIb"]), SP_SUM_COLS)
    M.write_tsv(cons / "fgfr2_IIIc_reference_agreement_summary_by_species.tsv",
                summary_by_species(all_agree_rows["IIIc"]), SP_SUM_COLS)

    pos_rows = []
    for iso, rows in (("IIIb", all_agree_rows["IIIb"]), ("IIIc", all_agree_rows["IIIc"])):
        by = defaultdict(list)
        for r in rows:
            by[r["human_reference_residue_index"]].append(r)
        for hidx, rs in sorted(by.items(), key=lambda x: int(x[0])):
            n = len(rs)
            n_id = sum(1 for r in rs if r["agreement_class"] == "identical_to_human")
            n_co = sum(1 for r in rs if r["agreement_class"] == "conservative_substitution")
            n_nc = sum(1 for r in rs if r["agreement_class"] == "nonconservative_substitution")
            n_gap = sum(1 for r in rs if r["agreement_class"] == "gap_or_missing")
            pid = round(n_id / n, 4) if n else 0.0
            pic = round((n_id + n_co) / n, 4) if n else 0.0
            pclass = ("position_highly_conserved" if pid >= 0.8 else
                      "position_conservative" if pic >= 0.8 else
                      "position_gap_rich" if (n_gap / n) >= 0.5 else "position_variable")
            pos_rows.append({"human_reference_isoform": iso, "human_reference_residue_index": hidx,
                             "alignment_col": rs[0]["alignment_col"],
                             "human_reference_aa": rs[0]["human_reference_aa"], "n_species": n,
                             "n_identical": n_id, "n_conservative": n_co, "n_nonconservative": n_nc,
                             "n_gap": n_gap, "percent_identical": pid,
                             "percent_identical_or_conservative": pic,
                             "position_agreement_class": pclass})
    M.write_tsv(cons / "fgfr2_reference_agreement_summary_by_position.tsv", pos_rows, POS_SUM_COLS)

    LH = {"IIIb": 51, "IIIc": 42}
    seg_rows = []
    for iso, agree_rows, map_rows in (("IIIb", all_agree_rows["IIIb"], all_map_rows["IIIb"]),
                                      ("IIIc", all_agree_rows["IIIc"], all_map_rows["IIIc"])):
        lh = max([M.to_int(r["human_reference_residue_index"], 0) for r in agree_rows] + [LH[iso]])
        by_sp = defaultdict(list)
        for r in agree_rows:
            by_sp[r["species"]].append(r)
        ins_by_sp = defaultdict(int)
        for r in map_rows:
            if r["mapping_status"] == "insertion_relative_to_human":
                ins_by_sp[r["species"]] += 1

        def seg_of(hidx: int) -> str:
            if hidx <= W:
                return "left_boundary_window"
            if hidx > lh - W:
                return "right_boundary_window"
            return "cassette_core"

        for sp, rs in by_sp.items():
            buckets = defaultdict(list)
            for r in rs:
                hidx = M.to_int(r["human_reference_residue_index"], 0)
                buckets[seg_of(hidx)].append(r)
                buckets["full_cassette"].append(r)
            for seg_type, srs in buckets.items():
                n = len(srs)
                n_id = sum(1 for r in srs if r["agreement_class"] == "identical_to_human")
                n_co = sum(1 for r in srs if r["agreement_class"] == "conservative_substitution")
                n_nc = sum(1 for r in srs if r["agreement_class"] == "nonconservative_substitution")
                n_gap = sum(1 for r in srs if r["agreement_class"] == "gap_or_missing")
                pid = round(n_id / n, 4) if n else 0.0
                pic = round((n_id + n_co) / n, 4) if n else 0.0
                pnc = round(n_nc / n, 4) if n else 0.0
                gapf = round(n_gap / n, 4) if n else 1.0
                if gapf >= 0.5:
                    st, wn = "segment_gap_rich_review", f"gap_fraction={gapf}"
                elif pic >= 0.9:
                    st, wn = "segment_well_conserved", ""
                elif pic >= 0.7:
                    st, wn = "segment_moderately_conserved", ""
                else:
                    st, wn = "segment_divergent_review", f"percent_id_or_cons={pic}"
                seg_rows.append({"species": sp, "isoform": iso, "segment_type": seg_type,
                                 "n_positions": n, "percent_identical": pid,
                                 "percent_identical_or_conservative": pic,
                                 "percent_nonconservative": pnc, "gap_fraction": gapf,
                                 "segment_conservation_status": st, "segment_warning": wn})
            if ins_by_sp.get(sp):
                seg_rows.append({"species": sp, "isoform": iso, "segment_type": "insertion_region",
                                 "n_positions": ins_by_sp[sp], "percent_identical": "",
                                 "percent_identical_or_conservative": "", "percent_nonconservative": "",
                                 "gap_fraction": "", "segment_conservation_status":
                                 "segment_insertion_relative_to_human",
                                 "segment_warning": "residues inserted relative to human reference"})
    M.write_tsv(cons / "fgfr2_cassette_segment_agreement_summary.tsv", seg_rows, SEG_COLS)

    main_ids = {sid for sid, _ in comb_aln
                if M.is_main_use(man_comb.get(sid, {}).get("recommended_use", ""))}

    def discriminating(items, restrict_main: bool):
        if not items:
            return []
        L = max(len(s) for _, s in items)
        b = [s for sid, s in items if sid.split("|")[1] == "IIIb"
             and (not restrict_main or sid in main_ids)]
        c = [s for sid, s in items if sid.split("|")[1] == "IIIc"
             and (not restrict_main or sid in main_ids)]
        out = []
        for col in range(L):
            def stat(group):
                chars = [s[col] if col < len(s) else "-" for s in group]
                res = [x for x in chars if x not in M.GAP_CHARS]
                gf = round((len(chars) - len(res)) / len(chars), 4) if chars else 1.0
                if not res:
                    return "", 0.0, gf
                cnt = Counter(res)
                mj, mn = cnt.most_common(1)[0]
                return mj, round(mn / len(res), 4), gf
            b_maj, b_frac, b_gap = stat(b)
            c_maj, c_frac, c_gap = stat(c)
            diff = 1 if (b_maj and c_maj and b_maj != c_maj) else 0
            disc = round(min(b_frac, c_frac) * diff * (1 - max(b_gap, c_gap)), 4)
            if b_gap >= 0.5 and c_gap >= 0.5:
                pclass = "gap_rich_review"
            elif b_gap >= 0.5 and c_frac >= 0.7:
                pclass = "IIIc_specific_conserved"
            elif c_gap >= 0.5 and b_frac >= 0.7:
                pclass = "IIIb_specific_conserved"
            elif b_frac >= 0.7 and c_frac >= 0.7:
                pclass = "isoform_discriminating_conserved" if b_maj != c_maj else "shared_conserved"
            else:
                pclass = "variable"
            informative = (b_gap < 0.5 and c_gap < 0.5)
            # IIIb and IIIc number their cassettes independently (46 vs 48 residues),
            # so each panel keeps its own reference index. Collapsing them onto one
            # shared axis is what truncated IIIc and destroyed its marker.
            b_idx, hb_aa = comb_proj["IIIb"].at(col)
            c_idx, hc_aa = comb_proj["IIIc"].at(col)
            href_idx = b_idx if b_idx is not None else (c_idx if c_idx is not None else "")
            resolved = comb_proj["IIIb"].available or comb_proj["IIIc"].available
            out.append({"human_reference_residue_index": href_idx, "combined_alignment_col": col + 1,
                        "alignment_col": col + 1, "human_IIIb_aa": hb_aa, "human_IIIc_aa": hc_aa,
                        "human_IIIb_aa_if_available": hb_aa, "human_IIIc_aa_if_available": hc_aa,
                        "human_IIIb_reference_index": b_idx if b_idx is not None else "",
                        "human_IIIc_reference_index": c_idx if c_idx is not None else "",
                        "human_reference_source": comb_proj["IIIb"].source,
                        "IIIb_major_aa": b_maj, "IIIc_major_aa": c_maj,
                        "IIIb_major_aa_fraction": b_frac, "IIIc_major_aa_fraction": c_frac,
                        "IIIb_gap_fraction": b_gap, "IIIc_gap_fraction": c_gap,
                        "between_isoform_difference": diff, "discriminating_score": disc,
                        "position_class": pclass, "informative_column": str(informative).lower(),
                        "gap_rich_excluded_from_main_plot": str(not informative).lower(),
                        "position_warning": "" if resolved else "human_reference_mapping_failed"})
        return out

    disc_main = discriminating(comb_aln, restrict_main=True)
    disc_all = discriminating(comb_aln, restrict_main=False)
    disc_inf = [r for r in disc_main if r["informative_column"] == "true"]
    M.write_tsv(cons / "fgfr2_IIIb_IIIc_discriminating_positions_main_only.tsv", disc_main, DISC_COLS)
    M.write_tsv(cons / "fgfr2_IIIb_IIIc_discriminating_positions_all_species.tsv", disc_all, DISC_COLS)
    M.write_tsv(cons / "fgfr2_IIIb_IIIc_discriminating_positions_informative.tsv", disc_inf, DISC_COLS)

    n_disc = sum(1 for r in disc_inf if r["position_class"] == "isoform_discriminating_conserved")
    print(f"[OK] reference agreement built (BLOSUM62={'yes' if _BLOSUM is not None else 'no'})")
    print(f"     maps: IIIb={len(all_map_rows['IIIb'])} IIIc={len(all_map_rows['IIIc'])} "
          f"combined={len(all_map_rows['combined'])} rows")
    print(f"     agreement: IIIb={len(all_agree_rows['IIIb'])} IIIc={len(all_agree_rows['IIIc'])} positions")
    print(f"     segments={len(seg_rows)}; discriminating informative cols={len(disc_inf)} "
          f"(isoform_discriminating_conserved={n_disc})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
