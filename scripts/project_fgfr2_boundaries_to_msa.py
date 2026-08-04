#!/usr/bin/env python3

from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402


ALN = {
    "full_length": ("fgfr2_full_length_protein_msa.aln.faa", "full_length_protein",
                    "fgfr2_full_length_msa_coordinate_map.tsv"),
    "IIIb_cassette": ("fgfr2_IIIb_cassette_msa.aln.faa", "IIIb_cassette",
                      "fgfr2_IIIb_cassette_msa_coordinate_map.tsv"),
    "IIIc_cassette": ("fgfr2_IIIc_cassette_msa.aln.faa", "IIIc_cassette",
                      "fgfr2_IIIc_cassette_msa_coordinate_map.tsv"),
    "IIIb_IIIc_combined_cassette": ("fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa",
                                    "IIIb_IIIc_combined_cassette",
                                    "fgfr2_IIIb_IIIc_combined_cassette_msa_coordinate_map.tsv"),
}
MAP_COLS = ["msa_name", "msa_input_id", "species", "isoform", "protein_id", "transcript_id",
            "alignment_col", "residue_index_ungapped", "residue_aa", "is_gap",
            "native_protein_aa_position", "sequence_type", "recommended_use"]
VAL_COLS = ["msa_name", "msa_input_id", "check", "status", "expected", "observed", "detail"]
PROJ_COLS = ["species", "display_species_name", "isoform", "protein_id", "transcript_id",
             "recommended_use", "native_cassette_start_aa", "native_cassette_end_aa",
             "full_length_msa_start_col", "full_length_msa_end_col",
             "cassette_msa_start_col", "cassette_msa_end_col",
             "left_boundary_gap_fraction_window", "right_boundary_gap_fraction_window",
             "internal_region_gap_fraction", "boundary_window_size",
             "left_boundary_gap_fraction_w3", "right_boundary_gap_fraction_w3",
             "left_boundary_gap_fraction_w5", "right_boundary_gap_fraction_w5",
             "boundary_projection_status", "boundary_projection_confidence",
             "boundary_projection_warning"]


def manifest_lookup(inp: Path, name: str) -> Dict[str, Dict[str, str]]:
    fn = {
        "full_length_protein": "fgfr2_full_length_protein_msa_input_manifest.tsv",
        "IIIb_cassette": "fgfr2_IIIb_cassette_msa_input_manifest.tsv",
        "IIIc_cassette": "fgfr2_IIIc_cassette_msa_input_manifest.tsv",
        "IIIb_IIIc_combined_cassette": "fgfr2_IIIb_IIIc_combined_cassette_msa_input_manifest.tsv",
    }[name]
    return {r["msa_input_id"]: r for r in M.read_tsv(inp / fn)}


def ungapped_to_col(aligned: str) -> Dict[int, int]:
    m: Dict[int, int] = {}
    idx = 0
    for col, ch in enumerate(aligned, start=1):
        if ch not in M.GAP_CHARS:
            idx += 1
            m[idx] = col
    return m


def col_gap_fractions(items: List[Tuple[str, str]]) -> List[float]:
    if not items:
        return []
    L = max(len(s) for _, s in items)
    n = len(items)
    out = []
    for c in range(L):
        g = sum(1 for _, s in items if c >= len(s) or s[c] in M.GAP_CHARS)
        out.append(g / n)
    return out


def window_gap(gaps: List[float], col: Optional[int], w: int) -> Optional[float]:
    if col is None or not gaps:
        return None
    lo, hi = max(0, col - 1 - w), min(len(gaps), col - 1 + w + 1)
    seg = gaps[lo:hi]
    return round(sum(seg) / len(seg), 4) if seg else None


def main() -> int:
    ap = argparse.ArgumentParser(description="MSA coordinate maps + boundary projection (Parts 4,5).")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    inp, maps = dirs["inputs"], dirs["maps"]

    master = {r["species"].lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}

    aligned_cache: Dict[str, List[Tuple[str, str]]] = {}
    u2c_cache: Dict[str, Dict[str, Dict[int, int]]] = {}
    gaps_cache: Dict[str, List[float]] = {}
    man_cache: Dict[str, Dict[str, Dict[str, str]]] = {}
    val_rows: List[Dict[str, object]] = []

    # Build and validate coordinate maps.
    for msa_name, (aln_file, stype_name, map_name) in ALN.items():
        aln_path = dirs["alignments"] / aln_file
        items = [(i, M.clean_alignment_seq(s)) for i, s in M.read_fasta(aln_path)]
        aligned_cache[msa_name] = items
        gaps_cache[msa_name] = col_gap_fractions(items)
        man = manifest_lookup(inp, stype_name)
        man_cache[msa_name] = man
        u2c_cache[msa_name] = {}
        rows: List[Dict[str, object]] = []
        for sid, aseq in items:
            mrow = man.get(sid, {})
            sp = mrow.get("species", sid.split("|")[0])
            iso = mrow.get("isoform", "")
            pid = mrow.get("protein_id", "")
            tx = mrow.get("transcript_id", "")
            ruse = mrow.get("recommended_use", "")
            native_start = M.to_int(mrow.get("native_start_aa"), 1) or 1
            u2c = ungapped_to_col(aseq)
            u2c_cache[msa_name][sid] = u2c
            ridx = 0
            for col, ch in enumerate(aseq, start=1):
                is_gap = ch in M.GAP_CHARS
                if is_gap:
                    rows.append({"msa_name": msa_name, "msa_input_id": sid, "species": sp,
                                 "isoform": iso, "protein_id": pid, "transcript_id": tx,
                                 "alignment_col": col, "residue_index_ungapped": "",
                                 "residue_aa": "-", "is_gap": "true",
                                 "native_protein_aa_position": "", "sequence_type": stype_name,
                                 "recommended_use": ruse})
                else:
                    ridx += 1
                    rows.append({"msa_name": msa_name, "msa_input_id": sid, "species": sp,
                                 "isoform": iso, "protein_id": pid, "transcript_id": tx,
                                 "alignment_col": col, "residue_index_ungapped": ridx,
                                 "residue_aa": ch, "is_gap": "false",
                                 "native_protein_aa_position": native_start + ridx - 1,
                                 "sequence_type": stype_name, "recommended_use": ruse})
            # validation: non-gap count == input length; last native pos == manifest native_end
            exp_len = M.to_int(mrow.get("sequence_length"))
            obs_len = ridx
            val_rows.append({"msa_name": msa_name, "msa_input_id": sid,
                             "check": "non_gap_equals_input_length",
                             "status": "ok" if (exp_len is None or exp_len == obs_len) else "fail",
                             "expected": exp_len if exp_len is not None else "", "observed": obs_len,
                             "detail": ""})
            exp_end = M.to_int(mrow.get("native_end_aa"))
            obs_end = native_start + obs_len - 1 if obs_len else ""
            val_rows.append({"msa_name": msa_name, "msa_input_id": sid,
                             "check": "native_end_matches_manifest",
                             "status": "ok" if (exp_end is None or exp_end == obs_end) else "fail",
                             "expected": exp_end if exp_end is not None else "",
                             "observed": obs_end, "detail": ""})
        M.write_tsv(maps / map_name, rows, MAP_COLS)

    M.write_tsv(maps / "msa_coordinate_map_validation.tsv", val_rows, VAL_COLS)

    # Project cassette boundaries.
    full_items = aligned_cache["full_length"]
    full_u2c = u2c_cache["full_length"]
    full_gaps = gaps_cache["full_length"]
    cas_man = {"IIIb": man_cache["IIIb_cassette"], "IIIc": man_cache["IIIc_cassette"]}
    cas_u2c = {"IIIb": u2c_cache["IIIb_cassette"], "IIIc": u2c_cache["IIIc_cassette"]}
    _cas_items = {"IIIb": dict(aligned_cache["IIIb_cassette"]),
                 "IIIc": dict(aligned_cache["IIIc_cassette"])}

    proj: List[Dict[str, object]] = []
    # collect full-length start cols per isoform among MAIN species for shift detection
    main_start_cols: Dict[str, List[int]] = {"IIIb": [], "IIIc": []}

    # first pass to compute start cols
    prelim: List[Dict[str, object]] = []
    for sid, _ in full_items:
        mrow = man_cache["full_length"].get(sid, {})
        sp, iso = mrow.get("species", ""), mrow.get("isoform", "")
        cas_m = cas_man.get(iso, {}).get(sid, {})
        cstart = M.to_int(cas_m.get("native_start_aa"))
        cend = M.to_int(cas_m.get("native_end_aa"))
        u2c = full_u2c.get(sid, {})
        fl_start = u2c.get(cstart) if cstart else None
        fl_end = u2c.get(cend) if cend else None
        prelim.append({"sid": sid, "sp": sp, "iso": iso, "cstart": cstart, "cend": cend,
                       "fl_start": fl_start, "fl_end": fl_end, "ruse": mrow.get("recommended_use", "")})
        if fl_start and M.is_main_use(mrow.get("recommended_use", "")):
            main_start_cols[iso].append(fl_start)
    med_start = {iso: (statistics.median(v) if v else None) for iso, v in main_start_cols.items()}

    for p in prelim:
        sid, sp, iso = p["sid"], p["sp"], p["iso"]
        mr = master.get(sp.lower(), {})
        cstart, cend = p["cstart"], p["cend"]
        fl_start, fl_end = p["fl_start"], p["fl_end"]
        ruse = p["ruse"]
        # cassette-only projection: first/last non-gap col of this id
        cu2c = cas_u2c.get(iso, {}).get(sid, {})
        cas_start_col = min(cu2c.values()) if cu2c else None
        cas_end_col = max(cu2c.values()) if cu2c else None

        lw5 = window_gap(full_gaps, fl_start, 5)
        rw5 = window_gap(full_gaps, fl_end, 5)
        lw3 = window_gap(full_gaps, fl_start, 3)
        rw3 = window_gap(full_gaps, fl_end, 3)
        internal = None
        if fl_start and fl_end and fl_end >= fl_start and full_gaps:
            seg = full_gaps[fl_start - 1:fl_end]
            internal = round(sum(seg) / len(seg), 4) if seg else None

        warn = ""
        if fl_start is None or fl_end is None:
            status, conf = "msa_boundary_unresolved", "none"
            warn = "cassette boundary could not be projected onto full-length MSA"
        else:
            local = max(lw5 or 0.0, rw5 or 0.0)
            shift = (med_start.get(iso) is not None and abs(fl_start - med_start[iso]) > 40)
            if local <= 0.2 and not shift:
                status, conf = "msa_boundary_projected_high_confidence", "high"
            elif local <= 0.5 and not shift:
                status, conf = "msa_boundary_projected_with_minor_gaps", "minor_gaps"
            elif shift:
                status, conf = "msa_boundary_shift_review", "shift_review"
                warn = f"full-length start col {fl_start} deviates from group median {med_start[iso]}"
            else:
                status, conf = "msa_boundary_projected_gap_rich_review", "gap_rich"
                warn = f"boundary window heavily gapped (max_window_gap={round(local,3)})"

        proj.append({
            "species": sp, "display_species_name": mr.get("display_species_name", sp),
            "isoform": iso, "protein_id": man_cache["full_length"].get(sid, {}).get("protein_id", ""),
            "transcript_id": man_cache["full_length"].get(sid, {}).get("transcript_id", ""),
            "recommended_use": ruse, "native_cassette_start_aa": cstart if cstart else "",
            "native_cassette_end_aa": cend if cend else "",
            "full_length_msa_start_col": fl_start if fl_start else "",
            "full_length_msa_end_col": fl_end if fl_end else "",
            "cassette_msa_start_col": cas_start_col if cas_start_col else "",
            "cassette_msa_end_col": cas_end_col if cas_end_col else "",
            "left_boundary_gap_fraction_window": lw5 if lw5 is not None else "",
            "right_boundary_gap_fraction_window": rw5 if rw5 is not None else "",
            "internal_region_gap_fraction": internal if internal is not None else "",
            "boundary_window_size": "5",
            "left_boundary_gap_fraction_w3": lw3 if lw3 is not None else "",
            "right_boundary_gap_fraction_w3": rw3 if rw3 is not None else "",
            "left_boundary_gap_fraction_w5": lw5 if lw5 is not None else "",
            "right_boundary_gap_fraction_w5": rw5 if rw5 is not None else "",
            "boundary_projection_status": status, "boundary_projection_confidence": conf,
            "boundary_projection_warning": warn,
        })

    M.write_tsv(maps / "fgfr2_exon_boundary_msa_projection.tsv", proj, PROJ_COLS)

    n_val_fail = sum(1 for v in val_rows if v["status"] == "fail")
    from collections import Counter
    print(f"[OK] coordinate maps written (validation fails={n_val_fail})")
    print(f"     boundary projection rows={len(proj)} "
          f"status={dict(Counter(p['boundary_projection_status'] for p in proj))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
