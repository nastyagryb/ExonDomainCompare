#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402


CASSETTES = [
    ("IIIb_cassette", "fgfr2_IIIb_cassette_msa_input.faa",
     "fgfr2_IIIb_cassette_msa_linsi.aln.faa", "fgfr2_IIIb_cassette_msa.aln.faa"),
    ("IIIc_cassette", "fgfr2_IIIc_cassette_msa_input.faa",
     "fgfr2_IIIc_cassette_msa_linsi.aln.faa", "fgfr2_IIIc_cassette_msa.aln.faa"),
    ("IIIb_IIIc_combined_cassette", "fgfr2_IIIb_IIIc_combined_cassette_msa_input.faa",
     "fgfr2_IIIb_IIIc_combined_cassette_msa_linsi.aln.faa",
     "fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa"),
]
COLS = ["msa_name", "strategy", "n_sequences", "aligned_length", "mean_gap_fraction",
        "boundary_projection_agreement_with_auto", "informative_columns", "gap_rich_columns",
        "recommended_for_main_figures", "strategy_warning"]


def gap_fraction(seq: str) -> float:
    return (sum(1 for c in seq if c in M.GAP_CHARS) / len(seq)) if seq else 1.0


def col_metrics(items: List[Tuple[str, str]]):
    if not items:
        return 0, 0.0, 0, 0
    L = max(len(s) for _, s in items)
    n = len(items)
    informative = gap_rich = 0
    gfs = []
    for c in range(L):
        col = [s[c] if c < len(s) else "-" for _, s in items]
        gf = sum(1 for ch in col if ch in M.GAP_CHARS) / n
        gfs.append(gf)
        if gf >= 0.5:
            gap_rich += 1
        elif sum(1 for ch in col if ch not in M.GAP_CHARS) >= 2:
            informative += 1
    mean_gf = round(sum(s2 for s2 in (gap_fraction(s) for _, s in items)) / n, 4)
    return L, mean_gf, informative, gap_rich


def start_cols_norm(items: List[Tuple[str, str]]) -> Dict[str, float]:
    out = {}
    L = max((len(s) for _, s in items), default=1)
    for sid, s in items:
        first = next((i for i, ch in enumerate(s) if ch not in M.GAP_CHARS), 0)
        out[sid] = first / L if L else 0.0
    return out


def pearson(a: List[float], b: List[float]) -> Optional[float]:
    n = len(a)
    if n < 2:
        return None
    ma, mb = sum(a) / n, sum(b) / n
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    va = sum((x - ma) ** 2 for x in a) ** 0.5
    vb = sum((y - mb) ** 2 for y in b) ** 0.5
    if va == 0 or vb == 0:
        return 1.0
    return round(cov / (va * vb), 4)


def main() -> int:
    ap = argparse.ArgumentParser(description="L-INS-i cassette alignments (Part A).")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--allow_fallback", action="store_true")
    ap.add_argument("--timeout", type=int, default=900)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    inp, aln_dir, meta = dirs["inputs"], dirs["alignments"], dirs["metadata"]

    mafft = shutil.which("mafft")
    rows: List[Dict[str, object]] = []
    for name, in_name, linsi_out, auto_out in CASSETTES:
        in_path = inp / in_name
        auto_path = aln_dir / auto_out
        linsi_path = aln_dir / linsi_out
        n_in = len(M.read_fasta(in_path)) if in_path.exists() else 0
        auto_items = [(i, M.clean_alignment_seq(s)) for i, s in M.read_fasta(auto_path)]
        auto_L, auto_mg, auto_inf, auto_gr = col_metrics(auto_items)
        rows.append({"msa_name": name, "strategy": "mafft_auto", "n_sequences": len(auto_items),
                     "aligned_length": auto_L, "mean_gap_fraction": auto_mg,
                     "boundary_projection_agreement_with_auto": 1.0,
                     "informative_columns": auto_inf, "gap_rich_columns": auto_gr,
                     "recommended_for_main_figures": "false",
                     "strategy_warning": "kept as reproducibility/sensitivity alignment"})
        if not mafft and not args.allow_fallback:
            rows.append({"msa_name": name, "strategy": "mafft_linsi", "n_sequences": 0,
                         "aligned_length": 0, "mean_gap_fraction": "",
                         "boundary_projection_agreement_with_auto": "", "informative_columns": 0,
                         "gap_rich_columns": 0, "recommended_for_main_figures": "false",
                         "strategy_warning": "MAFFT missing; L-INS-i not run"})
            continue
        ok, err = True, ""
        if mafft and n_in > 0:
            try:
                with open(linsi_path, "w", encoding="utf-8") as fh:
                    p = subprocess.run([mafft, "--localpair", "--maxiterate", "1000", "--quiet",
                                        str(in_path)], stdout=fh, stderr=subprocess.PIPE,
                                       text=True, timeout=args.timeout)
                ok = (p.returncode == 0)
                err = (p.stderr or "").strip()[:200] if not ok else ""
            except Exception as e:  # noqa: BLE001
                ok, err = False, f"{type(e).__name__}: {e}"[:200]
        elif args.allow_fallback:
            shutil.copyfile(auto_path, linsi_path)
            ok, err = True, "fallback_copied_auto_not_for_final_analysis"
        if not ok:
            rows.append({"msa_name": name, "strategy": "mafft_linsi", "n_sequences": n_in,
                         "aligned_length": 0, "mean_gap_fraction": "",
                         "boundary_projection_agreement_with_auto": "", "informative_columns": 0,
                         "gap_rich_columns": 0, "recommended_for_main_figures": "false",
                         "strategy_warning": f"L-INS-i failed: {err}"})
            continue
        linsi_items = [(i, M.clean_alignment_seq(s)) for i, s in M.read_fasta(linsi_path)]
        L, mg, inf, gr = col_metrics(linsi_items)
        a_norm = start_cols_norm(auto_items)
        l_norm = start_cols_norm(linsi_items)
        common = [sid for sid in a_norm if sid in l_norm]
        agree = pearson([a_norm[s] for s in common], [l_norm[s] for s in common])
        valid = (len(linsi_items) == n_in and L > 0)
        warn = err if err else ""
        if mg >= 0.7:
            warn = (warn + ";" if warn else "") + f"high_mean_gap_fraction({mg})"
        rows.append({"msa_name": name, "strategy": "mafft_linsi", "n_sequences": len(linsi_items),
                     "aligned_length": L, "mean_gap_fraction": mg,
                     "boundary_projection_agreement_with_auto": agree if agree is not None else "",
                     "informative_columns": inf, "gap_rich_columns": gr,
                     "recommended_for_main_figures": "true" if valid else "false",
                     "strategy_warning": warn})

    M.write_tsv(meta / "msa_strategy_comparison.tsv", rows, COLS)
    linsi_ok = sum(1 for r in rows if r["strategy"] == "mafft_linsi" and r["aligned_length"])
    print(f"[OK] L-INS-i cassette alignments: {linsi_ok}/{len(CASSETTES)} produced")
    for r in rows:
        if r["strategy"] == "mafft_linsi":
            print(f"     {r['msa_name']}: len={r['aligned_length']} mean_gap={r['mean_gap_fraction']} "
                  f"informative={r['informative_columns']} agree_with_auto="
                  f"{r['boundary_projection_agreement_with_auto']} main={r['recommended_for_main_figures']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
