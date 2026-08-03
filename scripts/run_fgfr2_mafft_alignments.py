#!/usr/bin/env python3
"""
Run FGFR2 MAFFT alignments.

Run MAFFT (--auto) on the four prepared MSA inputs and report alignment-quality
statistics. MAFFT is REQUIRED for final analysis: if it is missing the module fails
with a clear message and writes metadata/msa_dependency_check_failed.txt. A trivial
pad-to-length fallback exists ONLY for smoke tests and is always labelled
not_for_final_analysis (never used unless --allow_fallback is given).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402


ALIGN_SPECS = [
    ("full_length_protein", "fgfr2_full_length_protein_msa_input.faa",
     "fgfr2_full_length_protein_msa.aln.faa"),
    ("IIIb_cassette", "fgfr2_IIIb_cassette_msa_input.faa",
     "fgfr2_IIIb_cassette_msa.aln.faa"),
    ("IIIc_cassette", "fgfr2_IIIc_cassette_msa_input.faa",
     "fgfr2_IIIc_cassette_msa.aln.faa"),
    ("IIIb_IIIc_combined_cassette", "fgfr2_IIIb_IIIc_combined_cassette_msa_input.faa",
     "fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa"),
]
MANIFEST_COLS = [
    "msa_name", "input_file", "alignment_file", "aligner", "aligner_mode", "for_final_analysis",
    "n_sequences", "ungapped_min_length", "ungapped_max_length", "aligned_length",
    "mean_gap_fraction", "max_gap_fraction_sequence", "max_gap_fraction_value",
    "alignment_status", "alignment_warning",
]


def gap_fraction(seq: str) -> float:
    if not seq:
        return 1.0
    g = sum(1 for c in seq if c in M.GAP_CHARS)
    return g / len(seq)


def mafft_version(mafft_bin: str) -> str:
    try:
        p = subprocess.run([mafft_bin, "--version"], capture_output=True, text=True, timeout=30)
        return (p.stderr or p.stdout or "").strip().splitlines()[0] if (p.stderr or p.stdout) else "unknown"
    except Exception:  # noqa: BLE001
        return "unknown"


def run_mafft(mafft_bin: str, inp: Path, out: Path, timeout: int) -> Tuple[bool, str]:
    try:
        with open(out, "w", encoding="utf-8") as fh:
            p = subprocess.run([mafft_bin, "--auto", "--quiet", str(inp)],
                               stdout=fh, stderr=subprocess.PIPE, text=True, timeout=timeout)
        if p.returncode != 0:
            return False, (p.stderr or "mafft non-zero exit").strip()[:300]
        return True, ""
    except subprocess.TimeoutExpired:
        return False, "mafft timeout"
    except Exception as e:  # noqa: BLE001
        return False, f"{type(e).__name__}: {e}"[:300]


def fallback_align(inp: Path, out: Path) -> None:
    """Smoke-test only: right-pad sequences with gaps to equal length. NOT an alignment."""
    items = M.read_fasta(inp)
    maxlen = max((len(M.clean_alignment_seq(s)) for _, s in items), default=0)
    padded = [(i, M.clean_alignment_seq(s) + "-" * (maxlen - len(M.clean_alignment_seq(s))))
              for i, s in items]
    M.write_fasta(out, padded)


def stats_for(aln: Path) -> Dict[str, object]:
    items = [(i, M.clean_alignment_seq(s)) for i, s in M.read_fasta(aln)]
    if not items:
        return {"n_sequences": 0, "ungapped_min_length": "", "ungapped_max_length": "",
                "aligned_length": 0, "mean_gap_fraction": "", "max_gap_fraction_sequence": "",
                "max_gap_fraction_value": ""}
    aligned_len = max(len(s) for _, s in items)
    ungapped_lens = [len(M.ungapped(s)) for _, s in items]
    gfs = [(i, gap_fraction(s)) for i, s in items]
    mean_gf = sum(g for _, g in gfs) / len(gfs)
    max_id, max_gf = max(gfs, key=lambda x: x[1])
    return {
        "n_sequences": len(items), "ungapped_min_length": min(ungapped_lens),
        "ungapped_max_length": max(ungapped_lens), "aligned_length": aligned_len,
        "mean_gap_fraction": round(mean_gf, 4), "max_gap_fraction_sequence": max_id,
        "max_gap_fraction_value": round(max_gf, 4),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Run MAFFT alignments.")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--allow_fallback", action="store_true",
                    help="smoke-test only: pad-to-length pseudo-alignment if MAFFT missing")
    ap.add_argument("--timeout", type=int, default=600)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    aln_dir, meta = dirs["alignments"], dirs["metadata"]
    inp = dirs["inputs"]

    mafft_bin = shutil.which("mafft")
    use_fallback = False
    if not mafft_bin:
        msg = ("MAFFT not found on PATH. MAFFT is required for final MSA analysis.\n"
               "Install e.g. `conda install -c bioconda mafft` or `brew install mafft`,\n"
               "then re-run: python scripts/run_fgfr2_msa_boundary_module.py --base <base>\n")
        (meta / "msa_dependency_check_failed.txt").write_text(
            f"{M.now_iso()}\n{msg}", encoding="utf-8")
        if not args.allow_fallback:
            print("[FAIL] " + msg, file=sys.stderr)
            return 3
        use_fallback = True
        print("[WARN] MAFFT missing; using pad-to-length fallback (not_for_final_analysis).",
              file=sys.stderr)

    aligner = "mafft_fallback_pad" if use_fallback else "mafft"
    aligner_mode = "pad_to_length_NOT_AN_ALIGNMENT" if use_fallback else "--auto"
    ver = "fallback" if use_fallback else mafft_version(mafft_bin)

    rows: List[Dict[str, object]] = []
    for name, in_name, out_name in ALIGN_SPECS:
        in_path = inp / in_name
        out_path = aln_dir / out_name
        n_in = len(M.read_fasta(in_path)) if in_path.exists() else 0
        if not in_path.exists() or n_in == 0:
            rows.append({"msa_name": name, "input_file": in_name, "alignment_file": out_name,
                         "aligner": aligner, "aligner_mode": aligner_mode,
                         "for_final_analysis": "false" if use_fallback else "true",
                         "n_sequences": 0, "ungapped_min_length": "", "ungapped_max_length": "",
                         "aligned_length": 0, "mean_gap_fraction": "",
                         "max_gap_fraction_sequence": "", "max_gap_fraction_value": "",
                         "alignment_status": "msa_failed_missing_input",
                         "alignment_warning": f"input missing or empty ({in_name})"})
            continue
        if use_fallback:
            fallback_align(in_path, out_path)
            ok, err = True, ""
        else:
            ok, err = run_mafft(mafft_bin, in_path, out_path, args.timeout)
        if not ok:
            rows.append({"msa_name": name, "input_file": in_name, "alignment_file": out_name,
                         "aligner": aligner, "aligner_mode": aligner_mode,
                         "for_final_analysis": "true", "n_sequences": n_in,
                         "ungapped_min_length": "", "ungapped_max_length": "",
                         "aligned_length": 0, "mean_gap_fraction": "",
                         "max_gap_fraction_sequence": "", "max_gap_fraction_value": "",
                         "alignment_status": "msa_failed_mafft_error",
                         "alignment_warning": err})
            continue
        st = stats_for(out_path)
        # mark review-sequence presence (gap-rich) without failing
        warn, status = "", ("msa_pass" if not use_fallback else "msa_pass_with_review_sequences")
        if isinstance(st["max_gap_fraction_value"], float) and st["max_gap_fraction_value"] > 0.5:
            status = "msa_pass_with_review_sequences"
            warn = f"gap-rich sequence present (max_gap_fraction={st['max_gap_fraction_value']})"
        if use_fallback:
            warn = (warn + ";" if warn else "") + "not_for_final_analysis_fallback_used"
        rows.append({"msa_name": name, "input_file": in_name, "alignment_file": out_name,
                     "aligner": aligner, "aligner_mode": aligner_mode,
                     "for_final_analysis": "false" if use_fallback else "true",
                     "alignment_status": status, "alignment_warning": warn, **st})

    M.write_tsv(meta / "msa_run_manifest.tsv", rows, MANIFEST_COLS)
    M.write_tsv(meta / "msa_dependency_versions.tsv", [
        {"tool": "mafft", "version": ver, "path": mafft_bin or "", "required": "true"},
        {"tool": "python", "version": sys.version.split()[0], "path": sys.executable, "required": "true"},
    ], ["tool", "version", "path", "required"])
    n_pass = sum(1 for r in rows if str(r["alignment_status"]).startswith("msa_pass"))
    print(f"[OK] alignments: {n_pass}/{len(rows)} passed (aligner={aligner}, version={ver})")
    for r in rows:
        print(f"     {r['msa_name']}: {r['alignment_status']} n={r['n_sequences']} "
              f"len={r['aligned_length']} mean_gap={r['mean_gap_fraction']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
