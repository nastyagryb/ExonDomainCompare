#!/usr/bin/env python3
"""
run_fgfr2_msa_boundary_module.py  (MSA boundary-robustness sprint, Parts 1 + 12 + orchestration)

Central runner for the pre-InterPro MSA / boundary-robustness module. Performs the MAFFT
dependency check, orchestrates every step (input prep -> MAFFT -> coordinate maps/projection
-> conservation -> discriminating residues -> protein integrity -> splice QC -> robustness
-> figures -> reports), integrates summaries into species_qc_master.tsv (Part 12) and writes
the reproducibility manifests (Part 1).

This module is strictly pre-InterProScan: it never runs InterProScan and never invents
domain annotations. MSA is an independent robustness/QC layer that does not relabel IIIb/IIIc.

Acceptance command:
    python scripts/run_fgfr2_msa_boundary_module.py --base results/final_30_until_interpro_prepare
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fgfr2_msa_common as M  # noqa: E402

SCRIPT_VERSION = "1.0"
SCRIPTS_DIR = Path(__file__).resolve().parent

# data-generating + validation steps (run before figures)
DATA_STEPS = [
    ("reconcile_labels", "reconcile_fgfr2_exon_type_labels.py"),
    ("prepare_inputs", "prepare_fgfr2_msa_inputs.py"),
    ("correct_cassette_windows", "correct_fgfr2_cassette_windows.py"),
    ("mafft", "run_fgfr2_mafft_alignments.py"),
    ("cassette_linsi", "run_fgfr2_cassette_linsi_alignments.py"),
    ("project_boundaries", "project_fgfr2_boundaries_to_msa.py"),
    ("conservation", "analyze_fgfr2_msa_conservation.py"),
    ("discriminating", "analyze_fgfr2_isoform_discriminating_residues.py"),
    ("reference_agreement", "build_fgfr2_reference_agreement.py"),
    ("protein_integrity", "check_fgfr2_pre_interpro_protein_integrity.py"),
    ("splice_qc", "check_fgfr2_splice_site_qc.py"),
    ("robustness", "score_fgfr2_boundary_robustness.py"),
    ("validate_rescue", "validate_and_rescue_fgfr2_labels.py"),
]
# steps to re-run on pass 2 (consume rescue overrides); reconcile is NOT reset
RERUN_STEPS = [s for s in DATA_STEPS if s[0] != "reconcile_labels"]
MAXRESCUE_SCRIPT = "maximal_rescue_fgfr2_labels.py"
FIGURE_STEPS = [
    ("figures", "make_fgfr2_msa_boundary_figures.py"),
    ("reference_figures", "make_fgfr2_msa_reference_figures.py"),
]
# local synteny / gene-neighborhood validation runs ONCE after the final post-rescue state
# (it consumes the post-rescue truth table); MCScanX block-level synteny is intentionally omitted.
SYNTENY_STEPS = [
    ("synteny", "validate_fgfr2_local_synteny_neighborhood.py"),
    ("synteny_figures", "make_fgfr2_synteny_figures.py"),
]
STEP_SCRIPTS = (DATA_STEPS + [("maximal_rescue", MAXRESCUE_SCRIPT)] + FIGURE_STEPS + SYNTENY_STEPS)

EXPECTED_INPUTS = [
    "species_qc_master.tsv", "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv",
    "fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv", "fgfr2_pair_level_qc_summary.tsv",
    "fgfr2_III_pair_audit.tsv", "fgfr2_refined_uncertainty_classes.tsv",
    "fgfr2_interpro_clean_unique.fasta", "fgfr2_interpro_unique_mapping.tsv",
    "fgfr2_interpro_id_mapping.tsv", "selected_fgfr2_proteins.faa",
    "fgfr2_cassette_cds_block_map.tsv", "fgfr2_unique_cds_block_table.tsv",
    "cds_phase_boundary_audit.tsv", "fgfr2_transcript_cds_reconstruction_audit.tsv",
    "species_phylogenetic_order.tsv",
]

RCLASS_RANK = {"unresolved_or_annotation_dependent_boundary": 4, "review_boundary": 3,
               "supported_boundary_with_minor_flags": 2, "robust_boundary": 1}
INTEG_RANK = {"missing_sequence_fail": 5, "invalid_sequence_review": 4,
              "protein_length_outlier_review": 3, "protein_integrity_pass_with_minor_warning": 2,
              "protein_integrity_pass": 1}


def run_step(script: str, base: Path, extra: Optional[List[str]] = None) -> int:
    cmd = [sys.executable, str(SCRIPTS_DIR / script), "--base", str(base)] + (extra or [])
    print(f"\n>>> {script}", flush=True)
    p = subprocess.run(cmd)
    return p.returncode


def propagate_reconciliation(base: Path, dirs: Dict[str, Path]) -> None:
    """Add legacy-preserving reconciliation columns (legacy/upstream/final/validated/consistency)
    to the major upstream audit tables and species_qc_master.tsv. Existing columns are preserved;
    the upstream isoform column is never overwritten — final_isoform_label is added alongside."""
    rec = M.read_tsv(dirs["maps"] / "fgfr2_exon_type_label_reconciliation.tsv")
    if not rec:
        print("[WARN] reconciliation table missing; skipping propagation.", file=sys.stderr)
        return
    by_pair = {(r["species"].lower(), r["upstream_label"]): r for r in rec}
    add_cols = ["legacy_label", "upstream_label", "previous_pipeline_label", "final_isoform_label",
                "validated_exon_type", "label_consistency_status", "label_reconciliation_action",
                "validation_group", "rescue_status", "final_claim_status", "final_label_source",
                "final_claim_status_after_rescue", "maximal_rescue_decision"]

    def patch(path: Path, sp_key: str, iso_key: str) -> None:
        if not path or not path.exists():
            return
        rows = M.read_tsv(path)
        if not rows:
            return
        fields = list(rows[0].keys())
        for c in add_cols:
            if c not in fields:
                fields.append(c)
        for row in rows:
            sp = (row.get(sp_key) or "").lower()
            up = row.get(iso_key) or ""
            r = by_pair.get((sp, up), {})
            row["legacy_label"] = up
            row["upstream_label"] = up
            row["previous_pipeline_label"] = up
            row["final_isoform_label"] = r.get("final_isoform_label", up)
            row["validated_exon_type"] = r.get("validated_exon_type", up)
            row["label_consistency_status"] = r.get("label_consistency_status", "no_reconciliation")
            row["label_reconciliation_action"] = r.get("label_reconciliation_action", "keep_upstream_label")
            row["validation_group"] = r.get("validation_group", "")
            row["rescue_status"] = r.get("rescue_status", "")
            row["final_claim_status"] = r.get("final_claim_status", "")
            row["final_label_source"] = r.get("final_label_source", "")
            row["final_claim_status_after_rescue"] = r.get("final_claim_status_after_rescue",
                                                           r.get("final_claim_status", ""))
            row["maximal_rescue_decision"] = r.get("maximal_rescue_decision", "")
        M.write_tsv(path, rows, fields)

    for name, sp_key, iso_key in [
            ("fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv", "species_canonical", "inferred_isoform"),
            ("fgfr2_cassette_cds_block_map.tsv", "species", "isoform"),
            ("fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv", "species_canonical", "inferred_isoform")]:
        p = M.locate(base, name)
        patch(p, sp_key, iso_key)

    # species_qc_master: per-species joined reconciliation summary (legacy preserved)
    master_p = M.require(base, "species_qc_master.tsv", "11_pre_interpro_master")
    master = M.read_tsv(master_p)
    fields = list(master[0].keys()) if master else []
    msum = ["upstream_isoform_labels", "final_isoform_labels", "label_consistency_status",
            "label_reconciliation_action", "validation_group", "final_claim_status",
            "rescue_status_summary", "final_claim_status_after_rescue", "maximal_rescue_decision_summary",
            "recommended_use_pre_rescue", "recommended_use_post_rescue"]

    def _post_use(claim_set):
        if any(c.startswith("primary_claim") for c in claim_set):
            return "main_analysis" if all(c.startswith("primary_claim") for c in claim_set) \
                else "main_analysis_partial"
        if "supplement_review" in claim_set:
            return "supplement_only"
        if "excluded_from_primary_claim" in claim_set:
            return "exclude_from_primary_claim"
        return "review"
    for c in msum:
        if c not in fields:
            fields.append(c)
    by_sp: Dict[str, List[Dict[str, str]]] = {}
    for r in rec:
        by_sp.setdefault(r["species"].lower(), []).append(r)
    for row in master:
        rs = by_sp.get(row["species"].lower(), [])
        row["upstream_isoform_labels"] = ";".join(f"{r['upstream_label']}={r['upstream_label']}" for r in rs)
        row["final_isoform_labels"] = ";".join(f"{r['upstream_label']}->{r['final_isoform_label']}" for r in rs)
        row["label_consistency_status"] = ";".join(sorted({r["label_consistency_status"] for r in rs})) or "no_reconciliation"
        row["label_reconciliation_action"] = ";".join(sorted({r["label_reconciliation_action"] for r in rs})) or "keep_upstream_label"
        row["validation_group"] = ";".join(sorted({r.get("validation_group", "") for r in rs} - {""}))
        row["final_claim_status"] = ";".join(sorted({r.get("final_claim_status", "") for r in rs} - {""}))
        row["rescue_status_summary"] = ";".join(sorted({r.get("rescue_status", "") for r in rs} - {""}))
        row["final_claim_status_after_rescue"] = ";".join(sorted(
            {r.get("final_claim_status_after_rescue", "") for r in rs} - {""})) or \
            row["final_claim_status"]
        row["maximal_rescue_decision_summary"] = ";".join(sorted(
            {r.get("maximal_rescue_decision", "") for r in rs} - {""}))
        claim_set = {r.get("final_claim_status_after_rescue") or r.get("final_claim_status", "")
                     for r in rs} - {""}
        row["recommended_use_pre_rescue"] = row.get("recommended_use", "")
        row["recommended_use_post_rescue"] = _post_use(claim_set)
    M.write_tsv(master_p, master, fields)
    n_swap = sum(1 for r in rec if r["label_consistency_status"] == "swapped_relative_to_upstream")
    print(f"[OK] reconciliation propagated (legacy preserved); {n_swap} swapped cassettes flagged")


def update_master(base: Path, dirs: Dict[str, Path]) -> None:
    master_p = M.require(base, "species_qc_master.tsv", "11_pre_interpro_master")
    master = M.read_tsv(master_p)
    fields = list(master[0].keys()) if master else []

    run_man = {r["msa_name"]: r for r in M.read_tsv(dirs["metadata"] / "msa_run_manifest.tsv")}
    full_status = run_man.get("full_length_protein", {}).get("alignment_status", "")
    cass_status = ";".join(sorted({run_man.get(n, {}).get("alignment_status", "")
                                   for n in ("IIIb_cassette", "IIIc_cassette") if n in run_man}))
    proj = M.read_tsv(dirs["maps"] / "fgfr2_exon_boundary_msa_projection.tsv")
    scores = M.read_tsv(dirs["robustness"] / "fgfr2_boundary_robustness_scores.tsv")
    integ = M.read_tsv(dirs["protein_integrity"] / "fgfr2_pre_interpro_protein_integrity_qc.tsv")
    splice = M.read_tsv(dirs["splice_qc"] / "fgfr2_splice_site_boundary_qc.tsv")
    diag = {r["species"].lower() for r in
            M.read_tsv(dirs["review_diagnostics"] / "fgfr2_msa_review_case_diagnostics.tsv")}

    by_sp: Dict[str, Dict[str, list]] = {}
    for r in proj:
        by_sp.setdefault(r["species"].lower(), {}).setdefault("proj", []).append(r)
    for r in scores:
        by_sp.setdefault(r["species"].lower(), {}).setdefault("score", []).append(r)
    for r in integ:
        by_sp.setdefault(r["species"].lower(), {}).setdefault("integ", []).append(r)
    for r in splice:
        by_sp.setdefault(r["species"].lower(), {}).setdefault("splice", []).append(r)

    new_cols = ["full_length_msa_status", "cassette_msa_status", "msa_boundary_projection_summary",
                "msa_gap_review_count", "cassette_conservation_summary", "boundary_robustness_summary",
                "boundary_robustness_min_score", "boundary_robustness_class", "splice_site_qc_summary",
                "protein_integrity_summary", "review_case_msa_diagnostic_status"]
    for c in new_cols:
        if c not in fields:
            fields.append(c)

    for row in master:
        sp = row["species"].lower()
        d = by_sp.get(sp, {})
        projs = d.get("proj", [])
        scs = d.get("score", [])
        igs = d.get("integ", [])
        sps = d.get("splice", [])
        proj_sum = ";".join(sorted({p["boundary_projection_status"] for p in projs})) or "no_data"
        gap_review = sum(1 for p in projs if "review" in p["boundary_projection_status"])
        cass_cons_vals = [M.to_float(s.get("cassette_conservation_score")) for s in scs]
        cass_cons_vals = [v for v in cass_cons_vals if v is not None]
        cass_cons_sum = (f"mean={round(sum(cass_cons_vals)/len(cass_cons_vals),3)}"
                         if cass_cons_vals else "unavailable")
        rclasses = [s["boundary_robustness_class"] for s in scs]
        worst_class = max(rclasses, key=lambda x: RCLASS_RANK.get(x, 0)) if rclasses else "no_data"
        min_score = min((M.to_float(s.get("boundary_robustness_score"), 1.0) for s in scs), default="")
        splice_sum = ";".join(sorted({s["splice_site_qc_status"] for s in sps})) or "no_data"
        integ_worst = (max((i["protein_integrity_status"] for i in igs),
                           key=lambda x: INTEG_RANK.get(x, 0)) if igs else "no_data")
        row["full_length_msa_status"] = full_status
        row["cassette_msa_status"] = cass_status
        row["msa_boundary_projection_summary"] = proj_sum
        row["msa_gap_review_count"] = gap_review
        row["cassette_conservation_summary"] = cass_cons_sum
        row["boundary_robustness_summary"] = ";".join(sorted(set(rclasses))) or "no_data"
        row["boundary_robustness_min_score"] = min_score if min_score != "" else ""
        row["boundary_robustness_class"] = worst_class
        row["splice_site_qc_summary"] = splice_sum
        row["protein_integrity_summary"] = integ_worst
        row["review_case_msa_diagnostic_status"] = (
            "review_case_present" if sp in diag else "no_msa_review_flag")

    M.write_tsv(master_p, master, fields)
    print(f"[OK] species_qc_master.tsv updated with {len(new_cols)} MSA/robustness columns")


def write_manifests(base: Path, dirs: Dict[str, Path], cmd_used: str,
                    step_status: Dict[str, int], mafft_ver: str) -> None:
    meta = dirs["metadata"]
    # input file manifest + hashes
    in_rows = []
    for name in EXPECTED_INPUTS:
        p = M.locate(base, name)
        in_rows.append({"input_name": name, "found": "true" if p else "false",
                        "path": str(p.relative_to(base)) if p else "",
                        "sha256": M.sha256_file(p) if p else "",
                        "size_bytes": (p.stat().st_size if p else "")})
    M.write_tsv(meta / "msa_input_file_manifest.tsv", in_rows,
                ["input_name", "found", "path", "sha256", "size_bytes"])
    # output file manifest + hashes
    md = M.module_dir(base)
    out_rows = []
    for sd in ["inputs", "alignments", "maps", "conservation", "robustness", "splice_qc",
               "protein_integrity", "review_diagnostics", "figures", "tables", "captions"]:
        for f in sorted((md / sd).glob("*")):
            if f.is_file():
                out_rows.append({"subdir": sd, "output_name": f.name,
                                 "path": str(f.relative_to(base)), "sha256": M.sha256_file(f),
                                 "size_bytes": f.stat().st_size})
    M.write_tsv(meta / "msa_output_file_manifest.tsv", out_rows,
                ["subdir", "output_name", "path", "sha256", "size_bytes"])
    # script versions
    sv = []
    for _, scr in STEP_SCRIPTS + [("runner", "run_fgfr2_msa_boundary_module.py"),
                                  ("common", "_fgfr2_msa_common.py")]:
        sp = SCRIPTS_DIR / scr
        sv.append({"script": scr, "sha256": M.sha256_file(sp) if sp.exists() else "",
                   "exists": "true" if sp.exists() else "false"})
    M.write_tsv(meta / "msa_script_versions.tsv", sv, ["script", "sha256", "exists"])
    # dependency versions (mafft + python) — refreshed here too
    M.write_tsv(meta / "msa_dependency_versions.tsv", [
        {"tool": "mafft", "version": mafft_ver, "path": shutil.which("mafft") or "", "required": "true"},
        {"tool": "python", "version": sys.version.split()[0], "path": sys.executable, "required": "true"},
        {"tool": "matplotlib", "version": _try_ver("matplotlib"), "path": "", "required": "true"},
        {"tool": "numpy", "version": _try_ver("numpy"), "path": "", "required": "false"},
    ], ["tool", "version", "path", "required"])
    # run manifest json
    manifest = {
        "module": "12_msa_boundary_robustness_pre_interpro",
        "script_version": SCRIPT_VERSION, "run_timestamp_utc": M.now_iso(),
        "command_used": cmd_used, "python_version": sys.version.split()[0],
        "mafft_version": mafft_ver, "base": str(base),
        "step_status": step_status, "pre_interproscan": True,
        "msa_relabels_isoforms": False, "fake_interpro_domains": False,
    }
    (meta / "msa_module_run_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[OK] manifests written ({len(in_rows)} inputs, {len(out_rows)} outputs)")


def _try_ver(mod: str) -> str:
    try:
        m = __import__(mod)
        return getattr(m, "__version__", "unknown")
    except Exception:  # noqa: BLE001
        return "unavailable"


def write_reports(base: Path, dirs: Dict[str, Path]) -> None:
    cmd = [sys.executable, str(SCRIPTS_DIR / "make_fgfr2_msa_boundary_reports.py"),
           "--base", str(base)]
    if (SCRIPTS_DIR / "make_fgfr2_msa_boundary_reports.py").exists():
        subprocess.run(cmd)


def main() -> int:
    ap = argparse.ArgumentParser(description="Run the pre-InterPro MSA boundary module.")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--allow_fallback", action="store_true",
                    help="smoke-test only: allow pad-to-length fallback if MAFFT missing")
    ap.add_argument("--genome_dir", type=Path, default=None,
                    help="optional genome FASTA dir for splice-site QC")
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    cmd_used = "python " + " ".join([str(Path(sys.argv[0]).name)] + sys.argv[1:])

    # ---- dependency check (Part 1) ----
    mafft_bin = shutil.which("mafft")
    if not mafft_bin and not args.allow_fallback:
        msg = ("MAFFT not found on PATH. MAFFT is required for final MSA analysis.\n"
               "Install via `conda install -c bioconda mafft` or `brew install mafft` and re-run.\n"
               "Re-run with --allow_fallback ONLY for smoke tests (not_for_final_analysis).\n")
        (dirs["metadata"] / "msa_dependency_check_failed.txt").write_text(
            f"{M.now_iso()}\n{msg}", encoding="utf-8")
        print("[FAIL] " + msg, file=sys.stderr)
        return 3
    mafft_ver = "missing"
    if mafft_bin:
        try:
            p = subprocess.run([mafft_bin, "--version"], capture_output=True, text=True, timeout=30)
            mafft_ver = (p.stderr or p.stdout or "").strip().splitlines()[0]
        except Exception:  # noqa: BLE001
            mafft_ver = "unknown"

    def _extra(step: str) -> List[str]:
        if step in ("mafft", "cassette_linsi") and args.allow_fallback:
            return ["--allow_fallback"]
        if step == "splice_qc" and args.genome_dir:
            return ["--genome_dir", str(args.genome_dir)]
        return []

    step_status: Dict[str, int] = {}
    # ---- pass 1: data + validation steps ----
    for step, scr in DATA_STEPS:
        rc = run_step(scr, base, _extra(step))
        step_status[step] = rc
        if rc != 0 and step in ("reconcile_labels", "prepare_inputs", "mafft",
                                "project_boundaries", "validate_rescue"):
            print(f"[FAIL] critical step '{step}' failed (rc={rc}); aborting.", file=sys.stderr)
            write_manifests(base, dirs, cmd_used, step_status, mafft_ver)
            return rc

    # ---- maximal rescue (pass 1): emit sequence/provenance-validated overrides ----
    # gate is NOT enforced here (corrected data not yet propagated); --no_control_gate
    step_status["maximal_rescue"] = run_step(MAXRESCUE_SCRIPT, base, ["--no_control_gate"])

    # ---- pass 2: if rescue produced overrides, regenerate downstream once (reconcile NOT reset) ----
    ovr_flag = dirs["maps"] / "maximal_rescue_overrides.flag"
    sel_flag = dirs["maps"] / "rescue_selection_changed.flag"
    if ovr_flag.exists() or sel_flag.exists():
        for f in (ovr_flag, sel_flag):
            if f.exists():
                f.unlink()
        print("\n[INFO] rescue produced validated overrides; regenerating downstream once "
              "(consuming corrected candidates).")
        for step, scr in RERUN_STEPS:
            step_status[step + "_rerun"] = run_step(scr, base, _extra(step))

    # ---- maximal rescue (final pass): enforce hard validation gate on corrected data ----
    gate_rc = run_step(MAXRESCUE_SCRIPT, base, ["--final_pass"])
    step_status["maximal_rescue_final"] = gate_rc

    propagate_reconciliation(base, dirs)  # legacy-preserving final-label propagation
    update_master(base, dirs)           # Part 12

    if gate_rc == 4:
        print("[FAIL] maximal-rescue validation gate hard-failed; NOT generating primary figures.",
              file=sys.stderr)
        write_reports(base, dirs)
        write_manifests(base, dirs, cmd_used, step_status, mafft_ver)
        return 4

    # ---- local synteny / gene-neighborhood validation (once, on final post-rescue state) ----
    step_status["synteny"] = run_step(SYNTENY_STEPS[0][1], base, [])

    # Clear STALE cross-checked figure tables before the figure phase. The per-script
    # consistency gate cross-validates figure6 / figure8 / figure6C, but these tables are
    # produced by DIFFERENT figure scripts (figure6/figure8 by make_..._boundary_figures,
    # figure6C by make_..._reference_figures). On a re-run, a stale figure6C from a previous
    # rescue state would make the first script's gate fail spuriously. Removing them first
    # makes the gate treat not-yet-regenerated tables as absent (skipped) rather than stale,
    # so each script regenerates its own tables and the gate is order-independent.
    for _stale in ("figure6_msa_projected_boundary_map.tsv",
                   "figure8_boundary_robustness_evidence_stack.tsv",
                   "figure6C_human_referenced_residue_agreement_map.tsv"):
        try:
            (dirs["tables"] / _stale).unlink()
        except FileNotFoundError:
            pass

    # ---- figures (gate passed) ----
    for step, scr in FIGURE_STEPS:
        step_status[step] = run_step(scr, base, [])
    step_status["synteny_figures"] = run_step(SYNTENY_STEPS[1][1], base, [])

    write_reports(base, dirs)           # Part 14 + I
    write_manifests(base, dirs, cmd_used, step_status, mafft_ver)  # Part 1

    # The plotting steps below are PRESENTATION only. The closure (module 13) consumes
    # the MSA/rescue DATA outputs (post-rescue truth table, robustness, conservation),
    # not these figures. For a small custom run without human/mouse the figure scripts'
    # positive-control gates can legitimately fail (e.g. "0 reference neighbors"); that
    # must NOT strand the whole pipeline. Figure-step failures are therefore reported as
    # warnings but do not fail the module. For the full-30 panel these figures succeed,
    # so behaviour there is unchanged.
    OPTIONAL_FIGURE_STEPS = {"figures", "reference_figures", "synteny_figures"}
    failed = {k: v for k, v in step_status.items() if v != 0}
    critical_failed = {k: v for k, v in failed.items() if k not in OPTIONAL_FIGURE_STEPS}
    print("\n==================== MSA BOUNDARY MODULE SUMMARY ====================")
    for step, scr in STEP_SCRIPTS:
        print(f"  {step:20s} rc={step_status.get(step)}")
    if failed:
        print(f"[DONE WITH WARNINGS] failed steps: {failed}")
    if critical_failed:
        print(f"[FAIL] critical (non-figure) steps failed: {critical_failed}")
        return 1
    if failed:
        print("[NOTE] only optional figure steps failed (e.g. custom run without "
              "human/mouse controls). Core MSA/rescue DATA outputs are complete; "
              "treating module as successful for the downstream closure.")
        return 0
    print("[DONE] all MSA boundary-robustness steps completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
