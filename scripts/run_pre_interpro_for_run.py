#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.contracts import portable_runtime_record, stamp_payload  # noqa: E402
from exondomaincompare.runs.legacy import LegacyRunAdapter  # noqa: E402
from exondomaincompare.runs.registry import resolve_run_record  # noqa: E402
from exondomaincompare.config import load_config  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=Path(__file__).resolve().parent.parent)
REPO = RUNTIME_CONFIG.repository_root
RUNS_ROOT = RUNTIME_CONFIG.runs_root
# read-only validated example freeze — must never be written to by this wrapper
FREEZE_DIR = (REPO / "results" / "final_30_until_interpro_prepare").resolve()
CLOSURE_RUNNER = REPO / "run_fgfr2_pipeline_current_final_pre_interpro.sh"

CASE_STUDY = "FGFR2_IIIb_IIIc"
PRIMARY_FASTA_REL = "results/13_final_pre_interpro_closure/freeze/final_pre_interpro_proteins_primary.faa"
REVIEW_FASTA_REL = "results/13_final_pre_interpro_closure/freeze/final_pre_interpro_proteins_all_review_included.faa"

# full30 reference expectation (informational only — never a hard gate here)
FULL30_EXPECTED_PRIMARY = 58
FULL30_EXPECTED_REVIEW = 60

# Essential run-local outputs that MUST exist (Steps 1-11) before v3 may be
# skipped in cached mode. All paths are relative to runs/<run_id>/results/.
REQUIRED_V3_OUTPUTS = [
    "01_species_registry/species_registry.tsv",
    "02_models/genes.tsv",
    "02_models/transcripts.tsv",
    "02_models/cds_features.tsv",
    "04_isoform_evidence_v2_3_human_calibrated/fgfr2_isoform_evidence.tsv",
    "05b_selection_with_isoforms_v2_7_marker_validated/"
    "fgfr2_III_final_selected_protein_validation_summary.tsv",
    "06_protein_export_v2_7_marker_validated/selected_fgfr2_proteins.faa",
    "07_interpro_prepare_v2_7_marker_validated/fgfr2_interpro_clean_unique.fasta",
    "07_interpro_prepare_v2_7_marker_validated/fgfr2_interpro_id_mapping.tsv",
    "09_paper_ready_qc_v2_9/fgfr2_paper_ready_species_qc.tsv",
    "09_paper_ready_qc_v2_9/figures_v2_22_final_qc_display/"
    "fgfr2_pair_level_qc_summary.tsv",
    "11_pre_interpro_master/species_qc_master_pre_interpro.tsv",
]


def has_required_v3_outputs(results_dir: Path):
    return _check_required(results_dir, REQUIRED_V3_OUTPUTS)


# MSA reuse requires its own complete post-rescue cache.
REQUIRED_MSA_OUTPUTS = [
    "12_msa_boundary_robustness_pre_interpro/maps/fgfr2_post_rescue_final_truth_table.tsv",
    "12_msa_boundary_robustness_pre_interpro/robustness/fgfr2_boundary_robustness_scores.tsv",
]


def _check_required(results_dir: Path, rel_paths: List[str]):
    found: List[str] = []
    missing: List[str] = []
    for relp in rel_paths:
        p = results_dir / relp
        try:
            ok = p.is_file() and p.stat().st_size > 0
        except OSError:
            ok = False
        (found if ok else missing).append(relp)
    return (len(missing) == 0, found, missing)


def has_required_msa_outputs(results_dir: Path):
    ok, found, missing = _check_required(results_dir, REQUIRED_MSA_OUTPUTS)
    if not ok:
        return ok, found, missing
    gate_json = (results_dir / "12_msa_boundary_robustness_pre_interpro"
                 / "maps" / "fgfr2_maximal_rescue_validation_gate.json")
    gate_ok = False
    if gate_json.is_file():
        try:
            gate_ok = not bool(json.loads(gate_json.read_text(encoding="utf-8"))
                               .get("hard_fail", True))
        except (OSError, ValueError):
            gate_ok = False
    if not gate_ok:
        missing = missing + ["12_msa_boundary_robustness_pre_interpro/maps/"
                             "fgfr2_maximal_rescue_validation_gate.json (passed gate)"]
        return False, found, missing
    return True, found, missing


# --------------------------------------------------------------------------- #
# small utilities
# --------------------------------------------------------------------------- #
def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO))
    except ValueError:
        return str(path.resolve())


def read_json(path: Path, fallback: Any = None) -> Any:
    if not path.exists():
        return fallback
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return fallback


def write_json(path: Path, data: Any) -> None:
    if isinstance(data, dict):
        run_root = path
        while run_root.parent != RUNS_ROOT and run_root != run_root.parent:
            run_root = run_root.parent
        data = portable_runtime_record(
            data, repository_root=REPO,
            run_root=run_root if run_root.parent == RUNS_ROOT else None,
        )
        run_id = str(data.get("run_id") or (
            path.parent.name if path.name == "status.json" else ""))
        data = stamp_payload(
            data,
            payload_type="status" if path.name == "status.json" else path.stem,
            run_id=run_id,
            dataset_id=run_id,
            profile=RUNTIME_CONFIG.public_identity(),
            generator="scripts/run_pre_interpro_for_run.py",
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def detect_python() -> str:
    return RUNTIME_CONFIG.local_python().selected


def count_fasta(path: Path) -> Optional[int]:
    if not path.exists():
        return None
    n = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n




# Ensembl-style species identifier (mirrors scripts/create_new_run.py).
SPECIES_ID_RE = re.compile(r"^[a-z][a-z0-9]+_[a-z0-9_]+$")


def preflight_species_list(species: List[str]) -> None:
    for name in species:
        if name != name.strip() or " " in name or "\t" in name:
            suggestion = re.sub(r"\s+", "_", name.strip().lower())
            raise SystemExit(
                f'Invalid species identifier in species_list.txt: "{name}". '
                f'Use "{suggestion}".')
        if not SPECIES_ID_RE.match(name):
            suggestion = re.sub(r"\s+", "_", name.strip().lower())
            raise SystemExit(
                f'Invalid species identifier in species_list.txt: "{name}". '
                f'Use "{suggestion}" (lowercase underscore identifier).')


# --------------------------------------------------------------------------- #
# status / config / readme updates
# --------------------------------------------------------------------------- #
def update_status(run_dir: Path, **fields: Any) -> None:
    status_path = run_dir / "status.json"
    status = read_json(status_path, {}) or {}
    status.update(fields)
    status["last_updated"] = now_iso()
    write_json(status_path, status)


def mark_running(run_dir: Path) -> None:
    update_status(run_dir, status="running", current_step="pre_interpro_pipeline",
                  pre_interpro_status="running")


# Curated human reference FGFR2 IIIb/IIIc files (control layer; NOT analysed species).
HUMAN_REFERENCE_FILES = [
    "references/fgfr2_iii_segments/human_FGFR2_IIIb_segment.fasta",
    "references/fgfr2_iii_segments/human_FGFR2_IIIc_segment.fasta",
    "reference/human_FGFR2_IIIb_protein.faa",
    "reference/human_FGFR2_IIIc_protein.faa",
]


def record_human_reference_control(run_dir: Path, species: List[str],
                                   human_in_panel: bool) -> None:
    present = [p for p in HUMAN_REFERENCE_FILES if (REPO / p).is_file()]
    info: Dict[str, Any] = {
        "enabled": True,
        "source": ("run_panel_and_validated_example_dataset" if human_in_panel
                   else "validated_example_dataset"),
        "human_role": ("analysed_species_plus_reference_control" if human_in_panel
                       else "human_reference_control"),
        "homo_sapiens_in_panel": human_in_panel,
        "species_panel_unchanged": True,
        "human_reference_control_files": present,
        "note": ("homo_sapiens is part of the selected analysed panel and is also linked "
                 "to the curated human FGFR2 IIIb/IIIc reference."
                 if human_in_panel else
                 "Curated human FGFR2 IIIb/IIIc is reused only as a reference/control "
                 "layer and is not added as an analysed species."),
    }
    update_status(run_dir, human_reference=info)
    cfg_path = run_dir / "run_config.json"
    cfg = read_json(cfg_path, {}) or {}
    cfg["human_reference"] = info
    write_json(cfg_path, cfg)


def mark_failed(run_dir: Path, message: str, log_path: Path,
                err_path: Optional[Path] = None,
                failed_step: str = "", failed_reason: str = "") -> None:
    update_status(run_dir, status="failed", current_step="pre_interpro_failed",
                  pre_interpro_status="failed", error=message,
                  failed_step=failed_step, failed_reason=failed_reason,
                  log_file=rel(log_path),
                  err_file=rel(err_path) if err_path else "")


def _tail(text: str, n: int) -> str:
    lines = text.splitlines()
    return "\n".join(lines[-n:])


def dump_failure_context(run_dir: Path, results_dir: Path, env: Dict[str, str],
                         species: List[str], logs_dir: Path,
                         log_path: Path, err_path: Path) -> Dict[str, str]:
    closure_log = results_dir / "13_final_pre_interpro_closure" / "final_pre_interpro_run_log.txt"
    step_status = results_dir / "13_final_pre_interpro_closure" / "final_pre_interpro_step_status.tsv"
    ctx_path = logs_dir / "pre_interpro_failure_context.log"

    closure_text = ""
    if closure_log.exists():
        try:
            closure_text = closure_log.read_text(encoding="utf-8", errors="replace")
        except OSError:
            closure_text = ""

    failed_step, failed_cmd = "", ""
    for line in closure_text.splitlines():
        if ">>> STEP" in line:
            failed_step = line.split(">>> STEP", 1)[1].strip()
        if "CMD:" in line and failed_step:
            failed_cmd = line.split("CMD:", 1)[1].strip()
    # A concise reason line: the first [ERROR]/Traceback marker near the end.
    reason = ""
    for line in reversed(closure_text.splitlines()):
        low = line.lower()
        if "[error]" in low or "traceback" in low or "did not produce" in low:
            reason = line.strip()
            break

    tail100 = _tail(closure_text, 100)
    species_preview = ", ".join(species[:12]) + (" ..." if len(species) > 12 else "")

    parts: List[str] = []
    parts.append("=" * 70)
    parts.append("  PRE-INTERPRO FAILURE CONTEXT")
    parts.append("=" * 70)
    parts.append(f"  failed step   : {failed_step or '<unknown>'}")
    parts.append(f"  substep cmd   : {failed_cmd or '<unknown>'}")
    parts.append(f"  reason        : {reason or '<see substep output below>'}")
    parts.append(f"  BASE          : {env.get('BASE', '')}")
    parts.append(f"  SPECIES_LIST  : {env.get('SPECIES_LIST', '')}")
    parts.append(f"  SKIP_V3       : {env.get('SKIP_V3', '')}   SKIP_MSA: {env.get('SKIP_MSA', '')}   NO_ENSEMBL_REST: {env.get('NO_ENSEMBL_REST', '')}")
    parts.append(f"  species_count : {len(species)}")
    parts.append(f"  species       : {species_preview}")
    parts.append(f"  closure log   : {rel(closure_log)}")
    parts.append(f"  step status   : {rel(step_status)}")
    parts.append("-" * 70)
    parts.append("  last 100 lines of failed substep output (closure log):")
    parts.append("-" * 70)
    parts.append(tail100 if tail100 else "  <closure log empty or missing>")
    parts.append("=" * 70)
    ctx = "\n".join(parts)

    try:
        ctx_path.write_text(ctx + "\n", encoding="utf-8")
    except OSError:
        pass

    # Dedicated per-step log for the failing v3 step (A1) when identifiable.
    if failed_step.startswith("A1") and closure_text:
        try:
            (logs_dir / "step_A1_v3.log").write_text(closure_text, encoding="utf-8")
        except OSError:
            pass

    print(ctx)
    return {"failed_step": failed_step, "failed_reason": reason,
            "context_file": rel(ctx_path)}


def mark_complete(run_dir: Path, run_id: str, primary_n: Optional[int],
                  review_n: Optional[int]) -> None:
    update_status(
        run_dir,
        status="pre_interpro_complete",
        current_step="pre_interpro_complete",
        pre_interpro_status="complete",
        primary_fasta_status="available" if primary_n is not None else "not_available",
        review_fasta_status="available" if review_n is not None else "not_available",
        primary_fasta_count=primary_n if primary_n is not None else 0,
        review_fasta_count=review_n if review_n is not None else 0,
        error="",
        # Clear any stale failure context from a previous failed attempt so the UI does not
        # keep flagging a now-successful run as "failed".
        failed_reason="",
        failed_step="",
        next_action="run_cluster_roundtrip_command",
        next_actions=[
            "Submit the cluster analysis (InterProScan + pyTMHMM) from your local "
            "terminal: "
            f"python scripts/interpro_cluster/submit_cluster_analysis.py --run-id {run_id}",
            "Then check and fetch results: "
            f"python scripts/interpro_cluster/check_cluster_analysis.py --run-id {run_id}; "
            f"python scripts/interpro_cluster/fetch_cluster_analysis.py --run-id {run_id}",
        ],
    )


def update_config_fastas(run_dir: Path, primary: Path, review: Path,
                         primary_n: Optional[int], review_n: Optional[int]) -> None:
    cfg_path = run_dir / "run_config.json"
    cfg = read_json(cfg_path, {}) or {}
    if primary.exists():
        cfg["primary_fasta_path"] = rel(primary)
        cfg["primary_fasta_count"] = primary_n
    if review.exists():
        cfg["review_fasta_path"] = rel(review)
        cfg["review_fasta_count"] = review_n
    cfg["pre_interpro_completed_at"] = now_iso()
    write_json(cfg_path, cfg)


def update_readme(run_dir: Path, run_id: str, primary_n: Optional[int]) -> None:
    readme = run_dir / "00_README_NEXT_STEPS.md"
    if not readme.exists():
        return
    text = readme.read_text(encoding="utf-8")
    banner = (
        f"> **Pre-InterPro complete** ({now_iso()}). "
        f"Primary FASTA: {primary_n if primary_n is not None else '?'} sequences. "
        "Cluster submit/check/fetch commands are ready below.\n"
    )
    if "Pre-InterPro complete" not in text:
        # insert the banner right after the top read-only notice paragraph
        marker = "not touched.\n"
        idx = text.find(marker)
        if idx >= 0:
            cut = idx + len(marker)
            text = text[:cut] + "\n" + banner + text[cut:]
        else:
            text = banner + "\n" + text
    # flip the pre-InterPro step status label to complete
    for planned in ("### 1. Pre-InterPro pipeline — status: planned",
                    "### 1. Pre-InterPro pipeline — status: available"):
        text = text.replace(planned, "### 1. Pre-InterPro pipeline — status: complete")
    readme.write_text(text, encoding="utf-8")


# --------------------------------------------------------------------------- #
# core
# --------------------------------------------------------------------------- #
def build_env(run_id: str, run_dir: Path, results_dir: Path, species_path: Path,
              pre_interpro_dir: Path, env_mode: str, force: bool,
              skip_v3: bool, skip_msa: bool) -> Dict[str, str]:
    env = os.environ.copy()
    env["RUN_ID"] = run_id
    env["RUN_DIR"] = rel(run_dir)
    env["RESULTS_DIR"] = rel(results_dir)
    # BASE keeps every generated result inside the run folder.
    env["BASE"] = rel(results_dir)
    env["SPECIES_LIST"] = rel(species_path)
    env["PRE_INTERPRO_DIR"] = rel(pre_interpro_dir)
    env["PYTHON"] = detect_python()
    # The product name supports offline target/paralog identification.
    config = _run_config(run_dir).get("gene_config")
    if config and (REPO / config).is_file():
        env["GENE_CONFIG"] = config
    # Share a writable Matplotlib cache across figure subprocesses.
    if not env.get("MPLCONFIGDIR"):
        mpl_cache = REPO / ".cache" / "matplotlib"
        try:
            mpl_cache.mkdir(parents=True, exist_ok=True)
            env["MPLCONFIGDIR"] = str(mpl_cache)
        except OSError:
            pass
    # Small custom panels may legitimately lack full-panel paper figures.
    env["PAPER_FIGURES_OPTIONAL"] = "1"
    if force:
        env["FORCE"] = "1"
    # Reuse the early pipeline only when its complete run-local cache exists.
    env["SKIP_V3"] = "1" if skip_v3 else "0"
    if env_mode == "cached":
        # Protein-fetch and MSA reuse depend on separate validated caches.
        env["SKIP_MSA"] = "1" if skip_msa else "0"
        env["NO_ENSEMBL_REST"] = "1" if skip_v3 else "0"
    else:  # live
        env["SKIP_MSA"] = "0"
        env["NO_ENSEMBL_REST"] = "0"
    return env


ENV_KEYS = ["RUN_ID", "RUN_DIR", "RESULTS_DIR", "BASE", "SPECIES_LIST",
            "PRE_INTERPRO_DIR", "PYTHON", "GENE_CONFIG", "FORCE", "SKIP_V3",
            "SKIP_MSA", "NO_ENSEMBL_REST"]


def _run_config(run_dir: Path) -> Dict[str, str]:
    try:
        return json.loads((run_dir / "run_config.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def print_v3_decision(env_mode: str, cache_ok: bool, missing: List[str],
                      effective_skip_v3: bool) -> None:
    if env_mode == "cached" and cache_ok:
        print("cached mode: complete run-local v3 outputs found; SKIP_V3=1")
    elif env_mode == "cached":
        print("cached mode requested, but run-local v3 outputs are incomplete; "
              "running required v3 steps (SKIP_V3=0)")
        for m in missing:
            print(f"    missing: {m}")
    else:
        print(f"live mode: running full v3 (SKIP_V3={'1' if effective_skip_v3 else '0'})")


def write_manifest(setup_dir: Path, run_id: str, species: List[str],
                   species_path: Path, results_dir: Path, command: str,
                   start: str, end: Optional[str], exit_code: Optional[int],
                   primary: Path, review: Path, primary_n: Optional[int],
                   review_n: Optional[int], status: str,
                   extra: Optional[Dict[str, Any]] = None) -> Path:
    setup_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "run_id": run_id,
        "species_count": len(species),
        "species_list_path": rel(species_path),
        "results_dir": rel(results_dir),
        "start_time": start,
        "end_time": end,
        "command": command,
        "exit_code": exit_code,
        "primary_fasta_path": rel(primary) if primary.exists() else None,
        "review_fasta_path": rel(review) if review.exists() else None,
        "primary_fasta_count": primary_n,
        "review_fasta_count": review_n,
        "status": status,
    }
    if extra:
        manifest.update(extra)
    path = setup_dir / "pre_interpro_run_manifest.json"
    write_json(path, manifest)
    return path


def safety_check_results_dir(results_dir: Path) -> None:
    resolved = results_dir.resolve()
    if resolved == FREEZE_DIR or str(resolved).startswith(str(FREEZE_DIR) + os.sep):
        raise SystemExit(
            "ERROR: refusing to run — results dir resolves inside the read-only "
            f"example freeze ({rel(FREEZE_DIR)}). This wrapper only writes into runs/."
        )
    if not str(resolved).startswith(str(RUNS_ROOT.resolve()) + os.sep):
        raise SystemExit(
            f"ERROR: results dir must be inside runs/ (got {resolved}).")


def main(argv: Optional[List[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        description="Run the pre-InterPro FGFR2 pipeline inside a runs/<run_id>/ folder.")
    ap.add_argument("--run-id", required=True, help="run id under runs/")
    ap.add_argument("--force", action="store_true", help="pass FORCE=1 to the pipeline")
    ap.add_argument("--dry-run", action="store_true",
                    help="print command + environment and validate only")
    ap.add_argument("--skip-existing", action="store_true",
                    help="reuse an existing primary FASTA (no rerun)")
    ap.add_argument("--env", choices=["cached", "live"], default="cached",
                    help="cached (default): reuse run outputs, no download; live: full run")
    ap.add_argument("--config")
    ap.add_argument("--local-profile")
    ap.add_argument("--lrz-profile")
    args = ap.parse_args(argv)

    global RUNTIME_CONFIG, REPO, RUNS_ROOT, FREEZE_DIR, CLOSURE_RUNNER
    RUNTIME_CONFIG = load_config(
        config_path=args.config, repository_root=REPO,
        local_profile=args.local_profile, lrz_profile=args.lrz_profile,
    )
    REPO = RUNTIME_CONFIG.repository_root
    RUNS_ROOT = RUNTIME_CONFIG.runs_root
    _FREEZE_DIR = (REPO / "results" / "final_30_until_interpro_prepare").resolve()
    CLOSURE_RUNNER = REPO / "run_fgfr2_pipeline_current_final_pre_interpro.sh"

    # --- Part B: read + validate run config ---
    record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    run_dir = record.path if record else (RUNS_ROOT / args.run_id).resolve()
    if record and record.read_only and not args.dry_run:
        raise SystemExit("ERROR: selected legacy run is registered read-only; copy it before retry/resume.")
    if not run_dir.is_dir():
        raise SystemExit(f"ERROR: run folder not found: {rel(run_dir)}. "
                         "Create it first with scripts/create_new_run.py.")
    adapter = LegacyRunAdapter(run_dir, expected_run_id=args.run_id)
    cfg = adapter.config()
    species = adapter.species()
    if not args.dry_run:
        adapter.materialize_legacy_compatibility()
    species_path = run_dir / "species_list.txt"
    results_dir = run_dir / "results"

    problems: List[str] = []
    if len(species) < 1:
        problems.append("canonical species configuration has no species")
    cs = cfg.get("case_study", "")
    if cs and cs != CASE_STUDY:
        problems.append(f"case_study is '{cs}', expected {CASE_STUDY}")
    if not CLOSURE_RUNNER.exists():
        problems.append(f"pipeline runner not found: {rel(CLOSURE_RUNNER)}")
    if problems:
        raise SystemExit("ERROR: run validation failed:\n  - " + "\n  - ".join(problems))
    if not args.dry_run:
        results_dir.mkdir(parents=True, exist_ok=True)

    # Preflight: refuse invalid species identifiers before Step A1 (guards legacy
    # runs whose species_list.txt still contains 'genus species' with a space).
    preflight_species_list(species)

    # Human reference controls do not add human to the analysed species panel.
    human_in_panel = any((s or "").strip().lower() == "homo_sapiens" for s in species)
    record_human_reference_control(run_dir, species, human_in_panel)

    safety_check_results_dir(results_dir)

    pre_interpro_dir = results_dir / "13_final_pre_interpro_closure"
    primary = run_dir / PRIMARY_FASTA_REL
    review = run_dir / REVIEW_FASTA_REL
    logs_dir = run_dir / "logs"
    setup_dir = results_dir / "00_run_setup"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / "pre_interpro_pipeline.log"
    err_path = logs_dir / "pre_interpro_pipeline.err"

    # --- v3 cache decision: never blindly skip Steps 1-11 ---
    cache_ok, cache_found, cache_missing = has_required_v3_outputs(results_dir)
    effective_skip_v3 = cache_ok if args.env == "cached" else False
    # --- module-12 (MSA/rescue) cache decision: INDEPENDENT of the v3 cache. Only
    #     skip the MSA module when its post-rescue truth table already exists, else
    #     the closure step (A5) is stranded without its single source of truth. ---
    msa_ok, msa_found, msa_missing = has_required_msa_outputs(results_dir)
    effective_skip_msa = msa_ok if args.env == "cached" else False
    v3_info: Dict[str, Any] = {
        "requested_env": args.env,
        "effective_skip_v3": effective_skip_v3,
        "v3_cache_found": cache_found,
        "v3_cache_missing_files": cache_missing,
        "effective_skip_msa": effective_skip_msa,
        "msa_cache_missing_files": msa_missing,
    }

    env = build_env(args.run_id, run_dir, results_dir, species_path,
                    pre_interpro_dir, args.env, args.force,
                    effective_skip_v3, effective_skip_msa)
    command = f"bash {rel(CLOSURE_RUNNER)}"

    # --- Part G: dry-run ---
    if args.dry_run:
        print("=" * 64)
        print("  DRY RUN — pre-InterPro pipeline")
        print("=" * 64)
        print(f"  run_id        : {args.run_id}")
        print(f"  run_dir       : {rel(run_dir)}")
        print(f"  species_count : {len(species)}")
        print(f"  requested mode: {args.env}")
        print(f"  v3 cache      : {'complete' if cache_ok else 'incomplete'} "
              f"({len(cache_found)}/{len(REQUIRED_V3_OUTPUTS)} required outputs present)")
        print(f"  effective SKIP_V3 : {'1' if effective_skip_v3 else '0'}")
        print(f"  module-12 cache   : {'complete' if msa_ok else 'incomplete'} "
              f"({len(msa_found)}/{len(REQUIRED_MSA_OUTPUTS)} required MSA outputs present)")
        print(f"  effective SKIP_MSA: {'1' if effective_skip_msa else '0'}")
        if cache_missing:
            print("  missing v3 outputs:")
            for m in cache_missing:
                print(f"      - {m}")
        print(f"  command       : {command}")
        print("  environment   :")
        for k in ENV_KEYS:
            if k in env:
                print(f"      {k}={env[k]}")
        print(f"  primary FASTA : {rel(primary)} "
              f"({'exists' if primary.exists() else 'not yet'})")
        print(f"  freeze safe   : results dir is inside runs/ and not the example freeze ✓")
        print("  note          : dry run — nothing executed, no files changed except dry-run log")
        print("=" * 64)
        write_json(setup_dir / "pre_interpro_dry_run.json", {
            "run_id": args.run_id, "when": now_iso(), "command": command,
            "env": {k: env[k] for k in ENV_KEYS if k in env},
            "species_count": len(species),
            **v3_info,
        })
        return

    # --- Part H: skip-existing ---
    if args.skip_existing and primary.exists():
        primary_n = count_fasta(primary)
        review_n = count_fasta(review)
        print(f"[skip-existing] primary FASTA already present "
              f"({primary_n} sequences) — not rerunning the pipeline.")
        update_config_fastas(run_dir, primary, review, primary_n, review_n)
        mark_complete(run_dir, args.run_id, primary_n, review_n)
        write_manifest(setup_dir, args.run_id, species, species_path, results_dir,
                       command + "  [skipped: existing FASTA reused]", now_iso(),
                       now_iso(), 0, primary, review, primary_n, review_n,
                       "complete_skipped_existing", extra=v3_info)
        update_readme(run_dir, args.run_id, primary_n)
        _print_summary(args.run_id, run_dir, primary_n, review_n, primary)
        return

    # --- Part D/E: run the pipeline, capturing logs ---
    start = now_iso()
    mark_running(run_dir)
    update_status(run_dir, **v3_info)
    write_manifest(setup_dir, args.run_id, species, species_path, results_dir,
                   command, start, None, None, primary, review, None, None, "running",
                   extra=v3_info)

    print(f"[run] pre-InterPro pipeline for {args.run_id} (env={args.env})")
    print_v3_decision(args.env, cache_ok, cache_missing, effective_skip_v3)
    print(f"[run] BASE={env['BASE']}  logs -> {rel(log_path)}")
    exit_code = 1
    try:
        with log_path.open("w", encoding="utf-8") as out, \
                err_path.open("w", encoding="utf-8") as errf:
            proc = subprocess.run(["bash", str(CLOSURE_RUNNER)], cwd=str(REPO),
                                  env=env, stdout=out, stderr=errf, check=False)
            exit_code = proc.returncode
    except Exception as exc:  # pragma: no cover - unexpected launch failure
        mark_failed(run_dir, f"failed to launch pipeline: {exc}", log_path)
        write_manifest(setup_dir, args.run_id, species, species_path, results_dir,
                       command, start, now_iso(), exit_code, primary, review,
                       None, None, "failed", extra=v3_info)
        raise SystemExit(f"ERROR: could not launch pipeline: {exc}")

    end = now_iso()
    primary_n = count_fasta(primary)
    review_n = count_fasta(review)

    if exit_code != 0:
        ctx = dump_failure_context(run_dir, results_dir, env, species,
                                   logs_dir, log_path, err_path)
        mark_failed(run_dir,
                    f"pipeline exited with code {exit_code}", log_path, err_path,
                    failed_step=ctx.get("failed_step", ""),
                    failed_reason=ctx.get("failed_reason", ""))
        write_manifest(setup_dir, args.run_id, species, species_path, results_dir,
                       command, start, end, exit_code, primary, review,
                       primary_n, review_n, "failed", extra=v3_info)
        raise SystemExit(
            f"ERROR: pre-InterPro pipeline failed (exit {exit_code}).\n"
            f"  failed step   : {ctx.get('failed_step') or '<unknown>'}\n"
            f"  reason        : {ctx.get('failed_reason') or '<see context>'}\n"
            f"  context       : {ctx.get('context_file', '')}\n"
            f"  log: {rel(log_path)}\n  err: {rel(err_path)}")

    # success
    update_config_fastas(run_dir, primary, review, primary_n, review_n)
    mark_complete(run_dir, args.run_id, primary_n, review_n)
    write_manifest(setup_dir, args.run_id, species, species_path, results_dir,
                   command, start, end, exit_code, primary, review,
                   primary_n, review_n, "complete", extra=v3_info)
    update_readme(run_dir, args.run_id, primary_n)
    rebuilt = refresh_derived_layer(run_dir, args.run_id)
    _print_summary(args.run_id, run_dir, primary_n, review_n, primary, rebuilt)


def invalidate_derived_layer(run_dir: Path) -> None:
    update_status(
        run_dir,
        post_interpro_status="stale",
        website_indices_status="stale",
        derived_layer_invalidated_at=now_iso(),
    )


def refresh_derived_layer(run_dir: Path, run_id: str) -> Dict[str, Any]:
    invalidate_derived_layer(run_dir)
    info: Dict[str, Any] = {"cluster_outputs_reused": False, "indices_rebuilt": False}
    py = detect_python()

    if _cluster_outputs_valid(run_dir):
        info["cluster_outputs_reused"] = True
        cmd = [py, str(REPO / "scripts" / "run_post_interpro_for_run.py"),
               "--run-id", run_id]
    else:
        closure = run_dir / "results" / "13_final_pre_interpro_closure"
        cmd = [py, str(REPO / "scripts" / "build_website_indices.py"),
               "--run-dir", rel(closure),
               "--outdir", rel(run_dir / "website_indices")]

    log = run_dir / "logs" / "derived_layer_rebuild.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log.open("w", encoding="utf-8") as handle:
            handle.write(f"[{now_iso()}] $ {' '.join(cmd)}\n")
            handle.flush()
            proc = subprocess.run(cmd, cwd=str(REPO), stdout=handle,
                                  stderr=subprocess.STDOUT, check=False)
        info["exit_code"] = proc.returncode
        if proc.returncode == 0 and not info["cluster_outputs_reused"]:
            adapter_cmd = [
                py, "-m", "exondomaincompare.adapters.fgfr2_core_analysis_adapter",
                "--run-id", run_id,
            ]
            with log.open("a", encoding="utf-8") as handle:
                handle.write(f"[{now_iso()}] $ {' '.join(adapter_cmd)}\n")
                handle.flush()
                adapter_proc = subprocess.run(
                    adapter_cmd, cwd=str(REPO), stdout=handle,
                    stderr=subprocess.STDOUT, check=False)
            info["adapter_exit_code"] = adapter_proc.returncode
            if adapter_proc.returncode != 0:
                info["exit_code"] = adapter_proc.returncode
    except Exception as exc:  # pragma: no cover - launch failure
        info["error"] = f"{type(exc).__name__}: {exc}"
        update_status(run_dir, website_indices_status="failed")
        return info

    if info.get("exit_code") == 0:
        info["indices_rebuilt"] = True
        update_status(run_dir, website_indices_status="complete")
    else:
        # An index rebuild that fails must leave the run marked stale rather than
        # complete, so the readiness gate keeps the run out of "Results ready".
        update_status(run_dir, website_indices_status="failed")
    return info


def _cluster_outputs_valid(run_dir: Path) -> bool:
    ips = run_dir / "results" / "14_interproscan" / "primary" / "output"
    tm = (run_dir / "results" / "15_exon_domain_boundary_post_interpro"
          / "pytmhmm_primary" / "output")
    has_ips = ips.is_dir() and any(p.suffix == ".tsv" and p.stat().st_size > 0
                                   for p in ips.rglob("*") if p.is_file())
    has_tm = tm.is_dir() and any(p.is_file() and p.stat().st_size > 0
                                 for p in tm.rglob("*"))
    return has_ips and has_tm


def _print_summary(run_id: str, run_dir: Path, primary_n: Optional[int],
                   review_n: Optional[int], primary: Path,
                   rebuilt: Optional[Dict[str, Any]] = None) -> None:
    print("\n" + "=" * 64)
    print("  PRE-INTERPRO COMPLETE")
    print("=" * 64)
    print(f"  run_id        : {run_id}")
    print(f"  run_dir       : {rel(run_dir)}")
    print(f"  primary FASTA : {primary_n if primary_n is not None else 'missing'} sequences")
    print(f"  review FASTA  : {review_n if review_n is not None else 'missing'} sequences")
    if primary_n is not None and primary_n != FULL30_EXPECTED_PRIMARY:
        print(f"  note          : full30 reference freeze is {FULL30_EXPECTED_PRIMARY}/"
              f"{FULL30_EXPECTED_REVIEW}; this run differs (informational only).")
    if rebuilt is not None:
        how = ("post-cluster analysis re-run, returned cluster output reused"
               if rebuilt.get("cluster_outputs_reused") else "website indices rebuilt")
        state = "ok" if rebuilt.get("indices_rebuilt") else "FAILED — see log"
        print(f"  derived layer : {how} ({state})")
    if not (rebuilt or {}).get("cluster_outputs_reused"):
        print(f"  next command  : python scripts/interpro_cluster/submit_cluster_analysis.py --run-id {run_id}")
    print("=" * 64 + "\n")


if __name__ == "__main__":
    main()
