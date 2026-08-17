#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.contracts import portable_runtime_record, stamp_payload  # noqa: E402
from exondomaincompare.runs.legacy import LegacyRunAdapter  # noqa: E402
from exondomaincompare.runs.registry import resolve_run_record  # noqa: E402
from exondomaincompare.config import load_config  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=Path(__file__).resolve().parent.parent)
REPO = RUNTIME_CONFIG.repository_root
SCRIPTS = REPO / "scripts"
RUNS_ROOT = RUNTIME_CONFIG.runs_root
FREEZE_DIR = (REPO / "results" / "final_30_until_interpro_prepare").resolve()

CASE_STUDY = "FGFR2_IIIb_IIIc"
PRIMARY_FASTA_REL = "results/13_final_pre_interpro_closure/freeze/final_pre_interpro_proteins_primary.faa"

ALL_STEPS = ["architecture", "boundary", "audit", "indices"]


# --------------------------------------------------------------------------- #
# utilities
# --------------------------------------------------------------------------- #
def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def rel(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO))
    except ValueError:
        return str(Path(path).resolve())


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
            generator="scripts/run_post_interpro_for_run.py",
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def detect_python() -> str:
    return RUNTIME_CONFIG.local_python().selected


def first_match(folder: Path, suffixes: Tuple[str, ...]) -> Optional[Path]:
    if not folder.is_dir():
        return None
    for p in sorted(folder.rglob("*")):
        if p.is_file() and p.name.lower().endswith(suffixes):
            return p
    return None


# --------------------------------------------------------------------------- #
# status updates
# --------------------------------------------------------------------------- #
def _finalize_run_status(run_dir: Path) -> None:
    try:
        scripts_dir = str(Path(__file__).resolve().parent)
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        from exondomaincompare.shared_gene_analysis.finalize_run_status import finalize
        report = finalize(run_dir)
        print(f"  run status       : {report.get('status')} ({report.get('reason')})")
    except Exception as err:  # pragma: no cover - status must never break the run
        print(f"  WARN: canonical status finalisation skipped: {err}")


def update_status(run_dir: Path, **fields: Any) -> None:
    status_path = run_dir / "status.json"
    status = read_json(status_path, {}) or {}
    status.update(fields)
    status["last_updated"] = now_iso()
    write_json(status_path, status)


# --------------------------------------------------------------------------- #
# run paths + environment
# --------------------------------------------------------------------------- #
class RunPaths:
    def __init__(self, run_id: str):
        self.run_id = run_id
        record = resolve_run_record(RUNTIME_CONFIG, run_id)
        self.record = record
        self.run_dir = (
            record.path if record else (RUNS_ROOT / run_id).resolve())
        self.results = self.run_dir / "results"
        self.closure = self.results / "13_final_pre_interpro_closure"
        self.interpro_primary = self.results / "14_interproscan" / "primary"
        self.post = self.results / "15_exon_domain_boundary_post_interpro"
        self.pytmhmm_primary = self.post / "pytmhmm_primary"
        self.final_analysis = self.results / "16_final_thesis_analyses"
        self.website_indices = self.run_dir / "website_indices"
        self.logs = self.run_dir / "logs"
        self.setup = self.results / "00_run_setup"
        self.primary_fasta = self.run_dir / PRIMARY_FASTA_REL
        self.truth_table = self.closure / "final_pre_interpro_truth_table.tsv"
        self.interpro_out = self.interpro_primary / "output"
        self.pytmhmm_out = self.pytmhmm_primary / "output"


def build_env(rp: RunPaths) -> Dict[str, str]:
    env = os.environ.copy()
    env["RUN_ID"] = rp.run_id
    env["RUN_DIR"] = rel(rp.run_dir)
    env["BASE"] = rel(rp.results)
    env["RESULTS_DIR"] = rel(rp.results)
    env["FGFR2_RESULTS_DIR"] = rel(rp.results)
    env["PRE_INTERPRO_DIR"] = rel(rp.closure)
    env["INTERPROSCAN_PRIMARY_DIR"] = rel(rp.interpro_primary)
    env["PYTMHMM_PRIMARY_DIR"] = rel(rp.pytmhmm_primary)
    env["POST_INTERPRO_DIR"] = rel(rp.post)
    env["FINAL_ANALYSIS_DIR"] = rel(rp.final_analysis)
    env["WEBSITE_INDICES_DIR"] = rel(rp.website_indices)
    env["PYTHON"] = detect_python()
    return env


ENV_KEYS = ["RUN_ID", "RUN_DIR", "BASE", "RESULTS_DIR", "FGFR2_RESULTS_DIR",
            "PRE_INTERPRO_DIR", "INTERPROSCAN_PRIMARY_DIR", "PYTMHMM_PRIMARY_DIR",
            "POST_INTERPRO_DIR", "FINAL_ANALYSIS_DIR", "WEBSITE_INDICES_DIR", "PYTHON"]


# --------------------------------------------------------------------------- #
# step -> command list
# --------------------------------------------------------------------------- #
def step_commands(rp: RunPaths, py: str) -> Dict[str, List[List[str]]]:
    return {
        # exon-block reconstruction -> sanitation (merges) -> architecture figures.
        # Figures consume the merged override file, so a single figures pass AFTER
        # sanitation is correct (documented in the manifest).
        "architecture": [
            [py, str(SCRIPTS / "reconstruct_exon_blocks_post_interpro.py")],
            [py, str(SCRIPTS / "sanitize_exon_block_coordinates.py")],
            [py, str(SCRIPTS / "make_fgfr2_post_interpro_exon_domain_figures.py")],
        ],
        "boundary": [
            [py, str(SCRIPTS / "analyze_exon_domain_boundary_consistency.py")],
            [py, str(SCRIPTS / "build_human_reference_sanity_check.py")],
        ],
        "audit": [
            [py, str(SCRIPTS / "audit_review_and_minor_flag_cases.py")],
        ],
        "indices": [
            [py, str(SCRIPTS / "build_website_indices.py"),
             "--run-dir", rel(rp.closure), "--outdir", rel(rp.website_indices)],
            [py, "-m", "exondomaincompare.adapters.fgfr2_core_analysis_adapter",
             "--run-id", rp.run_id],
        ],
    }


# best-effort sub-commands that should not fail the whole run if they error
SOFT_COMMANDS = {"build_human_reference_sanity_check.py"}


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #
def validate_inputs(rp: RunPaths, cfg: Dict[str, Any]) -> Tuple[List[str], Dict[str, Any]]:
    problems: List[str] = []
    found: Dict[str, Any] = {}

    if not rp.run_dir.is_dir():
        return ([f"run folder not found: {rel(rp.run_dir)}"], found)

    cs = cfg.get("case_study", "")
    if cs and cs != CASE_STUDY:
        problems.append(f"case_study is '{cs}', expected {CASE_STUDY}")

    if rp.primary_fasta.exists():
        found["primary_fasta"] = rel(rp.primary_fasta)
    else:
        problems.append(f"primary FASTA missing: {rel(rp.primary_fasta)}")

    if rp.truth_table.exists():
        found["truth_table"] = rel(rp.truth_table)
    else:
        problems.append(f"truth table missing: {rel(rp.truth_table)}")

    # InterProScan TSV (required); GFF3 / JSON optional
    ips_tsv = rp.interpro_out / "input.fasta.tsv"
    if not ips_tsv.exists():
        ips_tsv = first_match(rp.interpro_out, (".tsv",)) or ips_tsv
    if ips_tsv.exists():
        found["interproscan_tsv"] = rel(ips_tsv)
    else:
        problems.append(
            f"InterProScan TSV missing: {rel(rp.interpro_out / 'input.fasta.tsv')} "
            "(fetch cluster results first)")
    for opt, suf in (("interproscan_gff3", ".gff3"), ("interproscan_json", ".json")):
        m = (rp.interpro_out / f"input.fasta{suf}")
        if not m.exists():
            m = first_match(rp.interpro_out, (suf,))
        if m and m.exists():
            found[opt] = rel(m)

    # pyTMHMM output (required): any summary/tsv/hits file under the output dir
    tm = first_match(rp.pytmhmm_out, (".summary", ".tsv", ".txt", ".json", ".gff3"))
    if tm:
        found["pytmhmm_output"] = rel(tm)
    else:
        problems.append(
            f"pyTMHMM output missing under {rel(rp.pytmhmm_out)} (fetch cluster results first)")

    return (problems, found)


def resolve_steps(args: argparse.Namespace) -> Tuple[List[str], List[str]]:
    if args.steps:
        requested = [s.strip() for s in args.steps.split(",") if s.strip()]
        unknown = [s for s in requested if s not in ALL_STEPS]
        if unknown:
            raise SystemExit(f"ERROR: unknown --steps values: {', '.join(unknown)}. "
                             f"Valid: {', '.join(ALL_STEPS)}")
        steps = [s for s in ALL_STEPS if s in requested]
    else:
        steps = list(ALL_STEPS)
    skipped: List[str] = []
    if args.skip_boundary and "boundary" in steps:
        steps.remove("boundary"); skipped.append("boundary")
    if args.skip_audit and "audit" in steps:
        steps.remove("audit"); skipped.append("audit")
    if args.skip_website_indices and "indices" in steps:
        steps.remove("indices"); skipped.append("indices")
    return steps, skipped


# --------------------------------------------------------------------------- #
# core
# --------------------------------------------------------------------------- #
def outputs_exist(rp: RunPaths) -> bool:
    idx = rp.website_indices / "boundary_consistency_index.json"
    bc = rp.final_analysis / "exon_domain_boundary_consistency"
    return idx.exists() and bc.is_dir()


def write_manifest(rp: RunPaths, start: str, end: Optional[str], status: str,
                   executed: List[Dict[str, Any]], planned_steps: List[str],
                   skipped: List[str], found: Dict[str, Any],
                   error: str = "") -> Path:
    rp.setup.mkdir(parents=True, exist_ok=True)
    manifest = {
        "run_id": rp.run_id,
        "start_time": start,
        "end_time": end,
        "status": status,
        "planned_steps": planned_steps,
        "executed_steps": [e["step"] for e in executed],
        "skipped_steps": skipped,
        "commands": executed,
        "input_files": found,
        "output_folders": {
            "post_interpro": rel(rp.post),
            "final_analyses": rel(rp.final_analysis),
            "website_indices": rel(rp.website_indices),
        },
        "notes": (
            "Order: exon-block reconstruction -> sanitation (merges overrides) -> "
            "domain-architecture figures (consume merged overrides) -> boundary "
            "consistency -> final audit -> run-local website indices."
        ),
        "error": error,
    }
    path = rp.setup / "post_interpro_run_manifest.json"
    write_json(path, manifest)
    return path


def run_command(cmd: List[str], env: Dict[str, str], log_path: Path,
                err_path: Path) -> int:
    header = f"\n{'=' * 70}\n[{now_iso()}] $ {' '.join(cmd)}\n{'=' * 70}\n"
    with log_path.open("a", encoding="utf-8") as out, \
            err_path.open("a", encoding="utf-8") as errf:
        out.write(header)
        out.flush()
        proc = subprocess.run(cmd, cwd=str(REPO), env=env,
                              stdout=out, stderr=errf, check=False)
        return proc.returncode


def main(argv: Optional[List[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        description="Run local post-InterPro FGFR2 analysis inside a runs/<run_id>/ folder.")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="run even if outputs already exist")
    ap.add_argument("--skip-existing", action="store_true",
                    help="reuse existing final analysis + indices (no rerun)")
    ap.add_argument("--steps", help="comma list: " + ",".join(ALL_STEPS))
    ap.add_argument("--skip-boundary", action="store_true")
    ap.add_argument("--skip-audit", action="store_true")
    ap.add_argument("--skip-website-indices", action="store_true")
    ap.add_argument("--config")
    ap.add_argument("--local-profile")
    ap.add_argument("--lrz-profile")
    args = ap.parse_args(argv)

    global RUNTIME_CONFIG, REPO, SCRIPTS, RUNS_ROOT, FREEZE_DIR
    RUNTIME_CONFIG = load_config(
        config_path=args.config, repository_root=REPO,
        local_profile=args.local_profile, lrz_profile=args.lrz_profile,
    )
    REPO = RUNTIME_CONFIG.repository_root
    _SCRIPTS = REPO / "scripts"
    _RUNS_ROOT = RUNTIME_CONFIG.runs_root
    FREEZE_DIR = (REPO / "results" / "final_30_until_interpro_prepare").resolve()

    rp = RunPaths(args.run_id)
    if not rp.run_dir.is_dir():
        raise SystemExit(f"ERROR: run folder not found: {rel(rp.run_dir)}. "
                         "Create it with scripts/create_new_run.py first.")
    if rp.record and rp.record.read_only and not args.dry_run:
        raise SystemExit(
            "ERROR: selected legacy run is registered read-only; copy it before retry/resume.")
    adapter = LegacyRunAdapter(rp.run_dir, expected_run_id=args.run_id)
    if not args.dry_run:
        adapter.materialize_legacy_compatibility()

    # safety: never operate on the example freeze
    if rp.results.resolve() == FREEZE_DIR or \
            str(rp.results.resolve()).startswith(str(FREEZE_DIR) + os.sep):
        raise SystemExit("ERROR: results dir resolves inside the read-only example freeze.")

    cfg = adapter.config()
    steps, skipped = resolve_steps(args)
    py = detect_python()
    env = build_env(rp)
    cmds = step_commands(rp, py)

    rp.logs.mkdir(parents=True, exist_ok=True)
    rp.setup.mkdir(parents=True, exist_ok=True)
    log_path = rp.logs / "post_interpro_pipeline.log"
    err_path = rp.logs / "post_interpro_pipeline.err"

    problems, found = validate_inputs(rp, cfg)

    # --- Part G: dry run ---
    if args.dry_run:
        print("=" * 66)
        print("  DRY RUN — post-InterPro analysis")
        print("=" * 66)
        print(f"  run_id        : {rp.run_id}")
        print(f"  run_dir       : {rel(rp.run_dir)}")
        print(f"  steps         : {', '.join(steps) or '(none)'}")
        if skipped:
            print(f"  skipped       : {', '.join(skipped)}")
        print("  environment   :")
        for k in ENV_KEYS:
            print(f"      {k}={env.get(k)}")
        print("  planned commands:")
        for s in steps:
            for c in cmds[s]:
                print(f"      [{s}] {' '.join(c)}")
        print("  input validation:")
        for k, v in found.items():
            print(f"      found   {k}: {v}")
        for p in problems:
            print(f"      MISSING {p}")
        print(f"  freeze safe   : results dir is inside runs/ and not the example freeze")
        print("  note          : dry run — nothing executed")
        print("=" * 66)
        write_json(rp.setup / "post_interpro_dry_run.json", {
            "run_id": rp.run_id, "when": now_iso(), "steps": steps,
            "skipped": skipped, "env": {k: env.get(k) for k in ENV_KEYS},
            "commands": {s: [" ".join(c) for c in cmds[s]] for s in steps},
            "found": found, "problems": problems,
        })
        return

    # --- Part H: skip existing ---
    if args.skip_existing and not args.force and outputs_exist(rp):
        print("[skip-existing] final analysis + website indices already present — not rerunning.")
        update_status(rp.run_dir, post_interpro_status="complete",
                      website_indices_status="complete",
                      current_step="analysis_complete", status="complete")
        write_manifest(rp, now_iso(), now_iso(), "complete_skipped_existing",
                       [], steps, skipped, found)
        return

    # --- Part B: hard input gate ---
    if problems:
        report = rp.setup / "post_interpro_missing_inputs.json"
        write_json(report, {"run_id": rp.run_id, "when": now_iso(),
                            "missing": problems, "found": found})
        update_status(rp.run_dir, post_interpro_status="failed",
                      current_step="post_interpro_failed",
                      error="missing required inputs: " + "; ".join(problems),
                      log_file=rel(report))
        raise SystemExit("ERROR: missing required post-InterPro inputs:\n  - "
                         + "\n  - ".join(problems)
                         + f"\nReport: {rel(report)}")

    # fresh logs for a real run
    log_path.write_text("", encoding="utf-8")
    err_path.write_text("", encoding="utf-8")

    start = now_iso()
    update_status(rp.run_dir, post_interpro_status="running",
                  current_step="post_interpro_analysis", error="")
    write_manifest(rp, start, None, "running", [], steps, skipped, found)

    executed: List[Dict[str, Any]] = []
    print(f"[run] post-InterPro analysis for {rp.run_id}  logs -> {rel(log_path)}")
    for step in steps:
        for cmd in cmds[step]:
            name = Path(cmd[1]).name if len(cmd) > 1 else cmd[0]
            print(f"  [{step}] {name} …")
            code = run_command(cmd, env, log_path, err_path)
            executed.append({"step": step, "command": " ".join(cmd), "exit_code": code})
            soft = name in SOFT_COMMANDS
            if code != 0 and not soft:
                update_status(rp.run_dir, post_interpro_status="failed",
                              current_step="post_interpro_failed",
                              error=f"step '{step}' command failed: {name} (exit {code})",
                              log_file=rel(log_path), err_file=rel(err_path))
                write_manifest(rp, start, now_iso(), "failed", executed, steps,
                               skipped, found,
                               error=f"{name} exit {code}")
                raise SystemExit(
                    f"ERROR: post-InterPro step '{step}' failed at {name} (exit {code}).\n"
                    f"  log: {rel(log_path)}\n  err: {rel(err_path)}")
            if code != 0 and soft:
                print(f"    (note) optional step {name} exited {code}; continuing.")

    end = now_iso()
    indices_done = "indices" in steps
    website_status = "complete" if indices_done else \
        (read_json(rp.run_dir / "status.json", {}) or {}).get("website_indices_status", "not_started")

    update_status(rp.run_dir, post_interpro_status="complete",
                  website_indices_status=website_status,
                  current_step="analysis_complete", error="")
    # Current run artefacts determine readiness after repaired earlier failures.
    _finalize_run_status(rp.run_dir)
    write_manifest(rp, start, end, "complete", executed, steps, skipped, found)

    print("\n" + "=" * 66)
    print("  POST-INTERPRO COMPLETE")
    print("=" * 66)
    print(f"  run_id           : {rp.run_id}")
    print(f"  executed steps   : {', '.join(dict.fromkeys(e['step'] for e in executed))}")
    print(f"  post-InterPro out: {rel(rp.post)}")
    print(f"  final analyses   : {rel(rp.final_analysis)}")
    print(f"  website indices  : {rel(rp.website_indices)}")
    print("=" * 66 + "\n")


if __name__ == "__main__":
    main()
