#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.cluster.ssh_common import configure  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from exondomaincompare.contracts import stamp_payload  # noqa: E402
from exondomaincompare.runs.registry import resolve_run_record  # noqa: E402
from exondomaincompare.config import (  # noqa: E402
    CONTROL_PATH_ENV,
    ConfigurationError,
    RuntimeConfig,
    load_config,
)

RUNTIME_CONFIG = load_config(repository_root=Path(__file__).resolve().parents[2])
REPO = RUNTIME_CONFIG.repository_root
CLUSTER_DIR = REPO / "scripts" / "interpro_cluster"
SUBMIT = CLUSTER_DIR / "submit_cluster_analysis.py"
CHECK = CLUSTER_DIR / "check_cluster_analysis.py"
FETCH = CLUSTER_DIR / "fetch_cluster_analysis.py"
POST = REPO / "scripts" / "run_post_interpro_for_run.py"
CORE_POST_MODULE = "exondomaincompare.framework.run_core_gene_analysis"


def _require_cluster_profile(config: RuntimeConfig) -> None:
    try:
        config.require_cluster()
    except ConfigurationError as exc:
        raise SystemExit(
            "Cluster profile is not ready: "
            f"{exc}\nConfigure it once with '.venv/bin/edc cluster configure', then run "
            "'.venv/bin/edc cluster doctor --redact-paths'."
        ) from None


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(path: Path, fallback: Any = None) -> Any:
    if not path.exists():
        return fallback
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return fallback


def write_json(path: Path, data: Any) -> None:
    if isinstance(data, dict):
        run_id = str(data.get("run_id") or path.parent.name)
        data = stamp_payload(
            data, payload_type="status" if path.name == "status.json" else path.stem,
            run_id=run_id, dataset_id=run_id,
            profile=RUNTIME_CONFIG.public_identity(),
            generator="scripts/interpro_cluster/run_cluster_roundtrip.py",
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


class Roundtrip:
    def __init__(self, run_id: str, config: RuntimeConfig | None = None):
        self.run_id = run_id
        config = config or RUNTIME_CONFIG
        self.config = config
        record = resolve_run_record(config, run_id)
        self.record = record
        self.run_dir = record.path if record else config.runs_root / run_id
        self.status_path = self.run_dir / "status.json"
        self.log_path = self.run_dir / "logs" / "cluster_roundtrip.log"
        self.py = config.local_python().selected
        self.control_path: Optional[Path] = None  # SSH ControlMaster socket, if enabled

    def cm_log(self, msg: str) -> None:
        p = self.run_dir / "logs" / "ssh_controlmaster.log"
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("a", encoding="utf-8") as fh:
            fh.write(f"[{now_iso()}] {msg}\n")

    def _control_socket_path(self) -> Path:
        h = hashlib.sha1(self.run_id.encode("utf-8")).hexdigest()[:8]
        return Path("/tmp") / f"fgfr2_{h}_cm.sock"

    def open_ssh_master(self) -> None:
        sock = self._control_socket_path()
        opts = ["-o", "ControlMaster=auto", "-o", f"ControlPath={sock}",
                "-o", "ControlPersist=30m", "-o", "ConnectTimeout=30"]

        if sock.exists():
            subprocess.run([self.config.executable_token("ssh"), "-O", "exit",
                            "-o", f"ControlPath={sock}", self.config.ssh_target],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
            try:
                sock.unlink()
            except OSError:
                pass

        self.cm_log(f"control_path={sock} (len={len(str(sock))})")
        if len(str(sock)) > 100:
            reason = f"control socket path too long ({len(str(sock))} chars): {sock}"
            self.cm_log(reason)
            self.log(f"SSH multiplexing could not be established: {reason}. "
                     "Multiple LRZ password/MFA prompts may be required.")
            return

        master_cmd = [self.config.executable_token("ssh"),
                      *self.config.lrz.get("ssh_options", []),
                      *opts, self.config.ssh_target, "true"]
        self.cm_log(f"opening master: {' '.join(master_cmd)}")
        self.log("Opening a shared SSH connection (aim: one login / MFA for the whole round-trip)…")
        try:
            proc = subprocess.run(master_cmd, cwd=str(REPO))
            open_rc = proc.returncode
        except Exception as exc:
            self.cm_log(f"master open raised: {exc!r}")
            self.log(f"SSH multiplexing could not be established: {exc}. "
                     "Multiple LRZ password/MFA prompts may be required.")
            return

        check = subprocess.run(
            [self.config.executable_token("ssh"), "-O", "check",
             "-o", f"ControlPath={sock}", self.config.ssh_target],
            capture_output=True, text=True, check=False,
        )
        check_out = (check.stderr or check.stdout or "").strip()
        self.cm_log(f"master open exit={open_rc}; 'ssh -O check' exit={check.returncode}; "
                    f"check_output={check_out!r}")

        if check.returncode == 0 and sock.exists():
            self.control_path = sock
            os.environ[CONTROL_PATH_ENV] = str(sock)
            self.log("Shared SSH connection established; submit/check/fetch reuse it "
                     "(no extra MFA prompts).")
            return

        reason = check_out or f"master connection did not persist (open exit {open_rc})"
        self.cm_log(f"falling back to plain ssh/scp; reason: {reason}")
        self.log(f"SSH multiplexing could not be established: {reason}. "
                 "Multiple LRZ password/MFA prompts may be required.")

    def close_ssh_master(self) -> None:
        if not self.control_path:
            return
        subprocess.run(
            [self.config.executable_token("ssh"), "-O", "exit",
             "-o", f"ControlPath={self.control_path}", self.config.ssh_target],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
        )
        self.cm_log("closed shared SSH master connection.")
        os.environ.pop(CONTROL_PATH_ENV, None)
        self.control_path = None

    def log(self, msg: str) -> None:
        line = f"[{now_iso()}] {msg}"
        print(line, flush=True)
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        with self.log_path.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")

    def set_phase(self, phase: str, **extra: Any) -> None:
        st = read_json(self.status_path, {}) or {}
        rt = st.get("cluster_roundtrip", {}) or {}
        rt.update({"phase": phase, "updated_at": now_iso(), **extra})
        st["cluster_roundtrip"] = rt
        st["cluster_profile"] = self.config.lrz_profile_name
        st["cluster_continuation_command"] = self.config.command([
            ".venv/bin/edc", "cluster", "roundtrip",
            "--run-id", self.run_id,
            "--local-profile", self.config.local_profile_name,
            "--lrz-profile", self.config.lrz_profile_name,
        ])
        st["cluster_continuation_compatibility_command"] = self.config.command([
            "python", "scripts/interpro_cluster/run_cluster_roundtrip.py",
            "--run-id", self.run_id,
            "--local-profile", self.config.local_profile_name,
            "--lrz-profile", self.config.lrz_profile_name,
        ])
        stage_map = {
            "submitting": "cluster_submitted", "submitted": "cluster_submitted",
            "polling": "cluster_submitted", "cluster_completed": "cluster_complete",
            "fetching": "cluster_complete", "fetched": "cluster_complete",
            "post_interpro": "post_domain_processing",
        }
        if phase in stage_map:
            st["status"] = stage_map[phase]
            st["current_step"] = stage_map[phase]
        elif phase.endswith("_failed") or phase in ("cluster_error", "poll_timeout"):
            st["status"] = "failed"
            st["current_step"] = phase
            st["failed_reason"] = extra.get("reason") or f"Cluster roundtrip phase failed: {phase}"
        st["last_updated"] = now_iso()
        write_json(self.status_path, st)

    def run(self, script: Path, *args: str) -> int:
        cmd = [
            self.py, str(script), "--run-id", self.run_id,
            "--local-profile", self.config.local_profile_name,
            "--lrz-profile", self.config.lrz_profile_name,
            *args,
        ]
        self.log(f"$ {' '.join(cmd)}")
        proc = subprocess.run(cmd, cwd=str(REPO))
        self.log(f"  -> exit {proc.returncode}")
        return proc.returncode

    def run_module(self, module: str, *args: str) -> int:
        cmd = [
            self.py, "-m", module, "--run-id", self.run_id,
            "--local-profile", self.config.local_profile_name,
            "--lrz-profile", self.config.lrz_profile_name,
            *args,
        ]
        self.log(f"$ {' '.join(cmd)}")
        proc = subprocess.run(cmd, cwd=str(REPO))
        self.log(f"  -> exit {proc.returncode}")
        return proc.returncode

    def cluster_status(self) -> str:
        return str((read_json(self.status_path, {}) or {}).get("cluster_analysis_status", ""))

    def has_job_ids(self) -> bool:
        jobs = (read_json(self.status_path, {}) or {}).get("cluster_jobs", {}) or {}
        return bool(jobs.get("interproscan_job_id") or jobs.get("pytmhmm_job_id"))

    def submit(self) -> None:
        if self.has_job_ids():
            self.log("Cluster job IDs already present in status.json — skipping submit.")
            return
        self.set_phase("submitting")
        rc = self.run(SUBMIT)
        if rc != 0:
            self.set_phase("submit_failed", exit_code=rc)
            raise SystemExit(f"submit failed (exit {rc}).")
        self.set_phase("submitted")

    def poll(self, poll_minutes: float, max_hours: float) -> None:
        interval = max(poll_minutes, 0.5) * 60.0
        deadline = time.time() + max_hours * 3600.0
        self.set_phase("polling", poll_minutes=poll_minutes, max_hours=max_hours)
        first = True
        while True:
            self.run(CHECK)
            state = self.cluster_status()
            self.log(f"cluster_analysis_status = {state or 'unknown'}")
            if state == "completed":
                self.set_phase("cluster_completed")
                self.log("Both cluster jobs completed.")
                return
            if state == "error":
                self.set_phase("cluster_error")
                raise SystemExit("A cluster job failed. Inspect the cluster logs / fetch and check.")
            if time.time() >= deadline:
                self.set_phase("poll_timeout")
                raise SystemExit(f"Timed out after {max_hours} h waiting for cluster completion.")
            if first:
                self.log(f"Waiting for cluster jobs; polling every {max(poll_minutes,0.5)} min "
                         f"(deadline {max_hours} h).")
                first = False
            time.sleep(interval)

    def fetch(self) -> None:
        self.set_phase("fetching")
        rc = self.run(FETCH)
        state = self.cluster_status()
        if rc != 0 or state not in ("fetched_complete", "fetched_incomplete"):
            self.set_phase("fetch_failed", exit_code=rc, cluster_analysis_status=state)
            raise SystemExit(f"fetch failed (exit {rc}, status {state or 'unknown'}).")
        if state == "fetched_incomplete":
            self.log("WARNING: fetch reported incomplete outputs; post-InterPro may fail.")
        self.set_phase("fetched", cluster_analysis_status=state)

    def _is_core_only(self) -> bool:
        rc = read_json(self.run_dir / "run_config.json", {}) or {}
        if str(rc.get("run_mode", "")).lower() == "core_only_pilot":
            return True
        return rc.get("has_event") is False

    def post_interpro(self) -> None:
        core_only = self._is_core_only()
        label = "core post-InterPro" if core_only else "post-InterPro"
        if not core_only and not POST.exists():
            self.log(f"{POST.name} not found — skipping {label}.")
            return
        self.set_phase("post_interpro")
        rc = (self.run_module(CORE_POST_MODULE, "--post")
              if core_only else self.run(POST))
        if rc != 0:
            self.set_phase("post_interpro_failed", exit_code=rc)
            raise SystemExit(f"{label} failed (exit {rc}).")
        self.finalize()

    def finalize(self) -> None:
        st = read_json(self.status_path, {}) or {}
        scientific_status, reason = self._end_state(st)

        if scientific_status is None:
            self.set_phase("post_interpro_failed", reason=reason)
            raise SystemExit(reason)

        st.update({
            "cluster_analysis_status": "complete",
            "cluster_fetch_status": "complete",
            "interproscan_status": "complete",
            "pytmhmm_status": "complete",
            "post_interpro_status": "complete",
            "status": scientific_status,
            "next_action": "open_results",
            "next_actions": [],
            "explorable": True,
            "current_step": scientific_status,
            "readiness_reason": reason,
            "last_updated": now_iso(),
        })
        for stale in ("failed_reason", "failed_step", "failed_species", "error", "detail"):
            st.pop(stale, None)
        write_json(self.status_path, st)
        self.set_phase("complete")
        self.log(f"status.json finalized: {scientific_status} · explorable · "
                 f"next_action=open_results ({reason}).")

    def _end_state(self, st: Dict[str, Any]) -> Tuple[Optional[str], str]:
        if st.get("post_interpro_status") == "failed" or st.get("error"):
            return None, (st.get("failed_reason") or st.get("error")
                          or "Post-cluster processing failed.")

        shared = self._shared_verdict()
        if shared is not None:
            return shared

        verdict = self._readiness(st)
        if verdict is not None:
            ready, reason = verdict
            if ready:
                return "results_ready", reason
            return None, reason

        core = self._core_milestone_state()
        if core is not None:
            return core

        legacy = st.get("status")
        if legacy in ("results_ready", "results_partial", "post_cluster_partial",
                      "complete", "post_interpro_complete"):
            return ("results_ready" if legacy in ("complete", "post_interpro_complete")
                    else legacy), f"post-cluster status was {legacy}"
        return None, ("Post-cluster processing did not produce a valid scientific "
                      "end state.")

    def _shared_verdict(self) -> Optional[Tuple[Optional[str], str]]:
        try:
            sys.path.insert(0, str(REPO / "scripts"))
            from exondomaincompare.shared_gene_analysis.finalize_run_status import evaluate_run
            report = evaluate_run(self.run_dir)
        except Exception:
            return None
        if not report.get("decided"):
            return None
        if report["status"] == "results_ready":
            return "results_ready", report["reason"]
        return None, report["reason"]

    def _core_milestone_state(self) -> Optional[Tuple[Optional[str], str]]:
        try:
            sys.path.insert(0, str(REPO / "scripts"))
            from exondomaincompare.framework.core_run_milestones import evaluate_core_run, is_core_only_run
        except Exception:
            return None
        try:
            if not is_core_only_run(self.run_dir):
                return None
            report = evaluate_core_run(self.run_dir)
        except Exception:
            return None
        inferred = report.get("inferred_status")
        if inferred == "results_ready":
            return "results_ready", "every core milestone is complete"
        return None, (report.get("detail")
                      or f"the core run is {str(inferred or 'incomplete').replace('_', ' ')}")

    def _readiness(self, st: Dict[str, Any]) -> Optional[Tuple[bool, str]]:
        try:
            sys.path.insert(0, str(REPO / "scripts"))
            from exondomaincompare.shared_gene_analysis.run_availability import models_run, readiness
        except Exception:
            return None
        if not models_run(self.run_dir):
            return None
        try:
            rc = read_json(self.run_dir / "run_config.json", {}) or {}
            has_event = (rc.get("gene_symbol") or "FGFR2").upper() == "FGFR2"
            verdict = readiness(self.run_dir,
                                n_species=int(st.get("species_count") or 0),
                                has_event=has_event)
            return verdict.ready, verdict.reason
        except Exception:
            return None


def main(argv: Optional[List[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        description="Local one-command cluster round-trip (submit → poll → fetch → post-InterPro).")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--poll-interval-minutes", "--poll-minutes", dest="poll_minutes",
                    type=float, default=2.0,
                    help="polling interval in minutes (default 2, floor 0.5)")
    ap.add_argument("--timeout-hours", "--max-hours", dest="max_hours",
                    type=float, default=12.0,
                    help="overall deadline in hours (default 12)")
    ap.add_argument("--no-post-interpro", action="store_true",
                    help="stop after fetch (skip local post-InterPro analysis)")
    ap.add_argument("--no-ssh-master", action="store_true",
                    help="disable SSH ControlMaster multiplexing (use plain ssh/scp)")
    ap.add_argument("--submit-only", action="store_true")
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--fetch-only", action="store_true")
    ap.add_argument("--finalize-only", action="store_true",
                    help="re-derive the end state from the run's artefacts and record "
                         "it, without contacting the cluster. Use after the local "
                         "analysis was rebuilt outside a full round-trip.")
    ap.add_argument("--config")
    ap.add_argument("--local-profile")
    ap.add_argument("--lrz-profile")
    ap.add_argument("--dry-run", action="store_true",
                    help="render the selected profile plan without network or writes")
    args = ap.parse_args(argv)

    global RUNTIME_CONFIG, REPO, CLUSTER_DIR, SUBMIT, CHECK, FETCH, POST
    RUNTIME_CONFIG = load_config(
        config_path=args.config, repository_root=REPO,
        local_profile=args.local_profile, lrz_profile=args.lrz_profile,
    )
    configure(RUNTIME_CONFIG)
    REPO = RUNTIME_CONFIG.repository_root
    CLUSTER_DIR = REPO / "scripts" / "interpro_cluster"
    SUBMIT = CLUSTER_DIR / "submit_cluster_analysis.py"
    CHECK = CLUSTER_DIR / "check_cluster_analysis.py"
    FETCH = CLUSTER_DIR / "fetch_cluster_analysis.py"
    POST = REPO / "scripts" / "run_post_interpro_for_run.py"
    rt = Roundtrip(args.run_id, RUNTIME_CONFIG)
    if not rt.run_dir.is_dir():
        raise SystemExit(f"Run folder not found: runs/{args.run_id}. "
                         "Create it with scripts/create_new_run.py first.")

    if args.dry_run:
        print(json.dumps({
            "schema_version": "1.0", "dry_run": True, "network_contacted": False,
            "run_id": args.run_id,
            "local_profile": RUNTIME_CONFIG.local_profile_name,
            "lrz_profile": RUNTIME_CONFIG.lrz_profile_name,
            "continuation_command": RUNTIME_CONFIG.command([
                ".venv/bin/edc", "cluster", "roundtrip",
                "--run-id", args.run_id,
                "--local-profile", RUNTIME_CONFIG.local_profile_name,
                "--lrz-profile", RUNTIME_CONFIG.lrz_profile_name,
            ]),
            "compatibility_command": RUNTIME_CONFIG.command([
                "python", "scripts/interpro_cluster/run_cluster_roundtrip.py",
                "--run-id", args.run_id,
                "--local-profile", RUNTIME_CONFIG.local_profile_name,
                "--lrz-profile", RUNTIME_CONFIG.lrz_profile_name,
            ]),
            "steps": ["submit", "poll", "fetch",
                      "post-interpro" if not args.no_post_interpro else "stop-after-fetch"],
        }, indent=2))
        return
    if rt.record and rt.record.read_only:
        raise SystemExit(
            "Registered legacy run is read-only; copy it before retry/resume/roundtrip.")

    if args.finalize_only:
        rt.log("=" * 66)
        rt.log(f"Re-deriving end state for run {args.run_id} (finalize only)")
        rt.finalize()
        return

    _require_cluster_profile(RUNTIME_CONFIG)

    for s in (SUBMIT, CHECK, FETCH):
        if not s.exists():
            raise SystemExit(f"Required cluster wrapper missing: {s}")

    exclusive = sum([args.submit_only, args.check_only, args.fetch_only])
    if exclusive > 1:
        raise SystemExit("Use at most one of --submit-only / --check-only / --fetch-only.")

    rt.log("=" * 66)
    rt.log(f"Cluster round-trip for run {args.run_id}")
    rt.log("Login / 2FA (if any) happen in THIS terminal; no credentials are stored.")
    rt.log("=" * 66)

    if not args.no_ssh_master:
        rt.open_ssh_master()

    try:
        if args.check_only:
            rt.run(CHECK)
            rt.log(f"cluster_analysis_status = {rt.cluster_status() or 'unknown'}")
            return
        if args.submit_only:
            rt.submit()
            rt.log("Submit-only done. Poll later with --check-only or re-run without --submit-only.")
            return
        if args.fetch_only:
            rt.fetch()
            if not args.no_post_interpro:
                rt.post_interpro()
            rt.log("Fetch-only phase done.")
            _final_message(rt, args)
            return

        rt.submit()
        rt.poll(args.poll_minutes, args.max_hours)
        rt.fetch()
        if args.no_post_interpro:
            rt.log("Skipping post-InterPro (--no-post-interpro). Run it later with "
                   f"python scripts/run_post_interpro_for_run.py --run-id {args.run_id}")
            rt.set_phase("fetched_no_post")
        else:
            rt.post_interpro()

        rt.log("=" * 66)
        rt.log("Cluster round-trip finished.")
        rt.log("=" * 66)
        _final_message(rt, args)
    finally:
        rt.close_ssh_master()


def _final_message(rt: "Roundtrip", args: argparse.Namespace) -> None:
    if args.no_post_interpro:
        rt.log("Cluster outputs fetched. Open the webapp and click "
               "\"Run Post-InterPro locally\" to finish.")
    else:
        rt.log("Results are ready. Return to the webapp and click Refresh / Open results.")


if __name__ == "__main__":
    main()
