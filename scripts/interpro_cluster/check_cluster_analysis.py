#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from exondomaincompare.cluster.ssh_common import clean_ssh_output, configure, ssh_cmd

SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
from exondomaincompare.contracts import stamp_payload  # noqa: E402
from exondomaincompare.runs.registry import resolve_run_record  # noqa: E402
from exondomaincompare.config import load_config  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=Path(__file__).resolve().parents[2])
PROJECT_ROOT = RUNTIME_CONFIG.repository_root


def run_cmd(cmd: list[str], capture: bool = True, check: bool = False) -> str:
    print("+", " ".join(cmd))
    completed = subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
        check=check,
    )
    return completed.stdout.strip() if capture and completed.stdout else ""


def ssh_run(remote_command: str) -> str:
    return clean_ssh_output(run_cmd(ssh_cmd(remote_command), capture=True))


def _exit_ok(exit_code: str) -> bool:
    code = (exit_code or "").strip()
    if not code:
        return True
    first = code.split(":", 1)[0]
    try:
        return int(first) == 0
    except ValueError:
        return True


def read_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing JSON file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, data: dict) -> None:
    run_id = str(data.get("run_id") or path.parent.name)
    data = stamp_payload(
        data, payload_type="status", run_id=run_id, dataset_id=run_id,
        profile=RUNTIME_CONFIG.public_identity(),
        generator="scripts/interpro_cluster/check_cluster_analysis.py",
    )
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _squeue_row(job_id: str, cleaned: str):
    for line in cleaned.splitlines():
        parts = line.split("|")
        if parts and parts[0].strip() == str(job_id):
            return parts
    return None


def check_job(job_id: str) -> dict:
    squeue_cmd = f"squeue -j {job_id} -h -o '%i|%T|%M|%R|%j'"
    squeue_clean = ssh_run(squeue_cmd)
    row = _squeue_row(job_id, squeue_clean)
    if row is not None:
        return {
            "job_id": job_id,
            "source": "squeue",
            "state": row[1] if len(row) > 1 else "UNKNOWN",
            "elapsed": row[2] if len(row) > 2 else "",
            "reason_or_node": row[3] if len(row) > 3 else "",
            "job_name": row[4] if len(row) > 4 else "",
            "raw": squeue_clean,
        }

    sacct_cmd = f"sacct -j {job_id} --format=JobID,JobName,State,Elapsed,ExitCode -P -n"
    sacct_clean = ssh_run(sacct_cmd)

    if sacct_clean:
        lines = sacct_clean.splitlines()
        main_line = None
        for line in lines:
            first = line.split("|", 1)[0].strip()
            if first == str(job_id):
                main_line = line
                break
        if main_line is not None:
            parts = main_line.split("|")
            return {
                "job_id": job_id,
                "source": "sacct",
                "state": parts[2] if len(parts) > 2 else "UNKNOWN",
                "elapsed": parts[3] if len(parts) > 3 else "",
                "exit_code": parts[4] if len(parts) > 4 else "",
                "job_name": parts[1] if len(parts) > 1 else "",
                "raw": sacct_clean,
            }

    return {
        "job_id": job_id,
        "source": "none",
        "state": "UNKNOWN_NOT_IN_SQUEUE_OR_SACCT",
        "raw": (sacct_clean or squeue_clean or ""),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Check InterProScan and pyTMHMM SLURM job status.")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--config")
    parser.add_argument("--local-profile")
    parser.add_argument("--lrz-profile")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    global RUNTIME_CONFIG, PROJECT_ROOT
    RUNTIME_CONFIG = load_config(
        config_path=args.config, repository_root=PROJECT_ROOT,
        local_profile=args.local_profile, lrz_profile=args.lrz_profile,
    )
    PROJECT_ROOT = RUNTIME_CONFIG.repository_root
    configure(RUNTIME_CONFIG)
    if args.dry_run:
        print(json.dumps({
            "schema_version": "1.0", "dry_run": True, "network_contacted": False,
            "run_id": args.run_id, "local_profile": RUNTIME_CONFIG.local_profile_name,
            "lrz_profile": RUNTIME_CONFIG.lrz_profile_name,
            "steps": ["read job identities", "query Slurm status"],
        }, indent=2))
        return

    record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    run_dir = record.path if record else RUNTIME_CONFIG.runs_root / args.run_id
    status_path = run_dir / "status.json"
    status = read_json(status_path)

    cluster_jobs = status.get("cluster_jobs", {})
    interpro_job_id = cluster_jobs.get("interproscan_job_id")
    pytmhmm_job_id = cluster_jobs.get("pytmhmm_job_id")

    if not interpro_job_id and not pytmhmm_job_id:
        raise RuntimeError("No cluster job IDs found in status.json.")

    result = {}

    print()
    print(f"Run ID: {args.run_id}")
    print("=" * 70)

    if interpro_job_id:
        print()
        print("InterProScan")
        print("-" * 70)
        interpro_status = check_job(interpro_job_id)
        result["interproscan"] = interpro_status
        for k, v in interpro_status.items():
            if k != "raw":
                print(f"{k}: {v}")

    if pytmhmm_job_id:
        print()
        print("pyTMHMM")
        print("-" * 70)
        pytmhmm_status = check_job(pytmhmm_job_id)
        result["pytmhmm"] = pytmhmm_status
        for k, v in pytmhmm_status.items():
            if k != "raw":
                print(f"{k}: {v}")

    print()
    print("=" * 70)

    states = [v.get("state", "") for v in result.values()]
    running_like = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "REQUEUED", "RESIZING"}
    failed_like = {"FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "BOOT_FAIL", "DEADLINE"}

    def is_running(s: str) -> bool:
        return any(tok in s for tok in running_like)

    def is_failed(s: str) -> bool:
        return any(tok in s for tok in failed_like)

    completed_with_bad_exit = any(
        "COMPLETED" in v.get("state", "") and not _exit_ok(v.get("exit_code", ""))
        for v in result.values()
    )

    if any(is_running(s) for s in states):
        overall = "running"
        print("Overall: at least one job is still running or pending.")
    elif any(is_failed(s) for s in states) or completed_with_bad_exit:
        overall = "error"
        print("Overall: at least one job appears to have failed. Check logs after fetching or on the cluster.")
    elif states and all("COMPLETED" in s for s in states):
        overall = "completed"
        print("Overall: jobs completed successfully. You can fetch results now.")
    else:
        overall = "unknown"
        print("Overall: status is unclear. This can happen if sacct is delayed or unavailable.")

    status["cluster_analysis_status"] = overall
    status["cluster_status_detail"] = result
    status["updated_at"] = datetime.now().isoformat(timespec="seconds")
    write_json(status_path, status)

    print()
    print("Next command when completed:")
    print(f"python scripts/interpro_cluster/fetch_cluster_analysis.py --run-id {args.run_id}")


if __name__ == "__main__":
    main()
