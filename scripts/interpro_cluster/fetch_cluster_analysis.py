#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional


SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
from exondomaincompare.contracts import file_sha256, stamp_payload  # noqa: E402
from exondomaincompare.runs.registry import resolve_run_record  # noqa: E402
from exondomaincompare.config import CONTROL_PATH_ENV, load_config  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=Path(__file__).resolve().parents[2])
PROJECT_ROOT = RUNTIME_CONFIG.repository_root


def display(cmd: List[str]) -> str:
    return " ".join(shlex.quote(part) for part in cmd)


def run_cmd(cmd: List[str], *, capture: bool = False, check: bool = True) -> subprocess.CompletedProcess:
    print("+", display(cmd))
    return subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
        check=check,
    )


def read_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing JSON file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def read_json_if_exists(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, data: dict) -> None:
    run_id = str(data.get("run_id") or (
        path.parent.name if path.name == "status.json" else ""))
    data = stamp_payload(
        data,
        payload_type="status" if path.name == "status.json" else path.stem,
        run_id=run_id, dataset_id=run_id,
        profile=RUNTIME_CONFIG.public_identity(),
        generator="scripts/interpro_cluster/fetch_cluster_analysis.py",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

def ssh_mux_opts(control_path: Path) -> List[str]:
    return [
        "-o", "ControlMaster=auto",
        "-o", f"ControlPath={control_path}",
        "-o", "ControlPersist=180",
        "-o", "ConnectTimeout=30",
    ]


def close_master(control_path: Path) -> None:
    subprocess.run(
        [RUNTIME_CONFIG.executable_token("ssh"), "-O", "exit",
         "-o", f"ControlPath={control_path}", RUNTIME_CONFIG.ssh_target],
        text=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
    )


def discover_remote_files(remote_root: str, control_path: Path) -> List[str]:
    cmd = [RUNTIME_CONFIG.executable_token("ssh"), *RUNTIME_CONFIG.lrz.get("ssh_options", []),
           *ssh_mux_opts(control_path), RUNTIME_CONFIG.ssh_target,
           "find", remote_root, "-maxdepth", "6", "-type", "f", "-print"]
    result = run_cmd(cmd, capture=True, check=True)
    files = [ln.strip() for ln in (result.stdout or "").splitlines() if ln.strip()]
    return files


def scp_files_into(remote_files: List[str], local_dir: Path, control_path: Path) -> List[str]:
    if not remote_files:
        return []
    local_dir.mkdir(parents=True, exist_ok=True)
    sources = [f"{RUNTIME_CONFIG.ssh_target}:{rf}" for rf in remote_files]
    cmd = [RUNTIME_CONFIG.executable_token("scp"),
           *RUNTIME_CONFIG.lrz.get("scp_options", []),
           *ssh_mux_opts(control_path), "-p", *sources, f"{local_dir}/"]
    run_cmd(cmd, capture=False, check=True)
    return [Path(rf).name for rf in remote_files]


def classify(remote_file: str, interpro_dir: str, pytmhmm_dir: str,
             targets: Dict[str, Path]) -> Optional[Path]:
    name = Path(remote_file).name

    def under(base: str, sub: str) -> bool:
        return remote_file.startswith(f"{base}/{sub}/") or remote_file == f"{base}/{sub}"

    # InterProScan
    if under(interpro_dir, "output"):
        return targets["interpro_output"]
    if remote_file.startswith(interpro_dir + "/"):
        if name.endswith(".sbatch"):
            return targets["interpro_slurm"]
        if name.endswith((".out", ".err", ".log")):
            return targets["interpro_logs"]
        return None
    # pyTMHMM
    if under(pytmhmm_dir, "output"):
        return targets["pytmhmm_output"]
    if remote_file.startswith(pytmhmm_dir + "/"):
        if name.endswith(".sbatch"):
            return targets["pytmhmm_slurm"]
        if name.endswith((".out", ".err", ".log")):
            return targets["pytmhmm_logs"]
        return None
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch InterProScan and pyTMHMM outputs from LRZ.")
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
    run_id = args.run_id
    record = resolve_run_record(RUNTIME_CONFIG, run_id)
    run_dir = record.path if record else RUNTIME_CONFIG.runs_root / run_id
    if record and record.read_only and not args.dry_run:
        raise RuntimeError("Registered legacy run is read-only; copy it before fetch.")
    status_path = run_dir / "status.json"
    status = read_json(status_path)

    cluster_jobs = status.get("cluster_jobs", {})
    stored_profile = str(cluster_jobs.get("profile") or status.get("cluster_profile") or "")
    if stored_profile and stored_profile != RUNTIME_CONFIG.lrz_profile_name:
        raise RuntimeError(
            f"Pending run requires LRZ profile {stored_profile!r}, "
            f"not {RUNTIME_CONFIG.lrz_profile_name!r}."
        )
    remote_root = cluster_jobs.get("remote_dir") or RUNTIME_CONFIG.remote_run_root(run_id)
    remote_interpro_dir = cluster_jobs.get("remote_interproscan_dir") or f"{remote_root}/interproscan"
    remote_pytmhmm_dir = cluster_jobs.get("remote_pytmhmm_dir") or f"{remote_root}/pytmhmm"
    interpro_job_id = cluster_jobs.get("interproscan_job_id")
    pytmhmm_job_id = cluster_jobs.get("pytmhmm_job_id")

    if not remote_interpro_dir or not remote_pytmhmm_dir:
        raise RuntimeError(
            "Remote directories are missing in status.json. Did you run submit_cluster_analysis.py?")
    if not remote_root:
        remote_root = str(Path(remote_interpro_dir).parent)

    if args.dry_run:
        print(json.dumps({
            "schema_version": "1.0", "dry_run": True, "network_contacted": False,
            "run_id": run_id, "local_profile": RUNTIME_CONFIG.local_profile_name,
            "lrz_profile": RUNTIME_CONFIG.lrz_profile_name,
            "remote_run": f"profile:{RUNTIME_CONFIG.lrz_profile_name}:{run_id}",
            "steps": ["discover remote outputs", "fetch exact files", "verify checksums"],
        }, indent=2))
        return

    interpro_base = run_dir / "results" / "14_interproscan" / "primary"
    pytmhmm_base = run_dir / "results" / "15_exon_domain_boundary_post_interpro" / "pytmhmm_primary"

    targets: Dict[str, Path] = {
        "interpro_output": interpro_base / "output",
        "interpro_logs": interpro_base / "logs",
        "interpro_slurm": interpro_base / "slurm",
        "pytmhmm_output": pytmhmm_base / "output",
        "pytmhmm_logs": pytmhmm_base / "logs",
        "pytmhmm_slurm": pytmhmm_base / "slurm",
    }
    for d in targets.values():
        d.mkdir(parents=True, exist_ok=True)

    interpro_manifest_path = interpro_base / "interproscan_manifest.json"
    pytmhmm_manifest_path = pytmhmm_base / "pytmhmm_manifest.json"

    shared_cp = os.environ.get(CONTROL_PATH_ENV, "").strip()
    owns_master = not shared_cp
    _private_sock = Path("/tmp") / f"fgfr2_{hashlib.sha1(run_id.encode()).hexdigest()[:8]}_cm.sock"
    try:
        control_path = Path(shared_cp) if shared_cp else _private_sock

        print(f"Discovering remote files under: {remote_root}")
        try:
            remote_files = discover_remote_files(remote_root, control_path)
        except subprocess.CalledProcessError as exc:
            stderr = (exc.stderr or "").strip()
            raise RuntimeError(
                f"Remote discovery failed (ssh find). stderr:\n{stderr}") from exc

        print(f"Discovered {len(remote_files)} remote file(s).")

        grouped: Dict[str, List[str]] = {key: [] for key in targets}
        skipped: List[str] = []
        for rf in remote_files:
            dest = classify(rf, remote_interpro_dir, remote_pytmhmm_dir, targets)
            if dest is None:
                skipped.append(rf)
                continue
            key = next(k for k, v in targets.items() if v == dest)
            grouped[key].append(rf)

        for key, files in grouped.items():
            if not files:
                continue
            print(f"Fetching {len(files)} file(s) -> {targets[key]}")
            try:
                scp_files_into(files, targets[key], control_path)
            except subprocess.CalledProcessError as exc:
                print(f"WARNING: scp into {targets[key]} failed: {exc}")

        if skipped:
            print(f"Skipped {len(skipped)} unrelated remote file(s) (not part of the result set).")
    finally:
        if owns_master:
            close_master(control_path)
            try:
                if control_path.exists():
                    control_path.unlink()
            except OSError:
                pass

    now = datetime.now().isoformat(timespec="seconds")

    interpro_output_dir = targets["interpro_output"]
    interpro_logs_dir = targets["interpro_logs"]
    pytmhmm_output_dir = targets["pytmhmm_output"]
    pytmhmm_logs_dir = targets["pytmhmm_logs"]

    interpro_files = sorted(p.name for p in interpro_output_dir.iterdir() if p.is_file())
    pytmhmm_files = sorted(p.name for p in pytmhmm_output_dir.iterdir() if p.is_file())

    old_interpro_manifest = read_json_if_exists(interpro_manifest_path)
    old_pytmhmm_manifest = read_json_if_exists(pytmhmm_manifest_path)

    interpro_manifest = {
        **old_interpro_manifest,
        "run_id": run_id,
        "tool": "InterProScan",
        "job_id": interpro_job_id,
        "remote_dir": f"profile:{RUNTIME_CONFIG.lrz_profile_name}:{run_id}/interproscan",
        "fetched_at": now,
        "local_output_dir": "run:results/14_interproscan/primary/output",
        "local_logs_dir": "run:results/14_interproscan/primary/logs",
        "files": interpro_files,
        "file_sha256": {p.name: file_sha256(p) for p in interpro_output_dir.iterdir()
                        if p.is_file()},
    }
    pytmhmm_manifest = {
        **old_pytmhmm_manifest,
        "run_id": run_id,
        "tool": "pyTMHMM",
        "job_id": pytmhmm_job_id,
        "remote_dir": f"profile:{RUNTIME_CONFIG.lrz_profile_name}:{run_id}/pytmhmm",
        "fetched_at": now,
        "local_output_dir": "run:results/15_exon_domain_boundary_post_interpro/"
                            "pytmhmm_primary/output",
        "local_logs_dir": "run:results/15_exon_domain_boundary_post_interpro/"
                          "pytmhmm_primary/logs",
        "files": pytmhmm_files,
        "file_sha256": {p.name: file_sha256(p) for p in pytmhmm_output_dir.iterdir()
                        if p.is_file()},
    }
    write_json(interpro_manifest_path, interpro_manifest)
    write_json(pytmhmm_manifest_path, pytmhmm_manifest)

    expected_interpro = [
        interpro_output_dir / "input.fasta.tsv",
        interpro_output_dir / "input.fasta.gff3",
        interpro_output_dir / "input.fasta.json",
    ]
    expected_pytmhmm = [
        pytmhmm_output_dir / "pytmhmm_summary_all.tsv",
        pytmhmm_output_dir / "pytmhmm_transmembrane_hits.tsv",
    ]
    interpro_ok = all(p.exists() and p.stat().st_size > 0 for p in expected_interpro)
    pytmhmm_ok = all(p.exists() and p.stat().st_size > 0 for p in expected_pytmhmm)

    status["cluster_analysis_status"] = "fetched_complete" if interpro_ok and pytmhmm_ok else "fetched_incomplete"
    status["cluster_fetch_status"] = "complete" if interpro_ok and pytmhmm_ok else "incomplete"
    status["fetched_at"] = now
    status["updated_at"] = now
    status["cluster_fetch_summary"] = {
        "interproscan_ok": interpro_ok,
        "pytmhmm_ok": pytmhmm_ok,
        "interproscan_files": interpro_files,
        "pytmhmm_files_count": len(pytmhmm_files),
        "pytmhmm_summary_exists": (pytmhmm_output_dir / "pytmhmm_summary_all.tsv").exists(),
        "pytmhmm_transmembrane_hits_exists": (pytmhmm_output_dir / "pytmhmm_transmembrane_hits.tsv").exists(),
    }
    write_json(status_path, status)

    print()
    print("Fetch finished.")
    print(f"InterProScan OK: {interpro_ok}")
    print(f"pyTMHMM OK: {pytmhmm_ok}")
    print(f"cluster_fetch_status: {status['cluster_fetch_status']}")
    print()
    print("InterProScan output:", interpro_output_dir)
    print("pyTMHMM output:", pytmhmm_output_dir)

    if not interpro_ok:
        print("\nWARNING: InterProScan expected output files are missing or empty.")
        print("Check logs here:", interpro_logs_dir)
    if not pytmhmm_ok:
        print("\nWARNING: pyTMHMM expected summary files are missing or empty.")
        print("Check logs here:", pytmhmm_logs_dir)


if __name__ == "__main__":
    main()
