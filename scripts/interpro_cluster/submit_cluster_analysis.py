#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from ssh_common import clean_ssh_output, configure, scp_cmd, ssh_cmd, ssh_target

SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
from framework.data_contract import (  # noqa: E402
    file_sha256, resolve_path_reference, stamp_payload,
)
from framework.local_registry import resolve_run_record  # noqa: E402
from framework.portable_config import (  # noqa: E402
    ConfigurationError, load_config, remote_shell_path,
)

RUNTIME_CONFIG = load_config(repository_root=Path(__file__).resolve().parents[2])
PROJECT_ROOT = RUNTIME_CONFIG.repository_root


def run_cmd(cmd: list[str], *, cwd: Path | None = None, capture: bool = False) -> str:
    print("+", " ".join(cmd))
    completed = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
        check=True,
    )
    return completed.stdout.strip() if capture and completed.stdout else ""


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def read_json(path: Path) -> dict:
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
        generator="scripts/interpro_cluster/submit_cluster_analysis.py",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def safe_run_id(run_id: str) -> str:
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")
    cleaned = "".join(ch if ch in allowed else "_" for ch in run_id)
    if not cleaned:
        raise ValueError("run_id is empty after sanitizing.")
    return cleaned


def remote_shell_quote(path: str) -> str:
    return remote_shell_path(path)


def submit_remote(remote_dir: str, slurm_path: str) -> str:
    cmd = f"cd {remote_shell_quote(remote_dir)} && sbatch {remote_shell_quote(slurm_path)}"
    output = clean_ssh_output(run_cmd(ssh_cmd(cmd), capture=True))
    print(output)

    # Expected: Submitted batch job 371278
    job_id = None
    for token in output.split():
        if token.isdigit():
            job_id = token
    if not job_id:
        raise RuntimeError(f"Could not parse job id from sbatch output: {output}")
    return job_id


def main() -> None:
    parser = argparse.ArgumentParser(description="Submit InterProScan and pyTMHMM jobs for one local run.")
    parser.add_argument("--run-id", required=True, help="Run ID, e.g. 2026-07-06_1530_fgfr2_custom")
    parser.add_argument("--force", action="store_true", help="Allow resubmission even if job IDs already exist.")
    parser.add_argument("--config")
    parser.add_argument("--local-profile")
    parser.add_argument("--lrz-profile")
    parser.add_argument("--dry-run", action="store_true",
                        help="validate and render a plan without filesystem or network changes")
    args = parser.parse_args()

    global RUNTIME_CONFIG, PROJECT_ROOT
    RUNTIME_CONFIG = load_config(
        config_path=args.config, repository_root=PROJECT_ROOT,
        local_profile=args.local_profile, lrz_profile=args.lrz_profile,
    )
    PROJECT_ROOT = RUNTIME_CONFIG.repository_root
    configure(RUNTIME_CONFIG)
    try:
        RUNTIME_CONFIG.require_cluster()
    except ConfigurationError as exc:
        raise SystemExit(f"Cluster profile is not ready: {exc}") from None

    run_id = safe_run_id(args.run_id)

    record = resolve_run_record(RUNTIME_CONFIG, run_id)
    run_dir = record.path if record else RUNTIME_CONFIG.runs_root / run_id
    if record and record.read_only and not args.dry_run:
        raise SystemExit("Registered legacy run is read-only; copy it before submission.")
    status_path = run_dir / "status.json"
    if args.dry_run:
        plan = {
            "schema_version": "1.0",
            "dry_run": True,
            "network_contacted": False,
            "run_id": run_id,
            "local_profile": RUNTIME_CONFIG.local_profile_name,
            "lrz_profile": RUNTIME_CONFIG.lrz_profile_name,
            "input": "run:results/14_interproscan/primary/input/"
                     "final_pre_interpro_proteins_primary.faa",
            "remote_run": f"profile:{RUNTIME_CONFIG.lrz_profile_name}:{run_id}",
            "steps": ["validate input", "prepare Slurm scripts", "upload", "submit"],
        }
        print(json.dumps(plan, indent=2))
        return

    interpro_base = run_dir / "results" / "14_interproscan" / "primary"
    interpro_input_dir = interpro_base / "input"
    interpro_output_dir = interpro_base / "output"
    interpro_logs_dir = interpro_base / "logs"
    interpro_slurm_dir = interpro_base / "slurm"
    interpro_manifest_path = interpro_base / "interproscan_manifest.json"

    pytmhmm_base = run_dir / "results" / "15_exon_domain_boundary_post_interpro" / "pytmhmm_primary"
    pytmhmm_input_dir = pytmhmm_base / "input"
    pytmhmm_output_dir = pytmhmm_base / "output"
    pytmhmm_logs_dir = pytmhmm_base / "logs"
    pytmhmm_slurm_dir = pytmhmm_base / "slurm"
    pytmhmm_manifest_path = pytmhmm_base / "pytmhmm_manifest.json"

    for d in [
        interpro_input_dir,
        interpro_output_dir,
        interpro_logs_dir,
        interpro_slurm_dir,
        pytmhmm_input_dir,
        pytmhmm_output_dir,
        pytmhmm_logs_dir,
        pytmhmm_slurm_dir,
        run_dir / "logs",
        run_dir / "slurm",
    ]:
        d.mkdir(parents=True, exist_ok=True)

    interpro_fasta = interpro_input_dir / "final_pre_interpro_proteins_primary.faa"
    freeze_fasta = run_dir / "results" / "13_final_pre_interpro_closure" / "freeze" / "final_pre_interpro_proteins_primary.faa"

    # Resolve the source FASTA in a gene-agnostic way. FGFR2 runs keep using the
    # freeze path; generic core-only runs may point run_config at a different
    # primary FASTA via `cluster_input_fasta` / `primary_fasta_path`.
    def _resolve_source_fasta() -> Path | None:
        run_config = read_json(run_dir / "run_config.json")
        for key in ("cluster_input_fasta", "primary_fasta_path", "primary_fasta_expected_path"):
            rel = run_config.get(key)
            if rel:
                cand = resolve_path_reference(
                    str(rel), repository_root=PROJECT_ROOT, run_root=run_dir)
                if cand.exists():
                    return cand
        if freeze_fasta.exists():
            return freeze_fasta
        return None

    if not interpro_fasta.exists():
        source_fasta = _resolve_source_fasta()
        if source_fasta is not None:
            print(f"InterProScan input FASTA not found. Copying from:\n  {source_fasta}\n→ {interpro_fasta}")
            shutil.copy2(source_fasta, interpro_fasta)
        else:
            raise FileNotFoundError(
                "Input FASTA not found.\n"
                f"Expected either:\n  {interpro_fasta}\n"
                f"or run_config cluster_input_fasta / primary_fasta_path,\n"
                f"or the freeze path:\n  {freeze_fasta}"
            )

    if interpro_fasta.stat().st_size == 0:
        raise RuntimeError(f"Input FASTA is empty: {interpro_fasta}")

    pytmhmm_fasta = pytmhmm_input_dir / "final_pre_interpro_proteins_primary.faa"
    shutil.copy2(interpro_fasta, pytmhmm_fasta)

    status = read_json(status_path)
    cluster_jobs = status.get("cluster_jobs", {})

    if cluster_jobs and not args.force:
        existing = {
            k: v for k, v in cluster_jobs.items()
            if k in {"interproscan_job_id", "pytmhmm_job_id"} and v
        }
        if existing:
            raise RuntimeError(
                "This run already has cluster job IDs in status.json. "
                "Use --force only if you intentionally want to resubmit.\n"
                f"Existing jobs: {existing}"
            )

    remote_dir = RUNTIME_CONFIG.remote_run_root(run_id)
    remote_interpro_dir = f"{remote_dir}/interproscan"
    remote_pytmhmm_dir = f"{remote_dir}/pytmhmm"

    print("Creating remote directories...")
    run_cmd(ssh_cmd(
        f"mkdir -p "
        f"{remote_shell_quote(remote_interpro_dir + '/output')} "
        f"{remote_shell_quote(remote_interpro_dir + '/temp')} "
        f"{remote_shell_quote(remote_pytmhmm_dir + '/output')}"
    ))

    print("Uploading FASTA...")
    run_cmd(scp_cmd([str(interpro_fasta), f"{ssh_target()}:{remote_dir}/input.fasta"]))

    interpro = RUNTIME_CONFIG.lrz.get("interproscan", {})
    pytmhmm = RUNTIME_CONFIG.lrz.get("pytmhmm", {})
    slurm_partition = str(RUNTIME_CONFIG.lrz.get("partition"))
    slurm_account = str(RUNTIME_CONFIG.lrz.get("account") or "")
    account_line = f"#SBATCH --account={slurm_account}\n" if slurm_account else ""
    interpro_launcher = str(interpro.get("launcher"))
    interpro_env = str(interpro.get("environment") or "")
    interpro_module = str(interpro.get("module") or "")
    interpro_module_line = f"module load {remote_shell_quote(interpro_module)}\n" if interpro_module else ""
    interpro_env_line = f'export PATH="{interpro_env}/bin:$PATH"\n' if interpro_env else ""
    pytmhmm_launcher = str(pytmhmm.get("launcher"))
    pytmhmm_env = str(pytmhmm.get("environment") or "")
    pytmhmm_module = str(pytmhmm.get("module") or "")
    pytmhmm_module_line = f"module load {remote_shell_quote(pytmhmm_module)}\n" if pytmhmm_module else ""
    pytmhmm_env_line = f'export PATH="{pytmhmm_env}/bin:$PATH"\n' if pytmhmm_env else ""
    remote_python = str(pytmhmm.get("python") or "python")

    interpro_slurm = f"""#!/bin/bash
#SBATCH --job-name=interpro_{run_id[:20]}
#SBATCH --partition={slurm_partition}
{account_line}#SBATCH --output=interpro_%j.out
#SBATCH --error=interpro_%j.err
#SBATCH --time={interpro.get("time", "08:00:00")}
#SBATCH --cpus-per-task={interpro.get("cpus", 4)}
#SBATCH --mem={interpro.get("memory", "24G")}

set -euo pipefail

{interpro_module_line}{interpro_env_line}
export PYTHONNOUSERSITE=1

cd "$SLURM_SUBMIT_DIR"

mkdir -p output temp

echo "Working directory:"
pwd

echo "Input FASTA:"
ls -lh ../input.fasta

echo "InterProScan version:"
{remote_shell_quote(interpro_launcher)} --version

echo "Starting InterProScan..."

{remote_shell_quote(interpro_launcher)} \\
  -i ../input.fasta \\
  -f TSV,GFF3,JSON \\
  -d output \\
  -cpu "$SLURM_CPUS_PER_TASK" \\
  -goterms \\
  -pa \\
  -dp \\
  -T "$PWD/temp"

echo "InterProScan finished."
echo "Output files:"
find output -maxdepth 3 -type f -print -exec ls -lh {{}} \\;
"""

    pytmhmm_slurm = f"""#!/bin/bash
#SBATCH --job-name=pytmhmm_{run_id[:20]}
#SBATCH --partition={slurm_partition}
{account_line}#SBATCH --output=pytmhmm_%j.out
#SBATCH --error=pytmhmm_%j.err
#SBATCH --time={pytmhmm.get("time", "02:00:00")}
#SBATCH --cpus-per-task={pytmhmm.get("cpus", 1)}
#SBATCH --mem={pytmhmm.get("memory", "8G")}

set -euo pipefail

{pytmhmm_module_line}{pytmhmm_env_line}
export PYTHONNOUSERSITE=1

cd "$SLURM_SUBMIT_DIR"

mkdir -p output

echo "Working directory:"
pwd

echo "Input FASTA:"
ls -lh ../input.fasta

echo "pyTMHMM:"
command -v {remote_shell_quote(pytmhmm_launcher)}
{remote_shell_quote(pytmhmm_launcher)} -h | head -40 || true

cd output
rm -f *

echo "Starting pyTMHMM..."
{remote_shell_quote(pytmhmm_launcher)} -f ../../input.fasta

echo "pyTMHMM finished."

echo "Creating pyTMHMM combined summary tables..."

{remote_shell_quote(remote_python)} - <<'PYINNER'
from pathlib import Path
import re

outdir = Path(".")
summary_files = sorted(outdir.glob("*.summary"))

with open("pytmhmm_summary_all.tsv", "w", encoding="utf-8") as w:
    w.write("sequence_id\\tline\\n")
    for f in summary_files:
        seq_id = f.name[:-len(".summary")]
        for line in f.read_text(errors="replace").splitlines():
            line = line.strip()
            if line:
                w.write(f"{{seq_id}}\\t{{line}}\\n")

with open("pytmhmm_transmembrane_hits.tsv", "w", encoding="utf-8") as w:
    w.write("sequence_id\\tline\\n")
    for f in summary_files:
        seq_id = f.name[:-len(".summary")]
        for line in f.read_text(errors="replace").splitlines():
            line_clean = line.strip()
            if re.search(r"transmembrane|tmhelix|inside|outside|TM", line_clean, re.IGNORECASE):
                w.write(f"{{seq_id}}\\t{{line_clean}}\\n")

print("summary files:", len(summary_files))
print("created: pytmhmm_summary_all.tsv")
print("created: pytmhmm_transmembrane_hits.tsv")
PYINNER

echo "Output files:"
find . -maxdepth 3 -type f | sort
"""

    local_interpro_slurm = interpro_slurm_dir / f"interproscan_{run_id}.sbatch"
    local_pytmhmm_slurm = pytmhmm_slurm_dir / f"pytmhmm_{run_id}.sbatch"

    write_text(local_interpro_slurm, interpro_slurm)
    write_text(local_pytmhmm_slurm, pytmhmm_slurm)

    tmp_interpro = PROJECT_ROOT / ".tmp_interproscan_submit.sbatch"
    tmp_pytmhmm = PROJECT_ROOT / ".tmp_pytmhmm_submit.sbatch"
    write_text(tmp_interpro, interpro_slurm)
    write_text(tmp_pytmhmm, pytmhmm_slurm)

    print("Uploading SLURM scripts...")
    run_cmd(scp_cmd([str(tmp_interpro), f"{ssh_target()}:{remote_interpro_dir}/run_interproscan.sbatch"]))
    run_cmd(scp_cmd([str(tmp_pytmhmm), f"{ssh_target()}:{remote_pytmhmm_dir}/run_pytmhmm.sbatch"]))

    try:
        tmp_interpro.unlink(missing_ok=True)
        tmp_pytmhmm.unlink(missing_ok=True)
    except Exception:
        pass

    print("Submitting InterProScan job...")
    interpro_job_id = submit_remote(remote_interpro_dir, "run_interproscan.sbatch")

    print("Submitting pyTMHMM job...")
    pytmhmm_job_id = submit_remote(remote_pytmhmm_dir, "run_pytmhmm.sbatch")

    now = datetime.now().isoformat(timespec="seconds")

    continuation = RUNTIME_CONFIG.command([
        ".venv/bin/edc", "cluster", "roundtrip",
        "--run-id", run_id,
        "--local-profile", RUNTIME_CONFIG.local_profile_name,
        "--lrz-profile", RUNTIME_CONFIG.lrz_profile_name,
    ])
    compatibility_continuation = RUNTIME_CONFIG.command([
        "python", "scripts/interpro_cluster/run_cluster_roundtrip.py",
        "--run-id", run_id,
        "--local-profile", RUNTIME_CONFIG.local_profile_name,
        "--lrz-profile", RUNTIME_CONFIG.lrz_profile_name,
    ])
    input_sha256 = file_sha256(interpro_fasta)
    cluster_jobs = {
        "profile": RUNTIME_CONFIG.lrz_profile_name,
        "remote_run_ref": f"profile:{RUNTIME_CONFIG.lrz_profile_name}:{run_id}",
        "interproscan_job_id": interpro_job_id,
        "pytmhmm_job_id": pytmhmm_job_id,
        "submitted_at": now,
        "input_sha256": input_sha256,
        "continuation_command": continuation,
        "compatibility_command": compatibility_continuation,
    }

    status.update({
        "run_id": run_id,
        "cluster_analysis_status": "submitted",
        "cluster_jobs": cluster_jobs,
        "cluster_profile": RUNTIME_CONFIG.lrz_profile_name,
        "cluster_continuation_command": continuation,
        "cluster_continuation_compatibility_command": compatibility_continuation,
        "updated_at": now,
    })
    write_json(status_path, status)

    write_json(interpro_manifest_path, {
        "run_id": run_id,
        "tool": "InterProScan",
        "tool_profile_ref": f"profile:{RUNTIME_CONFIG.lrz_profile_name}:interproscan",
        "remote_dir": f"profile:{RUNTIME_CONFIG.lrz_profile_name}:{run_id}/interproscan",
        "job_id": interpro_job_id,
        "input_fasta_local": "run:results/14_interproscan/primary/input/"
                             "final_pre_interpro_proteins_primary.faa",
        "input_sha256": input_sha256,
        "submitted_at": now,
        "expected_outputs": [
            "input.fasta.tsv",
            "input.fasta.gff3",
            "input.fasta.json",
        ],
    })

    write_json(pytmhmm_manifest_path, {
        "run_id": run_id,
        "tool": "pyTMHMM",
        "tool_profile_ref": f"profile:{RUNTIME_CONFIG.lrz_profile_name}:pytmhmm",
        "remote_dir": f"profile:{RUNTIME_CONFIG.lrz_profile_name}:{run_id}/pytmhmm",
        "job_id": pytmhmm_job_id,
        "input_fasta_local": "run:results/15_exon_domain_boundary_post_interpro/"
                             "pytmhmm_primary/input/final_pre_interpro_proteins_primary.faa",
        "input_sha256": input_sha256,
        "submitted_at": now,
        "expected_outputs": [
            "pytmhmm_summary_all.tsv",
            "pytmhmm_transmembrane_hits.tsv",
            "*.summary",
            "*.annotation",
            "*.plot",
        ],
    })

    print()
    print("Submitted successfully.")
    print(f"Run ID: {run_id}")
    print(f"Remote profile: {RUNTIME_CONFIG.lrz_profile_name}")
    print(f"InterProScan job ID: {interpro_job_id}")
    print(f"pyTMHMM job ID: {pytmhmm_job_id}")
    print()
    print("Check status with:")
    print(f"python scripts/interpro_cluster/check_cluster_analysis.py --run-id {run_id}")


if __name__ == "__main__":
    main()
