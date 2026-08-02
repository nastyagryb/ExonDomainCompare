#!/usr/bin/env python3
"""Check the files intended for the public repository."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tomllib
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
ALLOWED_TOP_LEVEL = {
    ".gitattributes", ".github", ".gitignore", "CHANGELOG.md", "CITATION.cff",
    "LICENSE", "Makefile", "README.md", "SECURITY.md", "config", "configs",
    "datasets", "docs", "pyproject.toml", "reference", "references",
    "requirements", "requirements.txt", "run_fgfr2_pipeline_current_final_pre_interpro.sh",
    "run_fgfr2_pipeline_current_v3.sh", "scripts", "species_list_final_30.txt",
    "src", "tests", "webapp",
}
FORBIDDEN_TOP_LEVEL = {
    ".venv", "venv", "runs", "results", "artifacts", "tmp", "cache", "logs",
    "packages",
}
FORBIDDEN_PARTS = {"node_modules", "__pycache__", ".pytest_cache"}


def release_files() -> list[Path]:
    process = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
        cwd=ROOT, text=True, capture_output=True, check=True,
    )
    return [Path(line) for line in process.stdout.splitlines() if line]


def check_dataset_hashes() -> None:
    sums = ROOT / "datasets" / "SHA256SUMS"
    for line in sums.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        path = ROOT / "datasets" / relative
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError(f"Dataset checksum mismatch: {relative}")


def main() -> int:
    files = release_files()
    if not files:
        raise RuntimeError("No release files found")
    for required in ("README.md", "LICENSE", "CITATION.cff", "CHANGELOG.md",
                     "datasets/DATA_NOTICE.md", ".github/workflows/ci.yml"):
        if Path(required) not in files:
            raise RuntimeError(f"Required release file missing: {required}")
    for relative in files:
        if relative.parts[0] not in ALLOWED_TOP_LEVEL:
            raise RuntimeError(f"Unexpected top-level path: {relative}")
        if relative.parts[0] in FORBIDDEN_TOP_LEVEL or \
                FORBIDDEN_PARTS.intersection(relative.parts):
            raise RuntimeError(f"Private or generated path included: {relative}")
        path = ROOT / relative
        if path.is_file() and path.stat().st_size >= 100_000_000:
            raise RuntimeError(f"File is too large for GitHub: {relative}")
    if Path("RELEASE_CANDIDATE_STATUS.md") in files:
        raise RuntimeError("Internal release status must not be published")
    tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    yaml.safe_load((ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8"))
    json.loads((ROOT / "datasets/registry.json").read_text(encoding="utf-8"))
    check_dataset_hashes()
    print(f"Public release check passed: {len(files)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
