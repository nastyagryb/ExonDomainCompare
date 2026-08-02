"""Portable projections for website and download payloads.

Raw run records deliberately retain the paths used by the analysis runtime.
Anything exposed by the website or its download catalogue passes through this
module so those implementation details never disclose a user's machine.
"""
from __future__ import annotations

import json
import re
import shlex
from pathlib import Path
from typing import Any

from exondomaincompare.config import discover_repository_root


REPO_ROOT = discover_repository_root(__file__)
_SEP = "/"
_PERSONAL_ABSOLUTE_RE = re.compile(
    rf"(?:{_SEP}Users{_SEP}[^{_SEP}\s]+{_SEP}|"
    rf"{_SEP}home{_SEP}[^{_SEP}\s]+{_SEP}|"
    rf"{_SEP}mnt{_SEP}data{_SEP}|"
    r"[A-Za-z]:\\Users\\[^\\\s]+\\)"
)
_WINDOWS_ABSOLUTE_RE = re.compile(r"^[A-Za-z]:[\\/]")
_ENV_ASSIGNMENT_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")
_PATH_KEYS = {
    "annotation_file",
    "cwd",
    "executable",
    "gff",
    "genomic_gff",
    "outdir",
    "project_root",
    "protein_faa",
    "python_executable",
    "repository_root",
    "root",
    "synteny_source_file",
    "source_file",
    "path",
    "venv",
    "virtual_environment",
}


def contains_personal_absolute_path(value: str) -> bool:
    return bool(_PERSONAL_ABSOLUTE_RE.search(value))


def is_machine_absolute_path(value: str) -> bool:
    """Recognise POSIX and Windows absolute paths without treating URLs as paths."""
    stripped = value.strip()
    return stripped.startswith("/") or bool(_WINDOWS_ABSOLUTE_RE.match(stripped))


def portable_path(value: str, project_root: Path | None = None) -> str:
    """Return a non-personal representation of one path-like value."""
    project_root = project_root or REPO_ROOT
    if not value or not (
        is_machine_absolute_path(value) or contains_personal_absolute_path(value)
    ):
        return value
    normalised = value.replace("\\", "/")
    root_normalised = str(project_root).replace("\\", "/").rstrip("/")
    if normalised.casefold().startswith(root_normalised.casefold() + "/"):
        return normalised[len(root_normalised) + 1:]
    path = Path(normalised)
    try:
        return path.resolve().relative_to(project_root.resolve()).as_posix()
    except (OSError, ValueError):
        return f"external_input/{path.name}" if path.name else "external_input"


def portable_command(value: str, project_root: Path | None = None) -> str:
    """Sanitise path-bearing command tokens without rewriting arbitrary text."""
    project_root = project_root or REPO_ROOT
    if not value:
        return value
    try:
        tokens = shlex.split(value)
    except ValueError:
        tokens = value.split()
    portable: list[str] = []
    for index, token in enumerate(tokens):
        option = ""
        candidate = token
        if "=" in token and (
            token.startswith("-") or _ENV_ASSIGNMENT_RE.match(token)
        ):
            option, candidate = token.split("=", 1)
            option += "="
        cleaned = portable_path(candidate, project_root)
        if (
            cleaned.endswith(("/python", "/python3"))
            and (index == 0 or "/.venv/bin/python" in candidate.replace("\\", "/"))
        ):
            cleaned = "python"
        portable.append(option + cleaned)
    return shlex.join(portable)


def _is_path_key(key: str) -> bool:
    return (
        key in _PATH_KEYS
        or key.endswith(("_path", "_file", "_dir"))
        or key in {"files", "source_files"}
    )


def sanitize_public_payload(
    value: Any, key: str = "", project_root: Path | None = None
) -> Any:
    """Copy a JSON-compatible payload while sanitising path-specific fields."""
    project_root = project_root or REPO_ROOT
    if isinstance(value, dict):
        return {
            item_key: sanitize_public_payload(item_value, item_key, project_root)
            for item_key, item_value in value.items()
        }
    if isinstance(value, list):
        return [sanitize_public_payload(item, key, project_root) for item in value]
    if not isinstance(value, str):
        return value
    if key == "command":
        return portable_command(value, project_root)
    if _is_path_key(key) and (
        is_machine_absolute_path(value) or contains_personal_absolute_path(value)
    ):
        return portable_path(value, project_root)
    return value


def write_public_download_projections(run_dir: Path) -> list[Path]:
    """Write sanitised copies of raw records that are intentionally downloadable."""
    run_dir = Path(run_dir)
    public_dir = run_dir / "website_indices" / "public"
    sources = {
        "run_config.json": run_dir / "run_config.json",
        "post_cluster_qc.json": (
            run_dir / "results" / "15_domain_architecture" / "post_cluster_qc.json"
        ),
    }
    written: list[Path] = []
    for output_name, source in sources.items():
        if not source.is_file():
            continue
        payload = json.loads(source.read_text(encoding="utf-8"))
        projected = sanitize_public_payload(payload)
        public_dir.mkdir(parents=True, exist_ok=True)
        output = public_dir / output_name
        output.write_text(
            json.dumps(projected, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        written.append(output)
    return written


def rebuild_existing_public_projections(run_dir: Path) -> list[Path]:
    """Sanitise existing website JSON and refresh downloadable projections only."""
    run_dir = Path(run_dir)
    indices_dir = run_dir / "website_indices"
    written: list[Path] = []
    if indices_dir.is_dir():
        for path in sorted(indices_dir.rglob("*.json")):
            if "public" in path.relative_to(indices_dir).parts:
                continue
            payload = json.loads(path.read_text(encoding="utf-8"))
            projected = sanitize_public_payload(payload)
            if projected != payload:
                path.write_text(
                    json.dumps(projected, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8",
                )
                written.append(path)
    written.extend(write_public_download_projections(run_dir))
    return written
