"""Checksum-backed canonical scientific output identities."""
from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

from exondomaincompare.runs.layout import RunLayout, RunLayoutVersion, validate_run_id

MANIFEST_VERSION = "1.0"


class CanonicalOutputError(ValueError):
    pass


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_entry(
    *,
    run_root: Path,
    result_type: str,
    path: Path,
    producer: str,
    schema_version: str = "1.0",
) -> dict[str, Any]:
    run_root = run_root.resolve()
    resolved = path.resolve()
    try:
        relative = resolved.relative_to(run_root)
    except ValueError as exc:
        raise CanonicalOutputError("Canonical output is outside its run.") from exc
    if not resolved.is_file() or resolved.is_symlink():
        raise CanonicalOutputError(f"Canonical output is unavailable: {path}")
    return {
        "result_type": str(result_type),
        "path": "run:" + relative.as_posix(),
        "schema_version": str(schema_version),
        "producer": str(producer),
        "size_bytes": resolved.stat().st_size,
        "sha256": file_sha256(resolved),
    }


def write_manifest(
    run_root: Path,
    *,
    run_id: str,
    dataset_id: str | None = None,
    entries: Iterable[Mapping[str, Any]],
    regenerated_derived: Iterable[Mapping[str, Any]] = (),
) -> Path:
    validate_run_id(run_id)
    layout = RunLayout(run_root, RunLayoutVersion.CANONICAL_V2)
    rows = [dict(row) for row in entries]
    types = [str(row.get("result_type") or "") for row in rows]
    duplicates = sorted({value for value in types if types.count(value) > 1})
    if not all(types) or duplicates:
        raise CanonicalOutputError(
            f"Every canonical result type must be unique; duplicates={duplicates}.")
    for row in rows:
        reference = str(row.get("path") or "")
        if not reference.startswith("run:"):
            raise CanonicalOutputError("Canonical output paths must be run: references.")
        path = run_root / reference[4:]
        if (
            not path.is_file()
            or path.is_symlink()
            or file_sha256(path) != row.get("sha256")
            or path.stat().st_size != int(row.get("size_bytes", -1))
        ):
            raise CanonicalOutputError(
                f"Canonical output identity mismatch: {reference}")
    document = {
        "manifest_version": MANIFEST_VERSION,
        "run_id": run_id,
        "dataset_id": dataset_id or run_id,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "outputs": sorted(rows, key=lambda row: str(row["result_type"])),
        "regenerated_derived": list(regenerated_derived),
    }
    layout.ensure_parent_for(layout.outputs_manifest)
    temporary = layout.outputs_manifest.with_name(".canonical_outputs.json.tmp")
    temporary.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, layout.outputs_manifest)
    return layout.outputs_manifest


def read_manifest(run_root: Path) -> dict[str, Any]:
    path = Path(run_root) / "outputs" / "canonical_outputs.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise CanonicalOutputError(f"Canonical output manifest unavailable: {exc}") from exc
    if not isinstance(data, dict) or data.get("manifest_version") != MANIFEST_VERSION:
        raise CanonicalOutputError("Unsupported canonical output manifest.")
    return data


def verify_manifest(run_root: Path) -> list[str]:
    run_root = Path(run_root).resolve()
    data = read_manifest(run_root)
    problems: list[str] = []
    seen: set[str] = set()
    for row in data.get("outputs", []):
        result_type = str(row.get("result_type") or "")
        if not result_type or result_type in seen:
            problems.append(f"duplicate_or_missing_result_type:{result_type}")
        seen.add(result_type)
        reference = str(row.get("path") or "")
        if not reference.startswith("run:") or ".." in Path(reference[4:]).parts:
            problems.append(f"unsafe_path:{reference}")
            continue
        path = (run_root / reference[4:]).resolve()
        try:
            path.relative_to(run_root)
        except ValueError:
            problems.append(f"path_escape:{reference}")
            continue
        if not path.is_file() or path.is_symlink():
            problems.append(f"missing:{reference}")
        elif path.stat().st_size != int(row.get("size_bytes", -1)):
            problems.append(f"size_mismatch:{reference}")
        elif file_sha256(path) != row.get("sha256"):
            problems.append(f"checksum_mismatch:{reference}")
    return problems
