from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

from exondomaincompare.config import RuntimeConfig
from exondomaincompare.contracts import stamp_payload
from exondomaincompare.runs.layout import validate_run_id

REGISTRY_VERSION = "2.0"
REGISTRY_FILE = "local_registry.json"
MANAGED_ROOT_IDS = {
    "configured-runs", "bundled-release-datasets", "repository-legacy-runs",
}


class RegistryError(RuntimeError):
    pass


class RunCollisionError(RegistryError):
    def __init__(self, run_id: str, candidates: Iterable["RunRecord"]):
        self.run_id = run_id
        self.candidates = tuple(candidates)
        roots = ", ".join(row.root_id for row in self.candidates)
        super().__init__(
            f"Run id {run_id!r} exists in multiple roots ({roots}); "
            "select or register an explicit run binding."
        )


@dataclass(frozen=True)
class RunRecord:
    run_id: str
    path: Path
    root_id: str
    kind: str
    read_only: bool
    explicit: bool = False

    @property
    def dataset_id(self) -> str:
        return f"run:{self.run_id}"


def registry_path(config: RuntimeConfig) -> Path:
    return config.paths.registry / REGISTRY_FILE


def _root_id(kind: str, path: Path) -> str:
    digest = hashlib.sha256(str(path.resolve()).encode()).hexdigest()[:12]
    return f"{kind}-{digest}"


def required_roots(config: RuntimeConfig) -> list[dict[str, Any]]:
    roots = [{
        "id": "configured-runs",
        "kind": "canonical",
        "path": str(config.paths.runs),
        "read_only": False,
    }]
    bundled = config.repository_root / "datasets" / "runs"
    if bundled.is_dir():
        roots.append({
            "id": "bundled-release-datasets",
            "kind": "bundled_example",
            "path": str(bundled),
            "read_only": True,
        })
    if config.paths.legacy_runs != config.paths.runs and config.paths.legacy_runs.is_dir():
        roots.append({
            "id": "repository-legacy-runs",
            "kind": "repository_legacy",
            "path": str(config.paths.legacy_runs),
            "read_only": False,
        })
    return roots


def initial_registry(config: RuntimeConfig) -> dict[str, Any]:
    return stamp_payload(
        {
            "registry_version": REGISTRY_VERSION,
            "roots": required_roots(config),
            "runs": [],
            "created_at": datetime.now(timezone.utc).isoformat(),
        },
        payload_type="local_registry",
        profile=config.public_identity(),
        generator="framework.local_registry",
    )


def _normalize(config: RuntimeConfig, value: Mapping[str, Any]) -> dict[str, Any]:
    data = dict(value)
    version = str(data.get("registry_version") or "1.0")
    if version not in {"1.0", REGISTRY_VERSION}:
        raise RegistryError(f"Unsupported local registry version: {version!r}.")
    data["registry_version"] = REGISTRY_VERSION
    roots = [dict(row) for row in data.get("roots", []) if isinstance(row, Mapping)]
    by_id = {
        str(row.get("id")): row for row in roots
        if row.get("id") and str(row.get("id")) not in MANAGED_ROOT_IDS
    }
    for row in required_roots(config):
        by_id[row["id"]] = row
    data["roots"] = list(by_id.values())
    data["runs"] = [
        dict(row) for row in data.get("runs", []) if isinstance(row, Mapping)
    ]
    data["hidden_runs"] = [
        dict(row) for row in data.get("hidden_runs", []) if isinstance(row, Mapping)
    ]
    return data


def read_registry(config: RuntimeConfig) -> dict[str, Any]:
    path = registry_path(config)
    if not path.is_file():
        return initial_registry(config)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RegistryError(f"Could not read private registry: {exc}") from exc
    if not isinstance(data, dict):
        raise RegistryError("Local registry must be a JSON object.")
    return _normalize(config, data)


def write_registry(config: RuntimeConfig, data: Mapping[str, Any]) -> Path:
    path = registry_path(config)
    path.parent.mkdir(parents=True, exist_ok=True)
    normalized = _normalize(config, data)
    normalized["updated_at"] = datetime.now(timezone.utc).isoformat()
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(normalized, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)
    return path


def write_initial_registry(config: RuntimeConfig) -> Path:
    return write_registry(config, read_registry(config))


def _safe_root(row: Mapping[str, Any]) -> Path | None:
    raw = str(row.get("path") or "")
    if not raw:
        return None
    try:
        return Path(raw).expanduser().resolve()
    except OSError:
        return None


def _candidate(root: Path, run_id: str) -> Path | None:
    candidate = (root / run_id).resolve()
    try:
        candidate.relative_to(root)
    except ValueError:
        return None
    return candidate if candidate.is_dir() else None


def discover_candidates(config: RuntimeConfig, run_id: str) -> list[RunRecord]:
    validate_run_id(run_id)
    registry = read_registry(config)
    hidden = {
        (str(row.get("run_id") or ""), str(row.get("root_id") or ""))
        for row in registry.get("hidden_runs", [])
    }
    result: list[RunRecord] = []
    for row in registry.get("roots", []):
        root = _safe_root(row)
        if root is None or not root.is_dir():
            continue
        candidate = _candidate(root, run_id)
        if candidate is None:
            continue
        root_id = str(row.get("id") or _root_id("root", root))
        if (run_id, root_id) in hidden:
            continue
        result.append(RunRecord(
            run_id=run_id,
            path=candidate,
            root_id=root_id,
            kind=str(row.get("kind") or "legacy_external"),
            read_only=bool(row.get("read_only", True)),
        ))
    return result


def _explicit_binding(config: RuntimeConfig, run_id: str) -> RunRecord | None:
    registry = read_registry(config)
    roots = {str(row.get("id")): row for row in registry.get("roots", [])}
    matches = [
        row for row in registry.get("runs", [])
        if str(row.get("run_id") or "") == run_id and bool(row.get("selected", True))
    ]
    if len(matches) > 1:
        raise RegistryError(f"Multiple selected registry bindings for {run_id!r}.")
    if not matches:
        return None
    row = matches[0]
    root_row = roots.get(str(row.get("root_id") or ""), {})
    raw_path = str(row.get("path") or "")
    path = Path(raw_path).expanduser().resolve() if raw_path else None
    if path is None or not path.is_dir() or path.name != run_id:
        raise RegistryError(f"Selected registry binding for {run_id!r} is unavailable.")
    return RunRecord(
        run_id=run_id, path=path,
        root_id=str(row.get("root_id") or ""),
        kind=str(root_row.get("kind") or row.get("kind") or "legacy_external"),
        read_only=bool(root_row.get("read_only", row.get("read_only", True))),
        explicit=True,
    )


def resolve_run_record(config: RuntimeConfig, run_id: str) -> RunRecord | None:
    validate_run_id(run_id)
    explicit = _explicit_binding(config, run_id)
    if explicit is not None:
        return explicit
    candidates = discover_candidates(config, run_id)
    if not candidates:
        return None
    user_candidates = [row for row in candidates if row.kind != "bundled_example"]
    if user_candidates:
        candidates = user_candidates
    if len(candidates) > 1:
        raise RunCollisionError(run_id, candidates)
    return candidates[0]


def discover_run(config: RuntimeConfig, run_id: str) -> Path | None:
    record = resolve_run_record(config, run_id)
    return record.path if record else None


def discover_runs(config: RuntimeConfig) -> tuple[list[RunRecord], dict[str, list[RunRecord]]]:
    registry = read_registry(config)
    names: set[str] = set()
    for row in registry.get("roots", []):
        root = _safe_root(row)
        if root is None or not root.is_dir():
            continue
        for child in root.iterdir():
            if (
                child.is_dir()
                and not child.name.startswith((".", "_"))
            ):
                try:
                    validate_run_id(child.name)
                except ValueError:
                    continue
                names.add(child.name)
    records: list[RunRecord] = []
    collisions: dict[str, list[RunRecord]] = {}
    for run_id in sorted(names):
        try:
            record = resolve_run_record(config, run_id)
        except RunCollisionError as exc:
            collisions[run_id] = list(exc.candidates)
            continue
        if record is not None:
            records.append(record)
    return records, collisions


def register_root(
    config: RuntimeConfig,
    root: Path,
    *,
    kind: str = "legacy_external",
    read_only: bool = True,
    root_id: str | None = None,
) -> str:
    source = root.expanduser().resolve()
    if not source.is_dir():
        raise RegistryError(f"Run root does not exist: {source}")
    registry = read_registry(config)
    for row in registry.get("roots", []):
        if _safe_root(row) == source:
            return str(row["id"])
    identifier = root_id or _root_id(kind, source)
    if any(str(row.get("id")) == identifier for row in registry.get("roots", [])):
        raise RegistryError(f"Registry root id collision: {identifier}")
    registry["roots"].append({
        "id": identifier, "kind": kind, "path": str(source),
        "read_only": bool(read_only),
    })
    write_registry(config, registry)
    return identifier


def register_run_binding(
    config: RuntimeConfig,
    *,
    run_id: str,
    path: Path,
    root_id: str,
    identity_sha256: str,
) -> None:
    validate_run_id(run_id)
    resolved = path.expanduser().resolve()
    if resolved.name != run_id or not resolved.is_dir():
        raise RegistryError("Run binding path/identity mismatch.")
    registry = read_registry(config)
    rows = [
        row for row in registry.get("runs", [])
        if str(row.get("run_id") or "") != run_id
    ]
    rows.append({
        "run_id": run_id,
        "path": str(resolved),
        "root_id": root_id,
        "identity_sha256": identity_sha256,
        "selected": True,
        "registered_at": datetime.now(timezone.utc).isoformat(),
    })
    registry["runs"] = rows
    registry["hidden_runs"] = [
        row for row in registry.get("hidden_runs", [])
        if not (
            str(row.get("run_id") or "") == run_id
            and str(row.get("root_id") or "") == root_id
        )
    ]
    write_registry(config, registry)


def unregister(
    config: RuntimeConfig, *, root_id: str | None = None, run_id: str | None = None
) -> dict[str, int]:
    registry = read_registry(config)
    removed_roots = removed_runs = 0
    if run_id:
        validate_run_id(run_id)
        matching = [
            row for row in registry.get("runs", [])
            if str(row.get("run_id") or "") == run_id
        ]
        before = len(registry.get("runs", []))
        registry["runs"] = [
            row for row in registry.get("runs", [])
            if str(row.get("run_id") or "") != run_id
        ]
        removed_runs = before - len(registry["runs"])
        hidden = list(registry.get("hidden_runs", []))
        for row in matching:
            marker = {
                "run_id": run_id,
                "root_id": str(row.get("root_id") or ""),
                "hidden_at": datetime.now(timezone.utc).isoformat(),
            }
            if not any(
                old.get("run_id") == marker["run_id"]
                and old.get("root_id") == marker["root_id"]
                for old in hidden
            ):
                hidden.append(marker)
        registry["hidden_runs"] = hidden
    if root_id:
        if root_id in {"configured-runs", "repository-legacy-runs"}:
            raise RegistryError("Required roots cannot be unregistered.")
        before = len(registry.get("roots", []))
        registry["roots"] = [
            row for row in registry.get("roots", [])
            if str(row.get("id") or "") != root_id
        ]
        removed_roots = before - len(registry["roots"])
        before_runs = len(registry.get("runs", []))
        registry["runs"] = [
            row for row in registry.get("runs", [])
            if str(row.get("root_id") or "") != root_id
        ]
        removed_runs += before_runs - len(registry["runs"])
        registry["hidden_runs"] = [
            row for row in registry.get("hidden_runs", [])
            if str(row.get("root_id") or "") != root_id
        ]
    write_registry(config, registry)
    return {"roots": removed_roots, "runs": removed_runs}


def hide_discovered_run(
    config: RuntimeConfig, *, run_id: str, root_id: str
) -> dict[str, int]:
    validate_run_id(run_id)
    if not root_id:
        raise RegistryError("A root id is required to hide a discovered run.")
    registry = read_registry(config)
    known_root_ids = {
        str(row.get("id") or "") for row in registry.get("roots", [])
    }
    if root_id not in known_root_ids:
        raise RegistryError(f"Unknown registry root: {root_id}")
    hidden = list(registry.get("hidden_runs", []))
    already_hidden = any(
        str(row.get("run_id") or "") == run_id
        and str(row.get("root_id") or "") == root_id
        for row in hidden
    )
    if not already_hidden:
        hidden.append({
            "run_id": run_id,
            "root_id": root_id,
            "hidden_at": datetime.now(timezone.utc).isoformat(),
        })
    registry["hidden_runs"] = hidden
    before = len(registry.get("runs", []))
    registry["runs"] = [
        row for row in registry.get("runs", [])
        if not (
            str(row.get("run_id") or "") == run_id
            and str(row.get("root_id") or "") == root_id
        )
    ]
    write_registry(config, registry)
    return {
        "hidden": 0 if already_hidden else 1,
        "bindings_removed": before - len(registry["runs"]),
    }
