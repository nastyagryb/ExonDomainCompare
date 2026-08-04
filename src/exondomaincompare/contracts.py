from __future__ import annotations

import copy
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

CONTRACT_VERSION = "1.0"
RESERVED_KEY = "_exondomain"
SELF_HASH_EXCLUDED_FIELDS = frozenset({"content_sha256", "generated_at"})


class ContractIdentityError(ValueError):
    pass


def file_sha256(path: Path | str, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _without_self_hash(value: Any) -> Any:
    cleaned = copy.deepcopy(value)
    if isinstance(cleaned, dict) and isinstance(cleaned.get(RESERVED_KEY), dict):
        for field in SELF_HASH_EXCLUDED_FIELDS:
            cleaned[RESERVED_KEY].pop(field, None)
    return cleaned


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        _without_self_hash(value), ensure_ascii=False, sort_keys=True,
        separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")


def canonical_json_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()




def tree_content_identity(root: Path | str, *,
                          suffixes: Iterable[str] | None = None,
                          ignored_parts: Iterable[str] = ()) -> dict[str, Any]:
    base = Path(root)
    allowed = {suffix.lower() for suffix in suffixes} if suffixes else None
    ignored = set(ignored_parts)
    files = []
    if base.is_file():
        candidates = [base]
        label_root = base.parent
    elif base.is_dir():
        candidates = sorted(base.rglob("*"))
        label_root = base
    else:
        candidates = []
        label_root = base
    for path in candidates:
        if not path.is_file() or any(part in ignored for part in path.parts):
            continue
        if allowed is not None and path.suffix.lower() not in allowed:
            continue
        files.append({
            "path": str(path.relative_to(label_root)),
            "sha256": file_sha256(path),
            "bytes": path.stat().st_size,
        })
    return {
        "algorithm": "sha256",
        "file_count": len(files),
        "files": files,
        "content_sha256": canonical_json_sha256(files),
    }


_SCIENTIFIC_SOURCE_SUFFIXES = {
    ".json", ".tsv", ".csv", ".txt", ".faa", ".fasta", ".fa", ".gff3", ".yaml", ".yml",
}
_SCIENTIFIC_SOURCE_IGNORES = {
    "_ncbi_datasets_cache", "00_run_setup", "website_indices", "figures", "plots",
    "packages",
}


def run_source_identity(run_dir: Path | str) -> dict[str, Any]:
    run = Path(run_dir)
    if "final_30_until_interpro_prepare" in run.parts:
        return {
            "algorithm": "trusted-release-manifest",
            "manifest": "artifacts/release_phase_a/fgfr2_freeze_manifest.tsv",
            "manifest_sha256": (
                "f757fedbdb28f7d4625735eb3066698f31575bb9a6d7f5adde73abcd301e88a4"
            ),
        }
    source = run / "results" if (run / "results").is_dir() else run
    return tree_content_identity(
        source,
        suffixes=_SCIENTIFIC_SOURCE_SUFFIXES,
        ignored_parts=_SCIENTIFIC_SOURCE_IGNORES,
    )


def write_freshness_contract(run_dir: Path | str, output_dir: Path | str, *,
                             generator: str) -> Path:
    run = Path(run_dir)
    out = Path(output_dir) / "_freshness.json"
    payload = stamp_payload(
        {
            "run_id": run.name,
            "dataset_id": run.name,
            "source_identity": run_source_identity(run),
        },
        payload_type="freshness",
        run_id=run.name,
        dataset_id=run.name,
        generator=generator,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return out


def write_payload_contracts(output_dir: Path | str, *, run_id: str,
                            dataset_id: str, generator: str) -> Path:
    directory = Path(output_dir)
    payloads = {}
    for path in sorted(directory.glob("*.json")):
        if path.name.startswith("_") or path.name == "payload_contracts.json":
            continue
        payloads[path.name] = {
            "sha256": file_sha256(path),
            "bytes": path.stat().st_size,
            "run_id": run_id,
            "dataset_id": dataset_id,
        }
    contract = stamp_payload(
        {
            "run_id": run_id,
            "dataset_id": dataset_id,
            "payloads": payloads,
        },
        payload_type="payload_contracts",
        run_id=run_id,
        dataset_id=dataset_id,
        generator=generator,
    )
    path = directory / "_payload_contracts.json"
    path.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    return path


def verify_payload_contract(path: Path | str, contract: Mapping[str, Any], *,
                            expected_run_id: str, expected_dataset_id: str) -> tuple[bool, str]:
    file_path = Path(path)
    entry = (contract.get("payloads") or {}).get(file_path.name)
    if not isinstance(entry, Mapping):
        return False, "payload is absent from the contract sidecar"
    if entry.get("run_id") != expected_run_id:
        return False, "payload run identity mismatch"
    if entry.get("dataset_id") != expected_dataset_id:
        return False, "payload dataset identity mismatch"
    if not file_path.is_file() or file_sha256(file_path) != entry.get("sha256"):
        return False, "payload checksum mismatch"
    return True, "payload identity and checksum match"


def freshness_verdict(run_dir: Path | str, record: Mapping[str, Any]) -> tuple[bool, str]:
    expected_run = str(record.get("run_id") or record.get(RESERVED_KEY, {}).get("run_id") or "")
    run = Path(run_dir)
    if expected_run and expected_run != run.name:
        return False, "run identity mismatch"
    current = run_source_identity(run)
    if current.get("content_sha256") != record.get("source_identity", {}).get("content_sha256"):
        if current != record.get("source_identity"):
            return False, "source content checksum changed"
    return True, "source content checksum matches"


def stamp_payload(payload: Mapping[str, Any], *, payload_type: str,
                  run_id: str = "", dataset_id: str = "",
                  profile: Mapping[str, str] | None = None,
                  sources: Iterable[Mapping[str, str]] = (),
                  generator: str = "") -> dict[str, Any]:
    result = copy.deepcopy(dict(payload))
    metadata = dict(result.get(RESERVED_KEY) or {})
    metadata.update({
        "contract_version": CONTRACT_VERSION,
        "payload_type": payload_type,
        "run_id": str(run_id or result.get("run_id") or ""),
        "dataset_id": str(dataset_id or result.get("dataset_id") or ""),
        "generator": generator,
        "profile": dict(profile or {}),
        "source_checksums": [dict(row) for row in sources],
        "generated_at": datetime.now(timezone.utc).isoformat(),
    })
    result[RESERVED_KEY] = metadata
    metadata["content_sha256"] = canonical_json_sha256(result)
    return result


def normalize_payload(payload: Any, *, payload_type: str,
                      expected_run_id: str = "", expected_dataset_id: str = "") -> Any:
    if not isinstance(payload, dict):
        return payload
    metadata = payload.get(RESERVED_KEY)
    if not isinstance(metadata, dict):
        result = copy.deepcopy(payload)
        result[RESERVED_KEY] = {
            "contract_version": "legacy-unversioned",
            "payload_type": payload_type,
            "run_id": str(payload.get("run_id") or expected_run_id or ""),
            "dataset_id": str(payload.get("dataset_id") or expected_dataset_id or ""),
            "compatibility_adapter": "phase_b_legacy_read",
        }
        return result
    actual_run = str(metadata.get("run_id") or payload.get("run_id") or "")
    actual_dataset = str(metadata.get("dataset_id") or payload.get("dataset_id") or "")
    if expected_run_id and actual_run and actual_run != expected_run_id:
        raise ContractIdentityError(
            f"Payload run identity mismatch: expected {expected_run_id!r}, got {actual_run!r}."
        )
    if expected_dataset_id and actual_dataset and actual_dataset != expected_dataset_id:
        raise ContractIdentityError(
            "Payload dataset identity mismatch: "
            f"expected {expected_dataset_id!r}, got {actual_dataset!r}."
        )
    expected_hash = metadata.get("content_sha256")
    if expected_hash and expected_hash != canonical_json_sha256(payload):
        raise ContractIdentityError("Payload content checksum does not match its contract metadata.")
    return copy.deepcopy(payload)


def portable_path_reference(path: Path | str, *, repository_root: Path,
                            run_root: Path | None = None) -> str:
    resolved = Path(path).expanduser().resolve()
    if run_root:
        try:
            return "run:" + str(resolved.relative_to(run_root.resolve()))
        except ValueError:
            pass
    try:
        return "repo:" + str(resolved.relative_to(repository_root.resolve()))
    except ValueError:
        return "external:" + resolved.name


def portable_runtime_record(value: Any, *, repository_root: Path,
                            run_root: Path | None = None) -> Any:
    if isinstance(value, dict):
        return {
            key: portable_runtime_record(
                child, repository_root=repository_root, run_root=run_root
            )
            for key, child in value.items()
        }
    if isinstance(value, list):
        return [
            portable_runtime_record(
                child, repository_root=repository_root, run_root=run_root
            )
            for child in value
        ]
    if not isinstance(value, str):
        return value

    text = value
    replacements: list[tuple[str, str]] = []
    if run_root:
        replacements.append((str(run_root.expanduser().resolve()), "run:."))
    replacements.append((str(repository_root.expanduser().resolve()), "repo:."))
    for absolute, logical in replacements:
        text = text.replace(absolute, logical)

    if text.startswith(("run:", "repo:", "external:", "PATH:")):
        return text
    if Path(text).expanduser().is_absolute():
        return portable_path_reference(
            text, repository_root=repository_root, run_root=run_root
        )

    text = re.sub(
        r"(?<!\S)/(?:[^/\s]+/)+(?P<tool>python(?:3(?:\.\d+)?)?|node|npm|"
        r"datasets|mafft|ssh|scp|rsync)(?=\s|$)",
        lambda match: f"PATH:{match.group('tool')}",
        text,
    )
    return text


def resolve_path_reference(reference: str, *, repository_root: Path,
                           run_root: Path | None = None) -> Path:
    if reference.startswith("run:") and run_root:
        return (run_root / reference[4:]).resolve()
    if reference.startswith("repo:"):
        return (repository_root / reference[5:]).resolve()
    if reference.startswith("external:"):
        raise ValueError("External portable references require an explicit execution-time path.")
    legacy = Path(reference).expanduser()
    if legacy.is_absolute():
        return legacy.resolve()
    return (repository_root / legacy).resolve()
