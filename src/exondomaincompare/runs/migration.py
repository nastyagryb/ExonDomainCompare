from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from exondomaincompare.config import RuntimeConfig
from exondomaincompare.runs.legacy import LegacyRunAdapter, LegacyRunError
from exondomaincompare.runs.outputs import build_entry, verify_manifest, write_manifest
from exondomaincompare.runs.registry import (
    RegistryError,
    register_root,
    register_run_binding,
    resolve_run_record,
    unregister,
)
from exondomaincompare.runs.layout import RunLayout, RunLayoutVersion, validate_run_id

JOURNAL_VERSION = "1.0"
JOURNAL_STATES = {
    "planned", "copying", "copied", "validating", "ready", "registered",
    "quarantined", "rolled_back", "failed",
}


class MigrationError(RuntimeError):
    pass


class MigrationSecurityError(MigrationError):
    pass


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _identity_sha256(value: Mapping[str, Any]) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _portable_legacy_value(value: Any, source_run: Path) -> Any:
    if isinstance(value, Mapping):
        return {
            str(key): _portable_legacy_value(child, source_run)
            for key, child in value.items()
        }
    if isinstance(value, list):
        return [_portable_legacy_value(child, source_run) for child in value]
    if isinstance(value, str) and os.path.isabs(value):
        path = Path(value).expanduser().resolve()
        try:
            return "run:" + path.relative_to(source_run).as_posix()
        except ValueError:
            return f"legacy-external:{path.name}"
    return value


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


class MigrationService:
    def __init__(
        self,
        config: RuntimeConfig,
        *,
        free_bytes: Callable[[Path], int] | None = None,
    ):
        self.config = config
        self.free_bytes = free_bytes or (
            lambda path: shutil.disk_usage(path).free)

    @property
    def journals(self) -> Path:
        return self.config.paths.registry / "migrations"

    def _journal_id(
        self, mode: str, source: Path, run_id: str, destination: Path | None
    ) -> str:
        value = "\0".join((
            JOURNAL_VERSION, mode, str(source.resolve()), run_id,
            str(destination.resolve()) if destination else "",
        ))
        return hashlib.sha256(value.encode()).hexdigest()[:20]

    def journal_path(self, journal_id: str) -> Path:
        if not journal_id or not all(c in "0123456789abcdef" for c in journal_id):
            raise MigrationSecurityError("Unsafe journal id.")
        return self.journals / f"{journal_id}.json"

    def read_journal(self, journal_id: str) -> dict[str, Any]:
        path = self.journal_path(journal_id)
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise MigrationError(f"Journal unavailable: {journal_id}") from exc
        if (
            not isinstance(data, dict)
            or data.get("journal_version") != JOURNAL_VERSION
            or data.get("state") not in JOURNAL_STATES
        ):
            raise MigrationError(f"Invalid migration journal: {journal_id}")
        return data

    def _write_journal(
        self, journal: dict[str, Any], state: str, **updates: Any
    ) -> dict[str, Any]:
        if state not in JOURNAL_STATES:
            raise MigrationError(f"Unknown journal state: {state}")
        journal.update(updates)
        journal["state"] = state
        journal["updated_at"] = _now()
        history = list(journal.get("history") or [])
        if not history or history[-1].get("state") != state:
            history.append({"state": state, "at": journal["updated_at"]})
        journal["history"] = history
        _atomic_json(self.journal_path(str(journal["journal_id"])), journal)
        return journal

    def _source_runs(
        self, root: Path, selected: Sequence[str] | None
    ) -> list[Path]:
        source = root.expanduser().resolve()
        if not source.is_dir():
            raise MigrationError(f"Source root does not exist: {source}")
        if source.is_symlink():
            raise MigrationSecurityError("Source root may not be a symlink.")
        requested = set(selected or ())
        for run_id in requested:
            validate_run_id(run_id)
        candidates = []
        if selected:
            candidates = [source / run_id for run_id in selected]
        else:
            candidates = [
                child for child in sorted(source.iterdir())
                if child.is_dir() and not child.name.startswith((".", "_"))
            ]
        result = []
        for candidate in candidates:
            resolved = candidate.resolve()
            if not _is_relative_to(resolved, source) or resolved.parent != source:
                raise MigrationSecurityError("Run source escapes its selected root.")
            try:
                LegacyRunAdapter(resolved, expected_run_id=resolved.name)
            except LegacyRunError as exc:
                if selected:
                    raise MigrationError(str(exc)) from exc
                continue
            result.append(resolved)
        if requested - {path.name for path in result}:
            missing = sorted(requested - {path.name for path in result})
            raise MigrationError(f"Selected runs are unavailable: {missing}")
        return result

    def plan(
        self,
        *,
        source_root: Path,
        mode: str,
        run_ids: Sequence[str] | None = None,
        destination_root: Path | None = None,
    ) -> list[dict[str, Any]]:
        if mode not in {"register", "copy", "move"}:
            raise MigrationError(f"Unsupported migration mode: {mode}")
        source = source_root.expanduser().resolve()
        rows = []
        for run in self._source_runs(source, run_ids):
            destination = (
                (destination_root or self.config.paths.runs).expanduser().resolve()
                / run.name
            ) if mode in {"copy", "move"} else None
            if destination and (
                _is_relative_to(destination, run)
                or _is_relative_to(run, destination)
                or destination == run
            ):
                raise MigrationSecurityError("Source and destination overlap.")
            rows.append({
                "run_id": run.name,
                "mode": mode,
                "source": str(run),
                "destination": str(destination) if destination else "",
                "source_writes": False,
                "network_contacted": False,
            })
        return rows

    def register(
        self,
        *,
        source_root: Path,
        run_ids: Sequence[str] | None = None,
        read_only: bool = True,
        dry_run: bool = False,
    ) -> list[dict[str, Any]]:
        source = source_root.expanduser().resolve()
        runs = self._source_runs(source, run_ids)
        if dry_run:
            return self.plan(
                source_root=source, mode="register", run_ids=run_ids)
        root_id = register_root(
            self.config, source, kind="legacy_external", read_only=read_only)
        results = []
        for run in runs:
            adapter = LegacyRunAdapter(run, expected_run_id=run.name)
            journal_id = self._journal_id("register", run, run.name, None)
            path = self.journal_path(journal_id)
            if path.is_file():
                existing = self.read_journal(journal_id)
                if existing["state"] == "registered":
                    results.append(existing)
                    continue
            journal = {
                "journal_version": JOURNAL_VERSION,
                "journal_id": journal_id,
                "mode": "register",
                "run_id": run.name,
                "source_root": str(source),
                "source": str(run),
                "destination": "",
                "created_at": _now(),
                "source_writes": False,
            }
            self._write_journal(journal, "planned")
            identity = _identity_sha256(adapter.config())
            register_run_binding(
                self.config, run_id=run.name, path=run, root_id=root_id,
                identity_sha256=identity)
            self._write_journal(
                journal, "registered", root_id=root_id,
                identity_sha256=identity,
                rollback={"operation": "unregister", "source_deleted": False},
            )
            results.append(journal)
        return results

    def _snapshot(self, run: Path) -> list[dict[str, Any]]:
        rows = []
        for path in sorted(run.rglob("*")):
            if path.is_symlink():
                raise MigrationSecurityError(
                    f"Symlink in migration source: {path.relative_to(run)}")
            if not path.is_file():
                continue
            relative = path.relative_to(run).as_posix()
            if ".." in Path(relative).parts:
                raise MigrationSecurityError("Unsafe source path.")
            rows.append({
                "source": relative,
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            })
        return rows

    def _destination_for(
        self, adapter: LegacyRunAdapter, relative: str
    ) -> tuple[str, str]:
        path = Path(relative)
        lower = relative.lower()
        name = path.name.lower()
        if name in {"run_config.json", "status.json", "species_list.txt",
                    "gene_config.yaml"}:
            return f"legacy_preserved/{relative}", "retain_raw_identity"
        if "/packages/" in f"/{lower}" or name.endswith(".zip"):
            return f"legacy_preserved/old_packages/{relative}", "retain_and_copy_package"
        if lower.startswith("website_indices/"):
            return "website/indices/" + relative[len("website_indices/"):], "map_derived"
        if lower.startswith("logs/"):
            return f"legacy_preserved/{relative}", "retain_log"
        if any(token in lower for token in (
            "/_ncbi_datasets_cache/", "/cache/", "/tmp/", "/temp/")):
            return f"legacy_preserved/cache_review/{relative}", "retain_uncertain_cache"
        if (
            "final_pre_interpro_proteins_primary.faa" in lower
            or (
                "14_interproscan" in lower and "/input/" in f"/{lower}"
                and name.endswith((".faa", ".fasta", ".fa"))
            )
        ):
            return f"inputs/cluster/legacy/{relative}", "map_cluster_input"
        species = adapter.species()
        scope = (
            f"species/{species[0]}" if len(species) == 1
            else "comparative" if len(species) > 1
            else "legacy_unresolved"
        )
        if "14_interproscan" in lower and "/output/" in f"/{lower}":
            return f"scientific/{scope}/annotations/interpro/{relative}", "map_raw_annotation"
        if "pytmhmm" in lower and "/output/" in f"/{lower}":
            return f"scientific/{scope}/annotations/pytmhmm/{relative}", "map_raw_annotation"
        if any(token in lower for token in ("/figures/", "/plots/")) or name.endswith(
                (".png", ".pdf", ".svg")):
            return f"scientific/{scope}/figures/{relative}", "map_figure"
        if lower.startswith("results/"):
            return f"scientific/{scope}/legacy_mapped/{relative}", "map_scientific"
        return f"legacy_preserved/{relative}", "retain_uncertain"

    def _dispositions(
        self, adapter: LegacyRunAdapter, snapshot: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        rows = []
        destinations: set[str] = set()
        for source in snapshot:
            destination, disposition = self._destination_for(
                adapter, str(source["source"]))
            if destination in destinations:
                raise MigrationError(f"Migration destination collision: {destination}")
            destinations.add(destination)
            rows.append({
                **source,
                "destination": destination,
                "disposition": disposition,
            })
        return rows

    def _same_snapshot(
        self, before: Sequence[Mapping[str, Any]], after: Sequence[Mapping[str, Any]]
    ) -> bool:
        keys = ("source", "size_bytes", "sha256")
        return [
            tuple(row.get(key) for key in keys) for row in before
        ] == [
            tuple(row.get(key) for key in keys) for row in after
        ]

    def _write_disposition(self, root: Path, rows: Sequence[Mapping[str, Any]]) -> Path:
        path = root / "migration" / "disposition.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        columns = [
            "source", "destination", "disposition", "size_bytes", "sha256"]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
                extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def _copy_one(self, source: Path, destination: Path, expected: Mapping[str, Any]) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if (
                destination.is_file()
                and destination.stat().st_size == int(expected["size_bytes"])
                and _sha256(destination) == expected["sha256"]
            ):
                return
            raise MigrationError(f"Partial-copy mismatch: {destination}")
        temporary = destination.with_name(f".{destination.name}.copying")
        with source.open("rb") as src, temporary.open("wb") as dst:
            for chunk in iter(lambda: src.read(1024 * 1024), b""):
                dst.write(chunk)
        if (
            temporary.stat().st_size != int(expected["size_bytes"])
            or _sha256(temporary) != expected["sha256"]
        ):
            recovery = (
                self.config.paths.quarantine / "checksum-failures"
                / f"{temporary.name}-{hashlib.sha256(str(temporary).encode()).hexdigest()[:8]}"
            )
            recovery.parent.mkdir(parents=True, exist_ok=True)
            temporary.replace(recovery)
            raise MigrationError(f"Checksum mismatch while copying {source}")
        os.replace(temporary, destination)

    def copy_run(
        self,
        *,
        source_run: Path,
        destination_root: Path | None = None,
        dry_run: bool = False,
        stop_after: str | None = None,
    ) -> dict[str, Any]:
        source = source_run.expanduser().resolve()
        adapter = LegacyRunAdapter(source, expected_run_id=source.name)
        run_id = validate_run_id(adapter.run_id)
        destination = (
            destination_root or self.config.paths.runs
        ).expanduser().resolve() / run_id
        if (
            destination == source
            or _is_relative_to(destination, source)
            or _is_relative_to(source, destination)
        ):
            raise MigrationSecurityError("Source and destination overlap.")
        journal_id = self._journal_id("copy", source, run_id, destination)
        if dry_run:
            return {
                "journal_version": JOURNAL_VERSION,
                "journal_id": journal_id,
                "mode": "copy",
                "run_id": run_id,
                "source": str(source),
                "destination": str(destination),
                "state": "planned",
                "dry_run": True,
                "source_writes": False,
            }
        journal_path = self.journal_path(journal_id)
        if journal_path.is_file():
            journal = self.read_journal(journal_id)
            if journal["state"] == "registered":
                return journal
        else:
            journal = {
                "journal_version": JOURNAL_VERSION,
                "journal_id": journal_id,
                "mode": "copy",
                "run_id": run_id,
                "source_root": str(source.parent),
                "source": str(source),
                "destination": str(destination),
                "created_at": _now(),
                "source_writes": False,
            }
            self._write_journal(journal, "planned")
        if destination.exists():
            raise MigrationError(f"Destination already exists: {destination}")

        snapshot = journal.get("source_manifest")
        if not snapshot:
            snapshot = self._snapshot(source)
            dispositions = self._dispositions(adapter, snapshot)
            journal["source_manifest"] = snapshot
            journal["dispositions"] = dispositions
            journal["source_manifest_sha256"] = _identity_sha256(
                {"files": snapshot})
            self._write_journal(journal, "planned")
        dispositions = list(journal["dispositions"])
        required = sum(int(row["size_bytes"]) for row in dispositions)
        destination.parent.mkdir(parents=True, exist_ok=True)
        free = self.free_bytes(destination.parent)
        if free < required + max(10 * 1024 * 1024, required // 20):
            self._write_journal(
                journal, "failed", failure="insufficient_space",
                required_bytes=required, free_bytes=free)
            raise MigrationError("Insufficient free space for verified staging copy.")

        stage_parent = self.config.paths.migration_staging / journal_id
        stage = stage_parent / run_id
        if not stage.exists():
            record = _portable_legacy_value({
                key: value for key, value in adapter.config().items()
                if key != "_legacy_adapter"
            }, source)
            record["migration"] = {
                "mode": "copy", "source_manifest_sha256":
                journal["source_manifest_sha256"],
            }
            status = _portable_legacy_value({
                key: value for key, value in adapter.status().items()
                if key != "_legacy_adapter"
            }, source)
            status["migration_state"] = "copying"
            gene = {
                "gene_symbol": record.get("gene_symbol", ""),
                "analysis_id": record.get("analysis_id", ""),
                "event_id": record.get("event_id", ""),
                "event_type": record.get("event_type", ""),
                "source": "legacy-normalized-without-biological-reinterpretation",
            }
            RunLayout(stage, RunLayoutVersion.CANONICAL_V2).initialize(
                run_record=record, status=status, gene=gene,
                species=adapter.species() or ["unknown_species"],
            )
        self._write_journal(
            journal, "copying", staging=str(stage),
            copied_files=sum(
                1 for row in dispositions
                if (stage / row["destination"]).is_file()))
        if stop_after == "copying":
            return journal
        for row in dispositions:
            self._copy_one(
                source / row["source"], stage / row["destination"], row)
        self._write_disposition(stage, dispositions)
        _atomic_json(stage / "migration" / "source_manifest.json", {
            "journal_version": JOURNAL_VERSION,
            "run_id": run_id,
            "source_manifest_sha256": journal["source_manifest_sha256"],
            "files": snapshot,
        })
        compatibility_links = []
        protected_initial = {
            "run.json", "status.json", "config/gene.json", "config/species.tsv",
            "run_config.json", "species_list.txt", "gene_config.yaml",
        }
        for row in dispositions:
            legacy_relative = str(row["source"])
            if legacy_relative in protected_initial:
                continue
            link = stage / legacy_relative
            target = stage / str(row["destination"])
            if link.exists() or link.is_symlink():
                continue
            link.parent.mkdir(parents=True, exist_ok=True)
            relative_target = os.path.relpath(target, start=link.parent)
            link.symlink_to(relative_target)
            if not link.is_file() or not _is_relative_to(link.resolve(), stage.resolve()):
                raise MigrationSecurityError("Compatibility projection escaped staging.")
            compatibility_links.append({
                "legacy_path": "run:" + legacy_relative,
                "canonical_path": "run:" + str(row["destination"]),
                "owner": "framework.run_migration",
                "type": "internal_symlink_projection",
            })
        LegacyRunAdapter(stage).materialize_legacy_compatibility()
        _atomic_json(stage / "migration" / "compatibility_projections.json", {
            "schema_version": "1.0",
            "run_id": run_id,
            "links": compatibility_links,
        })
        self._write_journal(
            journal, "copied", copied_files=len(dispositions),
            copied_bytes=required)
        if stop_after == "copied":
            return journal

        self._write_journal(journal, "validating")
        current = self._snapshot(source)
        if not self._same_snapshot(snapshot, current):
            self._write_journal(journal, "failed", failure="source_changed")
            raise MigrationError("Source changed during migration.")
        for row in dispositions:
            copied = stage / row["destination"]
            if (
                not copied.is_file()
                or copied.stat().st_size != int(row["size_bytes"])
                or _sha256(copied) != row["sha256"]
            ):
                self._write_journal(
                    journal, "failed", failure="destination_checksum_mismatch")
                raise MigrationError("Destination checksum validation failed.")

        candidate_map = adapter.canonical_candidates()
        by_source = {str(row["source"]): row for row in dispositions}
        entries = []
        alternates: dict[str, list[str]] = {}
        for result_type, candidates in candidate_map.items():
            retained = [
                path for path in candidates
                if path.relative_to(source).as_posix() in by_source
            ]
            if not retained:
                continue
            chosen = retained[0]
            if len(retained) > 1:
                alternates[result_type] = [
                    "run:" + by_source[path.relative_to(source).as_posix()]["destination"]
                    for path in retained[1:]
                ]
            mapped = by_source[chosen.relative_to(source).as_posix()]
            entries.append(build_entry(
                run_root=stage, result_type=result_type,
                path=stage / mapped["destination"],
                producer="legacy-mapped-without-recomputation",
            ))
        write_manifest(
            stage, run_id=run_id, dataset_id=run_id, entries=entries,
            regenerated_derived=[{
                "type": "canonical_manifest",
                "source_checksums_match": True,
                "alternates_retained": alternates,
            }],
        )
        if verify_manifest(stage):
            self._write_journal(
                journal, "failed", failure="canonical_output_validation")
            raise MigrationError("Canonical output manifest validation failed.")
        state_path = stage / "status.json"
        state = json.loads(state_path.read_text(encoding="utf-8"))
        state["migration_state"] = "ready"
        _atomic_json(state_path, state)
        self._write_journal(journal, "ready")
        if stop_after == "ready":
            return journal

        # Packages are user artifacts, not scientific run contents. Verify and
        # stage them before the destination or registry becomes visible.
        package_root = self.config.paths.packages / run_id
        package_stage = stage_parent / "packages"
        package_plan: list[tuple[Path, Path, dict[str, Any]]] = []
        for package in adapter.old_packages():
            target = package_root / package.name
            expected = {
                "size_bytes": package.stat().st_size,
                "sha256": _sha256(package),
            }
            if target.exists():
                if (
                    not target.is_file()
                    or target.stat().st_size != expected["size_bytes"]
                    or _sha256(target) != expected["sha256"]
                ):
                    self._write_journal(
                        journal, "failed", failure="package_destination_collision")
                    raise MigrationError(
                        f"Package destination collision: {target.name}")
            else:
                staged_package = package_stage / package.name
                self._copy_one(package, staged_package, expected)
                package_plan.append((staged_package, target, expected))

        os.replace(stage, destination)
        root_id = register_root(
            self.config, destination.parent, kind="canonical", read_only=False,
            root_id=(
                "configured-runs" if destination.parent == self.config.paths.runs
                else None),
        )
        try:
            register_run_binding(
                self.config, run_id=run_id, path=destination, root_id=root_id,
                identity_sha256=_identity_sha256(adapter.config()))
        except Exception as exc:
            rollback = (
                self.config.paths.quarantine / "registry-failures"
                / f"{run_id}-{journal_id}"
            )
            rollback.parent.mkdir(parents=True, exist_ok=True)
            os.replace(destination, rollback)
            self._write_journal(
                journal, "rolled_back", failure="registry_write_failed",
                recovery_path=str(rollback), error=type(exc).__name__)
            raise MigrationError("Registry failed after copy; destination quarantined.") from exc

        copied_packages = []
        try:
            for staged_package, target, expected in package_plan:
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.exists():
                    raise MigrationError(
                        f"Package destination changed during migration: {target.name}")
                os.replace(staged_package, target)
                copied_packages.append({
                    "path": f"package:{run_id}/{target.name}",
                    "sha256": expected["sha256"],
                })
            for package in adapter.old_packages():
                target = package_root / package.name
                if not any(row["path"].endswith("/" + target.name)
                           for row in copied_packages):
                    copied_packages.append({
                        "path": f"package:{run_id}/{target.name}",
                        "sha256": _sha256(target),
                    })
        except Exception as exc:
            try:
                unregister(self.config, run_id=run_id)
            except RegistryError:
                pass
            rollback = (
                self.config.paths.quarantine / "package-finalization-failures"
                / f"{run_id}-{journal_id}"
            )
            rollback.parent.mkdir(parents=True, exist_ok=True)
            if destination.is_dir():
                os.replace(destination, rollback)
            self._write_journal(
                journal, "rolled_back", failure="package_finalization_failed",
                recovery_path=str(rollback), error=type(exc).__name__)
            raise MigrationError(
                "Package finalization failed; registry and destination rolled back."
            ) from exc
        state = json.loads((destination / "status.json").read_text(encoding="utf-8"))
        state["migration_state"] = "registered"
        _atomic_json(destination / "status.json", state)
        self._write_journal(
            journal, "registered", root_id=root_id,
            destination_manifest_sha256=_sha256(
                destination / "migration" / "source_manifest.json"),
            packages=copied_packages,
            rollback={
                "source_untouched": True,
                "unregister_run_id": run_id,
                "destination_recovery": str(destination),
            },
        )
        return journal

    def copy(
        self,
        *,
        source_root: Path,
        run_ids: Sequence[str] | None = None,
        destination_root: Path | None = None,
        dry_run: bool = False,
    ) -> list[dict[str, Any]]:
        source = source_root.expanduser().resolve()
        if dry_run:
            return self.plan(
                source_root=source, mode="copy", run_ids=run_ids,
                destination_root=destination_root)
        return [
            self.copy_run(
                source_run=run, destination_root=destination_root)
            for run in self._source_runs(source, run_ids)
        ]

    def move(
        self,
        *,
        copy_journal_id: str,
        confirmed: bool,
        dry_run: bool = False,
    ) -> dict[str, Any]:
        if not confirmed:
            raise MigrationError(
                "Move requires explicit confirmation and a successful copy journal.")
        copy_journal = self.read_journal(copy_journal_id)
        if copy_journal.get("mode") != "copy" or copy_journal["state"] != "registered":
            raise MigrationError("Move requires a registered copy journal.")
        source = Path(copy_journal["source"]).resolve()
        destination = Path(copy_journal["destination"]).resolve()
        run_id = validate_run_id(str(copy_journal["run_id"]))
        journal_id = self._journal_id("move", source, run_id, destination)
        if self.journal_path(journal_id).is_file():
            existing = self.read_journal(journal_id)
            if existing["state"] == "quarantined":
                return existing
        journal = {
            "journal_version": JOURNAL_VERSION,
            "journal_id": journal_id,
            "mode": "move",
            "run_id": run_id,
            "source": str(source),
            "destination": str(destination),
            "copy_journal_id": copy_journal_id,
            "created_at": _now(),
        }
        if dry_run:
            return {**journal, "state": "planned", "dry_run": True}
        self._write_journal(journal, "planned")
        if not source.is_dir() or not destination.is_dir():
            raise MigrationError("Move source or verified destination is unavailable.")
        before = copy_journal.get("source_manifest") or []
        if not self._same_snapshot(before, self._snapshot(source)):
            self._write_journal(journal, "failed", failure="source_changed")
            raise MigrationError("Source changed after the verified copy.")
        if verify_manifest(destination):
            self._write_journal(journal, "failed", failure="destination_changed")
            raise MigrationError("Destination canonical outputs no longer verify.")
        selected = resolve_run_record(self.config, run_id)
        if selected is None or selected.path != destination:
            self._write_journal(journal, "failed", failure="registry_not_destination")
            raise MigrationError("Registry is not bound to the verified destination.")
        self._write_journal(journal, "registered")
        quarantine = (
            self.config.paths.quarantine / "moved-sources"
            / f"{run_id}-{journal_id}"
        )
        if quarantine.exists():
            raise MigrationError(f"Quarantine collision: {quarantine}")
        quarantine.parent.mkdir(parents=True, exist_ok=True)
        os.replace(source, quarantine)
        self._write_journal(
            journal, "quarantined", quarantine=str(quarantine),
            recovery={
                "operation": "move quarantine back to source after verifying destination",
                "source": str(source),
                "quarantine": str(quarantine),
                "registry_remains": str(destination),
            },
        )
        return journal

    def rollback(self, journal_id: str) -> dict[str, Any]:
        journal = self.read_journal(journal_id)
        mode = journal.get("mode")
        if mode == "register" and journal["state"] == "registered":
            unregister(
                self.config, run_id=str(journal["run_id"]),
                root_id=str(journal.get("root_id") or "") or None)
            return self._write_journal(
                journal, "rolled_back", source_deleted=False)
        if mode == "move" and journal["state"] == "quarantined":
            source = Path(journal["source"]).resolve()
            quarantine = Path(journal["quarantine"]).resolve()
            if source.exists() or not quarantine.is_dir():
                raise MigrationError("Move rollback path is not safe.")
            source.parent.mkdir(parents=True, exist_ok=True)
            os.replace(quarantine, source)
            return self._write_journal(
                journal, "rolled_back", source_restored=True)
        if mode == "copy" and journal["state"] in {"registered", "ready", "rolled_back"}:
            if journal["state"] == "registered":
                unregister(self.config, run_id=str(journal["run_id"]))
            destination = Path(journal["destination"]).resolve()
            recovery = self.config.paths.quarantine / "copy-rollbacks" / (
                f"{journal['run_id']}-{journal_id}")
            if destination.is_dir():
                recovery.parent.mkdir(parents=True, exist_ok=True)
                os.replace(destination, recovery)
            return self._write_journal(
                journal, "rolled_back", recovery_path=str(recovery),
                source_untouched=True)
        raise MigrationError(
            f"Journal {journal_id} cannot be rolled back from {journal['state']}.")
