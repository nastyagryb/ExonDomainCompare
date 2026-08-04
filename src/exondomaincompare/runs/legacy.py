from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from exondomaincompare.runs.layout import (
    CANONICAL_LAYOUT_VERSION,
    RunLayout,
    RunLayoutError,
    RunLayoutVersion,
    read_species_tsv,
    validate_run_id,
)


class LegacyRunError(ValueError):
    pass


class UnsupportedLegacyRun(LegacyRunError):
    pass


def _json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@dataclass(frozen=True)
class LegacyRunDescription:
    run_id: str
    layout_version: str
    supported: bool
    partial: bool
    numbered_fgfr2: bool
    generic_core_mirrors: bool
    old_packages: bool
    has_status: bool
    raw_provenance: tuple[str, ...]
    missing: tuple[str, ...]


class LegacyRunAdapter:
    def __init__(self, root: Path, *, expected_run_id: str | None = None):
        self.root = Path(root).expanduser().resolve()
        if not self.root.is_dir():
            raise UnsupportedLegacyRun(f"Run directory does not exist: {self.root}")
        self._canonical = (self.root / "run.json").is_file()
        self.layout = RunLayout(
            self.root,
            RunLayoutVersion.CANONICAL_V2 if self._canonical
            else RunLayoutVersion.LEGACY_V1,
        )
        identity = self.config()
        run_id = str(identity.get("run_id") or self.root.name)
        try:
            validate_run_id(run_id)
        except RunLayoutError as exc:
            raise LegacyRunError(str(exc)) from exc
        if run_id != self.root.name:
            raise LegacyRunError(
                f"Run identity {run_id!r} does not match directory {self.root.name!r}.")
        if expected_run_id is not None and run_id != expected_run_id:
            raise LegacyRunError(
                f"Expected run {expected_run_id!r}, found {run_id!r}.")
        self.run_id = run_id

    @property
    def is_canonical(self) -> bool:
        return self._canonical

    def config(self) -> dict[str, Any]:
        path = self.root / ("run.json" if self._canonical else "run_config.json")
        data = _json(path)
        if not data:
            if self._canonical:
                raise UnsupportedLegacyRun("Canonical run.json is missing or invalid.")
            data = {"run_id": self.root.name}
        normalized = dict(data)
        normalized.setdefault("run_id", self.root.name)
        normalized.setdefault("dataset_id", normalized["run_id"])
        normalized.setdefault(
            "layout_version",
            CANONICAL_LAYOUT_VERSION if self._canonical else "legacy-1",
        )
        normalized["_legacy_adapter"] = {
            "normalized_in_memory": True,
            "source": f"run:{path.name}",
            "raw_sha256": (
                _sha256_bytes(path.read_bytes()) if path.is_file() else ""
            ),
        }
        return normalized

    def status(self) -> dict[str, Any]:
        path = self.root / "status.json"
        raw = _json(path)
        state = dict(raw)
        state.setdefault("run_id", self.run_id)
        state.setdefault("dataset_id", self.run_id)
        state.setdefault("schema_version", "legacy-unversioned" if not raw else "1.0")
        if not raw:
            state.update({
                "status": "unknown",
                "current_step": "unknown",
                "next_action": "inspect_missing_status",
                "availability_reason": "status.json is missing or invalid",
            })
        state["_legacy_adapter"] = {
            "normalized_in_memory": True,
            "source": "run:status.json" if path.is_file() else "missing",
            "raw_sha256": _sha256_bytes(path.read_bytes()) if path.is_file() else "",
        }
        return state

    def species(self) -> list[str]:
        if self._canonical:
            path = self.root / "config" / "species.tsv"
            try:
                return read_species_tsv(path)
            except (OSError, RunLayoutError):
                return []
        cfg = self.config()
        values = cfg.get("species_ids")
        if isinstance(values, list):
            return [str(item) for item in values if str(item).strip()]
        path = self.root / "species_list.txt"
        if path.is_file():
            return [
                line.strip() for line in path.read_text(
                    encoding="utf-8", errors="replace").splitlines()
                if line.strip() and not line.lstrip().startswith("#")
            ]
        registry = self.root / "results" / "01_species_registry" / "species_registry.tsv"
        if registry.is_file():
            with registry.open(encoding="utf-8", errors="replace", newline="") as handle:
                rows = csv.DictReader(handle, delimiter="\t")
                return [
                    str(row.get("species_id") or row.get("species") or "").strip()
                    for row in rows
                    if str(row.get("species_id") or row.get("species") or "").strip()
                ]
        return []

    def describe(self) -> LegacyRunDescription:
        results = self.root / "results"
        numbered = any(
            child.is_dir() and child.name[:2].isdigit()
            for child in results.iterdir()
        ) if results.is_dir() else False
        mirrors = any(
            (results / name).is_dir()
            for name in ("core_gene_analysis", "generic_gene_analysis")
        )
        packages = (
            results / "generic_gene_analysis" / "packages"
        ).is_dir()
        missing = []
        if not ((self.root / "run.json").is_file()
                or (self.root / "run_config.json").is_file()):
            missing.append("run_identity")
        if not (self.root / "status.json").is_file():
            missing.append("status")
        if not self.species():
            missing.append("species")
        provenance = [
            item for item, present in (
                ("numbered_fgfr2_stages", numbered),
                ("generic_core_mirrors", mirrors),
                ("old_packages", packages),
                ("partial_or_missing_metadata", bool(missing)),
            ) if present
        ]
        return LegacyRunDescription(
            run_id=self.run_id,
            layout_version=(
                CANONICAL_LAYOUT_VERSION if self._canonical else "legacy-1"),
            supported=not self._canonical or (
                self.config().get("layout_version") == CANONICAL_LAYOUT_VERSION),
            partial=bool(missing),
            numbered_fgfr2=numbered,
            generic_core_mirrors=mirrors,
            old_packages=packages,
            has_status=(self.root / "status.json").is_file(),
            raw_provenance=tuple(provenance),
            missing=tuple(missing),
        )

    def logical(self, path: Path) -> str:
        return self.layout.logical(path)

    def resolve_logical(self, reference: str) -> Path:
        if not str(reference).startswith("run:"):
            raise LegacyRunError("Only run: references are accepted.")
        relative = str(reference)[4:]
        if relative.startswith("/") or ".." in Path(relative).parts:
            raise LegacyRunError("Unsafe logical run path.")
        candidate = (self.root / relative).resolve()
        try:
            candidate.relative_to(self.root)
        except ValueError as exc:
            raise LegacyRunError("Logical path escapes the run.") from exc
        return candidate

    def canonical_candidates(self) -> dict[str, list[Path]]:
        names = {
            "cluster_primary_fasta": {
                "final_pre_interpro_proteins_primary.faa", "proteins_primary.faa"},
            "gene_models": {"gene_model_index.tsv", "genes.tsv"},
            "primary_selection": {
                "primary_selection_evidence.tsv", "selected_transcripts.tsv"},
            "protein_sequences": {
                "proteins_all_isoforms.faa", "selected_fgfr2_proteins.faa"},
            "exon_structure": {"exon_protein_map.tsv"},
            "isoform_differences": {
                "event_region_evidence.tsv", "fgfr2_III_pair_difference_positions.tsv"},
            "comparative_msa": {
                "final_fgfr2_full_length_protein_msa.aln.faa",
                "comparative_protein_msa.aln.faa"},
            "synteny_neighbors": {
                "synteny_neighbors.tsv", "fgfr2_local_gene_neighborhood.tsv"},
            "interpro_raw": {"input.fasta.tsv"},
            "pytmhmm_raw": {"pytmhmm_summary_all.tsv"},
            "domain_features": {
                "domain_features.tsv", "interpro_domain_features_normalized.tsv"},
            "boundary_distances": {"exon_domain_boundary_distances.tsv"},
        }
        result = {key: [] for key in names}
        for path in self.root.rglob("*"):
            if not path.is_file() or path.is_symlink():
                continue
            for key, candidates in names.items():
                if path.name in candidates:
                    result[key].append(path)
        return {key: sorted(values) for key, values in result.items() if values}

    def materialize_legacy_compatibility(self) -> tuple[Path, ...]:
        if not self._canonical:
            return ()
        created: list[Path] = []
        cfg = self.config()
        species = self.species()
        gene = _json(self.root / "config" / "gene.json")
        projections: list[tuple[Path, bytes]] = []
        legacy_config = {
            key: value for key, value in cfg.items()
            if key not in {"_legacy_adapter", "paths"}
        }
        legacy_config.update({
            "run_id": self.run_id,
            "species_ids": species,
            "species_count": len(species),
            "species_list_path": "run:species_list.txt",
            "gene_config": "run:gene_config.yaml",
            "run_dir": "run:.",
            "results_dir": "run:results",
            "website_indices_dir": "run:website_indices",
            "_compatibility_projection": {
                "schema_version": "1.0",
                "source": "run:run.json",
                "source_sha256": _sha256_bytes((self.root / "run.json").read_bytes()),
                "owner": "framework.legacy_run_adapter",
            },
        })
        projections.append((
            self.root / "run_config.json",
            (json.dumps(legacy_config, indent=2, sort_keys=True) + "\n").encode(),
        ))
        projections.append((
            self.root / "species_list.txt",
            ("\n".join(species) + "\n").encode(),
        ))
        gene_projection = dict(gene)
        gene_projection["_compatibility_projection"] = {
            "source": "run:config/gene.json",
            "owner": "framework.legacy_run_adapter",
        }
        projections.append((
            self.root / "gene_config.yaml",
            (json.dumps(gene_projection, indent=2, sort_keys=True) + "\n").encode(),
        ))
        for path, content in projections:
            if path.exists():
                continue
            temporary = path.with_name(f".{path.name}.tmp")
            temporary.write_bytes(content)
            os.replace(temporary, path)
            created.append(path)
        return tuple(created)

    def old_packages(self) -> list[Path]:
        base = self.root / "results" / "generic_gene_analysis" / "packages"
        return sorted(path for path in base.glob("*") if path.is_file()) if base.is_dir() else []
