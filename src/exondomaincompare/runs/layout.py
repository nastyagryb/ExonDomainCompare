from __future__ import annotations

import csv
import json
import os
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Sequence


CANONICAL_LAYOUT_VERSION = "canonical-2.0"
STATUS_SCHEMA_VERSION = "2.0"
RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")


class RunLayoutError(ValueError):
    pass


class RunLayoutVersion(str, Enum):
    LEGACY_V1 = "legacy-1"
    CANONICAL_V2 = CANONICAL_LAYOUT_VERSION


def validate_run_id(run_id: str) -> str:
    value = str(run_id or "")
    if not RUN_ID_PATTERN.fullmatch(value) or value in {".", ".."}:
        raise RunLayoutError(f"Unsafe run id: {value!r}.")
    return value


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


@dataclass(frozen=True)
class RunLayout:
    root: Path
    version: RunLayoutVersion = RunLayoutVersion.LEGACY_V1

    @property
    def run(self) -> Path:
        return self.root / (
            "run.json" if self.version is RunLayoutVersion.CANONICAL_V2
            else "run_config.json"
        )

    @property
    def config(self) -> Path:
        return self.run

    @property
    def status(self) -> Path:
        return self.root / "status.json"

    @property
    def gene_config(self) -> Path:
        return (
            self.root / "config" / "gene.json"
            if self.version is RunLayoutVersion.CANONICAL_V2
            else self.root / "gene_config.yaml"
        )

    @property
    def species_config(self) -> Path:
        return (
            self.root / "config" / "species.tsv"
            if self.version is RunLayoutVersion.CANONICAL_V2
            else self.root / "species_list.txt"
        )

    @property
    def inputs(self) -> Path:
        return self.root / "inputs"

    @property
    def source_inputs(self) -> Path:
        return self.inputs / "sources"

    @property
    def cluster_input(self) -> Path:
        if self.version is RunLayoutVersion.LEGACY_V1:
            return self.root / "results" / "14_interproscan" / "primary" / "input"
        return self.inputs / "cluster"

    @property
    def cluster_primary_fasta(self) -> Path:
        if self.version is RunLayoutVersion.LEGACY_V1:
            return self.cluster_input / "input.fasta"
        return self.cluster_input / "primary.faa"

    @property
    def cluster_output(self) -> Path:
        if self.version is RunLayoutVersion.LEGACY_V1:
            return self.root / "results" / "14_interproscan" / "primary" / "output"
        return self.root / "scientific" / "annotations" / "interpro"

    @property
    def scientific(self) -> Path:
        return self.root / "scientific"

    def species(self, species_id: str) -> Path:
        validate_run_id(species_id)
        return self.scientific / "species" / species_id

    @property
    def comparative(self) -> Path:
        return self.scientific / "comparative"

    @property
    def website_indices(self) -> Path:
        return (
            self.root / "website" / "indices"
            if self.version is RunLayoutVersion.CANONICAL_V2
            else self.root / "website_indices"
        )

    @property
    def website_dataset(self) -> Path:
        return (
            self.root / "website" / "dataset"
            if self.version is RunLayoutVersion.CANONICAL_V2
            else self.root / "website_indices"
        )

    @property
    def outputs_manifest(self) -> Path:
        return self.root / "outputs" / "canonical_outputs.json"

    @property
    def migration(self) -> Path:
        return self.root / "migration"

    @property
    def legacy_preserved(self) -> Path:
        return self.root / "legacy_preserved"

    @property
    def logs(self) -> Path:
        return self.root / "logs"

    @property
    def packages(self) -> Path:
        if self.version is RunLayoutVersion.LEGACY_V1:
            return self.root / "results" / "generic_gene_analysis" / "packages"
        raise RunLayoutError("Canonical packages live under AppPaths.packages.")

    def logical(self, path: Path) -> str:
        resolved_root = self.root.resolve()
        resolved = path.resolve()
        try:
            relative = resolved.relative_to(resolved_root)
        except ValueError as exc:
            raise RunLayoutError("Path is outside the selected run.") from exc
        return "run:" + relative.as_posix()

    def ensure_parent_for(self, path: Path) -> Path:
        resolved_root = self.root.resolve()
        candidate = path.resolve(strict=False)
        try:
            candidate.relative_to(resolved_root)
        except ValueError as exc:
            raise RunLayoutError("Writer path escapes the selected run.") from exc
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def initialize(
        self,
        *,
        run_record: Mapping[str, Any],
        status: Mapping[str, Any],
        gene: Mapping[str, Any],
        species: Sequence[str],
    ) -> tuple[Path, ...]:
        if self.version is not RunLayoutVersion.CANONICAL_V2:
            raise RunLayoutError("Only canonical layouts can be initialized.")
        run_id = validate_run_id(str(run_record.get("run_id") or self.root.name))
        if self.root.name != run_id:
            raise RunLayoutError("Run identity does not match its directory name.")
        if self.root.exists() and any(self.root.iterdir()):
            raise FileExistsError(f"Run destination is not empty: {self.root}")
        normalized_species = []
        for value in species:
            token = str(value or "").strip()
            if not token or "\t" in token or "\n" in token:
                raise RunLayoutError(f"Invalid species identifier: {value!r}.")
            if token not in normalized_species:
                normalized_species.append(token)
        if not normalized_species:
            raise RunLayoutError("A run requires at least one species.")

        self.root.mkdir(parents=True, exist_ok=True)
        record = dict(run_record)
        record.update({
            "schema_version": "2.0",
            "layout_version": CANONICAL_LAYOUT_VERSION,
            "run_id": run_id,
            "dataset_id": str(record.get("dataset_id") or run_id),
            "paths": {
                "gene": "run:config/gene.json",
                "species": "run:config/species.tsv",
                "status": "run:status.json",
            },
        })
        state = dict(status)
        state.update({
            "schema_version": STATUS_SCHEMA_VERSION,
            "layout_version": CANONICAL_LAYOUT_VERSION,
            "run_id": run_id,
            "dataset_id": str(state.get("dataset_id") or record["dataset_id"]),
        })
        gene_record = dict(gene)
        gene_record.setdefault("schema_version", "1.0")
        gene_record.setdefault("run_id", run_id)

        _atomic_json(self.run, record)
        _atomic_json(self.status, state)
        _atomic_json(self.gene_config, gene_record)
        self.species_config.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.species_config.with_name(".species.tsv.tmp")
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["species_id"])
            writer.writerows((item,) for item in normalized_species)
        os.replace(temporary, self.species_config)
        return (self.run, self.status, self.gene_config, self.species_config)

    def initial_tree(self) -> tuple[str, ...]:
        return tuple(
            sorted(
                path.relative_to(self.root).as_posix()
                for path in self.root.rglob("*")
                if path.is_file()
            )
        )


RunPaths = RunLayout


def read_species_tsv(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if rows.fieldnames != ["species_id"]:
            raise RunLayoutError("config/species.tsv must have one species_id column.")
        values = [str(row.get("species_id") or "").strip() for row in rows]
    return [value for value in values if value]
