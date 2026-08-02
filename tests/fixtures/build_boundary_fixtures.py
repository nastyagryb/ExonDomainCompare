#!/usr/bin/env python3
"""Rebuild the compact Boundary-dashboard fixtures from a real run directory.

The dashboard tests used to read two live production runs, so deleting a run broke
the suite and editing one silently changed what was asserted. The fixtures below are
verbatim copies of the *scientific* Core tables of those runs — no value is rewritten
— reduced to the files ``protein_coordinate_model.build_models_for_run`` actually
reads, with the run config and status trimmed to their read keys and all absolute
personal paths dropped.

Usage (only needed when a fixture must be regenerated from a run that still exists):

    python tests/fixtures/build_boundary_fixtures.py <run_dir> <fixture_dir>
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

# Exactly the Core tables the coordinate model + boundary dashboard read.
CORE_TABLES = (
    "exon_protein_map.tsv",
    "domain_features.tsv",
    "interpro_annotations.tsv",
    "tm_features.tsv",
    "exon_domain_boundary_distances.tsv",
    "event_candidate_regions.tsv",
    "primary_selection_evidence.tsv",
    "protein_isoform_index.tsv",
    "proteins_primary.faa",
)
CONFIG_KEYS = ("run_id", "gene_symbol", "species_count", "species_ids", "species_taxids",
               "species_scientific_names", "workflow", "run_mode", "experimental")
STATUS_KEYS = ("run_id", "status", "gene_symbol", "species_count", "run_mode",
               "pre_interpro_status", "cluster_analysis_status", "post_interpro_status")
# `_transcript_models` reads only these; the rest of the served index is display state.
MODEL_KEYS = ("protein_id", "transcript_id", "protein_length", "curation_status",
              "is_primary", "role")
BLOCK_KEYS = ("id", "label", "start", "end", "exon_number", "transcript_exon_number",
              "transcript_id", "shared_exon_group_id", "genomic_start", "genomic_end",
              "cds_start", "cds_end", "phase", "strand", "source", "feature_type")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(run_dir: Path, dst: Path) -> None:
    core_src, core_dst = run_dir / "results/core_gene_analysis", dst / "results/core_gene_analysis"
    core_dst.mkdir(parents=True, exist_ok=True)
    (dst / "results/01_species_registry").mkdir(parents=True, exist_ok=True)

    copied = {}
    for name in CORE_TABLES:
        src = core_src / name
        if not src.is_file():
            continue
        shutil.copy2(src, core_dst / name)
        copied[f"results/core_gene_analysis/{name}"] = _sha256(src)

    reg = run_dir / "results/01_species_registry/species_registry.tsv"
    if reg.is_file():
        shutil.copy2(reg, dst / "results/01_species_registry/species_registry.tsv")
        copied["results/01_species_registry/species_registry.tsv"] = _sha256(reg)

    cfg = json.loads((run_dir / "run_config.json").read_text(encoding="utf-8"))
    trimmed = {k: cfg[k] for k in CONFIG_KEYS if k in cfg}
    trimmed["run_dir"] = dst.as_posix()
    trimmed["results_dir"] = (dst / "results").as_posix()
    (dst / "run_config.json").write_text(json.dumps(trimmed, indent=2) + "\n", encoding="utf-8")

    st = json.loads((run_dir / "status.json").read_text(encoding="utf-8"))
    (dst / "status.json").write_text(
        json.dumps({k: st[k] for k in STATUS_KEYS if k in st}, indent=2) + "\n", encoding="utf-8")

    src_idx = run_dir / "website_indices/coordinate_track_index.json"
    if not src_idx.is_file():
        src_idx = run_dir / "website_indices/generic/coordinate_track_index.json"
    if src_idx.is_file():
        idx = json.loads(src_idx.read_text(encoding="utf-8"))
        out = {"species": [{
            "species": sp.get("species"),
            "models": [{**{k: m.get(k) for k in MODEL_KEYS},
                        "blocks": [{k: b.get(k) for k in BLOCK_KEYS if k in b}
                                   for b in (m.get("blocks") or [])
                                   if b.get("feature_type") in (None, "coding_exon")]}
                       for m in (sp.get("models") or [])],
        } for sp in (idx.get("species") or [])]}
        (dst / "website_indices").mkdir(parents=True, exist_ok=True)
        (dst / "website_indices/coordinate_track_index.json").write_text(
            json.dumps(out, indent=1) + "\n", encoding="utf-8")
        copied["website_indices/coordinate_track_index.json"] = (
            f"reduced_from:{_sha256(src_idx)}")

    (dst / "fixture_provenance.json").write_text(json.dumps({
        "fixture_type": "boundary_dashboard_core_tables",
        "source_run_id": cfg.get("run_id"),
        "gene_symbol": cfg.get("gene_symbol"),
        "species_ids": cfg.get("species_ids"),
        "generated_by": "tests/fixtures/build_boundary_fixtures.py",
        "source_file_sha256": copied,
        "note": "Core tables copied verbatim; no scientific value is modified.",
    }, indent=2) + "\n", encoding="utf-8")

    blob = "\n".join(p.read_text(encoding="utf-8", errors="ignore")
                     for p in dst.rglob("*") if p.is_file())
    if "/Users/" in blob or "/home/" in blob:
        raise SystemExit(f"refusing to write {dst}: fixture still carries a personal path")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    build(Path(sys.argv[1]), Path(sys.argv[2]))
    print(f"wrote {sys.argv[2]}")
