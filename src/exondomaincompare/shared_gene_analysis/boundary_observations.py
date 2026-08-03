#!/usr/bin/env python3
"""Canonical exon–domain boundary observation table.

One row per species-specific Boundary observation, built from the canonical
per-species Boundary index inside ``protein_coordinate_model.json`` and — when
the run has more than one species — enriched with the comparative mapping from
the comparable-boundary groups.

The table therefore exists whenever the canonical JSON exists: it is never
required that somebody produced a TSV beforehand, and it is available for a
single-species run as well as for a multi-species run.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

# Column order is the download contract; do not reorder without updating the
# Data & Downloads report and the tests.
COLUMNS: Sequence[str] = (
    "species_id",
    "scientific_name",
    "primary_protein",
    "boundary_id",
    "comparable_boundary_group_id",
    "exon_transition",
    "native_aa_position",
    "MSA_column",
    "nearest_domain_instance_id",
    "nearest_domain_label",
    "nearest_edge",
    "signed_distance",
    "absolute_distance",
    "boundary_class",
    "mapping_method",
    "mapping_confidence",
    "status",
)

# Used wherever a value is genuinely not available, so an empty cell never has to
# be read as a biological statement.
UNAVAILABLE = "unavailable"


def _load(path: Path) -> Any:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def coordinate_model_path(run_dir: Path) -> Path:
    return Path(run_dir) / "website_indices" / "generic" / "protein_coordinate_model.json"


def _comparable_groups(run_dir: Path, model_doc: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Comparable-boundary groups from the canonical comparative sources.

    The comparative dataset index is preferred because it is the layer the UI and
    the packages read; the boundary dashboard inside the coordinate model is the
    fallback so the table also works before the comparative index was written.
    """
    comparative = _load(
        Path(run_dir) / "website_indices" / "generic"
        / "comparative_dataset_index.json") or {}
    groups = comparative.get("comparable_boundary_groups") or []
    if groups:
        return groups
    multi = (model_doc.get("boundary_dashboard") or {}).get("multi_species") or {}
    return multi.get("comparable_boundary_groups") or []


def _group_lookup(groups: Iterable[Dict[str, Any]]) -> Dict[tuple, Dict[str, Any]]:
    """(species_id, boundary_id) -> the group observation for that boundary."""
    out: Dict[tuple, Dict[str, Any]] = {}
    for g in groups:
        gid = g.get("comparable_boundary_group_id")
        members = (g.get("per_species_native_positions")
                   or g.get("members") or g.get("observations") or [])
        for obs in members:
            key = (obs.get("species_id"), obs.get("boundary_id"))
            if not key[1]:
                continue
            out[key] = {
                "comparable_boundary_group_id": gid,
                "msa_column": obs.get("msa_column") or g.get("msa_column"),
                "mapping_method": obs.get("mapping_method") or g.get("mapping_method"),
                "mapping_confidence": obs.get("mapping_confidence")
                if obs.get("mapping_confidence") is not None else g.get("confidence"),
                "mapping_status": obs.get("mapping_status") or g.get("mapping_status"),
            }
    return out


def build_rows(run_dir: Path,
               species_ids: Optional[Sequence[str]] = None) -> List[Dict[str, Any]]:
    """One row per Boundary observation of the requested species.

    ``species_ids`` of ``None`` means every species in the run.
    """
    doc = _load(coordinate_model_path(run_dir)) or {}
    models = doc.get("models") or doc.get("coordinate_models") or []
    if species_ids is not None:
        wanted = set(species_ids)
        models = [m for m in models if m.get("species_id") in wanted]

    lookup = _group_lookup(_comparable_groups(Path(run_dir), doc)) if len(
        doc.get("models") or []) > 1 else {}
    single_species = len(doc.get("models") or []) < 2

    rows: List[Dict[str, Any]] = []
    for model in models:
        sid = model.get("species_id") or ""
        protein = model.get("protein_id") or ""
        for b in model.get("exon_boundaries") or []:
            bid = b.get("boundary_id") or b.get("id") or ""
            mapped = lookup.get((sid, bid), {})
            signed = b.get("signed_distance")
            if signed is None:
                signed = b.get("signed_distance_aa")
            absolute = b.get("absolute_distance")
            if absolute is None:
                absolute = b.get("absolute_distance_aa")
            if absolute is None and signed is not None:
                absolute = abs(signed)
            rows.append({
                "species_id": sid,
                "scientific_name": model.get("scientific_name") or sid,
                "primary_protein": protein,
                "boundary_id": bid,
                "comparable_boundary_group_id":
                    mapped.get("comparable_boundary_group_id")
                    or ("not_applicable_single_species" if single_species
                        else "unmapped"),
                "exon_transition": b.get("label")
                or f"{b.get('left_exon_label') or '?'} → {b.get('right_exon_label') or '?'}",
                "native_aa_position": b.get("protein_position")
                if b.get("protein_position") is not None else b.get("boundary_position_aa"),
                "MSA_column": mapped.get("msa_column")
                if mapped.get("msa_column") is not None else UNAVAILABLE,
                "nearest_domain_instance_id":
                    b.get("nearest_domain_instance_id") or UNAVAILABLE,
                "nearest_domain_label": b.get("nearest_domain_short_label")
                or b.get("nearest_domain_label") or UNAVAILABLE,
                "nearest_edge": b.get("nearest_edge")
                or b.get("nearest_edge_type") or UNAVAILABLE,
                "signed_distance": signed if signed is not None else UNAVAILABLE,
                "absolute_distance": absolute if absolute is not None else UNAVAILABLE,
                "boundary_class": b.get("boundary_class") or b.get("class") or UNAVAILABLE,
                "mapping_method": mapped.get("mapping_method")
                or ("single_species_native" if single_species else "unmapped"),
                "mapping_confidence": mapped.get("mapping_confidence")
                if mapped.get("mapping_confidence") is not None else UNAVAILABLE,
                "status": mapped.get("mapping_status")
                or b.get("mapping_status") or "mapped",
            })
    return rows


def write_tsv(rows: Sequence[Dict[str, Any]], path: Path) -> Path:
    """Write the observation rows; the header is always the full contract."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(COLUMNS), delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: ("" if row.get(c) is None else row.get(c))
                             for c in COLUMNS})
    return path


def table_path(run_dir: Path, species_id: Optional[str] = None) -> Path:
    """Where the observation table for a run (or one of its species) lives."""
    base = Path(run_dir) / "results" / "generic_gene_analysis" / "boundaries"
    name = (f"exon_domain_boundaries__{species_id}.tsv" if species_id
            else "exon_domain_boundaries_long.tsv")
    return base / name


def ensure_table(run_dir: Path, species_id: Optional[str] = None) -> Optional[Path]:
    """Write the observation table from the canonical JSON and return its path.

    Returns ``None`` when the run has no Boundary observation at all, so a caller
    can report an exact reason instead of offering an empty download.
    """
    rows = build_rows(run_dir, [species_id] if species_id else None)
    if not rows:
        return None
    return write_tsv(rows, table_path(run_dir, species_id))


def label_for(multi_species: bool) -> str:
    """Return the visible download label."""
    return ("All species Boundary observations (TSV)" if multi_species
            else "Boundary observations (TSV)")


def main() -> int:  # pragma: no cover - operator entry point
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    rows = build_rows(args.run_dir)
    out = args.out or (args.run_dir / "results" / "generic_gene_analysis"
                       / "comparative" / "exon_domain_boundaries_long.tsv")
    write_tsv(rows, out)
    print(f"OK — {len(rows)} boundary observation(s) written to {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
