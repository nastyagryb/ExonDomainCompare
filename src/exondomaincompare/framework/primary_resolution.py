#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List, Optional


class PrimaryResolutionError(RuntimeError):
    pass


def _rows(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")]


def _true(value: Any) -> bool:
    return str(value or "").strip().lower() in ("true", "1", "yes")


def resolve_primaries(core_dir: Path) -> Dict[str, Dict[str, Any]]:
    evidence = _rows(core_dir / "primary_selection_evidence.tsv")
    isoforms = _rows(core_dir / "protein_isoform_index.tsv")

    species = sorted({str(r.get("species_id") or "") for r in isoforms
                      if r.get("species_id")})
    # A legacy evidence table has no species column, so fall back to whatever the
    # isoform index knows; if that is empty too there is nothing to resolve.
    if not species:
        species = sorted({str(r.get("species_id") or "") for r in evidence
                          if r.get("species_id")})

    by_species: Dict[str, Dict[str, Any]] = {}
    for sid in species:
        pick: Optional[Dict[str, str]] = None
        source = ""

        pick = next((r for r in evidence
                     if str(r.get("species_id") or "") == sid
                     and _true(r.get("species_primary"))), None)
        if pick:
            source = "primary_selection_evidence.species_primary"

        if pick is None:
            pick = next((r for r in isoforms
                         if str(r.get("species_id") or "") == sid
                         and str(r.get("primary_status") or "").lower() == "primary"), None)
            if pick:
                source = "protein_isoform_index.primary_status"

        if pick is None and len(species) == 1:
            pick = next((r for r in evidence if _true(r.get("selected_primary"))), None)
            if pick:
                source = "primary_selection_evidence.selected_primary(single_species)"

        if pick is None or not (pick.get("protein_id") or "").strip():
            raise PrimaryResolutionError(
                f"no primary protein could be resolved for species '{sid}' in "
                f"{core_dir}. Checked species_primary and selected_primary in "
                f"primary_selection_evidence.tsv and primary_status in "
                f"protein_isoform_index.tsv. Refusing to guess: an arbitrary choice "
                f"produces a coordinate model of the wrong protein, whose domain and "
                f"boundary lookups silently return nothing."
            )

        length: Optional[int] = None
        raw_len = pick.get("protein_length") or pick.get("length_aa")
        try:
            length = int(float(raw_len)) if raw_len not in (None, "") else None
        except (TypeError, ValueError):
            length = None

        by_species[sid] = {
            "species_id": sid,
            "protein_id": (pick.get("protein_id") or "").strip(),
            "transcript_id": (pick.get("transcript_id") or "").strip(),
            "protein_length": length,
            "resolved_from": source,
        }
    return by_species
