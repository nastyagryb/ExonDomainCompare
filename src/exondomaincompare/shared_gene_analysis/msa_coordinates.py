#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

GAP_CHARS = "-."


def read_aligned_fasta(path: Path) -> List[Tuple[str, str]]:
    entries: List[Tuple[str, List[str]]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                entries.append((line[1:].strip(), []))
            elif entries:
                entries[-1][1].append(line)
    return [(h, "".join(parts)) for h, parts in entries]


def parse_header(header: str) -> Dict[str, str]:
    head = header.strip()
    species = ""
    if "|" in head:
        head, _, species = head.rpartition("|")
    protein = head.split()[0] if head.split() else ""
    gene = ""
    rest = head.split()[1:] if len(head.split()) > 1 else []
    if rest:
        gene = rest[0]
    return {"protein_id": protein.strip(), "gene_symbol": gene.strip(),
            "species_id": species.strip()}


def column_map(aligned: str) -> Dict[int, int]:
    out: Dict[int, int] = {}
    native = 0
    for col, ch in enumerate(aligned, start=1):
        if ch in GAP_CHARS:
            continue
        native += 1
        out[native] = col
    return out


def build_msa_coordinate_map(alignment_path: Path) -> Dict[str, Any]:
    path = Path(alignment_path)
    if not path.is_file():
        # Only the file name, never the absolute path: this dict is published in the
        # run index, where a machine-local path would travel with the results.
        return {"available": False,
                "reason": (f"no cross-species primary alignment found ({path.name}); "
                           f"expected for a single-species run"),
                "n_columns": 0, "species": []}

    records = read_aligned_fasta(path)
    if len(records) < 2:
        return {"available": False,
                "reason": ("a cross-species mapping needs at least two aligned "
                           f"primaries; found {len(records)}"),
                "n_columns": 0, "species": []}

    widths = {len(seq) for _, seq in records}
    if len(widths) != 1:
        return {"available": False,
                "reason": f"aligned sequences differ in width: {sorted(widths)}",
                "n_columns": 0, "species": []}
    n_columns = widths.pop()

    species: List[Dict[str, Any]] = []
    for header, seq in records:
        meta = parse_header(header)
        cmap = column_map(seq)
        species.append({
            "species_id": meta["species_id"],
            "protein_id": meta["protein_id"],
            "protein_length": len(cmap),
            "n_gap_columns": sum(1 for ch in seq if ch in GAP_CHARS),
            "native_to_column": [[n, c] for n, c in sorted(cmap.items())],
        })

    duplicate = len({s["species_id"] for s in species}) != len(species)
    if duplicate or any(not s["species_id"] for s in species):
        return {"available": False,
                "reason": ("every aligned sequence must carry a distinct species id; "
                           f"got {[s['species_id'] for s in species]}"),
                "n_columns": n_columns, "species": []}

    return {
        "available": True,
        "reason": ("Cross-species primary-protein alignment; one primary per species. "
                   "A shared column means the residues were aligned, not that they are "
                   "functionally equivalent."),
        "alignment_file": path.name,
        "n_columns": n_columns,
        "n_species": len(species),
        "species": species,
    }


def lookup_tables(coord_map: Dict[str, Any]) -> Dict[str, Dict[int, int]]:
    return {
        s["species_id"]: {int(n): int(c) for n, c in s.get("native_to_column") or []}
        for s in coord_map.get("species") or []
    }


def annotate_boundaries_with_columns(models: Sequence[Dict[str, Any]],
                                     coord_map: Dict[str, Any]) -> Dict[str, Any]:
    tables = lookup_tables(coord_map)
    aligned_protein = {s["species_id"]: s["protein_id"]
                       for s in coord_map.get("species") or []}
    report: Dict[str, Any] = {}

    for m in models:
        sid = m.get("species_id") or ""
        table = tables.get(sid) or {}
        expected = aligned_protein.get(sid)
        # Mapping a boundary through an alignment of a *different* protein would put it
        # at a plausible but wrong column, which is the failure mode this whole phase
        # exists to remove.
        protein_matches = bool(expected) and m.get("protein_id") == expected
        mapped = 0
        total = 0
        for b in m.get("exon_boundaries") or []:
            total += 1
            pos = b.get("protein_position")
            if pos is None:
                pos = b.get("boundary_position_aa")
            col = table.get(int(pos)) if (protein_matches and pos is not None) else None
            b["msa_column"] = col
            b["msa_mapping_status"] = (
                "mapped" if col is not None else
                "unmapped_protein_not_aligned" if not protein_matches else
                "unmapped_position_outside_alignment"
            )
            if col is not None:
                mapped += 1
        report[sid] = {
            "species_id": sid,
            "protein_id": m.get("protein_id"),
            "aligned_protein_id": expected,
            "protein_matches_alignment": protein_matches,
            "boundaries_total": total,
            "boundaries_mapped": mapped,
            "mapping_status": (
                "available" if mapped else
                "unavailable_protein_not_aligned" if not protein_matches else
                "unavailable_no_positions_mapped"
            ),
        }
    return report
