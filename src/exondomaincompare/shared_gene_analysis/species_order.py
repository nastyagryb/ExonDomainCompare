"""One canonical species order for every comparative view.

Comparative figures are read by scanning down a column, so the row order carries
meaning: species next to each other should be related. Before this module each
view picked its own order — the boundary matrix used the reference list, the
gallery used whatever the index happened to yield, the selectors were
alphabetical — and the same two species swapped places between two figures of
the same dataset.

Two orderings exist, and the difference is not cosmetic:

* **taxonomic** — species are grouped by clade and ordered within a clade by a
  curated arrangement. This is what this module produces, and it is what the
  outputs are labelled.
* **phylogenetic** — species are ordered by the tips of an actual tree. It
  requires a tree, computed or supplied. No tree is used here, so nothing
  produced by this module may be described that way.

The validated 30-species panel keeps the approved order of
``reference/Species_list_final_30.txt`` exactly. Species outside that panel are
placed by clade and then alphabetically, which is a stated rule rather than an
accident of dictionary iteration.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from exondomaincompare.config import discover_repository_root

REPO = discover_repository_root(__file__)
REFERENCE_PANEL = REPO / "reference" / "Species_list_final_30.txt"

#: Clade order for species outside the reference panel: the same broad shape the
#: reference panel uses, from mammals outward.
CLADE_ORDER = ["mammal", "bird", "reptile", "amphibian", "fish", "invertebrate",
               "other"]

#: Display names for the clade filter.
CLADE_LABELS = {
    "mammal": "Mammals",
    "bird": "Birds",
    "reptile": "Reptiles",
    "amphibian": "Amphibians",
    "fish": "Fishes",
    "invertebrate": "Invertebrates",
    "other": "Other",
}

#: The coarser taxon grouping the FGFR2 figures already use for their sidebars.
_PRIMATES = {
    "homo_sapiens", "pan_troglodytes", "gorilla_gorilla_gorilla", "pongo_abelii",
    "macaca_mulatta", "callithrix_jacchus",
}
TAXON_GROUPS = {
    "mammal": "Other mammals",
    "bird": "Birds",
    "reptile": "Reptiles",
    "amphibian": "Amphibians",
    "fish": "Teleost fish",
    "invertebrate": "Invertebrates",
    "other": "Other",
}

#: The order those groups appear in, matching the reference panel.
TAXON_GROUP_ORDER = ["Primates", "Other mammals", "Birds", "Reptiles",
                     "Amphibians", "Teleost fish", "Invertebrates", "Other"]

#: Lineage strings, as far as the offline registry can state them honestly. Only
#: the ranks the clade actually determines are listed; nothing is guessed.
_LINEAGE = {
    "mammal": "Eukaryota; Metazoa; Chordata; Vertebrata; Mammalia",
    "bird": "Eukaryota; Metazoa; Chordata; Vertebrata; Aves",
    "reptile": "Eukaryota; Metazoa; Chordata; Vertebrata; Reptilia",
    "amphibian": "Eukaryota; Metazoa; Chordata; Vertebrata; Amphibia",
    "fish": "Eukaryota; Metazoa; Chordata; Vertebrata; Actinopterygii",
    "invertebrate": "Eukaryota; Metazoa",
    "other": "Eukaryota",
}

_REGISTRY: Optional[Dict[str, Dict[str, Any]]] = None


def _registry() -> Dict[str, Dict[str, Any]]:
    """The offline species registry, keyed by species_id."""
    global _REGISTRY
    if _REGISTRY is None:
        import sys
        if str(REPO / "scripts") not in sys.path:
            sys.path.insert(0, str(REPO / "scripts"))
        try:
            from build_species_registry_improved import KNOWN_SPECIES  # type: ignore
        except Exception:  # pragma: no cover - registry is optional
            KNOWN_SPECIES = {}
        _REGISTRY = {
            (entry.get("ensembl_species") or name.lower().replace(" ", "_")): {
                "scientific_name": entry.get("ncbi_species") or name,
                "common_name": entry.get("common_name") or "",
                "taxid": str(entry.get("taxid") or ""),
                "clade": entry.get("clade") or "other",
            }
            for name, entry in KNOWN_SPECIES.items()
        }
    return _REGISTRY


def reference_panel_order() -> Dict[str, int]:
    """The approved order of the validated 30-species panel."""
    order: Dict[str, int] = {}
    if REFERENCE_PANEL.is_file():
        for i, line in enumerate(REFERENCE_PANEL.read_text(encoding="utf-8").splitlines()):
            sid = line.strip().lower().replace(" ", "_")
            if sid:
                order[sid] = i
    return order


def scientific_name(species_id: str) -> str:
    """`gallus_gallus` -> `Gallus gallus`; only the genus is capitalised."""
    known = _registry().get(species_id)
    if known and known["scientific_name"]:
        return known["scientific_name"]
    parts = [p for p in str(species_id or "").replace(" ", "_").split("_") if p]
    if not parts:
        return species_id
    return " ".join([parts[0].capitalize(), *[p.lower() for p in parts[1:]]])


def clade_of(species_id: str) -> str:
    return _registry().get(species_id, {}).get("clade", "other")


def taxon_group(species_id: str) -> str:
    if species_id in _PRIMATES:
        return "Primates"
    return TAXON_GROUPS.get(clade_of(species_id), "Other")


def species_record(species_id: str, display_order: int,
                   ordering_method: str) -> Dict[str, Any]:
    """One row of the canonical order."""
    known = _registry().get(species_id, {})
    clade = known.get("clade", "other")
    return {
        "species_id": species_id,
        "scientific_name": scientific_name(species_id),
        "common_name": known.get("common_name", ""),
        "ncbi_taxonomy_id": known.get("taxid", ""),
        "taxonomic_lineage": _LINEAGE.get(clade, _LINEAGE["other"]),
        "major_clade": clade,
        "clade_label": CLADE_LABELS.get(clade, "Other"),
        "taxon_group": taxon_group(species_id),
        "display_order": display_order,
        # Populated only when a real tree is supplied; an empty value states that
        # no tree backs this ordering.
        "tree_tip_id": "",
        "ordering_method": ordering_method,
    }


def order_species(species_ids: Iterable[str]) -> List[str]:
    """Sort species into the canonical order.

    Members of the validated reference panel keep their approved positions and
    come first. Everything else follows, grouped by clade and alphabetical
    within a clade.
    """
    panel = reference_panel_order()
    unique = list(dict.fromkeys(s for s in species_ids if s))

    def key(species_id: str):
        if species_id in panel:
            return (0, panel[species_id], "")
        clade = clade_of(species_id)
        rank = CLADE_ORDER.index(clade) if clade in CLADE_ORDER else len(CLADE_ORDER)
        return (1, rank, scientific_name(species_id))

    return sorted(unique, key=key)


def build_species_order(species_ids: Sequence[str],
                        ordering_method: str = "taxonomic") -> Dict[str, Any]:
    """The canonical order for one dataset, as a serialisable document."""
    if ordering_method not in ("taxonomic", "phylogenetic"):
        raise ValueError("ordering_method must be 'taxonomic' or 'phylogenetic'")
    if ordering_method == "phylogenetic":
        raise ValueError(
            "A phylogenetic order requires a supplied or computed tree. This "
            "builder derives its order from taxonomy, so it may only be "
            "labelled 'taxonomic'.")
    ordered = order_species(species_ids)
    panel = reference_panel_order()
    rows = [species_record(sid, i, ordering_method) for i, sid in enumerate(ordered)]
    clades = list(dict.fromkeys(r["major_clade"] for r in rows))
    return {
        "contract": "species_order_v1",
        "ordering_method": ordering_method,
        "ordering_basis": (
            "Curated taxonomic arrangement. Species of the validated 30-species "
            "reference panel keep its approved order; any further species follow, "
            "grouped by clade and alphabetical within a clade. No phylogenetic "
            "tree is used, so this order is taxonomic and not phylogenetic."),
        "reference_panel_species": sum(1 for sid in ordered if sid in panel),
        "n_species": len(rows),
        "clades_present": [
            {"clade": c, "label": CLADE_LABELS.get(c, "Other"),
             "species": [r["species_id"] for r in rows if r["major_clade"] == c]}
            for c in CLADE_ORDER if c in clades
        ],
        "species": rows,
    }


TSV_COLUMNS = ["species_id", "scientific_name", "common_name", "ncbi_taxonomy_id",
               "taxonomic_lineage", "major_clade", "display_order", "tree_tip_id",
               "ordering_method"]


def write_species_order(species_ids: Sequence[str], outdir: Path,
                        ordering_method: str = "taxonomic") -> Dict[str, Path]:
    """Write ``species_order.json`` and ``species_order.tsv`` into ``outdir``."""
    doc = build_species_order(species_ids, ordering_method)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / "species_order.json"
    json_path.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n",
                         encoding="utf-8")
    tsv_path = outdir / "species_order.tsv"
    with tsv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(doc["species"])
    return {"json": json_path, "tsv": tsv_path}
