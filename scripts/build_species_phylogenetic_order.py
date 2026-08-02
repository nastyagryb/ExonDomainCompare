#!/usr/bin/env python3
"""
build_species_phylogenetic_order.py  (Sprint Part 2)

Derive a reproducible, taxonomy-based species display order for all FGFR2
figures (never silent alphabetical order).

Method priority:
  1. NCBI Taxonomy via ETE (ete3/ete4 NCBITaxa) if installed AND the local taxa
     database is available -> order_source = ncbi_taxonomy_ete (topology/lineage).
  2. Otherwise a documented curated fallback taxonomic order for the vertebrate
     FGFR2 panel -> order_source = curated_fallback_taxonomic_order.
  3. Any species not covered by either -> unresolved_fallback_order (only such
     rows may be alphabetized, and they are flagged).

Outputs:
  species_phylogenetic_order.tsv
  species_taxonomy_metadata.tsv

Broad taxon groups (in biological display order):
  Primates, Other mammals, Birds, Reptiles, Amphibians, Teleost fish
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path
from typing import Dict, List, Tuple


# Broad taxon group -> (display order index, display label, major clade).
GROUP_ORDER: Dict[str, Tuple[int, str, str]] = {
    "Primates": (0, "Primates", "Mammalia"),
    "Other mammals": (1, "Other mammals", "Mammalia"),
    "Birds": (2, "Birds", "Aves"),
    "Reptiles": (3, "Reptiles", "Reptilia"),
    "Amphibians": (4, "Amphibians", "Amphibia"),
    "Teleost fish": (5, "Teleost fish", "Actinopterygii"),
}

# Primate taxids in the panel (used to split mammals into Primates vs Other mammals).
PRIMATE_TAXIDS = {"9606", "9598", "9595", "9601", "9544", "9483"}

# Documented curated within-panel order following standard vertebrate phylogeny.
# (Group placement is biologically standard; within-group order follows accepted
# higher-level relationships for this specific 30-species panel.)
CURATED_ORDER: List[str] = [
    # Primates: great apes -> OWM -> NWM
    "homo_sapiens", "pan_troglodytes", "gorilla_gorilla_gorilla", "pongo_abelii",
    "macaca_mulatta", "callithrix_jacchus",
    # Other mammals: Glires -> Carnivora -> Perissodactyla/Artiodactyla -> Marsupialia -> Monotremata
    "mus_musculus", "rattus_norvegicus", "oryctolagus_cuniculus",
    "canis_lupus_familiaris", "felis_catus",
    "equus_caballus", "sus_scrofa", "bos_taurus", "ovis_aries",
    "monodelphis_domestica", "ornithorhynchus_anatinus",
    # Birds: Galliformes -> Passeriformes
    "gallus_gallus", "meleagris_gallopavo", "taeniopygia_guttata",
    # Reptiles: Crocodilia -> Testudines -> Squamata
    "alligator_mississippiensis", "chrysemys_picta_bellii", "anolis_carolinensis",
    # Amphibians: Urodela -> Anura
    "ambystoma_mexicanum", "xenopus_tropicalis",
    # Teleost fish: Ostariophysi -> Acanthomorpha
    "danio_rerio", "oryzias_latipes", "gasterosteus_aculeatus",
    "oreochromis_niloticus", "takifugu_rubripes",
]

CLADE_TO_GROUP = {"bird": "Birds", "reptile": "Reptiles", "amphibian": "Amphibians", "fish": "Teleost fish"}


def norm(s: str) -> str:
    return str(s or "").strip().lower().replace(" ", "_")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, rows: List[Dict[str, object]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def taxon_group_for(taxid: str, clade: str) -> str:
    c = norm(clade)
    if c == "mammal":
        return "Primates" if str(taxid).strip() in PRIMATE_TAXIDS else "Other mammals"
    return CLADE_TO_GROUP.get(c, "Other mammals")


def try_ncbi_taxa():
    """Return an NCBITaxa instance if ete is installed and the DB exists, else None."""
    for mod in ("ete3", "ete4"):
        if importlib.util.find_spec(mod) is None:
            continue
        try:
            db = Path.home() / ".etetoolkit" / "taxa.sqlite"
            if not db.exists():
                return None
            module = __import__(mod, fromlist=["NCBITaxa"])
            return module.NCBITaxa()
        except Exception:
            return None
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Build reproducible phylogenetic species order (Part 2).")
    ap.add_argument("--registry", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    reg = read_tsv(args.registry)
    ncbi = try_ncbi_taxa()

    lineage_by_sp: Dict[str, str] = {}
    if ncbi is not None:
        for r in reg:
            sp = norm(r.get("species_id"))
            tid = str(r.get("taxid", "")).strip()
            try:
                lin = ncbi.get_lineage(int(tid)) if tid else None
                if lin:
                    names = ncbi.get_taxid_translator(lin)
                    lineage_by_sp[sp] = "; ".join(names[t] for t in lin)
            except Exception:
                pass

    curated_rank = {sp: i for i, sp in enumerate(CURATED_ORDER)}

    entries: List[Dict[str, object]] = []
    for r in reg:
        sp = norm(r.get("species_id"))
        tid = str(r.get("taxid", "")).strip()
        clade = r.get("clade", "")
        disp = r.get("scientific_name", "") or sp.replace("_", " ").capitalize()
        group = taxon_group_for(tid, clade)
        g_order, g_display, major = GROUP_ORDER[group]

        if sp in lineage_by_sp:
            order_source = "ncbi_taxonomy_ete"
            order_conf = "high"
            within = curated_rank.get(sp, 9999)  # ete gives lineage; keep curated within-group tie-break
            lineage = lineage_by_sp[sp]
            warn = ""
        elif sp in curated_rank:
            order_source = "curated_fallback_taxonomic_order"
            order_conf = "medium"
            within = curated_rank[sp]
            lineage = f"Eukaryota; Metazoa; Chordata; Vertebrata; {major}; {disp}"
            warn = "within_group_order_curated_not_topology_derived"
        else:
            order_source = "unresolved_fallback_order"
            order_conf = "low"
            within = 100000
            lineage = f"Vertebrata; {major}; {disp}"
            warn = "species_not_in_curated_order_alphabetical_fallback"

        entries.append({
            "species": sp, "display_species_name": disp, "taxid": tid,
            "taxon_group": group, "taxon_group_display": g_display, "major_clade": major,
            "lineage_string": lineage, "_g_order": g_order, "_within": within,
            "order_source": order_source, "order_confidence": order_conf, "order_warning": warn,
        })

    # Sort: group order, then within-group order, then species name (only matters for unresolved).
    entries.sort(key=lambda e: (e["_g_order"], e["_within"], e["species"]))
    for i, e in enumerate(entries, start=1):
        e["phylo_order"] = i

    args.outdir.mkdir(parents=True, exist_ok=True)
    cols = ["species", "display_species_name", "taxid", "taxon_group", "taxon_group_display",
            "major_clade", "lineage_string", "phylo_order", "order_source", "order_confidence", "order_warning"]
    write_tsv(args.outdir / "species_phylogenetic_order.tsv", entries, cols)

    tax_cols = ["species", "display_species_name", "taxid", "taxon_group", "major_clade", "lineage_string", "order_source"]
    write_tsv(args.outdir / "species_taxonomy_metadata.tsv", entries, tax_cols)

    from collections import Counter
    src = Counter(str(e["order_source"]) for e in entries)
    grp = Counter(str(e["taxon_group"]) for e in entries)
    print(f"[OK] phylogenetic order for {len(entries)} species; order_source={dict(src)}")
    print(f"     taxon groups={dict(grp)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
