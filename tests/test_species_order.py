"""One canonical species order across every comparative view.

Each view used to order species for itself — alphabetically in the gallery
figures, by whatever the index returned in the matrices, by insertion order in
the selectors. The same two species therefore swapped places between two figures
of the same dataset, which makes a comparative reading unreliable: a reader
scanning down a column has to re-check the labels for every panel.

The order is taxonomic, not phylogenetic. No tree is computed or supplied, so
these tests also pin the wording: nothing may be presented as phylogenetic.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from shared_gene_analysis import species_order as so  # noqa: E402

FRONTEND = ROOT / "webapp" / "frontend" / "src"
JS_MODULE = FRONTEND / "pages" / "viewers" / "speciesOrder.js"

PANEL = ["homo_sapiens", "mus_musculus", "gallus_gallus", "anolis_carolinensis",
         "xenopus_tropicalis", "danio_rerio"]


# --------------------------------------------------------------------------- #
# The order itself
# --------------------------------------------------------------------------- #
def test_the_validated_reference_panel_keeps_its_approved_order():
    panel = so.reference_panel_order()
    assert len(panel) == 30
    ordered = so.order_species(list(panel))
    assert ordered == sorted(panel, key=panel.get)
    assert ordered[0] == "homo_sapiens"
    assert ordered[-1] == "oreochromis_niloticus"


def test_species_are_grouped_by_clade_not_alphabetically():
    """Alphabetically Danio rerio lands between Callithrix and Equus."""
    ids = ["equus_caballus", "danio_rerio", "callithrix_jacchus", "gallus_gallus"]
    ordered = so.order_species(ids)
    assert ordered == ["callithrix_jacchus", "equus_caballus", "gallus_gallus",
                       "danio_rerio"]
    assert ordered != sorted(ids)


def test_a_species_outside_the_panel_is_placed_by_clade_then_name():
    ids = ["mus_musculus", "vulpes_vulpes", "danio_rerio", "acipenser_ruthenus"]
    ordered = so.order_species(ids)
    # Panel members first in their approved order, then the rest by clade.
    assert ordered[:2] == ["mus_musculus", "danio_rerio"]
    assert set(ordered[2:]) == {"vulpes_vulpes", "acipenser_ruthenus"}


def test_the_order_is_deterministic_for_completely_unknown_species():
    ids = ["zzz_species", "aaa_species", "mmm_species"]
    assert so.order_species(ids) == so.order_species(list(reversed(ids)))


def test_a_single_species_dataset_orders_without_error():
    assert so.order_species(["gallus_gallus"]) == ["gallus_gallus"]
    doc = so.build_species_order(["gallus_gallus"])
    assert doc["n_species"] == 1
    assert doc["species"][0]["display_order"] == 0


def test_duplicates_collapse():
    assert so.order_species(["mus_musculus", "mus_musculus"]) == ["mus_musculus"]


# --------------------------------------------------------------------------- #
# What the order document has to say about itself
# --------------------------------------------------------------------------- #
def test_every_required_field_is_present():
    doc = so.build_species_order(PANEL)
    for row in doc["species"]:
        for field in ("species_id", "scientific_name", "common_name",
                      "ncbi_taxonomy_id", "taxonomic_lineage", "major_clade",
                      "display_order", "tree_tip_id", "ordering_method"):
            assert field in row, field
    assert [r["display_order"] for r in doc["species"]] == list(range(len(PANEL)))


def test_real_taxonomy_ids_are_carried_through():
    doc = so.build_species_order(["homo_sapiens", "danio_rerio"])
    ids = {r["species_id"]: r["ncbi_taxonomy_id"] for r in doc["species"]}
    assert ids["homo_sapiens"] == "9606"
    assert ids["danio_rerio"] == "7955"


def test_the_ordering_is_declared_taxonomic_and_never_phylogenetic():
    doc = so.build_species_order(PANEL)
    assert doc["ordering_method"] == "taxonomic"
    assert "not phylogenetic" in doc["ordering_basis"]
    # No tree backs this order, so no row may claim a tree tip.
    assert all(r["tree_tip_id"] == "" for r in doc["species"])


def test_claiming_a_phylogenetic_order_without_a_tree_is_refused():
    with pytest.raises(ValueError, match="tree"):
        so.build_species_order(PANEL, ordering_method="phylogenetic")


def test_the_document_lists_the_clades_present_in_order():
    doc = so.build_species_order(PANEL)
    assert [c["clade"] for c in doc["clades_present"]] == \
        ["mammal", "bird", "reptile", "amphibian", "fish"]


def test_writing_the_order_produces_both_files(tmp_path):
    paths = so.write_species_order(PANEL, tmp_path)
    assert paths["json"].is_file() and paths["tsv"].is_file()
    header = paths["tsv"].read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert header == so.TSV_COLUMNS
    assert len(paths["tsv"].read_text(encoding="utf-8").splitlines()) == len(PANEL) + 1


# --------------------------------------------------------------------------- #
# One definition, used everywhere
# --------------------------------------------------------------------------- #
def test_the_frontend_mirror_matches_the_backend_order():
    """The JS copy of the panel cannot drift from the reference list."""
    js = JS_MODULE.read_text(encoding="utf-8")
    block = js.split("const REFERENCE_PANEL = [", 1)[1].split("];", 1)[0]
    js_ids = re.findall(r'\["([a-z_]+)",\s*"([a-z]+)"\]', block)
    panel = so.reference_panel_order()
    assert [sid for sid, _ in js_ids] == sorted(panel, key=panel.get)
    for sid, clade in js_ids:
        assert clade == so.clade_of(sid), sid


def test_the_frontend_mirror_matches_the_backend_clade_order():
    js = JS_MODULE.read_text(encoding="utf-8")
    order = re.search(r"export const CLADE_ORDER = \[(.*?)\];", js, re.S).group(1)
    assert re.findall(r'"([a-z]+)"', order) == so.CLADE_ORDER


def test_no_comparative_view_sorts_species_alphabetically_any_more():
    viewers = FRONTEND / "pages" / "viewers"
    offenders = []
    for path in sorted(list(viewers.glob("*.js")) + list(viewers.glob("*.jsx"))):
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if "localeCompare" in line and re.search(r"species", line, re.I):
                offenders.append(f"{path.name}:{i}: {line.strip()}")
    assert not offenders, ("these views still order species alphabetically instead "
                           f"of canonically: {offenders}")


def test_the_shared_frontend_helper_is_the_one_used_by_the_gallery_figures():
    text = (FRONTEND / "pages/viewers/comparativeGalleryFigures.js").read_text(encoding="utf-8")
    assert 'from "./speciesOrder.js"' in text
    assert "return orderSpeciesRows(rows);" in text


def test_the_frontend_helper_states_the_ordering_is_taxonomic():
    js = JS_MODULE.read_text(encoding="utf-8")
    assert 'export const ORDERING_METHOD = "taxonomic"' in js
    assert "not a phylogenetic tree order" in js


def test_the_indices_carry_the_order_so_the_frontend_need_not_recompute_it():
    comparative = (ROOT / "src/exondomaincompare/shared_gene_analysis"
                   / "comparative_dataset.py").read_text(encoding="utf-8")
    assert '"species_order": order_doc' in comparative
    assert '"species_order.tsv"' in comparative
    synteny = (ROOT / "src/exondomaincompare/shared_gene_analysis/indices"
               / "synteny_locus.py").read_text(encoding="utf-8")
    assert "so.order_species(by_species.keys())" in synteny
    assert '"species_order": so.build_species_order' in synteny

    legacy = (ROOT / "scripts/shared_gene_analysis"
              / "comparative_dataset.py").read_text(encoding="utf-8")
    assert "exondomaincompare.shared_gene_analysis.comparative_dataset" in legacy


def test_the_validated_fgfr2_builder_uses_the_same_definition():
    text = (ROOT / "scripts" / "build_website_indices.py").read_text(encoding="utf-8")
    assert "_species_order.reference_panel_order()" in text
    assert "_species_order.TAXON_GROUP_ORDER" in text
    # The old private copy of the reference list is gone.
    assert 'reference" / "Species_list_final_30.txt"' not in text


def test_nothing_the_user_reads_calls_this_order_phylogenetic():
    """The word may be defined and disclaimed, but never used as the label."""
    doc = so.build_species_order(PANEL)
    emitted = [doc["ordering_method"],
               *(r["ordering_method"] for r in doc["species"]),
               *(c["label"] for c in doc["clades_present"])]
    assert all("phylogenetic" not in str(v).lower() for v in emitted)
    # The word appears in the basis text only to rule the claim out.
    basis = doc["ordering_basis"].lower()
    assert "no phylogenetic tree is used" in basis
    assert "taxonomic and not phylogenetic" in basis
