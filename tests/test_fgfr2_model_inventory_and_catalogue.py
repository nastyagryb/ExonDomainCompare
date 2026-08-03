"""The FGFR2 model inventory and the curated Gallery catalogue.

Two claims are under test here, and both were previously wrong in ways that a
reader could not have detected from the interface.

The first is arithmetic. The dataset was described as "58 models = 30 species ×
IIIb/IIIc", which cannot be true because 30 × 2 is 60. The two missing combinations
are not a rounding detail: they are the two proteins the freeze declined to admit to
the primary set, each with a recorded reason, and a description that multiplies them
away is a description that hides them.

The second is that a rendered file is not a Gallery card. 287 derived files became
90-odd cards, so a reader looking for the cassette evidence had to scroll past sixty
near-identical per-species thumbnails and four filter variants of one figure. The
tests below pin the curated shape: what a reader actually sees, per scope.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "shared_gene_analysis"))

DERIVED = ROOT / "results" / "derived" / "example"
INDICES = DERIVED / "website_indices"
MODEL_INDEX = INDICES / "protein_coordinate_model.json"
CATALOGUE = INDICES / "figure_catalogue.json"
GALLERY_INDEX = INDICES / "figure_index.json"
INVENTORY = INDICES / "tables" / "fgfr2_model_inventory.tsv"
FREEZE = ROOT / "results" / "final_30_until_interpro_prepare"

#: What the audit found in the freeze. These are the numbers the reports must use.
N_SPECIES = 30
N_ISOFORMS = 2
N_EXPECTED_COMBINATIONS = 60
N_MODELS = 58
N_IIIB_MODELS = 29
N_IIIC_MODELS = 29
N_SPECIES_WITH_BOTH = 28
N_SPECIES_WITH_ONE = 2
MISSING_COMBINATIONS = {("canis_lupus_familiaris", "IIIc"), ("pongo_abelii", "IIIb")}
CASSETTE_ONLY_MODELS = {"fgfr2:canis_lupus_familiaris:IIIb"}

#: Thirteen comparative main cards: the eleven curated scientific figures that
#: survived the withdrawal of the framework evidence stack and the synteny
#: neighbour-conservation matrix, plus the three shared whole-protein
#: exon–domain boundary figures.
N_COMPARATIVE_MAIN = 13
N_SUPPLEMENTS = 11
SPECIES_SCOPE_CARD_COUNTS = {6, 7}

#: Files that were separate Gallery cards and are now modes of one card. Each is
#: either byte-identical to the card's default view or a filter of it.
MERGED_SOURCE_FIGURES = {
    "Figure_9_FGFR2_local_synteny_neighborhood",       # byte-identical to Figure 9A
    "Supplement_full_length_MSA_outliers",             # byte-identical to QC histograms
    "Figure_10A_IIIb_exon_domain_architecture_primary",     # cassette filter
    "Figure_10B_IIIc_exon_domain_architecture_primary",     # cassette filter
    "Figure_10C_mammals_exon_domain_architecture_primary",  # taxon filter
    "Figure_10D_nonmammals_exon_domain_architecture_primary",
    "Figure_3C_exon_to_protein_cassette_coordinate_map",    # second view of Figure 3
    "Figure_6B_species_resolved_IIIb_IIIc_cassette_residue_map",
}


def _json(path: Path):
    if not path.is_file():
        pytest.skip(f"{path.relative_to(ROOT)} not built")
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def model_index():
    return _json(MODEL_INDEX)


@pytest.fixture(scope="module")
def catalogue():
    return _json(CATALOGUE)


@pytest.fixture(scope="module")
def inventory():
    if not INVENTORY.is_file():
        pytest.skip("model inventory not built")
    with INVENTORY.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --------------------------------------------------------------------------- #
# 1. the inventory
# --------------------------------------------------------------------------- #
def test_the_model_inventory_totals_are_exact(model_index):
    a = model_index["availability"]
    assert a["n_species"] == N_SPECIES
    assert a["n_models"] == N_MODELS
    assert a["n_expected_combinations"] == N_EXPECTED_COMBINATIONS
    assert a["n_models_per_isoform"] == {"IIIb": N_IIIB_MODELS, "IIIc": N_IIIC_MODELS}
    assert a["n_species_with_both_models"] == N_SPECIES_WITH_BOTH
    assert a["n_species_with_one_model"] == N_SPECIES_WITH_ONE


def test_the_58_models_are_not_explained_as_30_times_2(model_index):
    """58 is 60 expected combinations minus 2 the freeze did not admit."""
    a = model_index["availability"]
    assert a["n_species"] * N_ISOFORMS == N_EXPECTED_COMBINATIONS
    assert a["n_models"] + len(a["unavailable_combinations"]) == N_EXPECTED_COMBINATIONS
    explanation = a["explanation"]
    assert "60 expected" in explanation
    assert "58 have an architecture model" in explanation
    # The claim that would be false.
    assert "30 species × 2" not in explanation
    assert "30 x 2" not in explanation.lower()


def test_the_missing_combinations_are_named_with_their_reason(model_index):
    unavailable = model_index["availability"]["unavailable_combinations"]
    assert {(u["species_id"], u["isoform"]) for u in unavailable} == MISSING_COMBINATIONS
    for u in unavailable:
        assert u["omission_reason"], u
        assert u["readiness_class"] == "supplement_review_only"
        assert u["species_still_represented"] is True


def test_a_model_without_an_exon_series_says_which_layers_it_cannot_support(model_index):
    """One protein has a validated cassette and domains but no coding-exon series.

    Reporting it as a complete model would put an empty exon figure in the Gallery;
    dropping it would lose a validated architecture. It is a real model that cannot
    support two of the layers, and it says so.
    """
    cassette_only = [m for m in model_index["models"]
                     if m["availability_status"] == "cassette_only_no_exon_series"]
    assert {m["model_id"] for m in cassette_only} == CASSETTE_ONLY_MODELS
    for m in cassette_only:
        assert set(m["unavailable_layers"]) == {"exon_structure",
                                               "exon_domain_boundaries"}
        assert m["unavailable_reason"]
        assert m["representative_domains"], "the domain layer is still real"


def test_all_thirty_species_stay_in_the_dataset(model_index, catalogue, inventory):
    assert len(model_index["species_scope"]) == N_SPECIES
    assert len({r["species_id"] for r in inventory}) == N_SPECIES
    assert len(catalogue["species_scopes"]) == N_SPECIES
    # Including the two that are missing one isoform model.
    for species, _ in MISSING_COMBINATIONS:
        assert species in catalogue["species_scopes"]
        assert catalogue["species_scopes"][species]["n_models"] == 1


def test_the_inventory_table_has_a_row_per_expected_combination(inventory):
    assert len(inventory) == N_EXPECTED_COMBINATIONS
    required = ["species_id", "scientific_name", "model_id", "protein_id",
                "transcript_id", "isoform_label", "model_role", "is_primary_reference",
                "availability_status", "reconstruction_status", "review_status",
                "omission_reason"]
    assert list(inventory[0].keys()) == required
    absent = [r for r in inventory if r["availability_status"] == "no_architecture_model"]
    assert {(r["species_id"], r["isoform_label"]) for r in absent} == MISSING_COMBINATIONS
    for row in absent:
        assert row["omission_reason"], row
        assert not row["model_id"], "an absent combination has no model to identify"


def test_the_inventory_records_how_each_model_was_reconstructed(model_index):
    counts = model_index["availability"]["n_models_by_reconstruction_status"]
    assert sum(counts.values()) == N_MODELS
    assert counts["cassette_only_high_confidence"] == 1
    assert counts["native_exon_blocks_reconstructed"] == 2
    assert counts["coordinate_mapped"] == 55


def test_review_flags_are_carried_not_flattened(model_index):
    flagged = model_index["availability"]["review_flagged_models"]
    assert flagged, "the freeze raised review flags; they must survive into the model"
    for m in model_index["models"]:
        assert m["review_status"]
        if m["model_id"] in flagged:
            assert m["review_status"] != "no_review_flag"


# --------------------------------------------------------------------------- #
# 2. the model-role hierarchy
# --------------------------------------------------------------------------- #
def test_every_model_states_an_explicit_role_and_id(model_index):
    from exondomaincompare.shared_gene_analysis import model_roles

    for m in model_index["models"]:
        assert m["model_id"], m["protein_id"]
        assert model_roles.known_role(m["model_role"]), m["model_role"]
    assert not model_roles.role_errors(model_index["models"])


def test_exactly_one_primary_reference_per_species(model_index):
    by_species = {}
    for m in model_index["models"]:
        by_species.setdefault(m["species_id"], []).append(m)
    assert len(by_species) == N_SPECIES
    for species, models in by_species.items():
        primaries = [m for m in models if m["is_primary_reference"]]
        assert len(primaries) == 1, species
        # The choice follows the stated rule, not iteration order.
        available = {m["isoform"] for m in models}
        expected = "IIIc" if "IIIc" in available else "IIIb"
        assert primaries[0]["isoform"] == expected, species


def test_isoform_models_are_roles_not_two_unlabelled_primaries(model_index):
    roles = {m["model_role"] for m in model_index["models"]}
    assert roles == {"validated_isoform_IIIb", "validated_isoform_IIIc"}
    # is_primary_reference says what a model is *used as*; model_role says what it
    # *is*. Collapsing them would lose one of the two facts.
    assert any(m["is_primary_reference"] and m["model_role"] == "validated_isoform_IIIb"
               for m in model_index["models"])


def test_a_generic_gene_still_has_one_primary_reference_per_species():
    """The stricter rule must not have changed how a normal gene behaves."""
    from exondomaincompare.shared_gene_analysis import model_roles

    runs = sorted((ROOT / "runs").glob("*/website_indices/generic/"
                                       "protein_coordinate_model.json"))
    if not runs:
        pytest.skip("no generic run available")
    for path in runs:
        index = json.loads(path.read_text(encoding="utf-8"))
        models = index.get("models") or []
        assert not model_roles.role_errors(models), path
        for m in models:
            assert m["model_role"] == "primary_reference"
            assert m["is_primary_reference"] is True
        species = [m["species_id"] for m in models]
        assert len(set(species)) == len(species), "one model per species as before"


def test_the_renderer_refuses_a_model_without_an_explicit_identity(tmp_path):
    """Identity may not be inferred from array order or file name."""
    index = json.loads(MODEL_INDEX.read_text(encoding="utf-8")) \
        if MODEL_INDEX.is_file() else None
    if index is None:
        pytest.skip("model index not built")
    stripped = dict(index)
    stripped["models"] = [{k: v for k, v in index["models"][0].items()
                           if k not in ("model_id", "model_role")}]
    path = tmp_path / "model.json"
    path.write_text(json.dumps(stripped), encoding="utf-8")
    proc = subprocess.run(
        ["node", "scripts/plotting/render_main_figures.mjs", str(path),
         str(tmp_path / "out")],
        cwd=ROOT, capture_output=True, text=True)
    assert proc.returncode != 0
    assert "model_id" in proc.stderr


def test_rendered_files_record_the_model_they_belong_to():
    summary = DERIVED / "figures" / "main" / "render_summary.json"
    rows = _json(summary)
    assert len(rows) == 287
    for row in rows:
        assert row["model_id"], row["stem"]
        assert row["model_role"], row["stem"]
        assert row["figure_type"], row["stem"]
    # A species with two models produces two distinct sets of files rather than one
    # overwriting the other.
    human = [r for r in rows if r["species_id"] == "homo_sapiens"]
    assert len({r["model_id"] for r in human}) == 2
    assert len({r["stem"] for r in human}) == len(human)


# --------------------------------------------------------------------------- #
# 3. the curated catalogue
# --------------------------------------------------------------------------- #
def test_the_gallery_does_not_show_one_card_per_rendered_file(catalogue):
    rendered = list((DERIVED / "figures" / "main").glob("*.svg"))
    if not rendered:
        pytest.skip("derived figures not rendered")
    visible_comparative = catalogue["counts"]["n_comparative_main"]
    per_species = max(catalogue["counts"]["species_scope_card_counts"])
    assert visible_comparative == N_COMPARATIVE_MAIN
    assert per_species <= 8, "a species scope is a curated set, not a file listing"
    # 287 derived files; no view of the Gallery comes close to that many cards.
    assert visible_comparative + per_species < 25


def test_the_comparative_scope_is_the_curated_main_set(catalogue):
    cards = catalogue["comparative_cards"]
    assert len(cards) == N_COMPARATIVE_MAIN
    assert all(c["kind"] == "main" for c in cards)
    assert all(c["scope"] == "comparative" for c in cards)
    assert all(not c["species_id"] for c in cards), "a comparative card has no species"
    # One card per canonical figure type; no figure type twice.
    types = [c["figure_type"] for c in cards]
    assert len(set(types)) == len(types)


def test_a_species_scope_uses_the_accepted_species_pipeline(catalogue):
    expected_types = {
        "primary_exon_projection", "integrated_domain_architecture",
        "boundary_on_architecture", "signed_boundary_distances",
        "boundary_class_summary", "validated_exon_domain_architecture",
    }
    for species, scope in catalogue["species_scopes"].items():
        assert scope["n_cards"] in SPECIES_SCOPE_CARD_COUNTS, species
        types = {c["figure_type"] for c in scope["cards"]}
        assert expected_types <= types, species
        assert all(c["species_id"] == species for c in scope["cards"])
        assert all(c["scope"] == "species" for c in scope["cards"])
        assert scope["categories"], species


def test_a_species_with_both_isoforms_gets_a_model_selector_not_two_cards(catalogue):
    scope = catalogue["species_scopes"]["homo_sapiens"]
    assert scope["n_models"] == 2
    for card in scope["cards"]:
        options = (card["model_selection"] or {}).get("options") or []
        assert [o["label"] for o in options] == ["IIIb", "IIIc"], card["figure_id"]
        assert sum(1 for o in options if o["is_default"]) == 1
        assert card["n_modes"] == 2, "both models are views of one card"
    # And no card id is duplicated per isoform.
    ids = [c["figure_id"] for c in scope["cards"]]
    assert len(set(ids)) == len(ids)
    assert not any("IIIb" in i or "IIIc" in i for i in ids)


def test_both_isoform_models_are_preserved_and_never_merged(catalogue, model_index):
    scope = catalogue["species_scopes"]["homo_sapiens"]
    proteins = {m["protein_id"] for m in scope["models"]}
    assert len(proteins) == 2, "two real proteins, not one synthetic average"
    comparison = next(c for c in scope["cards"]
                      if c["figure_type"] == "isoform_model_comparison")
    assert comparison["n_modes"] == 2
    assert "not merged" in comparison["interpretation"]
    # The two models keep their own coordinates.
    human = [m for m in model_index["models"] if m["species_id"] == "homo_sapiens"]
    assert len({m["protein_length"] for m in human}) == 2


def test_a_species_missing_a_model_shows_an_exact_unavailable_status(catalogue):
    for species, isoform in MISSING_COMBINATIONS:
        scope = catalogue["species_scopes"][species]
        assert scope["isoforms_unavailable"] == [isoform], species
        entry = next(m for m in scope["models"] if m["isoform"] == isoform)
        assert entry["availability_status"] == "no_architecture_model"
        assert entry["unavailable_reason"]
        assert not entry["model_id"], "no model is fabricated for it"
        # The species is not hidden, and the model it does have is shown.
        assert any(m["model_id"] for m in scope["models"])
        stated = [c for c in scope["cards"]
                  if any(u["isoform"] == isoform
                         for u in (c["model_selection"] or {}).get("unavailable") or [])]
        assert stated, f"{species}: the absence is never stated on a card"


def test_a_card_with_no_figure_says_why_rather_than_disappearing(catalogue):
    scope = catalogue["species_scopes"]["canis_lupus_familiaris"]
    empty = [c for c in scope["cards"] if c["n_modes"] == 0]
    assert empty, "this species has no exon series, so those cards have no figure"
    for card in empty:
        assert card["availability"]["status"] == "unavailable"
        assert card["availability"]["unavailable_models"]
        assert all(u["reason"] for u in card["availability"]["unavailable_models"])
        assert not card["thumbnail"], "no placeholder image stands in for the figure"


def test_no_card_is_registered_per_export_format_or_filter(catalogue):
    cards = (catalogue["comparative_cards"] + catalogue["supplements"]
             + [c for s in catalogue["species_scopes"].values() for c in s["cards"]])
    for card in cards:
        ids = card["figure_id"]
        assert not any(ids.endswith(f".{fmt}") or f"_{fmt}_" in ids
                       for fmt in ("png", "svg", "pdf", "tsv")), ids
        # Formats live on the card, several per card.
        assert set(card["export_formats"]) <= {"png", "svg", "pdf", "tsv"}
    architecture = next(c for c in catalogue["comparative_cards"]
                        if c["figure_type"] == "all_species_exon_domain_architecture")
    labels = {m["label"] for m in architecture["modes"]}
    assert {"IIIb only", "IIIc only", "Mammals", "Non-mammals"} <= labels
    assert architecture["n_modes"] == 5, "the four filters are modes of one card"


def test_known_duplicates_are_merged_into_the_card_they_duplicate(catalogue):
    retired = set(catalogue["retired_source_figures"])
    assert retired == MERGED_SOURCE_FIGURES
    cards = catalogue["comparative_cards"] + catalogue["supplements"]
    ids = {c["figure_id"] for c in cards}
    for name in MERGED_SOURCE_FIGURES:
        assert f"fgfr2_cmp_{name}" not in ids
        assert f"fgfr2_supp_{name}" not in ids


def test_the_byte_identical_pairs_really_are_byte_identical():
    """The merge is justified by the files, not by their names."""
    figures = FREEZE / "13_final_pre_interpro_closure" / "figures"
    if not figures.is_dir():
        pytest.skip("freeze figures unavailable")
    pairs = [

        ("Figure_9A_FGFR2_local_synteny_5neighbor_paper",
         "Figure_9_FGFR2_local_synteny_neighborhood"),
        ("Supplement_full_length_MSA_QC_histograms",
         "Supplement_full_length_MSA_outliers"),
    ]
    for kept, merged in pairs:
        for ext in ("png", "svg", "pdf"):
            a, b = figures / f"{kept}.{ext}", figures / f"{merged}.{ext}"
            if a.is_file() and b.is_file():
                assert a.read_bytes() == b.read_bytes(), f"{kept} vs {merged} ({ext})"


def test_supplements_are_a_separate_hidden_level(catalogue):
    supplements = catalogue["supplements"]
    assert len(supplements) == N_SUPPLEMENTS
    assert all(c["kind"] == "supplement" for c in supplements)
    assert catalogue["filters"]["supplements_hidden_by_default"] is True
    # And they are not mixed into the main comparative set.
    main_ids = {c["figure_id"] for c in catalogue["comparative_cards"]}
    assert main_ids.isdisjoint({c["figure_id"] for c in supplements})


def test_a_multi_species_dataset_opens_on_the_comparative_scope(catalogue):
    assert catalogue["multi_species"] is True
    assert catalogue["default_scope"] == "comparative"


def test_species_scopes_follow_the_canonical_order(catalogue, model_index):
    from exondomaincompare.shared_gene_analysis import species_order

    listed = [s["species_id"] for s in catalogue["filters"]["species"]]
    assert listed == list(catalogue["species_scopes"].keys())
    assert listed == species_order.order_species(listed)
    assert listed != sorted(listed), "canonical order is not alphabetical"


def test_there_is_no_flat_per_species_appendix(catalogue):
    """The 62 per-species architecture files are one card per species, not 62."""
    validated = [c for s in catalogue["species_scopes"].values() for c in s["cards"]
                 if c["figure_type"] == "validated_exon_domain_architecture"]
    assert len(validated) == N_SPECIES
    assert sum(c["n_modes"] for c in validated) == N_MODELS
    # And none of them sits in the comparative scope or the supplements.
    assert all(c["scope"] == "species" for c in validated)
    supp = {c["figure_id"] for c in catalogue["supplements"]}
    assert supp.isdisjoint({c["figure_id"] for c in validated})


def test_every_card_carries_its_scientific_context(catalogue):
    cards = (catalogue["comparative_cards"] + catalogue["supplements"]
             + [c for s in catalogue["species_scopes"].values() for c in s["cards"]])
    for card in cards:
        for field in ("figure_type", "scientific_question", "interpretation",
                      "renderer", "category", "kind", "scope"):
            assert card[field], f"{card['figure_id']}: {field}"
        assert card["source_data"], card["figure_id"]


def test_no_card_is_a_synthetic_placeholder(catalogue):
    cards = (catalogue["comparative_cards"] + catalogue["supplements"]
             + [c for s in catalogue["species_scopes"].values() for c in s["cards"]])
    for card in cards:
        if card["n_modes"]:
            assert card["thumbnail"], card["figure_id"]
            for mode in card["modes"]:
                for rel in mode["formats"].values():
                    assert (ROOT / rel).is_file(), rel
        else:
            # No figure, and it says why instead of showing something invented.
            assert card["availability"]["unavailable_models"], card["figure_id"]


# --------------------------------------------------------------------------- #
# 4. the gallery index the frontend reads
# --------------------------------------------------------------------------- #
def test_the_gallery_index_keeps_scope_and_species_on_every_card():
    index = _json(GALLERY_INDEX)
    assert index["default_scope"] == "comparative"
    assert index["supplements_hidden_by_default"] is True
    for card in index["figures"]:
        assert card["scope"] in ("comparative", "species")
        if card["scope"] == "species":
            assert card["species_id"], card["figure_id"]
        else:
            assert not card["species_id"], card["figure_id"]


def test_the_visible_card_count_per_scope_is_curated():
    index = _json(GALLERY_INDEX)
    figures = index["figures"]
    comparative = [c for c in figures
                   if c["scope"] == "comparative" and c["kind"] == "main"]
    assert len(comparative) == N_COMPARATIVE_MAIN
    for species in index["species"]:
        visible = [c for c in figures if c["scope"] == "species"
                   and c["species_id"] == species["species_id"]
                   and c["kind"] == "main"]
        assert len(visible) in SPECIES_SCOPE_CARD_COUNTS, species["species_id"]


def test_the_gallery_frontend_orders_the_new_categories():
    source = (ROOT / "webapp" / "frontend" / "src" / "pages"
              / "FigureGallery.jsx").read_text(encoding="utf-8")
    index = _json(GALLERY_INDEX)
    used = {c["category"] for c in index["figures"]}
    for category in used:
        assert f'"{category}"' in source, f"{category} has no declared reading order"


# --------------------------------------------------------------------------- #
# 5. the freeze
# --------------------------------------------------------------------------- #
def test_no_derived_output_is_written_into_the_freeze(catalogue):
    for path in [MODEL_INDEX, CATALOGUE, GALLERY_INDEX, INVENTORY]:
        assert FREEZE not in path.parents, path
        assert str(path).startswith(str(DERIVED)), path
    for card in catalogue["comparative_cards"] + catalogue["supplements"]:
        for mode in card["modes"]:
            for rel in mode["formats"].values():
                # Freeze figures are *read* by the catalogue and never rewritten.
                assert (ROOT / rel).is_file(), rel


def test_the_freeze_figure_index_is_untouched():
    """The catalogue is served from the derived overlay, not by editing the freeze."""
    frozen = FREEZE / "13_final_pre_interpro_closure" / "website_indices" \
        / "figure_index.json"
    if not frozen.is_file():
        pytest.skip("frozen figure index unavailable")
    index = json.loads(frozen.read_text(encoding="utf-8"))
    ours = _json(GALLERY_INDEX)
    assert index != ours, "the freeze must keep its own index"
    assert "provenance" not in index or index.get("provenance", {}).get(
        "generated_by") != "scripts/fgfr2/gallery_catalogue.py"
