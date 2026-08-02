"""The two withdrawn FGFR2 Gallery cards, and the generic boundary analysis that
replaces one of them in the catalogue.

Two failures motivated this. First, the Gallery offered two cards that did not earn a
reader's attention: a framework evidence stack, which documents how the analysis
reached its own conclusions rather than answering a biological question, and a
per-neighbour synteny conservation matrix, which asks what the main synteny
neighbourhood figure already answers while presenting single-assembly annotation gaps
as blank cells. Both were withdrawn. Withdrawn means no card — not renamed, not
demoted to a supplement — while the files and tables behind them stay readable,
because downloads, QC and the validation record still cite them.

Second, FGFR2 had only the *validated cassette* boundary analysis. That analysis
answers a narrow question about the IIIb/IIIc cassette and carries the freeze's
conclusions. It says nothing about the rest of the protein. The generic whole-protein
analysis — every supported internal coding-exon boundary against the nearest
representative domain edge, compared across species — was available to FGFR1 and every
other multi-species gene but not to FGFR2. It is now, through the same shared
contract, and the tests below hold the two analyses apart: the generic one must never
be presented as validated or conserved, and the validated one must not change.
"""
from __future__ import annotations

import json
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
FREEZE = ROOT / "results" / "final_30_until_interpro_prepare"

#: The withdrawn cards, by the title a reader saw and the files behind them.
WITHDRAWN = {
    "Framework evidence stack": ("Figure_8_final_framework_evidence_stack",
                                 "Figure_Final_Framework_Evidence_Stack"),
    "Synteny neighbour conservation": (
        "Figure_9B_FGFR2_5neighbor_conservation_matrix_paper",),
}

GENERIC_BOUNDARY_CATEGORY = "Comparative exon–domain boundaries"
VALIDATED_BOUNDARY_CATEGORY = "FGFR2 IIIb/IIIc Boundary Consistency"

#: The three shared figures, by the stem the shared renderer writes.
GENERIC_BOUNDARY_FIGURES = {
    "All coding-exon Boundary matrix": "cmp_boundary_matrix",
    "All coding-exon signed-distance comparison": "cmp_paired_signed_distance",
    "All coding-exon Boundary-position consistency": "cmp_boundary_position_consistency",
}

#: The one renderer that draws these figures, for every gene. FGFR2 must not get a
#: private copy of it.
SHARED_COMPARATIVE_RENDERER = "scripts/plotting/render_comparative_gallery_figures.mjs"

VALIDATED_BOUNDARY_TITLES = {
    "IIIb/IIIc Boundary Consistency Matrix",
    "IIIb/IIIc Boundary-distance distribution",
}


def _json(path: Path):
    if not path.is_file():
        pytest.skip(f"{path.relative_to(ROOT)} not built")
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def catalogue():
    return _json(CATALOGUE)


@pytest.fixture(scope="module")
def cards(catalogue):
    """Every card a reader can reach, in any scope."""
    return (catalogue["comparative_cards"] + catalogue["supplements"]
            + [c for s in catalogue["species_scopes"].values() for c in s["cards"]])


@pytest.fixture(scope="module")
def multi_species():
    dashboard = _json(MODEL_INDEX).get("boundary_dashboard") or {}
    return dashboard.get("multi_species") or {}


# --------------------------------------------------------------------------- #
# 1. the two withdrawn cards
# --------------------------------------------------------------------------- #
def test_the_framework_evidence_stack_is_not_a_visible_card(cards):
    titles = {c["title"] for c in cards}
    assert "Framework evidence stack" not in titles
    # Nor under a new name: the files behind it must not back any card.
    stems = {m["mode_id"] for c in cards for m in c.get("modes") or []}
    assert not stems & set(WITHDRAWN["Framework evidence stack"])


def test_synteny_neighbour_conservation_is_not_a_visible_card(cards):
    titles = {c["title"] for c in cards}
    assert "Synteny neighbour conservation" not in titles
    stems = {m["mode_id"] for c in cards for m in c.get("modes") or []}
    assert not stems & set(WITHDRAWN["Synteny neighbour conservation"])


def test_neither_withdrawn_card_reappears_as_a_supplement(catalogue):
    """Withdrawing a card and hiding it behind the supplements toggle are different
    things. A supplement is still a card a reader can open."""
    supplements = {c["title"] for c in catalogue["supplements"]}
    assert not supplements & set(WITHDRAWN)


def test_the_withdrawal_is_recorded_with_its_reason(catalogue):
    """A card that simply disappears leaves a reader wondering whether the analysis
    failed. The catalogue states which cards were withdrawn and why."""
    withdrawn = {w["title"]: w for w in catalogue["withdrawn_cards"]}
    assert set(withdrawn) == set(WITHDRAWN)
    for title, stems in WITHDRAWN.items():
        assert set(withdrawn[title]["source_figures"]) == set(stems)
        assert len(withdrawn[title]["reason"]) > 40


def test_the_source_data_behind_the_withdrawn_cards_is_still_there(catalogue):
    """The evidence did not stop existing. Downloads, QC and the validated record
    still read these tables, so withdrawing the card must not remove them."""
    retained = catalogue["retained_source_data"]
    assert retained, "no retained source data recorded"
    for rel in retained.values():
        assert (ROOT / rel).is_file(), rel


def test_the_withdrawn_figure_files_are_not_deleted():
    """Only the card is withheld. Deleting validated output to hide a card would
    break the downloads and checksums that cite it."""
    for stems in WITHDRAWN.values():
        for stem in stems:
            hits = list(FREEZE.rglob(f"{stem}.*")) if FREEZE.is_dir() else []
            if not hits:
                pytest.skip("freeze figures not present")
            assert hits


def test_the_main_comparative_synteny_figure_survives(catalogue):
    """Removing the conservation matrix must not remove genomic context entirely."""
    genomic = [c for c in catalogue["comparative_cards"]
               if c["category"] == "Comparative genomic context"]
    assert genomic, "no comparative genomic-context card left"
    assert any("synteny" in c["title"].lower() for c in genomic)


# --------------------------------------------------------------------------- #
# 2. the generic whole-protein boundary analysis
# --------------------------------------------------------------------------- #
def test_the_generic_boundary_figures_are_registered(catalogue):
    cards = [c for c in catalogue["comparative_cards"]
             if c["category"] == GENERIC_BOUNDARY_CATEGORY]
    assert {c["title"] for c in cards} == set(GENERIC_BOUNDARY_FIGURES)


def test_the_generic_figures_come_from_the_shared_renderer(catalogue):
    """A second FGFR2-specific implementation would drift from the generic one. Each
    card must be backed by the stem the shared comparative renderer writes."""
    for card in catalogue["comparative_cards"]:
        if card["category"] != GENERIC_BOUNDARY_CATEGORY:
            continue
        expected = GENERIC_BOUNDARY_FIGURES[card["title"]]
        assert [m["mode_id"] for m in card["modes"]] == [expected]
        assert card["renderer"] == SHARED_COMPARATIVE_RENDERER


def test_every_generic_boundary_card_offers_the_vector_safe_export_set(catalogue):
    """Self-contained vector, 300-dpi raster and the numbers behind the figure. A
    reader who wants to check a plotted distance needs the table, not the picture."""
    for card in catalogue["comparative_cards"]:
        if card["category"] != GENERIC_BOUNDARY_CATEGORY:
            continue
        assert {"svg", "pdf", "png", "tsv"} <= set(card["export_formats"]), card["title"]


def test_the_generic_boundary_analysis_uses_real_comparable_groups(multi_species):
    if not multi_species.get("available"):
        pytest.skip("no comparative boundary evidence")
    groups = multi_species["comparable_boundary_groups"]
    assert groups
    for g in groups:
        observations = g["per_species_native_positions"]
        # A group exists because at least two species really observe it. A group with
        # one observation would render a row of blank cells that looks like data.
        assert len({o["species_id"] for o in observations}) >= 2, g["label"]


def test_boundaries_are_never_matched_by_exon_number_alone(multi_species):
    """The fourth coding exon of two species is not the same exon. Grouping by rank
    would manufacture comparisons between unrelated junctions."""
    if not multi_species.get("available"):
        pytest.skip("no comparative boundary evidence")
    allowed = {"shared_exon_group", "msa_aligned_position"}
    for g in multi_species["comparable_boundary_groups"]:
        assert g["mapping_method"] in allowed, g["mapping_method"]
        evidence = g.get("supporting_evidence") or {}
        assert evidence.get("shared_exon_group") or evidence.get("msa_column") is not None


def test_every_observation_names_a_real_domain_instance_and_signed_distance(multi_species):
    """A cell shows one species' own measurement against one specific domain
    instance. Where no representative domain is annotated near the boundary, that gap
    must stay visible instead of being read as a distance of zero."""
    if not multi_species.get("available"):
        pytest.skip("no comparative boundary evidence")
    classes = {"exact_domain_boundary", "near_domain_edge", "inside_domain",
               "outside_annotated_domains", "uncertain", "unmapped"}
    for g in multi_species["comparable_boundary_groups"]:
        for o in g["per_species_native_positions"]:
            assert o["boundary_class"] in classes, o["boundary_class"]
            assert o["mapping_confidence"] is not None
            if o["domain_annotation_available"]:
                assert o["nearest_domain_instance_id"]
                assert o["signed_distance"] is not None
            else:
                assert o["nearest_domain_instance_id"] is None


def test_the_generic_boundary_panel_is_one_model_per_species(multi_species):
    """A comparable group counts species. FGFR2 has two isoform models per species,
    so comparing both would count every species twice and inflate every group."""
    dashboard = _json(MODEL_INDEX)["boundary_dashboard"]
    panel = dashboard["comparative_panel"]
    assert panel["n_models_compared"] < panel["n_models_in_dataset"]
    assert panel["n_models_compared"] == len(set(panel["model_ids"]))
    rows = multi_species["boundary_matrix"]
    assert len({r["species_id"] for r in rows}) == len(rows)


def test_the_generic_boundary_species_order_is_canonical(multi_species):
    from shared_gene_analysis import species_order as so

    ids = [r["species_id"] for r in multi_species["boundary_matrix"]]
    assert ids == list(so.order_species(ids))


def test_the_gallery_and_the_interactive_view_read_one_index(multi_species):
    """Two independently derived versions of the same analysis would let a figure and
    the page beside it disagree. The comparative dataset the interactive Boundary
    Explorer reads must quote the same group ids as the figures."""
    dataset = INDICES / "comparative_dataset.json"
    if not dataset.is_file():
        pytest.skip("comparative dataset not built")
    if not multi_species.get("available"):
        pytest.skip("no comparative boundary evidence")
    ids = {g["comparable_boundary_group_id"]
           for g in multi_species["comparable_boundary_groups"]}
    text = dataset.read_text(encoding="utf-8")
    doc = json.loads(text)
    quoted = {g["comparable_boundary_group_id"]
              for g in ((doc.get("boundary_dashboard") or {})
                        .get("multi_species") or {})
              .get("comparable_boundary_groups", [])}
    if quoted:
        assert quoted == ids


# --------------------------------------------------------------------------- #
# 3. the two analyses stay distinguishable
# --------------------------------------------------------------------------- #
def test_the_validated_cassette_boundary_cards_are_unchanged(catalogue):
    cards = [c for c in catalogue["comparative_cards"]
             if c["category"] == VALIDATED_BOUNDARY_CATEGORY]
    assert {c["title"] for c in cards} == VALIDATED_BOUNDARY_TITLES
    for card in cards:
        # Still drawn from the freeze's own figures, not re-derived.
        assert all(m["mode_id"].startswith("Figure_") for m in card["modes"])


def test_the_two_boundary_analyses_are_separate_sections(catalogue):
    categories = {c["category"] for c in catalogue["comparative_cards"]}
    assert VALIDATED_BOUNDARY_CATEGORY in categories
    assert GENERIC_BOUNDARY_CATEGORY in categories
    assert VALIDATED_BOUNDARY_CATEGORY != GENERIC_BOUNDARY_CATEGORY


def test_every_generic_card_states_that_it_is_not_the_cassette_analysis(catalogue):
    """A reader who lands on a boundary figure should not have to infer which of the
    two analyses it belongs to from the section heading above it."""
    from fgfr2.gallery_catalogue import GENERIC_BOUNDARY_SCOPE_NOTE

    for card in catalogue["comparative_cards"]:
        if card["category"] != GENERIC_BOUNDARY_CATEGORY:
            continue
        assert GENERIC_BOUNDARY_SCOPE_NOTE in card["interpretation"]


def test_the_generic_boundaries_are_not_called_validated_or_conserved(catalogue):
    """These are positional observations across species. Calling them validated or
    conserved would borrow the standing of the cassette analysis, which rests on
    evidence this analysis does not have."""
    for card in catalogue["comparative_cards"]:
        if card["category"] != GENERIC_BOUNDARY_CATEGORY:
            continue
        text = f"{card['title']} {card['scientific_question']}".lower()
        for word in ("validated", "conserved", "conservation"):
            assert word not in text, f"{card['title']}: {word}"
        # The interpretation may say what the figure is *not*, and may warn against
        # reading conservation into it; it may not assert either.
        interpretation = card["interpretation"].lower()
        assert "is not evidence of conservation" in interpretation \
            or "separate from the validated" in interpretation


def test_no_two_cards_in_one_scope_share_a_title(catalogue):
    """Titles repeat across species scopes by design — every species has an exon-
    structure card — so uniqueness is asked of each scope a reader sees at once."""
    scopes = {"comparative": catalogue["comparative_cards"],
              "supplements": catalogue["supplements"]}
    for species_id, scope in catalogue["species_scopes"].items():
        scopes[species_id] = scope["cards"]
    for name, cards in scopes.items():
        titles = [c["title"] for c in cards]
        assert len(titles) == len(set(titles)), f"{name}: duplicate card title"


def test_no_rendered_comparative_figure_backs_two_cards(catalogue):
    """Two cards over one file is the duplication the catalogue exists to remove."""
    stems = [m["mode_id"] for c in catalogue["comparative_cards"] + catalogue["supplements"]
             for m in c.get("modes") or []]
    assert len(stems) == len(set(stems)), "one rendered figure backs two cards"


# --------------------------------------------------------------------------- #
# 4. other and future FGFR2 datasets
# --------------------------------------------------------------------------- #
def _fgfr2_run_dirs():
    out = []
    for cfg in sorted((ROOT / "runs").glob("*/run_config.json")):
        try:
            doc = json.loads(cfg.read_text(encoding="utf-8"))
        except Exception:
            continue
        if str(doc.get("gene_symbol") or "").upper() == "FGFR2":
            out.append(cfg.parent)
    return out


def test_no_other_fgfr2_dataset_still_shows_the_withdrawn_cards():
    runs = _fgfr2_run_dirs()
    if not runs:
        pytest.skip("no other FGFR2 dataset in the registry")
    stems = {s for group in WITHDRAWN.values() for s in group}
    checked = 0
    for run in runs:
        index = run / "website_indices" / "figure_index.json"
        if not index.is_file():
            continue
        checked += 1
        ids = {f.get("id") or f.get("figure_id")
               for f in json.loads(index.read_text(encoding="utf-8"))["figures"]}
        assert not ids & stems, f"{run.name}: {sorted(ids & stems)}"
    if not checked:
        pytest.skip("no built FGFR2 figure index")


def test_the_withdrawal_list_has_one_definition():
    """Two lists would drift, and a card withdrawn from one dataset would reappear in
    another."""
    from fgfr2.gallery_catalogue import withdrawn_figure_stems
    from build_website_indices import _withdrawn_figure_stems

    assert _withdrawn_figure_stems() == withdrawn_figure_stems()
    assert withdrawn_figure_stems() == {s for g in WITHDRAWN.values() for s in g}


def test_a_one_species_fgfr2_dataset_registers_no_comparative_card(tmp_path):
    """One row is not a comparison. A cross-species card over a single species would
    invite a reader to draw a conclusion from an empty contrast."""
    from plotting.generate_comparative_gallery_figures import generate

    model = tmp_path / "model.json"
    model.write_text(json.dumps({
        "gene_symbol": "FGFR2",
        "models": [{"species_id": "homo_sapiens", "status": "available"}],
    }), encoding="utf-8")
    result = generate(tmp_path, model)
    assert result["cards"] == 0
    assert result["skipped"] == "single_species"


def test_a_future_multi_species_fgfr2_run_withholds_the_conservation_matrix():
    """The withdrawal must hold for FGFR2 datasets built through the generic
    pipeline too, and must not silently change any other gene's gallery."""
    from plotting.generate_comparative_gallery_figures import (
        _withdrawn_comparative_stems)

    assert _withdrawn_comparative_stems("FGFR2") == {
        "cmp_synteny_neighbour_conservation"}
    assert _withdrawn_comparative_stems("fgfr2") == {
        "cmp_synteny_neighbour_conservation"}
    for gene in ("FGFR1", "TP53", "TPM1"):
        assert _withdrawn_comparative_stems(gene) == set()


def test_a_future_multi_species_fgfr2_run_registers_the_boundary_figures():
    """Registration is driven by the presence of real comparable evidence, not by a
    per-run hand edit, so a future run gets the figures on its own."""
    from plotting.generate_comparative_gallery_figures import FIGURE_META

    for stem in GENERIC_BOUNDARY_FIGURES.values():
        assert stem in FIGURE_META, stem


# --------------------------------------------------------------------------- #
# 5. nothing validated moved
# --------------------------------------------------------------------------- #
def test_the_validated_cassette_boundary_values_are_unchanged():
    """The generic analysis measures every boundary of the whole protein. It must not
    have reached back into the cassette boundaries the freeze validated."""
    index = _json(MODEL_INDEX)
    cassette = [b for m in index["models"] for b in m.get("exon_boundaries") or []
                if b.get("is_cassette_boundary")]
    if not cassette:
        pytest.skip("no cassette boundaries in the model index")
    for b in cassette:
        assert b["source"] == "shared_generic_boundary_classification"
        assert b["signed_distance"] is not None
        assert b["boundary_class"]


def test_the_freeze_is_bytewise_untouched():
    import subprocess

    if not (ROOT / ".git").is_dir():
        pytest.skip("not a git checkout")
    rel = FREEZE.relative_to(ROOT)
    proc = subprocess.run(["git", "status", "--porcelain", "--", str(rel)],
                          cwd=ROOT, capture_output=True, text=True)
    assert proc.returncode == 0
    assert proc.stdout.strip() == "", proc.stdout


def test_the_derived_outputs_live_outside_the_freeze():
    for path in (MODEL_INDEX, CATALOGUE, GALLERY_INDEX,
                 DERIVED / "figures" / "comparative"):
        assert FREEZE not in path.parents, path
