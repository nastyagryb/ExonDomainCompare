from __future__ import annotations

import json
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FRONTEND = PROJECT_ROOT / "webapp" / "frontend" / "src"

FGFR1_RUN_ID = "2026-07-23_1100_fgfr1_gallus_core_pilot"
TP53_RUN_ID = "2026-07-21_1436_custom_run"

SVG_NS = "{http://www.w3.org/2000/svg}"


def run_dir(run_id: str) -> Path:
    return PROJECT_ROOT / "runs" / run_id


def src(rel: str) -> str:
    return (FRONTEND / rel).read_text(encoding="utf-8")


_BLOCK_COMMENT = re.compile(r"/\*.*?\*/", re.DOTALL)
_LINE_COMMENT = re.compile(r"^\s*//.*$", re.MULTILINE)


def rendered(rel: str) -> str:
    return _LINE_COMMENT.sub("", _BLOCK_COMMENT.sub("", src(rel)))


def figures_index(run_id: str) -> dict:
    fp = run_dir(run_id) / "website_indices" / "figures_index.json"
    if not fp.exists():
        pytest.skip(f"no figures index for {run_id}")
    return json.loads(fp.read_text())


def cards(run_id: str) -> list[dict]:
    return figures_index(run_id).get("figures") or []


# --------------------------------------------------------------------------- #
# Part 1 — Candidate Evidence page layout
# --------------------------------------------------------------------------- #

def test_candidate_ranking_is_a_compact_strip_not_a_half_page_column():
    text = rendered("pages/GeneExplorer.jsx")
    css = (FRONTEND / "App.css").read_text(encoding="utf-8")

    assert "CandidateStrip" in text, "compact candidate strip is missing"
    assert "cand-strip" in text and "cand-chipcard" in text

    # The old two-column master-detail grid must be gone from both layers.
    assert "cmd-layout" not in text, "the master-detail grid is still rendered"
    assert ".cmd-layout" not in css, "the master-detail grid CSS still defines a half-page column"
    assert "grid-template-columns: minmax(340px, 1fr) minmax(360px, 1.3fr)" not in css

    # The strip is a horizontally scrolling row, not a grid column.
    strip = re.search(r"\.cand-strip \{[^}]*\}", css)
    assert strip, ".cand-strip is not styled"
    assert "overflow-x: auto" in strip.group(0)
    assert "display: flex" in strip.group(0)


def test_compact_candidate_item_shows_only_ranking_essentials():
    text = rendered("pages/GeneExplorer.jsx")
    strip = text.split("function CandidateStrip", 1)[1].split("\nfunction ", 1)[0]
    assert "rank_label" in strip
    assert "aa {c.aa_start}–{c.aa_end}" in strip
    assert "c.length" in strip
    assert "overall_score" in strip
    assert "affected" in strip
    # Full evidence detail belongs to the selected candidate, not the strip.
    for detail in ("score_components", "confidence_reason", "exon_evidence"):
        assert detail not in strip, f"{detail} should only appear after selection"


def test_selected_candidate_header_states_score_strength_and_validation():
    text = rendered("pages/GeneExplorer.jsx")
    head = text.split("cand-sel-head", 1)[1][:1400]
    assert "Evidence score:" in head
    assert "Evidence strength:" in head
    assert "Biological validation: not validated" in head
    assert "csh-rank" in head, "the candidate id is not given a larger header treatment"


def test_workspace_offers_four_full_width_analysis_views():
    text = rendered("pages/GeneExplorer.jsx")
    section = text.split("const CAND_SECTIONS", 1)[1].split("];", 1)[0]
    for label in ("Candidate overview", "Full Isoform Alignment",
                  "Candidate-focused alignment", "Exon & domain context"):
        assert label in section, f"missing workspace view: {label}"
    # The full alignment is a top-level view, not a subsection of isoform evidence.
    assert '"Isoform evidence"' not in section


def test_full_alignment_uses_the_whole_content_width():
    css = (FRONTEND / "App.css").read_text(encoding="utf-8")
    # A bare 1fr track keeps min-width:auto and lets the wide alignment grid push
    # the page sideways instead of scrolling inside its own container.
    assert "grid-template-columns: 320px minmax(0, 1fr)" in css
    assert ".iso-msa .msa-scroll { overflow-x: auto; }" in css
    assert re.search(r"\.iso-msa \{[^}]*min-width: 0", css)
    assert re.search(r"\.cand-alignment\.full-width[^{]*\{[^}]*width: 100%", css)


# --------------------------------------------------------------------------- #
# Part 2 — full within-species Isoform Alignment
# --------------------------------------------------------------------------- #

def test_alignment_offers_all_five_modes_and_does_not_default_to_the_candidate():
    text = rendered("pages/viewers/MsaExplorer.jsx")
    modes = text.split("const ALN_MODES", 1)[1].split("];", 1)[0]
    for label in ("Full alignment", "Differences to primary", "Variable regions",
                  "Conserved regions", "Candidate-focused region"):
        assert label in modes, f"missing alignment mode: {label}"
    # "candidate" may only be the initial mode of the dedicated focused view.
    assert 'useState(focusCandidate ? "candidate" : "full")' in text


def test_complete_alignment_is_not_cropped_to_the_candidate_interval():
    text = rendered("pages/viewers/MsaExplorer.jsx")
    # The auto-scroll effect must be gated on the focused view.
    effect = re.search(r"useEffect\(\(\) => \{\s*if \(focusCandidate && rawBand\)", text)
    assert effect, "the candidate auto-scroll is not restricted to the focused view"
    # Only the candidate mode narrows the visible column range.
    assert "const inCandidateMode = mode === \"candidate\" && rawBand" in text
    assert "const viewStart = inCandidateMode" in text


def test_alignment_navigation_covers_blocks_termini_and_zoom():
    text = rendered("pages/viewers/MsaExplorer.jsx")
    for control in ("var ▶", "◀ var", "gap ▶", "◀ gap", "N-term", "C-term", "Reset view"):
        assert control in text, f"missing navigation control: {control}"
    assert "ZOOM_STEPS" in text, "horizontal zoom is missing"
    assert "AlignmentMinimap" in text, "the full-alignment minimap is missing"
    assert "stepBlock" in text and "blocksOf" in text


def test_block_navigation_actually_advances_the_viewport():
    text = rendered("pages/viewers/MsaExplorer.jsx")
    assert "const LEAD = 8" in text
    assert "goTo = (col) => setOffset(clampOffset(col - LEAD))" in text
    assert "const here = offset + LEAD" in text
    assert "Math.floor(windowSize / 4)" not in text


def test_sequence_selection_offers_presets_and_defaults_to_all_isoforms():
    text = rendered("pages/viewers/MsaExplorer.jsx")
    for preset in ("Select all", "Primary + differing", "Curated only", "Reset selection"):
        assert preset in text, f"missing sequence preset: {preset}"
    # Nothing is hidden until the user hides it.
    assert "useState(() => new Set())" in text


def test_每_row_reports_the_full_protein_model_identity():
    text = rendered("pages/viewers/MsaExplorer.jsx")
    block = text.split("aln-name-meta", 1)[1][:900]
    assert "transcript_id" in block
    assert "protein_length" in block
    assert "curation_status" in block
    assert "primary" in block and "alternative" in block
    assert "aln-ident" in block, "identity to primary is not shown per row"


def test_alignment_annotations_use_generic_definitions():
    text = rendered("pages/viewers/MsaExplorer.jsx")
    for label in ("Variable columns", "Discriminating columns", "Gap boundaries",
                  "Conserved columns"):
        assert label in text
    # Generic vocabulary only — no FGFR2 cassette terminology in the generic view.
    iso_view = text.split("function IsoformMsaView", 1)[1].split("\nfunction ", 1)[0]
    for banned in ("IIIb", "IIIc", "cassette"):
        assert banned not in iso_view, f"FGFR2-specific term leaked into the generic view: {banned}"


def test_candidate_overlay_can_be_switched_off():
    text = rendered("pages/viewers/MsaExplorer.jsx")
    assert "showCandidate" in text
    assert "const band = showCandidate ? rawBand : null" in text
    assert "Highlight {selectedCandidate?.rank_label" in text


# --------------------------------------------------------------------------- #
# Part 3/4 — candidate-focused view and pairwise comparisons
# --------------------------------------------------------------------------- #

def test_candidate_focused_alignment_is_separate_from_the_full_alignment():
    text = rendered("pages/GeneExplorer.jsx")
    assert 'sect === "alignment"' in text
    assert 'sect === "focused"' in text
    # The focused view passes the flag; the full view must not.
    focused = text.split('sect === "focused"', 1)[1][:1200]
    assert "focusCandidate" in focused
    full = text.split('sect === "alignment"', 1)[1].split('sect === "focused"', 1)[0]
    assert "focusCandidate" not in full


def test_pairwise_comparisons_are_collapsed_with_a_compact_summary():
    text = rendered("pages/GeneExplorer.jsx")
    block = text.split("function PairwiseComparisonDetails", 1)[1].split("\nfunction ", 1)[0]
    assert "<details" in block, "pairwise comparisons are not collapsible"
    assert "open" not in re.search(r"<details[^>]*>", block).group(0), \
        "the pairwise section must be collapsed by default"
    assert "Pairwise comparison details" in block
    # Compact pre-expansion summary.
    assert "pairwise comparison" in block
    assert "supporting" in block and "exonAligned" in block


# --------------------------------------------------------------------------- #
# Part 5 — alignment exports (one menu, one figure source per format)
# --------------------------------------------------------------------------- #

def test_alignment_export_menu_offers_one_figure_in_three_formats_plus_data():
    text = rendered("pages/GeneExplorer.jsx")
    block = text.split("function AlignmentExportMenu", 1)[1].split("\nfunction ", 1)[0]
    assert "— SVG" in block and "— PDF" in block and "— PNG" in block
    assert "Alignment FASTA" in block
    assert "Alignment summary TSV" in block
    assert "Candidate-region TSV" in block
    # All three formats come from the same figure specification.
    assert block.count("buildFigure()") >= 3, "formats must share one figure source"
    # The residue-level alignment is offered as one document, not as loose pages.
    assert "multi-page PDF" in block
    # No technical extras in the normal menu.
    for banned in ("coordinate JSON", "Advanced data", "provenance"):
        assert banned not in block


# --------------------------------------------------------------------------- #
# Part 6 — multi-species behaviour
# --------------------------------------------------------------------------- #

def test_multi_species_keeps_a_species_selector_and_per_species_alignment():
    text = rendered("pages/GeneExplorer.jsx")
    panel = text.split("function CandidateEvidencePanel", 1)[1].split("\nfunction ", 1)[0]
    assert "speciesOptions" in panel and "setActiveSpecies" in panel
    # Within-species alignment is resolved per species, not once per dataset.
    assert "withinSpeciesAlignment(model, activeSpecies)" in text
    assert "Isoform alignment unavailable: only one protein model is available." in text


def test_cross_species_msa_stays_a_separate_top_level_view():
    text = rendered("pages/GeneExplorer.jsx")
    # A within-species alignment must never be presented as conservation.
    msa = rendered("pages/viewers/MsaExplorer.jsx")
    assert "Within-species protein isoform alignment — not cross-species" in msa
    assert "Cross-species primary-protein alignment — one primary isoform per species." not in msa
    assert "cross_species_msa" in msa
    assert '"msa"' in text


# --------------------------------------------------------------------------- #
# Parts 9 + 17 — figure sources render as self-contained vector output
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def rendered_alignment_figures(tmp_path_factory) -> dict:
    if shutil.which("node") is None:
        pytest.skip("node is not available")
    index = run_dir(FGFR1_RUN_ID) / "website_indices" / "isoform_alignment_index.json"
    if not index.exists():
        pytest.skip("no isoform alignment index for the reference run")
    out = tmp_path_factory.mktemp("alnfig")
    proc = subprocess.run(
        ["node", str(PROJECT_ROOT / "scripts" / "plotting" / "render_alignment_figures.mjs"),
         str(index), str(out)],
        capture_output=True, text=True, cwd=PROJECT_ROOT, check=False,
    )
    assert proc.returncode == 0, f"figure rendering failed: {proc.stderr}"
    return {"dir": out, "summary": json.loads((out / "summary.json").read_text())}


def test_alignment_figures_use_the_real_reference_dataset(rendered_alignment_figures):
    s = rendered_alignment_figures["summary"]
    assert s["gene"] == "FGFR1"
    assert s["species"] == "Gallus gallus"
    assert s["primary"] == "NP_990841.2"
    assert s["n_rows"] == 8, "all eight FGFR1 protein models must be drawn"
    assert s["n_columns"] == 823, "the complete alignment must be represented"
    # The candidate maps onto the alignment, not onto raw residue indices.
    assert s["candidate_columns"] == [31, 118]


def test_both_alignment_figures_are_produced(rendered_alignment_figures):
    out = rendered_alignment_figures["dir"]
    assert (out / "full_isoform_alignment.svg").exists()
    assert (out / "candidate_alignment_detail.svg").exists()


@pytest.mark.parametrize("name", ["full_isoform_alignment.svg", "candidate_alignment_detail.svg"])
def test_exported_svg_is_valid_and_self_contained(rendered_alignment_figures, name):
    path = rendered_alignment_figures["dir"] / name
    text = path.read_text()

    root = ET.parse(path).getroot()          # must parse as real XML
    assert root.tag == f"{SVG_NS}svg"
    # An intrinsic size is required by Preview and by PDF converters.
    assert root.get("width") and root.get("height")
    assert root.get("viewBox")

    # Nothing may depend on the application stylesheet or the browser DOM.
    assert "var(--" not in text, "CSS custom properties do not resolve in a standalone SVG"
    assert "foreignObject" not in text
    assert "class=" not in text
    # Explicit typography and paint.
    assert "font-family=" in text and "font-size=" in text
    assert "fill=" in text and "stroke=" in text


def test_full_alignment_figure_is_not_uniform_bars(rendered_alignment_figures):
    text = (rendered_alignment_figures["dir"] / "full_isoform_alignment.svg").read_text()
    assert text.count("<rect") > 100, "too few marks to show per-column structure"
    for label in ("alignment column", "identity", "Variable columns",
                  "residue differing from primary", "gap in this isoform",
                  "Major variable blocks"):
        assert label in text, f"missing legend/axis element: {label}"
    assert "823 alignment columns" in text
    assert "Exploratory analysis" in text, "the exploratory status must be stated"


def test_candidate_detail_figure_carries_both_coordinate_systems(rendered_alignment_figures):
    text = (rendered_alignment_figures["dir"] / "candidate_alignment_detail.svg").read_text()
    assert "alignment column" in text
    assert "aa" in text
    assert "affected" in text and "unaffected" in text
    assert "reference" in text, "the primary protein must be labelled as the reference"
    assert "not a biologically validated splicing event" in text


def test_interactive_svg_components_carry_literal_colours():
    offenders = []
    for path in sorted(FRONTEND.rglob("*.jsx")):
        text = path.read_text(encoding="utf-8")
        if "<svg" not in text:
            continue
        if "var(--" in text:
            offenders.append(str(path.relative_to(FRONTEND)))
    assert not offenders, f"CSS variables inside SVG components: {offenders}"


def test_png_export_targets_publication_resolution():
    text = rendered("pages/viewers/plotExport.js")
    assert "scale = 4" in text, "PNG export must rasterise at publication resolution"
    assert "serializeSvg(svgEl, { width: W, height: H })" in text
    # The serialiser must produce a standalone document.
    ser = text.split("export function serializeSvg", 1)[1].split("\nexport ", 1)[0]
    assert 'setAttribute("width"' in ser and 'setAttribute("height"' in ser
    assert 'removeAttribute("class")' in ser
    assert "#ffffff" in ser, "a white paper background is required for raster output"


# --------------------------------------------------------------------------- #
# Part 7 — Figure Gallery information architecture
# --------------------------------------------------------------------------- #

CATEGORIES = ("Exon structure", "Isoform analysis", "Domain architecture",
              "Exon–domain boundaries", "Genomic context", "Exploratory candidates",
              "Supplements")


def test_gallery_declares_the_seven_categories_in_reading_order():
    text = rendered("pages/FigureGallery.jsx")
    block = text.split("const CATEGORY_ORDER", 1)[1].split("];", 1)[0]
    # "Supplements" is referenced through a constant, so resolve it first.
    supplements = re.search(r'const SUPPLEMENTS = "([^"]+)"', text)
    assert supplements and supplements.group(1) == "Supplements"
    resolved = block.replace("SUPPLEMENTS", f'"{supplements.group(1)}"')
    for cat in CATEGORIES:
        assert cat in resolved, f"missing gallery category: {cat}"
    positions = [resolved.index(c) for c in CATEGORIES]
    assert positions == sorted(positions), "categories are not in reading order"


def test_gallery_card_shows_identity_and_stage_but_not_provenance_text():
    text = rendered("pages/FigureGallery.jsx")
    card = text.split("function FigCard", 1)[1].split("\nfunction ", 1)[0]
    assert "scientific_question" in card
    assert "gene" in card and "species" in card and "protein_id" in card
    assert "STAGE_LABEL" in card, "the analysis stage is not shown"
    for fmt in ("SVG", "PDF", "PNG", "TSV"):
        assert f">{fmt}</a>" in card, f"missing {fmt} download on the card"
    # Reproducibility provenance moved to the detail view.
    assert "feature_sources" not in card, "provenance text is still on the card face"
    assert "Sources:" not in card


def test_one_card_per_figure_with_several_formats():
    text = rendered("pages/FigureGallery.jsx")
    assert "const seen = new Set()" in text
    assert "if (seen.has(key)) continue" in text
    # Formats are attributes of one card, not separate cards.
    assert "formats: {" in text


def test_validated_fgfr2_gallery_keeps_its_curated_set_and_group_order():
    text = rendered("pages/FigureGallery.jsx")
    assert 'f.kind !== "main" || f.category === SUPPLEMENTS' in text, \
        "the per-figure supplement kind is no longer honoured"
    # The index's own category order still reaches the rendered order; the
    # declared reading order only picks the categories that are actually present.
    ordering = text.split("const orderedCategories", 1)[1][:700]
    assert "index.categories" in ordering, \
        "the index's own category order is no longer honoured"
    assert "preferred.filter((c) => present.has(c))" in ordering
    # The FGFR2 "Supplement" group is folded into the shared Supplements category.
    assert 'group === "Supplement" ? SUPPLEMENTS : group' in text

    fp = (PROJECT_ROOT / "results" / "final_30_until_interpro_prepare"
          / "13_final_pre_interpro_closure" / "website_indices" / "figure_index.json")
    if not fp.exists():
        pytest.skip("FGFR2 example figure index not present")
    doc = json.loads(fp.read_text())
    kinds = {f.get("kind") for f in doc["figures"]}
    assert "supplement" in kinds and "main" in kinds, \
        "the FGFR2 index no longer distinguishes main figures from supplements"
    assert doc["groups"][0] == "Framework", "the curated FGFR2 group order changed"


def test_generic_indices_may_also_mark_a_card_as_a_supplement():
    text = rendered("pages/FigureGallery.jsx")
    assert 'kind: f.kind === "supplement" ? "supplement" : "main"' in text, \
        "the generic normaliser drops the per-card supplement marking"
    assert 'kind: "main",' not in text, "a card's kind is still hardcoded"


@pytest.mark.parametrize("run_id", [FGFR1_RUN_ID, TP53_RUN_ID])
def test_the_signature_supplement_is_marked_as_a_supplement(run_id):
    supplements = [c for c in cards(run_id)
                   if "member_signature" in (c.get("figure_id") or "")]
    assert supplements, "the member-database signature card is missing"
    for card in supplements:
        assert card.get("kind") == "supplement", \
            f"{card['figure_id']} would be shown as a main figure"


def test_gallery_offers_a_scope_selector_only_for_multi_species_datasets():
    text = rendered("pages/FigureGallery.jsx")
    assert "multiSpecies" in text
    # Comparative is the default Scope; a Scope remembered in the URL or the
    # session wins over it, which is what keeps a linked species view stable.
    initial = text.split("function readInitialScope", 1)[1].split("\n}", 1)[0]
    assert 'if (!multiSpecies) return "comparative";' in initial
    assert initial.rstrip().endswith('return "comparative";')
    assert "useState(() => readInitialScope(multiSpecies))" in text
    assert "{multiSpecies && (" in text, "the scope selector must be conditional"
    assert '<option value="comparative">Comparative</option>' in text


def test_captions_are_editable_and_downloadable_but_not_permanent_page_furniture():
    text = rendered("pages/FigureGallery.jsx")
    assert "function autoCaption" in text
    editor = text.split("function CaptionEditor", 1)[1].split("\nconst STAGE_LABEL", 1)[0]
    assert "<textarea" in editor, "the caption must be editable"
    assert ".caption.txt" in editor, "a downloadable caption file is missing"
    # Captions live in the detail drawer, not on the card face.
    card = text.split("function FigCard", 1)[1].split("\nfunction ", 1)[0]
    assert "CaptionEditor" not in card
    # The caption states the analysis stage and the exploratory status.
    auto = text.split("function autoCaption", 1)[1].split("\nfunction ", 1)[0]
    assert "analysis stage" in auto
    assert "Exploratory analysis unless explicitly marked as validated." in auto


# --------------------------------------------------------------------------- #
# Parts 16 + 19 — the final single-species figure set
# --------------------------------------------------------------------------- #

REMOVED_CARDS = (
    "generic_domain_architecture",
    "generic_exon_domain_boundary_distribution",
    "domain_exon_projection",
    "domain_candidate_overlay",
    "domain_boundary_overlay",
    "boundary_evidence_supplement",
    "selected_boundary_detail",
)


@pytest.mark.parametrize("run_id", [FGFR1_RUN_ID, TP53_RUN_ID])
def test_removed_and_redundant_figure_cards_are_gone(run_id):
    ids = [c.get("figure_id", "") for c in cards(run_id)]
    offenders = [i for i in ids if any(bad in i for bad in REMOVED_CARDS)]
    assert not offenders, f"redundant cards still registered: {offenders}"


@pytest.mark.parametrize("run_id", [FGFR1_RUN_ID, TP53_RUN_ID])
def test_every_card_declares_one_canonical_category(run_id):
    missing = [c.get("figure_id") for c in cards(run_id)
               if c.get("category") not in CATEGORIES]
    assert not missing, f"cards without a canonical category: {missing}"


@pytest.mark.parametrize("run_id", [FGFR1_RUN_ID, TP53_RUN_ID])
def test_no_duplicate_figure_cards(run_id):
    ids = [c.get("figure_id") for c in cards(run_id)]
    dupes = {i for i in ids if ids.count(i) > 1}
    assert not dupes, f"duplicate figure cards: {dupes}"


@pytest.mark.parametrize("run_id", [FGFR1_RUN_ID, TP53_RUN_ID])
def test_each_card_offers_all_three_download_formats(run_id):
    thin = [c.get("figure_id") for c in cards(run_id)
            if c.get("status") == "available"
            and not (c.get("svg_url") and c.get("pdf_url") and c.get("png_url"))]
    assert not thin, f"cards missing a download format: {thin}"


@pytest.mark.parametrize("run_id", [FGFR1_RUN_ID, TP53_RUN_ID])
def test_every_card_carries_a_caption(run_id):
    missing = [c.get("figure_id") for c in cards(run_id)
               if c.get("status") == "available" and not c.get("caption")]
    assert not missing, f"cards without a caption: {missing}"


# The final single-species figure set: one card per analysis, named by the analysis
# it answers rather than by the file it happens to produce. Each entry is matched
# against the card ids, so a renamed or dropped figure fails here by name.
REQUIRED_ANALYSES = {
    "exon-to-protein projection": "primary_exon_projection",
    "transcript and protein structure": "transcript_exon_structure",
    "transcript-model comparison": "transcript_model_comparison",
    "differences-only comparison": "transcript_model_comparison_differences",
    "isoform alignment overview": "full_isoform_alignment",
    "wrapped residue-level alignment": "wrapped_alignment",
    "candidate alignment detail": "candidate_alignment_detail",
    "integrated domain architecture": "integrated_domain_architecture",
    "member-database signatures": "member_signature_supplement",
    "boundaries on the architecture": "boundary_on_architecture",
    "signed boundary distances": "signed_boundary_distances",
    "boundary-class summary": "boundary_class_summary",
    "local genomic neighbourhood": "neighbourhood",
    "candidate ranking": "candidate_ranking",
    "candidate exon/domain context": "candidate_domain_context",
}


def test_single_species_figure_set_covers_the_required_analyses():
    entries = [c for c in cards(FGFR1_RUN_ID) if c.get("status") == "available"]
    by_cat: dict[str, list[str]] = {}
    for c in entries:
        by_cat.setdefault(c.get("category", "?"), []).append(c.get("figure_id", ""))

    for cat in ("Exon structure", "Isoform analysis", "Domain architecture",
                "Exon–domain boundaries", "Genomic context", "Exploratory candidates"):
        assert by_cat.get(cat), f"no figure in category {cat}"

    ids = [c.get("figure_id", "") for c in entries]
    for analysis, token in REQUIRED_ANALYSES.items():
        assert any(token in i for i in ids), \
            f"the final set answers no figure for: {analysis}"

    # One card per analysis: a second card matching the same token would mean the
    # same question is answered twice, or a format was registered as its own card.
    for analysis, token in REQUIRED_ANALYSES.items():
        hits = [i for i in ids if token in i]
        # "transcript_model_comparison" is a prefix of its differences-only sibling.
        if token == "transcript_model_comparison":
            hits = [i for i in hits if not i.endswith("_differences")]
        assert len(hits) == 1, f"{analysis} is covered by {len(hits)} cards: {hits}"

    main = [c for c in entries if c.get("category") != "Supplements"]
    assert len(main) <= 16, f"too many main gallery cards: {len(main)}"


def test_candidate_figures_use_cautious_wording():
    entries = cards(FGFR1_RUN_ID)
    cand = [c for c in entries if c.get("category") == "Exploratory candidates"]
    assert cand, "no candidate figures registered"
    for c in cand:
        blob = " ".join(str(c.get(k, "")) for k in
                        ("title", "scientific_question", "interpretation", "caption")).lower()
        assert "exploratory" in blob or "not biologically validated" in blob or "not validated" in blob, \
            f"{c.get('figure_id')} does not mark the candidate as exploratory"
        assert "validated event" not in blob.replace("not validated event", "")


def test_no_generated_index_leaks_personal_absolute_paths():
    for run_id in (FGFR1_RUN_ID, TP53_RUN_ID):
        idx = run_dir(run_id) / "website_indices"
        if not idx.exists():
            continue
        for fp in idx.rglob("*.json"):
            assert "/Users/" not in fp.read_text(), f"personal path leaked into {fp}"


# --------------------------------------------------------------------------- #
# Gene Explorer keeps only its compact export menu
# --------------------------------------------------------------------------- #

def test_gene_explorer_modules_do_not_embed_generated_figure_lists():
    for rel in ("pages/viewers/ExonMap.jsx", "pages/viewers/ProteinArchitecture.jsx",
                "pages/viewers/BoundaryExplorer.jsx"):
        text = rendered(rel)
        assert "publication_figures" not in text, f"{rel} still embeds a static figure list"
        assert "figures_index" not in text
