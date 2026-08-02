"""Export validation for the single-species publication figures.

These tests inspect the actual exported artefacts — a standalone SVG, a vector
PDF and a 300 dpi PNG — rather than the source code that produces them, because
the defect being guarded against was invisible in the source: the figures looked
correct on screen and exported as a single full-page JPEG.

The reference dataset is the real post-cluster FGFR1 / Gallus gallus run.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pdf_probe import probe_pdf, probe_png, probe_svg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FGFR1_RUN = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"
TP53_RUN = ROOT / "runs" / "2026-07-21_1436_custom_run"
# The production renderer, deliberately: these tests validate the artefacts that
# actually ship, produced by the same code the Gene Explorer exports through.
RENDERER = ROOT / "scripts" / "plotting" / "render_main_figures.mjs"

# The figures this phase is responsible for, in Gallery reading order.
MAIN_FIGURES = [
    "primary_exon_projection",
    "integrated_domain_architecture",
    "boundary_on_architecture",
    "signed_boundary_distances",
    "boundary_class_summary",
]

# Available on demand (e.g. from the Gene Explorer with a boundary selected), but
# deliberately not registered as permanent Gallery cards.
ON_DEMAND_FIGURES = ["selected_boundary_detail", "selected_signed_distance"]

# Publication page geometry: a figure narrower than a single column or wider than
# a landscape page is a layout defect. The broken pipeline produced a 55-inch page.
MIN_PAGE_IN = 2.0
MAX_PAGE_IN = 14.0


def _coordinate_model(run_dir: Path) -> Path:
    return run_dir / "website_indices" / "generic" / "protein_coordinate_model.json"


def _render(run_dir: Path, out_dir: Path, selected: str = "E4 → E5") -> Path:
    """Render the main figures and normalise the per-species stems for the tests."""
    if shutil.which("node") is None:
        pytest.skip("node is required to render the figure specifications")
    model = _coordinate_model(run_dir)
    if not model.exists():
        pytest.skip(f"coordinate model missing: {model}")
    result = subprocess.run(
        ["node", str(RENDERER), str(model), str(out_dir),
         # exercise the on-demand selected-boundary export as well
         f"--selected-boundary={selected}"],
        capture_output=True, text=True, cwd=ROOT,
    )
    if result.returncode != 0:
        pytest.fail(f"figure renderer failed:\n{result.stdout}\n{result.stderr}")
    # The renderer namespaces its output by species; the assertions below address
    # figures by kind, so link the primary species' figures to the bare kind.
    doc = json.loads(model.read_text())
    primary = next((m for m in doc["models"] if m.get("role") == "primary"), doc["models"][0])
    sp = primary["species_id"]
    for kind in MAIN_FIGURES + ON_DEMAND_FIGURES:
        for ext in ("svg", "pdf", "tsv"):
            for prefix in ("main", "ondemand"):
                src = out_dir / f"{prefix}_{sp}_{kind}.{ext}"
                if src.exists():
                    (out_dir / f"{kind}.{ext}").write_bytes(src.read_bytes())
    return out_dir


def _context(run_dir: Path) -> dict:
    """The real dataset facts the figures are expected to depict."""
    doc = json.loads(_coordinate_model(run_dir).read_text())
    primary = next((m for m in doc["models"] if m.get("role") == "primary"), doc["models"][0])
    boundaries = [b for b in primary.get("exon_boundaries") or []
                  if b.get("signed_distance") is not None]
    class_counts: dict[str, int] = {}
    for b in boundaries:
        key = b.get("classification") or b.get("category") or "unknown"
        class_counts[key] = class_counts.get(key, 0) + 1
    return {
        "gene": primary.get("gene_symbol") or doc.get("gene_symbol"),
        "species": primary.get("species_id"),
        "proteinId": primary.get("protein_id"),
        "transcriptId": primary.get("transcript_id"),
        "proteinLength": primary.get("protein_length"),
        # `models` is one entry per species; the protein models of that species are
        # its transcript-derived isoform models.
        "n_models": primary.get("n_transcript_models") or len(doc["models"]),
        "n_coding_exons": len(primary.get("exons") or []),
        "n_boundaries": len(boundaries),
        "near_edge_threshold_aa": primary.get("near_edge_threshold_aa"),
        "class_counts": class_counts,
        "domain_instances": [
            {"domain_instance_id": d.get("domain_instance_id"),
             "full_label": d.get("full_label"),
             "start": d.get("start"), "end": d.get("end"),
             "instance_number": d.get("instance_number")}
            for d in primary.get("representative_domains") or []
        ],
        "tm": primary.get("tm_regions") or [],
    }


@pytest.fixture(scope="module")
def figures(tmp_path_factory) -> Path:
    if not FGFR1_RUN.exists():
        pytest.skip(f"reference run missing: {FGFR1_RUN}")
    return _render(FGFR1_RUN, tmp_path_factory.mktemp("fgfr1_figures"))


@pytest.fixture(scope="module")
def context() -> dict:
    if not FGFR1_RUN.exists():
        pytest.skip(f"reference run missing: {FGFR1_RUN}")
    return _context(FGFR1_RUN)


# --------------------------------------------------------------------------- #
# The reference dataset really is the dataset we claim to plot
# --------------------------------------------------------------------------- #

def test_reference_dataset_matches_the_expected_real_results(context):
    assert context["gene"] == "FGFR1"
    assert context["species"] == "gallus_gallus"
    assert context["proteinId"] == "NP_990841.2"
    assert context["transcriptId"] == "NM_205510.2"
    assert context["proteinLength"] == 817
    assert context["n_models"] == 8
    assert context["n_coding_exons"] == 17
    assert context["n_boundaries"] == 16
    assert context["near_edge_threshold_aa"] == 5
    counts = context["class_counts"]
    # 6 near a domain edge, 8 inside a domain, 2 outside any annotated domain
    assert sum(counts.values()) == 16
    assert counts.get("near_domain_edge", counts.get("near_edge")) == 6
    assert counts.get("inside_domain") == 8
    assert counts.get("outside_annotated_domains", counts.get("outside_domain")) == 2


def test_the_four_domain_instances_are_distinct_and_correctly_placed(context):
    instances = context["domain_instances"]
    assert len(instances) == 4
    # Every instance is uniquely addressable, so three Ig-like domains sharing one
    # InterPro accession can never collapse onto each other.
    ids = [d["domain_instance_id"] for d in instances]
    assert len(set(ids)) == 4
    assert ids == [
        "IPR007110:33-118",
        "IPR007110:145-244",
        "IPR007110:253-355",
        "IPR001245:476-750",
    ]
    labels = [d["full_label"] for d in instances]
    assert labels == [
        "Ig-like domain 1 · aa 33–118",
        "Ig-like domain 2 · aa 145–244",
        "Ig-like domain 3 · aa 253–355",
        "Ser-Thr/Tyr kinase domain · aa 476–750",
    ]


def test_one_real_pytmhmm_region_is_present(context):
    assert len(context["tm"]) == 1
    assert context["tm"][0]["source"] == "pytmhmm"


# --------------------------------------------------------------------------- #
# PDF: true vector output at a sensible physical size
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_pdf_opens_and_is_a_single_well_formed_page(figures, name):
    info = probe_pdf(figures / f"{name}.pdf")
    assert info.ok, info.error
    assert info.header.startswith("%PDF-")
    assert info.n_pages == 1


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_pdf_contains_no_raster_image(figures, name):
    """The defining defect: the page was one embedded JPEG."""
    info = probe_pdf(figures / f"{name}.pdf")
    assert info.images == [], f"{name}.pdf embeds raster data: {info.images}"
    assert not info.is_single_raster_page


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_pdf_has_selectable_text_in_an_embedded_font(figures, name):
    info = probe_pdf(figures / f"{name}.pdf")
    assert info.has_vector_text, f"{name}.pdf carries no text objects"
    assert info.fonts, f"{name}.pdf references no font"
    assert info.n_text_ops >= 5


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_pdf_page_has_publication_dimensions(figures, name):
    info = probe_pdf(figures / f"{name}.pdf")
    assert MIN_PAGE_IN <= info.width_in <= MAX_PAGE_IN, \
        f"{name}.pdf is {info.width_in:.1f}in wide"
    assert MIN_PAGE_IN * 0.5 <= info.height_in <= MAX_PAGE_IN, \
        f"{name}.pdf is {info.height_in:.1f}in tall"
    # A figure taller than it is wide, or absurdly letterboxed, indicates the page
    # was sized from pixels rather than from the layout preset.
    assert 0.4 <= info.aspect <= 12.0, f"{name}.pdf aspect ratio is {info.aspect:.1f}"


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_pdf_draws_vector_geometry_and_an_explicit_white_background(figures, name):
    info = probe_pdf(figures / f"{name}.pdf")
    assert info.n_path_ops >= 5, f"{name}.pdf has almost no vector geometry"
    assert info.n_fill_ops >= 5
    # The paper rectangle is the first fill and is painted pure white.
    raw = (figures / f"{name}.pdf").read_bytes().decode("latin-1")
    assert "1 1 1 rg" in raw, f"{name}.pdf has no explicit white background"


def test_exon_pdf_carries_readable_exon_labels(figures):
    info = probe_pdf(figures / "primary_exon_projection.pdf")
    for label in ("E1", "E5", "E10", "E17"):
        assert info.contains(label), f"exon label {label} missing from the PDF text layer"
    assert info.contains("Primary protein position")
    assert info.contains("817")


def test_domain_pdf_names_every_track_separately(figures):
    info = probe_pdf(figures / "integrated_domain_architecture.pdf")
    for track in ("Representative domains", "Family / superfamily",
                  "Membrane topology", "Coding exons", "Exon boundaries",
                  "Candidate regions"):
        assert info.contains(track), f"track label {track!r} missing"
    assert info.contains("TM helix")


def test_boundary_pdf_carries_a_class_legend(figures):
    info = probe_pdf(figures / "boundary_on_architecture.pdf")
    assert info.contains("Near domain edge")
    assert info.contains("Inside domain")
    assert info.contains("Outside annotated domains")
    assert info.contains("threshold")


def test_signed_distance_pdf_shows_signed_values_and_the_zero_reference(figures):
    info = probe_pdf(figures / "signed_boundary_distances.pdf")
    assert info.contains("0 = domain edge")
    assert info.contains("E1 -> E2") or info.contains("E1")
    # The sign is the whole point of replacing the absolute-distance histogram.
    assert "-2" in info.text
    assert "+1" in info.text
    assert info.contains("start edge") and info.contains("end edge")


def test_no_exported_pdf_is_a_full_page_seventy_two_dpi_raster(figures):
    """Regression guard for the exact failure mode that was replaced."""
    for name in MAIN_FIGURES:
        info = probe_pdf(figures / f"{name}.pdf")
        assert not info.has_raster_image
        assert info.n_text_ops > 0
        assert info.width_pt < 2000, \
            f"{name}.pdf page is {info.width_pt:.0f}pt wide, i.e. sized from pixels"


# --------------------------------------------------------------------------- #
# SVG: standalone, explicitly painted, no stylesheet dependency
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_svg_is_valid_xml_with_explicit_dimensions(figures, name):
    svg = probe_svg(figures / f"{name}.svg")
    assert svg["width"] and svg["height"] and svg["viewBox"]
    assert svg["n_marks"] > 5
    assert svg["n_text"] > 3


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_svg_carries_no_stylesheet_dependency(figures, name):
    """Every mark must paint itself; a CSS-only fill renders black when exported."""
    svg = probe_svg(figures / f"{name}.svg")
    assert not svg["has_css_var"], f"{name}.svg contains an unresolved var(--…)"
    assert not svg["has_foreign_object"]
    assert not svg["has_external_css"]
    assert not svg["has_class_attr"], f"{name}.svg relies on CSS classes"
    assert not svg["has_style_element"]
    assert svg["unpainted"] == [], \
        f"{name}.svg has {len(svg['unpainted'])} unpainted marks: {set(svg['unpainted'])}"


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_svg_text_declares_its_own_font(figures, name):
    raw = (figures / f"{name}.svg").read_text(encoding="utf-8")
    root = ET.fromstring(raw)
    texts = [el for el in root.iter() if el.tag.split("}")[-1] == "text"]
    assert texts
    for el in texts:
        assert el.get("font-family"), "text element without an explicit font-family"
        assert el.get("font-size"), "text element without an explicit font-size"
        assert el.get("fill"), "text element without an explicit fill"


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_svg_has_no_black_fallback_blocks(figures, name):
    """A black fill on a feature block is the signature of the old export bug."""
    raw = (figures / f"{name}.svg").read_text(encoding="utf-8")
    root = ET.fromstring(raw)
    for el in root.iter():
        if el.tag.split("}")[-1] != "rect":
            continue
        fill = (el.get("fill") or "").lower()
        area = float(el.get("width") or 0) * float(el.get("height") or 0)
        # Small black marks are legitimate (ticks, rules); a black *block* is not.
        if area > 60:
            assert fill not in ("#000", "#000000", "black"), \
                f"{name}.svg has a black feature block of area {area:.0f}"


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_svg_paints_an_explicit_white_background(figures, name):
    raw = (figures / f"{name}.svg").read_text(encoding="utf-8")
    root = ET.fromstring(raw)
    first_rect = next(el for el in root.iter() if el.tag.split("}")[-1] == "rect")
    assert (first_rect.get("fill") or "").lower() in ("#ffffff", "#fff", "white")


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_all_svg_coordinates_lie_inside_the_viewbox(figures, name):
    raw = (figures / f"{name}.svg").read_text(encoding="utf-8")
    root = ET.fromstring(raw)
    _, _, vw, vh = (float(v) for v in root.get("viewBox").split())
    tol = 1.0
    for el in root.iter():
        tag = el.tag.split("}")[-1]
        if tag == "rect":
            x, y = float(el.get("x")), float(el.get("y"))
            w, h = float(el.get("width")), float(el.get("height"))
            assert -tol <= x and x + w <= vw + tol, f"{name}.svg rect overflows in x"
            assert -tol <= y and y + h <= vh + tol, f"{name}.svg rect overflows in y"
        elif tag == "line":
            for a, lim in ((("x1", "x2"), vw), (("y1", "y2"), vh)):
                for attr in a:
                    v = float(el.get(attr))
                    assert -tol <= v <= lim + tol, f"{name}.svg line overflows ({attr}={v})"
        elif tag == "text":
            x, y = float(el.get("x")), float(el.get("y"))
            assert -tol <= x <= vw + tol, f"{name}.svg text overflows in x"
            assert -tol <= y <= vh + tol, f"{name}.svg text overflows in y"


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_no_figure_text_falls_below_the_publication_minimum(figures, name):
    raw = (figures / f"{name}.svg").read_text(encoding="utf-8")
    root = ET.fromstring(raw)
    sizes = [float(el.get("font-size")) for el in root.iter()
             if el.tag.split("}")[-1] == "text"]
    assert sizes
    assert min(sizes) >= 5.5, f"{name}.svg has {min(sizes)}pt text"


def test_the_renderer_reports_no_layout_warnings(figures):
    summary = json.loads((figures / "render_summary.json").read_text())
    offenders = {f["stem"]: f["warnings"] for f in summary if f["warnings"]}
    assert not offenders, f"layout warnings: {offenders}"


# --------------------------------------------------------------------------- #
# Semantic content: the figure answers its scientific question
# --------------------------------------------------------------------------- #

def _svg_text(path: Path) -> str:
    root = ET.fromstring(path.read_text(encoding="utf-8"))
    return "\n".join(el.text or "" for el in root.iter()
                     if el.tag.split("}")[-1] == "text")


def test_exon_map_carries_no_domain_track(figures):
    """The exon map answers the exon-to-protein question only."""
    text = _svg_text(figures / "primary_exon_projection.svg")
    assert "Coding exons" in text
    for forbidden in ("Ig-like", "kinase", "Representative domains", "TM helix"):
        assert forbidden not in text, f"exon map should not show domains ({forbidden})"


def test_exon_map_states_gene_species_ids_and_termini(figures):
    text = _svg_text(figures / "primary_exon_projection.svg")
    assert "FGFR1" in text
    assert "Gallus gallus" in text
    assert "NP_990841.2" in text
    assert "NM_205510.2" in text
    assert "817 aa" in text
    assert "\n1\n" in f"\n{text}\n"  # explicit start residue
    assert "817" in text


def test_exon_map_labels_every_coding_exon(figures):
    text = _svg_text(figures / "primary_exon_projection.svg")
    labels = {line.strip() for line in text.splitlines()}
    for i in range(1, 18):
        assert f"E{i}" in labels, f"exon label E{i} is missing"


def test_exon_map_labels_the_candidate_outside_the_blocks(figures):
    text = _svg_text(figures / "primary_exon_projection.svg")
    assert "C1 · aa 31–118" in text
    assert "not validated" in text


def test_domain_figure_names_each_instance_in_the_legend(figures):
    text = _svg_text(figures / "integrated_domain_architecture.svg")
    for label in ("Ig-like domain 1 · aa 33–118",
                  "Ig-like domain 2 · aa 145–244",
                  "Ig-like domain 3 · aa 253–355",
                  "Ser-Thr/Tyr kinase domain · aa 476–750"):
        assert label in text, f"legend entry missing: {label}"


def test_domain_figure_keeps_family_and_topology_as_separate_rows(figures):
    text = _svg_text(figures / "integrated_domain_architecture.svg")
    assert "Family / superfamily" in text
    assert "Membrane topology" in text
    assert "TM helix · aa 372–394" in text
    assert "Family / homologous superfamily" in text


def test_domain_figure_uses_cautious_candidate_wording(figures):
    text = _svg_text(figures / "integrated_domain_architecture.svg")
    assert "biological validation: not validated" in text


def test_boundary_figure_keeps_interpretation_out_of_the_architecture(figures):
    """Distances and identifiers belong to the subtitle, not the plot area."""
    root = ET.fromstring((figures / "selected_boundary_detail.svg")
                         .read_text(encoding="utf-8"))
    texts = [(float(el.get("y")), el.text or "") for el in root.iter()
             if el.tag.split("}")[-1] == "text"]
    selected = [(y, t) for y, t in texts if t.startswith("Selected:")]
    assert len(selected) == 1, "expected exactly one selected-boundary annotation"
    header_bottom = selected[0][0]
    # Nothing below the header block may restate the signed distance.
    for y, t in texts:
        if y > header_bottom:
            assert "signed distance" not in t.lower(), \
                f"interpretation text inside the plot: {t!r}"


def test_boundary_figure_annotation_names_the_domain_instance(figures):
    text = _svg_text(figures / "selected_boundary_detail.svg")
    assert "Selected:" in text
    assert "signed distance" in text
    assert "→" in text  # E4 → E5 style transition label


def test_signed_distance_figure_lists_every_boundary_transition(figures):
    text = _svg_text(figures / "signed_boundary_distances.svg")
    for i in range(1, 17):
        assert f"E{i} → E{i + 1}" in text, f"missing row E{i} → E{i + 1}"


def test_signed_distance_figure_shows_the_near_edge_band_and_edge_symbols(figures):
    text = _svg_text(figures / "signed_boundary_distances.svg")
    assert "±5 aa" in text
    assert "0 = domain edge" in text
    assert "domain start edge" in text
    assert "domain end edge" in text


def test_signed_distance_figure_reproduces_the_real_signed_values(figures):
    """The seven documented FGFR1 boundaries, read back out of the figure."""
    text = _svg_text(figures / "signed_boundary_distances.svg")
    for value in ("-2", "+1", "+2", "-39", "+3", "-45", "+4"):
        assert value in text.splitlines() or value in text, \
            f"signed distance {value} not shown"


def test_boundary_class_summary_reports_the_real_distribution(figures):
    text = _svg_text(figures / "boundary_class_summary.svg")
    assert "16 internal coding-exon boundaries" in text
    assert "6 (38%)" in text
    assert "8 (50%)" in text
    assert "2 (13%)" in text


def test_absolute_distance_histogram_is_not_among_the_main_figures(figures):
    """The signed plot replaced it; keeping both would be redundant."""
    names = {p.stem for p in figures.glob("*.svg")}
    assert not any("histogram" in n for n in names)
    assert "signed_boundary_distances" in names


# --------------------------------------------------------------------------- #
# PNG: rasterised from the same specification at publication resolution
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def pngs(figures) -> Path:
    """Rasterise the exported SVGs at 300 dpi, as the download path does."""
    if shutil.which("rsvg-convert") is None:
        pytest.skip("rsvg-convert is required to rasterise the SVG exports")
    out = figures / "png"
    out.mkdir(exist_ok=True)
    for name in MAIN_FIGURES:
        subprocess.run(
            ["rsvg-convert", "--dpi-x", "300", "--dpi-y", "300", "--zoom", "4.1667",
             str(figures / f"{name}.svg"), "-o", str(out / f"{name}.png")],
            check=True, capture_output=True,
        )
    return out


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_png_has_the_expected_size_and_aspect_ratio(pngs, figures, name):
    png = probe_png(pngs / f"{name}.png")
    svg = probe_svg(figures / f"{name}.svg")
    expected = float(svg["width"]) / float(svg["height"])
    assert abs(png["aspect"] - expected) < 0.05, \
        f"{name}.png aspect {png['aspect']:.2f} != figure aspect {expected:.2f}"
    # 300 dpi over a physical width in points is width/72*300 pixels.
    expected_px = float(svg["width"]) / 72 * 300
    assert abs(png["width"] - expected_px) / expected_px < 0.02, \
        f"{name}.png is {png['width']}px wide, expected ~{expected_px:.0f}px at 300 dpi"


@pytest.mark.parametrize("name", MAIN_FIGURES)
def test_png_is_not_blank_and_has_readable_contrast(pngs, name):
    png = probe_png(pngs / f"{name}.png")
    assert not png["is_blank"]
    assert png["contrast"] > 0.5, f"{name}.png contrast is only {png['contrast']:.2f}"
    # A mostly dark image would mean the black-fallback bug returned.
    assert png["frac_light"] > 0.5, f"{name}.png is only {png['frac_light']:.0%} light"
    assert png["frac_dark"] < 0.25, f"{name}.png is {png['frac_dark']:.0%} dark"


# --------------------------------------------------------------------------- #
# Source tables: every figure ships its numbers
# --------------------------------------------------------------------------- #

def test_each_figure_group_ships_a_source_table(figures):
    for name in ("primary_exon_projection", "integrated_domain_architecture",
                 "boundary_on_architecture"):
        path = figures / f"{name}.tsv"
        assert path.exists(), f"missing source table {name}.tsv"
        lines = path.read_text().strip().splitlines()
        assert len(lines) > 1
        assert "\t" in lines[0]


def test_exon_source_table_lists_all_seventeen_coding_exons(figures):
    lines = (figures / "primary_exon_projection.tsv").read_text().strip().splitlines()
    assert len(lines) == 18  # header + 17 exons
    assert lines[0].split("\t")[:4] == [
        "exon_label", "exon_number", "protein_start_aa", "protein_end_aa"]


def test_boundary_source_table_carries_the_domain_instance_column(figures):
    header = (figures / "boundary_on_architecture.tsv").read_text().splitlines()[0].split("\t")
    for column in ("transition", "nearest_domain_instance_id", "nearest_domain_start",
                   "nearest_domain_end", "signed_distance_aa", "boundary_class"):
        assert column in header, f"boundary table lacks {column}"


# --------------------------------------------------------------------------- #
# TP53 regression: the pipeline is not FGFR1-specific
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def tp53_figures(tmp_path_factory) -> Path:
    if not TP53_RUN.exists():
        pytest.skip(f"TP53 regression run missing: {TP53_RUN}")
    return _render(TP53_RUN, tmp_path_factory.mktemp("tp53_figures"))


def test_tp53_renders_the_same_figure_set_as_vector_output(tp53_figures):
    for name in MAIN_FIGURES:
        info = probe_pdf(tp53_figures / f"{name}.pdf")
        assert info.ok, f"{name}.pdf: {info.error}"
        assert info.images == []
        assert info.has_vector_text
        assert MIN_PAGE_IN <= info.width_in <= MAX_PAGE_IN


def test_tp53_svgs_are_standalone_and_fully_painted(tp53_figures):
    for name in MAIN_FIGURES:
        svg = probe_svg(tp53_figures / f"{name}.svg")
        assert not svg["has_css_var"]
        assert not svg["has_class_attr"]
        assert svg["unpainted"] == []


def test_tp53_figures_name_their_own_gene_and_species(tp53_figures):
    context = _context(TP53_RUN)
    text = _svg_text(tp53_figures / "primary_exon_projection.svg")
    assert context["gene"] in text
    assert context["proteinId"] in text
