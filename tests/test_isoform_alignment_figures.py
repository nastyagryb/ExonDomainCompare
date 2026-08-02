"""Export validation for the within-species isoform alignment figures.

Three figures are validated as artefacts rather than as source code, because the
defects being guarded against are invisible in the source: an alignment figure
can look plausible on screen and still export as a rasterised page, drop the
residue information it claims to show, or degenerate into one flat coverage bar
per protein model.

  * full_isoform_alignment      — the complete alignment at column resolution
  * wrapped_alignment           — the same alignment at residue resolution,
                                  wrapped Jalview-style over several PDF pages
  * candidate_alignment_detail  — one candidate interval at residue resolution

Reference dataset: the real FGFR1 / Gallus gallus post-cluster run, with the
TP53 / Danio rerio run as a second, structurally different regression case.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pdf_probe import probe_pdf, probe_png, probe_svg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
# The production renderer, so the tests validate the artefacts the Gallery ships.
RENDERER = ROOT / "scripts" / "plotting" / "render_alignment_figures.mjs"
FGFR1_RUN = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"
TP53_RUN = ROOT / "runs" / "2026-07-21_1436_custom_run"

# Publication page geometry. A figure narrower than a single journal column or
# wider than a landscape page was sized from pixels, not from a layout preset.
MIN_PAGE_IN = 2.0
MAX_PAGE_IN = 14.0
# Smallest type any of these figures may set, in points.
MIN_FONT_PT = 5.5

SINGLE_PAGE_FIGURES = ["full_isoform_alignment", "candidate_alignment_detail"]


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #

def _render(run_dir: Path, out_dir: Path) -> dict:
    if shutil.which("node") is None:
        pytest.skip("node is required to render the figure specifications")
    index = run_dir / "website_indices" / "isoform_alignment_index.json"
    if not index.exists():
        pytest.skip(f"no isoform alignment index: {index}")
    proc = subprocess.run(
        ["node", str(RENDERER), str(index), str(out_dir)],
        capture_output=True, text=True, cwd=ROOT,
    )
    if proc.returncode != 0:
        pytest.fail(f"alignment figure rendering failed:\n{proc.stdout}\n{proc.stderr}")
    summary = json.loads((out_dir / "summary.json").read_text())
    return {"dir": out_dir, "summary": summary}


@pytest.fixture(scope="module")
def fgfr1(tmp_path_factory) -> dict:
    if not FGFR1_RUN.exists():
        pytest.skip(f"reference run missing: {FGFR1_RUN}")
    return _render(FGFR1_RUN, tmp_path_factory.mktemp("fgfr1_alignment_figures"))


@pytest.fixture(scope="module")
def tp53(tmp_path_factory) -> dict:
    if not TP53_RUN.exists():
        pytest.skip(f"regression run missing: {TP53_RUN}")
    return _render(TP53_RUN, tmp_path_factory.mktemp("tp53_alignment_figures"))


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

@lru_cache(maxsize=None)
def _pdf(path: Path):
    """Probe a PDF once; the wrapped document is large and probed repeatedly."""
    return probe_pdf(path)


def _texts(path: Path) -> list[str]:
    root = ET.fromstring(path.read_text(encoding="utf-8"))
    return [el.text or "" for el in root.iter() if el.tag.split("}")[-1] == "text"]


def _text(path: Path) -> str:
    return "\n".join(_texts(path))


def _fills(path: Path) -> set[str]:
    root = ET.fromstring(path.read_text(encoding="utf-8"))
    return {(el.get("fill") or "").lower() for el in root.iter()
            if el.tag.split("}")[-1] == "rect"} - {"", "none"}


def _svgs(out: Path) -> list[Path]:
    return sorted(out.glob("*.svg"))


def _wrapped_pages(out: Path) -> list[Path]:
    return sorted(out.glob("wrapped_alignment_p*.svg"),
                  key=lambda p: int(re.search(r"_p(\d+)", p.name).group(1)))


# --------------------------------------------------------------------------- #
# The reference dataset really is the dataset we claim to draw
# --------------------------------------------------------------------------- #

def test_reference_dataset_is_the_real_fgfr1_run(fgfr1):
    s = fgfr1["summary"]
    assert s["gene"] == "FGFR1"
    assert s["species"] == "Gallus gallus"
    assert s["primary"] == "NP_990841.2"
    assert s["transcript"] == "NM_205510.2"
    assert s["n_rows"] == 8, "all eight protein models must be drawn"
    assert s["n_columns"] == 823, "the complete alignment must be represented"
    assert s["n_exons"] == 17
    # The candidate is mapped onto alignment columns, not used as a residue index.
    assert s["candidate"] == "C1"
    assert s["candidate_aa"] == [31, 118]
    assert s["candidate_columns"] == [31, 118]


def test_every_expected_artefact_is_written(fgfr1):
    out = fgfr1["dir"]
    for name in SINGLE_PAGE_FIGURES:
        assert (out / f"{name}.svg").exists(), f"{name}.svg missing"
        assert (out / f"{name}.pdf").exists(), f"{name}.pdf missing"
    assert (out / "wrapped_alignment.pdf").exists()
    assert _wrapped_pages(out), "the wrapped alignment produced no pages"
    assert (out / "alignment.fasta").exists()
    assert (out / "alignment_summary.tsv").exists()


def test_the_renderer_reports_no_layout_warnings(fgfr1):
    offenders = {f["name"]: f["warnings"] for f in fgfr1["summary"]["render"]
                 if f["warnings"]}
    assert not offenders, f"layout warnings: {offenders}"


# --------------------------------------------------------------------------- #
# PDF: true vector output at a publication page size
# --------------------------------------------------------------------------- #

PDF_NAMES = [f"{n}.pdf" for n in SINGLE_PAGE_FIGURES] + ["wrapped_alignment.pdf"]


@pytest.mark.parametrize("name", PDF_NAMES)
def test_pdf_is_well_formed(fgfr1, name):
    info = _pdf(fgfr1["dir"] / name)
    assert info.ok, info.error
    assert info.header.startswith("%PDF-")
    assert info.n_pages >= 1


@pytest.mark.parametrize("name", PDF_NAMES)
def test_pdf_contains_no_raster_image(fgfr1, name):
    """The defect being replaced: an alignment exported as a rasterised page."""
    info = _pdf(fgfr1["dir"] / name)
    assert info.images == [], f"{name} embeds raster data: {info.images}"
    assert not info.is_single_raster_page


@pytest.mark.parametrize("name", PDF_NAMES)
def test_pdf_text_is_selectable_and_uses_a_referenced_font(fgfr1, name):
    info = _pdf(fgfr1["dir"] / name)
    assert info.has_vector_text, f"{name} carries no text objects"
    assert info.fonts, f"{name} references no font"
    assert info.n_text_ops >= 10


@pytest.mark.parametrize("name", PDF_NAMES)
def test_pdf_page_has_publication_dimensions(fgfr1, name):
    info = _pdf(fgfr1["dir"] / name)
    assert MIN_PAGE_IN <= info.width_in <= MAX_PAGE_IN, \
        f"{name} is {info.width_in:.1f}in wide"
    assert MIN_PAGE_IN <= info.height_in <= MAX_PAGE_IN, \
        f"{name} is {info.height_in:.1f}in tall"
    assert info.width_pt < 2000, f"{name} page was sized from pixels"


@pytest.mark.parametrize("name", PDF_NAMES)
def test_pdf_draws_vector_geometry_on_explicit_white_paper(fgfr1, name):
    info = _pdf(fgfr1["dir"] / name)
    assert info.n_path_ops >= 20, f"{name} has almost no vector geometry"
    assert info.n_fill_ops >= 20
    raw = (fgfr1["dir"] / name).read_bytes().decode("latin-1")
    assert "1 1 1 rg" in raw, f"{name} has no explicit white background"


def test_single_figure_pdfs_are_one_page_each(fgfr1):
    for name in SINGLE_PAGE_FIGURES:
        assert _pdf(fgfr1["dir"] / f"{name}.pdf").n_pages == 1


def test_wrapped_pdf_is_one_document_with_several_pages(fgfr1):
    info = _pdf(fgfr1["dir"] / "wrapped_alignment.pdf")
    expected = fgfr1["summary"]["wrapped"]["nPages"]
    assert expected > 1, "an 823-column alignment cannot be one legible page"
    assert info.n_pages == expected, \
        f"wrapped PDF has {info.n_pages} pages, expected {expected}"
    # Every page must carry its own text, not just the first one.
    assert info.n_text_ops > 100 * expected


# --------------------------------------------------------------------------- #
# SVG: standalone, explicitly painted, no stylesheet dependency
# --------------------------------------------------------------------------- #

def test_all_svg_exports_are_valid_standalone_documents(fgfr1):
    paths = _svgs(fgfr1["dir"])
    assert len(paths) >= 4, "expected the overview, the candidate figure and wrapped pages"
    for path in paths:
        svg = probe_svg(path)                       # raises on invalid XML
        assert svg["width"] and svg["height"] and svg["viewBox"], f"{path.name}: no size"
        assert svg["n_marks"] > 20, f"{path.name}: almost nothing drawn"
        assert svg["n_text"] > 10, f"{path.name}: almost no text"
        assert not svg["has_css_var"], f"{path.name}: unresolved var(--…)"
        assert not svg["has_foreign_object"], f"{path.name}: foreignObject"
        assert not svg["has_class_attr"], f"{path.name}: relies on CSS classes"
        assert not svg["has_style_element"], f"{path.name}: embedded stylesheet"
        assert not svg["has_external_css"], f"{path.name}: external stylesheet"
        assert svg["unpainted"] == [], \
            f"{path.name}: {len(svg['unpainted'])} unpainted marks"


def test_every_svg_text_element_declares_its_own_typography(fgfr1):
    for path in _svgs(fgfr1["dir"]):
        root = ET.fromstring(path.read_text(encoding="utf-8"))
        elements = [el for el in root.iter() if el.tag.split("}")[-1] == "text"]
        assert elements, f"{path.name}: no text"
        for el in elements:
            assert el.get("font-family"), f"{path.name}: text without font-family"
            assert el.get("font-size"), f"{path.name}: text without font-size"
            assert el.get("fill"), f"{path.name}: text without fill"


def test_no_svg_sets_type_below_the_publication_minimum(fgfr1):
    for path in _svgs(fgfr1["dir"]):
        root = ET.fromstring(path.read_text(encoding="utf-8"))
        sizes = [float(el.get("font-size")) for el in root.iter()
                 if el.tag.split("}")[-1] == "text"]
        assert sizes
        assert min(sizes) >= MIN_FONT_PT, f"{path.name} has {min(sizes)}pt text"


def test_every_svg_paints_explicit_white_paper(fgfr1):
    for path in _svgs(fgfr1["dir"]):
        root = ET.fromstring(path.read_text(encoding="utf-8"))
        rects = [el for el in root.iter() if el.tag.split("}")[-1] == "rect"]
        assert rects
        first = rects[0]
        assert (first.get("fill") or "").lower() in ("#ffffff", "#fff", "white"), \
            f"{path.name}: the first rectangle is not white paper"


def test_no_svg_mark_falls_outside_its_viewbox(fgfr1):
    tol = 1.0
    for path in _svgs(fgfr1["dir"]):
        root = ET.fromstring(path.read_text(encoding="utf-8"))
        _, _, vw, vh = (float(v) for v in root.get("viewBox").split())
        for el in root.iter():
            tag = el.tag.split("}")[-1]
            if tag == "rect":
                x, y = float(el.get("x")), float(el.get("y"))
                w, h = float(el.get("width")), float(el.get("height"))
                assert -tol <= x and x + w <= vw + tol, f"{path.name}: rect overflows in x"
                assert -tol <= y and y + h <= vh + tol, f"{path.name}: rect overflows in y"
            elif tag == "line":
                for attrs, lim in ((("x1", "x2"), vw), (("y1", "y2"), vh)):
                    for attr in attrs:
                        v = float(el.get(attr))
                        assert -tol <= v <= lim + tol, \
                            f"{path.name}: line overflows ({attr}={v})"
            elif tag == "text":
                x, y = float(el.get("x")), float(el.get("y"))
                assert -tol <= x <= vw + tol, f"{path.name}: text overflows in x"
                assert -tol <= y <= vh + tol, f"{path.name}: text overflows in y"


# --------------------------------------------------------------------------- #
# Figure 1 — the overview is a real alignment, not eight coverage bars
# --------------------------------------------------------------------------- #

def test_overview_names_every_protein_model_with_its_transcript_and_length(fgfr1):
    text = _text(fgfr1["dir"] / "full_isoform_alignment.svg")
    for protein_id in fgfr1["summary"]["protein_ids"]:
        assert protein_id in text, f"protein model {protein_id} is not labelled"
    for transcript_id in ("NM_205510.2", "XM_015297361.4", "XM_046903352.1"):
        assert transcript_id in text, f"transcript {transcript_id} is not labelled"
    assert "817 aa" in text
    assert "curated" in text and "predicted" in text


def test_overview_reports_identity_to_the_primary_for_every_row(fgfr1):
    lines = _texts(fgfr1["dir"] / "full_isoform_alignment.svg")
    identities = [line for line in lines if re.fullmatch(r"\d{1,3}(\.\d)?%", line.strip())]
    assert len(identities) >= fgfr1["summary"]["n_rows"], \
        f"only {len(identities)} identity values for 8 models"
    assert any(v.startswith("100") for v in identities), "the primary must read 100%"
    assert any(not v.startswith("100") for v in identities), \
        "every model at 100% means identity is not being computed"


def test_overview_carries_conservation_difference_and_gap_tracks(fgfr1):
    text = _text(fgfr1["dir"] / "full_isoform_alignment.svg")
    for label in ("Conservation", "Difference density", "Gap density",
                  "Variable columns", "Primary aa mapping"):
        assert label in text, f"aggregate track missing: {label}"
    # The tracks must report a real observed range, not a hard-coded axis.
    assert re.search(r"0–1, max 0\.\d\d", text), "no observed maximum on the density tracks"


def test_overview_states_both_coordinate_systems(fgfr1):
    text = _text(fgfr1["dir"] / "full_isoform_alignment.svg")
    assert "Alignment column" in text
    assert "NP_990841.2 residue position (aa)" in text
    assert "823 alignment columns" in text


def test_overview_labels_the_candidate_intervals(fgfr1):
    text = _text(fgfr1["dir"] / "full_isoform_alignment.svg")
    assert "C1 · aa 31–118" in text
    assert "Candidates" in text
    assert "Exploratory analysis" in text, "the exploratory status must be stated"


def test_overview_legend_explains_the_alignment_colours(fgfr1):
    text = _text(fgfr1["dir"] / "full_isoform_alignment.svg")
    for entry in ("primary protein (aligned residues)",
                  "alternative model (aligned residues)",
                  "residue differing from primary",
                  "gap in this isoform",
                  "variable column (models disagree)",
                  "conservation (consensus fraction)",
                  "exploratory candidate interval"):
        assert entry in text, f"legend entry missing: {entry}"


def test_overview_is_not_a_row_of_uniform_coverage_bars(fgfr1):
    """Per-column structure, not one flat rectangle per protein model."""
    path = fgfr1["dir"] / "full_isoform_alignment.svg"
    svg = probe_svg(path)
    assert svg["n_marks"] > 300, \
        f"only {svg['n_marks']} marks: this cannot resolve 823 columns"
    text = _text(path)
    # Difference and gap structure is quantified, not merely mentioned.
    assert re.search(r"\d+ variable and \d+ gap-containing columns of 823", text), \
        "the overview does not report its difference and gap structure"
    assert re.search(r"Major variable blocks: \d+ · major gap blocks: \d+", text)


def test_overview_uses_the_alignment_palette_and_no_colour_per_isoform(fgfr1):
    """Colour encodes alignment state; a colour per model would encode nothing."""
    fills = _fills(fgfr1["dir"] / "full_isoform_alignment.svg")
    # The semantic colours that must all be present.
    for colour in ("#39536e",   # primary protein
                   "#d55e00",   # residue differing from the primary
                   "#eef1f5"):  # gap in this isoform
        assert colour in fills, f"semantic colour {colour} is not used"
    # Eight arbitrary categorical hues would push the palette far wider than the
    # handful of roles the figure actually distinguishes.
    assert len(fills) <= 16, f"{len(fills)} distinct rectangle fills: palette is categorical"


# --------------------------------------------------------------------------- #
# Figure 2 — wrapped residue-level alignment
# --------------------------------------------------------------------------- #

def test_wrapped_export_is_paginated_into_blocks(fgfr1):
    layout = fgfr1["summary"]["wrapped"]
    assert 60 <= layout["colsPerBlock"] <= 100, \
        f"{layout['colsPerBlock']} columns per block is outside the legible range"
    assert layout["nBlocks"] * layout["colsPerBlock"] >= fgfr1["summary"]["n_columns"]
    assert layout["n_pages_written"] == layout["nPages"]
    assert len(_wrapped_pages(fgfr1["dir"])) == layout["nPages"]


def test_wrapped_export_shows_residue_letters(fgfr1):
    letters = 0
    for path in _wrapped_pages(fgfr1["dir"]):
        letters += sum(1 for t in _texts(path)
                       if len(t.strip()) == 1 and t.strip().isalpha())
    # Eight models over 823 columns cannot be shown with a few hundred glyphs.
    assert letters > 4000, f"only {letters} residue letters across the wrapped export"


def test_wrapped_export_repeats_the_sequence_labels_in_every_block(fgfr1):
    for path in _wrapped_pages(fgfr1["dir"]):
        texts = _texts(path)
        blocks = [t for t in texts if t.startswith("cols ")]
        assert blocks, f"{path.name}: no block header"
        for protein_id in fgfr1["summary"]["protein_ids"]:
            occurrences = sum(1 for t in texts if t.startswith(protein_id))
            assert occurrences == len(blocks), \
                f"{path.name}: {protein_id} labelled {occurrences}x for {len(blocks)} blocks"


def test_wrapped_export_covers_every_alignment_column(fgfr1):
    covered: set[int] = set()
    for path in _wrapped_pages(fgfr1["dir"]):
        for t in _texts(path):
            m = re.fullmatch(r"cols (\d+)–(\d+)", t.strip())
            if m:
                covered.update(range(int(m.group(1)), int(m.group(2)) + 1))
    n_cols = fgfr1["summary"]["n_columns"]
    assert covered == set(range(1, n_cols + 1)), \
        f"the wrapped export covers {len(covered)} of {n_cols} alignment columns"


def test_wrapped_export_annotates_every_block(fgfr1):
    for path in _wrapped_pages(fgfr1["dir"]):
        texts = _texts(path)
        n_blocks = sum(1 for t in texts if t.startswith("cols "))
        for row in ("variable", "conservation", "primary aa"):
            # One occurrence per block, plus at most one in the legend.
            occurrences = sum(1 for t in texts if t == row)
            assert n_blocks <= occurrences <= n_blocks + 1, \
                f"{path.name}: the {row!r} row is not repeated per block"
        assert "page 1 of" in "\n".join(texts) or "page " in "\n".join(texts)


def test_wrapped_export_marks_the_candidate_intervals(fgfr1):
    text = "\n".join(_text(p) for p in _wrapped_pages(fgfr1["dir"]))
    assert "C1" in text
    assert "exploratory candidate interval" in text
    assert "no isoform difference shown here is a validated splicing event" in text


# --------------------------------------------------------------------------- #
# Figure 3 — candidate-focused detail
# --------------------------------------------------------------------------- #

def test_candidate_figure_shows_affected_and_unaffected_models(fgfr1):
    text = _text(fgfr1["dir"] / "candidate_alignment_detail.svg")
    summary = fgfr1["summary"]
    assert summary["affected"] and summary["unaffected"], \
        "the reference candidate must have both affected and unaffected models"
    for protein_id in summary["affected"]:
        assert f"{protein_id} · affected" in text, f"{protein_id} not marked affected"
    for protein_id in summary["unaffected"]:
        assert f"{protein_id} · unaffected" in text, f"{protein_id} not marked unaffected"
    assert f"{summary['primary']} ★ · reference" in text, \
        "the primary protein must be labelled as the reference"


def test_candidate_figure_carries_both_coordinate_systems(fgfr1):
    text = _text(fgfr1["dir"] / "candidate_alignment_detail.svg")
    assert "primary aa 31–118 on NP_990841.2 = alignment columns 31–118" in text
    assert "alignment column" in text
    assert "primary aa" in text
    # Both rulers must carry real numbers from the two systems.
    assert re.search(r"^cols \d+–\d+$", text, re.MULTILINE)


def test_candidate_figure_shows_residues_and_reports_local_identity(fgfr1):
    path = fgfr1["dir"] / "candidate_alignment_detail.svg"
    letters = sum(1 for t in _texts(path)
                  if len(t.strip()) == 1 and t.strip().isalpha())
    assert letters > 500, f"only {letters} residue letters in the candidate figure"
    text = _text(path)
    assert "identity within the interval" in text
    # A model with no residue at all inside the interval is not "0% identical".
    assert "no aligned residues" in text
    assert re.search(r"\d+ variable and \d+ gap-containing columns inside the candidate",
                     text)


def test_candidate_figure_states_the_exon_association(fgfr1):
    text = _text(fgfr1["dir"] / "candidate_alignment_detail.svg")
    assert "primary exons" in text
    assert "E1" in text and "E2" in text
    assert "coding exon of the primary protein" in text


def test_candidate_figure_does_not_replace_the_overview(fgfr1):
    """Both figures must exist, and the detail must be the narrower one."""
    detail = probe_svg(fgfr1["dir"] / "candidate_alignment_detail.svg")
    overview = probe_svg(fgfr1["dir"] / "full_isoform_alignment.svg")
    assert "823 alignment columns" in overview["raw"]
    assert "823" not in _text(fgfr1["dir"] / "candidate_alignment_detail.svg")
    assert detail["raw"] != overview["raw"]


# --------------------------------------------------------------------------- #
# Rasterised output
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def fgfr1_pngs(fgfr1, tmp_path_factory) -> Path:
    if shutil.which("rsvg-convert") is None:
        pytest.skip("rsvg-convert is required to rasterise the SVG exports")
    out = tmp_path_factory.mktemp("fgfr1_alignment_pngs")
    for name in SINGLE_PAGE_FIGURES:
        subprocess.run(
            ["rsvg-convert", "--dpi-x", "300", "--dpi-y", "300", "--zoom", "4.1667",
             str(fgfr1["dir"] / f"{name}.svg"), "-o", str(out / f"{name}.png")],
            check=True, capture_output=True,
        )
    return out


@pytest.mark.parametrize("name", SINGLE_PAGE_FIGURES)
def test_rasterised_figure_has_publication_resolution_and_real_content(fgfr1_pngs, name):
    png = probe_png(fgfr1_pngs / f"{name}.png")
    # 300 dpi over a page of at least two inches.
    assert png["width"] >= 2 * 300
    assert png["height"] >= 300
    assert not png["is_blank"], f"{name}.png is a blank page"
    assert png["contrast"] > 0.5, f"{name}.png has almost no contrast"
    # Predominantly white paper with ink on it, not a solid block of colour.
    assert png["frac_light"] > 0.5, f"{name}.png is not mostly paper"
    assert png["frac_dark"] < 0.5


# --------------------------------------------------------------------------- #
# Source data shipped alongside the figures
# --------------------------------------------------------------------------- #

def test_alignment_fasta_holds_every_aligned_sequence(fgfr1):
    lines = (fgfr1["dir"] / "alignment.fasta").read_text().splitlines()
    headers = [line for line in lines if line.startswith(">")]
    assert len(headers) == fgfr1["summary"]["n_rows"]
    assert headers[0].startswith(f">{fgfr1['summary']['primary']}")
    assert "primary" in headers[0] and "curated" in headers[0]
    sequences: dict[str, str] = {}
    current = None
    for line in lines:
        if line.startswith(">"):
            current = line[1:].split()[0]
            sequences[current] = ""
        elif current:
            sequences[current] += line.strip()
    assert set(sequences) == set(fgfr1["summary"]["protein_ids"])
    for protein_id, seq in sequences.items():
        assert len(seq) == fgfr1["summary"]["n_columns"], \
            f"{protein_id} is not aligned to {fgfr1['summary']['n_columns']} columns"


def test_alignment_summary_tsv_is_the_figure_source_table(fgfr1):
    rows = [line.split("\t") for line in
            (fgfr1["dir"] / "alignment_summary.tsv").read_text().splitlines() if line]
    header, body = rows[0], rows[1:]
    for column in ("protein_id", "transcript_id", "curation_status", "protein_length",
                   "alignment_columns", "identity_to_primary_pct", "gap_columns",
                   "differing_columns", "candidate_id", "candidate_alignment_columns",
                   "candidate_role"):
        assert column in header, f"missing TSV column: {column}"
    assert len(body) == fgfr1["summary"]["n_rows"]
    assert all(len(r) == len(header) for r in body), "ragged TSV"
    idx = header.index("differing_columns")
    assert any(int(r[idx]) > 0 for r in body), \
        "no model differs from the primary anywhere: the table carries no difference data"


# --------------------------------------------------------------------------- #
# Second run: a structurally different alignment renders the same figure set
# --------------------------------------------------------------------------- #

def test_tp53_run_renders_the_same_figure_set(tp53):
    s = tp53["summary"]
    assert s["gene"] == "TP53"
    assert s["species"] == "Danio rerio"
    assert s["n_rows"] == 12, "every TP53 protein model must be drawn"
    assert s["n_columns"] == 374
    for name in SINGLE_PAGE_FIGURES:
        assert (tp53["dir"] / f"{name}.svg").exists()
        assert (tp53["dir"] / f"{name}.pdf").exists()
    assert (tp53["dir"] / "wrapped_alignment.pdf").exists()
    assert _wrapped_pages(tp53["dir"])


def test_tp53_pdfs_are_vector_with_selectable_text(tp53):
    for name in PDF_NAMES:
        info = _pdf(tp53["dir"] / name)
        assert info.ok, info.error
        assert info.images == [], f"{name} embeds raster data"
        assert info.has_vector_text
        assert MIN_PAGE_IN <= info.width_in <= MAX_PAGE_IN
        assert MIN_PAGE_IN <= info.height_in <= MAX_PAGE_IN


def test_tp53_svgs_are_standalone_and_legible(tp53):
    for path in _svgs(tp53["dir"]):
        svg = probe_svg(path)
        assert not svg["has_css_var"] and not svg["has_class_attr"]
        assert svg["unpainted"] == []
        root = ET.fromstring(path.read_text(encoding="utf-8"))
        sizes = [float(el.get("font-size")) for el in root.iter()
                 if el.tag.split("}")[-1] == "text"]
        assert min(sizes) >= MIN_FONT_PT, f"{path.name} has {min(sizes)}pt text"


def test_tp53_overview_shows_all_twelve_models_with_identities(tp53):
    text = _text(tp53["dir"] / "full_isoform_alignment.svg")
    for protein_id in tp53["summary"]["protein_ids"]:
        assert protein_id in text, f"{protein_id} is missing from the TP53 overview"
    assert "Conservation" in text and "Gap density" in text
    # TP53 models differing by a single residue must not be reported as identical.
    assert re.search(r"\d\d\.\d%", text), "sub-percent identity differences are rounded away"
