from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(__file__).resolve().parents[1]
FGFR1_RUN = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"
TP53_RUN = ROOT / "runs" / "2026-07-21_1436_custom_run"
FIGURE_DIR = Path("results") / "generic_gene_analysis" / "figures" / "main"

RENDER_DPI = 150
MAX_BBOX_SLACK_PX = 12


def _require_tools():
    if not shutil.which("rsvg-convert"):
        pytest.skip("rsvg-convert is required to rasterise the SVG")
    if not (shutil.which("sips") or shutil.which("pdftoppm")):
        pytest.skip("no PDF rasteriser available (sips or pdftoppm)")
    try:
        import PIL  # noqa: F401
    except ImportError:
        pytest.skip("Pillow is required for the pixel comparison")


def _svg_to_png(svg: Path, out: Path, dpi: int = RENDER_DPI) -> None:
    zoom = dpi / 72.0
    subprocess.run(["rsvg-convert", "--zoom", f"{zoom:.4f}", str(svg), "-o", str(out)],
                   capture_output=True, check=True)


def _pdf_to_png(pdf: Path, out: Path, dpi: int = RENDER_DPI) -> None:
    if shutil.which("pdftoppm"):
        stem = out.with_suffix("")
        subprocess.run(["pdftoppm", "-r", str(dpi), "-png", "-singlefile",
                        str(pdf), str(stem)], capture_output=True, check=True)
        return
    subprocess.run(["sips", "-s", "format", "png",
                    "-s", "dpiHeight", str(dpi), "-s", "dpiWidth", str(dpi),
                    str(pdf), "--out", str(out)], capture_output=True, check=True)


def _ink(img) -> float:
    grey = img.convert("L")
    dark = sum(c for v, c in enumerate(grey.histogram()) if v < 235)
    return dark / float(grey.width * grey.height)


def _ink_bbox(img):
    # After thresholding, ink is non-zero and paper is zero, so getbbox() returns
    # the extent of the drawn content.
    mask = img.convert("L").point(lambda v: 0 if v > 235 else 255)
    box = mask.getbbox()
    if box is None:
        return None
    x0, y0, x1, y1 = box
    return (x0 / mask.width, y0 / mask.height, x1 / mask.width, y1 / mask.height)


def _figure_pairs(run_dir: Path) -> list[tuple[str, Path, Path]]:
    out = []
    for svg in sorted((run_dir / FIGURE_DIR).glob("main_*.svg")):
        pdf = svg.with_suffix(".pdf")
        if pdf.exists():
            out.append((svg.stem, svg, pdf))
    return out


def _cases() -> list[tuple[str, Path, Path]]:
    cases: list[tuple[str, Path, Path]] = []
    for run in (FGFR1_RUN, TP53_RUN):
        if run.exists():
            cases.extend(_figure_pairs(run))
    return cases


CASES = _cases()
IDS = [c[0] for c in CASES]


PARITY_HARNESS = Path(__file__).with_name("check_backend_parity.mjs")
MODEL = Path("website_indices") / "generic" / "protein_coordinate_model.json"

RUNS = [r for r in (FGFR1_RUN, TP53_RUN) if (r / MODEL).exists()]


@pytest.mark.skipif(not RUNS, reason="no coordinate model available")
@pytest.mark.parametrize("run", RUNS, ids=[r.name for r in RUNS])
def test_both_backends_draw_every_mark_of_every_figure(run):
    if not shutil.which("node"):
        pytest.skip("node is required to exercise the figure backends")
    proc = subprocess.run(["node", str(PARITY_HARNESS), str(run / MODEL)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, \
        f"backend parity failed for {run.name}:\n{proc.stdout}\n{proc.stderr}"
    # A harness that silently checked nothing would also exit zero.
    assert "every mark drawn by both backends" in proc.stdout
    assert proc.stdout.count("both backends agree") >= 5, \
        f"too few figures checked for {run.name}:\n{proc.stdout}"


@pytest.mark.skipif(not CASES, reason="no rendered main figures found")
@pytest.mark.parametrize("name,svg,pdf", CASES, ids=IDS or None)
def test_neither_format_rasterises_to_a_blank_page(name, svg, pdf, tmp_path):
    _require_tools()
    from PIL import Image

    a, b = tmp_path / "svg.png", tmp_path / "pdf.png"
    _svg_to_png(svg, a)
    _pdf_to_png(pdf, b)
    with Image.open(a) as si, Image.open(b) as pi:
        assert _ink(si) > 0.005, f"{name}: the SVG rasterises to an almost blank page"
        assert _ink(pi) > 0.005, f"{name}: the PDF rasterises to an almost blank page"


@pytest.mark.skipif(not CASES, reason="no rendered main figures found")
@pytest.mark.parametrize("name,svg,pdf", CASES, ids=IDS or None)
def test_pdf_and_svg_place_their_content_in_the_same_area(name, svg, pdf, tmp_path):
    _require_tools()
    from PIL import Image

    a, b = tmp_path / "svg.png", tmp_path / "pdf.png"
    _svg_to_png(svg, a)
    _pdf_to_png(pdf, b)
    with Image.open(a) as si, Image.open(b) as pi:
        # Compared as a shape, not in pixels: the two rasterisers do not always
        # honour the same resolution, but the page proportions must still match.
        s_aspect, p_aspect = si.width / si.height, pi.width / pi.height
        assert abs(s_aspect - p_aspect) / s_aspect <= 0.02, (
            f"{name}: page proportions differ (SVG {s_aspect:.3f} vs "
            f"PDF {p_aspect:.3f}) — one format is distorted")
        sb, pb = _ink_bbox(si), _ink_bbox(pi)
        assert sb and pb, f"{name}: one format has no inked content"
        for i, edge in enumerate(("left", "top", "right", "bottom")):
            reference = min(si.width, pi.width) if i % 2 == 0 \
                else min(si.height, pi.height)
            slack = MAX_BBOX_SLACK_PX / reference
            assert abs(sb[i] - pb[i]) <= slack, (
                f"{name}: the {edge} edge of the drawn content differs between "
                f"SVG and PDF (SVG {sb[i]:.3f} vs PDF {pb[i]:.3f})")


@pytest.mark.skipif(not CASES, reason="no rendered main figures found")
@pytest.mark.parametrize("name,svg,pdf", CASES, ids=IDS or None)
def test_neither_format_paints_a_dark_page(name, svg, pdf, tmp_path):
    _require_tools()
    from PIL import Image

    for label, src, conv in (("SVG", svg, _svg_to_png), ("PDF", pdf, _pdf_to_png)):
        out = tmp_path / f"{label}.png"
        conv(src, out)
        with Image.open(out) as img:
            grey = img.convert("L")
            hist = grey.histogram()
            total = grey.width * grey.height
            very_dark = sum(hist[:60]) / total
            assert very_dark < 0.25, \
                f"{name}: {label} is {very_dark:.0%} near-black — paint fell back to black"
