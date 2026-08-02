"""Minimal PDF inspector for export validation, using only the standard library.

Enough of the PDF structure is parsed to answer the questions the export tests
ask: is this a genuine vector figure at a sensible physical size, or is it a
single full-page raster screenshot?

Deliberately dependency-free so the export tests run anywhere the project runs.
"""

from __future__ import annotations

import re
import zlib
from dataclasses import dataclass, field
from pathlib import Path

PT_PER_INCH = 72.0


@dataclass
class PdfInfo:
    path: Path
    ok: bool
    header: str = ""
    n_pages: int = 0
    # page geometry in PDF points (1/72 inch)
    width_pt: float = 0.0
    height_pt: float = 0.0
    fonts: list[str] = field(default_factory=list)
    images: list[str] = field(default_factory=list)
    text: str = ""
    n_text_ops: int = 0
    n_path_ops: int = 0
    n_fill_ops: int = 0
    error: str = ""

    @property
    def width_in(self) -> float:
        return self.width_pt / PT_PER_INCH

    @property
    def height_in(self) -> float:
        return self.height_pt / PT_PER_INCH

    @property
    def width_mm(self) -> float:
        return self.width_in * 25.4

    @property
    def aspect(self) -> float:
        return self.width_pt / self.height_pt if self.height_pt else 0.0

    @property
    def has_raster_image(self) -> bool:
        return bool(self.images)

    @property
    def is_single_raster_page(self) -> bool:
        """True for a page whose visible content is one embedded bitmap."""
        return bool(self.images) and self.n_text_ops == 0

    @property
    def has_vector_text(self) -> bool:
        return self.n_text_ops > 0 and bool(self.fonts)

    def contains(self, needle: str) -> bool:
        """Case-insensitive search of the extracted text layer."""
        return needle.lower() in self.text.lower()


def _decode_streams(raw: bytes) -> str:
    """Concatenate every content stream we can decode, as latin-1 text."""
    out = []
    for match in re.finditer(rb"stream\r?\n(.*?)\r?\nendstream", raw, re.DOTALL):
        body = match.group(1)
        for candidate in (body,):
            try:
                out.append(zlib.decompress(candidate).decode("latin-1", "replace"))
                break
            except zlib.error:
                # Uncompressed content stream (our own writer emits these).
                out.append(candidate.decode("latin-1", "replace"))
                break
    return "\n".join(out)


_TJ_SIMPLE = re.compile(r"\((?:\\.|[^\\()])*\)\s*Tj")
_TJ_ARRAY = re.compile(r"\[(.*?)\]\s*TJ", re.DOTALL)
_STR_IN_ARRAY = re.compile(r"\((?:\\.|[^\\()])*\)")
# matplotlib writes glyph-by-glyph shows inside a single text object
_TEXT_OBJ = re.compile(r"\bBT\b")


def _unescape(s: str) -> str:
    s = s[1:-1]  # strip the parentheses
    s = re.sub(r"\\([()\\])", r"\1", s)
    s = s.replace(r"\n", " ").replace(r"\r", " ").replace(r"\t", " ")
    return re.sub(r"\\[0-7]{1,3}", "?", s)


def probe_pdf(path) -> PdfInfo:
    path = Path(path)
    if not path.exists():
        return PdfInfo(path=path, ok=False, error="file does not exist")
    raw = path.read_bytes()
    info = PdfInfo(path=path, ok=False, header=raw[:8].decode("latin-1", "replace"))
    if not raw.startswith(b"%PDF-"):
        info.error = "missing %PDF- header"
        return info
    if b"%%EOF" not in raw[-2048:]:
        info.error = "missing %%EOF trailer"
        return info

    info.n_pages = len(re.findall(rb"/Type\s*/Page[^s]", raw))
    # A malformed writer may omit /Type /Page; fall back to the page tree count.
    if info.n_pages == 0:
        m = re.search(rb"/Count\s+(\d+)", raw)
        info.n_pages = int(m.group(1)) if m else 0

    box = re.search(rb"/MediaBox\s*\[\s*([\d.+-]+)\s+([\d.+-]+)\s+"
                    rb"([\d.+-]+)\s+([\d.+-]+)\s*\]", raw)
    if box:
        x0, y0, x1, y1 = (float(v) for v in box.groups())
        info.width_pt = abs(x1 - x0)
        info.height_pt = abs(y1 - y0)
    else:
        info.error = "no /MediaBox"
        return info

    info.fonts = sorted({m.decode("latin-1") for m in
                         re.findall(rb"/BaseFont\s*/([#\w+.-]+)", raw)})
    info.images = sorted({m.decode("latin-1") for m in
                          re.findall(rb"/(DCTDecode|JPXDecode|CCITTFaxDecode)", raw)})
    # An /Image XObject with a pixel-array filter is also a raster.
    if re.search(rb"/Subtype\s*/Image", raw) and not info.images:
        if re.search(rb"/Filter\s*/FlateDecode", raw):
            info.images.append("FlateImage")

    content = _decode_streams(raw)
    info.n_text_ops = len(_TEXT_OBJ.findall(content))
    # Vector geometry is rectangles (re), Bezier curves (c) and line segments (l).
    info.n_path_ops = sum(len(re.findall(rf"(?<![\w.]){op}(?![\w.])", content))
                          for op in ("re", "c", "l"))
    info.n_fill_ops = len(re.findall(r"(?<![\w.])[fFbBSs]\*?(?![\w.])", content))

    chunks = []
    for m in _TJ_SIMPLE.finditer(content):
        chunks.append(_unescape(m.group(0).rsplit("Tj", 1)[0].strip()))
    for m in _TJ_ARRAY.finditer(content):
        chunks.append("".join(_unescape(s) for s in _STR_IN_ARRAY.findall(m.group(1))))
    info.text = "\n".join(c for c in chunks if c.strip())

    info.ok = True
    return info


def probe_svg(path) -> dict:
    """Structural checks for an exported standalone SVG."""
    import xml.etree.ElementTree as ET

    path = Path(path)
    text = path.read_text(encoding="utf-8")
    root = ET.fromstring(text)  # raises on invalid XML
    tags = [el.tag.split("}")[-1] for el in root.iter()]

    def attrs_of(tag):
        return [el.attrib for el in root.iter() if el.tag.split("}")[-1] == tag]

    shape_tags = ("rect", "circle", "path", "polygon", "ellipse", "line")
    unpainted = []
    for tag in shape_tags:
        for a in attrs_of(tag):
            has_fill = "fill" in a
            has_stroke = "stroke" in a
            if tag == "line":
                if not has_stroke:
                    unpainted.append(tag)
            elif not has_fill and not has_stroke:
                unpainted.append(tag)

    texts = [(el.text or "") for el in root.iter() if el.tag.split("}")[-1] == "text"]
    return {
        "width": root.get("width"),
        "height": root.get("height"),
        "viewBox": root.get("viewBox"),
        "n_marks": sum(tags.count(t) for t in shape_tags),
        "n_text": tags.count("text"),
        "texts": texts,
        "unpainted": unpainted,
        "has_foreign_object": "foreignObject" in tags,
        "has_style_element": "style" in tags,
        "has_class_attr": any("class" in el.attrib for el in root.iter()),
        "has_css_var": "var(--" in text,
        "has_external_css": "<?xml-stylesheet" in text or "@import" in text,
        "raw": text,
    }


def probe_png(path) -> dict:
    """Dimensions, dpi (pHYs) and simple contrast statistics for a PNG."""
    from PIL import Image

    path = Path(path)
    with Image.open(path) as im:
        dpi = im.info.get("dpi")
        grey = im.convert("L")
        histogram = grey.histogram()
        total = sum(histogram) or 1
        dark = sum(histogram[:96]) / total
        light = sum(histogram[192:]) / total
        extrema = grey.getextrema()
        return {
            "width": im.width,
            "height": im.height,
            "dpi": tuple(round(v) for v in dpi) if dpi else None,
            "mode": im.mode,
            "aspect": im.width / im.height if im.height else 0,
            "frac_dark": dark,
            "frac_light": light,
            "contrast": (extrema[1] - extrema[0]) / 255,
            "is_blank": extrema[0] == extrema[1],
        }


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        p = Path(arg)
        if p.suffix.lower() == ".pdf":
            i = probe_pdf(p)
            print(f"{p.name}: ok={i.ok} pages={i.n_pages} "
                  f"{i.width_pt:.0f}x{i.height_pt:.0f}pt "
                  f"({i.width_in:.1f}x{i.height_in:.1f}in) fonts={i.fonts} "
                  f"images={i.images} text_objs={i.n_text_ops} paths={i.n_path_ops}")
            if i.text:
                print("   text:", i.text[:300].replace("\n", " | "))
            if i.error:
                print("   error:", i.error)
        elif p.suffix.lower() == ".svg":
            print(f"{p.name}:", {k: v for k, v in probe_svg(p).items() if k != "raw"})
        elif p.suffix.lower() == ".png":
            print(f"{p.name}:", probe_png(p))
