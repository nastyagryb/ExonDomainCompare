// Shared, dependency-free SVG export helpers for the scientific viewers
// (Exon Map style). Extracted so the Boundary explorer can reuse the exact same
// SVG/PNG/PDF export behaviour without duplicating logic or touching the
// accepted Domain Architecture component.

export function niceStep(span) {
  const raw = span / 9;
  const pow = Math.pow(10, Math.floor(Math.log10(Math.max(1, raw))));
  const n = raw / pow;
  const m = n >= 5 ? 5 : n >= 2 ? 2 : 1;
  return Math.max(5, m * pow);
}

export function overlaps(a, b) {
  return a && b && a.start != null && b.start != null && a.start <= b.end && a.end >= b.start;
}

// Explicit ink/line colours for SVG presentation attributes. Exported SVG files
// are standalone documents, so a `var(--line-strong)` reference there resolves to
// nothing and the stroke silently disappears — figures must carry literal values.
export const INK = "#1c2433";
export const LINE = "#e3e8f0";
export const LINE_STRONG = "#d2d9e6";
export const FIG_FONT = "Helvetica, Arial, 'DejaVu Sans', sans-serif";

export function downloadBlob(blob, name) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

// Serialize an <svg> element into a standalone, self-contained SVG document:
// explicit width/height, a literal font stack, a white paper background and no
// dependency on the page's stylesheet, custom properties or foreignObject.
export function serializeSvg(svgEl, opts = {}) {
  const clone = svgEl.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
  clone.removeAttribute("class");
  // A viewBox-only SVG has no intrinsic size in Preview or a PDF converter.
  const vb = (clone.getAttribute("viewBox") || "").trim().split(/[\s,]+/).map(Number);
  const w = opts.width ?? (vb.length === 4 ? vb[2] : svgEl.clientWidth || 1000);
  const h = opts.height ?? (vb.length === 4 ? vb[3] : svgEl.clientHeight || 400);
  clone.setAttribute("width", String(w));
  clone.setAttribute("height", String(h));
  if (vb.length !== 4) clone.setAttribute("viewBox", `0 0 ${w} ${h}`);

  const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
  style.textContent = `text{font-family:${FIG_FONT};fill:${INK}}`;
  const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  bg.setAttribute("x", String(vb.length === 4 ? vb[0] : 0));
  bg.setAttribute("y", String(vb.length === 4 ? vb[1] : 0));
  bg.setAttribute("width", String(w));
  bg.setAttribute("height", String(h));
  bg.setAttribute("fill", "#ffffff");
  clone.insertBefore(bg, clone.firstChild);
  clone.insertBefore(style, clone.firstChild);
  return `<?xml version="1.0" encoding="UTF-8"?>\n` + new XMLSerializer().serializeToString(clone);
}

// Rasterize an <svg> element to a canvas and hand the blob to cb(blob, w, h).
// scale 4 on a ~1000 px figure box yields ~300 dpi at typical journal widths.
export function rasterizeSvg(svgEl, W, H, cb, type = "image/png", quality, scale = 4) {
  if (!svgEl) return;
  const str = serializeSvg(svgEl, { width: W, height: H });
  const img = new Image();
  const url = URL.createObjectURL(new Blob([str], { type: "image/svg+xml;charset=utf-8" }));
  img.onload = () => {
    const canvas = document.createElement("canvas");
    canvas.width = W * scale; canvas.height = H * scale;
    const ctx = canvas.getContext("2d");
    ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    canvas.toBlob((b) => cb(b, canvas.width, canvas.height), type, quality);
    URL.revokeObjectURL(url);
  };
  img.src = url;
}

// PDF export deliberately has no raster path any more.
//
// The previous implementation rasterised the live SVG to JPEG and wrapped that
// single image in a PDF page sized from canvas pixels. The result was a
// full-page 72 dpi photograph of the interface — no vector geometry, no
// selectable text, a page several feet wide, and solid black blocks wherever a
// fill came from the application stylesheet rather than from the mark itself.
//
// Publication PDFs are now produced from a declarative figure specification
// (figureSpec.js) via figureExport.js, which emits real path and text operators.
// Anything that needs a PDF must build a figure specification.
