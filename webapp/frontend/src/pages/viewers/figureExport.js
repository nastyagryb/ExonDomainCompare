// Download helpers for the publication figure specifications.
//
// All three formats come from one figure specification (see figureSpec.js), so a
// downloaded SVG, PDF and PNG are views of the same figure rather than three
// independent pipelines. Nothing here screenshots the interactive component.

import { createFigure, renderPdfPages } from "./figureSpec";
import { downloadBlob } from "./plotExport";

/** Publication raster resolution. */
export const EXPORT_DPI = 300;

export function downloadFigureSvg(fig, stem) {
  downloadBlob(new Blob([fig.toSvg()], { type: "image/svg+xml;charset=utf-8" }),
    `${stem}.svg`);
}

export function downloadFigurePdf(fig, stem) {
  downloadBlob(new Blob([fig.toPdf()], { type: "application/pdf" }), `${stem}.pdf`);
}

/** One multi-page PDF from a sequence of figure specifications. */
export function downloadFigurePdfPages(figures, stem) {
  downloadBlob(new Blob([renderPdfPages(figures)], { type: "application/pdf" }),
    `${stem}.pdf`);
}

/**
 * Rasterise the figure's own SVG at 300 dpi.
 *
 * The figure is specified in points, so the pixel size follows from the physical
 * geometry: width_pt / 72 * dpi. The result therefore has the correct physical
 * aspect ratio and carries no interface chrome.
 */
export function downloadFigurePng(fig, stem, { dpi = EXPORT_DPI } = {}) {
  const scale = dpi / 72;
  const pxW = Math.round(fig.width * scale);
  const pxH = Math.round(fig.height * scale);
  const svg = fig.toSvg();
  const url = URL.createObjectURL(new Blob([svg], { type: "image/svg+xml;charset=utf-8" }));
  const img = new Image();
  img.onload = () => {
    const canvas = document.createElement("canvas");
    canvas.width = pxW;
    canvas.height = pxH;
    const ctx = canvas.getContext("2d");
    // The SVG already paints its own white paper, but an opaque canvas keeps
    // viewers that ignore the background rect consistent with the PDF.
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, pxW, pxH);
    ctx.drawImage(img, 0, 0, pxW, pxH);
    canvas.toBlob((blob) => {
      URL.revokeObjectURL(url);
      if (blob) downloadBlob(blob, `${stem}.png`);
    }, "image/png");
  };
  img.onerror = () => URL.revokeObjectURL(url);
  img.src = url;
}

export function downloadFigureTsv(text, stem) {
  downloadBlob(new Blob([text], { type: "text/tab-separated-values;charset=utf-8" }),
    `${stem}.tsv`);
}

/**
 * Replay a self-contained SVG into a figure specification, so it can be written
 * as a real vector PDF.
 *
 * This is lossless for the figures this project generates, because they are built
 * from the same four primitives with explicit presentation attributes and no CSS.
 * A figure that relied on a stylesheet could not be converted — and would already
 * be broken as an SVG export.
 */
export function svgToFigure(svgText) {
  const doc = new DOMParser().parseFromString(svgText, "image/svg+xml");
  const root = doc.documentElement;
  if (root.tagName === "parsererror") throw new Error("figure SVG is not valid XML");
  const [, , vbW, vbH] = (root.getAttribute("viewBox") || "").split(/\s+/).map(Number);
  const width = Number(root.getAttribute("width")) || vbW || 600;
  const height = Number(root.getAttribute("height")) || vbH || 400;

  // A converted figure keeps its own geometry rather than a preset's, so pass the
  // width explicitly and lift the height clamp by using the widest preset.
  const fig = createFigure({ preset: "full", width, height });
  fig.resize(height);

  const num = (el, name, dflt = 0) => {
    const v = parseFloat(el.getAttribute(name));
    return Number.isFinite(v) ? v : dflt;
  };
  const paint = (el) => {
    const fill = el.getAttribute("fill");
    const stroke = el.getAttribute("stroke");
    const opacity = el.getAttribute("opacity");
    return {
      fill: fill && fill !== "none" ? fill : "none",
      stroke: stroke && stroke !== "none" ? stroke : undefined,
      lw: num(el, "stroke-width", 0.5),
      opacity: opacity != null ? Number(opacity) : undefined,
    };
  };

  for (const el of root.querySelectorAll("rect, line, circle, text")) {
    const tag = el.tagName.toLowerCase();
    if (tag === "rect") {
      fig.rect(num(el, "x"), num(el, "y"), num(el, "width"), num(el, "height"), paint(el));
    } else if (tag === "line") {
      const p = paint(el);
      fig.line(num(el, "x1"), num(el, "y1"), num(el, "x2"), num(el, "y2"), {
        stroke: p.stroke || "#000000", lw: p.lw, opacity: p.opacity,
        dash: el.getAttribute("stroke-dasharray") || undefined,
      });
    } else if (tag === "circle") {
      fig.circle(num(el, "cx"), num(el, "cy"), num(el, "r"), paint(el));
    } else if (tag === "text") {
      const anchorMap = { middle: "middle", end: "end", start: "start" };
      fig.text(num(el, "x"), num(el, "y"), el.textContent || "", {
        size: num(el, "font-size", 7),
        fill: el.getAttribute("fill") || "#000000",
        anchor: anchorMap[el.getAttribute("text-anchor")] || "start",
        weight: el.getAttribute("font-weight") === "bold" ? "bold" : "normal",
        italic: el.getAttribute("font-style") === "italic",
      });
    }
  }
  return fig;
}

/** Vector PDF for a figure that is only available as a self-contained SVG string. */
export function downloadSvgAsPdf(svgText, stem) {
  downloadFigurePdf(svgToFigure(svgText), stem);
}

/**
 * Build the standard export menu entries for a figure.
 *
 * `build` is called lazily on click, so the exported figure always reflects the
 * current selection and filter state.
 */
export function figureExportItems(build, stem, { tsv } = {}) {
  const items = [
    ["Main figure — SVG (vector)", () => downloadFigureSvg(build(), stem)],
    ["Main figure — PDF (vector)", () => downloadFigurePdf(build(), stem)],
    [`Main figure — PNG (${EXPORT_DPI} dpi)`, () => downloadFigurePng(build(), stem)],
  ];
  if (tsv) items.push(["Source table — TSV", () => downloadFigureTsv(tsv(), stem)]);
  return items;
}
