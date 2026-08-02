// Declarative figure specification with two vector backends.
//
// A figure is built once from canonical coordinate data using a small set of
// primitives, then rendered either as a standalone SVG document or as a real
// vector PDF. PNG is rasterised from the SVG at publication resolution.
//
// This replaces the previous export path, which drew the live React SVG into a
// canvas and wrapped the resulting JPEG in a PDF page. That produced a
// single-image page sized from canvas pixels (a 4000-point-wide "72 dpi" page),
// with no selectable text, and — because the interactive marks take their colour
// from the application stylesheet — solid black blocks wherever a `fill` came
// only from CSS.
//
// Invariants this module guarantees for every exported figure:
//   * every mark carries literal fill / stroke / font attributes
//   * no CSS custom properties, no external stylesheet, no foreignObject
//   * an explicit white paper background
//   * PDF text is emitted as real text objects in a base-14 font, so it stays
//     selectable and searchable without an embedding step
//   * PDF page geometry comes from a physical preset, not from pixel counts
//
// The module is intentionally dependency-free so figures can also be rendered and
// validated headlessly in Node, outside any bundler.

// --------------------------------------------------------------------------- //
// Publication layout presets (physical geometry and typography)
// --------------------------------------------------------------------------- //

// Widths follow common journal column measures. Font sizes are points at the
// preset's physical width; MIN_FONT_PT is the smallest type we allow to be
// exported, so a figure can never ship with illegible labels.
export const MIN_FONT_PT = 5.5;

export const PRESETS = {
  // single-column figure, ~89 mm
  compact: {
    id: "compact",
    widthPt: 252,
    minHeightPt: 90,
    maxHeightPt: 620,
    font: { title: 8.5, subtitle: 6.8, label: 6.2, tick: 5.8, legend: 6, small: 5.5 },
    lw: { thin: 0.4, rule: 0.6, outline: 0.9 },
    marker: 2.2,
    margin: { top: 8, right: 8, bottom: 8, left: 8 },
    legend: "stacked",
  },
  // double-column figure, ~183 mm
  double: {
    id: "double",
    widthPt: 522,
    minHeightPt: 120,
    maxHeightPt: 720,
    font: { title: 10.5, subtitle: 8, label: 7.4, tick: 6.8, legend: 7, small: 6.2 },
    lw: { thin: 0.5, rule: 0.7, outline: 1.1 },
    marker: 2.8,
    margin: { top: 10, right: 12, bottom: 10, left: 12 },
    legend: "inline",
  },
  // full page width in landscape, for wide coordinate tracks
  full: {
    id: "full",
    widthPt: 756,
    minHeightPt: 140,
    maxHeightPt: 900,
    font: { title: 12, subtitle: 9, label: 8.2, tick: 7.4, legend: 7.6, small: 6.6 },
    lw: { thin: 0.6, rule: 0.8, outline: 1.3 },
    marker: 3.2,
    margin: { top: 12, right: 14, bottom: 12, left: 14 },
    legend: "inline",
  },
};

export function preset(name) {
  return PRESETS[name] || PRESETS.double;
}

// Restrained scientific palette, colour-blind safe (Okabe-Ito derived). Shared by
// the interactive views and the exported figures so the two cannot drift apart.
export const PALETTE = {
  paper: "#ffffff",
  ink: "#1c2433",
  muted: "#6B7280",
  axis: "#9aa4b2",
  grid: "#e3e8f0",
  exon: "#A9BED4",
  exonEdge: "#7B93AE",
  exonAlt: "#E69F00",
  exonPrimary: "#7FA0C0",
  domain: "#0072B2",
  domainAlt: "#56B4E9",
  family: "#B9C0C8",
  tm: "#CC79A7",
  candidate: "#F0C070",
  candidateEdge: "#B8801A",
  boundary: "#4A5568",
  identity: "#117733",
  // boundary classes
  exact_domain_edge: "#005B8F",
  near_domain_edge: "#0072B2",
  inside_domain: "#009E73",
  outside_annotated_domains: "#D55E00",
  unavailable_or_uncertain: "#9AA0A6",
};

// A base-14 PDF font needs no embedding and keeps text searchable. The SVG side
// names the matching real families so both look the same.
const PDF_FONTS = { normal: "Helvetica", bold: "Helvetica-Bold", italic: "Helvetica-Oblique" };
export const SVG_FONT = "Helvetica, Arial, 'Nimbus Sans', 'DejaVu Sans', sans-serif";

// --------------------------------------------------------------------------- //
// Figure builder
// --------------------------------------------------------------------------- //

const n = (v) => {
  const x = Number(v);
  return Number.isFinite(x) ? Math.round(x * 100) / 100 : 0;
};

const escXml = (s) => String(s ?? "")
  .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
  .replace(/"/g, "&quot;");

// PDF strings are Latin-1 with escaped delimiters. Characters outside that range
// (en dashes, middle dots, Greek) are transliterated so the PDF never contains a
// broken glyph; the SVG keeps the original typography.
// Values are the WinAnsiEncoding code points the PDF text operators must carry:
// en dash is 0x96 and middle dot 0xB7 in WinAnsi, while arrows have no slot and
// are spelled out. The SVG keeps the original typography.
const PDF_TRANSLITERATE = {
  "–": "\u0096", "—": "\u0097", "·": "\u00b7", "→": "->", "←": "<-",
  "★": "*", "≈": "~", "≤": "<=", "≥": ">=", "±": "\u00b1",
  "⇥": "->", "◀": "<", "▶": ">", "＋": "+", "－": "-",
};

function pdfString(s) {
  let out = "";
  for (const ch of String(s ?? "")) {
    const mapped = PDF_TRANSLITERATE[ch] ?? ch;
    for (const c of mapped) {
      const code = c.codePointAt(0);
      if (c === "(" || c === ")" || c === "\\") out += `\\${c}`;
      else if (code < 32) out += " ";
      else if (code < 256) out += c;
      else out += "?";
    }
  }
  return out;
}

function hexToRgb(hex) {
  const h = String(hex || "").replace("#", "");
  if (h.length === 3) {
    return [parseInt(h[0] + h[0], 16) / 255, parseInt(h[1] + h[1], 16) / 255,
      parseInt(h[2] + h[2], 16) / 255];
  }
  if (h.length !== 6) return [0, 0, 0];
  return [parseInt(h.slice(0, 2), 16) / 255, parseInt(h.slice(2, 4), 16) / 255,
    parseInt(h.slice(4, 6), 16) / 255];
}

/**
 * Create a figure specification.
 *
 * @param {object} opts
 * @param {string} opts.preset  layout preset id: "compact" | "double" | "full"
 * @param {number} opts.height  figure height in points
 * @param {number} [opts.width] overrides the preset width (rarely needed)
 */
export function createFigure({ preset: presetName = "double", height, width } = {}) {
  const P = preset(presetName);
  const W = width ?? P.widthPt;
  const H = Math.max(P.minHeightPt, Math.min(P.maxHeightPt, height ?? P.minHeightPt));
  const marks = [];
  const warnings = [];

  const api = {
    preset: P,
    palette: PALETTE,
    width: W,
    height: H,
    get warnings() { return warnings; },

    /** Filled and/or stroked rectangle. */
    rect(x, y, w, h, { fill = "none", stroke, lw, opacity, rx } = {}) {
      marks.push({ t: "rect", x: n(x), y: n(y), w: n(Math.max(0.2, w)), h: n(h),
        fill, stroke, lw: lw ?? P.lw.thin, opacity, rx });
      return api;
    },

    line(x1, y1, x2, y2, { stroke = PALETTE.axis, lw, dash, opacity } = {}) {
      marks.push({ t: "line", x1: n(x1), y1: n(y1), x2: n(x2), y2: n(y2),
        stroke, lw: lw ?? P.lw.rule, dash, opacity });
      return api;
    },

    circle(cx, cy, r, { fill = PALETTE.ink, stroke, lw, opacity } = {}) {
      marks.push({ t: "circle", cx: n(cx), cy: n(cy), r: n(r), fill, stroke,
        lw: lw ?? P.lw.thin, opacity });
      return api;
    },

    /**
     * Text mark. `size` is a key of the preset font table or an explicit point
     * size; anything below the publication minimum is raised and recorded.
     */
    text(x, y, s, { size = "label", fill = PALETTE.ink, anchor = "start",
      weight = "normal", italic = false, opacity } = {}) {
      let pt = typeof size === "number" ? size : (P.font[size] ?? P.font.label);
      if (pt < MIN_FONT_PT) {
        warnings.push(`font size ${pt}pt raised to the ${MIN_FONT_PT}pt publication minimum`);
        pt = MIN_FONT_PT;
      }
      marks.push({ t: "text", x: n(x), y: n(y), s: String(s ?? ""), size: pt,
        fill, anchor, weight, italic, opacity });
      return api;
    },

    /**
     * Legend of swatch + label pairs; returns the y below the last entry.
     *
     * The inline layout wraps onto further lines instead of running past the
     * right margin, so a long entry can never be clipped.
     */
    legend(x, y, items, { size = "legend", swatch = 7, maxX } = {}) {
      const pt = typeof size === "number" ? size : (P.font[size] ?? P.font.legend);
      const right = maxX ?? (W - P.margin.right);
      if (P.legend === "stacked") {
        let ly = y;
        for (const [colour, label] of items) {
          api.rect(x, ly - swatch + 1, swatch, swatch, { fill: colour, stroke: PALETTE.axis });
          api.text(x + swatch + 4, ly, label, { size });
          ly += pt + 3.5;
        }
        return ly;
      }
      let lx = x;
      let ly = y;
      for (const [colour, label] of items) {
        const entryW = swatch + 4 + textWidth(label, pt);
        if (lx > x && lx + entryW > right) {
          lx = x;
          ly += pt + 4.5;
        }
        api.rect(lx, ly - swatch + 1, swatch, swatch, { fill: colour, stroke: PALETTE.axis });
        api.text(lx + swatch + 4, ly, label, { size });
        lx += entryW + 10;
      }
      return ly + pt + 4;
    },

    /** Marks in insertion order — used by tests and by the PDF/SVG backends. */
    get marks() { return marks; },

    /**
     * Trim (or extend) the canvas after drawing, so a figure carries no dead
     * whitespace. Clamped to the preset's height range.
     */
    resize(height) {
      api.height = Math.max(P.minHeightPt, Math.min(P.maxHeightPt, height));
      return api;
    },

    toSvg() { return renderSvg(api, marks); },
    toPdf() { return renderPdf(api, marks); },
  };
  return api;
}

/** Helvetica advance-width estimate, good enough for collision avoidance. */
export function textWidth(s, pt) {
  let units = 0;
  for (const ch of String(s ?? "")) {
    if ("iljI.,:;'|!".includes(ch)) units += 0.28;
    else if ("frt()[]-/ ".includes(ch)) units += 0.36;
    else if ("mwMW".includes(ch)) units += 0.86;
    else if (ch >= "A" && ch <= "Z") units += 0.68;
    else if (ch >= "0" && ch <= "9") units += 0.556;
    else units += 0.53;
  }
  return units * pt;
}

// --------------------------------------------------------------------------- //
// SVG backend
// --------------------------------------------------------------------------- //

function svgPaint(fill, stroke, lw, opacity) {
  let out = ` fill="${fill ?? "none"}"`;
  if (stroke) out += ` stroke="${stroke}" stroke-width="${n(lw)}"`;
  else out += ` stroke="none"`;
  if (opacity != null) out += ` opacity="${opacity}"`;
  return out;
}

function renderSvg(fig, marks) {
  const parts = [];
  // Explicit paper background: exported figures are viewed on many backgrounds.
  parts.push(`<rect x="0" y="0" width="${n(fig.width)}" height="${n(fig.height)}" `
    + `fill="${PALETTE.paper}" stroke="none" />`);
  for (const m of marks) {
    if (m.t === "rect") {
      parts.push(`<rect x="${m.x}" y="${m.y}" width="${m.w}" height="${m.h}"`
        + (m.rx ? ` rx="${n(m.rx)}"` : "")
        + svgPaint(m.fill, m.stroke, m.lw, m.opacity) + " />");
    } else if (m.t === "line") {
      parts.push(`<line x1="${m.x1}" y1="${m.y1}" x2="${m.x2}" y2="${m.y2}" `
        + `stroke="${m.stroke}" stroke-width="${n(m.lw)}"`
        + (m.dash ? ` stroke-dasharray="${m.dash}"` : "")
        + (m.opacity != null ? ` opacity="${m.opacity}"` : "") + " />");
    } else if (m.t === "circle") {
      parts.push(`<circle cx="${m.cx}" cy="${m.cy}" r="${m.r}"`
        + svgPaint(m.fill, m.stroke, m.lw, m.opacity) + " />");
    } else if (m.t === "text") {
      const anchor = m.anchor === "middle" ? "middle" : m.anchor === "end" ? "end" : "start";
      parts.push(`<text x="${m.x}" y="${m.y}" font-family="${SVG_FONT}" `
        + `font-size="${n(m.size)}" fill="${m.fill}" text-anchor="${anchor}"`
        + (m.weight === "bold" ? ` font-weight="bold"` : ` font-weight="normal"`)
        + (m.italic ? ` font-style="italic"` : "")
        + (m.opacity != null ? ` opacity="${m.opacity}"` : "")
        + `>${escXml(m.s)}</text>`);
    }
  }
  return `<?xml version="1.0" encoding="UTF-8"?>\n`
    + `<svg xmlns="http://www.w3.org/2000/svg" version="1.1" `
    + `width="${n(fig.width)}" height="${n(fig.height)}" `
    + `viewBox="0 0 ${n(fig.width)} ${n(fig.height)}">\n`
    + parts.join("\n") + `\n</svg>\n`;
}

// --------------------------------------------------------------------------- //
// PDF backend — real vector content, real text objects
// --------------------------------------------------------------------------- //

// PDF space has its origin bottom-left and y growing upwards, while the figure
// specification uses screen coordinates. One flip at emit time keeps every
// builder in a single coordinate system.
function pdfPage(fig, marks) {
  const H = fig.height;
  const fy = (y) => n(H - y);
  const ops = [];
  const usedFonts = new Set();

  const setFill = (hex) => {
    const [r, g, b] = hexToRgb(hex);
    ops.push(`${n(r)} ${n(g)} ${n(b)} rg`);
  };
  const setStroke = (hex) => {
    const [r, g, b] = hexToRgb(hex);
    ops.push(`${n(r)} ${n(g)} ${n(b)} RG`);
  };

  // Opaque paper.
  setFill(PALETTE.paper);
  ops.push(`0 0 ${n(fig.width)} ${n(H)} re f`);

  for (const m of marks) {
    const alpha = m.opacity != null ? Number(m.opacity) : 1;
    // A minimal PDF has no soft-mask resources, so translucency is approximated
    // by blending the mark towards the paper. Geometry stays vector either way.
    const blend = (hex) => {
      if (alpha >= 0.999) return hex;
      const [r, g, b] = hexToRgb(hex);
      const mix = (c) => Math.round((c * alpha + 1 * (1 - alpha)) * 255);
      return `#${[mix(r), mix(g), mix(b)].map((v) => v.toString(16).padStart(2, "0")).join("")}`;
    };

    if (m.t === "rect") {
      const hasFill = m.fill && m.fill !== "none";
      if (hasFill) setFill(blend(m.fill));
      if (m.stroke) { setStroke(blend(m.stroke)); ops.push(`${n(m.lw)} w`); }
      ops.push(`${m.x} ${fy(m.y + m.h)} ${m.w} ${m.h} re`);
      ops.push(hasFill && m.stroke ? "B" : hasFill ? "f" : "S");
    } else if (m.t === "line") {
      setStroke(blend(m.stroke));
      ops.push(`${n(m.lw)} w`);
      if (m.dash) {
        const d = String(m.dash).trim().split(/[\s,]+/).map(n).join(" ");
        ops.push(`[${d}] 0 d`);
      } else {
        ops.push("[] 0 d");
      }
      ops.push(`${m.x1} ${fy(m.y1)} m ${m.x2} ${fy(m.y2)} l S`);
    } else if (m.t === "circle") {
      // Four Bézier arcs; k is the standard circle-to-cubic constant.
      const k = 0.5523 * m.r;
      const cx = m.cx, cy = fy(m.cy), r = m.r;
      const hasFill = m.fill && m.fill !== "none";
      if (hasFill) setFill(blend(m.fill));
      if (m.stroke) { setStroke(blend(m.stroke)); ops.push(`${n(m.lw)} w`); }
      ops.push(`${n(cx - r)} ${n(cy)} m`);
      ops.push(`${n(cx - r)} ${n(cy + k)} ${n(cx - k)} ${n(cy + r)} ${n(cx)} ${n(cy + r)} c`);
      ops.push(`${n(cx + k)} ${n(cy + r)} ${n(cx + r)} ${n(cy + k)} ${n(cx + r)} ${n(cy)} c`);
      ops.push(`${n(cx + r)} ${n(cy - k)} ${n(cx + k)} ${n(cy - r)} ${n(cx)} ${n(cy - r)} c`);
      ops.push(`${n(cx - k)} ${n(cy - r)} ${n(cx - r)} ${n(cy - k)} ${n(cx - r)} ${n(cy)} c`);
      ops.push(hasFill && m.stroke ? "b" : hasFill ? "f" : "s");
    } else if (m.t === "text") {
      const key = m.weight === "bold" ? "bold" : m.italic ? "italic" : "normal";
      usedFonts.add(key);
      const str = pdfString(m.s);
      if (!str) continue;
      // PDF has no text-anchor, so the advance width is applied here.
      let x = m.x;
      if (m.anchor === "middle") x = m.x - textWidth(m.s, m.size) / 2;
      else if (m.anchor === "end") x = m.x - textWidth(m.s, m.size);
      setFill(blend(m.fill));
      ops.push("BT");
      ops.push(`/F_${key} ${n(m.size)} Tf`);
      ops.push(`${n(x)} ${fy(m.y)} Td`);
      ops.push(`(${str}) Tj`);
      ops.push("ET");
    }
  }

  return { width: fig.width, height: H, content: ops.join("\n"), usedFonts };
}

function renderPdf(fig, marks) {
  return assemblePdf([pdfPage(fig, marks)]);
}

/**
 * One PDF holding one page per figure specification.
 *
 * A wrapped residue-level alignment does not fit on a single page, and splitting
 * it into separate files would make it unreadable as one document.
 */
export function renderPdfPages(figures) {
  const pages = figures.map((f) => pdfPage(f, f.marks));
  if (!pages.length) throw new Error("a PDF needs at least one page");
  return assemblePdf(pages);
}

function assemblePdf(pages) {
  const usedFonts = new Set();
  for (const p of pages) for (const k of p.usedFonts) usedFonts.add(k);
  const fonts = usedFonts.size ? [...usedFonts] : ["normal"];

  // Object layout: 1 catalog, 2 page tree, then per page a page object and its
  // content stream, then the shared font objects.
  const objects = [];
  const pageObjNum = (i) => 3 + i * 2;
  const contentObjNum = (i) => 4 + i * 2;
  const fontObjStart = 3 + pages.length * 2;
  const fontRefs = fonts
    .map((k, i) => `/F_${k} ${fontObjStart + i} 0 R`).join(" ");

  objects.push("<< /Type /Catalog /Pages 2 0 R >>");
  objects.push(`<< /Type /Pages /Kids [`
    + pages.map((_, i) => `${pageObjNum(i)} 0 R`).join(" ")
    + `] /Count ${pages.length} >>`);
  pages.forEach((p, i) => {
    objects.push(`<< /Type /Page /Parent 2 0 R `
      + `/MediaBox [0 0 ${n(p.width)} ${n(p.height)}] `
      + `/Resources << /Font << ${fontRefs} >> >> /Contents ${contentObjNum(i)} 0 R >>`);
    objects.push({ stream: p.content });
  });
  for (const k of fonts) {
    objects.push(`<< /Type /Font /Subtype /Type1 /BaseFont /${PDF_FONTS[k]} `
      + `/Encoding /WinAnsiEncoding >>`);
  }

  // A PDF byte stream is not UTF-8. Text shown with WinAnsiEncoding must be one
  // byte per character, so encoding via TextEncoder would turn a middle dot into
  // the two bytes that render as "Â·".
  const enc = (s) => {
    const out = new Uint8Array(s.length);
    for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i) & 0xff;
    return out;
  };
  const chunks = [];
  let length = 0;
  const push = (bytes) => { chunks.push(bytes); length += bytes.length; };

  push(enc("%PDF-1.4\n%\u00e2\u00e3\u00cf\u00d3\n"));
  const offsets = [0];
  objects.forEach((obj, i) => {
    offsets.push(length);
    if (typeof obj === "object" && obj.stream != null) {
      const body = enc(obj.stream);
      push(enc(`${i + 1} 0 obj\n<< /Length ${body.length} >>\nstream\n`));
      push(body);
      push(enc("\nendstream\nendobj\n"));
    } else {
      push(enc(`${i + 1} 0 obj\n${obj}\nendobj\n`));
    }
  });

  const xref = length;
  let table = `xref\n0 ${objects.length + 1}\n0000000000 65535 f \n`;
  for (let i = 1; i <= objects.length; i++) {
    table += `${String(offsets[i]).padStart(10, "0")} 00000 n \n`;
  }
  push(enc(table));
  push(enc(`trailer\n<< /Size ${objects.length + 1} /Root 1 0 R >>\n`
    + `startxref\n${xref}\n%%EOF\n`));

  const out = new Uint8Array(length);
  let o = 0;
  for (const c of chunks) { out.set(c, o); o += c.length; }
  return out;
}

// --------------------------------------------------------------------------- //
// Axis helper — readable major/minor ticks for a coordinate range
// --------------------------------------------------------------------------- //

/**
 * Major and minor tick positions for [lo, hi] at a round step.
 * Returns { major, minor, step } with values in data units.
 */
export function axisTicks(lo, hi, targetMajor = 8) {
  const span = Math.max(1, hi - lo);
  const raw = span / Math.max(1, targetMajor);
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) || mag * 10;
  const minorStep = step / (step / mag === 2.5 ? 2.5 : 2);
  const major = [];
  const minor = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-9; t += step) major.push(Math.round(t));
  for (let t = Math.ceil(lo / minorStep) * minorStep; t <= hi + 1e-9; t += minorStep) {
    const v = Math.round(t);
    if (!major.includes(v)) minor.push(v);
  }
  return { major, minor, step };
}

/**
 * Place labels for a series of blocks without overlap.
 *
 * Returns one entry per block: either an inside label (the block is wide enough),
 * or a label on one of two alternating rows below the track with a leader line,
 * or `null` when even the alternating rows would collide — in which case the
 * caller should fall back to a side key rather than draw unreadable text.
 */
export function placeBlockLabels(blocks, { size, minInsidePad = 2, rows = 2 } = {}) {
  const out = [];
  const rowEnd = new Array(rows).fill(-Infinity);
  for (const b of blocks) {
    const w = textWidth(b.label, size);
    if (b.x1 - b.x0 >= w + minInsidePad * 2) {
      out.push({ ...b, mode: "inside", labelX: (b.x0 + b.x1) / 2, width: w });
      continue;
    }
    const cx = (b.x0 + b.x1) / 2;
    const start = cx - w / 2;
    let placed = false;
    for (let r = 0; r < rows; r++) {
      if (start > rowEnd[r] + 2) {
        rowEnd[r] = start + w;
        out.push({ ...b, mode: "below", row: r, labelX: cx, width: w });
        placed = true;
        break;
      }
    }
    if (!placed) out.push({ ...b, mode: "none", width: w });
  }
  return out;
}
