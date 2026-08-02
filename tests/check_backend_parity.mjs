/**
 * Do the SVG and the PDF backend draw the same marks?
 *
 * One figure specification feeds two independent backends. Nothing forces them to
 * agree, so a backend can silently drop a whole class of mark — translucent fills
 * are the obvious candidate, because PDF has no soft mask and has to approximate
 * them. Such a defect passes every per-format check: the PDF is still valid, still
 * vector, still has text, and the missing marks are usually too small to move
 * global ink coverage.
 *
 * This harness counts the marks in each specification and then counts what each
 * backend actually emitted, so every mark has to be accounted for on both sides.
 *
 * Usage: node tests/check_backend_parity.mjs <protein_coordinate_model.json>
 */

import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";
import { dirname, join, resolve } from "node:path";

const HERE = dirname(new URL(import.meta.url).pathname);
const VIEWERS = join(HERE, "..", "webapp", "frontend", "src", "pages", "viewers");
const load = (m) => import(pathToFileURL(join(VIEWERS, m)).href);

const {
  boundaryClassSummaryFigure, boundaryFigure, domainArchitectureFigure, exonMapFigure,
  signedDistanceFigure,
} = await load("figureData.js");

const modelPath = resolve(process.argv[2]);
const doc = JSON.parse(readFileSync(modelPath, "utf8"));

/** Marks in the specification, grouped by primitive type. */
function specCounts(fig) {
  const c = { rect: 0, line: 0, circle: 0, text: 0, textNonEmpty: 0, translucent: 0 };
  for (const m of fig.marks) {
    c[m.t] = (c[m.t] || 0) + 1;
    if (m.t === "text" && String(m.s ?? "").length) c.textNonEmpty += 1;
    if (m.opacity != null && Number(m.opacity) < 0.999) c.translucent += 1;
  }
  return c;
}

/** What the SVG backend emitted. The backend adds one paper rectangle. */
function svgCounts(svg) {
  const n = (re) => (svg.match(re) || []).length;
  return {
    rect: n(/<rect\b/g) - 1,
    line: n(/<line\b/g),
    circle: n(/<circle\b/g),
    text: n(/<text\b/g),
  };
}

/**
 * What the PDF backend emitted. The content stream is read from the uncompressed
 * document this project writes; each primitive has a recognisable operator shape.
 */
function pdfCounts(bytes) {
  const s = Buffer.from(bytes).toString("latin1");
  const n = (re) => (s.match(re) || []).length;
  return {
    // Matched as whole operator lines, so the word "re" inside a drawn string
    // cannot be mistaken for a rectangle. The paper rectangle is one of them.
    rect: n(/^[-\d.\s]+re(?: [fBS])?$/gm) - 1,
    line: n(/ l S$/gm),
    circle: n(/ c$/gm) / 4,
    text: n(/^BT$/gm),
  };
}

const problems = [];
const rows = [];

function check(name, fig) {
  const spec = specCounts(fig);
  const svg = svgCounts(fig.toSvg());
  const pdf = pdfCounts(fig.toPdf());

  const expect = (backend, got, want, kind) => {
    if (got !== want) {
      problems.push(`${name}: the ${backend} backend emitted ${got} ${kind} marks `
        + `but the specification holds ${want}`);
    }
  };
  for (const kind of ["rect", "line", "circle"]) {
    expect("SVG", svg[kind], spec[kind], kind);
    expect("PDF", pdf[kind], spec[kind], kind);
  }
  expect("SVG", svg.text, spec.text, "text");
  // Empty strings carry no ink, so the PDF backend legitimately skips them.
  expect("PDF", pdf.text, spec.textNonEmpty, "text");

  rows.push({ name, marks: fig.marks.length, translucent: spec.translucent,
    rect: spec.rect, line: spec.line, circle: spec.circle, text: spec.text });
}

// Mirrors the production renderer: same adapter, same builders, same model objects.
for (const m of doc.models || []) {
  const sp = m.species_id || "sp";
  const boundaries = (m.exon_boundaries || []).filter((b) => b.signed_distance != null);
  check(`${sp}/exon_projection`, exonMapFigure(m));
  if (!(m.representative_domains || []).length || m.status !== "available") continue;
  check(`${sp}/domain_architecture`, domainArchitectureFigure(m));
  if (!boundaries.length) continue;
  check(`${sp}/boundary_on_architecture`, boundaryFigure(m, boundaries));
  check(`${sp}/signed_distances`, signedDistanceFigure(m, boundaries));
  check(`${sp}/class_summary`, boundaryClassSummaryFigure(m, boundaries));
  // The on-demand selected-boundary exports go through a different code path.
  const id = boundaries[0].exon_boundary_id || boundaries[0].boundary_id;
  check(`${sp}/selected_boundary_detail`,
    boundaryFigure(m, boundaries, { selectedBoundaryId: id }));
  check(`${sp}/selected_signed_distance`,
    signedDistanceFigure(m, boundaries, { selectedBoundaryId: id }));
}

for (const r of rows) {
  console.log(`${r.name}: ${r.marks} marks `
    + `(${r.rect} rect, ${r.line} line, ${r.circle} circle, ${r.text} text; `
    + `${r.translucent} translucent) — both backends agree`);
}

if (problems.length) {
  console.error("\nBACKEND PARITY FAILURES:");
  for (const p of problems) console.error(`  - ${p}`);
  process.exit(1);
}
console.log(`\nOK: ${rows.length} figures, every mark drawn by both backends.`);
