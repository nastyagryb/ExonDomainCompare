// Publication figure builders for the single-species main figures.
//
// Each builder turns canonical coordinate-model data into a figure
// specification (see figureSpec.js), which renders to a standalone SVG, a real
// vector PDF and a 300 dpi PNG. The interactive React views and these builders
// read the same records, so coordinate ranges, feature order, labels, colours
// and legends cannot drift apart.
//
// Every figure answers exactly one scientific question, stated in its subtitle.

import {
  PALETTE, axisTicks, createFigure, placeBlockLabels, preset, textWidth,
  // The explicit extension keeps this module importable by plain Node, so the
  // figures can be rendered and validated without a browser or bundler.
} from "./figureSpec.js";
import { canonClass } from "./boundaryClasses.js";
import { candidateDisplayLayout, candidateLabelFits } from "./candidateDisplay.js";
import {
  BOUNDARY_CLASS_LABEL, FEATURE_STYLES, boundaryStyleKey, domainInstanceFill,
} from "./semanticStyles.js";

// --------------------------------------------------------------------------- //
// Domain instance identity
// --------------------------------------------------------------------------- //

import { prettyDomainName } from "./domainNames.js";

export { prettyDomainName };

/**
 * Normalise representative domains into instances with stable identity.
 *
 * Three FGFR1 Ig-like domains share the accession IPR007110 but are distinct
 * feature instances at different coordinates. Resolving a feature by accession
 * alone therefore collapses them onto each other, so every consumer works with
 * `domain_instance_id` = `<accession>:<start>-<end>` instead. A backend-supplied
 * id always wins; the derivation here keeps older data usable.
 */
export function domainInstances(domains) {
  const list = (domains || [])
    .map((d) => ({
      accession: d.interpro_accession || d.accession || d.domain_id || "",
      rawName: d.interpro_name || d.domain_name || d.label || d.name || "",
      start: Number(d.start_aa ?? d.start),
      end: Number(d.end_aa ?? d.end),
      featureType: d.interpro_type || d.feature_type || "DOMAIN",
      source: d.domain_source || d.source || "",
      memberSignatures: d.member_databases || d.member_signatures || [],
      // Identity and labels established by the coordinate model always win; the
      // derivation below only covers data written before that model existed.
      backendId: d.domain_instance_id || null,
      backendShort: d.short_label || null,
      backendFull: d.full_label || null,
      backendNumber: d.instance_number ?? null,
      backendCount: d.n_instances_of_accession ?? null,
    }))
    .filter((d) => Number.isFinite(d.start) && Number.isFinite(d.end))
    .sort((a, b) => a.start - b.start || a.end - b.end);

  // Instance numbers run in coordinate order within one accession.
  const counts = {};
  for (const d of list) counts[d.accession] = (counts[d.accession] || 0) + 1;
  const seen = {};

  return list.map((d, i) => {
    seen[d.accession] = (seen[d.accession] || 0) + 1;
    const instanceNumber = d.backendNumber ?? seen[d.accession];
    const instanceCount = d.backendCount ?? counts[d.accession];
    const base = prettyDomainName(d.rawName);
    // Only number a label when the accession really does repeat.
    const shortLabel = d.backendShort
      || (instanceCount > 1 ? `${base} ${instanceNumber}` : base);
    return {
      domain_instance_id: d.backendId || `${d.accession}:${d.start}-${d.end}`,
      interpro_accession: d.accession,
      short_label: shortLabel,
      full_label: d.backendFull || `${shortLabel} · aa ${d.start}–${d.end}`,
      instance_number: instanceNumber,
      instance_count: instanceCount,
      start: d.start,
      end: d.end,
      feature_type: d.featureType,
      source: d.source,
      member_signatures: d.memberSignatures,
      display_order: i + 1,
      raw_name: d.rawName,
    };
  });
}

/** Colour ramp for representative domain instances, stable by display order. */
export const domainColour = (inst) => domainInstanceFill(inst?.display_order);

// --------------------------------------------------------------------------- //
// Shared scaffolding
// --------------------------------------------------------------------------- //

/** Species names are set in italics in scientific figures. */
export function speciesLabel(species) {
  if (!species) return "";
  return String(species).replace(/_/g, " ")
    .replace(/^([a-z])/, (m) => m.toUpperCase());
}

export function headerBlock(fig, { title, species, subtitle, question, isoform = "" }) {
  const P = fig.preset;
  let y = P.margin.top + P.font.title;
  if (title) {
    // The isoform belongs in the title, not the subtitle. Where a gene's analysis
    // is about two mutually exclusive isoforms of one species, a reader comparing
    // two figures has to be able to tell them apart at a glance, and a protein
    // accession buried in the subtitle does not do that.
    const head = isoform ? `${title} ${isoform}` : title;
    fig.text(P.margin.left, y, head, { size: "title", weight: "bold" });
    if (species) {
      const w = textWidth(head, P.font.title);
      fig.text(P.margin.left + w + 4, y, `· ${speciesLabel(species)}`,
        { size: "title", italic: true, fill: PALETTE.ink });
    }
    y += P.font.subtitle + 4;
  }
  if (subtitle) {
    fig.text(P.margin.left, y, subtitle, { size: "subtitle", fill: PALETTE.muted });
    y += P.font.subtitle + 3;
  }
  if (question) {
    fig.text(P.margin.left, y, question, { size: "small", fill: PALETTE.muted, italic: true });
    y += P.font.small + 4;
  }
  return y;
}

/**
 * Amino-acid coordinate axis with a strong major grid and subtle minor ticks.
 * The first and last residue are always labelled explicitly.
 */
function aaAxis(fig, { x0, x1, y, lo, hi, scale, label }) {
  const P = fig.preset;
  const { major, minor } = axisTicks(lo, hi, Math.max(5, Math.round((x1 - x0) / 70)));
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of minor) {
    if (t < lo || t > hi) continue;
    fig.line(scale(t), y, scale(t), y + 2, { stroke: PALETTE.axis, lw: P.lw.thin, opacity: 0.7 });
  }
  const drawn = [];
  const place = (value, text) => {
    const x = scale(value);
    const w = textWidth(text, P.font.tick);
    if (drawn.some(([a, b]) => x - w / 2 < b + 3 && x + w / 2 > a - 3)) return;
    drawn.push([x - w / 2, x + w / 2]);
    fig.line(x, y, x, y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(x, y + 3.5 + P.font.tick + 1, text,
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  };
  // Termini first so a crowded interior tick yields to them.
  place(lo, String(lo));
  place(hi, String(hi));
  for (const t of major) {
    if (t <= lo || t >= hi) continue;
    place(t, String(t));
  }
  if (label) {
    fig.text((x0 + x1) / 2, y + 3.5 + P.font.tick * 2 + 6, label,
      { size: "label", anchor: "middle", fill: PALETTE.ink });
  }
  return y + 3.5 + P.font.tick * 2 + (label ? 8 : 2);
}

/** Left-hand track label in a reserved gutter, so it can never overlap a mark. */
function trackLabel(fig, x, yCentre, text, { bold = false } = {}) {
  fig.text(x, yCentre + fig.preset.font.label * 0.35, text,
    { size: "label", anchor: "end", fill: PALETTE.ink, weight: bold ? "bold" : "normal" });
}

/**
 * Exon blocks with readable identifiers.
 *
 * Labels go inside the block when it is wide enough, otherwise onto one of two
 * alternating rows underneath with a leader line. Dark text is never placed on a
 * dark fill: the in-block colour is chosen from the block's own fill.
 */
function exonTrack(fig, { exons, y, h, scale, palette = PALETTE.exon,
  edge = PALETTE.exonEdge, showLabels = true, labelRows = 2, highlight = null,
  eventExon = null }) {
  const P = fig.preset;
  const blocks = exons.map((e) => ({
    x0: scale(e.start), x1: scale(e.end + 1), label: e.label, exon: e,
  }));
  for (const b of blocks) {
    const isHi = highlight && highlight(b.exon);
    const isEvent = !isHi && eventExon && eventExon(b.exon);
    const ev = FEATURE_STYLES.validated_event;
    fig.rect(b.x0, y, Math.max(0.6, b.x1 - b.x0), h, {
      fill: isHi ? PALETTE.exonPrimary : (isEvent ? ev.fill : palette),
      stroke: isHi ? PALETTE.ink : (isEvent ? ev.stroke : edge),
      lw: isHi ? P.lw.outline : (isEvent ? ev.strokeWidth : P.lw.thin),
    });
  }
  if (!showLabels) return { bottom: y + h, placed: [] };

  const placed = placeBlockLabels(blocks, { size: P.font.small, rows: labelRows });
  let bottom = y + h;
  const rowY = (r) => y + h + 3 + P.font.small + r * (P.font.small + 2.5);
  for (const p of placed) {
    if (p.mode === "inside") {
      // Light fill, so dark ink stays legible inside the block.
      fig.text(p.labelX, y + h / 2 + P.font.small * 0.36, p.label,
        { size: "small", anchor: "middle", fill: PALETTE.ink });
    } else if (p.mode === "below") {
      const ly = rowY(p.row);
      fig.text(p.labelX, ly, p.label,
        { size: "small", anchor: "middle", fill: PALETTE.muted });
      fig.line((p.x0 + p.x1) / 2, y + h, p.labelX, ly - P.font.small,
        { stroke: PALETTE.axis, lw: P.lw.thin, opacity: 0.75 });
      bottom = Math.max(bottom, ly + 2);
    }
    // mode "none": the block is too narrow even for an offset row; the caller
    // provides a side key rather than drawing unreadable text.
  }
  return { bottom, placed };
}

// The boundary vocabulary is canonicalised through the shared module, because the
// coordinate model and the boundary index spell the same classes differently
// (`near_domain_edge` vs `near_edge`). Resolving both to one canonical key keeps a
// raw database string from ever reaching a figure label.
export const BOUNDARY_CLASS_ORDER = [
  "exact_domain_edge", "near_domain_edge", "inside_domain",
  "outside_annotated_domains", "unavailable_or_uncertain",
];

const CLASS_LABELS = BOUNDARY_CLASS_LABEL;

/** Canonical class of a boundary record, whichever vocabulary it was written in. */
export const boundaryClassOf = (b) =>
  canonClass(b?.boundary_class || b?.classification || b?.category);

export const boundaryClassColour = (c) =>
  FEATURE_STYLES[boundaryStyleKey(canonClass(c))].fill;
export const boundaryClassLabel = (c) => CLASS_LABELS[canonClass(c)];

/** Whether an exon carries an established, validated alternative-splicing event. */
export function isEventExon(exon) {
  return Boolean(exon?.is_cassette_exon || exon?.is_event_exon);
}

/**
 * The validated-event track: an established event, on the same axis, in its own row.
 *
 * Only a real, externally established event belongs here. An exploratory candidate
 * has its own track and its own style, and the two must never be drawn alike — the
 * distinction between "this is known" and "this is a positional guess" is the most
 * consequential thing a reader takes from these figures.
 */
function eventTrack(fig, { events, y, h, scale, x1 }) {
  const P = fig.preset;
  const style = FEATURE_STYLES.validated_event;
  for (const e of events) {
    const s = Number(e.start ?? e.start_aa);
    const t = Number(e.end ?? e.end_aa);
    if (!Number.isFinite(s) || !Number.isFinite(t)) continue;
    const ex0 = scale(s);
    const ew = Math.max(1.5, scale(t + 1) - ex0);
    fig.rect(ex0, y, ew, h, {
      fill: style.fill, stroke: style.stroke, lw: style.strokeWidth,
    });
    const label = `${e.label || e.cassette_id || "event"} · aa ${s}–${t}`;
    const lx = Math.min(ex0 + ew + 4, x1 - textWidth(label, P.font.small));
    fig.text(lx, y + h / 2 + P.font.small * 0.36, label,
      { size: "small", fill: PALETTE.muted });
  }
  return y + h;
}

// --------------------------------------------------------------------------- //
// Figure 1 — Exon map / primary exon-to-protein projection
// --------------------------------------------------------------------------- //

/**
 * Which coding exons produce which regions of the selected primary protein?
 *
 * Deliberately carries no domain track: the exon-to-protein projection is the
 * single question this figure answers.
 */
export function exonMapFigureSpec({
  gene, species, proteinId, transcriptId, proteinLength,
  exons = [], candidates = [], selectedCandidateId = null, presetName = "full",
  isoform = "", events = [],
}) {
  const P = preset(presetName);
  const gutter = 62;
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 10;
  const len = Math.max(1, proteinLength);
  const scale = (aa) => x0 + ((Math.min(Math.max(aa, 1), len + 1) - 1) / len) * (x1 - x0);

  const coding = exons.filter((e) => (e.feature_type || "coding_exon") === "coding_exon");
  const nCand = candidates.length;
  const fig = createFigure({ preset: presetName, height: 300 });

  let y = headerBlock(fig, {
    title: gene, species, isoform,
    subtitle: `Primary protein ${proteinId} · transcript ${transcriptId} · ${len} aa · `
      + `${coding.length} coding exons`,
    question: "Which coding exons produce which regions of the primary protein?",
  });
  y += 6;

  // --- exon track -----------------------------------------------------------
  // A validated event exon is drawn in the exon series it belongs to, not beside
  // it. FGFR2's cassette *is* one of the coding exons, and lifting it out would
  // hide the very context — all the other exons — the reader needs.
  const eventExonIds = new Set(exons.filter(isEventExon).map((e) => e.id));
  const exonH = 15;
  trackLabel(fig, x0 - 8, y + exonH / 2, "Coding exons");
  const { bottom: exonBottom, placed } = exonTrack(fig, {
    exons: coding, y, h: exonH, scale,
    eventExon: (e) => eventExonIds.has(e.id),
  });
  y = exonBottom + 4;

  // --- internal boundary markers -------------------------------------------
  const boundaries = coding.slice(0, -1).map((e, i) => ({
    position: e.end, from: e.label, to: coding[i + 1].label,
  }));
  const bndY = y;
  const bndH = 5;
  trackLabel(fig, x0 - 8, bndY + bndH / 2, "Exon boundaries");
  for (const b of boundaries) {
    fig.line(scale(b.position + 1), bndY, scale(b.position + 1), bndY + bndH,
      { stroke: PALETTE.boundary, lw: P.lw.rule });
  }
  y = bndY + bndH + 8;

  // --- validated event track ------------------------------------------------
  if (events.length) {
    const evH = 8;
    trackLabel(fig, x0 - 8, y + evH / 2, "Validated event");
    y = eventTrack(fig, { events, y, h: evH, scale, x1 }) + 8;
  }

  // --- candidate track ------------------------------------------------------
  if (nCand) {
    const candH = 8;
    trackLabel(fig, x0 - 8, y + candH / 2, "Candidates");
    const selected = candidates.find((c) => candidateId(c) === selectedCandidateId)
      || candidates[0];
    for (const c of candidates) {
      const isSel = candidateId(c) === candidateId(selected);
      const cx0 = scale(c.aa_start);
      const cw = Math.max(1, scale(c.aa_end + 1) - cx0);
      // Non-selected candidates stay visible but clearly subordinate.
      fig.rect(cx0, y, cw, candH, {
        fill: PALETTE.candidate,
        stroke: isSel ? PALETTE.candidateEdge : "none",
        lw: isSel ? P.lw.outline : 0,
        opacity: isSel ? 1 : 0.35,
      });
    }
    // The candidate identifier and interval sit outside the blocks, so no long
    // text is ever drawn into the plot area.
    const sLabel = `${candidateLabel(selected)} · aa ${selected.aa_start}–${selected.aa_end}`;
    const sx = scale((selected.aa_start + selected.aa_end) / 2);
    const lw = textWidth(sLabel, P.font.small);
    const lx = Math.min(Math.max(sx, x0 + lw / 2), x1 - lw / 2);
    y += candH + 3 + P.font.small;
    fig.text(lx, y, sLabel, { size: "small", anchor: "middle", fill: PALETTE.candidateEdge });
    y += 6;
  }

  // --- axis and legend ------------------------------------------------------
  y = aaAxis(fig, { x0, x1, y, lo: 1, hi: len, scale, label: "Primary protein position (aa)" });
  y += 2;

  const legend = [[PALETTE.exon, "Coding exon"], [PALETTE.boundary, "Internal exon boundary"]];
  if (events.length) {
    legend.push([FEATURE_STYLES.validated_event.fill,
      `${events[0].event_label || "Validated event exon"} (established, not a candidate)`]);
  }
  if (nCand) legend.push([PALETTE.candidate, "Exploratory candidate region (not validated)"]);
  y = fig.legend(x0, y + P.font.legend, legend);

  // A side key rescues exons too narrow to label in place.
  const unlabelled = placed.filter((p) => p.mode === "none").map((p) => p.label);
  if (unlabelled.length) {
    fig.text(x0, y + P.font.small + 1,
      `Exons too narrow to label in place: ${unlabelled.join(", ")}`,
      { size: "small", fill: PALETTE.muted });
    y += P.font.small + 4;
  }

  return finalise(fig, y);
}

// --------------------------------------------------------------------------- //
// Figure 2 — Integrated domain architecture
// --------------------------------------------------------------------------- //

/**
 * How do annotated protein features, the coding-exon structure and the
 * exploratory candidates align along the primary protein?
 *
 * Track order is fixed: representative domains, family/superfamily, membrane
 * topology, coding exons, exon boundaries, candidates.
 */
export function domainArchitectureFigureSpec({
  gene, species, proteinId, transcriptId, proteinLength,
  domains = [], families = [], tm = [], exons = [], candidates = [],
  selectedDomainInstanceId = null, presetName = "full",
  isoform = "", events = [], showAllCandidates = false,
}) {
  const P = preset(presetName);
  const gutter = 96;
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 10;
  const len = Math.max(1, proteinLength);
  const scale = (aa) => x0 + ((Math.min(Math.max(aa, 1), len + 1) - 1) / len) * (x1 - x0);

  const instances = domainInstances(domains);
  const coding = exons.filter((e) => (e.feature_type || "coding_exon") === "coding_exon");
  const fig = createFigure({ preset: presetName, height: 420 });

  let y = headerBlock(fig, {
    title: gene, species, isoform,
    subtitle: `Primary protein ${proteinId} · transcript ${transcriptId} · ${len} aa · `
      + `${instances.length} representative domain instances`,
    question: "How do annotated domains, membrane topology and coding exons align "
      + "along the primary protein?",
  });
  y += 8;
  const eventExonIds = new Set(exons.filter(isEventExon).map((e) => e.id));

  const GAP = 13;   // vertical separation between tracks
  const H = 14;

  // --- 1. representative domain instances -----------------------------------
  trackLabel(fig, x0 - 10, y + H / 2, "Representative domains", { bold: true });
  for (const d of instances) {
    const dx0 = scale(d.start);
    const dw = Math.max(1.2, scale(d.end + 1) - dx0);
    const isSel = selectedDomainInstanceId && d.domain_instance_id === selectedDomainInstanceId;
    fig.rect(dx0, y, dw, H, {
      fill: domainColour(d),
      stroke: isSel ? PALETTE.ink : "#ffffff",
      // Selection is an outline, never a translucent rectangle across tracks.
      lw: isSel ? P.lw.outline * 1.6 : P.lw.thin,
    });
    // A repeated instance is numbered inside the block; the full name is in the
    // legend, so no long text is drawn into the plot. When the label already ends
    // in a number the digit is omitted: a block reading "1" under a legend entry
    // reading "Ig-like domain 2" invites exactly the wrong reading, and that
    // happens whenever the display numbering and the per-accession instance
    // numbering differ (FGFR2's Ig1 and Ig3 share one accession, Ig2 does not).
    const labelled = /\d+$/.test(String(d.short_label || ""));
    const inner = (d.instance_count > 1 && !labelled) ? String(d.instance_number) : "";
    if (inner && dw >= textWidth(inner, P.font.small) + 5) {
      fig.text(dx0 + dw / 2, y + H / 2 + P.font.small * 0.36, inner,
        { size: "small", anchor: "middle", fill: "#ffffff", weight: "bold" });
    }
  }
  y += H + GAP;

  // --- 2. family / superfamily ---------------------------------------------
  if (families.length) {
    const fH = 8;
    trackLabel(fig, x0 - 10, y + fH / 2, "Family / superfamily");
    // Neutral grey and a thinner row: this is not another structural domain.
    const rows = packIntervals(families.map((f) => ({
      start: Number(f.start_aa ?? f.start), end: Number(f.end_aa ?? f.end),
      label: prettyDomainName(f.interpro_name || f.domain_name),
    })));
    for (const f of rows) {
      const fx0 = scale(f.start);
      fig.rect(fx0, y + f.row * (fH + 2.5), Math.max(1, scale(f.end + 1) - fx0), fH,
        { fill: PALETTE.family, stroke: "#ffffff", lw: P.lw.thin });
    }
    const nRows = Math.max(...rows.map((r) => r.row)) + 1;
    y += nRows * (fH + 2.5) + GAP - 2;
  }

  // --- 3. membrane topology -------------------------------------------------
  if (tm.length) {
    const tH = 9;
    trackLabel(fig, x0 - 10, y + tH / 2, "Membrane topology");
    for (const t of tm) {
      const s = Number(t.start_aa ?? t.start);
      const e = Number(t.end_aa ?? t.end);
      if (!Number.isFinite(s) || !Number.isFinite(e)) continue;
      const tx0 = scale(s);
      const tw = Math.max(1.5, scale(e + 1) - tx0);
      fig.rect(tx0, y, tw, tH, { fill: PALETTE.tm, stroke: PALETTE.ink, lw: P.lw.thin });
      const lbl = `TM helix · aa ${s}–${e}`;
      const lx = Math.min(tx0 + tw + 4, x1 - textWidth(lbl, P.font.small));
      fig.text(lx, y + tH / 2 + P.font.small * 0.36, lbl,
        { size: "small", fill: PALETTE.muted });
    }
    y += tH + GAP;
  }

  // --- 4. coding exons ------------------------------------------------------
  trackLabel(fig, x0 - 10, y + H / 2, "Coding exons");
  const { bottom } = exonTrack(fig, {
    exons: coding, y, h: H - 2, scale, labelRows: 2,
    eventExon: (e) => eventExonIds.has(e.id),
  });
  y = bottom + 6;

  // --- 5. exon boundaries ---------------------------------------------------
  const bH = 5;
  trackLabel(fig, x0 - 10, y + bH / 2, "Exon boundaries");
  for (const e of coding.slice(0, -1)) {
    fig.line(scale(e.end + 1), y, scale(e.end + 1), y + bH,
      { stroke: PALETTE.boundary, lw: P.lw.thin, opacity: 0.85 });
  }
  y += bH + GAP - 3;

  // --- 6. validated event ---------------------------------------------------
  if (events.length) {
    const evH = 8;
    trackLabel(fig, x0 - 10, y + evH / 2, "Validated event");
    y = eventTrack(fig, { events, y, h: evH, scale, x1 }) + GAP - 3;
  }

  // --- 7. candidates --------------------------------------------------------
  if (candidates.length) {
    // Same display clusters and lane assignment as the interactive architecture
    // track: overlapping candidates are packed, never drawn on top of each other.
    const layout = candidateDisplayLayout(candidates,
      { selectedId: null, showAll: showAllCandidates });
    const cH = 8, cGap = 2;
    trackLabel(fig, x0 - 10, y + cH / 2, "Candidate regions");
    for (const c of layout.visible) {
      const cx0 = scale(c.aa_start);
      const w = Math.max(1, scale(c.aa_end + 1) - cx0);
      const ly = y + layout.laneOf(c) * (cH + cGap);
      fig.rect(cx0, ly, w, cH,
        { fill: PALETTE.candidate, stroke: PALETTE.candidateEdge, lw: P.lw.thin, opacity: 0.8 });
      if (candidateLabelFits(w)) {
        fig.text(cx0 + w / 2, ly + cH - 2, candidateLabel(c),
          { size: "small", anchor: "middle", fill: PALETTE.candidateEdge });
      }
    }
    y += layout.laneCount * (cH + cGap) + 3 + P.font.small;
    const hidden = layout.hiddenCount
      ? ` · ${layout.hiddenCount} lower-ranked cluster(s) not shown` : "";
    fig.text(x0, y, `Exploratory candidates · ${layout.total} display cluster(s) in `
      + `${layout.laneCount} lane(s)${hidden} · biological validation: not validated`,
      { size: "small", fill: PALETTE.muted, italic: true });
    y += 6;
  }

  // --- axis and legend ------------------------------------------------------
  y = aaAxis(fig, { x0, x1, y, lo: 1, hi: len, scale, label: "Primary protein position (aa)" });

  const legend = instances.map((d) => [domainColour(d), d.full_label]);
  if (families.length) legend.push([PALETTE.family, "Family / homologous superfamily"]);
  if (tm.length) legend.push([PALETTE.tm, "Predicted TM helix (pyTMHMM)"]);
  legend.push([PALETTE.exon, "Coding exon"]);
  if (events.length) {
    legend.push([FEATURE_STYLES.validated_event.fill,
      `${events[0].event_label || "Validated event exon"} (established, not a candidate)`]);
  }
  if (candidates.length) legend.push([PALETTE.candidate, "Exploratory candidate"]);
  y = fig.legend(x0, y + P.font.legend + 2, legend, { size: "legend" });

  return finalise(fig, y);
}

/** Greedy interval packing so overlapping features never occupy the same row. */
function packIntervals(items) {
  const sorted = [...items].sort((a, b) => a.start - b.start || b.end - a.end);
  const rowEnds = [];
  return sorted.map((it) => {
    let row = rowEnds.findIndex((end) => it.start > end);
    if (row === -1) { row = rowEnds.length; rowEnds.push(it.end); } else { rowEnds[row] = it.end; }
    return { ...it, row };
  });
}

// --------------------------------------------------------------------------- //
// Figure 3 — Boundaries on the domain architecture
// --------------------------------------------------------------------------- //

/**
 * Where do internal coding-exon boundaries fall relative to annotated domains?
 *
 * Interpretation text, distances and identifiers belong to the subtitle, the
 * side annotation and the caption — never inside the architecture.
 */
export function boundaryFigureSpec({
  gene, species, proteinId, proteinLength, domains = [], exons = [],
  boundaries = [], candidates = [], selectedBoundaryId = null,
  nearEdgeThreshold = 5, showCandidates = true, presetName = "full",
  isoform = "", events = [], selectedCandidateId = null, showAllCandidates = false,
}) {
  const P = preset(presetName);
  const gutter = 92;
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 10;
  const len = Math.max(1, proteinLength);
  const scale = (aa) => x0 + ((Math.min(Math.max(aa, 1), len + 1) - 1) / len) * (x1 - x0);

  const instances = domainInstances(domains);
  const coding = exons.filter((e) => (e.feature_type || "coding_exon") === "coding_exon");
  const selected = boundaries.find((b) => boundaryId(b) === selectedBoundaryId) || null;
  const fig = createFigure({ preset: presetName, height: 340 });

  const counts = {};
  for (const b of boundaries) {
    const c = boundaryClassOf(b);
    counts[c] = (counts[c] || 0) + 1;
  }
  const countText = BOUNDARY_CLASS_ORDER.filter((c) => counts[c])
    .map((c) => [c, counts[c]])
    .map(([c, n]) => `${n} ${boundaryClassLabel(c).toLowerCase()}`).join(" · ");

  let y = headerBlock(fig, {
    title: gene, species, isoform,
    subtitle: `${proteinId} · ${boundaries.length} internal coding-exon boundaries · ${countText}`,
    question: "Where do internal coding-exon boundaries fall relative to "
      + "annotated domain edges?",
  });
  // Selected-boundary facts go here, as a concise annotation outside the plot.
  if (selected) {
    fig.text(P.margin.left, y, selectedAnnotation(selected),
      { size: "small", fill: PALETTE.ink });
    y += P.font.small + 5;
  }
  y += 5;

  // --- domains --------------------------------------------------------------
  const H = 14;
  trackLabel(fig, x0 - 10, y + H / 2, "Representative domains", { bold: true });
  for (const d of instances) {
    const dx0 = scale(d.start);
    const dw = Math.max(1.2, scale(d.end + 1) - dx0);
    fig.rect(dx0, y, dw, H, { fill: domainColour(d), stroke: "#ffffff", lw: P.lw.thin });
    const inner = (d.instance_count > 1 && !/\d+$/.test(String(d.short_label || "")))
      ? String(d.instance_number) : "";
    if (inner && dw >= textWidth(inner, P.font.small) + 5) {
      fig.text(dx0 + dw / 2, y + H / 2 + P.font.small * 0.36, inner,
        { size: "small", anchor: "middle", fill: "#ffffff", weight: "bold" });
    }
  }
  const domainTop = y;
  y += H + 10;

  // --- exons ----------------------------------------------------------------
  const eventExonIds = new Set(exons.filter(isEventExon).map((e) => e.id));
  trackLabel(fig, x0 - 10, y + (H - 2) / 2, "Coding exons");
  const { bottom } = exonTrack(fig, {
    exons: coding, y, h: H - 2, scale,
    eventExon: (e) => eventExonIds.has(e.id),
  });
  const exonBottom = y + H - 2;
  y = bottom + 6;

  // --- boundary markers, coloured by class ---------------------------------
  const mH = 9;
  trackLabel(fig, x0 - 10, y + mH / 2, "Boundary class");
  for (const b of boundaries) {
    const pos = Number(b.boundary_position_aa ?? b.position);
    if (!Number.isFinite(pos)) continue;
    const bx = scale(pos);
    const cls = boundaryClassOf(b);
    const isSel = selected && boundaryId(b) === boundaryId(selected);
    // A marker spans domain and exon tracks so the relation is readable, while
    // remaining visually subordinate to the domain blocks.
    fig.line(bx, domainTop - 3, bx, exonBottom + 3,
      { stroke: boundaryClassColour(cls), lw: isSel ? P.lw.outline * 1.4 : P.lw.thin,
        opacity: isSel ? 1 : 0.55 });
    fig.line(bx, y, bx, y + mH,
      { stroke: boundaryClassColour(cls), lw: isSel ? P.lw.outline * 1.6 : P.lw.rule });
    if (isSel) {
      fig.circle(bx, y - 3, P.marker * 0.8,
        { fill: boundaryClassColour(cls), stroke: PALETTE.ink, lw: P.lw.thin });
    }
  }
  y += mH + 8;

  // --- validated event ------------------------------------------------------
  // The event's own boundaries are the ones a reader of this figure is most likely
  // to be looking for, so the interval is drawn on the same axis as the classes.
  if (events.length) {
    const evH = 7;
    trackLabel(fig, x0 - 10, y + evH / 2, "Validated event");
    y = eventTrack(fig, { events, y, h: evH, scale, x1 }) + 8;
  }

  // --- candidate context ----------------------------------------------------
  // The same display clusters and lanes as the Domain Architecture track and the
  // Boundary Explorer: one module decides grouping and lane, so a candidate cannot
  // sit in a different lane depending on which figure a reader is looking at.
  if (showCandidates && candidates.length) {
    const layout = candidateDisplayLayout(candidates,
      { selectedId: selectedCandidateId, showAll: showAllCandidates });
    const cH = 7, cGap = 2;
    trackLabel(fig, x0 - 10, y + cH / 2, "Candidates");
    for (const c of layout.visible) {
      const cx0 = scale(c.aa_start);
      const w = Math.max(1, scale(c.aa_end + 1) - cx0);
      const ly = y + layout.laneOf(c) * (cH + cGap);
      fig.rect(cx0, ly, w, cH,
        { fill: PALETTE.candidate, stroke: PALETTE.candidateEdge, lw: P.lw.thin, opacity: 0.7 });
      if (candidateLabelFits(w)) {
        fig.text(cx0 + w / 2, ly + cH - 2, candidateLabel(c),
          { size: "small", anchor: "middle", fill: PALETTE.candidateEdge });
      }
    }
    y += layout.laneCount * (cH + cGap) + 8;
  }

  y = aaAxis(fig, { x0, x1, y, lo: 1, hi: len, scale, label: "Primary protein position (aa)" });

  const present = BOUNDARY_CLASS_ORDER.filter((c) => counts[c]);
  const legend = present.map((c) => [boundaryClassColour(c), boundaryClassLabel(c)]);
  legend.push([PALETTE.exon, "Coding exon"]);
  if (events.length) {
    legend.push([FEATURE_STYLES.validated_event.fill,
      `${events[0].event_label || "Validated event exon"} (established, not a candidate)`]);
  }
  y = fig.legend(x0, y + P.font.legend + 2, legend);
  fig.text(x0, y + P.font.small + 1,
    `Near-domain-edge threshold: ${nearEdgeThreshold} aa`,
    { size: "small", fill: PALETTE.muted });
  y += P.font.small + 4;

  return finalise(fig, y);
}

function selectedAnnotation(b) {
  const pos = b.boundary_position_aa ?? b.position;
  const label = boundaryTransitionLabel(b);
  const dom = b.nearest_domain_label
    || (b.nearest_domain_start != null
      ? `${prettyDomainName(b.nearest_domain_name)} · aa ${b.nearest_domain_start}–${b.nearest_domain_end}`
      : prettyDomainName(b.nearest_domain_name));
  const signed = b.signed_distance_aa ?? b.signed_distance;
  const sign = signed > 0 ? `+${signed}` : String(signed);
  const edge = b.nearest_edge || b.domain_edge_type || "";
  return `Selected: ${label} at aa ${pos} · nearest ${dom} `
    + `· ${edge} edge · signed distance ${sign} aa `
    + `· ${boundaryClassLabel(boundaryClassOf(b))}`;
}

/** "E4 → E5" style transition label, derived from the boundary identifier. */
export function boundaryTransitionLabel(b) {
  if (b.transition_label) return b.transition_label;
  if (b.exon_from_label && b.exon_to_label) return `${b.exon_from_label} → ${b.exon_to_label}`;
  const m = String(b.exon_boundary_id || "").match(/cds(\d+)_end/);
  if (m) {
    const i = Number(m[1]);
    return `E${i} → E${i + 1}`;
  }
  return String(b.exon_boundary_id || "boundary");
}

const boundaryId = (b) => b?.exon_boundary_id || b?.boundary_id
  || `${b?.boundary_position_aa ?? b?.position}`;
const candidateId = (c) => c?.candidate_id || c?.id || c?.rank_label;
const candidateLabel = (c) => c?.rank_label || c?.label || "C1";

// --------------------------------------------------------------------------- //
// Figure 4 — Signed distance to the nearest representative-domain edge
// --------------------------------------------------------------------------- //

/**
 * How far, and in which direction, does each exon boundary sit from the nearest
 * domain edge?
 *
 * This replaces an absolute-distance histogram, which discarded the sign and
 * hid which domain instance was involved.
 */
export function signedDistanceFigureSpec({
  gene, species, proteinId, boundaries = [], nearEdgeThreshold = 5,
  selectedBoundaryId = null, groupByDomain = false, presetName = "double", isoform = "",
}) {
  const P = preset(presetName);
  const rows = [...boundaries]
    .map((b) => ({
      label: boundaryTransitionLabel(b),
      position: Number(b.boundary_position_aa ?? b.position),
      signed: Number(b.signed_distance_aa ?? b.signed_distance ?? 0),
      cls: boundaryClassOf(b),
      edge: b.nearest_edge || b.domain_edge_type || "",
      domain: b.nearest_domain_label || prettyDomainName(b.nearest_domain_name),
      instanceId: b.nearest_domain_instance_id || b.nearest_domain_accession,
      id: boundaryId(b),
    }))
    // Exon order along the protein is the meaningful order here.
    .sort((a, b) => a.position - b.position);

  const gutter = Math.max(58, ...rows.map((r) => textWidth(r.label, P.font.label) + 10));
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 8;
  const rowH = Math.max(11, P.font.label + 5);

  const maxAbs = Math.max(10, ...rows.map((r) => Math.abs(r.signed)));
  const lim = Math.ceil(maxAbs / 10) * 10;
  const scale = (v) => x0 + ((v + lim) / (2 * lim)) * (x1 - x0);

  // Optional grouping inserts a heading row per domain instance.
  const groups = [];
  if (groupByDomain) {
    const seen = new Map();
    for (const r of rows) {
      if (!seen.has(r.instanceId)) seen.set(r.instanceId, { name: r.domain, rows: [] });
      seen.get(r.instanceId).rows.push(r);
    }
    groups.push(...seen.values());
  } else {
    groups.push({ name: null, rows });
  }

  const bodyH = groups.reduce(
    (acc, g) => acc + g.rows.length * rowH + (g.name ? P.font.label + 6 : 0), 0);
  const fig = createFigure({ preset: presetName, height: bodyH + 170 });

  let y = headerBlock(fig, {
    title: gene, species, isoform,
    subtitle: `${proteinId} · signed distance from each internal coding-exon boundary `
      + `to the nearest representative-domain edge`,
    question: "How far, and on which side, does each exon boundary sit relative to "
      + "the nearest domain edge?",
  });
  y += 8;

  const top = y;
  const bottom = top + bodyH;

  // Near-edge band and the domain-edge zero line.
  fig.rect(scale(-nearEdgeThreshold), top - 3,
    scale(nearEdgeThreshold) - scale(-nearEdgeThreshold), bodyH + 6,
    { fill: PALETTE.grid, stroke: "none", opacity: 0.9 });
  fig.text(scale(0), top - 6, `±${nearEdgeThreshold} aa`,
    { size: "small", anchor: "middle", fill: PALETTE.muted });
  fig.line(scale(0), top - 3, scale(0), bottom + 3,
    { stroke: PALETTE.ink, lw: P.lw.rule });

  let ry = top;
  for (const g of groups) {
    if (g.name) {
      fig.text(P.margin.left, ry + P.font.label, g.name,
        { size: "label", weight: "bold", fill: PALETTE.ink });
      ry += P.font.label + 6;
    }
    for (const r of g.rows) {
      const cy = ry + rowH / 2;
      const isSel = selectedBoundaryId && r.id === selectedBoundaryId;
      fig.text(x0 - 6, cy + P.font.label * 0.35, r.label,
        { size: "label", anchor: "end",
          fill: isSel ? PALETTE.ink : PALETTE.muted,
          weight: isSel ? "bold" : "normal" });
      // Stem from the domain edge to the boundary: the lollipop makes both the
      // sign and the magnitude readable at a glance.
      fig.line(scale(0), cy, scale(r.signed), cy,
        { stroke: boundaryClassColour(r.cls), lw: P.lw.rule, opacity: 0.75 });
      const cx = scale(r.signed);
      const rad = P.marker * (isSel ? 1.35 : 1);
      if (r.edge === "start") {
        // Open marker = distance measured to a domain start edge.
        fig.circle(cx, cy, rad, { fill: PALETTE.paper,
          stroke: boundaryClassColour(r.cls), lw: P.lw.outline });
      } else {
        // Filled marker = distance measured to a domain end edge.
        fig.circle(cx, cy, rad, { fill: boundaryClassColour(r.cls),
          stroke: isSel ? PALETTE.ink : "none", lw: P.lw.thin });
      }
      // The value sits on the far side of the marker so it cannot collide with
      // the stem, and flips to the near side when the marker is close enough to
      // an edge that the label would run into the row label or the margin.
      const vLabel = r.signed > 0 ? `+${r.signed}` : String(r.signed);
      const vw = textWidth(vLabel, P.font.small);
      let outward = r.signed >= 0 ? 1 : -1;
      if (outward < 0 && cx - rad - 3 - vw < x0) outward = 1;
      else if (outward > 0 && cx + rad + 3 + vw > x1) outward = -1;
      fig.text(cx + outward * (rad + 3), cy + P.font.small * 0.35, vLabel,
        { size: "small", anchor: outward > 0 ? "start" : "end", fill: PALETTE.muted });
      ry += rowH;
    }
  }

  y = bottom + 8;
  const { major } = axisTicks(-lim, lim, 8);
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < -lim || t > lim) continue;
    const x = scale(t);
    fig.line(x, y, x, y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(x, y + 4 + P.font.tick, t > 0 ? `+${t}` : String(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 4;
  fig.text((x0 + x1) / 2, y, "Signed distance to nearest domain edge (aa) · 0 = domain edge",
    { size: "label", anchor: "middle" });
  y += P.font.label + 6;

  const classes = [...new Set(rows.map((r) => r.cls))];
  y = fig.legend(x0, y, classes.map((c) => [boundaryClassColour(c), boundaryClassLabel(c)]));
  fig.text(x0, y + P.font.small + 1,
    "Open marker: distance to a domain start edge · "
    + "filled marker: distance to a domain end edge · "
    + "negative: upstream of that edge",
    { size: "small", fill: PALETTE.muted });
  y += P.font.small + 4;

  return finalise(fig, y);
}

// --------------------------------------------------------------------------- //
// Figure 5 — Boundary-class summary (compact second panel)
// --------------------------------------------------------------------------- //

/** How are the internal coding-exon boundaries distributed across classes? */
export function boundaryClassSummarySpec({
  gene, species, proteinId, boundaries = [], nearEdgeThreshold = 5, presetName = "compact", isoform = "",
}) {
  const P = preset(presetName);
  const counts = {};
  for (const b of boundaries) {
    const c = boundaryClassOf(b);
    counts[c] = (counts[c] || 0) + 1;
  }
  const rows = BOUNDARY_CLASS_ORDER.filter((c) => counts[c])
    .map((c) => ({ cls: c, n: counts[c] }));
  const total = boundaries.length || 1;

  const gutter = Math.max(...rows.map((r) => textWidth(boundaryClassLabel(r.cls), P.font.label)))
    + 12;
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 26;
  const rowH = P.font.label + 8;
  const fig = createFigure({ preset: presetName, height: rows.length * rowH + 120 });

  let y = headerBlock(fig, {
    title: gene, species, isoform,
    subtitle: `${proteinId} · ${boundaries.length} internal coding-exon boundaries`,
    question: "How are exon boundaries distributed across domain-relation classes?",
  });
  y += 6;

  // Round the count axis up to a whole tick so the bars can be read off it.
  const rawMax = Math.max(...rows.map((r) => r.n));
  const { major, minor } = axisTicks(0, rawMax, 4);
  const axisMax = Math.max(rawMax, major[major.length - 1] || rawMax) || 1;
  const cx = (v) => x0 + ((x1 - x0) * v) / axisMax;
  const top = y;

  // Gridlines sit behind the bars so a count can be estimated without the label.
  for (const t of major) {
    if (t <= 0) continue;
    fig.line(cx(t), top, cx(t), top + rows.length * rowH,
      { stroke: PALETTE.grid, lw: P.lw.thin });
  }

  for (const r of rows) {
    const cy = y + rowH / 2;
    fig.text(x0 - 6, cy + P.font.label * 0.35, boundaryClassLabel(r.cls),
      { size: "label", anchor: "end" });
    const w = cx(r.n) - x0;
    fig.rect(x0, y + 2, Math.max(0.8, w), rowH - 6,
      { fill: boundaryClassColour(r.cls), stroke: "none" });
    fig.text(x0 + w + 4, cy + P.font.label * 0.35,
      `${r.n} (${Math.round((r.n / total) * 100)}%)`, { size: "small", fill: PALETTE.ink });
    y += rowH;
  }

  // Count axis.
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of minor) {
    if (t <= 0 || t > axisMax) continue;
    fig.line(cx(t), y, cx(t), y + 2, { stroke: PALETTE.axis, lw: P.lw.thin, opacity: 0.7 });
  }
  for (const t of major) {
    if (t < 0 || t > axisMax) continue;
    fig.line(cx(t), y, cx(t), y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(cx(t), y + 4 + P.font.tick, String(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick + 3;
  fig.text((x0 + x1) / 2, y + P.font.label, "Number of internal exon boundaries",
    { size: "label", anchor: "middle" });
  y += P.font.label + 8;
  fig.text(P.margin.left, y,
    `Near-domain-edge threshold: ${nearEdgeThreshold} aa · `
    + "classes are mutually exclusive", { size: "small", fill: PALETTE.muted });
  y += P.font.small + 4;

  return finalise(fig, y);
}

// --------------------------------------------------------------------------- //
// Output helpers
// --------------------------------------------------------------------------- //

/**
 * Trim the canvas to what was actually drawn and expose the render backends.
 * Builders reserve a generous height up front because track heights depend on
 * label placement, which is only known once the marks exist.
 */
export function finalise(fig, contentBottom) {
  fig.resize(contentBottom + fig.preset.margin.bottom);
  return fig;
}

/** Tab-separated source table for a figure, so every plot ships its numbers. */
export function tsv(rows, columns) {
  const head = columns.join("\t");
  const body = rows.map((r) => columns.map((c) => {
    const v = r[c];
    return v == null ? "" : String(v).replace(/[\t\n\r]/g, " ");
  }).join("\t"));
  return [head, ...body].join("\n") + "\n";
}
