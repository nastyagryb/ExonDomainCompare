import { useEffect, useMemo, useState } from "react";
import { fileUrl } from "../../api";
import { Badge, Spinner, Empty, AvailabilityState } from "../../ui";
import { useIndex, unavailableState } from "./common";
import { useScientificSelection } from "../../components/ScientificSelectionContext";
// The interactive alignment and the exported alignment figures share one data
// model, one row order, one set of column definitions and one palette, so the
// screen and the figure cannot describe the same alignment differently.
import {
  ALN_COLOURS, aaToColumn, alignmentProfile, identityLabel, identityPct, orderRows,
} from "./alignmentFigure";

const TAXA = ["all", "Primates", "Other mammals", "Birds", "Reptiles", "Amphibians", "Teleost fish"];

/** Legend swatch in a figure colour, outlined so pale fills stay visible. */
const swatch = (fill) => ({ background: fill, border: `1px solid ${ALN_COLOURS.gapEdge}` });

export default function MsaExplorer({ preloaded, species, isoformMode, focusCandidate }) {
  const { data, loading } = useIndex((client) => client.msa(), preloaded);
  const isIsoform = isoformMode || data?.mode === "isoform_alignment"
    || (data?.tabs?.length === 1 && data.tabs[0]?.key === "isoform");
  const [tab, setTab] = useState(isIsoform ? "isoform" : "combined_cassette");
  const [search, setSearch] = useState("");
  const [taxon, setTaxon] = useState("all");
  const [includeReview, setIncludeReview] = useState(true);
  const [highlight, setHighlight] = useState(species || null);
  const [offset, setOffset] = useState(0);
  const WINDOW = 110;
  useEffect(() => {
    // eslint-disable-next-line react-hooks/set-state-in-effect
    if (species) setHighlight(species);
  }, [species]);
  useEffect(() => {
    // eslint-disable-next-line react-hooks/set-state-in-effect
    if (isIsoform) setTab("isoform");
  }, [isIsoform]);

  if (loading) return <Spinner label="Loading alignments…" />;
  if (!data?.available) {
    const why = unavailableState(data, "Alignment");
    return <AvailabilityState why={why} />;
  }

  if (isIsoform && data.alignments?.isoform?.available) {
    return <IsoformMsaView data={data} focusCandidate={focusCandidate} />;
  }

  const tabs = data.tabs || [];
  const isDisc = tab === "discriminating";
  const aln = data.alignments?.[tab];

  const allRows = (aln?.rows || []);
  const rows = allRows.filter((r) => {
    if (taxon !== "all" && r.taxon_group !== taxon) return false;
    if (!includeReview && r.is_review) return false;
    if (search && !`${r.display_species_name} ${r.species}`.toLowerCase().includes(search.toLowerCase())) return false;
    return true;
  });

  const nCols = aln?.n_columns || 0;
  const windowed = tab === "full_length" && nCols > WINDOW;
  const start = windowed ? offset : 0;
  const end = windowed ? Math.min(nCols, offset + WINDOW) : nCols;
  const discCols = tab === "combined_cassette" ? new Set(data.discriminating_columns_combined || []) : new Set();

  return (
    <div className="viewer msa-viewer">
      <div className="seg tabs-seg">
        {tabs.map((t) => (
          <button key={t.key} className={`seg-btn${tab === t.key ? " on" : ""}`}
            onClick={() => { setTab(t.key); setOffset(0); }}
            disabled={t.key !== "discriminating" && !data.alignments?.[t.key]?.available}>{t.label}</button>
        ))}
      </div>

      {isDisc ? (
        <DiscriminatingMsa data={data} />
      ) : !aln?.available ? (
        <Empty title="Alignment not available" />
      ) : (
        <>
          <div className="viewer-controls">
            <input className="search" placeholder="Find species…" value={search} onChange={(e) => setSearch(e.target.value)} />
            <select value={taxon} onChange={(e) => setTaxon(e.target.value)}>
              {TAXA.map((t) => <option key={t} value={t}>{t === "all" ? "All taxa" : t}</option>)}
            </select>
            <select value={highlight || ""} onChange={(e) => setHighlight(e.target.value)}>
              <option value="">Highlight species…</option>
              {allRows.map((r) => <option key={`${r.species}-${r.isoform}`} value={r.species}>{r.display_species_name} {r.isoform}</option>)}
            </select>
            <label className="check inline"><input type="checkbox" checked={includeReview} onChange={(e) => setIncludeReview(e.target.checked)} /><span>Include review</span></label>
            <span className="spacer" />
            <span className="muted small">{rows.length} rows · {nCols} cols</span>
            {aln.file && <a className="btn ghost sm" href={fileUrl(aln.file)}>FASTA</a>}
          </div>

          {windowed && (
            <div className="col-pager">
              <button className="btn ghost sm" disabled={offset === 0} onClick={() => setOffset(Math.max(0, offset - WINDOW))}>← prev</button>
              <span className="muted small">columns {start + 1}–{end} of {nCols}</span>
              <button className="btn ghost sm" disabled={end >= nCols} onClick={() => setOffset(offset + WINDOW)}>next →</button>
            </div>
          )}

          <MsaGrid rows={rows} start={start} end={end} highlight={highlight} discCols={discCols} />
          <MsaConservation conservation={data.conservation} species={highlight} />
          <div className="legend res-legend">
            <span className="legend-item"><span className="msa-swatch hl" />Highlighted</span>
            <span className="legend-item"><span className="msa-swatch human" />Human reference</span>
            <span className="legend-item"><span className="msa-swatch gapsw" />Gap</span>
            {tab === "combined_cassette" && <span className="legend-item"><span className="gold-tick" />Discriminating column</span>}
          </div>
        </>
      )}
    </div>
  );
}

// Generic (gene-agnostic) per-column annotations for any protein alignment.
//   conserved   — all displayed rows share one residue and no gaps
//   gapContaining — at least one gap in the column
//   substitution — >1 distinct non-gap residue
//   variable    — rows differ by residue OR gap state
//   discriminating — rows differ by residue or gap state (groups/models disagree)
// gap boundaries are transitions into/out of a contiguous gap block;
// contiguous variable columns are merged into variable blocks.
function computeColumnStats(rows, nCols) {
  // One shared computation with the exported figures (alignmentProfile), so a
  // column called variable on screen is variable in the figure as well.
  const primary = rows.find((r) => r.is_primary) || rows[0];
  const profile = alignmentProfile(rows, nCols, primary?.seq || "");
  const stats = [];
  for (let c = 0; c < nCols; c++) {
    stats.push({
      conserved: Boolean(profile.conserved[c]),
      gapContaining: Boolean(profile.gapped[c]),
      substitution: Boolean(profile.substitution[c]),
      variable: Boolean(profile.variable[c]),
      // For a within-species isoform alignment these are the same predicate: the
      // models either agree at a column or they do not.
      discriminating: Boolean(profile.variable[c]),
      conservation: profile.conservation[c],
      gapFraction: profile.gapFraction[c],
    });
  }
  // gap boundaries + variable blocks
  const gapStart = new Set();
  const gapEnd = new Set();
  const variableBlocks = [];
  let runStart = -1;
  for (let c = 0; c < nCols; c++) {
    const g = stats[c].gapContaining;
    const prevG = c > 0 && stats[c - 1].gapContaining;
    const nextG = c < nCols - 1 && stats[c + 1].gapContaining;
    if (g && !prevG) gapStart.add(c);
    if (g && !nextG) gapEnd.add(c);
    if (stats[c].variable) { if (runStart < 0) runStart = c; }
    else if (runStart >= 0) { variableBlocks.push([runStart, c - 1]); runStart = -1; }
  }
  if (runStart >= 0) variableBlocks.push([runStart, nCols - 1]);
  const counts = {
    conserved: stats.filter((s) => s.conserved).length,
    variable: stats.filter((s) => s.variable).length,
    substitution: stats.filter((s) => s.substitution).length,
    gap: stats.filter((s) => s.gapContaining).length,
    discriminating: stats.filter((s) => s.discriminating).length,
  };
  return { stats, gapStart, gapEnd, variableBlocks, counts };
}

// Contiguous runs of columns satisfying `pred`, used for block navigation.
function blocksOf(nCols, pred) {
  const out = [];
  let start = -1;
  for (let c = 0; c < nCols; c++) {
    if (pred(c)) { if (start < 0) start = c; }
    else if (start >= 0) { out.push([start, c - 1]); start = -1; }
  }
  if (start >= 0) out.push([start, nCols - 1]);
  return out;
}

// Alignment modes. "Candidate-focused region" is one mode among several and is
// never the default, so the complete alignment stays the primary object of study.
const ALN_MODES = [
  ["full", "Full alignment", "Every column of the alignment."],
  ["diff", "Differences to primary", "Residues differing from the primary protein are marked."],
  ["variable", "Variable regions", "Jumps between columns where the shown models differ."],
  ["conserved", "Conserved regions", "Columns where every shown model shares one residue."],
  ["candidate", "Candidate-focused region", "Restricted to the selected candidate interval."],
];

const ZOOM_STEPS = [60, 120, 200, 320, 500];

// Full-alignment minimap: one thin column per alignment column, coloured by
// annotation, with the current viewport as a draggable window.
function AlignmentMinimap({ nCols, stats, band, offset, windowSize, onJump }) {
  const W = 1000, H = 26;
  const scale = W / Math.max(1, nCols);
  const bar = Math.max(0.6, scale);
  const jump = (e) => {
    const box = e.currentTarget.getBoundingClientRect();
    const frac = (e.clientX - box.left) / box.width;
    onJump(Math.round(frac * nCols - windowSize / 2));
  };
  return (
    <div className="aln-minimap" title="Click to jump to an alignment position">
      <svg viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="none" onClick={jump}
        role="img" aria-label="Alignment overview minimap">
        <rect x="0" y="0" width={W} height={H} fill={ALN_COLOURS.trackBase} />
        {band && (
          // Pale fill with an amber edge: the band must locate the candidate
          // without hiding the variable and gap ticks drawn over it.
          <rect x={band[0] * scale} y="0" width={Math.max(1, (band[1] - band[0] + 1) * scale)}
            height={H} fill={ALN_COLOURS.matchBandCell}
            stroke={ALN_COLOURS.candidateEdge} strokeWidth="0.8" />
        )}
        {stats.map((s, c) => {
          if (!s.variable && !s.gapContaining) return null;
          const fill = s.gapContaining ? ALN_COLOURS.gapDensity : ALN_COLOURS.diff;
          return <rect key={c} x={c * scale} y={s.gapContaining ? 13 : 3} width={bar}
            height={10} fill={fill} />;
        })}
        <rect x={offset * scale} y="0" width={Math.max(2, windowSize * scale)} height={H}
          fill="none" stroke={ALN_COLOURS.reference} strokeWidth="1.4" />
      </svg>
      <div className="aln-minimap-key">
        <span><i className="mm-var" style={{ background: ALN_COLOURS.diff }} />variable column</span>
        <span><i className="mm-gap" style={{ background: ALN_COLOURS.gapDensity }} />gap column</span>
        {band && <span><i className="mm-band" style={{ background: ALN_COLOURS.candidate }} />
          exploratory candidate interval</span>}
        <span className="muted small">viewport {offset + 1}–{Math.min(nCols, offset + windowSize)} of {nCols}</span>
      </div>
    </div>
  );
}

/** Within-species isoform mode — same MsaExplorer shell, isoform_alignment dataset. */
function IsoformMsaView({ data, focusCandidate }) {
  const selection = useScientificSelection();
  const [offset, setOffset] = useState(0);
  // "candidate" is only the initial mode for the dedicated candidate-focused
  // view; the full alignment never opens pre-cropped to the candidate.
  const [mode, setMode] = useState(focusCandidate ? "candidate" : "full");
  const [zoom, setZoom] = useState(1);
  const [showConserved, setShowConserved] = useState(false);
  const [showVariable, setShowVariable] = useState(true);
  const [showGaps, setShowGaps] = useState(true);
  const [showDiscriminating, setShowDiscriminating] = useState(false);
  const [showCandidate, setShowCandidate] = useState(true);
  const [hidden, setHidden] = useState(() => new Set());
  const windowSize = ZOOM_STEPS[zoom];
  const diffMode = mode === "diff";
  const aln = data.alignments?.isoform || {};
  // Primary first, then accession — the row order of the exported figures.
  const allRows = useMemo(() => orderRows(aln.rows || []), [aln.rows]);
  const rows = useMemo(() => allRows.filter((r) => !hidden.has(r.protein_id)), [allRows, hidden]);
  const primaryRow = rows.find((r) => r.is_primary) || rows[0];
  const primaryAligned = primaryRow?.seq || "";
  const nCols = aln.n_columns || primaryAligned.length;

  const colInfo = useMemo(() => computeColumnStats(rows, nCols), [rows, nCols]);
  const identity = useMemo(() => {
    const m = {};
    for (const r of rows) {
      m[r.protein_id] = r.is_primary ? 100 : identityPct(r.seq || "", primaryAligned);
    }
    return m;
  }, [rows, primaryAligned]);
  // Displayed identity keeps one decimal just below 100%, so a model differing by
  // a single residue is not reported as identical — as in the exported figures.
  const identityText = useMemo(() => {
    const m = {};
    for (const r of rows) {
      m[r.protein_id] = r.is_primary ? "100%" : identityLabel(r.seq || "", primaryAligned);
    }
    return m;
  }, [rows, primaryAligned]);

  const selectedCandidate = selection?.selectedCandidate;
  const affected = selectedCandidate
    ? selection.affectedProteinsFor(selectedCandidate.candidate_id) : null;
  const rawBand = useMemo(() => {
    if (!selectedCandidate || !primaryAligned) return null;
    const c0 = aaToColumn(primaryAligned, selectedCandidate.aa_start);
    const c1 = aaToColumn(primaryAligned, selectedCandidate.aa_end);
    if (c0 == null || c1 == null) return null;
    return [Math.min(c0, c1), Math.max(c0, c1)];
  }, [selectedCandidate, primaryAligned]);
  const band = showCandidate ? rawBand : null;

  // Only the dedicated candidate-focused view scrolls to the candidate; the full
  // alignment stays anchored at the N-terminus so nothing is hidden by default.
  useEffect(() => {
    if (focusCandidate && rawBand) {
      // eslint-disable-next-line react-hooks/set-state-in-effect
      setOffset(Math.max(0, rawBand[0] - 5));
    }
  }, [focusCandidate, rawBand?.[0]]); // eslint-disable-line react-hooks/exhaustive-deps

  // Navigation targets derived from the same column statistics as the display.
  const gapBlocks = useMemo(
    () => blocksOf(nCols, (c) => colInfo.stats[c]?.gapContaining),
    [nCols, colInfo],
  );
  const clampOffset = (v) => Math.max(0, Math.min(Math.max(0, nCols - 1), v));
  // A small fixed lead keeps the target block visible near the left edge. Using a
  // fraction of the window instead would silently swallow short jumps, because
  // subtracting a third of a 120-column window from column 30 lands back at 0.
  const LEAD = 8;
  const goTo = (col) => setOffset(clampOffset(col - LEAD));
  const stepBlock = (blocks, dir) => {
    if (!blocks.length) return;
    // Compare against the anchor the viewport is actually pinned to, so every
    // click advances to the next block rather than re-selecting the current one.
    const here = offset + LEAD;
    const next = dir > 0
      ? blocks.find(([s]) => s > here)
      : [...blocks].reverse().find(([s]) => s < here);
    if (next) goTo(next[0]);
    else if (dir > 0) goTo(blocks[0][0]);            // wrap to the first block
    else goTo(blocks[blocks.length - 1][0]);         // wrap to the last block
  };

  const inCandidateMode = mode === "candidate" && rawBand;
  const viewStart = inCandidateMode ? Math.max(0, rawBand[0] - 10) : offset;
  const viewEnd = inCandidateMode
    ? Math.min(nCols, rawBand[1] + 11)
    : Math.min(nCols, offset + windowSize);
  const cols = [];
  for (let c = viewStart; c < viewEnd; c++) cols.push(c);
  const end = viewEnd;
  const cnt = colInfo.counts;

  const colClass = (c) => {
    const s = colInfo.stats[c] || {};
    let cls = "";
    if (showConserved && s.conserved) cls += " col-conserved";
    if (showVariable && s.variable) cls += " col-variable";
    if (showDiscriminating && s.discriminating) cls += " col-disc";
    if (showGaps && (colInfo.gapStart.has(c) || colInfo.gapEnd.has(c))) cls += " col-gapbound";
    // The "variable regions" / "conserved regions" modes de-emphasise everything
    // outside the region of interest instead of hiding columns, so the coordinate
    // system stays continuous.
    if (mode === "variable" && !s.variable) cls += " col-dim";
    if (mode === "conserved" && !s.conserved) cls += " col-dim";
    return cls;
  };

  // Residue colours come from the shared figure palette rather than from the
  // stylesheet, so the screen and the exported figure encode the same meaning:
  // primary dark, match neutral, difference highlighted, gap light, and the
  // candidate interval as a warmer version of the same cell — never a colour per
  // isoform, which would encode nothing the row label does not already say.
  const cellStyle = (ch, { isPri, inBand, differs }) => {
    if (ch === "-") {
      return { background: ALN_COLOURS.gap, color: ALN_COLOURS.gapEdge };
    }
    if (differs) {
      return { background: ALN_COLOURS.diffCell, color: ALN_COLOURS.diffInk };
    }
    const bg = isPri
      ? (inBand ? ALN_COLOURS.primaryBandCell : ALN_COLOURS.primaryCell)
      : (inBand ? ALN_COLOURS.matchBandCell : ALN_COLOURS.matchCell);
    return { background: bg, color: ALN_COLOURS.ink };
  };

  const isCrossSpecies = data.mode === "cross_species_msa";
  return (
    <div className="viewer msa-viewer iso-msa">
      {!isCrossSpecies && (
        <div className="arch-note info">
          <b>Within-species protein isoform alignment — not cross-species conservation.</b>
        </div>
      )}
      {/* Mode selector — the alignment view, not a set of independent toggles. */}
      <div className="seg aln-mode-seg">
        {ALN_MODES.map(([id, label, hint]) => (
          <button key={id} className={`seg-btn${mode === id ? " on" : ""}`} title={hint}
            disabled={id === "candidate" && !rawBand}
            onClick={() => setMode(id)}>{label}</button>
        ))}
      </div>

      <div className="viewer-controls aln-toolbar">
        <Badge cls="accepted" soft>{aln.tool || "MAFFT"}</Badge>
        <span className="muted small">{rows.length} of {allRows.length} proteins · {nCols} columns</span>
        <span className="spacer" />
        {/* Navigation between the scientifically interesting regions. */}
        <span className="aln-navgroup" title="Jump between contiguous variable blocks">
          <button className="btn ghost sm" disabled={inCandidateMode}
            onClick={() => stepBlock(colInfo.variableBlocks, -1)}>◀ var</button>
          <button className="btn ghost sm" disabled={inCandidateMode}
            onClick={() => stepBlock(colInfo.variableBlocks, 1)}>var ▶</button>
        </span>
        <span className="aln-navgroup" title="Jump between contiguous gap blocks">
          <button className="btn ghost sm" disabled={inCandidateMode}
            onClick={() => stepBlock(gapBlocks, -1)}>◀ gap</button>
          <button className="btn ghost sm" disabled={inCandidateMode}
            onClick={() => stepBlock(gapBlocks, 1)}>gap ▶</button>
        </span>
        <span className="aln-navgroup">
          <button className="btn ghost sm" disabled={inCandidateMode}
            onClick={() => setOffset(0)} title="Jump to the N-terminus">N-term</button>
          <button className="btn ghost sm" disabled={inCandidateMode}
            onClick={() => setOffset(Math.max(0, nCols - windowSize))}
            title="Jump to the C-terminus">C-term</button>
          {rawBand && (
            <button className="btn ghost sm" disabled={inCandidateMode}
              onClick={() => goTo(rawBand[0])}
              title="Scroll the full alignment to the selected candidate">
              ⇥ {selectedCandidate?.rank_label || "candidate"}</button>
          )}
        </span>
        <span className="aln-navgroup" title="Horizontal zoom (alignment columns per view)">
          <button className="btn ghost sm" disabled={zoom === 0 || inCandidateMode}
            onClick={() => setZoom((z) => Math.max(0, z - 1))}>＋</button>
          <button className="btn ghost sm" disabled={zoom >= ZOOM_STEPS.length - 1 || inCandidateMode}
            onClick={() => setZoom((z) => Math.min(ZOOM_STEPS.length - 1, z + 1))}>－</button>
          <button className="btn ghost sm"
            onClick={() => { setZoom(1); setOffset(0); setMode("full"); setHidden(new Set()); }}
            title="Reset mode, zoom, position and sequence selection">Reset view</button>
        </span>
        <button className="btn ghost sm" disabled={inCandidateMode || !offset}
          onClick={() => setOffset(Math.max(0, offset - windowSize))}>←</button>
        <span className="muted small">cols {viewStart + 1}–{end}</span>
        <button className="btn ghost sm" disabled={inCandidateMode || end >= nCols}
          onClick={() => setOffset(offset + windowSize)}>→</button>
      </div>

      {/* Column-annotation toggles, separate from the view mode. */}
      <div className="viewer-controls aln-annot-controls">
        <label className="check inline"><input type="checkbox" checked={showVariable}
          onChange={(e) => setShowVariable(e.target.checked)} /><span>Variable columns</span></label>
        <label className="check inline"><input type="checkbox" checked={showDiscriminating}
          onChange={(e) => setShowDiscriminating(e.target.checked)} /><span>Discriminating columns</span></label>
        <label className="check inline"><input type="checkbox" checked={showGaps}
          onChange={(e) => setShowGaps(e.target.checked)} /><span>Gap boundaries</span></label>
        <label className="check inline"><input type="checkbox" checked={showConserved}
          onChange={(e) => setShowConserved(e.target.checked)} /><span>Conserved columns</span></label>
        {rawBand && (
          <label className="check inline"><input type="checkbox" checked={showCandidate}
            onChange={(e) => setShowCandidate(e.target.checked)} />
            <span>Highlight {selectedCandidate?.rank_label || "candidate"}</span></label>
        )}
      </div>

      <div className="msa-annot-summary">
        <span className="mc-chip">Conserved: <b>{cnt.conserved}</b></span>
        <span className="mc-chip">Variable: <b>{cnt.variable}</b></span>
        <span className="mc-chip">Substitutions: <b>{cnt.substitution}</b></span>
        <span className="mc-chip">Gap columns: <b>{cnt.gap}</b></span>
        <span className="mc-chip">Discriminating: <b>{cnt.discriminating}</b></span>
        <span className="mc-chip">Variable blocks: <b>{colInfo.variableBlocks.length}</b></span>
      </div>

      {selectedCandidate && rawBand && (
        <p className="muted sm">{selectedCandidate.rank_label} · primary aa
          {" "}{selectedCandidate.aa_start}–{selectedCandidate.aa_end}
          {" "}= alignment columns {rawBand[0] + 1}–{rawBand[1] + 1}
          {affected?.size ? ` · ${affected.size} affected isoform(s), marked below` : ""}
          {showCandidate ? "." : " (highlight off)."}</p>
      )}

      <div className="msa-iso-filter">
        <span className="muted small">{isCrossSpecies ? "Species:" : "Isoforms:"}</span>
        {allRows.map((r) => (
          <button key={r.protein_id}
            className={`chip sm${hidden.has(r.protein_id) ? "" : " sel"}`}
            onClick={() => setHidden((h) => {
              const n = new Set(h);
              if (n.has(r.protein_id)) n.delete(r.protein_id); else n.add(r.protein_id);
              return n;
            })}
            disabled={r.is_primary}
            title={r.is_primary ? "reference is always shown" : "toggle row"}>
            {isCrossSpecies ? (r.display_species_name || r.species) : r.protein_id}
            {r.is_primary ? " ★" : ""}</button>
        ))}
        <span className="spacer" />
        {/* Selection presets (Part 2B). All isoforms are shown by default. */}
        <button className="btn ghost sm" disabled={!hidden.size}
          onClick={() => setHidden(new Set())}>Select all</button>
        {/* "Differs" means any column differs by residue or by gap state: a model
            that is 100% identical over its aligned residues but carries a
            deletion is still a difference, and must not be hidden here. */}
        <button className="btn ghost sm" title="Keep the primary plus every model that differs from it"
          onClick={() => setHidden(new Set(allRows
            .filter((r) => !r.is_primary && (r.seq || "") === primaryAligned)
            .map((r) => r.protein_id)))}>Primary + differing</button>
        {allRows.some((r) => r.curation_status) && (
          <button className="btn ghost sm" title="Keep only curated protein models"
            onClick={() => setHidden(new Set(allRows
              .filter((r) => !r.is_primary && r.curation_status && r.curation_status !== "curated")
              .map((r) => r.protein_id)))}>Curated only</button>
        )}
        <button className="btn ghost sm" disabled={!hidden.size}
          onClick={() => setHidden(new Set())}>Reset selection</button>
      </div>

      {/* Full-alignment minimap — makes the complete alignment navigable without
          uncontrolled horizontal scrolling (Part 2D). */}
      {!inCandidateMode && nCols > windowSize && (
        <AlignmentMinimap nCols={nCols} stats={colInfo.stats} gapStart={colInfo.gapStart}
          band={band} offset={offset} windowSize={windowSize}
          onJump={(c) => setOffset(clampOffset(c))} />
      )}

      <div className="msa-scroll">
        <table className="msa-grid">
          <tbody>
            <tr className="msa-ruler">
              <td className="msa-name">Column</td>
              {cols.map((c) => <td key={c} className={`msa-cell msa-tick${colClass(c)}`}>
                {(c + 1) % 10 === 0 ? c + 1 : ""}</td>)}
            </tr>
            {rows.map((r) => {
              const seq = r.seq || "";
              const isAff = affected?.has(r.protein_id);
              const isPri = r.is_primary;
              const meta = [
                r.transcript_id,
                isPri ? "primary" : "alternative",
                r.curation_status,
                r.protein_length != null ? `${r.protein_length} aa` : null,
                identity[r.protein_id] != null
                  ? `${identityText[r.protein_id]} id to ${isCrossSpecies ? "reference" : "primary"}` : null,
              ].filter(Boolean);
              return (
                <tr key={r.protein_id} className={`${isPri ? "msa-primary" : ""}${isAff ? " msa-affected" : ""}`}>
                  {/* Sticky row label carrying the full model identity (Part 2B). */}
                  <td className="msa-name aln-name"
                    title={`${isCrossSpecies ? `${r.display_species_name} · ` : ""}${r.protein_id} · ${meta.join(" · ")}`}>
                    <span className="aln-name-id">
                      {isCrossSpecies
                        ? <>{r.display_species_name || r.species} <span className="msa-iso">{r.protein_id}</span></>
                        : <>{r.protein_id}{isPri ? " ★" : ""}</>}
                      {isCrossSpecies && isPri ? " ★" : ""}
                    </span>
                    <span className="aln-name-meta">
                      {r.transcript_id && <span className="msa-iso">{r.transcript_id}</span>}
                      {r.protein_length != null && <span className="msa-iso">{r.protein_length} aa</span>}
                      {r.curation_status && (
                        <span className={`aln-tag${r.curation_status === "curated" ? " curated" : ""}`}>
                          {r.curation_status}</span>)}
                      <span className="aln-tag">{isPri ? "primary" : "alternative"}</span>
                      {identity[r.protein_id] != null && (
                        <b className="aln-ident">{identityText[r.protein_id]}</b>)}
                    </span>
                  </td>
                  {cols.map((c) => {
                    const ch = seq[c] || "-";
                    const inBand = Boolean(band && c >= band[0] && c <= band[1]);
                    const differs = !isPri && ch !== "-" && primaryAligned[c]
                      && primaryAligned[c] !== "-" && ch !== primaryAligned[c];
                    const gap = ch === "-";
                    return (
                      <td key={c}
                        className={`msa-cell${gap ? " gap" : ""}${inBand ? " msa-band" : ""}${diffMode && differs ? " msa-diff" : ""}${colClass(c)}`}
                        style={cellStyle(ch, { isPri, inBand, differs })}
                        title={differs
                          ? `column ${c + 1}: ${primaryAligned[c]} → ${ch} (differs from the primary)`
                          : undefined}>{ch}</td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
      {/* Same colour semantics, same wording as the exported figure legends. A
          thin border keeps the palest swatches visible on white. */}
      <div className="legend res-legend">
        <span className="legend-item">
          <span className="msa-swatch" style={swatch(ALN_COLOURS.primaryCell)} />
          primary protein residue (reference)</span>
        <span className="legend-item">
          <span className="msa-swatch" style={swatch(ALN_COLOURS.matchCell)} />
          residue matching the primary</span>
        <span className="legend-item">
          <span className="msa-swatch" style={swatch(ALN_COLOURS.diffCell)} />
          residue differing from primary</span>
        <span className="legend-item">
          <span className="msa-swatch" style={swatch(ALN_COLOURS.gap)} />
          gap in this isoform (insertion / deletion)</span>
        {band && <span className="legend-item">
          <span className="msa-swatch" style={swatch(ALN_COLOURS.matchBandCell)} />
          exploratory candidate interval</span>}
        <span className="legend-item"><span className="msa-swatch hl" />affected isoform</span>
        {showVariable && <span className="legend-item"><span className="col-swatch variable" />variable column (models disagree)</span>}
        {showDiscriminating && <span className="legend-item"><span className="col-swatch disc" />discriminating column</span>}
        {showGaps && <span className="legend-item"><span className="col-swatch gapbound" />gap boundary</span>}
        {showConserved && <span className="legend-item"><span className="col-swatch conserved" />conserved column</span>}
      </div>
    </div>
  );
}

function MsaGrid({ rows, start, end, highlight, discCols }) {
  const cols = [];
  for (let c = start; c < end; c++) cols.push(c);
  return (
    <div className="msa-scroll">
      <table className="msa-grid">
        <tbody>
          {rows.map((r) => {
            const cls = `${r.is_human ? "msa-human " : ""}${r.species === highlight ? "msa-hl " : ""}${r.is_review ? "msa-review " : ""}`;
            return (
              <tr key={`${r.species}-${r.isoform}`} className={cls}>
                <td className="msa-name" title={`${r.display_species_name} · ${r.isoform}${r.is_review ? " · review" : ""}`}>
                  <span className={`iso-mini-dot iso-${r.isoform.toLowerCase()}`} />{r.display_species_name} <span className="msa-iso">{r.isoform}</span>
                </td>
                {cols.map((c) => {
                  const ch = r.seq[c] || "-";
                  const gap = ch === "-";
                  return <td key={c} className={`msa-cell${gap ? " gap" : ""}${discCols.has(c + 1) ? " disc" : ""}`}>{ch}</td>;
                })}
              </tr>
            );
          })}
        </tbody>
      </table>
      {rows.length === 0 && <p className="muted pad">No rows match.</p>}
    </div>
  );
}

function MsaConservation({ conservation, species }) {
  const rows = (conservation || []).filter((c) => c.species === species && ["cassette", "left_boundary_window", "right_boundary_window"].includes(c.region_type));
  if (!species || !rows.length) return null;
  return (
    <div className="msa-cons">
      <span className="mc-title">Conservation · {species.replaceAll("_", " ")}</span>
      {rows.map((c, i) => (
        <span key={i} className="mc-chip" title={c.conservation_status}>
          {c.region_type.replace("_window", "")}: <b>{c.mean_conservation_score ?? "—"}</b> <em>(gap {c.mean_gap_fraction ?? "—"})</em>
        </span>
      ))}
    </div>
  );
}

function DiscriminatingMsa({ data }) {
  const items = (data.discriminating || []).filter((d) => d.is_discriminating).sort((a, b) => (a.combined_alignment_col || 0) - (b.combined_alignment_col || 0));
  const combined = data.alignments?.combined_cassette;
  if (!items.length) return <Empty title="No discriminating residues" />;
  return (
    <div className="disc-msa">
      <table className="mini-tbl">
        <thead><tr><th>Human pos</th><th>Combined col</th><th>IIIb major</th><th>IIIc major</th><th>Class</th><th>Score</th></tr></thead>
        <tbody>
          {items.map((d, i) => (
            <tr key={i}>
              <td>{d.human_reference_residue_index ?? "—"}</td>
              <td><span className="gold-tick inline" /> {d.combined_alignment_col ?? "—"}</td>
              <td><b className="aa">{d.IIIb_major_aa}</b></td>
              <td><b className="aa">{d.IIIc_major_aa}</b></td>
              <td>{d.position_class}</td>
              <td>{d.discriminating_score ?? "—"}</td>
            </tr>
          ))}
        </tbody>
      </table>
      {combined?.available && (
        <>
          <p className="muted small">Combined cassette alignment · discriminating columns marked in gold:</p>
          <MsaGrid rows={combined.rows} start={0} end={combined.n_columns} highlight="homo_sapiens"
            discCols={new Set(data.discriminating_columns_combined || [])} />
        </>
      )}
    </div>
  );
}
