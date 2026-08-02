// The shared local-synteny renderer.
//
// One component draws the neighbourhood for every gene and every dataset: the
// generic single-species view, the generic comparative view, FGFR2 and every
// future run. It replaces the flex-row of fixed-width buttons that used to
// scroll horizontally, which is what silently hid the outermost locus whenever
// the container was narrower than the row.
//
// Geometry rules that make the view trustworthy:
//   * Loci sit on a slot grid whose half-width is the widest displayed side, so
//     the target is at the exact centre and a human reference row aligns
//     slot-for-slot with the species row.
//   * The drawing is one `viewBox`-scaled SVG. It shrinks with the container
//     instead of overflowing it, so no locus can be clipped at either edge.
//   * An empty slot is drawn as empty. Nothing is invented to fill the grid.

import { useMemo } from "react";
import { PALETTE } from "./figureSpec.js";
import { orthologyStyle, slotGrid } from "./syntenyModel.js";

const VIEW_W = 1000;
const ROW_H = 74;
const HEAD_H = 20;
const GENE_H = 22;
const LABEL_GAP = 13;

/** Shorten a label to the available width, keeping the full text in the tooltip. */
function fit(text, width, size) {
  const per = size * 0.55;
  const max = Math.max(2, Math.floor(width / per));
  const s = String(text || "");
  return s.length <= max ? s : `${s.slice(0, Math.max(1, max - 1))}…`;
}

function GeneArrow({ x, y, w, h, locus, selected, onSelect, speciesId }) {
  const style = orthologyStyle(locus.orthology_class);
  const dir = locus.strand === "-" ? -1 : 1;
  const tip = Math.min(9, w * 0.28);
  // A locus with no annotated strand is drawn as a plain box rather than being
  // given a direction it does not have.
  const points = locus.strand
    ? (dir > 0
      ? `${x},${y} ${x + w - tip},${y} ${x + w},${y + h / 2} ${x + w - tip},${y + h} ${x},${y + h}`
      : `${x + w},${y} ${x + tip},${y} ${x},${y + h / 2} ${x + tip},${y + h} ${x + w},${y + h}`)
    : `${x},${y} ${x + w},${y} ${x + w},${y + h} ${x},${y + h}`;
  const parts = [
    locus.is_target ? `${locus.symbol} — target gene` : locus.symbol,
    locus.orthology_definition,
    locus.strand ? `Transcription direction ${locus.strand}` : "Strand not annotated",
    locus.is_target ? "Central slot" : `${locus.side} ${locus.rank}`,
    locus.distance != null && !locus.is_target
      ? `${Number(locus.distance).toLocaleString()} bp from the target` : "",
    locus.source_symbol && locus.source_symbol !== locus.symbol
      ? `Raw annotation symbol: ${locus.source_symbol}` : "",
  ].filter(Boolean);
  return (
    <g className={`syn-locus${selected ? " is-selected" : ""}`}
      role="button" tabIndex={0}
      onClick={() => onSelect?.(locus, speciesId)}
      onKeyDown={(e) => {
        if (e.key === "Enter" || e.key === " ") onSelect?.(locus, speciesId);
      }}>
      <polygon points={points} fill={style.fill} stroke={style.stroke}
        strokeWidth={locus.is_target ? 1.6 : 1}
        strokeDasharray={style.dashed ? "3 2" : undefined} />
      {selected && (
        <polygon points={points} fill="none" stroke={PALETTE.ink} strokeWidth={2.2}
          strokeOpacity={0.75} />
      )}
      <title>{parts.join("\n")}</title>
    </g>
  );
}

function SpeciesRow({ row, grid, y, colW, onSelect, selectedId, subtitle }) {
  const labelSize = 11;
  const geneY = y + HEAD_H;
  // Long symbols in narrow slots alternate between two label rows so adjacent
  // labels cannot run into each other.
  const needsStagger = row.loci.some(
    (n) => String(n.symbol).length * labelSize * 0.55 > colW - 4);
  return (
    <g>
      <text x={0} y={y + 11} fontSize={12} fontWeight={600} fill={PALETTE.ink}>
        {row.displayName}
      </text>
      {subtitle && (
        <text x={VIEW_W} y={y + 11} fontSize={10.5} textAnchor="end" fill={PALETTE.muted}>
          {subtitle}
        </text>
      )}
      <line x1={0} y1={geneY + GENE_H / 2} x2={VIEW_W} y2={geneY + GENE_H / 2}
        stroke={PALETTE.grid} strokeWidth={1.4} />
      {row.loci.map((locus) => {
        const col = grid.columnOf(locus);
        const x = col * colW + 2;
        const w = colW - 4;
        const id = `${row.speciesId}:${locus.slot_x}`;
        const stagger = needsStagger && col % 2 === 1 ? LABEL_GAP : 0;
        const style = orthologyStyle(locus.orthology_class);
        return (
          <g key={id}>
            <GeneArrow x={x} y={geneY} w={w} h={GENE_H} locus={locus}
              speciesId={row.speciesId}
              selected={selectedId === id} onSelect={onSelect} />
            {stagger > 0 && (
              <line x1={x + w / 2} y1={geneY + GENE_H + 1} x2={x + w / 2}
                y2={geneY + GENE_H + stagger - 1}
                stroke={PALETTE.axis} strokeWidth={0.7} />
            )}
            <text x={x + w / 2} y={geneY + GENE_H + LABEL_GAP + stagger}
              fontSize={labelSize} textAnchor="middle"
              fontWeight={locus.is_target ? 700 : 500}
              fontStyle={style.italic ? "italic" : "normal"}
              fill={locus.is_target ? PALETTE.ink : style.text}>
              {fit(locus.symbol, w + 6, labelSize)}
              <title>{locus.symbol}</title>
            </text>
          </g>
        );
      })}
    </g>
  );
}

/**
 * Target-centred local neighbourhood for one species, optionally above a
 * reference row drawn on the same slot grid.
 */
export default function SyntenyNeighbourhood({
  row, referenceRow, onSelect, selectedId, showCounts = true,
}) {
  const rows = useMemo(
    () => [referenceRow, row].filter(Boolean), [row, referenceRow]);
  const grid = useMemo(() => slotGrid(rows), [rows]);
  if (!row) return null;

  const colW = VIEW_W / grid.columns;
  const height = rows.length * ROW_H + 26;
  const targetX = grid.targetColumn * colW;

  return (
    <div className="synteny-neighbourhood">
      <svg viewBox={`0 0 ${VIEW_W} ${height}`} width="100%" height={height}
        preserveAspectRatio="xMidYMid meet" role="img"
        aria-label={`Local gene neighbourhood around ${row.target.symbol}`}>
        {/* The central target column, so the eye finds the anchor immediately. */}
        <rect x={targetX} y={0} width={colW} height={rows.length * ROW_H}
          fill={PALETTE.domain} fillOpacity={0.06} />
        {rows.map((r, i) => (
          <SpeciesRow key={r.speciesId || i} row={r} grid={grid} colW={colW}
            y={i * ROW_H} onSelect={onSelect} selectedId={selectedId}
            subtitle={r.isHumanReference || (referenceRow && r === referenceRow)
              ? "reference" : ""} />
        ))}
        <text x={targetX + colW / 2} y={rows.length * ROW_H + 14} fontSize={10}
          textAnchor="middle" fill={PALETTE.muted}>target locus</text>
        <text x={0} y={rows.length * ROW_H + 14} fontSize={10} fill={PALETTE.muted}>
          upstream
        </text>
        <text x={VIEW_W} y={rows.length * ROW_H + 14} fontSize={10} textAnchor="end"
          fill={PALETTE.muted}>downstream</text>
      </svg>
      {showCounts && (
        <p className="muted sm syn-counts">
          {row.countsLabel}
          {row.omissionReason ? ` — ${row.omissionReason}` : ""}
        </p>
      )}
    </div>
  );
}
