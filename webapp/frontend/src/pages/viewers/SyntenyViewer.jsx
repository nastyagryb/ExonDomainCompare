import { useEffect, useMemo, useState } from "react";
import { fileUrl } from "../../api";
import { Badge, Spinner, Empty, Drawer, Field, Menu } from "../../ui";
import { useIndex } from "./common";
import SyntenyNeighbourhood from "./SyntenyNeighbourhood";
import { legendEntries, normaliseSyntenyIndex, orthologyStyle } from "./syntenyModel";
import {
  syntenyNeighbourhoodFigureSpec, neighbourConservationMatrixFigureSpec,
  syntenyRowsTsv,
} from "./syntenyFigures";
import {
  downloadFigureSvg, downloadFigurePdf, downloadFigurePng, downloadFigureTsv,
} from "./figureExport";

// Local synteny for every gene and every dataset.
//
// The rows come from the backend contract `shared_synteny_v1` and are drawn by
// the shared renderer, so the generic single-species view, the generic
// comparative view and FGFR2 are the same view with different data. The old
// implementation laid the loci out as a horizontally scrolling row of
// fixed-width buttons, which quietly pushed the outermost locus out of sight
// whenever the container was narrower than the row — the target then looked
// centred while a real neighbour had disappeared.

export default function SyntenyViewer({
  preloaded, species, embedded, generic, anchorSymbol, tsvPath: tsvPathProp,
}) {
  const { data, loading } = useIndex((client) => client.synteny(), preloaded);
  const [sel, setSel] = useState(species || null);
  const [compareHuman, setCompareHuman] = useState(false);
  const [picked, setPicked] = useState(null);
  const gene = picked?.locus || null;

  useEffect(() => {
    // eslint-disable-next-line react-hooks/set-state-in-effect
    if (species) setSel(species);
  }, [species]);

  const rows = useMemo(
    () => normaliseSyntenyIndex(data, { gene: anchorSymbol }), [data, anchorSymbol]);
  const current = useMemo(
    () => rows.find((r) => r.speciesId === (sel || species)) || rows[0],
    [rows, sel, species]);
  const humanRow = useMemo(
    () => rows.find((r) => r.speciesId === "homo_sapiens")
      || (data?.human_reference
        ? normaliseSyntenyIndex({ species: [data.human_reference] })[0] : null),
    [rows, data]);

  // The comparison control belongs to the selected species: a stale "on" state
  // carried over from another species would claim a reference row that this
  // species may not have.
  const canCompare = Boolean(humanRow) && current?.speciesId !== humanRow?.speciesId;
  useEffect(() => {
    // eslint-disable-next-line react-hooks/set-state-in-effect
    setCompareHuman(false);
  }, [current?.speciesId]);

  if (loading) return <Spinner label="Loading synteny…" />;

  if (!rows.length) {
    const status = data?.synteny_status;
    const title = status === "not_applicable"
      ? "Synteny not applicable for this run"
      : status === "not_computed"
        ? "Synteny not computed for this run"
        : "Local gene neighbourhood unavailable";
    return (
      <div className="viewer synteny-viewer">
        <Empty title={title} hint={data?.synteny_reason
          || "No local gene-neighbourhood table was found for this run."} />
        <p className="muted small">
          Local synteny provides supporting evidence for locus identity and genomic
          context, but does not by itself establish isoform or event identity.
        </p>
        {(data?.extraction_warning || data?.synteny_missing_source) && (
          <details className="tech-prov">
            <summary>Technical details</summary>
            {data.extraction_warning && (
              <Field label="Extractor status" wide>{data.extraction_warning}</Field>
            )}
            {data.synteny_missing_source && (
              <Field label="Missing run-local raw annotation" wide>
                <code>{data.synteny_missing_source}</code>
              </Field>
            )}
          </details>
        )}
      </div>
    );
  }

  const targetGene = data?.gene_symbol || data?.target_symbol || anchorSymbol
    || current?.target?.symbol || "target gene";
  const multiSpecies = rows.length > 1;
  const referenceRow = compareHuman && canCompare ? humanRow : null;
  const shown = [referenceRow, current].filter(Boolean);
  const legend = legendEntries(shown);
  const stem = `synteny_${targetGene}_${current?.speciesId || "species"}`.toLowerCase();

  const buildLocal = () => syntenyNeighbourhoodFigureSpec({
    gene: targetGene, rows: [current],
  });
  const buildCompared = () => syntenyNeighbourhoodFigureSpec({
    gene: targetGene, rows: shown,
    subtitle: "Selected species against the human reference on one slot grid",
  });
  const buildComparative = () => syntenyNeighbourhoodFigureSpec({
    gene: targetGene, rows,
    subtitle: "One row per species, target locus centred on a shared slot grid",
  });
  const buildMatrix = () => neighbourConservationMatrixFigureSpec({
    gene: targetGene, rows,
  });

  return (
    <div className="viewer synteny-viewer">
      <div className="viewer-head">
        <div>
          <b>{multiSpecies ? "Comparative local synteny" : "Local genomic neighbourhood"}</b>
          <p className="muted sm">
            {multiSpecies
              ? `One row per species, ${targetGene} centred. Locus context, not a `
                + "whole-genome conservation claim."
              : `Single-species locus context around ${targetGene}; no conservation claim.`}
          </p>
        </div>
      </div>

      <div className="viewer-controls">
        {!embedded && rows.length > 1 && (
          <select value={current?.speciesId || ""} onChange={(e) => setSel(e.target.value)}
            aria-label="Species">
            {rows.map((r) => (
              <option key={r.speciesId} value={r.speciesId}>{r.displayName}</option>
            ))}
          </select>
        )}
        {canCompare && (
          <label className="check inline">
            <input type="checkbox" checked={compareHuman}
              onChange={(e) => setCompareHuman(e.target.checked)} />
            <span>Compare to human</span>
          </label>
        )}
        {current?.statusLabel && !generic && (
          <Badge soft title={current.statusDefinition}>{current.statusLabel}</Badge>
        )}
        <span className="spacer" />
        <Menu label="Export figure">
          <button className="menu-item" onClick={() => downloadFigureSvg(
            referenceRow ? buildCompared() : buildLocal(), stem)}>SVG (vector)</button>
          <button className="menu-item" onClick={() => downloadFigurePdf(
            referenceRow ? buildCompared() : buildLocal(), stem)}>PDF (vector)</button>
          <button className="menu-item" onClick={() => downloadFigurePng(
            referenceRow ? buildCompared() : buildLocal(), stem)}>PNG (300 dpi)</button>
          {multiSpecies && (
            <>
              <button className="menu-item" onClick={() => downloadFigureSvg(
                buildComparative(), `${stem}_all_species`)}>All species · SVG</button>
              <button className="menu-item" onClick={() => downloadFigurePdf(
                buildComparative(), `${stem}_all_species`)}>All species · PDF</button>
              <button className="menu-item" onClick={() => downloadFigurePng(
                buildComparative(), `${stem}_all_species`)}>All species · PNG</button>
              <button className="menu-item" onClick={() => downloadFigureSvg(
                buildMatrix(), `${stem}_conservation_matrix`)}>Conservation matrix · SVG</button>
              <button className="menu-item" onClick={() => downloadFigurePdf(
                buildMatrix(), `${stem}_conservation_matrix`)}>Conservation matrix · PDF</button>
              <button className="menu-item" onClick={() => downloadFigurePng(
                buildMatrix(), `${stem}_conservation_matrix`)}>Conservation matrix · PNG</button>
            </>
          )}
          <button className="menu-item" onClick={() => downloadFigureTsv(
            syntenyRowsTsv(multiSpecies ? rows : shown), `${stem}_loci`)}>
            Displayed loci (TSV)
          </button>
        </Menu>
        {(tsvPathProp || data?.source_tables?.neighbours
          || data?.source_tables?.resolved_5neighbor) && (
          <a className="btn ghost sm" href={fileUrl(tsvPathProp
            || data.source_tables.neighbours
            || data.source_tables.resolved_5neighbor)}>Source TSV</a>
        )}
      </div>

      <div className={`synteny-track${current?.isReview ? " is-review" : ""}`}>
        <SyntenyNeighbourhood row={current} referenceRow={referenceRow}
          onSelect={(locus, speciesId) => setPicked({ locus, speciesId })}
          selectedId={picked ? `${picked.speciesId}:${picked.locus.slot_x}` : null} />
      </div>

      {current?.targetPosition && (
        <p className="muted small">
          {targetGene} · {current.targetPosition}
          {current.targetCoordinateSource === "derived_from_neighbour_offsets"
            && " (span reconstructed from the recorded neighbour offsets)"}
          {current.target?.gene_id ? ` · ${current.target.gene_id}` : ""}
        </p>
      )}

      <div className="legend res-legend">
        {legend.map((e) => (
          <span className="legend-item" key={e.cls} title={e.definition}>
            <span className="gene-dot" style={{ background: e.fill, borderColor: e.stroke }} />
            {e.cls === "target" ? `${targetGene} (target, centred)` : e.label}
          </span>
        ))}
        <span className="legend-item">
          <span className="gene-arrow right legend-arrow" />Transcription direction
        </span>
      </div>

      <p className="muted small">
        Local synteny provides supporting evidence for locus identity and genomic
        context, but does not by itself establish isoform or event identity.
      </p>

      <Drawer open={Boolean(gene)} onClose={() => setPicked(null)}
        title={gene?.symbol || ""}
        subtitle={gene ? `${rows.find((r) => r.speciesId === picked.speciesId)?.displayName
          || picked.speciesId} · ${gene.is_target
            ? "target locus" : `${gene.side} ${gene.rank}`}` : ""}>
        {gene && (
          <>
            <div className="drawer-badges">
              <Badge cls={badgeClass(gene.orthology_class)} soft>
                {gene.orthology_label || orthologyStyle(gene.orthology_class).label}
              </Badge>
              {gene.mapping_confidence && gene.mapping_confidence !== "none" && (
                <Badge soft>confidence: {gene.mapping_confidence}</Badge>
              )}
            </div>
            <p className="muted small">
              {gene.orthology_definition || orthologyStyle(gene.orthology_class).definition}
            </p>
            <Field label="Displayed symbol">{gene.symbol || "—"}</Field>
            <Field label="Raw annotation symbol">{gene.source_symbol || "—"}</Field>
            {gene.gene_id && <Field label="Gene ID"><code>{gene.gene_id}</code></Field>}
            {gene.protein_id && <Field label="Protein ID"><code>{gene.protein_id}</code></Field>}
            <Field label="Transcription direction">{gene.strand || "not annotated"}</Field>
            <Field label="Position">{gene.seqid
              ? `${gene.seqid}:${(gene.genomic_start ?? "?").toLocaleString?.()
                  ?? gene.genomic_start}–${(gene.genomic_end ?? "?").toLocaleString?.()
                  ?? gene.genomic_end}`
              : "—"}</Field>
            <Field label={`Slot vs ${targetGene}`}>
              {`slot ${gene.slot_x}`}
              {gene.distance != null && !gene.is_target
                ? ` · ${Number(gene.distance).toLocaleString()} bp` : ""}
            </Field>
            {gene.percent_identity != null && (
              <Field label="% identity / coverage">
                {gene.percent_identity}% / {gene.coverage ?? "—"}
              </Field>
            )}
          </>
        )}
      </Drawer>
    </div>
  );
}

function badgeClass(cls) {
  if (cls === "target" || cls === "exact" || cls === "curated") return "accepted";
  if (cls === "rbh") return "minor";
  if (cls === "weak" || cls === "ambiguous") return "review";
  return "neutral";
}
