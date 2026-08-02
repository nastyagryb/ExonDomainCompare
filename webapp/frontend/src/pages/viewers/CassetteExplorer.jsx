import { useMemo, useState } from "react";
import { fileUrl } from "../../api";
import { Badge, Drawer, Field, Spinner, Empty, AvailabilityState } from "../../ui";
import { useIndex, unavailableState, RESIDUE_LEGEND, RESIDUE_LABEL } from "./common";

const MODES = [
  ["compare", "Human comparison"],
  ["heatmap", "Species heatmap"],
  ["discriminating", "IIIb/IIIc columns"],
];

export default function CassetteExplorer({ preloaded, species, embedded }) {
  const { data, loading } = useIndex((client) => client.cassette(), preloaded);
  const [panel, setPanel] = useState("IIIb");
  const [mode, setMode] = useState("compare");
  const [onlyDiff, setOnlyDiff] = useState(false);
  const [includeReview, setIncludeReview] = useState(false);
  const [sel, setSel] = useState(species || null);
  const [cell, setCell] = useState(null);
  // A species handed down from the page (deep link, cross-view jump) wins over the
  // local dropdown choice. React's documented way to adjust state on a changed prop
  // is during render, not from an effect.
  const [syncedSpecies, setSyncedSpecies] = useState(species);
  if (species && species !== syncedSpecies) {
    setSyncedSpecies(species);
    setSel(species);
  }

  const speciesList = useMemo(() => (data?.species || []), [data]);
  const current = useMemo(
    () => speciesList.find((s) => s.species === (sel || species)) || speciesList[0],
    [speciesList, sel, species]
  );

  if (loading) return <Spinner label="Loading cassette evidence…" />;
  if (!data?.available) {
    // The reason comes from the index, not from this component: a cassette table that
    // was never written must not be reported as a species without a cassette.
    const why = unavailableState(data, "Cassette evidence");
    return <AvailabilityState why={why} />;
  }
  if (data.evidence_level === "sequence_marker") {
    return <SequenceMarkerCassette data={data} species={species} />;
  }

  const humanRef = data.human_reference?.[panel] || [];
  const discriminating = data.discriminating || {};
  const panelData = current?.panels?.[panel];
  const isControlComparison = data.human_reference_role === "reference_control_only";
  const sourceHref = data.source_tables?.residue_map && !isControlComparison
    ? data.source_tables.residue_map
    : (data.source_tables?.comparison || data.source_tables?.motif_map || data.source_tables?.residue_map);

  return (
    <div className="viewer cassette-viewer">
      {isControlComparison && data.note && <div className="info-banner">{data.note}</div>}
      <div className="viewer-controls">
        <div className="seg">
          {["IIIb", "IIIc"].map((p) => (
            <button key={p} className={`seg-btn iso-tint-${p.toLowerCase()}${panel === p ? " on" : ""}`}
              onClick={() => setPanel(p)}>{p}</button>
          ))}
        </div>
        <div className="seg">
          {MODES.map(([id, l]) => (
            <button key={id} className={`seg-btn${mode === id ? " on" : ""}`} onClick={() => setMode(id)}>{l}</button>
          ))}
        </div>
        {mode !== "discriminating" && !embedded && (
          <select value={current?.species || ""} onChange={(e) => setSel(e.target.value)}>
            {speciesList.map((s) => (
              <option key={s.species} value={s.species}>{s.display_species_name}</option>
            ))}
          </select>
        )}
        {mode === "compare" && (
          <label className="check inline"><input type="checkbox" checked={onlyDiff} onChange={(e) => setOnlyDiff(e.target.checked)} /><span>Only differences</span></label>
        )}
        {mode === "heatmap" && (
          <label className="check inline"><input type="checkbox" checked={includeReview} onChange={(e) => setIncludeReview(e.target.checked)} /><span>Include review</span></label>
        )}
        <span className="spacer" />
        <a className="btn ghost sm" href={fileUrl(sourceHref)}>Source TSV</a>
      </div>

      {mode === "compare" && (
        <CompareView humanRef={humanRef} panelData={panelData} discriminating={discriminating}
          panel={panel} onlyDiff={onlyDiff} onCell={setCell} species={current} />
      )}
      {mode === "heatmap" && (
        <HeatmapView speciesList={speciesList} humanRef={humanRef} panel={panel}
          discriminating={discriminating} includeReview={includeReview} selected={current?.species} onCell={setCell} />
      )}
      {mode === "discriminating" && <DiscriminatingView discriminating={discriminating} />}

      <ResidueLegend />

      <Drawer open={Boolean(cell)} onClose={() => setCell(null)}
        title={cell ? `${cell.species} · ${cell.panel} · position ${cell.i}` : ""}
        subtitle={cell ? RESIDUE_LABEL[cell.cls] : ""}>
        {cell && (() => {
          const disc = discriminating[String(cell.i)];
          return (
            <>
              <div className="drawer-badges">
                <Badge cls={residueBadgeCls(cell.cls)}>{RESIDUE_LABEL[cell.cls]}</Badge>
                {cell.is_discriminating && <Badge cls="info" soft>IIIb/IIIc-discriminating column</Badge>}
              </div>
              {cell.is_discriminating && (
                <p className="drawer-lead">
                  This is an IIIb/IIIc-discriminating column — a position where human IIIb and human IIIc differ.
                  The gold marker flags the column; the cell colour describes the selected species residue relative
                  to the selected human {cell.panel} reference (it does not mean the species residue is wrong).
                </p>
              )}
              <Field label="Species / isoform">{cell.species} · {cell.panel}</Field>
              {cell.is_discriminating && disc && (
                <>
                  <Field label="Human IIIb residue">{disc.IIIb_aa || "—"}</Field>
                  <Field label="Human IIIc residue">{disc.IIIc_aa || "—"}</Field>
                </>
              )}
              <Field label={`Human ${cell.panel} reference residue`}>{cell.h_aa} @ position {cell.i}</Field>
              <Field label="Selected species residue">{cell.sp_aa}</Field>
              <Field label="Agreement with reference">{RESIDUE_LABEL[cell.cls]}{cell.agreement_class ? ` (${cell.agreement_class})` : ""}</Field>
              <Field label="Substitution class">{cell.substitution_class || "—"}</Field>
              <Field label="BLOSUM62">{cell.blosum ?? "—"}</Field>
              <Field label="MSA column">{cell.msa_column ?? "—"}</Field>
              <Field label="IIIb/IIIc-discriminating column">{cell.is_discriminating ? "yes — IIIb and IIIc differ here" : "no"}</Field>
            </>
          );
        })()}
      </Drawer>
    </div>
  );
}

function residueBadgeCls(cls) {
  return cls === "identical" ? "accepted" : cls === "gap" ? "neutral" : cls === "conservative" ? "minor" : "review";
}

// Sequence / marker level cassette view for runs without residue-level human-referenced
// agreement (typical custom runs). Shows cassette coordinates, length, IDs and marker
// status per species/isoform, plus the IIIb/IIIc-discriminating columns from the MSA layer.
function SequenceMarkerCassette({ data, species }) {
  const all = data.sequence_evidence || [];
  const [showAll, setShowAll] = useState(false);
  const forSpecies = species ? all.filter((e) => e.species === species) : all;
  const list = (!showAll && forSpecies.length) ? forSpecies : all;
  const discriminating = data.discriminating || {};
  const hasDisc = Object.keys(discriminating).length > 0;

  return (
    <div className="viewer cassette-viewer">
      <div className="info-banner">
        {data.note || "Cassette evidence available from sequence/MSA marker layer."}
      </div>

      {species && forSpecies.length > 0 && all.length > forSpecies.length && (
        <label className="check inline"><input type="checkbox" checked={showAll} onChange={(e) => setShowAll(e.target.checked)} /><span>Show all species in this run</span></label>
      )}

      <div className="seq-ev-grid">
        {list.map((e) => <SeqEvCard key={`${e.species}-${e.isoform}`} e={e} />)}
        {list.length === 0 && <Empty title="No cassette rows for this species" hint="No IIIb/IIIc cassette evidence was emitted for this species in this run." />}
      </div>

      {hasDisc && (
        <div className="seq-ev-disc">
          <h4 className="seq-ev-head">IIIb / IIIc-discriminating columns
            {data.residue_index_basis === "msa_column" && <span className="muted small"> (indexed by cassette MSA column)</span>}
          </h4>
          <DiscriminatingView discriminating={discriminating} />
        </div>
      )}

      <div className="viewer-controls">
        <span className="spacer" />
        {data.source_tables?.coordinate_map && (
          <a className="btn ghost sm" href={fileUrl(data.source_tables.coordinate_map)}>Coordinate TSV</a>
        )}
        {data.source_tables?.cassette_zoom && (
          <a className="btn ghost sm" href={fileUrl(data.source_tables.cassette_zoom)}>Cassette zoom TSV</a>
        )}
      </div>
    </div>
  );
}

function SeqEvCard({ e }) {
  const coords = e.cassette_available && e.cassette_start_aa != null && e.cassette_end_aa != null
    ? `${e.cassette_start_aa}–${e.cassette_end_aa} aa`
    : null;
  const marker = e.visual_review_flag && e.visual_review_flag !== "false"
    ? "flagged for review"
    : "validated (no review flag)";
  return (
    <div className={`seq-ev-card iso-frame-${(e.isoform || "").toLowerCase()}`}>
      <div className="seq-ev-card-head">
        <b>{e.display_species_name}</b>
        <span className={`iso iso-${(e.isoform || "").toLowerCase()}`}>{e.isoform}</span>
        {e.is_review
          ? <Badge cls="review" soft>review</Badge>
          : <Badge cls={e.claim_class || "neutral"} soft>{e.final_claim_status_after_rescue || "—"}</Badge>}
      </div>
      <Field label="Final isoform label">{e.final_isoform_label || "—"}</Field>
      <Field label="Validated exon type">{e.validated_exon_type || "—"}</Field>
      <Field label="Transcript ID"><code>{e.transcript_id || "—"}</code></Field>
      <Field label="Protein ID"><code>{e.protein_id || "—"}</code></Field>
      <Field label="Protein length">{e.protein_length ? `${e.protein_length} aa` : "—"}</Field>
      <Field label="Cassette coordinates">
        {coords || <span className="muted">Sequence-level evidence available; exact genomic exon coordinates not available.</span>}
      </Field>
      <Field label="Cassette length">{e.cassette_length_aa != null ? `${e.cassette_length_aa} aa` : "—"}
        {e.cassette_msa_start_col != null && e.cassette_msa_end_col != null
          ? ` · MSA cols ${e.cassette_msa_start_col}–${e.cassette_msa_end_col}` : ""}</Field>
      <Field label="Marker / validation">{marker}</Field>
      {e.cassette_exons?.length > 0 && (
        <Field label="Cassette exon(s)">{e.cassette_exons.map((x) => x.label || x.id).filter(Boolean).join(", ") || "—"}</Field>
      )}
    </div>
  );
}

function Res({ pos, panel, species, onCell, dim }) {
  // Cell colour ALWAYS reflects residue agreement (human vs species). The
  // discriminating flag is a column-level annotation shown only as a gold marker.
  const disc = pos.is_discriminating;
  const tip = [
    `${species} · ${panel} · position ${pos.i}`,
    `human ${pos.h_aa} → species ${pos.sp_aa}`,
    RESIDUE_LABEL[pos.cls] + (pos.substitution_class ? ` (${pos.substitution_class})` : ""),
    `IIIb/IIIc-discriminating column: ${disc ? "yes — IIIb and IIIc differ here" : "no"}`,
    pos.msa_column ? `MSA column ${pos.msa_column}` : null,
  ].filter(Boolean).join("\n");
  return (
    <button className={`res res-${pos.cls}${disc ? " disc" : ""}${dim ? " dim" : ""}`} title={tip}
      onClick={() => onCell({ ...pos, panel, species })}>
      {pos.sp_aa || "-"}
    </button>
  );
}

function CompareView({ humanRef, panelData, discriminating, panel, onlyDiff, onCell, species }) {
  if (!panelData?.available) {
    return <Empty title={`No ${panel} cassette for ${species?.display_species_name || "this species"}`} hint="No residue row was emitted for this isoform." />;
  }
  const byIndex = new Map(panelData.positions.map((p) => [p.i, p]));
  let cols = humanRef.map((h) => ({ h, sp: byIndex.get(h.i) }));
  if (onlyDiff) cols = cols.filter((c) => c.sp && c.sp.cls !== "identical");
  return (
    <div className="cassette-compare">
      <div className="cc-head">
        <span className="cc-meta">{panelData.n_identical} identical · {panelData.n_diff} substitutions · {panelData.n_gap} gaps</span>
        {panelData.is_review && <Badge cls="review" soft>review row</Badge>}
        {!panelData.is_review && <Badge cls={panelData.claim_class} soft>{panelData.final_claim_status_after_rescue}</Badge>}
      </div>
      <div className="res-rows">
        <div className="res-row">
          <span className="res-row-label">Human {panel}</span>
          <div className="res-track">
            {cols.map(({ h }, idx) => (
              <div key={idx} className={`res ref${discriminating[String(h.i)] ? " disc" : ""}`} title={`human ${panel} ${h.aa} @ ${h.i}${h.property ? "\n" + h.property : ""}`}>{h.aa}</div>
            ))}
          </div>
        </div>
        <div className="res-row">
          <span className="res-row-label">{species?.display_species_name}</span>
          <div className="res-track">
            {cols.map(({ h, sp }, idx) => sp
              ? <Res key={idx} pos={sp} panel={panel} species={species.display_species_name} onCell={onCell} />
              : <div key={idx} className="res res-gap" title={`pos ${h.i}: not available`}>·</div>)}
          </div>
        </div>
        <div className="res-row ticks">
          <span className="res-row-label" />
          <div className="res-track">
            {cols.map(({ h }, idx) => {
              const d = discriminating[String(h.i)];
              return (
                <div key={idx} className="res-tick">
                  {d ? <span className="gold-tick"
                    title={`IIIb/IIIc-discriminating column ${h.i}\nhuman IIIb ${d.IIIb_aa || "?"} vs human IIIc ${d.IIIc_aa || "?"}`} /> : null}
                </div>
              );
            })}
          </div>
        </div>
      </div>
    </div>
  );
}

function HeatmapView({ speciesList, humanRef, panel, discriminating, includeReview, selected, onCell }) {
  const rows = speciesList
    .map((s) => ({ s, p: s.panels?.[panel] }))
    .filter((r) => r.p?.available && (includeReview || !r.p.is_review));
  const indices = humanRef.map((h) => h.i);
  return (
    <div className="cassette-heatmap">
      <div className="chm-scroll">
        <table className="chm">
          <thead>
            <tr>
              <th className="chm-sticky">Species</th>
              {indices.map((i, idx) => (
                <th key={idx} className={discriminating[String(i)] ? "disc-col" : ""} title={discriminating[String(i)] ? `IIIb/IIIc-discriminating column ${i}` : `position ${i}`}>{i}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            <tr className="chm-ref">
              <td className="chm-sticky">Human {panel}</td>
              {humanRef.map((h, idx) => <td key={idx} className="res-mini ref">{h.aa}</td>)}
            </tr>
            {rows.map(({ s, p }) => {
              const byIndex = new Map(p.positions.map((x) => [x.i, x]));
              return (
                <tr key={s.species} className={s.species === selected ? "chm-sel" : ""}>
                  <td className="chm-sticky">{s.display_species_name}{p.is_review && <span className="rev-dot" title="review" />}</td>
                  {indices.map((i, idx) => {
                    const pos = byIndex.get(i);
                    if (!pos) return <td key={idx} className="res-mini res-gap" title={`pos ${i}: n/a`} />;
                    return (
                      <td key={idx} className={`res-mini res-${pos.cls}${pos.is_discriminating ? " disc" : ""}`}
                        title={`${s.display_species_name} · ${panel} · position ${i}\nhuman ${pos.h_aa} → species ${pos.sp_aa}\n${RESIDUE_LABEL[pos.cls]}\nIIIb/IIIc-discriminating column: ${pos.is_discriminating ? "yes — IIIb and IIIc differ here" : "no"}`}
                        onClick={() => onCell({ ...pos, panel, species: s.display_species_name })}>{pos.sp_aa || "-"}</td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
      {rows.length === 0 && <p className="muted pad">No species match (toggle “include review”).</p>}
    </div>
  );
}

function DiscriminatingView({ discriminating }) {
  const items = Object.values(discriminating).sort((a, b) => a.i - b.i);
  if (!items.length) return <Empty title="No IIIb/IIIc-discriminating columns" />;
  return (
    <div className="disc-table">
      <p className="cassette-legend-note">
        Columns where human IIIb and human IIIc differ. These are the positions that distinguish the two isoforms;
        species residues are compared to the selected human reference elsewhere in this explorer.
      </p>
      <table className="mini-tbl">
        <thead><tr><th>Pos</th><th>MSA col</th><th>IIIb</th><th>IIIc</th><th>IIIb property</th><th>IIIc property</th><th>Substitution</th><th>Score</th></tr></thead>
        <tbody>
          {items.map((d) => (
            <tr key={d.i}>
              <td><span className="gold-tick inline" /> {d.i}</td>
              <td>{d.msa_column ?? "—"}</td>
              <td><b className="aa">{d.IIIb_aa}</b></td>
              <td><b className="aa">{d.IIIc_aa}</b></td>
              <td>{d.IIIb_property || "—"}</td>
              <td>{d.IIIc_property || "—"}</td>
              <td>{d.substitution_class || "—"}</td>
              <td>{d.discriminating_score ?? "—"}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function ResidueLegend() {
  return (
    <div className="cassette-legend-wrap">
      <div className="legend res-legend">
        {RESIDUE_LEGEND.map(([c, l]) => <span key={c} className="legend-item"><span className={`res-swatch res-${c}`} />{l}</span>)}
        <span className="legend-item"><span className="gold-tick" />IIIb/IIIc-discriminating column</span>
      </div>
      <p className="cassette-legend-note">
        Gold markers indicate columns that distinguish human IIIb from human IIIc. The cell colour still
        describes the selected species residue relative to the selected human isoform reference — a gold marker
        does not mean the species residue is wrong.
      </p>
    </div>
  );
}
