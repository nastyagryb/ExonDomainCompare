import { useMemo, useState } from "react";
import { Badge, Menu } from "../../ui";
import { useScientificSelection } from "../../components/ScientificSelectionContext";
import { CANON_CLASS_COLOR, CANON_CLASS_LABEL, canonClass } from "./boundaryClasses";
import ComparativeBoundaryMatrix from "./ComparativeBoundaryMatrix";
import ComparativePairedPlot from "./ComparativePairedPlot";
import {
  comparativeMatrixFigureSpec, pairedSignedDistanceFigureSpec,
  consistencySummaryFigureSpec, comparativeArchitectureFigureSpec,
  comparativeLongTsv, comparativeMatrixTsv, comparableMappingTsv,
  isSupported, speciesTag,
} from "./comparativeFigures";
import {
  EMPTY_FILTERS, activeFilterCount, filterComparativeDataset,
} from "./comparativeFilters";
// Structural colours are read from the shared figure palette so the on-screen
// architecture panel and its exported counterpart cannot diverge.
import { PALETTE } from "./figureSpec";
import { speciesCompare } from "./speciesOrder";
import {
  downloadFigureSvg, downloadFigurePdf, downloadFigurePng, downloadFigureTsv,
} from "./figureExport";

// Comparative multi-species Exon–Domain-Boundary Explorer.
//
// Everything shown here comes from the canonical comparative index published by
// src/exondomaincompare/shared_gene_analysis/boundary_dashboard.py. The browser does not decide which
// boundaries are comparable, how they were matched, or how confident the mapping is; a
// second implementation of that logic would be free to disagree with the backend and
// with the exported tables, and a reader could not tell which answer the figures rest
// on. The frontend's job is filtering, layout and export.
//
// With no real comparable evidence the page shows an honest species inventory and says
// why — it never fabricates placeholder cells.

const MODE_LABEL = {
  signed: "Signed distance",
  absolute: "Absolute distance",
  class: "Class only",
};

const CASE_SEVERITY_CLS = { caution: "warn", review: "neutral" };

/** Multi-select chip row; empty selection means "all", which is the honest default. */
function ChipFilter({ label, options, selected, onChange, format }) {
  if (!options.length) return null;
  const toggle = (v) => onChange(selected.includes(v)
    ? selected.filter((x) => x !== v)
    : [...selected, v]);
  return (
    <div className="cbe-filter-group">
      <span className="cbe-filter-label">{label}</span>
      <div className="cbe-chips">
        {options.map((v) => (
          <button key={v} type="button"
            className={`cbe-chip${selected.includes(v) ? " on" : ""}`}
            onClick={() => toggle(v)}>
            {format ? format(v) : v}
          </button>
        ))}
      </div>
    </div>
  );
}

function SpeciesInventory({ rows }) {
  return (
    <div className="card">
      <div className="card-head">
        <h3>Species in this comparative run</h3>
        <Badge cls="neutral" soft>{rows.length} species</Badge>
      </div>
      <div className="table-scroll"><table className="mini-tbl">
        <thead><tr><th>Species</th><th>Taxonomic group</th><th>Primary protein</th>
          <th>Transcript</th><th>Length</th><th>Analysis status</th></tr></thead>
        <tbody>{rows.map((r) => (
          <tr key={r.species_id}>
            <td><i>{r.scientific_name || r.species_id}</i></td>
            <td>{r.taxonomic_group || "—"}</td>
            <td><code>{r.primary_protein || "—"}</code></td>
            <td><code>{r.transcript || "—"}</code></td>
            <td>{r.protein_length ?? "—"} aa</td>
            <td><Badge cls={r.analysis_status === "available" ? "accepted" : "neutral"} soft>
              {r.analysis_status || "—"}</Badge></td>
          </tr>
        ))}</tbody>
      </table></div>
    </div>
  );
}

/** One real row per species for the selected comparable-boundary group. */
function SelectedGroupDetail({ group, stat, openGene, setPage, onOpenSpecies }) {
  if (!group) {
    return (
      <div className="arch-note info">
        Select a matrix cell, a plot row or an inspection case to see the per-species
        detail for one comparable boundary.
      </div>
    );
  }
  const obs = group.per_species_native_positions || [];
  const tentative = !isSupported(group.mapping_status);
  return (
    <>
      <div className="cbe-detail-head">
        <b>{group.comparable_boundary_group_id}</b>
        <Badge cls={tentative ? "warn" : "accepted"} soft>
          {tentative ? "tentative mapping" : "supported mapping"}
        </Badge>
        <span className="muted sm">
          matched by {group.mapping_method === "msa_aligned_position"
            ? "MSA-aligned protein position" : group.mapping_method}
          {group.msa_column != null ? ` · alignment column ${group.msa_column}` : ""}
          {" · mapped in "}{Math.round((group.confidence ?? 0) * 100)}% of analysed species
        </span>
      </div>
      {(openGene || setPage) && obs.length > 0 && (
        <div className="cbe-linked-nav">
          {obs.map((o) => (
            <span key={o.species_id} className="cbe-linked-nav-sp">
              <button type="button" className="btn ghost sm"
                onClick={() => {
                  onOpenSpecies?.(o.species_id, "architecture");
                  openGene?.({ species: o.species_id, tab: "architecture",
                    boundaryId: o.boundary_id, domainId: o.nearest_domain_id });
                }}>
                Open {o.scientific_name || o.species_id} Domain Architecture
              </button>
              <button type="button" className="btn ghost sm"
                onClick={() => {
                  onOpenSpecies?.(o.species_id, "boundary");
                  openGene?.({ species: o.species_id, tab: "boundary",
                    boundaryId: o.boundary_id });
                }}>
                Open {o.scientific_name || o.species_id} Boundary
              </button>
            </span>
          ))}
          <button type="button" className="btn ghost sm"
            onClick={() => {
              const o = obs[0];
              onOpenSpecies?.(o.species_id, "exon");
              openGene?.({ species: o.species_id, tab: "exon",
                boundaryId: o.boundary_id });
            }}>
            Open selected Boundary in Exon Map
          </button>
        </div>
      )}
      {tentative && (
        <div className="arch-note warn">
          The observations below were matched across a small alignment-column offset.
          They describe closely positioned boundaries; this does not establish that they
          are the same junction, and they are not presented as equivalent.
        </div>
      )}
      <div className="table-scroll"><table className="mini-tbl">
        <thead><tr>
          <th>Species</th><th>Primary protein</th><th>Exon transition</th>
          <th>Native pos.</th><th>Column</th><th>Nearest domain instance</th>
          <th>Domain coords</th><th>Edge</th><th>Signed dist.</th>
          <th>Boundary class</th><th>Mapping method</th><th>Mapping confidence</th>
        </tr></thead>
        <tbody>{obs.map((o) => (
          <tr key={o.species_id}>
            <td><i>{o.scientific_name || o.species_id}</i></td>
            <td><code>{o.protein_id}</code></td>
            <td>{o.exon_transition || "—"}</td>
            <td>{o.native_position ?? "—"}</td>
            <td>{o.msa_column ?? "—"}</td>
            <td>{o.nearest_domain_label || <span className="muted">not annotated</span>}</td>
            <td>{o.nearest_domain_start != null
              ? `${o.nearest_domain_start}–${o.nearest_domain_end}` : "—"}</td>
            <td>{o.nearest_edge || "—"}</td>
            <td>{o.signed_distance == null ? "—"
              : `${o.signed_distance > 0 ? "+" : ""}${o.signed_distance} aa`}</td>
            <td>{CANON_CLASS_LABEL[canonClass(o.boundary_class)]}</td>
            <td className="sm">{o.mapping_method}</td>
            <td className="sm">{isSupported(o.mapping_status) ? "supported" : "tentative"}
              <br /><span className="muted">
                {Math.round((o.mapping_confidence ?? 0) * 100)}% species coverage</span></td>
          </tr>
        ))}</tbody>
      </table></div>
      {stat && (
        <div className="muted sm cbe-detail-stat">
          {(stat.raw_signed_distances || []).length <= 2
            ? `Cross-species difference: ${stat.cross_species_difference ?? "—"} aa `
              + "(difference between the two raw observations)"
            : `Observed range: ${(stat.distance_range || []).join(" to ")} aa · median `
              + `${stat.median_signed_distance ?? "—"} aa`}
          {stat.classes_differ && " · the species assign different boundary classes"}
          {!stat.domain_annotation_available_in_all
            && " · at least one species has no representative domain annotated nearby"}
        </div>
      )}
    </>
  );
}

/** Compact comparative local architecture, drawn on a boundary-anchored axis. */
function ComparativeArchitecture({ group, models }) {
  if (!group) return null;
  const obs = (group.per_species_native_positions || []).slice()
    .sort((a, b) => speciesCompare(a.species_id, b.species_id));
  const modelBySpecies = new Map((models || []).map((m) => [m.species_id, m]));
  const columns = [...new Set(obs.map((o) => o.msa_column).filter((c) => c != null))];
  const sharedColumn = columns.length === 1 ? columns[0] : null;
  const HALF = 70;
  const W = 780;
  const LEFT = 96;
  const innerW = W - LEFT - 20;
  const x = (delta) => LEFT + ((delta + HALF) / (2 * HALF)) * innerW;
  const LANE = 42;

  return (
    <div className="cbe-arch">
      <div className="muted sm">
        {sharedColumn != null && isSupported(group.mapping_status)
          ? `Tracks are anchored on the comparable boundary, which both species map to `
            + `alignment column ${sharedColumn}.`
          : "Tracks are anchored on each species' own boundary. The observations do not "
            + "share a single alignment column, so no common coordinate axis is implied."}
      </div>
      <svg viewBox={`0 0 ${W} ${obs.length * LANE + 54}`} className="cbe-arch-svg"
        role="img" aria-label="Comparative local domain architecture around the selected boundary">
        {obs.map((o, i) => {
          const ly = 8 + i * LANE;
          const pos = Number(o.native_position);
          const model = modelBySpecies.get(o.species_id);
          const colour = CANON_CLASS_COLOR[canonClass(o.boundary_class)];
          return (
            <g key={o.species_id}>
              <text x={LEFT - 8} y={ly + 13} textAnchor="end" fontSize="10"
                fontStyle="italic" fill={PALETTE.ink}>
                {speciesTag(o.species_id, o.scientific_name)}
              </text>
              <rect x={LEFT} y={ly + 8} width={innerW} height={2.4} fill={PALETTE.grid} />
              {(model?.representative_domains || []).map((d) => {
                const ds = Number(d.start) - pos;
                const de = Number(d.end) - pos;
                if (de < -HALF || ds > HALF) return null;
                const cs = Math.max(ds, -HALF);
                const ce = Math.min(de, HALF);
                const w = x(ce) - x(cs);
                return (
                  <g key={d.id}>
                    <title>{d.full_label || d.label}</title>
                    <rect x={x(cs)} y={ly} width={w} height={18} rx="2"
                      fill={PALETTE.domain} stroke={PALETTE.ink} strokeWidth="0.5" opacity="0.85" />
                    {w > 54 && (
                      <text x={x(cs) + w / 2} y={ly + 12.5} textAnchor="middle"
                        fontSize="8.5" fill="#fff">{d.short_label || d.label}</text>
                    )}
                  </g>
                );
              })}
              {o.nearest_edge_position != null
                && Math.abs(Number(o.nearest_edge_position) - pos) <= HALF && (
                <line x1={x(Number(o.nearest_edge_position) - pos)} y1={ly - 2}
                  x2={x(Number(o.nearest_edge_position) - pos)} y2={ly + 20}
                  stroke={PALETTE.ink} strokeWidth="0.9" strokeDasharray="2,2" />
              )}
              <line x1={x(0)} y1={ly - 4} x2={x(0)} y2={ly + 22}
                stroke={colour} strokeWidth="1.8" />
              <text x={LEFT} y={ly + 32} fontSize="8.5" fill={PALETTE.muted}>
                {o.exon_transition} · native aa {o.native_position}
                {o.msa_column != null ? ` · column ${o.msa_column}` : ""} ·
                {" "}{o.signed_distance > 0 ? "+" : ""}{o.signed_distance} aa to the
                {" "}{o.nearest_edge} edge of {o.nearest_domain_label || "no annotated domain"}
              </text>
            </g>
          );
        })}
        <line x1={LEFT} y1={obs.length * LANE + 14} x2={LEFT + innerW}
          y2={obs.length * LANE + 14} stroke={PALETTE.axis} strokeWidth="1" />
        {[-HALF, -35, 0, 35, HALF].map((t) => (
          <g key={t}>
            <line x1={x(t)} y1={obs.length * LANE + 14} x2={x(t)}
              y2={obs.length * LANE + 18} stroke={PALETTE.axis} strokeWidth="1" />
            <text x={x(t)} y={obs.length * LANE + 28} textAnchor="middle" fontSize="8.5"
              fill={PALETTE.muted}>{t > 0 ? `+${t}` : t}</text>
          </g>
        ))}
        <text x={LEFT + innerW / 2} y={obs.length * LANE + 44} textAnchor="middle"
          fontSize="9" fill={PALETTE.ink}>
          Amino-acid offset from the comparable boundary · 0 = boundary in each species
        </text>
      </svg>
    </div>
  );
}

export default function ComparativeBoundaryDashboard({
  multi, gene, threshold = 5, models = [], setPage, openGene,
}) {
  const [filters, setFilters] = useState(EMPTY_FILTERS);
  const [mode, setMode] = useState("signed");
  // Outside a ScientificSelectionProvider the context is undefined rather than an
  // error, so the explorer stays usable on its own; the local id then keeps the matrix,
  // the plot, the detail panel and the inspection cases linked to one selection.
  const selection = useScientificSelection();
  const [localGroupId, setLocalGroupId] = useState(null);
  const selectedGroupId = selection?.selectedComparableGroupId ?? localGroupId;
  const selectedSpeciesId = selection?.selectedComparativeSpeciesId ?? null;

  const selectGroup = (groupOrId, observation = null) => {
    const id = typeof groupOrId === "string"
      ? groupOrId : groupOrId?.comparable_boundary_group_id;
    setLocalGroupId(id || null);
    selection?.selectComparableGroup?.(id, observation);
  };

  const opts = multi?.filter_options || {};
  const view = useMemo(() => filterComparativeDataset(multi, filters), [multi, filters]);
  const groups = view.comparable_boundary_groups;
  const stats = view.distance_statistics;
  const cases = view.inspection_cases;
  const counts = view.counts;

  const selectedGroup = useMemo(
    () => groups.find((g) => g.comparable_boundary_group_id === selectedGroupId) || null,
    [groups, selectedGroupId]);
  const selectedStat = useMemo(
    () => stats.find((s) => s.comparable_boundary_group_id === selectedGroupId) || null,
    [stats, selectedGroupId]);

  const available = Boolean(multi?.available
    && (multi.comparable_boundary_groups || []).length);
  const stem = `${(gene || "gene").toLowerCase()}_comparative_boundaries`;

  // ---- exports: one builder per figure, all from the filtered dataset ------- //
  const buildMatrix = () => comparativeMatrixFigureSpec({
    gene, matrix: view.boundary_matrix, groups, selectedGroupId, mode,
    nearEdgeThreshold: threshold,
  });
  const buildPaired = () => pairedSignedDistanceFigureSpec({
    gene, groups, stats, selectedGroupId, nearEdgeThreshold: threshold,
  });
  const buildConsistency = () => consistencySummaryFigureSpec({
    gene, stats, groups, selectedGroupId,
  });
  const buildArchitecture = () => comparativeArchitectureFigureSpec({
    gene, group: selectedGroup, models,
  });

  const figureExports = [
    ["Comparative matrix", buildMatrix, `${stem}_matrix`],
    ["Paired signed distances", buildPaired, `${stem}_paired_signed_distance`],
    ["Consistency summary", buildConsistency, `${stem}_consistency`],
  ];

  const resetFilters = () => setFilters(EMPTY_FILTERS);
  const nActive = activeFilterCount({ ...EMPTY_FILTERS, ...filters });
  const set = (patch) => setFilters((f) => ({ ...f, ...patch }));

  return (
    <div className="boundary-dashboard-body cbe-explorer">
      <SpeciesInventory rows={multi?.species_rows || []} />

      {!available && (
        <div className="arch-note info">
          <b>Comparative exon–domain-boundary analysis is not available yet.</b> A
          cross-species boundary comparison requires at least two species with
          post-cluster domain annotations and mutually comparable boundaries (an
          MSA-aligned protein position, or a coding-exon group that demonstrably occurs
          in more than one species). Boundaries are never compared only because they
          share an ordinal name such as <i>E3→E4</i>. The species above are shown for
          inventory; no comparative cells are fabricated.
        </div>
      )}

      {available && (
        <>
          {/* ---- filters: one central filtered dataset for every view ------- */}
          <div className="card cbe-filter-card">
            <div className="card-head">
              <h3>Filters</h3>
              <div className="cbe-filter-meta">
                <Badge cls={nActive ? "warn" : "neutral"} soft>
                  {nActive} active filter{nActive === 1 ? "" : "s"}
                </Badge>
                <Badge cls="neutral" soft>{counts.visible_species} species visible</Badge>
                <Badge cls="neutral" soft>
                  {counts.visible_groups} of {counts.total_groups} comparable groups
                </Badge>
                <button className="btn ghost sm" onClick={resetFilters} disabled={!nActive}>
                  Reset all filters
                </button>
              </div>
            </div>
            <div className="cbe-filters">
              <ChipFilter label="Species" selected={filters.species}
                options={(opts.species || []).map((s) => s.species_id)}
                format={(v) => speciesTag(v, (opts.species || [])
                  .find((s) => s.species_id === v)?.scientific_name)}
                onChange={(v) => set({ species: v })} />
              <ChipFilter label="Taxonomic group" selected={filters.taxonomicGroups}
                options={opts.taxonomic_groups || []}
                onChange={(v) => set({ taxonomicGroups: v })} />
              <ChipFilter label="Boundary class" selected={filters.boundaryClasses}
                options={opts.boundary_classes || []}
                format={(v) => CANON_CLASS_LABEL[canonClass(v)]}
                onChange={(v) => set({ boundaryClasses: v })} />
              <ChipFilter label="Representative domain" selected={filters.domainGroups}
                options={(opts.representative_domain_groups || [])
                  .map((d) => d.interpro_accession)}
                format={(v) => (opts.representative_domain_groups || [])
                  .find((d) => d.interpro_accession === v)?.label || v}
                onChange={(v) => set({ domainGroups: v })} />
              <ChipFilter label="Mapping confidence" selected={filters.mappingStatuses}
                options={opts.mapping_statuses || []}
                format={(v) => (v === "tentative" ? "tentative" : "supported")}
                onChange={(v) => set({ mappingStatuses: v })} />
              <ChipFilter label="Domain edge" selected={filters.edges}
                options={opts.edges || []}
                format={(v) => `${v} edge`}
                onChange={(v) => set({ edges: v })} />
              <div className="cbe-filter-group">
                <span className="cbe-filter-label">Scope</span>
                <div className="cbe-chips">
                  <button type="button"
                    className={`cbe-chip${filters.exactNearOnly ? " on" : ""}`}
                    onClick={() => set({ exactNearOnly: !filters.exactNearOnly })}>
                    Exact / near edge only
                  </button>
                  <button type="button"
                    className={`cbe-chip${filters.inspectionOnly ? " on" : ""}`}
                    onClick={() => set({ inspectionOnly: !filters.inspectionOnly })}>
                    Inspection cases only
                  </button>
                  <button type="button"
                    className={`cbe-chip${filters.showUnmapped ? " on" : ""}`}
                    onClick={() => set({ showUnmapped: !filters.showUnmapped })}>
                    Show missing / unmapped
                  </button>
                </div>
              </div>
            </div>
          </div>

          {/* ---- matrix ------------------------------------------------------ */}
          <div className="card">
            <div className="card-head">
              <h3>Species × comparable-boundary groups</h3>
              <div className="cbe-head-actions">
                <div className="seg">
                  {Object.keys(MODE_LABEL).map((m) => (
                    <button key={m} className={`seg-btn${mode === m ? " on" : ""}`}
                      onClick={() => setMode(m)}>{MODE_LABEL[m]}</button>
                  ))}
                </div>
                <Menu label="Export" align="right">
                  <button className="menu-item"
                    onClick={() => downloadFigureSvg(buildMatrix(), `${stem}_matrix`)}>
                    Matrix — SVG (vector)</button>
                  <button className="menu-item"
                    onClick={() => downloadFigurePdf(buildMatrix(), `${stem}_matrix`)}>
                    Matrix — PDF (vector)</button>
                  <button className="menu-item"
                    onClick={() => downloadFigurePng(buildMatrix(), `${stem}_matrix`)}>
                    Matrix — PNG (300 dpi)</button>
                  <button className="menu-item" onClick={() => downloadFigureTsv(
                    comparativeMatrixTsv(view.boundary_matrix, groups, mode),
                    `${stem}_matrix`)}>
                    Matrix — TSV</button>
                </Menu>
              </div>
            </div>
            <ComparativeBoundaryMatrix matrix={view.boundary_matrix} groups={groups}
              mode={mode} selectedGroupId={selectedGroupId}
              selectedSpeciesId={selectedSpeciesId} threshold={threshold}
              onSelectCell={(observation, cell) => selectGroup(
                cell?.comparable_boundary_group_id, observation)} />
          </div>

          {/* ---- linked quantitative plot ------------------------------------ */}
          <div className="card">
            <div className="card-head">
              <h3>Signed distance to nearest domain edge, per species</h3>
              <div className="cbe-head-actions">
                <span className="muted sm">
                  {counts.visible_observations} of {counts.total_observations} observations
                </span>
                <Menu label="Export" align="right">
                  <button className="menu-item" onClick={() => downloadFigureSvg(
                    buildPaired(), `${stem}_paired_signed_distance`)}>
                    Paired signed distance — SVG (vector)</button>
                  <button className="menu-item" onClick={() => downloadFigurePdf(
                    buildPaired(), `${stem}_paired_signed_distance`)}>
                    Paired signed distance — PDF (vector)</button>
                  <button className="menu-item" onClick={() => downloadFigurePng(
                    buildPaired(), `${stem}_paired_signed_distance`)}>
                    Paired signed distance — PNG (300 dpi)</button>
                  <button className="menu-item" onClick={() => downloadFigureTsv(
                    comparativeLongTsv(groups), `${stem}_long`)}>
                    Long boundary table — TSV</button>
                </Menu>
              </div>
            </div>
            <ComparativePairedPlot groups={groups} stats={stats} threshold={threshold}
              selectedGroupId={selectedGroupId}
              onSelectGroup={(id, obs) => selectGroup(id, obs)} />
          </div>

          {/* ---- selected group detail + comparative architecture ----------- */}
          <div className="card">
            <div className="card-head">
              <h3>Selected comparable boundary</h3>
              {selectedGroup && (
                <Menu label="Export" align="right">
                  <button className="menu-item" onClick={() => downloadFigureSvg(
                    buildArchitecture(), `${stem}_architecture_${selectedGroupId}`)}>
                    Comparative architecture — SVG (vector)</button>
                  <button className="menu-item" onClick={() => downloadFigurePdf(
                    buildArchitecture(), `${stem}_architecture_${selectedGroupId}`)}>
                    Comparative architecture — PDF (vector)</button>
                  <button className="menu-item" onClick={() => downloadFigurePng(
                    buildArchitecture(), `${stem}_architecture_${selectedGroupId}`)}>
                    Comparative architecture — PNG (300 dpi)</button>
                </Menu>
              )}
            </div>
            <SelectedGroupDetail group={selectedGroup} stat={selectedStat}
              openGene={openGene} setPage={setPage}
              onOpenSpecies={(sid, hint) => selection?.openSpeciesView?.(sid, hint)} />
            <ComparativeArchitecture group={selectedGroup} models={models} />
          </div>

          {/* ---- consistency summary ---------------------------------------- */}
          <div className="card">
            <div className="card-head">
              <h3>Boundary-position consistency</h3>
              <div className="cbe-head-actions">
                <span className="muted sm">
                  annotation-supported comparability · cautious summary
                </span>
                <Menu label="Export" align="right">
                  <button className="menu-item" onClick={() => downloadFigureSvg(
                    buildConsistency(), `${stem}_consistency`)}>
                    Consistency summary — SVG (vector)</button>
                  <button className="menu-item" onClick={() => downloadFigurePdf(
                    buildConsistency(), `${stem}_consistency`)}>
                    Consistency summary — PDF (vector)</button>
                  <button className="menu-item" onClick={() => downloadFigurePng(
                    buildConsistency(), `${stem}_consistency`)}>
                    Consistency summary — PNG (300 dpi)</button>
                  <button className="menu-item" onClick={() => downloadFigureTsv(
                    comparableMappingTsv(groups, stats), `${stem}_mapping`)}>
                    Comparable-boundary mapping — TSV</button>
                </Menu>
              </div>
            </div>
            <div className="table-scroll"><table className="mini-tbl">
              <thead><tr>
                <th>Group</th><th>Species coverage</th><th>Mapping coverage</th>
                <th>Exact / near</th><th>Raw signed distances</th>
                <th>Cross-species difference</th><th>Dominant class</th>
                <th>Domain annotation</th><th>Mapping</th>
              </tr></thead>
              <tbody>{stats.map((s) => (
                <tr key={s.comparable_boundary_group_id}
                  className={s.comparable_boundary_group_id === selectedGroupId
                    ? "row-selected" : ""}
                  onClick={() => selectGroup(s.comparable_boundary_group_id)}
                  style={{ cursor: "pointer" }}>
                  <td>{s.comparable_boundary_group_id}</td>
                  <td>{s.species_with_mapped_boundary} / {s.n_species_available}</td>
                  <td>{Math.round((s.mapping_coverage || 0) * 100)}%</td>
                  <td>{Math.round((s.exact_or_near_proportion || 0) * 100)}%</td>
                  <td className="sm">{(s.raw_signed_distances || []).map((r) => (
                    `${speciesTag(r.species_id).split(" ")[0]} `
                    + `${r.signed_distance > 0 ? "+" : ""}${r.signed_distance}`
                  )).join(" · ")}</td>
                  <td>{s.cross_species_difference ?? "—"} aa</td>
                  <td>{s.dominant_class
                    ? CANON_CLASS_LABEL[canonClass(s.dominant_class)] : "—"}
                    {s.classes_differ && <span className="muted sm"> · differs</span>}</td>
                  <td>{s.domain_annotation_available_in_all
                    ? "all species" : <span className="muted">incomplete</span>}</td>
                  <td className="sm">{isSupported(s.mapping_status)
                    ? "supported" : "tentative"}</td>
                </tr>
              ))}</tbody>
            </table></div>
            <div className="muted sm">
              {stats.some((s) => s.primary_statistic === "raw_pair")
                ? "With two species the two raw observations and their difference are the "
                  + "result; median and range are reported for completeness but are not "
                  + "distribution statistics at n = 2."
                : "Median and range summarise the species' observations; the raw values "
                  + "remain in the exported long table."}
              {" "}A consistent boundary position is not by itself evidence of
              evolutionary conservation.
            </div>
          </div>

          {/* ---- inspection cases ------------------------------------------- */}
          {Boolean(cases.length) && (
            <div className="card">
              <div className="card-head">
                <h3>Inspection cases</h3>
                <Badge cls="neutral" soft>{cases.length} case{cases.length === 1 ? "" : "s"}</Badge>
              </div>
              <ul className="cbe-cases">
                {cases.map((c) => (
                  <li key={c.case_id}
                    className={c.comparable_boundary_group_id === selectedGroupId
                      ? "sel" : ""}>
                    <button type="button" className="cbe-case-btn"
                      onClick={() => selectGroup(c.comparable_boundary_group_id)}>
                      <span className="cbe-case-top">
                        <Badge cls={CASE_SEVERITY_CLS[c.severity] || "neutral"} soft>
                          {c.label}
                        </Badge>
                        <code>{c.comparable_boundary_group_id}</code>
                      </span>
                      <span className="cbe-case-detail">{c.detail}</span>
                    </button>
                  </li>
                ))}
              </ul>
              <div className="muted sm">
                A discrepancy between species may be biological, an annotation gap, or an
                alignment artefact. These cases mark positions worth checking; they are
                not labelled as errors.
              </div>
            </div>
          )}

          {/* ---- all exports in one place ----------------------------------- */}
          <div className="em-export">
            <span className="muted small">Export all:</span>
            {figureExports.map(([label, build, name]) => (
              <Menu key={name} label={label} align="left">
                <button className="menu-item"
                  onClick={() => downloadFigureSvg(build(), name)}>SVG (vector)</button>
                <button className="menu-item"
                  onClick={() => downloadFigurePdf(build(), name)}>PDF (vector)</button>
                <button className="menu-item"
                  onClick={() => downloadFigurePng(build(), name)}>PNG (300 dpi)</button>
              </Menu>
            ))}
            <button className="btn ghost sm" onClick={() => downloadFigureTsv(
              comparativeLongTsv(groups), `${stem}_long`)}>Long boundary TSV</button>
            <button className="btn ghost sm" onClick={() => downloadFigureTsv(
              comparativeMatrixTsv(view.boundary_matrix, groups, mode), `${stem}_matrix`)}>
              Matrix TSV</button>
            <button className="btn ghost sm" onClick={() => downloadFigureTsv(
              comparableMappingTsv(groups, stats), `${stem}_mapping`)}>
              Comparable-boundary mapping TSV</button>
          </div>
        </>
      )}
    </div>
  );
}
