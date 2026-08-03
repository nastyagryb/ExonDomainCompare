// One filtered comparative dataset for the whole Comparative Boundary Explorer.
//
// The matrix, the paired plot, the detail panel, the consistency summary and every
// exported table read the result of this single function. That is deliberate: when each
// view filtered its own copy, a figure and its own source table could disagree about
// which observations were included, and the reader had no way to notice.
//
// Filtering only ever removes real observations published by the backend. Nothing is
// re-derived here — group membership, mapping method, mapping status and confidence all
// come from src/exondomaincompare/shared_gene_analysis/boundary_dashboard.py.

import { canonClass } from "./boundaryClasses.js";
import { isSupported } from "./comparativeFigures.js";

export const EMPTY_FILTERS = {
  species: [],            // species_id allow-list; [] = all
  taxonomicGroups: [],    // taxonomic group allow-list; [] = all
  boundaryClasses: [],    // canonical class allow-list; [] = all
  domainGroups: [],       // InterPro accessions of the nearest representative domain
  mappingStatuses: [],    // supported_comparable / high_confidence_comparable / tentative
  edges: [],              // "start" / "end"
  exactNearOnly: false,   // keep only exact_domain_edge / near_domain_edge observations
  inspectionOnly: false,  // keep only groups that raised an inspection case
  showUnmapped: true,     // show species rows/cells without an observation
};

const has = (list, v) => !list || list.length === 0 || list.includes(v);

/** Number of filters the user has actively changed, for the "n active" badge. */
export function activeFilterCount(f) {
  let n = 0;
  if (f.species.length) n += 1;
  if (f.taxonomicGroups.length) n += 1;
  if (f.boundaryClasses.length) n += 1;
  if (f.domainGroups.length) n += 1;
  if (f.mappingStatuses.length) n += 1;
  if (f.edges.length) n += 1;
  if (f.exactNearOnly) n += 1;
  if (f.inspectionOnly) n += 1;
  if (!f.showUnmapped) n += 1;
  return n;
}

const observationPasses = (o, f) => has(f.species, o.species_id)
  && has(f.taxonomicGroups, o.taxonomic_group)
  && has(f.boundaryClasses, canonClass(o.boundary_class))
  && has(f.domainGroups, o.nearest_domain_accession)
  && has(f.edges, o.nearest_edge)
  && (!f.exactNearOnly
    || ["exact_domain_edge", "near_domain_edge"].includes(canonClass(o.boundary_class)));

/**
 * Apply the filter state to the canonical comparative index.
 *
 * @param multi the backend's boundary_dashboard.multi_species object
 * @param filters see EMPTY_FILTERS
 * @returns a filtered view with the same shape as `multi`, plus counts
 *
 * A group survives only if at least one of its observations passes. Groups that keep a
 * single observation are still shown — losing sight of them would hide exactly the
 * asymmetric cases that matter — but they are marked so no view can present them as a
 * confirmed cross-species pair.
 */
export function filterComparativeDataset(multi, filters = EMPTY_FILTERS) {
  const f = { ...EMPTY_FILTERS, ...filters };
  const allGroups = multi?.comparable_boundary_groups || [];
  const allStats = multi?.distance_statistics || [];
  const allCases = multi?.inspection_cases || [];
  const allMatrix = multi?.boundary_matrix || [];
  const speciesRows = multi?.species_rows || [];

  const caseGroupIds = new Set(allCases.map((c) => c.comparable_boundary_group_id));

  const groups = [];
  for (const g of allGroups) {
    if (!has(f.mappingStatuses, g.mapping_status)) continue;
    if (f.inspectionOnly && !caseGroupIds.has(g.comparable_boundary_group_id)) continue;
    const kept = (g.per_species_native_positions || []).filter((o) => observationPasses(o, f));
    if (!kept.length) continue;
    groups.push({
      ...g,
      per_species_native_positions: kept,
      n_observations_filtered_out:
        (g.per_species_native_positions || []).length - kept.length,
      // A pair claim requires two surviving observations *and* a supported mapping.
      // Filtering down to one observation must never leave a view free to draw a
      // connector, so the answer is computed once, here.
      connectable: kept.length >= 2 && isSupported(g.mapping_status),
    });
  }

  const keptIds = new Set(groups.map((g) => g.comparable_boundary_group_id));
  const keptObs = new Map();
  for (const g of groups) {
    for (const o of g.per_species_native_positions) {
      keptObs.set(`${g.comparable_boundary_group_id}::${o.species_id}`, o);
    }
  }

  const visibleSpecies = new Set();
  for (const g of groups) {
    for (const o of g.per_species_native_positions) visibleSpecies.add(o.species_id);
  }
  // An explicit species filter defines the rows even where a species has no surviving
  // observation, so an empty row stays visible as an absence rather than vanishing.
  const rowSpecies = f.species.length
    ? f.species
    : speciesRows.filter((r) => has(f.taxonomicGroups, r.taxonomic_group))
      .map((r) => r.species_id);

  const matrix = allMatrix
    .filter((r) => rowSpecies.includes(r.species_id))
    .map((r) => ({
      ...r,
      cells: (r.cells || [])
        .filter((c) => keptIds.has(c.comparable_boundary_group_id))
        .map((c) => {
          const survivor = keptObs.get(
            `${c.comparable_boundary_group_id}::${r.species_id}`);
          if (survivor) return c;
          // The observation existed but was filtered out. Reporting it as
          // "boundary_absent_or_unmapped" would claim the backend found nothing, so
          // it gets its own state and no numbers.
          if (c.observed) {
            return {
              ...c, observed: false, state: "filtered_out",
              native_position: null, signed_distance: null,
              absolute_distance: null, observation: null,
            };
          }
          return c;
        })
        .filter((c) => f.showUnmapped || c.observed),
    }))
    .filter((r) => f.showUnmapped || r.cells.some((c) => c.observed));

  const stats = allStats.filter((s) => keptIds.has(s.comparable_boundary_group_id));
  const cases = allCases.filter((c) => keptIds.has(c.comparable_boundary_group_id));

  return {
    ...multi,
    comparable_boundary_groups: groups,
    boundary_matrix: matrix,
    distance_statistics: stats,
    inspection_cases: cases,
    species_rows: speciesRows.filter((r) => rowSpecies.includes(r.species_id)),
    counts: {
      visible_species: (f.showUnmapped ? rowSpecies.length : visibleSpecies.size),
      visible_groups: groups.length,
      total_groups: allGroups.length,
      visible_observations: groups.reduce(
        (n, g) => n + g.per_species_native_positions.length, 0),
      total_observations: allGroups.reduce(
        (n, g) => n + (g.per_species_native_positions || []).length, 0),
      visible_cases: cases.length,
      active_filters: activeFilterCount(f),
    },
  };
}
