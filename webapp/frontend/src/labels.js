// Gene/event label abstraction (additive, low-risk).
//
// The app is validated for FGFR2 IIIb/IIIc. This helper lets display labels be
// driven by the active dataset's gene/event config (backend `ui_labels`, or a
// generic dataset_summary), while ALWAYS falling back to the current FGFR2
// wording. For FGFR2 the UI therefore looks exactly as before.
//
// Internally, components can think in generic terms (event_region, boundary_relation),
// but render whatever label the config provides.

export const DEFAULT_DATASET_LABELS = {
  geneExplorer: "Gene Explorer",
  eventRegion: "Cassette",
  eventRegionFull: "IIIb/IIIc cassette",
  eventDiscriminatingColumns: "IIIb/IIIc-discriminating columns",
  boundaryRelation: "Boundary Consistency",
  domainRelationDescription: "Cassette-to-domain boundary consistency",
  referenceComparison: "Human comparison",
  geneSymbol: "FGFR2",
  eventDisplayName: "IIIb/IIIc cassette",
};

// Accepts a backend dataset status (with `ui_labels` + gene/event fields) or a
// generic dataset_summary. Returns a stable camelCase label object.
export function getDatasetLabels(source) {
  const s = source || {};
  const ui = s.ui_labels || {};
  const pick = (snakeKey, fallback) =>
    (ui[snakeKey] != null && ui[snakeKey] !== "" ? ui[snakeKey] : fallback);
  return {
    geneExplorer: pick("gene_explorer", DEFAULT_DATASET_LABELS.geneExplorer),
    eventRegion: pick("event_region", DEFAULT_DATASET_LABELS.eventRegion),
    eventRegionFull: pick("event_region_full", DEFAULT_DATASET_LABELS.eventRegionFull),
    eventDiscriminatingColumns: pick(
      "event_discriminating_columns", DEFAULT_DATASET_LABELS.eventDiscriminatingColumns),
    boundaryRelation: pick("boundary_relation", DEFAULT_DATASET_LABELS.boundaryRelation),
    domainRelationDescription: pick(
      "domain_relation_description", DEFAULT_DATASET_LABELS.domainRelationDescription),
    referenceComparison: pick("reference_comparison", DEFAULT_DATASET_LABELS.referenceComparison),
    geneSymbol: s.gene_symbol || DEFAULT_DATASET_LABELS.geneSymbol,
    eventDisplayName: s.event_display_name || DEFAULT_DATASET_LABELS.eventDisplayName,
  };
}
