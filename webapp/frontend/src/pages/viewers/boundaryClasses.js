// Canonical (generic) exon-boundary classification vocabulary.
//
// Mirrors scripts/shared_gene_analysis/boundary_classification.py. Colours and
// labels are NOT defined here: they come from the shared scientific visual
// specification in semanticStyles.js, so the interactive Boundary views and the
// exported publication figures cannot drift apart. They previously did: this file
// used to carry its own green/orange/violet palette, which the exported figures
// never used, so a class had one colour on screen and another in the paper.
//
// Kept separate from the React components so fast-refresh's
// only-export-components rule stays satisfied. Never used for the frozen FGFR2
// Boundary Consistency vocabulary, which is a different scientific vocabulary with
// different thresholds and lives in boundary.js.

import { BOUNDARY_CLASS_LABEL, FEATURE_STYLES, boundaryStyleKey } from "./semanticStyles.js";

export const CANON_CLASS_ORDER = [
  "exact_domain_edge", "near_domain_edge", "inside_domain",
  "outside_annotated_domains", "unavailable_or_uncertain",
];

export const CANON_CLASS_COLOR = Object.fromEntries(
  CANON_CLASS_ORDER.map((c) => [c, FEATURE_STYLES[boundaryStyleKey(c)].fill]));

export const CANON_CLASS_LABEL = { ...BOUNDARY_CLASS_LABEL };

const LEGACY_ALIAS = {
  exact_edge: "exact_domain_edge",
  near_edge: "near_domain_edge",
  outside_domain: "outside_annotated_domains",
  unknown: "unavailable_or_uncertain",
};

export function canonClass(c) {
  if (!c) return "unavailable_or_uncertain";
  if (CANON_CLASS_COLOR[c]) return c;
  return LEGACY_ALIAS[c] || "unavailable_or_uncertain";
}
