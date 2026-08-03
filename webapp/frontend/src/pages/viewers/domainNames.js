/**
 * Readable display names for real InterPro short names.
 *
 * Single source of truth for the interactive views and the publication figure
 * renderer. Kept free of React and of any other import so the headless Node
 * renderer can load it too.
 *
 * The Python counterpart is `_DOMAIN_DISPLAY_NAMES` in
 * `src/exondomaincompare/shared_gene_analysis/protein_coordinate_model.py`; both must agree, so
 * that a figure, a table and an exported TSV name the same feature identically.
 *
 * This never invents a name: the accession and the raw short name travel with
 * every feature, and only the text shown to a reader is made readable.
 */
const DISPLAY_NAMES = {
  "Ig-like_dom": "Ig-like domain",
  "Ig-like_fold": "Ig-like fold",
  "Ig-like_dom_sf": "Ig-like domain superfamily",
  "Ser-Thr/Tyr_kinase_cat_dom": "Ser-Thr/Tyr kinase domain",
  "Kinase-like_dom_sf": "Kinase-like domain superfamily",
  Protein_kinase_ATP_BS: "Protein kinase ATP binding site",
  Tyr_kinase_AS: "Tyrosine kinase active site",
  FGF_rcpt_fam: "FGF receptor family",
  RTK: "Receptor tyrosine kinase family",
  disorder_prediction: "Predicted disorder",
};

export function prettyDomainName(name) {
  if (!name) return "Unnamed feature";
  const raw = String(name);
  if (DISPLAY_NAMES[raw]) return DISPLAY_NAMES[raw];
  const words = raw.replace(/_/g, " ").replace(/\bdom\b/g, "domain").trim();
  return words.charAt(0).toUpperCase() + words.slice(1);
}
