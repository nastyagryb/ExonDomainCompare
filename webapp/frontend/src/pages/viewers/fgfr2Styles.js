// The frozen FGFR2 visual vocabulary.
//
// FGFR2 is the immutable scientific regression target of this project, and its
// figures encode a different set of features than the generic single-species views:
// named architecture classes (Ig-like, kinase, signal peptide, FGFR family), the
// pyTMHMM receptor topology, and the IIIb/IIIc cassette distinction. None of these
// exist in the generic vocabulary, so they are NOT semantic keys in
// semanticStyles.js — folding them in would either lose the Ig-like/kinase
// distinction or silently repaint the frozen figures.
//
// They live here rather than inside the components for the same reason as the
// generic spec: a colour that is written twice is a colour that will drift. This
// module previously existed as duplicated constants in DomainArchitecture.jsx and
// BoundaryDetailTrack.jsx.
//
// Known and accepted collision: some of these hues are the same Okabe–Ito values the
// generic spec uses for other meanings — #009E73 is "kinase domain" here and
// "inside domain" there. The two never appear in one figure or one legend, and the
// FGFR2 side is frozen, so the collision is documented rather than resolved. See
// docs/architecture/figure_parity_contract.md.

/** Architecture classes of the FGFR2 domain track. */
export const FGFR2_DOMAIN_FILL = {
  ig_like_domain: "#56B4E9",
  kinase_domain: "#009E73",
  signal_peptide: "#CC79A7",
  other_domain: "#D9DCE1",
};

export const FGFR2_DOMAIN_LABEL = {
  ig_like_domain: "Ig-like",
  kinase_domain: "Kinase",
  signal_peptide: "Signal peptide",
  other_domain: "FGFR family",
};

/** A domain block whose class is not one of the named FGFR2 classes. */
export const FGFR2_DOMAIN_FALLBACK = "#B0B4BB";

/** pyTMHMM topology states of the receptor. */
export const FGFR2_TM_FILL = {
  receptor_tm: "#E69F00",
  n_terminal_signal_anchor: "#F6C36B",
};

/** The IIIb / IIIc alternative-cassette distinction. Mirrors --iiib / --iiic. */
export const FGFR2_CASSETTE_FILL = { iiib: "#138a9c", iiic: "#d2622a" };

/** Functional-site and disorder markers of the generic architecture track. */
export const SITE_MARKER_FILL = "#D55E00";

/** Ink for a label drawn on top of a filled architecture block. */
export const FGFR2_BLOCK_LABEL_INK = "#12151a";

/** Hairline outline that separates adjacent architecture blocks. */
export const FGFR2_BLOCK_OUTLINE = "#2b2f36";
