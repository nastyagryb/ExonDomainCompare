# Figure parity contract

The same scientific figure exists twice: once as an interactive React viewer, once as
an exported SVG/PDF/PNG for publication. This document states what the two are
required to share, what they are allowed to differ in, and which of those
requirements are enforced by a test rather than by care.

The contract exists because the two sides had genuinely drifted. Before it, the
interactive exon block was `#c9d3e2` while the exported one was `#A9BED4`, and the
five generic boundary classes carried one palette on screen and a different one in
the figures. A reader comparing the screen with the paper was comparing two
different colour systems, and nothing failed to warn anyone.

## The single source of truth

`webapp/frontend/src/pages/viewers/semanticStyles.js` defines one semantic style per
scientific feature kind. Every scientific colour, stroke, marker and label priority
comes from there — the interactive viewers and the publication renderer both read the
same table. Colours are the colour-blind-safe publication palette from
`figureSpec.js`; the interactive side adopted them, not the reverse, so the accepted
exported figures stayed byte-identical when the contract was introduced.

The 19 required semantic keys are: `coding_exon`, `alternative_exon`, `shared_exon`,
`shifted_boundary`, `representative_domain`, `family_superfamily`, `tm_helix`,
`candidate_region`, `boundary_exact`, `boundary_near`, `boundary_inside`,
`boundary_outside`, `boundary_uncertain`, `selected_feature`, `primary_sequence`,
`alternative_sequence`, `gap`, `variable_region`, `conserved_region`.

Five more were added because the viewers draw them and leaving them as literals would
have reopened the drift this module closes: `protein_backbone`, `exon_boundary_tick`,
`member_signature`, `functional_site`, `disorder_region`. A `CHROME` group covers
paper, grid and rules, which carry no scientific meaning but still have to be explicit,
because a grid line that vanishes in an export costs the reader the coordinate frame.

Each provides `fill`, `stroke`, `strokeWidth`, `opacity`, `text`, `marker`,
`labelPriority` and `printFallback`. Asking for an unknown key raises instead of
returning a default, because a silent default is precisely how a feature ends up
painted black in an exported SVG.

## What both sides must share

| Shared property | How it is held |
|---|---|
| Semantic style identity | both call `featureProps` / `featureStyle` from `semanticStyles.js`; no component may hardcode a publication colour |
| Boundary-class identity | `boundaryClasses.js` derives its colours and labels from the shared spec; `mainFigures.js` resolves classes through the same `boundaryStyleKey` |
| Domain-instance identity | both resolve an instance by `domain_instance_id` and colour it by `display_order` through `domainInstanceFill` |
| Candidate identity | C-labels come from `candidate_regions` in the coordinate model, ordered by `aa_start`, never from a per-figure ordering |
| Label text | domain names from `domainNames.js`, boundary-class names from `BOUNDARY_CLASS_LABEL`, exon names as `E<n>`, boundary names as `E<n> → E<n+1>` |
| Legend vocabulary | the legend of a class figure lists exactly the `BOUNDARY_CLASS_LABEL` values |
| Feature ordering | both read the coordinate model's order: representative domains, family/superfamily, TM topology, coding exons, boundaries, candidates |
| Coordinate transformation | both map amino-acid position 1…protein_length onto the drawing width; feature coordinates are taken from the model, never recomputed |
| Selected-feature meaning | a selection is an outline, never a recolouring and never a translucent slab across tracks, so the feature keeps the class colour the reader is judging |

## What they are allowed to differ in

Interaction is not a figure property. The interactive side may add, and the export
side may omit: hover state, focus outlines, selection handles, tooltips, responsive
scrolling and zoom, control widgets, and CSS-driven layout, spacing and transitions.

The two also legitimately differ in rendering technology: the interactive views draw
React SVG elements, the export path builds a declarative figure spec that is rendered
to SVG and to PDF. They are not required to produce identical DOM trees, and a test
that demanded that would fail for reasons that do not matter.

## Deliberately outside the contract

FGFR2 is the immutable scientific regression target, and its figures encode features
the generic views do not have. Two modules therefore keep their own palettes on
purpose:

- `boundary.js` — the **frozen FGFR2 Boundary Consistency vocabulary**
  (`aligned_to_domain_boundary`, `near_domain_boundary`, `within_domain`,
  `between_domains`, `review_or_missing`) with its own thresholds of 0–3 aa and
  4–15 aa. That is a different scientific classification from the generic five-class
  one, not a second styling of the same thing.
- `fgfr2Styles.js` — the **frozen FGFR2 architecture vocabulary**: the named domain
  classes (Ig-like, kinase, signal peptide, FGFR family), the pyTMHMM receptor
  topology states, and the IIIb/IIIc cassette distinction. Folding these into
  `representative_domain` would make Ig-like and kinase indistinguishable.

Being outside the shared spec does not license scattering: both are single modules,
and a test asserts that no component redefines the cassette colours locally.

One collision is documented rather than resolved. Some FGFR2 hues are the same
Okabe–Ito values the generic spec uses for other meanings — `#009E73` is "kinase
domain" in the FGFR2 architecture and "inside domain" in the generic boundary views.
The two never share a figure or a legend, and the FGFR2 side is frozen, so renaming
either would change an accepted result for no reader's benefit.

## What is enforced, and what is still convention

Enforced by `tests/test_figure_parity_contract.py`:

- the shared spec defines all 19 semantic keys, each with all eight properties;
- no viewer component hardcodes a colour that the shared spec owns;
- `boundaryClasses.js` carries no palette of its own, and its labels are the shared
  ones;
- `mainFigures.js` resolves domain and boundary colours through the shared spec;
- the count of CSS-dependent scientific marks is zero, so no scientific mark can go
  black when the stylesheet is absent;
- for each main figure, the exported SVG's boundary-class colours, class legend,
  domain-instance labels and candidate labels equal the values the shared spec and
  the coordinate model prescribe.

Still convention, and the honest limit of this contract: the **pixel geometry** of
the interactive views is computed by the components themselves, not by the shared
figure spec. Both sides read the same coordinates from the model, so a feature cannot
sit at a different amino-acid position in the two, but the mapping from amino acids to
pixels is written twice. Porting the interactive views onto `figureSpec.js` would make
that last piece structural too; it is the recommended next step and is not done here.
