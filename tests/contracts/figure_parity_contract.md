# Figure parity contract

The same scientific figure exists twice: once as an interactive React viewer, once as
an exported SVG/PDF/PNG for publication. This contract states what both sides must
share, what they are allowed to differ in, and which requirements are enforced by
tests.

## The single source of truth

`webapp/frontend/src/pages/viewers/semanticStyles.js` defines one semantic style per
scientific feature kind. Interactive viewers and publication renderers both read this
table. Unknown style keys raise an error instead of silently using a default.

The required semantic keys cover coding and alternative exons, domain layers,
transmembrane regions, Candidate regions, boundary classes, sequence roles, gaps,
variable regions, conserved regions and selection outlines. Each style provides its
fill, stroke, stroke width, opacity, text colour, marker, label priority and print
fallback.

## What both sides must share

| Shared property | Requirement |
| --- | --- |
| Semantic style identity | Both use `semanticStyles.js`; components do not redefine publication colours. |
| Boundary-class identity | Interactive and exported figures resolve the same class labels and colours. |
| Domain-instance identity | Both resolve domains by `domain_instance_id` and `display_order`. |
| Candidate identity | Candidate labels come from the coordinate model and use the same ordering. |
| Label text | Domain, boundary, exon and Candidate labels use the same vocabulary. |
| Legend vocabulary | Class legends list the shared boundary-class labels. |
| Feature ordering | Both read the order encoded in the coordinate model. |
| Coordinates | Both use the model's amino-acid coordinates without recomputing biology. |
| Selection | Selection is an outline and does not replace the scientific class colour. |

## What they are allowed to differ in

The interactive side may add hover state, focus outlines, tooltips, responsive
scrolling, zoom and controls. The export side may omit these interaction details.
React SVG and publication exporters may use different markup and layout technology;
identical DOM trees and pixel geometry are not required.

## Deliberately outside the contract

FGFR2 is the immutable scientific regression target and contains features the generic
views do not have. `boundary.js` keeps the frozen FGFR2 Boundary Consistency
vocabulary, while `fgfr2Styles.js` keeps the frozen FGFR2 domain, topology and
IIIb/IIIc cassette vocabulary. These remain centralised exceptions rather than local
component styles.

## What is enforced, and what is still convention

The test suite enforces complete semantic styles, shared class vocabularies, central
colour ownership, export-safe marks, matching labels and matching Candidate identity.
Pixel geometry remains convention: both sides read the same biological coordinates,
but their amino-acid-to-pixel mapping is implemented by their respective renderers.
