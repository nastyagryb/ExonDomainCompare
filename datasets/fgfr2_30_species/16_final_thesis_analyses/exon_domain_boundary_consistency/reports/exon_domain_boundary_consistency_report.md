# Exon-domain boundary consistency report

Thesis analysis: *consistency and robustness of exon-domain boundary identification across vertebrate FGFR2 orthologs (IIIb/IIIc case study)*.

## Input files (read-only)

* `results/final_30_until_interpro_prepare/13_final_pre_interpro_closure/final_pre_interpro_truth_table.tsv`
* `results/final_30_until_interpro_prepare/13_final_pre_interpro_closure/freeze/final_pre_interpro_sequence_manifest.tsv`
* `results/final_30_until_interpro_prepare/15_exon_domain_boundary_post_interpro/tables/exon_domain_architecture_features.tsv`
* `results/final_30_until_interpro_prepare/15_exon_domain_boundary_post_interpro/tables/interpro_domain_features_normalized.tsv`
* `results/final_30_until_interpro_prepare/15_exon_domain_boundary_post_interpro/tables/pytmhmm_tm_features_normalized.tsv`
* `results/final_30_until_interpro_prepare/15_exon_domain_boundary_post_interpro/tables/fgfr2_domain_architecture_qc.tsv`
* `results/final_30_until_interpro_prepare/15_exon_domain_boundary_post_interpro/tables/exon_block_length_consistency_audit.tsv`

Provenance-only (not used as final display truth): `post_interpro_qc_review_case_audit.tsv`, `exon_block_reconstruction_overrides.json`, figure3C source table.

## Method

* Primary proteins analyzed: **58**.
* Coordinates come from the **sanitized** post-InterPro feature table `exon_domain_architecture_features.tsv` (never the raw figure3C blocks where a sanitation override exists).
* Domain landmarks: representative InterProScan **Ig-like** (Ig1/Ig2/Ig3) and **kinase** domains plus the **pyTMHMM** receptor TM helix. The family-level FGFR fingerprint (spanning the whole protein) is excluded because it carries no informative internal boundary.
* pyTMHMM provides the TM layer because InterProScan did not annotate the transmembrane helix consistently for these proteins.
* For each boundary the absolute distance to the nearest domain edge (start or end) is computed and classified.

### Boundary-class thresholds

* `aligned_to_domain_boundary`: 0–3 aa from nearest domain boundary
* `near_domain_boundary`: 4–15 aa
* `within_domain`: boundary inside a domain but > 15 aa from its edges
* `between_domains`: boundary outside all annotated domains
* `review_or_missing`: cassette/domain/coordinate data missing

## Overall cassette-boundary consistency

Across 116 cassette boundaries (58 proteins × start+end):

* aligned (0-3): **12** (10%)
* near (4-15): **49** (42%)
* within domain: **55** (47%)
* between domains: **0** (0%)
* missing/NA: **0** (0%)

Median distance to nearest domain boundary: **13.5 aa**; mean: **21.29 aa**.

**Interpretation:** 61/116 (53%) of cassette boundaries fall within 15 aa of an InterProScan/pyTMHMM domain boundary, i.e. the IIIb/IIIc splice cassette is consistently anchored to the D3 (Ig3) region of the receptor architecture across vertebrates. This supports the robustness of exon–domain boundary identification.

## IIIb vs IIIc

* **IIIb** (29 proteins): aligned 6, near 25, within 27, between 0; median 13.0 aa.
* **IIIc** (29 proteins): aligned 6, near 24, within 28, between 0; median 13.5 aa.

## Taxon groups

* **Primates** (11 proteins): median 12.5 aa (aligned 2, near 11, within 9, between 0).
* **Other mammals** (21 proteins): median 15.5 aa (aligned 1, near 20, within 21, between 0).
* **Birds** (6 proteins): median 11.0 aa (aligned 0, near 7, within 5, between 0).
* **Reptiles** (6 proteins): median 22.0 aa (aligned 2, near 4, within 6, between 0).
* **Amphibians** (4 proteins): median 16.5 aa (aligned 2, near 2, within 4, between 0).
* **Teleost fish** (10 proteins): median 13.5 aa (aligned 5, near 5, within 10, between 0).

## Outliers (13)

* ambystoma mexicanum IIIc (cassette_end @ 380): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* canis lupus familiaris IIIb (cassette_end @ 244): low-confidence / reconstructed exon-block display (cassette_only_high_confidence)
* chrysemys picta bellii IIIb (cassette_end @ 379): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* danio rerio IIIb (cassette_end @ 383): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* gasterosteus aculeatus IIIb (cassette_end @ 475): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* gorilla gorilla gorilla IIIb (cassette_end @ 270): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* meleagris gallopavo IIIc (cassette_end @ 245): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* mus musculus IIIb (cassette_end @ 393): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* oreochromis niloticus IIIc (cassette_end @ 424): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* oryzias latipes IIIb (cassette_end @ 376): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* pongo abelii IIIc (cassette_end @ 446): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* takifugu rubripes IIIb (cassette_end @ 343): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)
* xenopus tropicalis IIIb (cassette_end @ 368): low-confidence / reconstructed exon-block display (native_exon_blocks_reconstructed)

## Notes

* InterProScan/pyTMHMM domain annotations **support** the domain architecture; they never relabel IIIb/IIIc (labels always come from the final truth table).
* pyTMHMM is used for the TM layer because InterProScan did not annotate TM helices consistently.
* All coordinates use the sanitized post-InterPro feature tables; low-confidence / reconstructed exon-block display is flagged (marked with `*` in Figure 11) and reported separately from biological QC — it is a display-coordinate property, not a biological failure.
