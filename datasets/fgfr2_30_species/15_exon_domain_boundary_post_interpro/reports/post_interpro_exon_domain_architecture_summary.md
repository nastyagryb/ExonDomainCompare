# Post-InterProScan + pyTMHMM FGFR2 exon-domain architecture — summary

This is a **post-InterProScan / pyTMHMM visualization + QC** step. It does **not** modify the pre-InterPro truth table, FASTA files, or primary/review membership. InterProScan provides domain annotations (Ig-like, kinase, …), pyTMHMM provides the transmembrane-helix prediction, and the final IIIb/IIIc labels still come from the final truth table. InterProScan/pyTMHMM do **not** validate or relabel IIIb/IIIc.

## Input files

- InterProScan TSV: `results/final_30_until_interpro_prepare/14_interproscan_primary/input.fasta.tsv`
- InterProScan folder: `results/final_30_until_interpro_prepare/14_interproscan_primary`
- pyTMHMM predictions: `results/final_30_until_interpro_prepare/15_pytmhmm_primary/pytmhmm_transmembrane_hits.tsv`
- Truth table: `results/final_30_until_interpro_prepare/13_final_pre_interpro_closure/final_pre_interpro_truth_table.tsv`
- Sequence manifest: `results/final_30_until_interpro_prepare/13_final_pre_interpro_closure/freeze/final_pre_interpro_sequence_manifest.tsv`
- Coordinate / cassette table: `results/final_30_until_interpro_prepare/13_final_pre_interpro_closure/tables/figure3C_exon_to_protein_cassette_coordinate_map.tsv`
- Member databases present: CDD, FunFam, Gene3D, MobiDBLite, PANTHER, PIRSF, PRINTS, Pfam, ProSitePatterns, ProSiteProfiles, SMART, SUPERFAMILY

## Counts (primary proteins)

- Proteins parsed (primary, in truth table): **58**
- With a pyTMHMM TM prediction (≥1 helix): **58**
- With a resolved receptor TM helix (used for QC/plots): **58**
- With ≥1 InterProScan domain match: **58**
- With a kinase domain detected: **58**
- With ≥1 Ig-like domain detected: **58**

## QC status distribution

| final_qc_status | count |
| --- | --- |
| architecture_supported | 31 |
| architecture_supported_with_minor_flags | 27 |
| interpro_partial_tmhmm_supported | 0 |
| review_unusual_domain_order | 0 |
| insufficient_domain_or_tm_data | 0 |
| failed_coordinate_mapping | 0 |

- architecture_supported: **31**
- architecture_supported_with_minor_flags: **27**
- review / unusual architecture: **0**

## Proteins with unusual / review domain order

_None — no protein showed a biologically impossible domain order._

## Notes

- **InterProScan** annotates protein domains (Ig-like, kinase, …); **pyTMHMM** predicts the transmembrane helix; the **final truth table** owns the IIIb/IIIc labels. Neither InterProScan nor pyTMHMM relabels IIIb/IIIc.
- pyTMHMM is used as the authoritative TM layer (InterProScan did not annotate the TM helix). Where InterProScan also reports a TM, agreement is recorded in the QC table (`tm_agreement`).
- FGFR2 is a single-pass type-I receptor. pyTMHMM often reports a second, N-terminal helix; this is treated as a signal-anchor (N-terminal signal region) and shown separately, not as the membrane-spanning receptor TM.
- Missing InterProScan TM is **not** a failure when pyTMHMM predicts the TM helix. A protein is only flagged when the pyTMHMM TM is missing or the domain order is biologically impossible (e.g. cassette slot overlapping the kinase, TM downstream of the kinase).
- Ig1/Ig2/Ig3 are numbered positionally (N→C). Where CDD FGFR-specific Ig signatures (IgI_1/2/3_FGFR) are absent, the numbering is positional only and flagged as a minor note.
- Coding exons and IIIb/IIIc cassette slots come from the existing coordinate table (numbered N→C), not from InterProScan/pyTMHMM. For major-native-offset species the reference `cassette_start/end_aa` can lie outside the species protein coordinate space; the drawn cassette slot therefore uses the protein-mapped cassette CDS block coordinates.