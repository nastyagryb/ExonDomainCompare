# Changelog

## 0.2.1 - 2026-08-22

- Corrected structure-matched permutation binning and added a second common-bin sensitivity specification.
- Added a deterministic, outcome-independent audit for selecting the 17 structure-mapped positions from the 25 FDR-significant columns.
- Recomputed and separated the valid leave-one-clade-out portability analysis from an invalid ancestral reconstruction.
- Quarantined the non-identifiable ancestral reconstruction and removed it from active claims and release checks.
- Recorded fixed random seeds for all IQ-TREE topology and AU-test stages.
- Replaced dependent protein-level binomial intervals with descriptive protein- and species-level topology counts.
- Documented the implemented lexicographic NCBI/Ensembl candidate-ranking rule.
- Corrected citation metadata and synchronized the software version.

## 0.2.0 - 2026-08-11

- Added portable local setup and startup commands.
- Added private per-user configuration and LRZ capability checks.
- Added user-scoped InterProScan and pyTMHMM installation support.
- Added checksummed FGFR2 and PTPN11 release datasets.
- Consolidated the source layout and retained supported compatibility commands.
- Unified Candidate display behavior across interactive views and exports.
- Preserved short protein-isoform indels as exploratory evidence.
