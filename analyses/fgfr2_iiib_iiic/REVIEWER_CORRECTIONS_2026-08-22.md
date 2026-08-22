# Reviewer corrections and thesis-impact register

## 1. Ancestral reconstruction and LOCO separation

The earlier ancestral reconstruction was not repaired by repeating the same numerical calculation. Its target IIIb/IIIc ancestral split was not identifiable from the unrooted tree, so a repeated run would reproduce an unjustified node label. The implementation and outputs are preserved under `quarantined/invalid_ancestral_reconstruction` and excluded from active verification.

The LOCO component was extracted into `scripts/analysis/09_loco_validation.py` and independently rerun with seed 20260804, 10,000 paired label permutations and 1,000 training bootstraps per fold. All six valid LOCO tables reproduce the former LOCO tables exactly. The active result is 56/56 correct held-out proteins with empirical *p* = 9.999e-05. It must be described as cross-clade portability of a sequence-calibrated partition, not independent validation of isoform truth.

Thesis impact: remove the 15/17 ancestral-difference and 11-site ancestral-core claims, figures and conclusion language. Retain LOCO with the narrower interpretation. A valid future ASR requires justified rooting or species-tree reconciliation and root/topology sensitivity analyses.

## 2. Structure-matched permutation correction

The original code derived accessibility bins from control residues but assigned target bins with different, full-universe cut points. The corrected implementation applies the same bin edges to target and control residues within each structure. The primary specification uses control-derived quartiles and the sensitivity specification uses full-universe quartiles.

Direct-contact enrichment is retained in both structures under both specifications. For 1NUN, both direct-contact *p* values are 9.999e-05. For 2FDB, the direct-contact values are 9.999e-05 and 0.005599. The 2FDB delta-SASA value changes from 0.045995 to 0.247475 and is therefore specification-sensitive. Counts are unchanged: 8/17 direct contacts in 1NUN, 10/17 in 2FDB, 13/17 in at least one complex and 15/17 within 8 Angstrom in at least one complex.

Thesis impact: retain direct-contact enrichment as the shared structural conclusion; do not present 2FDB delta-SASA as robust.

## 3. Deterministic 25-to-17 structural barcode

The new script `03b_select_structure_barcode.py` starts from the 25 FDR-significant columns and applies four outcome-independent criteria: different non-gap major residues, discriminating score at least 0.70, informative column and no gap-rich exclusion. It selects alignment columns 1, 2, 4, 6, 7, 8, 11, 16, 21, 28, 30, 31, 32, 46, 50, 51 and 52. Structural contact and delta-SASA values are not inspected during selection.

Thesis impact: state the formula and threshold in Methods and point to `barcode_selection_audit.tsv`; replace vague wording such as "additional conservation and mapping criteria".

## 4. IQ-TREE random seeds

The executable commands and manifest now record seeds 107252, 995853, 460166 and 160107 for the unconstrained, isoform-constrained, species-pair-constrained and AU-test runs. The three-row AU parse and numerical results are unchanged.

Thesis impact: add the four seeds to reproducibility information; no result claim changes.

## 5. Dependence-aware topology reporting

The caller-robustness output no longer attaches exact binomial intervals to 58 protein observations. It reports 58/58 proteins and 30/30 represented species descriptively for every leave-one-database-out scenario. Caller-offset uncertainty remains based on resampling species.

Thesis impact: remove any protein-level binomial confidence interval and report both descriptive levels. Explicitly acknowledge within-species and phylogenetic dependence.

## 6. Cross-annotation candidate selection

The documentation now matches the implementation. Within each species, source and predicted isoform, candidates are ranked lexicographically by reference coverage, score margin, best local-alignment score and protein length, all descending. The code does not rank by identity multiplied by coverage.

Thesis impact: replace any identity-times-coverage description with this exact lexicographic rule.

## 7. Status of analysis rules

The release describes rules as fixed and executable, not as preregistered or temporally pre-specified unless contemporaneous evidence exists. The structural-barcode rule is now explicit and machine-auditable.

Thesis impact: replace unsupported "pre-specified" claims with "fixed analysis rule" or "defined in the analysis workflow".

## 8. Versioning and citation metadata

The software version is synchronized at 0.2.1 and the case-study package at 1.0.1, both dated 22 August 2026. The author spelling is consistently Anastasiya Gryb. A version tag can provide a stable commit-level citation. No DOI is claimed; a DOI should be added only after an actual archival deposit such as Zenodo.

Thesis impact: cite the tagged release or archival DOI once available and record the exact commit. Do not invent or anticipate a DOI in the submitted text.
