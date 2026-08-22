# Methods

## Framework-derived evidence

The FGFR2 run was organized around 30 vertebrate species and one selected IIIb and IIIc protein model per species. Transcript and CDS annotations were retained with the selected protein identifiers. Coding-exon boundaries were projected into amino-acid coordinates before they were compared with protein-domain intervals.

The framework produced a full-length primary-protein alignment, separate IIIb and IIIc cassette alignments, a combined cassette alignment, a boundary projection and alignment QC tables. Protein annotations were normalized from InterPro member-database calls and pyTMHMM output. Comparative tables joined those records with exon architecture, isoform labels and synteny evidence. The SVG files in `results/00_framework_outputs/website_figures` are exports of those comparative records.

The final model table contains 60 proteins from 30 species. The paired sequence analyses use 28 complete main-analysis species pairs (56 sequences). Two species retained review status and were not used for the paired sequence tests. The InterPro ensemble used for the boundary analyses contains 58 proteins across the 30-species set.

## Domain-caller analysis

Multiple signatures from the same member database and protein were collapsed to a within-protein median. This prevents a database with redundant signatures from receiving extra weight. Systematic offsets were estimated after accounting for species and isoform. Uncertainty was calculated by resampling species rather than individual domain calls. A leave-one-database-out analysis repeated the topology classification after removing each member database in turn.

## Robustness surface

Boundary classification was evaluated over all combinations of an absolute D3-end tolerance from 0 to 25 amino acids and a required member-database consensus from 50% to 100%. The surface records the fraction of proteins passing at each point and the smallest tolerance that recovers all proteins for each consensus requirement.

## Isoform-associated sequence signal

The 28 complete IIIb/IIIc species pairs were taken from the combined cassette MSA. Henikoff weights reduced the contribution of closely related sequences. Weighted Jensen–Shannon divergence was calculated per alignment column. Significance was assessed by 10,000 paired permutations that swapped IIIb and IIIc labels within species. Position-wise values were corrected with the Benjamini–Hochberg procedure. Structural analysis started from the FDR-significant set and retained positions with different non-gap major residues, an unweighted discriminating score of at least 0.70, an informative alignment column and no gap-rich exclusion. The score is the smaller within-isoform major-residue fraction multiplied by one minus the larger within-isoform gap fraction. This deterministic rule selects 17 of the 25 FDR-significant positions. Structural contact and solvent-accessibility outcomes are not used for selection.

## Phylogenetic topology test

Three maximum-likelihood hypotheses were compared: an unconstrained topology, an isoform-monophyly constraint and a constraint pairing IIIb and IIIc within each species. IQ-TREE used LG+G4 for the AU comparison with 10,000 RELL replicates. The released commands record seeds 107252, 995853, 460166 and 160107 for the unconstrained, isoform-constrained, species-pair-constrained and AU-test runs, respectively. The cassette alignment is short, so this test is interpreted as evidence about the competing topologies rather than a resolved species phylogeny.

## Leave-one-clade-out validation

One vertebrate clade was removed at a time. Discriminating sites were selected using only the remaining clades and the held-out proteins were then classified. Training-species bootstraps measured whether perfect held-out classification depended on a particular training sample. A paired label-permutation null preserved species pairing while removing isoform identity. Because both labels and the feature definition originate within the analysed cohort, this analysis tests cross-clade portability of the sequence-calibrated partition rather than independently validating the biological labels.

An earlier ancestral reconstruction is not part of the active analysis. Its target IIIb/IIIc ancestral split was not identifiable from the unrooted topology used. The implementation and outputs are preserved under `quarantined/invalid_ancestral_reconstruction` for audit, but no ancestral-state result is interpreted or release-verified.

## Structural interface mapping

The 17 conservative barcode positions were mapped to FGFR2b–FGF10 (PDB 1NUN) and FGFR2c–FGF8b (PDB 2FDB). Direct contact was defined at 4.5 Å and near-interface placement at 8 Å using minimum heavy-atom distance. Receptor-alone and complex solvent-accessible surface area were calculated with Biopython’s Shrake–Rupley implementation. Matched permutations controlled for receptor-alone accessibility and local sequence position. Within each structure, target and control residues were assigned with the same quartile boundaries. The primary analysis used control-derived accessibility quartiles; a sensitivity analysis used quartiles from the full cassette-residue universe.

## Cross-annotation replication

Matched NCBI and Ensembl models were evaluated for eight species and both isoforms. For each species, source and predicted isoform, candidates were ranked lexicographically by reference coverage, score margin, best local-alignment score and protein length, all in descending order. This is the implemented selection rule; it is not an identity-times-coverage score. Candidate selection and cassette localization were followed by the same InterPro-based topology rule. This is a separate FGFR2 validation layer; the website does not automatically compare NCBI and Ensembl annotations.

Topology success fractions are reported descriptively at both protein and species level. Protein-level exact binomial intervals are not used because two isoforms from one species are not independent observations and the species themselves share phylogenetic history. Species-clustered bootstrap uncertainty is used for caller-offset estimates.
