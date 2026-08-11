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

The 28 complete IIIb/IIIc species pairs were taken from the combined cassette MSA. Henikoff weights reduced the contribution of closely related sequences. Weighted Jensen–Shannon divergence was calculated per alignment column. Significance was assessed by 10,000 paired permutations that swapped IIIb and IIIc labels within species. Position-wise values were corrected with the Benjamini–Hochberg procedure. The 17-site barcode used in the structural and ancestral analyses is a conservative subset with additional within-isoform conservation and mapping requirements; it is not the same object as the full set of FDR-significant columns.

## Phylogenetic topology test

Three maximum-likelihood hypotheses were compared: an unconstrained topology, an isoform-monophyly constraint and a constraint pairing IIIb and IIIc within each species. IQ-TREE used LG+G4 for the AU comparison with 10,000 RELL replicates. The cassette alignment is short, so this test is interpreted as evidence about the competing topologies rather than a resolved species phylogeny.

## Leave-one-clade-out validation

One vertebrate clade was removed at a time. Discriminating sites were selected using only the remaining clades and the held-out proteins were then classified. Training-species bootstraps measured whether perfect held-out classification depended on a particular training sample. A paired label-permutation null preserved species pairing while removing isoform identity.

## Ancestral reconstruction

Marginal empirical-Bayes ancestral states were reconstructed under LG+G4 on the existing maximum-likelihood tree. All-gap artifact columns were removed and gaps were treated as missing amino-acid observations. Clade jackknifes were used to identify states sensitive to taxon sampling. The ancestral sequences are conditional reconstructions from a short cassette alignment and should not be read as directly observed sequences.

## Structural interface mapping

The 17 conservative barcode positions were mapped to FGFR2b–FGF10 (PDB 1NUN) and FGFR2c–FGF8b (PDB 2FDB). Direct contact was defined at 4.5 Å and near-interface placement at 8 Å using minimum heavy-atom distance. Receptor-alone and complex solvent-accessible surface area were calculated with Biopython’s Shrake–Rupley implementation. Matched permutations controlled for receptor-alone accessibility and local sequence position.

## Cross-annotation replication

Matched NCBI and Ensembl models were evaluated for eight species and both isoforms. Candidate selection and cassette localization were followed by the same InterPro-based topology rule. This is a separate FGFR2 validation layer; the website does not automatically compare NCBI and Ensembl annotations.
