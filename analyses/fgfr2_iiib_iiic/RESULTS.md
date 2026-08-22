# Results

## Framework output

The final framework snapshot contains one IIIb and one IIIc model for each of 30 vertebrate species. It links transcript and protein identifiers, coding-exon projections, cassette coordinates, full-length and cassette MSAs, normalized domain calls, transmembrane predictions, synteny evidence and QC status. The comparative website exports make these layers inspectable from the same selected model set.

The website layer is evidence organization and exploration. The tests below were used to decide which patterns were robust enough to support thesis claims.

## Boundary consistency depends on context, but the topology is stable

The domain-caller analysis contained 290 protein-by-member-database observations from 58 proteins, 30 species and five InterPro member databases. Member-database identity explained a partial R² of 0.622 after accounting for species and isoform. Despite those systematic endpoint differences, every leave-one-database-out scenario retained the topology classification for all 58 proteins and for all observed models in each of the 30 represented species. These 58/58 and 30/30 fractions are descriptive; no protein-level binomial confidence interval is attached to dependent isoform observations.

The threshold analysis covered 1,326 combinations. At 80% member-database consensus, a tolerance of 12 amino acids was sufficient for all 58 proteins. In total, 55.1% of the tested parameter grid yielded 58/58 passing proteins. The result is therefore not tied to one isolated threshold, although exact endpoint coordinates remain caller-dependent.

## The alternative cassette carries an isoform-associated sequence signal

The weighted sequence analysis used 28 complete species pairs, 56 sequences and 69 alignment columns. The global paired-label permutation test gave *p* = 9.999 × 10⁻⁵. Twenty-five columns passed the 5% FDR threshold. The deterministic conservation rule retained 17 positions for structural mapping; the released audit table records every inclusion criterion for all 25 positions and confirms that structural outcomes did not enter selection.

In the AU test, the unconstrained and isoform-monophyly trees had almost identical likelihoods and were not rejected. The species-pair topology was 913.81 log-likelihood units worse and was rejected (*p*AU = 6.36 × 10⁻⁸). These values are read directly from the original IQ-TREE report; the released parser is restricted to the three-row AU table.

When each vertebrate clade was held out in turn, all 56 proteins were assigned to their source isoform partition. All 1,000 training-species bootstrap replicates per fold retained 100% held-out accuracy. The paired label-permutation null had a mean accuracy of 0.499 and a maximum of 0.821; the observed accuracy was 1.0 (*p* = 9.999 × 10⁻⁵). This demonstrates cross-clade portability within the analysed cohort, not independent validation of labels that helped define the partition.

## The conserved sequence signal is positioned at the ligand interface

All 17 barcode positions could be mapped in both structural complexes. Thirteen contacted a ligand directly in at least one complex and 15 lay within 8 Å in at least one complex. Direct-contact enrichment was retained with both common-bin specifications: *p* = 9.999 × 10⁻⁵ under both schemes for 1NUN and *p* = 9.999 × 10⁻⁵ versus 0.00560 for 2FDB. The 1NUN delta-SASA result was retained under both specifications, whereas the 2FDB value changed from *p* = 0.0460 with control-derived quartiles to *p* = 0.247 with universe-derived quartiles. Direct contact is therefore the robust shared structural result; the 2FDB delta-SASA inference is specification-sensitive.

No ancestral-state conclusion is reported. The earlier reconstruction was quarantined because an unrooted tree did not identify the directional ancestral split required by the claim.

## Annotation changes coordinates more readily than topology

The cross-annotation analysis compared 16 matched NCBI–Ensembl isoform pairs. All 16 retained the same topology class, although several absolute cassette and D3-end coordinates shifted. This replication supports a relative cassette–D3 relationship for the tested FGFR2 models. It does not establish that every annotation source or every gene will behave in the same way.

## Conclusion

For FGFR2 IIIb/IIIc, exact domain endpoints vary with the caller and some absolute protein coordinates vary with the selected annotation model. The relation between the alternative cassette and its D3 context is nevertheless stable across the tested callers, parameter ranges and matched annotations. The cassette also carries a clade-portable isoform-associated sequence signature concentrated at the ligand-facing surface. The present analyses do not establish the ancestral state of that signature.

The reusable contribution of ExonDomainCompare is the traceable path from model selection through exon-to-protein projection, alignment, protein annotation, comparative views and exportable evidence records. The FGFR2 analyses test one configured biological event built from those records; they do not turn the application into an automatic evolutionary or cross-provider inference system.
