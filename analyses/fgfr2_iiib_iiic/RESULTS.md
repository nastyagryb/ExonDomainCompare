# Results

## Framework output

The final framework snapshot contains one IIIb and one IIIc model for each of 30 vertebrate species. It links transcript and protein identifiers, coding-exon projections, cassette coordinates, full-length and cassette MSAs, normalized domain calls, transmembrane predictions, synteny evidence and QC status. The comparative website exports make these layers inspectable from the same selected model set.

The website layer is evidence organization and exploration. The tests below were used to decide which patterns were robust enough to support thesis claims.

## Boundary consistency depends on context, but the topology is stable

The domain-caller analysis contained 290 protein-by-member-database observations from 58 proteins, 30 species and five InterPro member databases. Member-database identity explained a partial R² of 0.622 after accounting for species and isoform. Despite those systematic endpoint differences, every leave-one-database-out scenario retained the topology classification for all 58 proteins.

The threshold analysis covered 1,326 combinations. At 80% member-database consensus, a tolerance of 12 amino acids was sufficient for all 58 proteins. In total, 55.1% of the tested parameter grid yielded 58/58 passing proteins. The result is therefore not tied to one isolated threshold, although exact endpoint coordinates remain caller-dependent.

## The alternative cassette carries an isoform-associated sequence signal

The weighted sequence analysis used 28 complete species pairs, 56 sequences and 69 alignment columns. The global paired-label permutation test gave *p* = 9.999 × 10⁻⁵. Twenty-five columns passed the 5% FDR threshold. The later 17-site barcode is a deliberately narrower set selected for conservation and biological interpretation.

In the AU test, the unconstrained and isoform-monophyly trees had almost identical likelihoods and were not rejected. The species-pair topology was 913.81 log-likelihood units worse and was rejected (*p*AU = 6.36 × 10⁻⁸). These values are read directly from the original IQ-TREE report; the released parser is restricted to the three-row AU table.

When each vertebrate clade was held out in turn, all 56 proteins were assigned to the correct isoform. All 1,000 training-species bootstrap replicates per fold retained 100% held-out accuracy. The paired label-permutation null had a mean accuracy of 0.499 and a maximum of 0.821; the observed accuracy was 1.0 (*p* = 9.999 × 10⁻⁵).

## The signal is evolutionarily old and positioned at the ligand interface

Fifteen of the 17 conservative barcode positions were reconstructed as different at the IIIb/IIIc ancestral split. Eleven positions met the stricter criterion of posterior probability at least 0.95 and unchanged maximum-posterior state in every clade jackknife. This supports an ancient core, while leaving the taxon-sensitive positions explicitly uncertain.

All 17 barcode positions could be mapped in both structural complexes. Thirteen contacted a ligand directly in at least one complex and 15 lay within 8 Å in at least one complex. The matched direct-contact permutation was significant for both 1NUN and 2FDB. The delta-SASA result was significant for 1NUN but not for 2FDB, so contact enrichment is the stronger shared structural result.

## Annotation changes coordinates more readily than topology

The cross-annotation analysis compared 16 matched NCBI–Ensembl isoform pairs. All 16 retained the same topology class, although several absolute cassette and D3-end coordinates shifted. This replication supports a relative cassette–D3 relationship for the tested FGFR2 models. It does not establish that every annotation source or every gene will behave in the same way.

## Thesis-level conclusion

For FGFR2 IIIb/IIIc, exact domain endpoints vary with the caller and some absolute protein coordinates vary with the selected annotation model. The relation between the alternative cassette and its D3 context is nevertheless stable across the tested callers, parameter ranges and matched annotations. The cassette also carries a clade-generalizing, largely ancestral sequence signature concentrated at the ligand-facing surface.

The reusable contribution of ExonDomainCompare is the traceable path from model selection through exon-to-protein projection, alignment, protein annotation, comparative views and exportable evidence records. The FGFR2 analyses test one configured biological event built from those records; they do not turn the application into an automatic evolutionary or cross-provider inference system.
