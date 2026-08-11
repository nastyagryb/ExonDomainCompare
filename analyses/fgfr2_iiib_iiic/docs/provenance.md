# Provenance

This snapshot was assembled on 11 August 2026 from the local thesis working repository and added to the `main` branch of ExonDomainCompare. The source repository parent revision was `69738d44b58a32fe862818bed0b81e1344c1bae7`; the destination parent revision was `418a81495c249d7029393914f868fa3426c1bd5a`. The source working tree also contained unrelated development work, which is not included in this directory.

Framework material was copied from three existing run locations:

- `results/final_30_until_interpro_prepare/13_final_pre_interpro_closure`
- `results/final_30_until_interpro_prepare/15_exon_domain_boundary_post_interpro`
- `results/derived/example`

The downstream analysis scripts and frozen result tables came from `FGFR2_Analysis_Suite_Modern`. No file outside `analyses/fgfr2_iiib_iiic` is changed by this release.

The combined cassette alignment in the framework closure and the alignment used by the downstream suite have the same SHA-256 checksum:

```text
60faedde39e6fc462bd5c61e7370c2f3df61f73249c3a6f31d239d6c9dcf3cc2
```

The AU table in `results/05_iqtree/parsed_AU_topology_test.tsv` was regenerated from `topology_AU_test.iqtree`. The earlier four-row parse was not retained. The release check compares all three rows with the original report.

Only compact, analysis-relevant artifacts were selected. The package excludes transient work directories, caches, virtual environments, repeated raster/PDF encodings, API download caches and raw 24 GB run material. Every published file is listed in `checks/release_manifest.tsv` with its size and SHA-256 checksum.
