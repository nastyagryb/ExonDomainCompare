# FGFR2 IIIb/IIIc case study

This directory contains the material used for the FGFR2 case study in the bachelor’s thesis *An Annotation-Aware Framework for Assessing the Consistency and Robustness of Exon–Domain Boundary Identification Across Vertebrate Orthologs: FGFR2 IIIb/IIIc as a Case Study*.

It starts with the evidence produced by ExonDomainCompare and then follows the analyses used to test boundary consistency, sequence-level isoform identity, evolutionary depth and structural placement. The website exports are included together with their source tables. They are not screenshots and they are not treated as statistical tests.

## Analysis path

```text
gene and species set
        |
        v
model selection and exon-to-protein projection
        |
        v
MSA, synteny, InterPro and pyTMHMM evidence
        |
        v
comparative website views
        |
        v
frozen FGFR2 framework snapshot
        |
        +-- domain-caller and threshold sensitivity
        +-- sequence divergence and phylogenetic tests
        +-- leave-one-clade-out validation and ancestral reconstruction
        +-- structural interface mapping
        +-- cross-annotation replication
```

The framework layer records what was selected, mapped and annotated. The downstream layer asks whether the resulting observations survive changes in caller, threshold, clade and annotation source.

## What is included

| Location | Contents |
|---|---|
| `data/framework_snapshot/cohort` | Final 30-species model set, primary proteins and run status |
| `data/framework_snapshot/exon_protein_mapping` | Exon-to-protein projections and coordinate audits |
| `data/framework_snapshot/alignments` | Cassette and full-length MSAs, boundary projections and MSA QC |
| `data/framework_snapshot/protein_annotations` | Normalized InterPro, pyTMHMM and exon–domain records |
| `data/framework_snapshot/comparative` | Tables behind comparative website views |
| `results/00_framework_outputs` | Website figure exports and selected website indices |
| `results/01`–`08` | Frozen outputs of the downstream FGFR2 analyses |
| `scripts/analysis` | Analysis and plotting programs |
| `docs` | Methods, result summary, file map and data notes |

`results/00_framework_outputs/website_figures` contains SVG figures and TSV exports. PDF and PNG copies were omitted because they duplicate the same figures. The original 24 GB working run, caches, downloaded API responses and virtual environment are not part of this release.

## Quick check

The frozen release can be checked without rerunning the analyses:

```bash
cd analyses/fgfr2_iiib_iiic
python scripts/verify_release.py
```

The check covers the cohort size, key result values, AU-test parsing, required files, file sizes and the absence of local absolute paths in the published scripts.

## Reproducing the core analyses

Python 3.13.1 was used for the frozen run.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
make core
```

New outputs are written to `results/reproduced`; the frozen reference results are not overwritten. IQ-TREE and cross-provider acquisition are separate targets because they require an external executable or network access:

```bash
make iqtree IQTREE=iqtree3
make cross-annotation
```

The website views were rendered by the main application from the comparative data set. Their generator remains part of the application at `scripts/plotting/generate_comparative_gallery_figures.py`. This case-study directory preserves the exact exports, their input tables and the website indices used to identify them.

## Reading the results

- [METHODS.md](METHODS.md) describes the cohort and tests.
- [RESULTS.md](RESULTS.md) gives the numerical results and their limits.
- [docs/website_view_to_file_map.tsv](docs/website_view_to_file_map.tsv) links each website view to its scientific question and exported table.
- [docs/data_dictionary.md](docs/data_dictionary.md) explains the role of each data group.
- [docs/provenance.md](docs/provenance.md) records where the frozen files came from.
- [checks/release_manifest.tsv](checks/release_manifest.tsv) records SHA-256 checksums for the published files.

The release distinguishes exploratory views from confirmatory tests. A website plot can reveal a pattern or a case worth inspecting; the statistical claim belongs to the corresponding downstream analysis.
