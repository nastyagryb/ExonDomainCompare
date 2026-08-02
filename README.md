# ExonDomainCompare

ExonDomainCompare links exon structure, protein isoforms, domain annotations and
exon-domain boundaries across species.

The repository contains two analysis modes:

- a validated, read-only FGFR2 IIIb/IIIc thesis dataset;
- a generic exploratory workflow for other genes and species.

Generic Candidates are protein-isoform difference evidence. They are not
validated splicing events.

## Quick start

Requirements: Git, Python 3.13, Node.js 20.19+ or 22.12+, npm, MAFFT and the
NCBI Datasets CLI.

```bash
git clone https://github.com/nastyagryb/ExonDomainCompare.git
cd ExonDomainCompare
./scripts/setup_local.sh
./scripts/start_local.sh
```

Open <http://127.0.0.1:5173>. The bundled FGFR2 and BCL2L1 datasets work without
an LRZ account.

## New analyses

1. Open `New Analysis` in the website.
2. Enter a gene and one or more species.
3. Start the local pre-cluster analysis.
4. If domain annotation is needed, run the command shown on `My Runs`.
5. Refresh the run after the command finishes.

LRZ setup is required only for the cluster round trip. See
[docs/LRZ.md](docs/LRZ.md).

## Documentation

- [Installation](docs/INSTALLATION.md)
- [LRZ setup](docs/LRZ.md)
- [Troubleshooting](docs/TROUBLESHOOTING.md)
- [Release contents](docs/RELEASE_CONTENTS.md)
- [Architecture](docs/architecture/README.md)
- [Bundled datasets](datasets/README.md)

## Verification

```bash
make doctor
make release-check
make test
make lint
make build
```

## Citation

Citation metadata is provided in [CITATION.cff](CITATION.cff).

## License

The software is available under the [MIT License](LICENSE). Bundled research
data are covered by the separate [data notice](datasets/DATA_NOTICE.md).

## Scientific scope

The generic workflow produces exploratory evidence and confidence scores. It
does not establish a biologically validated splice event. The bundled FGFR2
dataset is a fixed thesis result and is loaded read-only.
