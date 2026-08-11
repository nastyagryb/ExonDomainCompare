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

Windows users run ExonDomainCompare in WSL2. Follow the dedicated
[Windows setup](docs/WINDOWS.md) instead of the commands below.

```bash
git clone https://github.com/nastyagryb/ExonDomainCompare.git
cd ExonDomainCompare
./scripts/setup_local.sh
./scripts/start_local.sh
```

Open <http://127.0.0.1:5173>. The bundled FGFR2 and PTPN11 datasets work without
an LRZ account.

## One-time LRZ setup

This is required only before the first cluster annotation of a new run. It is
not required for the bundled datasets.

1. Configure the LRZ username, login host and Slurm partition.
2. Run the read-only cluster doctor.
3. If a required tool is missing, install it and run the doctor again.
4. Continue only when `ready_for_cluster_runs` is `true`.

The exact commands are in [docs/LRZ.md](docs/LRZ.md). Configuration is normally
done once per computer and LRZ account, before creating the first cluster run.

## New analyses

1. Open `New Analysis` in the website.
2. Enter a gene and one or more species.
3. Start the local pre-cluster analysis.
4. If domain annotation is needed, run the command shown on `My Runs`:
   `.venv/bin/edc cluster roundtrip --run-id RUN_ID`.
5. Refresh the run after the command finishes.

LRZ setup is required only for the cluster round trip. See
[docs/LRZ.md](docs/LRZ.md).

## Local run storage

Run outputs are stored outside the Git checkout so that a fresh clone stays
clean. Run `.venv/bin/edc doctor` to display the exact `runs_root` used by the
current installation.

- Linux and Windows/WSL2: `~/.local/share/ExonDomainCompare/runs`
- macOS: `~/Library/Application Support/ExonDomainCompare/runs`

These local results are not pushed to GitHub. Cluster results are downloaded
into the same local run folder after the round trip finishes.

## Documentation

- [Installation](docs/INSTALLATION.md)
- [Windows setup](docs/WINDOWS.md)
- [LRZ setup](docs/LRZ.md)
- [Troubleshooting](docs/TROUBLESHOOTING.md)
- [Release contents](docs/RELEASE_CONTENTS.md)
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
