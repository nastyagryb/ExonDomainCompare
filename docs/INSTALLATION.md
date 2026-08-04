# Installation

Windows users should install and run the application inside WSL2. Follow
[WINDOWS.md](WINDOWS.md) for the complete Windows procedure.

## 1. Install prerequisites

Required:

- Git
- Python 3.13
- Node.js 20.19+ or 22.12+
- npm
- MAFFT
- NCBI Datasets CLI

The official download pages are:

- [Python](https://www.python.org/downloads/)
- [Node.js](https://nodejs.org/en/download)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)

On macOS, Python, Node.js and MAFFT can also be installed with Homebrew:

```bash
brew install python@3.13 node@22 mafft
```

Install the NCBI Datasets CLI from its official page. Make sure the `datasets`
command is available in the terminal.

## 2. Clone the repository

```bash
git clone https://github.com/nastyagryb/ExonDomainCompare.git
cd ExonDomainCompare
```

## 3. Run setup

```bash
./scripts/setup_local.sh
```

Setup creates `.venv`, installs the tested Python and frontend dependencies,
creates private user-data folders and runs a local capability check.

## 4. Start the application

```bash
./scripts/start_local.sh
```

Open <http://127.0.0.1:5173>. Stop both local servers with `Ctrl+C`.

## 5. Check the installation

```bash
.venv/bin/edc doctor --redact-paths
```

The bundled datasets can be viewed without an LRZ account. Continue with
[LRZ.md](LRZ.md) only if new runs should use cluster annotation.

## Local run folder

Run outputs are kept outside the repository and are never committed. Use this
command to display the exact folder selected for the current installation:

```bash
.venv/bin/edc doctor
```

The defaults are `~/.local/share/ExonDomainCompare/runs` on Linux and WSL2,
and `~/Library/Application Support/ExonDomainCompare/runs` on macOS.
