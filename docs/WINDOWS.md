# Windows setup

ExonDomainCompare runs on Windows through WSL2 with Ubuntu. All commands after
the first step are entered in the Ubuntu terminal, not PowerShell.

## 1. Install WSL2

Open PowerShell as Administrator:

```powershell
wsl --install -d Ubuntu
```

Restart Windows, open Ubuntu and create the requested Linux username and
password.

## 2. Install system tools

```bash
sudo apt update
sudo apt install -y git curl ca-certificates build-essential mafft
```

Install Python 3.13 with uv:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source "$HOME/.local/bin/env"
uv python install 3.13
```

Install Node.js 22 with nvm:

```bash
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.4/install.sh | bash
source "$HOME/.bashrc"
nvm install 22
nvm alias default 22
```

Install the NCBI Datasets CLI:

```bash
mkdir -p "$HOME/.local/bin"
curl -o "$HOME/.local/bin/datasets" https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets
curl -o "$HOME/.local/bin/dataformat" https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/dataformat
chmod +x "$HOME/.local/bin/datasets" "$HOME/.local/bin/dataformat"
export PATH="$HOME/.local/bin:$PATH"
```

Add the PATH line permanently:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
```

## 3. Clone and start ExonDomainCompare

```bash
git clone https://github.com/nastyagryb/ExonDomainCompare.git
cd ExonDomainCompare
./scripts/setup_local.sh
./scripts/start_local.sh
```

Open <http://127.0.0.1:5173> in the Windows browser. Stop the application with
`Ctrl+C` in the Ubuntu terminal.

## 4. Optional LRZ setup

The bundled datasets work without LRZ. For cluster annotation of new runs,
follow [LRZ.md](LRZ.md) from the same Ubuntu terminal.

## Run location

New runs are stored inside WSL2 at:

```text
~/.local/share/ExonDomainCompare/runs
```

Show the exact path with:

```bash
.venv/bin/edc doctor
```
