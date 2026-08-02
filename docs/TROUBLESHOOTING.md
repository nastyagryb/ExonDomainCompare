# Troubleshooting

## `.venv/bin/python` not found

Run commands from the repository root and complete setup first:

```bash
./scripts/setup_local.sh
```

## `start_frontend.sh` not found

Do not run the command from inside `frontend` or `backend`. From the repository
root, use:

```bash
./scripts/start_local.sh
```

## Python or Node.js version is rejected

Use Python 3.13 and Node.js 20.19+ or 22.12+. Check them with:

```bash
python3.13 --version
node --version
```

## `datasets` or `mafft` is missing

Install the missing program and open a new terminal. Then run:

```bash
.venv/bin/edc doctor --redact-paths
```

## Port 5173 or 8000 is already in use

Stop the older ExonDomainCompare terminal with `Ctrl+C`. Then start the
application again.

## LRZ login fails

Check the username and hostname. Test the same SSH target directly. Do not add a
password to the configuration file.

## The LRZ partition check fails

Use a partition assigned to the user's LRZ project. Re-run `cluster configure`
with the correct value.

## Cluster output was fetched but the page is not ready

Open `My Runs` and refresh the status. If Post-InterPro is still offered, use the
local Post-InterPro button once. Do not resubmit completed cluster jobs.

## A bundled dataset is missing

Run:

```bash
cd datasets
shasum -a 256 -c SHA256SUMS
```

Re-clone the repository if a checksum fails.
