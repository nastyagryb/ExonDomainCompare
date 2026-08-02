# LRZ setup

An LRZ account is optional. It is needed only for InterProScan and pyTMHMM
cluster annotation of new runs.

## Information required

Ask the LRZ project administrator for:

- your LRZ username;
- the login hostname;
- an allowed Slurm partition;
- the Slurm account, if the project requires one.

## Configure the profile

Run this from the repository root:

```bash
.venv/bin/edc cluster configure \
  --user LRZ_USER \
  --host LRZ_LOGIN_HOST \
  --partition SLURM_PARTITION
```

Add `--account SLURM_ACCOUNT` only if required. The generated configuration is
private. Passwords and MFA tokens are never written to it.

## Check LRZ

```bash
.venv/bin/edc cluster doctor --redact-paths
```

Enter the normal LRZ password and MFA response in the terminal. The check is
read-only. It creates no directory and submits no job.

Continue only when the output contains:

```json
"ready_for_cluster_runs": true,
"missing": []
```

## Missing external tools

First display the installation plan:

```bash
.venv/bin/edc cluster tools install --tool all
```

Install only if the doctor reports InterProScan or pyTMHMM as missing:

```bash
.venv/bin/edc cluster tools install --tool all --confirm
```

The installation is user-scoped. It uses no `sudo`, changes no shell startup
file and submits no cluster job.

## Run annotation

Create the run in the website. Then copy its command from `My Runs`, or run:

```bash
.venv/bin/edc cluster roundtrip --run-id RUN_ID
```

Keep the terminal open until fetch and Post-InterPro analysis finish. Then use
`Refresh` on `My Runs`.

