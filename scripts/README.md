# Script inventory

Application and scientific library code lives in `src/exondomaincompare`.
This directory contains executable workflow and maintenance entry points only.

| Files | Classification | Purpose |
| --- | --- | --- |
| `edc.py`, `setup_local.sh`, `start_local.sh` | entry point | Local setup, diagnostics and application startup |
| `create_new_run.py`, `run_pre_interpro_for_run.py`, `run_post_interpro_for_run.py` | entry point | Run creation and FGFR2 workflow orchestration |
| Top-level `*.py` files not listed above | active workflow | Current FGFR2 pipeline stages called by the two root pipeline launchers |
| `fgfr2/*.py` | active workflow | FGFR2 projection, gallery and coordinate-model orchestration |
| `generic_gene/*.py` | active workflow | Generic exploratory gene-analysis stages |
| `interpro_cluster/*.py` | active workflow | Cluster submission, polling, retrieval and round-trip orchestration |
| `plotting/*.py`, `plotting/*.mjs` | active workflow | Figure registration, generation and shared rendering |
| `rebuild_dataset_indices.py` | maintenance | Rebuild published run indices without changing source evidence |
| `release/*.py` | maintenance | Public dataset assembly and release validation |

Compatibility aliases, historical one-command InterPro wrappers and unused
development scripts are not part of the public release.
