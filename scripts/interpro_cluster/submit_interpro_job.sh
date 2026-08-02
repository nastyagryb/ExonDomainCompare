#!/usr/bin/env bash
set -euo pipefail
echo "This legacy FASTA-only wrapper cannot preserve the versioned run/profile contract." >&2
echo "Create a run, then use: python scripts/interpro_cluster/run_cluster_roundtrip.py --run-id <run_id>" >&2
exit 2
