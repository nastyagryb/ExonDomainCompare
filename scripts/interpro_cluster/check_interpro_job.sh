#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
if [[ $# -ne 1 ]]; then
  echo "Usage: scripts/interpro_cluster/check_interpro_job.sh <run_id>" >&2
  exit 2
fi
echo "Legacy wrapper: delegating to the versioned configuration contract." >&2
exec "${PYTHON:-python}" scripts/interpro_cluster/check_cluster_analysis.py --run-id "$1"
