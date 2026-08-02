#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/frontend"
if [[ ! -d node_modules ]]; then
  echo "Run ./scripts/setup_local.sh first." >&2
  exit 2
fi
exec npm run dev -- --port 5173
