#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

if [[ ! -x "$ROOT/.venv/bin/python" ]]; then
  echo "Run ./scripts/setup_local.sh first." >&2
  exit 2
fi
PY="$ROOT/.venv/bin/python"

cd "$HERE/backend"

if ! "$PY" -c 'import fastapi, openpyxl' >/dev/null 2>&1; then
  echo "The project environment is incomplete; run ./scripts/setup_local.sh." >&2
  exit 2
fi

exec "$PY" -m uvicorn main:app --reload --port 8000
