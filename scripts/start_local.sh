#!/usr/bin/env bash
set -euo pipefail

EDC_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EDC_REPOSITORY_ROOT="$(cd "$EDC_SCRIPT_DIR/.." && pwd)"
EDC_FRONTEND="$EDC_REPOSITORY_ROOT/webapp/frontend"

if [[ ! -x "$EDC_REPOSITORY_ROOT/.venv/bin/python" || ! -d "$EDC_FRONTEND/node_modules" ]]; then
  echo "Run ./scripts/setup_local.sh first." >&2
  exit 2
fi

cd "$EDC_REPOSITORY_ROOT"
"$EDC_REPOSITORY_ROOT/.venv/bin/python" -m uvicorn \
  webapp.backend.main:app --host 127.0.0.1 --port 8000 &
EDC_BACKEND_PID=$!

cleanup() {
  kill "$EDC_BACKEND_PID" 2>/dev/null || true
  wait "$EDC_BACKEND_PID" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

cd "$EDC_FRONTEND"
npm run dev -- --host 127.0.0.1 --port 5173
