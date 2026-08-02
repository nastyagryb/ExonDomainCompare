#!/usr/bin/env bash
set -euo pipefail

EDC_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EDC_REPOSITORY_ROOT="$(cd "$EDC_SCRIPT_DIR/.." && pwd)"
cd "$EDC_REPOSITORY_ROOT"

if [[ -n "${EDC_PYTHON_BIN:-}" ]]; then
  EDC_SETUP_PYTHON="$EDC_PYTHON_BIN"
else
  EDC_SETUP_PYTHON="$(command -v python3.13 || true)"
fi

if [[ -z "$EDC_SETUP_PYTHON" ]]; then
  echo "Python 3.13 is required. Install it, then run this setup again." >&2
  exit 2
fi

"$EDC_SETUP_PYTHON" -c 'import sys; raise SystemExit(sys.version_info[:2] != (3, 13))' || {
  echo "The selected interpreter must be Python 3.13." >&2
  exit 2
}

command -v node >/dev/null 2>&1 || {
  echo "Node.js 20.19+ or 22.12+ is required." >&2
  exit 2
}
command -v npm >/dev/null 2>&1 || {
  echo "npm is required." >&2
  exit 2
}
node -e 'const [a,b]=process.versions.node.split(".").map(Number); process.exit((a===20&&b>=19)||a>=22?0:1)' || {
  echo "Node.js 20.19+ or 22.12+ is required." >&2
  exit 2
}

if [[ ! -x .venv/bin/python ]]; then
  "$EDC_SETUP_PYTHON" -m venv .venv
fi

.venv/bin/python -m pip install \
  -c requirements/constraints-py313-tested.txt \
  -e '.[test,render,synteny]'
npm ci --prefix webapp/frontend
.venv/bin/edc setup
.venv/bin/edc doctor --redact-paths

echo
echo "Setup complete. Start ExonDomainCompare with:"
echo "  ./scripts/start_local.sh"
