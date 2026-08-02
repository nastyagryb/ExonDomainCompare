#!/usr/bin/env python3
"""Deprecated alias for the gene-agnostic rebuild command.

Rebuilding indices was never FGFR2-specific, so the command now lives at
``scripts/rebuild_dataset_indices.py`` and covers every gene. This wrapper keeps
the old invocation working:

    python scripts/fgfr2/rebuild_fgfr2_gallery.py --dataset example
    python scripts/rebuild_dataset_indices.py --dataset example   # the same thing
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from rebuild_dataset_indices import (  # noqa: E402,F401
    fgfr2_run_ids_compat as fgfr2_run_ids, main, rebuild_freeze_dataset, rebuild_run,
)

if __name__ == "__main__":
    raise SystemExit(main())
