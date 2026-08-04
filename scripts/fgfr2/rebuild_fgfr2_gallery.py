#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from rebuild_dataset_indices import (
    fgfr2_run_ids_compat as fgfr2_run_ids, main, rebuild_freeze_dataset, rebuild_run,
)

if __name__ == "__main__":
    raise SystemExit(main())
