"""Compatibility entrypoint for :mod:`exondomaincompare.framework.gene_config`."""
from __future__ import annotations

import sys
from pathlib import Path

_src = Path(__file__).resolve().parents[2] / "src"
if _src.is_dir() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

if __name__ == "__main__":
    from exondomaincompare.framework.gene_config import _main

    raise SystemExit(_main())

else:
    from exondomaincompare.framework import gene_config as _canonical

    sys.modules[__name__] = _canonical
