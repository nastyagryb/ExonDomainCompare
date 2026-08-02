"""Compatibility launcher and import alias for :mod:`exondomaincompare.framework.run_core_gene_analysis`."""
from __future__ import annotations

import sys
from pathlib import Path

_src = Path(__file__).resolve().parents[2] / "src"
if _src.is_dir() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

if __name__ == "__main__":
    import runpy

    runpy.run_module("exondomaincompare.framework.run_core_gene_analysis", run_name="__main__")
else:
    from exondomaincompare.framework import run_core_gene_analysis as _canonical

    sys.modules[__name__] = _canonical
