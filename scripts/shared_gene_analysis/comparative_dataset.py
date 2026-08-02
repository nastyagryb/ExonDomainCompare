"""Compatibility entrypoint for :mod:`exondomaincompare.shared_gene_analysis.comparative_dataset`."""
from __future__ import annotations

import sys
from pathlib import Path

_src = Path(__file__).resolve().parents[2] / "src"
if _src.is_dir() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

if __name__ == "__main__":
    import runpy

    runpy.run_module("exondomaincompare.shared_gene_analysis.comparative_dataset", run_name="__main__")
else:
    from exondomaincompare.shared_gene_analysis import comparative_dataset as _canonical

    sys.modules[__name__] = _canonical
