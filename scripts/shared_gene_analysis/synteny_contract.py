"""Compatibility entrypoint for :mod:`exondomaincompare.shared_gene_analysis.synteny_contract`."""
from __future__ import annotations

import sys
from pathlib import Path

_src = Path(__file__).resolve().parents[2] / "src"
if _src.is_dir() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

from exondomaincompare.shared_gene_analysis import synteny_contract as _canonical

sys.modules[__name__] = _canonical
