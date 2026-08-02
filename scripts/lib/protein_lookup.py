"""Compatibility import alias for :mod:`exondomaincompare.scientific.protein_lookup`."""
from __future__ import annotations

import sys
from pathlib import Path

_src = Path(__file__).resolve().parents[2] / "src"
if _src.is_dir() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

from exondomaincompare.scientific import protein_lookup as _canonical

sys.modules[__name__] = _canonical
