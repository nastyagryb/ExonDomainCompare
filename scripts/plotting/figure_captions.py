"""Compatibility import alias for :mod:`exondomaincompare.presentation.figure_captions`."""
from __future__ import annotations

import sys
from pathlib import Path

_src = Path(__file__).resolve().parents[2] / "src"
if _src.is_dir() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

from exondomaincompare.presentation import figure_captions as _canonical

sys.modules[__name__] = _canonical
