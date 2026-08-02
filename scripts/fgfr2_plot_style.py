"""Compatibility import alias for :mod:`exondomaincompare.presentation.fgfr2_plot_style`."""
from __future__ import annotations

import sys
from pathlib import Path

_src = Path(__file__).resolve().parents[1] / "src"
if _src.is_dir() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

from exondomaincompare.presentation import fgfr2_plot_style as _canonical

sys.modules[__name__] = _canonical
