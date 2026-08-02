"""Compatibility alias for :mod:`exondomaincompare.runs.registry`."""
from pathlib import Path
import sys

_src = Path(__file__).resolve().parents[2] / "src"
if _src.is_dir() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

from exondomaincompare.runs import registry as _implementation

sys.modules[__name__] = _implementation
