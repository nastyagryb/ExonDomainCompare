#!/usr/bin/env python3
"""Compatibility launcher for the installed :command:`edc` CLI."""
from pathlib import Path
import sys

_src = Path(__file__).resolve().parents[1] / "src"
if _src.is_dir() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

from exondomaincompare.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
