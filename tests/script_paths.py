"""Repository-relative resolution of pipeline scripts for the test suite.

Tests must never depend on the current working directory or on absolute paths
from a particular machine or session. Everything is anchored on this file's
location inside the repository.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"


def script_path(name: str) -> Path:
    """Return the path of ``scripts/<name>`` and fail loudly if it is missing."""
    candidate = SCRIPTS_DIR / name
    if not candidate.exists():
        available = sorted(p.name for p in SCRIPTS_DIR.glob("*.py"))
        raise FileNotFoundError(
            f"{name} does not exist in {SCRIPTS_DIR}. Available top-level scripts: "
            + ", ".join(available)
        )
    return candidate


def load_script_module(name: str, module_name: str) -> ModuleType:
    """Import ``scripts/<name>`` as a module under ``module_name``."""
    path = script_path(name)
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot build an import spec for {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module
