"""Guard maintained source code against hardcoded machine-local absolute paths.

Code that points at the throwaway sandbox mount used by an earlier authoring
session, or at somebody's home directory, only runs in one checkout on one
machine. This test scans the Python and JS/JSX source trees and fails if such a
path reappears. See ``FORBIDDEN_PATTERNS`` below for the exact shapes covered.

Scope, and why it is drawn this way
-----------------------------------
Only maintained *source code* is checked: ``.py``, ``.js``, ``.jsx``, ``.mjs``
(and ``.ts``/``.tsx`` should any appear) under ``tests/``, ``scripts/`` and
``webapp/``. Deliberately not checked:

* ``artifacts/*.md`` reports, ``docs/`` and ``QC_migration_report_*.md`` — prose
  that quotes the path a historical run used;
* ``runs/``, ``results/`` (including the immutable freeze under
  ``results/final_30_until_interpro_prepare/``) — run logs and provenance
  metadata whose whole purpose is to record where a run happened.

Rewriting those would falsify recorded provenance, so they stay as they are.

The forbidden prefixes below are assembled from fragments on purpose: that way
this file contains no literal machine-local path and does not flag itself.
"""
from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

SOURCE_TREES = ("tests", "scripts", "webapp")
SOURCE_SUFFIXES = {".py", ".js", ".jsx", ".mjs", ".ts", ".tsx"}

# Caches, dependencies, build output and downloaded browsers are not our source.
EXCLUDED_DIR_NAMES = {
    ".pw-browsers",
    ".pytest_cache",
    ".venv",
    "__pycache__",
    "build",
    "dist",
    "node_modules",
    "venv",
}

_SEP = "/"
_USER = r"[A-Za-z0-9._-]+"
FORBIDDEN_PATTERNS = {
    "sandbox_mount": re.compile(_SEP + "mnt" + _SEP + "data"),
    "macos_home_directory": re.compile(_SEP + "Users" + _SEP + _USER + _SEP),
    "linux_home_directory": re.compile(_SEP + "home" + _SEP + _USER + _SEP),
}


def iter_source_files():
    for tree in SOURCE_TREES:
        root = REPO_ROOT / tree
        assert root.is_dir(), f"Expected source tree '{tree}' below {REPO_ROOT}"
        for path in sorted(root.rglob("*")):
            if path.suffix not in SOURCE_SUFFIXES or not path.is_file():
                continue
            if EXCLUDED_DIR_NAMES.intersection(path.parts):
                continue
            yield path


def find_machine_local_paths(text: str):
    hits = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        for label, pattern in FORBIDDEN_PATTERNS.items():
            if pattern.search(line):
                hits.append((lineno, label, line.strip()))
    return hits


def test_source_trees_contain_no_machine_local_paths():
    offenders = []
    for path in iter_source_files():
        for lineno, label, line in find_machine_local_paths(
            path.read_text(encoding="utf-8", errors="replace")
        ):
            rel = path.relative_to(REPO_ROOT)
            offenders.append(f"{rel}:{lineno}: {label}: {line[:160]}")
    assert not offenders, (
        "Machine-local absolute paths found in maintained source code. Resolve "
        "paths relative to the repository root instead, e.g. via "
        "Path(__file__).resolve().parents[1]:\n" + "\n".join(offenders)
    )


def test_scan_covers_the_expected_source_trees():
    """A filter bug must not turn the scan above into a vacuous pass."""
    by_tree = {tree: 0 for tree in SOURCE_TREES}
    for path in iter_source_files():
        by_tree[path.relative_to(REPO_ROOT).parts[0]] += 1
    assert all(count > 0 for count in by_tree.values()), by_tree
    assert sum(by_tree.values()) > 100, by_tree


def test_detector_recognises_each_forbidden_path_shape():
    """The detector must react to all three path shapes it claims to cover."""
    samples = {
        "sandbox_mount": _SEP.join(["", "mnt", "data", "script.py"]),
        "macos_home_directory": _SEP.join(["", "Users", "someone", "project"]),
        "linux_home_directory": _SEP.join(["", "home", "someone", "project"]),
    }
    for expected_label, sample in samples.items():
        hits = find_machine_local_paths(f'SOME_PATH = "{sample}"')
        assert [label for _, label, _ in hits] == [expected_label], (sample, hits)

    clean = 'assert "' + _SEP + "Users" + _SEP + '" not in generated_output'
    assert find_machine_local_paths(clean) == [], (
        "A bare /Users/ needle used by other tests must not count as a violation"
    )
