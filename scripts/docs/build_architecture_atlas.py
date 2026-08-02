#!/usr/bin/env python3
"""Reproducible build for the ExonDomainCompare Architecture Atlas.

Renders every D2 diagram source in ``docs/architecture/sources/`` to a vector
SVG (``docs/architecture/svg/``) and a vector PDF (``docs/architecture/``),
then compiles the Typst Architecture Guide to PDF.

Toolchain (all open-source, installable via Homebrew):
    * d2            - diagram compiler        (brew install d2)
    * rsvg-convert  - SVG -> PDF converter    (brew install librsvg)
    * typst         - guide typesetting       (brew install typst)

Usage:
    python scripts/docs/build_architecture_atlas.py            # build everything
    python scripts/docs/build_architecture_atlas.py --check    # verify tools only
    python scripts/docs/build_architecture_atlas.py --diagrams # diagrams only
    python scripts/docs/build_architecture_atlas.py --guide    # guide only

Equivalent Makefile target: ``make architecture``.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from framework.portable_config import load_config  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)
ARCH = PROJECT_ROOT / "docs" / "architecture"
SOURCES = ARCH / "sources"
SVG = ARCH / "svg"
GUIDE_SRC = ARCH / "guide" / "ExonDomainCompare_Architecture_Guide.typ"
GUIDE_PDF = ARCH / "ExonDomainCompare_Architecture_Guide.pdf"

D2_LAYOUT = "dagre"


def _which(tool: str) -> str | None:
    config_name = {"rsvg-convert": "rsvg_convert"}.get(tool, tool)
    return RUNTIME_CONFIG.executable(config_name)


def _run(cmd: list[str]) -> None:
    env = dict(os.environ)
    subprocess.run(cmd, check=True, env=env)


def check_tools(required: tuple[str, ...]) -> dict[str, str | None]:
    found: dict[str, str | None] = {}
    print("== toolchain ==")
    for tool in required:
        path = _which(tool)
        found[tool] = path
        print(f"  {tool:<14}: {path or 'NOT FOUND'}")
    missing = [t for t, p in found.items() if p is None]
    if missing:
        print(
            "\nMissing tools: " + ", ".join(missing) +
            "\nInstall with: brew install " +
            " ".join({"d2": "d2", "rsvg-convert": "librsvg", "typst": "typst"}.get(m, m) for m in missing)
        )
    return found


def diagram_sources() -> list[Path]:
    return sorted(p for p in SOURCES.glob("*.d2") if not p.name.startswith("_"))


def build_diagrams() -> int:
    d2 = _which("d2")
    rsvg = _which("rsvg-convert")
    if not d2 or not rsvg:
        print("Cannot build diagrams: missing d2 and/or rsvg-convert.")
        return 1
    SVG.mkdir(parents=True, exist_ok=True)
    sources = diagram_sources()
    if not sources:
        print(f"No diagram sources found in {SOURCES}")
        return 1
    for src in sources:
        stem = src.stem
        svg_out = SVG / f"{stem}.svg"
        pdf_out = ARCH / f"{stem}.pdf"
        print(f"-- {stem}")
        _run([d2, "--layout", D2_LAYOUT, str(src), str(svg_out)])
        _run([rsvg, "-f", "pdf", "-o", str(pdf_out), str(svg_out)])
    print(f"Built {len(sources)} diagrams (SVG + PDF).")
    return 0


def build_guide() -> int:
    typst = _which("typst")
    if not typst:
        print("Cannot build guide: typst not found.")
        return 1
    if not GUIDE_SRC.exists():
        print(f"Guide source not found: {GUIDE_SRC}")
        return 1
    print("-- Architecture Guide (Typst)")
    _run([typst, "compile", "--root", str(PROJECT_ROOT), str(GUIDE_SRC), str(GUIDE_PDF)])
    print(f"Built guide: {GUIDE_PDF.relative_to(PROJECT_ROOT)}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="only verify the toolchain")
    ap.add_argument("--diagrams", action="store_true", help="build diagrams only")
    ap.add_argument("--guide", action="store_true", help="build the guide only")
    args = ap.parse_args()

    if args.check:
        found = check_tools(("d2", "rsvg-convert", "typst"))
        return 0 if all(found.values()) else 1

    rc = 0
    do_all = not (args.diagrams or args.guide)
    if args.diagrams or do_all:
        check_tools(("d2", "rsvg-convert"))
        rc |= build_diagrams()
    if args.guide or do_all:
        check_tools(("typst",))
        rc |= build_guide()
    return rc


if __name__ == "__main__":
    sys.exit(main())
