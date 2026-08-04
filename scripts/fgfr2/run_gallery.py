#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from fgfr2 import coordinate_model as cm
from fgfr2 import gallery_catalogue as gc
from exondomaincompare.runs.registry import RegistryError, resolve_run_record
from exondomaincompare.config import load_config

RUNTIME_CONFIG = load_config(repository_root=ROOT)

MAIN_FIGURES = ("figures", "main")
COMPARATIVE_FIGURES = ("figures", "comparative")

MAIN_RENDERER = "scripts/plotting/render_main_figures.mjs"
COMPARATIVE_RENDERER = "scripts/plotting/render_comparative_gallery_figures.mjs"


def is_fgfr2_closure_run(run_dir: Path) -> bool:
    closure = Path(run_dir) / "results" / "13_final_pre_interpro_closure"
    return (closure / "final_pre_interpro_truth_table.tsv").is_file()


def _figures_dir(run_dir: Path, parts: tuple) -> Optional[Path]:
    path = Path(run_dir).joinpath(*parts)
    return path if path.is_dir() else None


def catalogue_for(run_dir: Path, model_index: Dict[str, Any],
                  cluster_ready: bool) -> Dict[str, Any]:
    run_dir = Path(run_dir)
    return gc.build_catalogue(
        model_index,
        main_dir=_figures_dir(run_dir, MAIN_FIGURES),
        comparative_dir=_figures_dir(run_dir, COMPARATIVE_FIGURES),
        dataset=gc.run_dataset(run_dir, cluster_ready=cluster_ready))


def build(run_dir: Path) -> Dict[str, Any]:
    run_dir = Path(run_dir)
    sources = cm.run_sources(run_dir)
    model_index = cm.build_index(sources=sources)
    catalogue = catalogue_for(run_dir, model_index, sources.post_cluster_available)
    return {"model_index": model_index, "catalogue": catalogue,
            "figure_index": gc.flatten_for_gallery(catalogue)}


def render_figures(run_dir: Path, outdir: Path) -> Dict[str, str]:
    import subprocess

    from fgfr2 import comparative_bridge as cb

    problems: Dict[str, str] = {}
    jobs = [("main", [MAIN_RENDERER, str(outdir / "protein_coordinate_model.json"),
                      str(Path(run_dir).joinpath(*MAIN_FIGURES))])]
    try:
        cb.build(run_dir)
        jobs.append(("comparative",
                     [COMPARATIVE_RENDERER, str(outdir / "comparative_model_index.json"),
                      str(outdir / "comparative_dataset.json"),
                      str(Path(run_dir).joinpath(*COMPARATIVE_FIGURES))]))
    except Exception as exc:
        problems["comparative_dataset"] = f"{type(exc).__name__}: {exc}"

    for name, cmd in jobs:
        try:
            node = RUNTIME_CONFIG.executable("node")
            if not node:
                problems.append("Node is unavailable in the selected local profile.")
                continue
            proc = subprocess.run([node, *cmd], cwd=ROOT, capture_output=True, text=True)
        except OSError as exc:
            problems[name] = f"{type(exc).__name__}: {exc}"
            continue
        if proc.returncode != 0:
            message = (proc.stderr or proc.stdout or "").strip()
            if "at least two species" not in message:
                problems[name] = message[-400:]

    try:
        if str(ROOT / "scripts" / "plotting") not in sys.path:
            sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
        from generate_shared_main_figures import _rasterise_png

        failed = 0
        for svg in sorted(Path(run_dir).joinpath("figures").rglob("*.svg")):
            png = svg.with_suffix(".png")
            if png.is_file() and png.stat().st_mtime >= svg.stat().st_mtime:
                continue
            if not _rasterise_png(svg, png):
                failed += 1
        if failed:
            problems["png"] = f"{failed} figure(s) could not be rasterised"
    except Exception as exc:
        problems["png"] = f"{type(exc).__name__}: {exc}"
    return problems


def write(run_dir: Path, outdir: Optional[Path] = None) -> Dict[str, Any]:
    run_dir = Path(run_dir)
    outdir = Path(outdir) if outdir else (run_dir / "website_indices")
    outdir.mkdir(parents=True, exist_ok=True)

    sources = cm.run_sources(run_dir)
    model_index = cm.build_index(sources=sources)

    def dump(name: str, payload: Any) -> str:
        (outdir / name).write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return name

    written = [dump("protein_coordinate_model.json", model_index)]
    problems = (render_figures(run_dir, outdir) if sources.post_cluster_available else {})

    catalogue = catalogue_for(run_dir, model_index, sources.post_cluster_available)
    written += [dump("figure_catalogue.json", catalogue),
                dump("figure_index.json", gc.flatten_for_gallery(catalogue))]
    registration: Dict[str, Any] = {}
    try:
        from plotting.figure_registration import normalise_run
        registration = normalise_run(run_dir)
    except Exception as err:
        registration = {"error": str(err)}
    return {
        "registration": registration,
        "written": written,
        "dataset": catalogue["dataset"],
        "cluster_ready": catalogue["cluster_ready"],
        "default_scope": catalogue["default_scope"],
        "n_models": model_index["n_models"],
        "species_scopes": sorted(catalogue["species_scopes"]),
        "render_problems": problems,
        **catalogue["counts"],
    }


def main(argv: Optional[list] = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description='The modern FGFR2 Figure Gallery for an individual run.',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True, help="Run directory name under runs/")
    args = ap.parse_args(argv)
    try:
        record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    except RegistryError as exc:
        raise SystemExit(str(exc)) from exc
    if record is None:
        raise SystemExit(f"Run not found: {args.run_id}")
    if record.read_only:
        raise SystemExit("Run is registered read-only; copy it before rebuilding Gallery.")
    run_dir = record.path
    if not is_fgfr2_closure_run(run_dir):
        raise SystemExit(f"Not an FGFR2 closure run: {args.run_id}")
    print(json.dumps(write(run_dir), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
