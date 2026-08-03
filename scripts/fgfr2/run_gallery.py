#!/usr/bin/env python3
"""The modern FGFR2 Figure Gallery for an individual run.

The curated catalogue in ``gallery_catalogue`` was written for, and could only describe, the
validated 30-species dataset: its asset roots were module constants under
``results/final_30_until_interpro_prepare`` and its model inventory came from a post-cluster
table. A newly created FGFR2 run therefore fell back to ``build_website_indices``'s
directory listing — the legacy ``Framework`` / numbered ``Figure 1–4`` catalogue — even
though its closure carries exactly the same table and figure names.

This module supplies the two things the catalogue needed to become a function of a dataset
rather than a description of one:

* the run's own asset roots (``gallery_catalogue.run_dataset``), and
* the run's own model inventory (``coordinate_model.run_sources``), which before the cluster
  round-trip is built from the closure truth table with the cluster-derived layers stated as
  pending.

It is called from the per-run production path, so a run gets the modern Gallery when it is
created and again when its round-trip completes. No figure is rendered here and no
post-cluster value is invented: a card whose figure does not exist yet says so.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from fgfr2 import coordinate_model as cm  # noqa: E402
from fgfr2 import gallery_catalogue as gc  # noqa: E402
from exondomaincompare.runs.registry import RegistryError, resolve_run_record  # noqa: E402
from exondomaincompare.config import load_config  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=ROOT)

#: Where a run keeps figures rendered by the shared renderers, if it has any. Absent for a
#: run before its round-trip, which is why the catalogue tolerates a missing directory.
MAIN_FIGURES = ("figures", "main")
COMPARATIVE_FIGURES = ("figures", "comparative")

MAIN_RENDERER = "scripts/plotting/render_main_figures.mjs"
COMPARATIVE_RENDERER = "scripts/plotting/render_comparative_gallery_figures.mjs"


def is_fgfr2_closure_run(run_dir: Path) -> bool:
    """Whether this run was produced by the FGFR2 event pipeline.

    The closure's own truth table is the test, not the directory: the shared orchestrator
    materialises the numbered stage folders for every gene, so an empty
    ``13_final_pre_interpro_closure/`` says nothing about which pipeline ran.
    """
    closure = Path(run_dir) / "results" / "13_final_pre_interpro_closure"
    return (closure / "final_pre_interpro_truth_table.tsv").is_file()


def _figures_dir(run_dir: Path, parts: tuple) -> Optional[Path]:
    path = Path(run_dir).joinpath(*parts)
    return path if path.is_dir() else None


def catalogue_for(run_dir: Path, model_index: Dict[str, Any],
                  cluster_ready: bool) -> Dict[str, Any]:
    """The run's catalogue over whatever figures the run currently has on disk."""
    run_dir = Path(run_dir)
    return gc.build_catalogue(
        model_index,
        main_dir=_figures_dir(run_dir, MAIN_FIGURES),
        comparative_dir=_figures_dir(run_dir, COMPARATIVE_FIGURES),
        dataset=gc.run_dataset(run_dir, cluster_ready=cluster_ready))


def build(run_dir: Path) -> Dict[str, Any]:
    """The run's catalogue and its flattened Gallery index, without writing anything."""
    run_dir = Path(run_dir)
    sources = cm.run_sources(run_dir)
    model_index = cm.build_index(sources=sources)
    catalogue = catalogue_for(run_dir, model_index, sources.post_cluster_available)
    return {"model_index": model_index, "catalogue": catalogue,
            "figure_index": gc.flatten_for_gallery(catalogue)}


def render_figures(run_dir: Path, outdir: Path) -> Dict[str, str]:
    """Draw the run's own shared-renderer figures from its own coordinate model.

    Returns one message per renderer that did not run, keyed by renderer, so a missing
    Node or a renderer error is reported rather than raised: rebuilding the indices must
    still finish, and the catalogue then reports the undrawn figures as missing instead of
    claiming them.

    The run renders the same figures as the validated freeze, from the same renderers,
    because the alternative was a Gallery whose every post-cluster card was permanently
    missing for every run that is not the freeze.
    """
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
            # A one-species run has no comparative figures to draw. The renderer saying so
            # is the expected answer, not a failure of this rebuild.
            if "at least two species" not in message:
                problems[name] = message[-400:]

    # The renderers emit vector and table; the 300-dpi raster the web view shows is
    # derived from the SVG here, so it can never fall out of step with it.
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
    except Exception as exc:  # pragma: no cover - optional layer
        problems["png"] = f"{type(exc).__name__}: {exc}"
    return problems


def write(run_dir: Path, outdir: Optional[Path] = None) -> Dict[str, Any]:
    """Write the run's coordinate model, render its figures, then write the catalogue.

    The order is the point. The catalogue reports what exists on disk, so the figures have
    to be drawn from the freshly written coordinate model before it is taken.

    ``figure_index.json`` is written last and deliberately replaces the file the generic
    per-file builder produced. Two production catalogues for one gene is what let a run show
    the legacy one; there is now a single answer to what a reader sees.
    """
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
    # Only a run whose annotation is back has anything post-cluster to draw; before that
    # the figures do not exist because the science does not exist yet.
    problems = (render_figures(run_dir, outdir) if sources.post_cluster_available else {})

    catalogue = catalogue_for(run_dir, model_index, sources.post_cluster_available)
    written += [dump("figure_catalogue.json", catalogue),
                dump("figure_index.json", gc.flatten_for_gallery(catalogue))]
    # The same registration gate the generic Gallery goes through, so an FGFR2 run
    # cannot show a card whose outputs are missing or belong to another run.
    registration: Dict[str, Any] = {}
    try:
        from plotting.figure_registration import normalise_run
        registration = normalise_run(run_dir)
    except Exception as err:  # pragma: no cover - registration must not lose figures
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
    ap = argparse.ArgumentParser(description=__doc__,
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
