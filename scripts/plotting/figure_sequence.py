#!/usr/bin/env python3
"""The one Figure Gallery build sequence.

Every Gallery card in every run comes from this list, in this order, whether the
run has one species or five. The main-figure stage runs first because it owns
the Gallery's main cards and retires the ones it supersedes; the later stages
then register their supplements into an already-cleaned card set. The
comparative stage is last and skips cleanly on a single-species run.

Keeping the sequence here rather than inline in the pipeline is what makes a
species Scope of a multi-species run reproduce a standalone single-species
Gallery: both go through the same stages with the same arguments, so a species
Scope cannot drift into a separate, reduced implementation.

Usage for an existing run::

    python -m plotting.figure_sequence --run-id <run_id>
"""
from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from exondomaincompare.runs.registry import RegistryError, resolve_run_record  # noqa: E402
from exondomaincompare.config import load_config  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=ROOT)

# (human label, module). Each module exposes generate(run_dir, model_json).
FIGURE_STAGES: Tuple[Tuple[str, str], ...] = (
    ("shared main", "plotting.generate_shared_main_figures"),
    ("Exon Map", "plotting.generate_exon_map_figures"),
    ("Domain Architecture", "plotting.generate_domain_figures"),
    ("Exon–Domain Boundary", "plotting.generate_boundary_figures"),
    ("Isoform Alignment", "plotting.generate_alignment_figures"),
    ("Comparative Gallery", "plotting.generate_comparative_gallery_figures"),
    # Last, and the only owner of the availability record: earlier stages append
    # their cards, this one decides which of them a reader may actually see.
    ("figure registration", "plotting.figure_registration"),
)

MODEL_REL = Path("website_indices") / "generic" / "protein_coordinate_model.json"


def model_path(run_dir: Path) -> Path:
    return run_dir / MODEL_REL


def run_figure_stages(run_dir: Path, model_json: Path,
                      logline: Callable[[str], None] = print,
                      stages: Sequence[Tuple[str, str]] = FIGURE_STAGES,
                      ) -> List[Dict[str, Any]]:
    """Run the Gallery stages best-effort and return one result per stage.

    A stage that raises is reported and skipped: a missing optional input must
    not cost the run every other figure.
    """
    results: List[Dict[str, Any]] = []
    for label, module_name in stages:
        entry: Dict[str, Any] = {"stage": label, "module": module_name}
        try:
            gen = importlib.import_module(module_name).generate
            res = gen(run_dir, model_json) or {}
            entry.update({"ok": True, **res})
            logline(f"Wrote {res.get('figures', 0)} publication {label} figure file(s).")
        except Exception as err:  # pragma: no cover - plotting optional
            entry.update({"ok": False, "error": str(err)})
            logline(f"WARN: {label} figure generation skipped: {err}")
        results.append(entry)
    return results


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--rebuild-model", action="store_true",
                    help="rebuild the protein-coordinate model before the figures")
    args = ap.parse_args(argv)

    try:
        record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    except RegistryError as exc:
        ap.error(str(exc))
    if record is None:
        ap.error(f"no such run: {args.run_id}")
    if record.read_only:
        ap.error("run is registered read-only; copy it before rebuilding figures")
    run_dir = record.path
    out = model_path(run_dir)
    if args.rebuild_model:
        from exondomaincompare.shared_gene_analysis.protein_coordinate_model import build_models_for_run
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(build_models_for_run(run_dir), indent=2))
    if not out.exists():
        ap.error(f"no protein-coordinate model at {out}; pass --rebuild-model")

    results = run_figure_stages(run_dir, out)
    print(json.dumps({"run_id": args.run_id, "stages": results}, indent=2))
    return 0 if all(r.get("ok") for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
