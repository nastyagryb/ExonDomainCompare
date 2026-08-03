#!/usr/bin/env python3
"""Regenerate the website indices of a dataset from its stored results.

    python scripts/rebuild_dataset_indices.py --dataset example
    python scripts/rebuild_dataset_indices.py --run-id <run_id>
    python scripts/rebuild_dataset_indices.py --gene FGFR2
    python scripts/rebuild_dataset_indices.py --all

One command for every dataset, whatever its gene. Rebuilding is how a change to
an index builder — a new figure scope, a canonical species order — reaches
datasets that were computed before it existed, without re-running any analysis:
only the derived JSON is written, never a result table or a FASTA.

The validated 30-species FGFR2 dataset is a read-only scientific record. Its
indices are read from the freeze and written to ``results/derived/example/``,
which the backend prefers when it exists, and the freeze is then verified
byte-for-byte. Ordinary runs are rebuilt in place.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import build_website_indices as bwi  # noqa: E402
from exondomaincompare.runs.registry import (  # noqa: E402
    RegistryError, discover_runs, resolve_run_record,
)
from exondomaincompare.config import load_config  # noqa: E402
from exondomaincompare.shared_gene_analysis.public_paths import sanitize_public_payload  # noqa: E402

FREEZE_CLOSURE = (ROOT / "results" / "final_30_until_interpro_prepare"
                  / "13_final_pre_interpro_closure")
DERIVED_ROOT = ROOT / "results" / "derived"
RUNTIME_CONFIG = load_config(repository_root=ROOT)

#: Indices this command regenerates for a validated (event-pipeline) dataset.
#: Everything else in the dataset stays exactly as computed.
REBUILT = ("synteny_locus_index.json", "figure_index.json",
           # The Figure 6B website view. Its IIIb/IIIc discriminating overlay is
           # derived per panel now, so an existing dataset has to be re-derived or its
           # interactive view keeps an overlay the figure beside it no longer draws.
           "cassette_residue_index.json")


def _digest(directory: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for path in sorted(directory.rglob("*")):
        if path.is_file():
            out[str(path.relative_to(directory))] = hashlib.sha256(
                path.read_bytes()).hexdigest()
    return out


def _rebuild(closure: Path, outdir: Path) -> Tuple[List[str], Dict[str, str]]:
    outdir.mkdir(parents=True, exist_ok=True)
    written: List[str] = []
    errors: Dict[str, str] = {}
    for name in REBUILT:
        builder = bwi.INDEX_BUILDERS.get(name)
        if builder is None:
            continue
        try:
            payload = sanitize_public_payload(builder(closure))
        except Exception as exc:  # a partial dataset simply keeps its old index
            errors[name] = f"{type(exc).__name__}: {exc}"
            continue
        (outdir / name).write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        written.append(name)
    return written, errors


def _render_derived_figures(indices: Path, figures: Path) -> Dict[str, object]:
    """Render the FGFR2 dataset's shared-renderer figures into the derived tree.

    The same vector-safe renderer every other gene uses: self-contained SVG, vector
    PDF, a 300-dpi raster for the web, and the source table beside each figure. No
    export depends on page CSS and none is a screenshot.

    The comparative figures are rendered from the reduced index — one model per
    species — because a comparative row can only be one protein; the per-species
    figures are rendered from the full index, where both isoform models are real.
    """
    import subprocess

    out: Dict[str, object] = {}
    model = indices / "protein_coordinate_model.json"
    reduced = indices / "comparative_model_index.json"
    dataset = indices / "comparative_dataset.json"

    jobs = [("main", ["node", "scripts/plotting/render_main_figures.mjs",
                      str(model), str(figures / "main")])]
    if reduced.is_file() and dataset.is_file():
        jobs.append(("comparative",
                     ["node", "scripts/plotting/render_comparative_gallery_figures.mjs",
                      str(reduced), str(dataset), str(figures / "comparative")]))

    # A renderer that is missing or fails is reported, not raised. Rebuilding the
    # indices is this command's job and it can finish without redrawing; aborting
    # here would leave the dataset with no rebuilt index because Node was absent.
    warnings: Dict[str, str] = {}
    for name, cmd in jobs:
        try:
            proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        except OSError as exc:
            warnings[name] = f"{type(exc).__name__}: {exc}"
            continue
        if proc.returncode != 0:
            warnings[name] = (proc.stderr or "").strip()[-400:]

    # The raster is derived from the vector, so it is produced here rather than by a
    # separate pass that could fall out of step with the SVG beside it.
    try:
        sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
        from generate_shared_main_figures import _rasterise_png
        made = failed = 0
        for svg in sorted(figures.rglob("*.svg")):
            png = svg.with_suffix(".png")
            if png.is_file() and png.stat().st_mtime >= svg.stat().st_mtime:
                continue
            if _rasterise_png(svg, png):
                made += 1
            else:
                failed += 1
        out["png_rasterised"] = made
        if failed:
            warnings["png"] = f"{failed} figure(s) could not be rasterised"
    except Exception as exc:  # pragma: no cover - optional layer
        warnings["png"] = f"{type(exc).__name__}: {exc}"
    if warnings:
        out["render_warnings"] = warnings
    return out


def rebuild_freeze_dataset(derived_root: Path = DERIVED_ROOT) -> Dict[str, object]:
    """Rebuild the validated dataset's indices outside the freeze and verify it."""
    before = _digest(FREEZE_CLOSURE)
    derived_root = Path(derived_root).expanduser().resolve()
    outdir = derived_root / "example" / "website_indices"
    written, errors = _rebuild(FREEZE_CLOSURE, outdir)

    # Rebuild the FGFR2 catalogue last so its curated figure index is authoritative.
    derived: Dict[str, object] = {}
    try:
        from fgfr2 import comparative_bridge, coordinate_model, gallery_catalogue
        _, inventory = coordinate_model.write_index(outdir)
        comparative_bridge.build(derived_root / "example")
        figures = derived_root / "example" / "figures"
        # Re-render before cataloguing. A card is only registered when its figure
        # exists, so cataloguing against stale files would silently drop the figures
        # a change to the renderer or the boundary contract has just added.
        render = _render_derived_figures(outdir, figures)
        catalogue = gallery_catalogue.write_catalogue(
            outdir,
            main_dir=(figures / "main") if (figures / "main").is_dir() else None,
            comparative_dir=(figures / "comparative")
            if (figures / "comparative").is_dir() else None)
        doc = json.loads(catalogue.read_text(encoding="utf-8"))
        written += ["protein_coordinate_model.json", "comparative_dataset.json",
                    "figure_catalogue.json", "figure_index.json",
                    str(inventory.relative_to(outdir))]
        derived = {"n_models": doc["availability"]["n_models"], **doc["counts"],
                   **render}
    except Exception as exc:
        errors["fgfr2_catalogue"] = f"{type(exc).__name__}: {exc}"

    after = _digest(FREEZE_CLOSURE)
    changed = sorted(k for k in set(before) | set(after) if before.get(k) != after.get(k))
    try:
        public_output = str(outdir.relative_to(ROOT))
    except ValueError:
        public_output = "<DERIVED_ROOT>/example/website_indices"
    return {
        "dataset": "example",
        "gene_symbol": "FGFR2",
        "output": public_output,
        "written": written,
        "catalogue": derived,
        "errors": errors,
        "freeze_unchanged": not changed,
        "freeze_changed_files": changed,
    }


def run_gene(run_dir: Path) -> str:
    cfg = run_dir / "run_config.json"
    if not cfg.is_file():
        return ""
    try:
        return str(json.loads(cfg.read_text(encoding="utf-8")).get("gene_symbol", ""))
    except (OSError, ValueError):
        return ""


def rebuild_coordinate_model(run_dir: Path) -> Dict[str, object]:
    """Rebuild a run's protein coordinate model from its stored result tables.

    The coordinate model is otherwise only written while a run executes, so a
    change to its contract — a new required field, a stricter rule — would reach
    new runs and leave existing ones behind. Rebuilding it here is what makes such
    a change a migration rather than a fork. No analysis is re-run: the model is
    derived from tables that already exist.
    """
    model_path = run_dir / "website_indices" / "generic" / "protein_coordinate_model.json"
    if not model_path.is_file():
        return {"rebuilt": False, "reason": "run has no coordinate model"}
    try:
        from exondomaincompare.shared_gene_analysis.protein_coordinate_model import build_models_for_run
        from exondomaincompare.shared_gene_analysis.validate_protein_coordinate_model import validate_index
        index = build_models_for_run(run_dir)
        errors = validate_index(index,
                                core_dir=run_dir / "results" / "core_gene_analysis")
        if errors:
            return {"rebuilt": False, "reason": f"{len(errors)} validation issue(s)",
                    "errors": errors[:5]}
        model_path.write_text(json.dumps(index, indent=2), encoding="utf-8")
        return {"rebuilt": True, "n_models": index.get("n_models", 0)}
    except Exception as exc:  # pragma: no cover - optional layer
        return {"rebuilt": False, "reason": f"{type(exc).__name__}: {exc}"}


def rebuild_exploratory_candidates(run_dir: Path) -> Dict[str, object]:
    """Re-derive a run's exploratory isoform-difference layer from its own sequences.

    The scan itself is a contract too: it once dropped any difference block shorter
    than six residues, which hid the real four-residue block between the two curated
    PTPN11 isoforms of mouse and rat. Re-deriving here is what makes removing that
    threshold a migration for existing runs instead of a rule that only new runs get.
    No analysis is re-run — the isoform sequences and the exon map already exist.
    """
    core = run_dir / "results" / "core_gene_analysis"
    if not (core / "event_candidate_regions.tsv").is_file():
        return {"rebuilt": False, "reason": "run has no exploratory candidate layer"}
    try:
        from exondomaincompare.framework.scan_isoform_event_candidates import scan, COLUMNS
        from exondomaincompare.framework.build_event_region_evidence import build_evidence
        from exondomaincompare.framework.cluster_event_region_evidence import build_clusters
        cfg = json.loads((run_dir / "run_config.json").read_text(encoding="utf-8")) \
            if (run_dir / "run_config.json").is_file() else {}
        rows = scan(run_dir)
        with (core / "event_candidate_regions.tsv").open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        ev = build_evidence(core, cfg.get("analysis_id", ""), cfg.get("gene_symbol", ""))
        cl = build_clusters(core, gap=5)
        mirrored = _mirror_exploratory_layer(run_dir)
        return {"rebuilt": True, "n_candidates": len(rows),
                "n_evidence": ev["n_evidence"], "n_clusters": cl["n_clusters"],
                **mirrored}
    except Exception as exc:  # pragma: no cover - optional layer
        return {"rebuilt": False, "reason": f"{type(exc).__name__}: {exc}"}


def _mirror_exploratory_layer(run_dir: Path) -> Dict[str, object]:
    """Carry a re-derived core candidate layer into the layers the pages read.

    ``core_gene_analysis/`` is where the scan writes, but a generic run serves its
    pages from ``generic_gene_analysis/`` and from the numbered stage folders, which
    are copies. Re-deriving only the core table would leave a run whose raw evidence
    holds the four-residue PTPN11 block while its Candidate page still says the
    isoforms are identical. Nothing is re-analysed here: the copies and the
    single-species indices are recomposed from the tables that were just written.
    """
    generic = run_dir / "results" / "generic_gene_analysis"
    if not generic.is_dir():
        return {"generic_layer": "absent"}
    from generic_gene import build_event_evidence
    from exondomaincompare.generic_gene import build_single_species_explorer
    from exondomaincompare.generic_gene.common import load_context
    ctx = load_context(run_dir.name)
    evidence = build_event_evidence.build(ctx)
    explorer = build_single_species_explorer.build(ctx)
    core = run_dir / "results" / "core_gene_analysis"
    stage_copies = {
        "04_event_evidence/event_region_evidence.tsv": "event_region_evidence.tsv",
        "05_event_region_detection/event_region_candidate_clusters.tsv":
            "event_region_candidate_clusters.tsv",
    }
    for rel, name in stage_copies.items():
        source = core / name
        if source.is_file():
            target = run_dir / "results" / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read_bytes())
    return {"generic_layer": "rebuilt",
            "n_generic_clusters": evidence["n_clusters"],
            "n_explorer_candidates": explorer.get("n_candidates", 0)}


def rebuild_run(run_id: str) -> Dict[str, object]:
    try:
        record = resolve_run_record(RUNTIME_CONFIG, run_id)
    except RegistryError as exc:
        raise SystemExit(str(exc)) from exc
    if record is None:
        raise SystemExit(f"Run not found: {run_id}")
    if record.read_only:
        raise SystemExit("Run is registered read-only; copy it before rebuilding indices.")
    run_dir = record.path
    # The exploratory layer feeds the coordinate model, so it is re-derived first.
    candidates = rebuild_exploratory_candidates(run_dir)
    base = {"dataset": run_id, "gene_symbol": run_gene(run_dir).upper(),
            "output": f"run:{run_id}/website_indices",
            "freeze_unchanged": True, "freeze_changed_files": [],
            "exploratory_candidates": candidates,
            "coordinate_model": rebuild_coordinate_model(run_dir)}
    closure = run_dir / "results" / "13_final_pre_interpro_closure"
    if not (closure / "final_pre_interpro_truth_table.tsv").is_file():
        # A run from the shared exploratory pipeline uses the shared builders.
        from exondomaincompare.shared_gene_analysis.build_fgfr2_compatible_indices import (  # noqa: WPS433
            build_fgfr2_compatible_indices)
        try:
            result = build_fgfr2_compatible_indices(run_dir)
        except Exception as exc:
            return {**base, "written": [], "errors": {"shared": f"{type(exc).__name__}: {exc}"}}
        return {**base, "written": [f"{s}.json" for s in result["written"]], "errors": {}}
    written, errors = _rebuild(closure, run_dir / "website_indices")
    # Write the curated FGFR2 Gallery after the generic index rebuild.
    gallery = _run_gallery_summary(run_dir)
    if gallery:
        written += ["figure_catalogue.json", "figure_index.json",
                    "protein_coordinate_model.json"]
    return {**base, "written": written, "errors": errors, "fgfr2_gallery": gallery}


def _run_gallery_summary(run_dir: Path) -> Dict[str, object]:
    """Write the run's FGFR2 Gallery and report it, or {} when the run has none."""
    try:
        from fgfr2 import run_gallery
        if not run_gallery.is_fgfr2_closure_run(run_dir):
            return {}
        built = run_gallery.write(run_dir)
        return {k: v for k, v in built.items() if k not in ("registration", "written")}
    except Exception as exc:  # pragma: no cover - optional layer
        return {"error": f"{type(exc).__name__}: {exc}"}


def run_ids(gene: str = "") -> List[str]:
    """Run directories, optionally restricted to one gene."""
    wanted = gene.strip().upper()
    out: List[str] = []
    records, _collisions = discover_runs(RUNTIME_CONFIG)
    for record in records:
        run_dir = record.path
        if run_dir.name.startswith((".", "_")):
            continue
        symbol = run_gene(run_dir).upper()
        if not symbol:
            continue
        if not wanted or symbol == wanted:
            out.append(run_dir.name)
    return out


def fgfr2_run_ids_compat() -> List[str]:
    """Kept for the deprecated scripts/fgfr2/rebuild_fgfr2_gallery.py wrapper."""
    return run_ids("FGFR2")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--run-id", help="Run directory name under runs/")
    group.add_argument("--dataset", choices=["example"],
                       help="The validated 30-species freeze dataset")
    group.add_argument("--gene", help="Every run of one gene, e.g. FGFR2")
    group.add_argument("--all", action="store_true",
                       help="The freeze dataset plus every run in the registry")
    ap.add_argument(
        "--derived-root", type=Path, default=DERIVED_ROOT,
        help="Derived example output root (tests should use a disposable directory).")
    args = ap.parse_args(argv)

    results: List[Dict[str, object]] = []
    if args.dataset == "example":
        results.append(rebuild_freeze_dataset(args.derived_root))
    elif args.run_id:
        results.append(rebuild_run(args.run_id))
    elif args.gene:
        if args.gene.strip().upper() == "FGFR2":
            results.append(rebuild_freeze_dataset(args.derived_root))
        results.extend(rebuild_run(rid) for rid in run_ids(args.gene))
    else:
        results.append(rebuild_freeze_dataset(args.derived_root))
        results.extend(rebuild_run(rid) for rid in run_ids())

    print(json.dumps(results, indent=2))
    failed = [r for r in results if r["errors"] or not r["freeze_unchanged"]]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
