"""Shared generic gene-analysis orchestrator.

Runs the shared, gene-agnostic modules for a run and materializes the SAME stage
folder structure used by every gene (``04_event_evidence`` … ``16_final_analyses``),
the canonical ``generic_gene_analysis/`` layer, and the shared ``website_indices/``.

The event stage runs a *generic exploratory event-region candidate search*; it
never claims validated events, invents markers, or uses IIIb/IIIc / cassette
terminology. Domain (15) and boundary (16) stages are materialized as
``pending_cluster`` until InterProScan/pyTMHMM outputs exist.

Usage:
    PYTHONPATH=scripts python -m generic_gene.run_generic_gene_analysis --run-id <id>
"""
from __future__ import annotations

import argparse
import datetime as _dt
import shutil
from pathlib import Path
from typing import Any, Dict, List

import sys as _sys

from . import (build_event_evidence, build_exon_protein_architecture,
               build_gene_model_summary, build_generic_msa_index,
               build_generic_precluster_figures, build_generic_website_indices,
               build_synteny_neighbourhood, select_primary_protein)
from exondomaincompare.generic_gene import build_single_species_explorer
from exondomaincompare.generic_gene.common import GenericContext, load_context, read_json, read_tsv, write_json, write_tsv
from exondomaincompare.generic_gene.stages import STAGES, event_layer_for_gene

# Real generated pre-cluster figure files, keyed by the figures_index item id used
# by the shared/rich figures index. These are enriched into the shared root
# figures_index.json so the Figure Gallery shows the actual files (PART 4).
_FIGURE_STEM_FOR_ID = {
    "transcript_exon_structure": "Figure_transcript_exon_structure",
    "primary_protein_exon_projection": "Figure_primary_protein_exon_projection",
    "isoform_alignment": "Figure_isoform_alignment",
    "local_gene_neighbourhood": "Figure_local_gene_neighbourhood",
    "exploratory_candidate_ranking": "Figure_exploratory_candidate_ranking",
}


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def _copy(src: Path, dst: Path) -> bool:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        return True
    return False


def _ensure_stage_dirs(ctx: GenericContext) -> None:
    results = ctx.run_dir / "results"
    for stage in STAGES:
        (results / stage).mkdir(parents=True, exist_ok=True)
    ctx.generic_dir.mkdir(parents=True, exist_ok=True)
    (ctx.run_dir / "website_indices").mkdir(parents=True, exist_ok=True)


def _materialize_stages(ctx: GenericContext, results_summary: Dict[str, Any]) -> Dict[str, Any]:
    """Copy canonical products into the shared numbered stage folders + reports."""
    R = ctx.run_dir / "results"
    stage_status: Dict[str, str] = {}

    # 04 event evidence
    s = R / "04_event_evidence"
    _copy(ctx.out("event_region_evidence.tsv"), s / "event_region_evidence.tsv")
    write_json(s / "event_region_evidence_summary.json", {
        "event_layer_type": event_layer_for_gene(ctx.gene_symbol)["event_layer_type"],
        "event_status": "exploratory_candidates_only",
        "n_evidence_rows": results_summary.get("n_evidence", 0),
        "disclaimer": "Exploratory isoform-difference evidence only; not validated events.",
        "generated_at": _now(),
    })
    stage_status["04_event_evidence"] = "available"

    # 05 event region detection
    s = R / "05_event_region_detection"
    _copy(ctx.out("event_region_candidate_clusters.tsv"), s / "event_region_candidate_clusters.tsv")
    write_json(s / "event_region_detection_report.json", {
        "detector": "generic_exploratory_event_candidate_search",
        "event_status": "exploratory_candidates_only",
        "n_candidate_clusters": results_summary.get("n_clusters", 0),
        "method": ("Protein-isoform pairwise difference regions, clustered; "
                   "checked against exon/protein boundaries. No markers/labels invented."),
        "generated_at": _now(),
    })
    stage_status["05_event_region_detection"] = "available"

    # 06 coordinate mapping
    s = R / "06_coordinate_mapping"
    _copy(ctx.out("exon_protein_architecture.tsv"), s / "exon_protein_architecture.tsv")
    write_json(s / "coordinate_mapping_report.json", {
        "n_exon_blocks": results_summary.get("n_exon_blocks", 0),
        "n_proteins": results_summary.get("n_proteins", 0),
        "source": "GFF CDS features mapped to protein aa coordinates.",
        "generated_at": _now(),
    })
    stage_status["06_coordinate_mapping"] = "available"

    # 07 msa
    s = R / "07_msa"
    msa_rows = read_tsv(ctx.out("msa_index.tsv"))
    msa_status = msa_rows[0].get("msa_status") if msa_rows else "unavailable_single_sequence"
    aln = ctx.generic_dir / "msa" / "isoform_msa.aln.faa"
    if aln.exists():
        _copy(aln, s / "protein_alignment.faa")
    write_json(s / "msa_status.json", {
        "msa_status": msa_status,
        "n_sequences": int(msa_rows[0].get("n_sequences", 0)) if msa_rows else 0,
        "tool": msa_rows[0].get("tool", "none") if msa_rows else "none",
        "reason": msa_rows[0].get("reason", "") if msa_rows else "",
        "alignment_file": "results/07_msa/protein_alignment.faa" if aln.exists() else "",
        "generated_at": _now(),
    })
    stage_status["07_msa"] = msa_status

    # 08 synteny
    s = R / "08_synteny"
    _copy(ctx.out("synteny_neighbourhood.tsv"), s / "synteny_neighbourhood.tsv")
    write_json(s / "synteny_report.json", {
        "n_neighbours": results_summary.get("n_neighbors", 0),
        "n_resolved": results_summary.get("n_resolved_neighbors", 0),
        "generated_at": _now(),
    })
    stage_status["08_synteny"] = "available"

    # 09 qc
    s = R / "09_qc"
    write_json(s / "run_qc_summary.json", {
        "gene_symbol": ctx.gene_symbol,
        "n_species": results_summary.get("n_species", 0),
        "n_gene_models": results_summary.get("n_gene_models", 0),
        "n_protein_isoforms": results_summary.get("n_protein_isoforms", 0),
        "n_synteny_neighbours": results_summary.get("n_neighbors", 0),
        "n_event_candidates": results_summary.get("n_clusters", 0),
        "msa_status": msa_status,
        "cluster_status": ctx.cluster_status,
        "generated_at": _now(),
    })
    stage_status["09_qc"] = "available"

    # 10 figures pre-domain (copy canonical figure files to stage-named files)
    s = R / "10_figures_pre_domain"
    fig_map = {
        "Figure_transcript_exon_structure": "transcript_exon_structure",
        "Figure_primary_protein_exon_projection": "primary_protein_exon_projection",
        "Figure_isoform_alignment": "isoform_alignment",
        "Figure_local_gene_neighbourhood": "local_gene_neighbourhood",
        "Figure_exploratory_candidate_ranking": "exploratory_candidate_ranking",
    }
    n_fig = 0
    for stem, out_stem in fig_map.items():
        for ext in ("svg", "pdf", "png"):
            if _copy(ctx.figures_dir / f"{stem}.{ext}", s / f"{out_stem}.{ext}"):
                n_fig += 1
    stage_status["10_figures_pre_domain"] = "available" if n_fig else "unavailable"

    # 15 domain architecture (pending until cluster)
    cluster_complete = ctx.cluster_status == "complete"
    s = R / "15_domain_architecture"
    dom_status = "available" if cluster_complete else "pending_cluster"
    _write_pending_domain(ctx, s, dom_status)
    stage_status["15_domain_architecture"] = dom_status

    # 16 final analyses (boundary; pending until cluster)
    s = R / "16_final_analyses"
    bnd_status = "available" if cluster_complete else "pending_cluster"
    _write_pending_boundary(ctx, s, bnd_status)
    stage_status["16_final_analyses"] = bnd_status

    return stage_status


def _write_pending_domain(ctx: GenericContext, stage_dir: Path, status: str) -> None:
    dom_cols = ["species_id", "protein_id", "domain_source", "domain_id", "domain_name",
                "start_aa", "end_aa", "score", "match_status", "interpro_accession",
                "interpro_description", "go_terms", "pathways", "feature_type",
                "domain_class_simplified", "source_file"]
    tm_cols = ["species_id", "protein_id", "start_aa", "end_aa", "source", "topology"]
    dom_rows: List[Dict[str, Any]] = []
    tm_rows: List[Dict[str, Any]] = []
    if status == "available":
        dom_rows = read_tsv(ctx.core("domain_features.tsv"))
        tm_rows = read_tsv(ctx.core("tm_features.tsv"))
    write_tsv(stage_dir / "domain_architecture.tsv", dom_rows, dom_cols)
    write_tsv(stage_dir / "tm_features.tsv", tm_rows, tm_cols)
    write_json(stage_dir / "domain_architecture_report.json", {
        "status": status,
        "reason": ("InterProScan/pyTMHMM domain architecture." if status == "available"
                   else "Pending the InterProScan/pyTMHMM cluster step."),
        "n_domains": len(dom_rows), "n_tm": len(tm_rows), "generated_at": _now(),
    })
    # canonical layer copies (PART 1 items 9)
    write_tsv(ctx.out("domain_architecture.tsv"), dom_rows, dom_cols)
    write_json(ctx.out("domain_architecture_status.json"), {"status": status})


def _write_pending_boundary(ctx: GenericContext, stage_dir: Path, status: str) -> None:
    cols = ["analysis_id", "gene_symbol", "species_id", "protein_id", "transcript_id",
            "exon_boundary_id", "boundary_position_aa", "nearest_domain_id",
            "nearest_domain_instance_id", "nearest_domain_start_aa", "nearest_domain_end_aa",
            "nearest_domain_name", "nearest_domain_boundary_type", "domain_edge_type",
            "signed_distance_aa", "absolute_distance_aa", "distance_aa", "category",
            "source"]
    rows: List[Dict[str, Any]] = []
    if status == "available":
        rows = read_tsv(ctx.core("exon_domain_boundary_distances.tsv"))
    write_tsv(stage_dir / "exon_domain_boundary_analysis.tsv", rows, cols)
    write_tsv(stage_dir / "exon_domain_boundary_distances.tsv", rows, cols)
    summary_rows = read_tsv(ctx.core("exon_domain_boundary_summary.tsv")) if status == "available" else []
    summary_cols = ["category", "count", "near_boundary_threshold_aa"]
    write_tsv(stage_dir / "exon_domain_boundary_summary.tsv", summary_rows, summary_cols)
    write_json(stage_dir / "boundary_analysis_report.json", {
        "status": status,
        "mode": "all_exon_boundaries",
        "reason": ("All-exon distances to domain boundaries." if status == "available"
                   else "Pending the InterProScan/pyTMHMM cluster step."),
        "n_rows": len(rows),
        "category_counts": {r.get("category", ""): r.get("count", "0")
                            for r in summary_rows if r.get("category")},
        "near_boundary_threshold_aa": 5,
        "source": "real fetched InterProScan coordinates" if status == "available" else "",
        "generated_at": _now(),
    })
    # canonical layer copies (PART 1 item 10)
    write_tsv(ctx.out("exon_domain_boundary_analysis.tsv"), rows, cols)
    write_json(ctx.out("exon_domain_boundary_analysis_status.json"), {"status": status})


def _update_status(ctx: GenericContext, stage_status: Dict[str, str], idx_result: Dict[str, Any]) -> None:
    routing = event_layer_for_gene(ctx.gene_symbol)
    status_path = ctx.run_dir / "status.json"
    st = read_json(status_path, {}) or {}
    st.update({
        "gene_symbol": ctx.gene_symbol,
        "analysis_id": ctx.analysis_id,
        **routing,
        "shared_pipeline": True,
        "stage_status": stage_status,
        "cluster_status": ctx.cluster_status,
        "website_indices": idx_result.get("website_indices", []),
        "generic_pipeline_generated_at": _now(),
    })
    if ctx.cluster_status == "pending" and st.get("cluster_analysis_status") == "not_started":
        st["current_step"] = "models_ready"
        st["cluster_roundtrip"] = {"phase": "pending_cluster", "updated_at": _now()}
    write_json(status_path, st)


def _write_shared_root_indices(ctx: GenericContext) -> Dict[str, Any]:
    """Write the RICH canonical indices to the shared ``website_indices/`` root.

    The frontend's shared renderer consumes ONE contract for every gene. Rather
    than a second flat schema, the shared root gets the same rich, UI-ready index
    shapes that FGFR2-style pages already understand (produced by the core index
    builder), plus real pre-cluster figure files enriched into figures_index.json.
    """
    scripts_dir = str(Path(__file__).resolve().parents[1])
    if scripts_dir not in _sys.path:
        _sys.path.insert(0, scripts_dir)
    from exondomaincompare.framework import build_core_gene_indices as bcgi  # noqa: E402

    cfg = bcgi._load_cfg(None, ctx.run_dir)
    src = bcgi.CoreSource(ctx.core_dir, f"run:{ctx.run_id}", cfg)
    root = ctx.run_dir / "website_indices"
    result = bcgi.build_core_indices(src, root)
    _enrich_root_figures(ctx, root)
    return {"website_indices_root": result.get("files", [])}


def _enrich_root_figures(ctx: GenericContext, root: Path) -> None:
    """Build a verified run-local figure contract; never return a URL to a missing file."""
    fig_path = root / "figures_index.json"
    manifest = read_tsv(ctx.out("figure_manifest.tsv"))
    figures = []
    for item in manifest:
        fid = item.get("figure_id", "")
        stem = _FIGURE_STEM_FOR_ID.get(fid, f"Figure_{fid}")
        requested = item.get("status", "")
        files = {ext: ctx.figures_dir / f"{stem}.{ext}" for ext in ("png", "svg", "pdf")}
        missing = [ext for ext, path in files.items() if not path.is_file() or path.stat().st_size == 0]
        status = requested
        error = item.get("error", "")
        if requested == "available" and missing:
            status = "failed"
            error = f"Generated figure file(s) missing or empty: {', '.join(missing)}"
        rel = f"results/generic_gene_analysis/figures/{stem}"
        base = f"/api/runs/{ctx.run_id}/files?path="
        def url(ext: str, inline: bool = False) -> str:
            if ext in missing or requested != "available":
                return ""
            from urllib.parse import quote
            return f"{base}{quote(rel + '.' + ext)}" + ("&inline=true" if inline else "")
        figures.append({
            "figure_id": fid, "title": item.get("title", ""),
            "scientific_question": item.get("scientific_question", ""),
            "interpretation": item.get("interpretation", ""),
            "stage": item.get("stage", ""), "status": status,
            "png_url": url("png", inline=True), "svg_url": url("svg"),
            "pdf_url": url("pdf"), "source_files": [
                x for x in item.get("source_files", "").split(";") if x],
            "error": error,
        })
    idx = {
        "schema_version": 2, "dataset_id": f"run:{ctx.run_id}",
        "analysis_id": ctx.analysis_id, "gene_symbol": ctx.gene_symbol,
        "figures": figures,
        "available": [{
            "id": f["figure_id"], "title": f["title"], "caption": f["interpretation"],
            "scientific_question": f["scientific_question"], "stage": f["stage"],
            "status": f["status"], "png_url": f["png_url"], "svg_url": f["svg_url"],
            "pdf_url": f["pdf_url"], "source_files": f["source_files"], "error": f["error"],
        } for f in figures if f["status"] == "available"],
        "pending": [{
            "id": f["figure_id"], "title": f["title"], "stage": f["stage"],
            "status": f["status"], "error": f["error"],
        } for f in figures if f["status"] != "available"],
    }
    _keep_publication_cards(fig_path, idx)
    write_json(fig_path, idx)


def _keep_publication_cards(fig_path: Path, idx: Dict[str, Any]) -> None:
    """Carry over the cards that the publication figure stages own.

    This function rebuilds the whole index from the pre-cluster figure manifest, so
    on its own it would delete the cards the publication figure stages registered
    into the same file. Those stages run after this orchestrator and retire the
    pre-cluster cards they replace, so the only thing needed here is not to throw
    their work away. Cards are recognised as belonging to another stage by their
    canonical Gallery category, which manifest-derived cards never carry.
    """
    previous = (read_json(fig_path, {}) or {}).get("figures") or []
    mine = {f.get("figure_id") for f in idx["figures"]}
    foreign = [f for f in previous
               if f.get("category") and f.get("figure_id") not in mine]
    if not foreign:
        return

    idx["figures"] = idx["figures"] + foreign
    idx["available"] = idx["available"] + [
        {"id": f.get("figure_id"), "title": f.get("title", ""),
         "caption": f.get("caption") or f.get("interpretation", ""),
         "scientific_question": f.get("scientific_question", ""),
         "stage": f.get("stage", ""), "status": f.get("status", "available"),
         "png_url": f.get("png_url", ""), "svg_url": f.get("svg_url", ""),
         "pdf_url": f.get("pdf_url", ""),
         "source_files": f.get("source_files") or [], "error": f.get("error", "")}
        for f in foreign if f.get("status", "available") == "available"]


def run(run_id: str) -> Dict[str, Any]:
    ctx = load_context(run_id)
    ctx.assert_not_freeze()
    _ensure_stage_dirs(ctx)

    summary: Dict[str, Any] = {}
    for mod in (build_gene_model_summary, select_primary_protein,
                build_exon_protein_architecture, build_synteny_neighbourhood,
                build_generic_msa_index, build_event_evidence):
        summary.update(mod.build(ctx))
    # Compose the scientifically richer single-species indices from the existing
    # model/parser/MSA/synteny products. This adds no alternative pipeline.
    summary.update(build_single_species_explorer.build(ctx))
    # figures depend on the canonical TSVs written above
    summary.update(build_generic_precluster_figures.build(ctx))
    # analysis-layer stack + canonical products (does not shadow the rich root)
    idx_result = build_generic_website_indices.build(ctx)
    # rich, UI-ready shared root indices (single contract for the frontend) —
    # written last so they are authoritative over any flat root files.
    try:
        root_result = _write_shared_root_indices(ctx)
        idx_result["website_indices_root"] = root_result.get("website_indices_root", [])
        # FGFR2-compatible coordinate_track / msa / synteny_locus for shared renderer parity.
        from exondomaincompare.shared_gene_analysis.build_fgfr2_compatible_indices import (  # noqa: E402
            build_fgfr2_compatible_indices,
        )
        compat = build_fgfr2_compatible_indices(ctx.run_dir)
        idx_result["fgfr2_compatible_indices"] = compat.get("written", [])
    except Exception as exc:  # pragma: no cover - non-blocking
        idx_result["shared_root_error"] = str(exc)
    # The legacy-compatible root builder above intentionally remains unchanged.
    # Re-apply the additive scientific single-species indices it does not know.
    summary.update(build_single_species_explorer.build(ctx))

    stage_status = _materialize_stages(ctx, summary)
    _update_status(ctx, stage_status, idx_result)

    summary["stage_status"] = stage_status
    summary["event_layer_type"] = idx_result.get("event_layer_type")
    summary["shared_root_error"] = idx_result.get("shared_root_error", "")
    return summary


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Shared generic gene-analysis orchestrator.")
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args(argv)
    res = run(args.run_id)
    print(f"OK generic pipeline  gene run={args.run_id}")
    print(f"   event_layer={res.get('event_layer_type')}  clusters={res.get('n_clusters')}  "
          f"figures={res.get('figures_generated')}  msa={res.get('msa_status')}")
    print(f"   stages: {res.get('stage_status')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
