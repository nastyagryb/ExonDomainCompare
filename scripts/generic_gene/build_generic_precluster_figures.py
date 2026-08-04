from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List

from exondomaincompare.generic_gene.common import GenericContext, load_context, read_json, read_tsv, write_tsv  # noqa: E402

# Import only the canonical plotting API (built on the FGFR2 primitives).
_SCRIPTS = Path(__file__).resolve().parents[1]
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
import plotting as plots  # noqa: E402

MANIFEST_COLUMNS = [
    "figure_id", "title", "scientific_question", "interpretation", "stage", "status",
    "svg", "pdf", "png", "source_files", "error",
]


def _int(v, default=0):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def _fig_transcript_structure(ctx: GenericContext, stem: str) -> bool:
    idx = read_json(ctx.run_dir / "website_indices" / "transcript_exon_structure_index.json", {}) or {}
    return plots.plot_transcript_exon_structure(
        ctx.figures_dir, stem, gene_symbol=ctx.gene_symbol,
        transcripts=idx.get("transcripts", []))


def _fig_exon_protein_architecture(ctx: GenericContext, stem: str) -> bool:
    arch = read_tsv(ctx.out("exon_protein_architecture.tsv"))
    sel = read_tsv(ctx.out("primary_selection_evidence.tsv"))
    primary_id = next((r.get("protein_id") for r in sel
                       if str(r.get("selected_primary", "")).lower() == "true"), "")
    if not primary_id and arch:
        primary_id = arch[0].get("protein_id", "")
    blocks = [r for r in arch if r.get("protein_id") == primary_id]
    clusters = read_tsv(ctx.out("event_region_candidate_clusters.tsv"))
    candidate_regions = [{"start_aa": c.get("representative_start_aa"),
                          "end_aa": c.get("representative_end_aa")} for c in clusters]
    return plots.plot_protein_exon_architecture(
        ctx.figures_dir, stem, gene_symbol=ctx.gene_symbol, primary_id=primary_id,
        exon_blocks=blocks, candidate_regions=candidate_regions)


def _fig_synteny(ctx: GenericContext, stem: str) -> bool:
    syn = read_tsv(ctx.out("synteny_neighbourhood.tsv"))
    return plots.plot_synteny_neighbourhood(ctx.figures_dir, stem,
                                            gene_symbol=ctx.gene_symbol, neighbours=syn)


def _fig_event_candidates(ctx: GenericContext, stem: str) -> bool:
    clusters = read_tsv(ctx.out("event_region_candidate_clusters.tsv"))
    arch = read_tsv(ctx.out("exon_protein_architecture.tsv"))
    max_aa = max((_int(r.get("protein_end_aa")) for r in arch), default=0) or None
    rows = [{"start_aa": c.get("representative_start_aa"),
             "end_aa": c.get("representative_end_aa"),
             "confidence": c.get("confidence"),
             "support_count": c.get("support_count")} for c in clusters]
    return plots.plot_evidence_regions_on_protein(
        ctx.figures_dir, stem, gene_symbol=ctx.gene_symbol, clusters=rows, max_aa=max_aa)


def _fig_isoform_alignment(ctx: GenericContext, stem: str) -> bool:
    idx = read_json(ctx.run_dir / "website_indices" / "isoform_alignment_index.json", {}) or {}
    candidates = read_json(
        ctx.run_dir / "website_indices" / "event_candidate_evidence_index.json", {}) or {}
    return plots.plot_isoform_alignment_overview(
        ctx.figures_dir, stem, gene_symbol=ctx.gene_symbol,
        sequences=idx.get("sequences", []), candidates=candidates.get("candidates", []))


def _fig_domain_architecture(ctx: GenericContext, stem: str) -> bool:
    domains = read_tsv(ctx.core("domain_features.tsv"))
    tm = read_tsv(ctx.core("tm_features.tsv"))
    exons = read_tsv(ctx.core("exon_protein_map.tsv"))
    selection = read_json(ctx.out("primary_selection_report.json"), {}) or {}
    candidates = read_json(
        ctx.run_dir / "website_indices" / "event_candidate_evidence_index.json", {}) or {}
    pid = selection.get("primary_protein_id", "")
    length = _int(selection.get("primary_length_aa"))
    return plots.plot_domain_architecture(
        ctx.figures_dir, stem, gene_symbol=ctx.gene_symbol, protein_id=pid,
        protein_length=length, domains=[d for d in domains if d.get("protein_id") == pid],
        tm_regions=[t for t in tm if t.get("protein_id") == pid],
        exon_boundaries=[_int(e.get("protein_end_aa")) for e in exons if e.get("protein_id") == pid],
        candidates=candidates.get("candidates", []))


def _fig_boundary_distribution(ctx: GenericContext, stem: str) -> bool:
    return plots.plot_exon_domain_boundary_distribution(
        ctx.figures_dir, stem, gene_symbol=ctx.gene_symbol,
        boundaries=read_tsv(ctx.core("exon_domain_boundary_distances.tsv")))


def _fig_candidate_domain_context(ctx: GenericContext, stem: str) -> bool:
    return plots.plot_candidate_domain_context(
        ctx.figures_dir, stem, gene_symbol=ctx.gene_symbol,
        contexts=read_tsv(ctx.out("candidate_domain_context.tsv")))


FIGURES = [
    ("transcript_exon_structure", "Transcript & Exon Structure", _fig_transcript_structure,
     "How do annotated transcripts differ at exon/CDS level?",
     "Transcript-relative exon structures reveal shared and alternative starts, ends and coding segments.",
     "website_indices/transcript_exon_structure_index.json"),
    ("primary_protein_exon_projection", "Primary Protein Exon Projection",
     _fig_exon_protein_architecture,
     "Where do coding-exon boundaries project onto the selected protein?",
     "Coding-exon projections and exploratory candidate regions are shown on the primary sequence.",
     "generic_gene_analysis/exon_protein_architecture.tsv;generic_gene_analysis/primary_selection_evidence.tsv"),
    ("isoform_alignment", "Isoform Alignment", _fig_isoform_alignment,
     "Where do protein models contain gaps or divergent sequence segments?",
     "The real within-species MAFFT alignment shows occupancy differences across protein models.",
     "results/07_msa/protein_alignment.faa;website_indices/isoform_alignment_index.json"),
    ("local_gene_neighbourhood", "Local Gene Neighbourhood", _fig_synteny,
     "Which annotated genes flank the target locus in this assembly?",
     "The target is centred between upstream and downstream annotated neighbours.",
     "generic_gene_analysis/synteny_neighbourhood.tsv"),
    ("exploratory_candidate_ranking", "Exploratory Candidate Evidence", _fig_event_candidates,
     "Which isoform-difference regions have the strongest transparent evidence?",
     "Candidates are ranked exploratory evidence and are not validated splicing events.",
     "generic_gene_analysis/event_candidate_ranking.tsv"),
]


def build(ctx: GenericContext) -> Dict[str, Any]:
    ctx.assert_not_freeze()
    ctx.figures_dir.mkdir(parents=True, exist_ok=True)
    plots.apply_style()
    manifest: List[Dict[str, Any]] = []
    made = 0
    for key, title, fn, question, interpretation, sources in FIGURES:
        stem = f"Figure_{key}"
        ok = False
        error = ""
        try:
            ok = fn(ctx, stem)
        except Exception as exc:
            ok = False
            error = str(exc)
        manifest.append({
            "figure_id": key,
            "title": title,
            "scientific_question": question,
            "interpretation": interpretation,
            "stage": "pre_cluster",
            "status": "available" if ok else "unavailable",
            "svg": f"figures/{stem}.svg" if ok else "",
            "pdf": f"figures/{stem}.pdf" if ok else "",
            "png": f"figures/{stem}.png" if ok else "",
            "source_files": sources,
            "error": error or ("" if ok else "Required source data were unavailable."),
        })
        made += int(ok)

    # Post-cluster figures remain pending unless real parsed cluster products exist.
    for key, title, fn, question, interpretation in [
        ("generic_domain_architecture", "Domain architecture", _fig_domain_architecture,
         "Which real InterPro/pyTMHMM features occur on the selected protein?",
         "Pending real InterProScan/pyTMHMM outputs."),
        ("generic_exon_domain_boundary_distribution", "Exon–domain boundary distribution",
         _fig_boundary_distribution,
         "How close are coding-exon boundaries to real domain edges?",
         "Pending real InterProScan domain coordinates."),
        ("generic_candidate_domain_context", "Candidate–Domain Context",
         _fig_candidate_domain_context,
         "How do exploratory candidates relate to real annotated domains?",
         "Pending real InterProScan domain coordinates."),
    ]:
        ok, error = False, ""
        if read_tsv(ctx.core("domain_features.tsv")):
            try:
                ok = fn(ctx, f"Figure_{key}")
            except Exception as exc:
                error = str(exc)
        manifest.append({
            "figure_id": key, "title": title, "scientific_question": question,
            "interpretation": (interpretation if not ok else
                               "Generated from real fetched InterProScan/pyTMHMM products."),
            "stage": "post_cluster", "status": "available" if ok else
                ("failed" if error else "pending_cluster"),
            "svg": f"figures/Figure_{key}.svg" if ok else "",
            "pdf": f"figures/Figure_{key}.pdf" if ok else "",
            "png": f"figures/Figure_{key}.png" if ok else "",
            "source_files": "domain_features.tsv;tm_features.tsv;exon_domain_boundary_distances.tsv",
            "error": error,
        })
        made += int(ok)

    write_tsv(ctx.out("figure_manifest.tsv"), manifest, MANIFEST_COLUMNS)
    return {"figure_manifest.tsv": len(manifest), "figures_generated": made}


def main() -> int:
    ap = argparse.ArgumentParser(description='Generic gene-agnostic pre-cluster figure files.')
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    ctx = load_context(args.run_id)
    res = build(ctx)
    print(f"OK figures  generated={res['figures_generated']}  manifest_rows={res['figure_manifest.tsv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
