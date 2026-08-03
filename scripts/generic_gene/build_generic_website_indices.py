"""Shared, gene-agnostic website indices.

Reads the canonical ``generic_gene_analysis/`` products and writes the shared
conceptual indices into ``website_indices/`` (same names for every gene):

  overview_index.json, evidence_stack.json, gene_explorer_index.json,
  protein_architecture_index.json, synteny_index.json, event_evidence_index.json,
  domain_architecture_index.json, exon_domain_boundaries_index.json,
  figures_index.json, available_views.json

Also writes the richer ``analysis_evidence_stack.json`` into the canonical layer.
No FGFR2 / IIIb / IIIc terminology.
"""
from __future__ import annotations

import argparse
from typing import Any, Dict, List

from exondomaincompare.framework import production_contract
from exondomaincompare.contracts import write_freshness_contract, write_payload_contracts
from exondomaincompare.shared_gene_analysis.public_paths import (
    sanitize_public_payload,
    write_public_download_projections,
)

from exondomaincompare.generic_gene.common import GenericContext, display_species, load_context, read_tsv, write_json
from exondomaincompare.generic_gene.stages import event_layer_for_gene

SCHEMA_VERSION = 2


def _gene_identity(ctx: GenericContext) -> Dict[str, Any]:
    """The requested-versus-source symbol record the core runner wrote, if any."""
    import json
    cfg = ctx.run_dir / "run_config.json"
    if not cfg.is_file():
        return {}
    try:
        data = json.loads(cfg.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    identity = data.get("gene_identity")
    if not isinstance(identity, dict) or not identity:
        return {}
    per_species = data.get("gene_identity_by_species")
    return {**identity,
            "by_species": per_species if isinstance(per_species, dict) else {}}


def _int(v, d=0):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return d


def _primary(sel: List[Dict[str, str]]) -> Dict[str, str]:
    for r in sel:
        if str(r.get("selected_primary", "")).lower() == "true":
            return r
    return sel[0] if sel else {}


def build(ctx: GenericContext) -> Dict[str, Any]:
    ctx.assert_not_freeze()
    routing = event_layer_for_gene(ctx.gene_symbol)

    gm = read_tsv(ctx.out("gene_model_summary.tsv"))
    iso = read_tsv(ctx.out("protein_isoform_summary.tsv"))
    sel = read_tsv(ctx.out("primary_selection_evidence.tsv"))
    arch = read_tsv(ctx.out("exon_protein_architecture.tsv"))
    syn = read_tsv(ctx.out("synteny_neighbourhood.tsv"))
    evid = read_tsv(ctx.out("event_region_evidence.tsv"))
    clusters = read_tsv(ctx.out("event_region_candidate_clusters.tsv"))
    msa = read_tsv(ctx.out("msa_index.tsv"))
    figman = read_tsv(ctx.out("figure_manifest.tsv"))

    species = sorted({r.get("species_id", "") for r in gm if r.get("species_id")})
    primary = _primary(sel)
    cluster_complete = ctx.cluster_status == "complete"
    domain_status = "available" if cluster_complete else "pending_cluster"
    boundary_status = "available" if cluster_complete else "pending_cluster"
    msa_status = (msa[0].get("msa_status") if msa else "unavailable_single_sequence")

    wi = ctx.run_dir / "website_indices"

    base_meta = {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": ctx.run_id,
        "analysis_id": ctx.analysis_id,
        "gene_symbol": ctx.gene_symbol,
        # How the requested symbol reached the assembly locus, so a gene annotated as
        # LOC… is explained rather than silently renamed.
        "gene_identity": _gene_identity(ctx),
        **routing,
        # Schema identity on every index (and therefore on every Gallery card built
        # from one), so a card can be checked against the run that produced it instead
        # of being trusted because its file happens to exist.
        **production_contract.resolve(ctx.gene_symbol).identity(),
    }

    # ---- protein_architecture_index ----
    by_protein: Dict[str, List[Dict[str, Any]]] = {}
    for b in arch:
        by_protein.setdefault(b.get("protein_id", ""), []).append({
            "exon_number": _int(b.get("exon_number")),
            "start_aa": _int(b.get("protein_start_aa")),
            "end_aa": _int(b.get("protein_end_aa")),
            "length_aa": _int(b.get("length_aa")),
            "phase": b.get("phase", ""),
        })
    protein_arch_index = {**base_meta, "primary_protein_id": primary.get("protein_id", ""),
                          "proteins": [{"protein_id": pid, "n_blocks": len(v), "blocks": v}
                                       for pid, v in by_protein.items()]}

    # ---- synteny_index ----
    synteny_index = {**base_meta,
                     "n_neighbours": len(syn),
                     "n_resolved": len([r for r in syn if r.get("status") == "resolved"]),
                     "neighbours": [{
                         "neighbor_symbol": r.get("neighbor_symbol", ""),
                         "side": r.get("side", ""), "order": _int(r.get("order")),
                         "orientation": r.get("orientation", ""),
                         "classification": r.get("classification", ""),
                         "status": r.get("status", ""),
                     } for r in syn]}

    # ---- event_evidence_index (per-source breakdown per cluster) ----
    def _cluster_sources(c: Dict[str, str]) -> List[Dict[str, str]]:
        out = []
        srcs = (c.get("evidence_sources", "") or "").split(";")
        for s in [x for x in srcs if x]:
            out.append({"source": s, "kind": (
                "exon_aligned" if "exon_aligned" in s else "protein_isoform_difference")})
        return out
    event_evidence_index = {**base_meta,
                            "event_status": routing["event_status"],
                            "n_candidate_clusters": len(clusters),
                            "n_evidence_rows": len(evid),
                            "disclaimer": "Exploratory isoform-difference candidates only; not validated events.",
                            "candidate_clusters": [{
                                "candidate_cluster_id": c.get("candidate_cluster_id", ""),
                                "start_aa": _int(c.get("representative_start_aa")),
                                "end_aa": _int(c.get("representative_end_aa")),
                                "length_aa": _int(c.get("representative_length_aa")),
                                "support_count": _int(c.get("support_count")),
                                "confidence": c.get("confidence", ""),
                                "confidence_reason": c.get("confidence_reason", ""),
                                "sources": _cluster_sources(c),
                            } for c in clusters]}

    # ---- figures_index ----
    figures_index = {**base_meta, "figures": [{
        "figure_id": f.get("figure_id", ""), "title": f.get("title", ""),
        "stage": f.get("stage", ""), "status": f.get("status", ""),
        "svg": f.get("svg", ""), "png": f.get("png", ""),
        "caption": f.get("caption", ""),
    } for f in figman]}

    # ---- domain + boundary indices (pending until cluster) ----
    domain_index = {**base_meta, "status": domain_status,
                    "reason": ("Domain architecture from InterProScan/pyTMHMM."
                               if cluster_complete else
                               "Pending the InterProScan/pyTMHMM cluster step.")}
    boundary_index = {**base_meta, "status": boundary_status,
                      "reason": ("All-exon distances to domain boundaries."
                                 if cluster_complete else
                                 "Pending the InterProScan/pyTMHMM cluster step.")}

    # Rich evidence stack.
    stack_items = _evidence_stack(ctx, routing, gm, iso, sel, arch, syn, clusters,
                                  msa_status, domain_status, boundary_status)
    evidence_stack = {**base_meta, "items": stack_items}

    # ---- gene_explorer_index ----
    gene_explorer_index = {**base_meta,
        "n_species": len(species), "n_gene_models": len(gm),
        "n_protein_isoforms": len(iso), "n_primary_proteins": len([r for r in iso
            if str(r.get("primary_status", "")).lower() == "primary"]),
        "species": [{"species_id": s, "display_name": display_species(s)} for s in species],
        "primary_protein_id": primary.get("protein_id", ""),
        "primary_selection_rule": primary.get("selection_rule", ""),
        "isoforms": [{"protein_id": r.get("protein_id", ""),
                      "transcript_id": r.get("transcript_id", ""),
                      "length_aa": _int(r.get("protein_length")),
                      "primary_status": r.get("primary_status", ""),
                      "source_kind": r.get("source_kind", "")} for r in iso]}

    # ---- available_views (shared concept; event-specific FGFR2 views off) ----
    views = {
        "overview": True,
        "gene_explorer": bool(gm),
        "protein_architecture": bool(arch),
        "synteny": bool(syn),
        "event_evidence": True,          # always present (exploratory or validated)
        "msa": msa_status == "available",
        "figures": any(f.get("status") == "available" for f in figman),
        "domain_architecture": domain_status == "available",
        "exon_domain_boundaries": boundary_status == "available",
        # FGFR2-only event-specific views are never enabled for generic genes:
        "cassette": False,
        "iiib_iiic_markers": False,
        "human_comparison": False,
        "boundary_consistency": False,
    }
    available_views = {**base_meta, "views": views,
                       "pending": {"domain_architecture": domain_status,
                                   "exon_domain_boundaries": boundary_status}}

    # ---- overview_index ----
    overview_index = {**base_meta,
        "kpis": {
            "species_analysed": len(species),
            "gene_models": len(gm),
            "protein_isoforms": len(iso),
            "primary_proteins": len([r for r in iso
                if str(r.get("primary_status", "")).lower() == "primary"]),
            "synteny_neighbours": synteny_index["n_resolved"],
            "exploratory_event_candidates": len(clusters),
        },
        "msa_status": msa_status,
        "domain_annotation": domain_status,
        "exon_domain_boundaries": boundary_status,
        "evidence_summary": [{"layer_id": it["layer_id"], "display_name": it["display_name"],
                              "status": it["status"], "confidence": it["confidence"]}
                             for it in stack_items],
        "available_views": views,
    }

    outputs = {
        "overview_index.json": overview_index,
        "evidence_stack.json": evidence_stack,
        "gene_explorer_index.json": gene_explorer_index,
        "protein_architecture_index.json": protein_arch_index,
        "synteny_index.json": synteny_index,
        "event_evidence_index.json": event_evidence_index,
        "domain_architecture_index.json": domain_index,
        "exon_domain_boundaries_index.json": boundary_index,
        "figures_index.json": figures_index,
        "available_views.json": available_views,
    }
    for name, data in outputs.items():
        write_json(wi / name, sanitize_public_payload(data))

    # The canonical rich stack also lives in the analysis layer.
    write_json(ctx.out("analysis_evidence_stack.json"), evidence_stack)
    write_public_download_projections(ctx.run_dir)
    write_freshness_contract(
        ctx.run_dir, wi,
        generator="scripts/generic_gene/build_generic_website_indices.py")
    write_payload_contracts(
        wi, run_id=ctx.run_id, dataset_id=ctx.run_id,
        generator="scripts/generic_gene/build_generic_website_indices.py")

    return {"website_indices": list(outputs.keys()),
            "analysis_evidence_stack.json": len(stack_items),
            "event_layer_type": routing["event_layer_type"]}


def _evidence_stack(ctx, routing, gm, iso, sel, arch, syn, clusters, msa_status,
                    domain_status, boundary_status) -> List[Dict[str, Any]]:
    primary = _primary(sel)
    n_primary = len([r for r in iso if str(r.get("primary_status", "")).lower() == "primary"])
    items = [
        {"layer_id": "gene_model", "display_name": "Gene model",
         "status": "available" if gm else "unavailable", "source": "NCBI GFF gene model",
         "confidence": "high" if gm else "n/a",
         "explanation": f"{len(gm)} gene/transcript model row(s) for {ctx.gene_symbol}.",
         "source_files": ["generic_gene_analysis/gene_model_summary.tsv"],
         "linked_view": "gene_explorer", "linked_figure": "generic_gene_model_overview"},
        {"layer_id": "protein_isoform", "display_name": "Protein isoform",
         "status": "available" if iso else "unavailable", "source": "Translated CDS isoforms",
         "confidence": "high" if iso else "n/a",
         "explanation": f"{len(iso)} protein isoform(s); {n_primary} selected as primary.",
         "source_files": ["generic_gene_analysis/protein_isoform_summary.tsv"],
         "linked_view": "gene_explorer", "linked_figure": "generic_gene_model_overview"},
        {"layer_id": "primary_selection", "display_name": "Primary selection",
         "status": "available" if sel else "unavailable", "source": "Documented selection hierarchy",
         "confidence": primary.get("confidence", "medium") or "medium",
         "explanation": f"Primary = {primary.get('protein_id','?')} via {primary.get('selection_rule','?')}.",
         "source_files": ["generic_gene_analysis/primary_selection_evidence.tsv"],
         "linked_view": "gene_explorer", "linked_figure": ""},
        {"layer_id": "exon_protein_architecture", "display_name": "Exon/protein architecture",
         "status": "available" if arch else "unavailable", "source": "GFF CDS → protein aa mapping",
         "confidence": "high" if arch else "n/a",
         "explanation": f"{len(arch)} coding-exon block(s) mapped to protein coordinates.",
         "source_files": ["generic_gene_analysis/exon_protein_architecture.tsv"],
         "linked_view": "protein_architecture", "linked_figure": "generic_exon_protein_architecture_primary"},
        {"layer_id": "synteny", "display_name": "Synteny",
         "status": "available" if syn else "unavailable", "source": "Genomic neighbourhood (GFF)",
         "confidence": "medium" if syn else "n/a",
         "explanation": f"{len([r for r in syn if r.get('status')=='resolved'])}/{len(syn)} neighbours resolved.",
         "source_files": ["generic_gene_analysis/synteny_neighbourhood.tsv"],
         "linked_view": "synteny", "linked_figure": "generic_synteny_neighbourhood"},
        {"layer_id": "exploratory_event", "display_name": "Exploratory event evidence",
         "status": "available" if clusters else "unavailable",
         "source": "Isoform-difference candidate search",
         "confidence": "exploratory",
         "explanation": (f"{len(clusters)} exploratory candidate cluster(s). "
                         "Not validated events; no markers or event labels invented."),
         "source_files": ["generic_gene_analysis/event_region_candidate_clusters.tsv"],
         "linked_view": "event_evidence", "linked_figure": "generic_exploratory_event_candidates"},
        {"layer_id": "external_uniprot", "display_name": "External UniProt evidence",
         "status": "optional", "source": "UniProt alternative-sequence (optional collector)",
         "confidence": "n/a",
         "explanation": "Optional external evidence; collected only when enabled.",
         "source_files": [], "linked_view": "event_evidence", "linked_figure": ""},
        {"layer_id": "msa", "display_name": "Multiple sequence alignment",
         "status": ("available" if msa_status == "available" else msa_status),
         "source": "MAFFT --auto (isoform-level)", "confidence": "n/a",
         "explanation": ("Isoform MSA built." if msa_status == "available"
                         else f"MSA {msa_status}."),
         "source_files": ["generic_gene_analysis/msa_index.tsv"],
         "linked_view": "msa", "linked_figure": ""},
        {"layer_id": "domain_annotation", "display_name": "Domain annotation",
         "status": domain_status, "source": "InterProScan / pyTMHMM",
         "confidence": "n/a",
         "explanation": ("Domain architecture available." if domain_status == "available"
                         else "Pending the InterProScan/pyTMHMM cluster step."),
         "source_files": ["generic_gene_analysis/domain_architecture.tsv"],
         "linked_view": "domain_architecture", "linked_figure": "generic_domain_architecture"},
        {"layer_id": "boundary_analysis", "display_name": "Exon–domain boundary analysis",
         "status": boundary_status, "source": "All-exon boundary distances",
         "confidence": "n/a",
         "explanation": ("Boundary analysis available." if boundary_status == "available"
                         else "Pending the InterProScan/pyTMHMM cluster step."),
         "source_files": ["generic_gene_analysis/exon_domain_boundary_analysis.tsv"],
         "linked_view": "exon_domain_boundaries", "linked_figure": "generic_exon_domain_boundary_distribution"},
    ]
    return items


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    ctx = load_context(args.run_id)
    res = build(ctx)
    print(f"OK website_indices  wrote={len(res['website_indices'])}  event_layer={res['event_layer_type']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
