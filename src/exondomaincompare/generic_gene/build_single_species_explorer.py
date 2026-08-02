"""Build the stage-aware scientific explorer products (species-aware).

This module composes existing shared pipeline outputs into the frontend indices
for the exploratory Gene Explorer. It is fully generic and species-aware: every
product (transcript structure, exploratory candidates, isoform alignment,
protein architecture) is built PER SPECIES from the per-species core tables
(``exon_protein_map.tsv``, ``protein_isoform_index.tsv``, the candidate cluster
TSV and the per-species MAFFT alignments). It never mixes isoforms from
different species into one alignment and never assumes a single species.

Each index carries a top-level ``species: [...]`` array (one entry per analysed
species) plus reference-species fields at the top level for backward
compatibility. It never infers domain/TM evidence before real cluster outputs
have been parsed and contains no gene- or species-specific hard-coding.
"""
from __future__ import annotations

import shutil
from datetime import datetime, timezone
from typing import Any, Dict, List, Mapping

from exondomaincompare.framework.primary_selection import classify_accession
from scripts.run_fgfr2_mafft_alignments import mafft_version

from exondomaincompare.shared_gene_analysis.common import shared_exon_group_id

from .common import GenericContext, read_fasta, read_json, read_tsv, write_json, write_tsv


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _int(value: Any, default: int | None = None) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


_shared_group = shared_exon_group_id


def _mafft_version() -> str:
    binary = shutil.which("mafft")
    return mafft_version(binary) if binary else "unavailable"


def _display_species(sid: str) -> str:
    return sid.replace("_", " ").strip().capitalize() if sid else ""


def _species_order(ctx: GenericContext, iso_rows: List[Mapping[str, Any]]) -> List[str]:
    """Analysed species in a stable order (species_list.txt is authoritative)."""
    slist = ctx.run_dir / "species_list.txt"
    order: List[str] = []
    if slist.is_file():
        order = [ln.strip() for ln in slist.read_text(encoding="utf-8").splitlines() if ln.strip()]
    for r in iso_rows:
        sid = str(r.get("species_id") or "")
        if sid and sid not in order:
            order.append(sid)
    return order or [""]


def _species_primaries(ctx: GenericContext) -> Dict[str, Dict[str, Any]]:
    """Map species_id -> {primary_protein_id, primary_transcript_id, primary_length_aa}.

    Prefers the multi-species ``primary_selection_index.json`` (species_primaries);
    falls back to the single primary report so a one-species run still works.
    """
    out: Dict[str, Dict[str, Any]] = {}
    idx = read_json(ctx.run_dir / "website_indices" / "primary_selection_index.json", {}) or {}
    for sp in idx.get("species_primaries") or []:
        sid = str(sp.get("species_id") or "")
        if sid:
            out[sid] = {
                "primary_protein_id": str(sp.get("primary_protein_id") or ""),
                "primary_transcript_id": str(sp.get("primary_transcript_id") or ""),
                "primary_length_aa": _int(sp.get("primary_length_aa")),
            }
    if not out:
        rep = read_json(ctx.out("primary_selection_report.json"), {}) or {}
        sid = str(idx.get("species_id") or rep.get("species_id") or "")
        out[sid] = {
            "primary_protein_id": str(rep.get("primary_protein_id")
                                      or idx.get("primary_protein_id") or ""),
            "primary_transcript_id": str(rep.get("primary_transcript_id")
                                         or idx.get("primary_transcript_id") or ""),
            "primary_length_aa": _int(rep.get("primary_length_aa")
                                      or idx.get("primary_length_aa")),
        }
    return out


# --------------------------------------------------------------------------- #
# per-species transcript / exon structure (from the per-species core tables)
# --------------------------------------------------------------------------- #
def _transcript_section(sid: str, primary_pid: str, primary_tid: str,
                        exon_rows: List[Mapping[str, Any]],
                        iso_rows: List[Mapping[str, Any]]) -> Dict[str, Any]:
    len_by_pid = {str(r.get("protein_id")): _int(r.get("protein_length"))
                  for r in iso_rows if r.get("protein_id")}
    tid_by_pid = {str(r.get("protein_id")): str(r.get("transcript_id") or "")
                  for r in iso_rows if r.get("protein_id")}
    # group coding exons by transcript
    by_tx: Dict[str, List[Mapping[str, Any]]] = {}
    for e in exon_rows:
        tid = str(e.get("transcript_id") or "")
        by_tx.setdefault(tid, []).append(e)
    transcripts: List[Dict[str, Any]] = []
    for tid, exons in by_tx.items():
        exons_sorted = sorted(exons, key=lambda x: _int(x.get("exon_number"), 0) or 0)
        pid = str(exons_sorted[0].get("protein_id") or "") if exons_sorted else ""
        tx_rows: List[Dict[str, Any]] = []
        for e in exons_sorted:
            n = _int(e.get("exon_number"), 0) or 0
            aa_s, aa_e = _int(e.get("protein_start_aa")), _int(e.get("protein_end_aa"))
            g_s, g_e = _int(e.get("cds_start")), _int(e.get("cds_end"))
            strand = str(e.get("normalized_strand") or e.get("strand") or "")
            length_aa = (aa_e - aa_s + 1) if (aa_s is not None and aa_e is not None) else None
            tx_rows.append({
                "transcript_id": tid, "protein_id": pid,
                "transcript_exon_number": n, "exon_label": f"E{n}",
                "exon_id": str(e.get("exon_id") or ""),
                "shared_exon_group_id": _shared_group(g_s, g_e, strand),
                "genomic_start": g_s, "genomic_end": g_e,
                "cds_start": g_s, "cds_end": g_e,
                "protein_start_aa": aa_s, "protein_end_aa": aa_e,
                "length_aa": length_aa,
                "phase": e.get("phase", ""), "strand": strand,
                "source": str(e.get("source") or "NCBI genomic GFF3"),
                "confidence": str(e.get("confidence") or "source_annotated"),
                "coding_status": "coding",
            })
        acc = classify_accession(pid, tid_by_pid.get(pid, tid))
        transcripts.append({
            "transcript_id": tid or tid_by_pid.get(pid, ""),
            "protein_id": pid,
            "is_primary": (pid == primary_pid) if primary_pid else (tid == primary_tid),
            "curation_status": ("curated" if acc["curated"] == "yes"
                                else "predicted" if acc["curated"] == "no" else "unknown"),
            "annotation_source": acc["source_label"],
            "protein_length": len_by_pid.get(pid),
            "exon_count": len(tx_rows),
            "exons": tx_rows,
        })
    # primary first, then by protein id
    transcripts.sort(key=lambda t: (not t["is_primary"], t.get("protein_id") or ""))
    return {
        "species_id": sid,
        "display_species_name": _display_species(sid),
        "primary_protein_id": primary_pid,
        "primary_transcript_id": primary_tid,
        "transcript_count": len(transcripts),
        "transcripts": transcripts,
    }


# --------------------------------------------------------------------------- #
# per-species exploratory candidates
# --------------------------------------------------------------------------- #
def _candidates_section(ctx: GenericContext, sid: str, primary_pid: str,
                        clusters: List[Mapping[str, Any]],
                        pair_rows: List[Mapping[str, Any]],
                        exon_rows: List[Mapping[str, Any]],
                        uniprot: Mapping[str, Any],
                        domain_context: Mapping[str, Any],
                        cluster_ready: bool) -> Dict[str, Any]:
    primary_exons = [e for e in exon_rows if str(e.get("protein_id")) == primary_pid]
    out: List[Dict[str, Any]] = []
    for i, c in enumerate(clusters, 1):
        cid = c.get("candidate_cluster_id") or f"CAND_{sid}_{i:03d}"
        start = _int(c.get("representative_start_aa"), 0) or 0
        end = _int(c.get("representative_end_aa"), start) or start
        supports = [r for r in pair_rows
                    if str(r.get("species_id") or sid) == sid and (
                        not r.get("event_candidate_id")
                        or r.get("event_candidate_id") in cid
                        or (_int(r.get("region_start_aa")) == start
                            and _int(r.get("region_end_aa")) == end))]
        start_match = any(_int(e.get("protein_start_aa")) == start for e in primary_exons)
        end_match = any(_int(e.get("protein_end_aa")) == end for e in primary_exons)
        support_count = _int(c.get("support_count"), len(supports)) or len(supports)
        iso_score = min(40, support_count * 5)
        exon_score = 20 if start_match and end_match else 10 if start_match or end_match else 0
        alignment_score = min(20, support_count * 3)
        external_status = "found" if uniprot.get("features") else (
            "not_found" if uniprot.get("status") in ("not_found", "complete") else "unavailable")
        external_score = 10 if external_status == "found" else 0
        domain = domain_context.get(cid, {})
        domain_score = _int(domain.get("domain_context_score"), 0) if cluster_ready else None
        penalty = 15 if end - start + 1 <= 1 else 0
        score = iso_score + exon_score + alignment_score + external_score + (domain_score or 0) - penalty
        confidence = "high" if score >= 65 else "medium" if score >= 35 else "low"
        out.append({
            "candidate_id": cid, "species_id": sid, "reference_protein": primary_pid,
            "aa_start": start, "aa_end": end, "length": end - start + 1,
            "candidate_class": "isoform_difference_region",
            "status": "exploratory_not_validated", "overall_score": score,
            "confidence_class": confidence,
            "confidence_reason": c.get("confidence_reason")
                or f"{support_count} supporting comparisons; exon-boundary score {exon_score}.",
            "protein_isoform_evidence": {
                "supporting_isoform_pairs": supports,
                "affected_proteins": sorted({x for r in supports for x in
                                             (r.get("protein_a", ""), r.get("protein_b", "")) if x}),
                "affected_transcripts": [x for x in c.get("transcripts_involved", "").split(";") if x],
                "classification": "present_absent_or_divergent",
                "sequence_difference_summary": f"{support_count} pairwise isoform comparison(s) support this region.",
                "evidence_confidence": confidence,
            },
            "exon_evidence": {
                "exon_aligned": start_match or end_match,
                "transcript_exon_numbers": [e.get("exon_number") for e in primary_exons
                                            if (_int(e.get("protein_start_aa"), 0) or 0) <= end
                                            and (_int(e.get("protein_end_aa"), 0) or 0) >= start],
                "exon_ids": [e.get("exon_id") for e in primary_exons
                             if (_int(e.get("protein_start_aa"), 0) or 0) <= end
                             and (_int(e.get("protein_end_aa"), 0) or 0) >= start],
                "projected_aa_coordinates": [start, end],
                "start_boundary_matches_exon_boundary": start_match,
                "end_boundary_matches_exon_boundary": end_match,
                "evidence_confidence": "high" if start_match and end_match else "medium",
                "explanation": "Transcript-relative CDS projection; shared genomic identity is recorded separately.",
            },
            "external_evidence": {"uniprot_status": external_status,
                                  "features": uniprot.get("features", []),
                                  "note": "Absence is not evidence against the candidate."},
            "alignment_evidence": {"alignment_start": start, "alignment_end": end,
                                   "supporting_aligned_isoforms": support_count,
                                   "gap_divergence_pattern": "isoform_difference_region",
                                   "alignment_confidence": confidence},
            "domain_evidence": (domain if cluster_ready else {
                "status": "pending_cluster", "overlap": None, "distance": None}),
            "score_components": {
                "isoform_support_score": iso_score, "exon_boundary_score": exon_score,
                "external_annotation_score": external_score, "alignment_score": alignment_score,
                "domain_context_score": domain_score, "penalty": penalty},
            "interpretation": (
                f"This region differs across multiple {ctx.gene_symbol} protein models in "
                f"{_display_species(sid)}"
                + (" and its boundaries coincide with projected coding-exon boundaries"
                   if start_match and end_match else "")
                + ". No curated UniProt alternative-sequence feature was found. "
                + ("Domain context uses fetched InterProScan output. " if cluster_ready
                   else "Domain context is pending InterProScan. ")
                + "It is not a validated splicing event."),
        })
    out.sort(key=lambda x: (-x["overall_score"], x["aa_start"]))
    return {
        "species_id": sid, "display_species_name": _display_species(sid),
        "reference_protein": primary_pid, "candidates": out,
    }


# --------------------------------------------------------------------------- #
# per-species within-species isoform alignment
# --------------------------------------------------------------------------- #
def _alignment_section(ctx: GenericContext, sid: str, primary_pid: str,
                       iso_rows: List[Mapping[str, Any]],
                       candidates: List[Mapping[str, Any]],
                       per_species_aln: Mapping[str, str]) -> Dict[str, Any]:
    tag = sid or "unknown"
    rel = per_species_aln.get(sid) or f"results/generic_gene_analysis/msa/isoform_msa__{tag}.aln.faa"
    path = ctx.run_dir / rel
    seqs = read_fasta(path) if path.is_file() else {}
    primary = primary_pid if primary_pid in seqs else next(iter(seqs), "")
    ordered_ids = ([primary] if primary in seqs else []) + [k for k in seqs if k != primary]
    alignment_length = len(next(iter(seqs.values()))) if seqs else 0
    valid = bool(seqs) and len({len(v) for v in seqs.values()}) == 1 and len(seqs) >= 2
    columns: List[str] = []
    if valid:
        for i in range(alignment_length):
            residues = [seqs[k][i] for k in ordered_ids if seqs[k][i] != "-"]
            columns.append("*" if residues and len(set(residues)) == 1 else ":" if residues else " ")
    tid_by_pid = {str(r.get("protein_id")): str(r.get("transcript_id") or "")
                  for r in iso_rows if r.get("protein_id")}
    return {
        "species_id": sid, "display_species_name": _display_species(sid),
        "status": "available" if valid else "unavailable",
        "reference_sequence": primary,
        "sequence_count": len(seqs), "alignment_length": alignment_length,
        "alignment_file": rel if path.is_file() else "",
        "sequences": [{"protein_id": pid, "transcript_id": tid_by_pid.get(pid, ""),
                       "is_primary": pid == primary, "aligned_sequence": seqs[pid]}
                      for pid in ordered_ids],
        "consensus": "".join(columns),
        "candidates": [{"candidate_id": c["candidate_id"], "aa_start": c["aa_start"],
                        "aa_end": c["aa_end"]} for c in candidates],
    }


# --------------------------------------------------------------------------- #
# per-species protein architecture (pre-cluster: exon boundaries only)
# --------------------------------------------------------------------------- #
def _architecture_section(sid: str, primary_pid: str, iso_rows: List[Mapping[str, Any]],
                          exon_rows: List[Mapping[str, Any]],
                          candidates: List[Mapping[str, Any]],
                          domains: List[Mapping[str, Any]], tm: List[Mapping[str, Any]],
                          cluster_ready: bool) -> Dict[str, Any]:
    proteins: List[Dict[str, Any]] = []
    for row in iso_rows:
        pid = str(row.get("protein_id") or "")
        exons = [e for e in exon_rows if str(e.get("protein_id")) == pid]
        proteins.append({
            "protein_id": pid, "transcript_id": row.get("transcript_id", ""),
            "length_aa": _int(row.get("protein_length"), 0),
            "is_primary": (pid == primary_pid) or row.get("primary_status") == "primary",
            "exon_boundaries": sorted({_int(e.get("protein_end_aa"), 0) for e in exons
                                       if (_int(e.get("protein_end_aa"), 0) or 0) > 0}),
            "candidate_overlays": [{
                "candidate_id": c["candidate_id"], "start_aa": c["aa_start"],
                "end_aa": c["aa_end"], "status": c["status"],
            } for c in candidates if c.get("reference_protein") == pid],
            "domains": [{
                "source": d.get("domain_source", ""), "signature_accession": d.get("domain_id", ""),
                "description": d.get("domain_name", ""), "start_aa": _int(d.get("start_aa")),
                "end_aa": _int(d.get("end_aa")), "score": d.get("score", ""),
                "interpro_accession": d.get("interpro_accession", ""),
                "interpro_description": d.get("interpro_description", ""),
                "feature_type": d.get("feature_type", "domain"),
            } for d in domains if d.get("protein_id") == pid],
            "tm_regions": [{
                "start_aa": _int(t.get("start_aa")), "end_aa": _int(t.get("end_aa")),
                "source": t.get("source", "pytmhmm"), "topology": t.get("topology", ""),
            } for t in tm if t.get("protein_id") == pid],
        })
    return {"species_id": sid, "display_species_name": _display_species(sid),
            "proteins": proteins}


# --------------------------------------------------------------------------- #
# reference-species synteny back-compat index (multi-species uses synteny_locus)
# --------------------------------------------------------------------------- #
def _build_synteny_backcompat(ctx: GenericContext, ref_sid: str) -> None:
    rows = [r for r in read_tsv(ctx.core("synteny_neighbors.tsv"))
            if not r.get("species_id") or str(r.get("species_id")) == ref_sid]
    cfg = read_json(ctx.run_dir / "run_config.json", {}) or {}
    neighbours = [{**r, "resolved": r.get("status") == "resolved",
                   "distance_to_target": _int(r.get("distance_to_target"))} for r in rows]
    idx = {
        "schema_version": 2, "label": "Local Gene Neighbourhood",
        "scope": "single_species_local_neighbourhood", "comparative_ready": True,
        "gene_symbol": ctx.gene_symbol, "species_id": ref_sid,
        "assembly_accession": (cfg.get("annotation_provenance") or {}).get("assembly_accession", ""),
        "annotation_source": "NCBI genomic GFF3",
        "neighbours": neighbours, "generated_at": _now(),
    }
    write_json(ctx.run_dir / "website_indices" / "synteny_index.json", idx)


def build(ctx: GenericContext) -> Dict[str, Any]:
    ctx.assert_not_freeze()
    wi = ctx.run_dir / "website_indices"
    iso_all = read_tsv(ctx.core("protein_isoform_index.tsv"))
    exon_all = read_tsv(ctx.core("exon_protein_map.tsv"))
    clusters_all = read_tsv(ctx.out("event_region_candidate_clusters.tsv"))
    pair_all = read_tsv(ctx.out("event_region_evidence.tsv"))
    uniprot = read_json(ctx.core("uniprot_event_evidence_report.json"), {}) or {}
    cluster_ready = ctx.cluster_status == "complete" and bool(read_tsv(ctx.core("domain_features.tsv")))
    domains_all = read_tsv(ctx.core("domain_features.tsv")) if cluster_ready else []
    tm_all = read_tsv(ctx.core("tm_features.tsv")) if cluster_ready else []
    domain_context = {r.get("candidate_id", ""): r
                      for r in read_tsv(ctx.out("candidate_domain_context.tsv"))} if cluster_ready else {}
    msa_status = read_json(ctx.out("msa_index.tsv").with_suffix(".status.json"), {}) or {}
    per_species_aln = (read_json(ctx.run_dir / "results" / "07_msa" / "msa_status.json", {}) or {}).get(
        "per_species_isoform_alignments", {})

    species_ids = _species_order(ctx, iso_all)
    primaries = _species_primaries(ctx)

    tx_sections: List[Dict[str, Any]] = []
    cand_sections: List[Dict[str, Any]] = []
    aln_sections: List[Dict[str, Any]] = []
    arch_sections: List[Dict[str, Any]] = []
    total_candidates = 0
    total_transcripts = 0

    for sid in species_ids:
        prim = primaries.get(sid, {})
        primary_pid = prim.get("primary_protein_id", "")
        primary_tid = prim.get("primary_transcript_id", "")
        iso_rows = [r for r in iso_all if str(r.get("species_id") or "") == sid] or (
            iso_all if len(species_ids) == 1 else [])
        exon_rows = [e for e in exon_all if str(e.get("species_id") or "") == sid] or (
            exon_all if len(species_ids) == 1 else [])
        clusters = [c for c in clusters_all if str(c.get("species_id") or "") == sid] or (
            clusters_all if len(species_ids) == 1 else [])

        tx = _transcript_section(sid, primary_pid, primary_tid, exon_rows, iso_rows)
        cand = _candidates_section(ctx, sid, primary_pid, clusters, pair_all, exon_rows,
                                   uniprot, domain_context, cluster_ready)
        aln = _alignment_section(ctx, sid, primary_pid, iso_rows, cand["candidates"], per_species_aln)
        arch = _architecture_section(sid, primary_pid, iso_rows, exon_rows, cand["candidates"],
                                     domains_all, tm_all, cluster_ready)
        tx_sections.append(tx)
        cand_sections.append(cand)
        aln_sections.append(aln)
        arch_sections.append(arch)
        total_candidates += len(cand["candidates"])
        total_transcripts += len(tx["transcripts"])

    ref = species_ids[0] if species_ids else ""
    ref_tx = tx_sections[0] if tx_sections else {"transcripts": []}
    ref_cand = cand_sections[0] if cand_sections else {"candidates": [], "reference_protein": ""}
    ref_aln = aln_sections[0] if aln_sections else {}
    ref_arch = arch_sections[0] if arch_sections else {"proteins": []}

    # ---- transcript_exon_structure_index.json (species-aware) --------------- #
    write_json(wi / "transcript_exon_structure_index.json", {
        "schema_version": 2,
        "status": "available" if total_transcripts else "failed",
        "scientific_question": "How do annotated transcripts differ at exon/CDS level?",
        "gene_symbol": ctx.gene_symbol,
        "n_species": len(species_ids),
        "species": tx_sections,
        "exon_numbering": "Transcript-relative E1..En; genomic equivalence requires shared_exon_group_id.",
        "generated_at": _now(),
        # reference-species back-compat
        "species_id": ref, "primary_transcript_id": ref_tx.get("primary_transcript_id", ""),
        "transcripts": ref_tx.get("transcripts", []),
    })

    # ---- event_candidate_evidence_index.json (species-aware) ---------------- #
    ranking_rows = [{
        "candidate_id": c["candidate_id"], "species_id": c["species_id"],
        "reference_protein": c["reference_protein"], "aa_start": c["aa_start"],
        "aa_end": c["aa_end"], "length": c["length"], "status": c["status"],
        "overall_score": c["overall_score"], "confidence_class": c["confidence_class"],
        **c["score_components"], "confidence_reason": c["confidence_reason"],
    } for sec in cand_sections for c in sec["candidates"]]
    ranking_rows.sort(key=lambda x: (x["species_id"], -x["overall_score"], x["aa_start"]))
    if ranking_rows:
        write_tsv(ctx.out("event_candidate_ranking.tsv"), ranking_rows, list(ranking_rows[0]))
    write_json(wi / "event_candidate_evidence_index.json", {
        "schema_version": 2, "status": "available", "gene_symbol": ctx.gene_symbol,
        "n_species": len(species_ids),
        "domain_context_status": "available" if cluster_ready else "pending_cluster",
        "species": cand_sections,
        "generated_at": _now(),
        # reference-species back-compat
        "reference_protein": ref_cand.get("reference_protein", ""),
        "candidates": ref_cand.get("candidates", []),
    })

    # ---- isoform_alignment_index.json (species-aware) ----------------------- #
    any_aln = any(s.get("status") == "available" for s in aln_sections)
    write_json(wi / "isoform_alignment_index.json", {
        "schema_version": 2,
        "status": "available" if any_aln else "unavailable",
        "kind": "within_species_protein_isoform_alignment",
        "disclaimer": "This is a protein-isoform alignment within one species, not a "
                      "cross-species conservation analysis.",
        "tool": "MAFFT", "version": _mafft_version(),
        "command": "mafft --auto --quiet <per_species_isoform_fasta>",
        "gene_symbol": ctx.gene_symbol, "n_species": len(species_ids),
        "species": aln_sections,
        "generated_at": msa_status.get("generated_at") or _now(),
        # reference-species back-compat
        "reference_sequence": ref_aln.get("reference_sequence", ""),
        "sequence_count": ref_aln.get("sequence_count", 0),
        "alignment_length": ref_aln.get("alignment_length", 0),
        "alignment_file": ref_aln.get("alignment_file", ""),
        "sequences": ref_aln.get("sequences", []),
        "consensus": ref_aln.get("consensus", ""),
        "candidates": ref_aln.get("candidates", []),
    })

    # ---- protein_architecture_index.json (species-aware) -------------------- #
    write_json(wi / "protein_architecture_index.json", {
        "schema_version": 2, "status": "available",
        "stage": "post_cluster" if cluster_ready else "pre_cluster",
        "scientific_question": "Which protein-level features occur along the amino-acid sequence?",
        "domain_annotation_status": "available" if cluster_ready else "pending_cluster",
        "message": "" if cluster_ready else "Domain and TM annotation pending cluster.",
        "gene_symbol": ctx.gene_symbol, "n_species": len(species_ids),
        "species": arch_sections,
        "generated_at": _now(),
        # reference-species back-compat
        "proteins": ref_arch.get("proteins", []),
    })

    # ---- per-species protein/transcript model evidence TSVs ----------------- #
    _write_model_evidence(ctx, species_ids, primaries, iso_all, tx_sections)
    _build_synteny_backcompat(ctx, ref)

    return {
        "n_transcripts": total_transcripts,
        "n_candidates": total_candidates,
        "n_species_explorer": len(species_ids),
        "alignment_status": "available" if any_aln else "unavailable",
        "cluster_ready": cluster_ready,
    }


def _write_model_evidence(ctx: GenericContext, species_ids: List[str],
                          primaries: Mapping[str, Dict[str, Any]],
                          iso_all: List[Mapping[str, Any]],
                          tx_sections: List[Mapping[str, Any]]) -> None:
    """Per-species protein/transcript model evidence tables (generic)."""
    tx_len_by_pid: Dict[str, Any] = {}
    for sec in tx_sections:
        for t in sec.get("transcripts", []):
            tx_len_by_pid[str(t.get("protein_id"))] = t.get("exon_count")
    protein_rows: List[Dict[str, Any]] = []
    transcript_rows: List[Dict[str, Any]] = []
    for sid in species_ids:
        primary_pid = primaries.get(sid, {}).get("primary_protein_id", "")
        iso_rows = [r for r in iso_all if str(r.get("species_id") or "") == sid] or (
            iso_all if len(species_ids) == 1 else [])
        for row in iso_rows:
            pid, tid = str(row.get("protein_id") or ""), str(row.get("transcript_id") or "")
            acc = classify_accession(pid, tid)
            selected = pid == primary_pid
            common = {
                "species_id": sid,
                "annotation_source": acc["source_label"],
                "curation_status": ("curated" if acc["curated"] == "yes"
                                    else "predicted" if acc["curated"] == "no" else "unknown"),
                "protein_length": _int(row.get("protein_length")),
                "exon_count": tx_len_by_pid.get(pid, ""),
                "selected_primary": str(selected).lower(),
                "alternative_model_type": "" if selected else "alternative_protein_model",
            }
            protein_rows.append({"protein_id": pid, "transcript_id": tid, **common})
            transcript_rows.append({"transcript_id": tid, "protein_id": pid, **common})
    if protein_rows:
        write_tsv(ctx.out("protein_model_evidence.tsv"), protein_rows, list(protein_rows[0]))
    if transcript_rows:
        write_tsv(ctx.out("transcript_model_evidence.tsv"), transcript_rows, list(transcript_rows[0]))
