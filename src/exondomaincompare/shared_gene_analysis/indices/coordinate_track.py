"""Build coordinate_track_index.json from core exon→protein map (FGFR2 CoordinateTrack contract).

The index is enriched with per-exon metadata (genomic/CDS coordinates, phase, strand,
shared exon group) and a per-model track list so the Exon Map can show readable
E1..En labels, full tooltips and an inline "Compare transcript models" panel.
Gene-agnostic: no cassette / IIIb / IIIc assumptions.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from ..common import (
    SharedRunContext, display_species, read_json, read_tsv, rel, shared_exon_group_id, to_int,
)


def _candidate_regions(candidates_blob: Dict[str, Any]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for i, c in enumerate(candidates_blob.get("candidates") or [], start=1):
        if not isinstance(c, dict):
            continue
        out.append({
            "candidate_id": c.get("candidate_id", ""),
            "rank_label": f"C{i}",
            "aa_start": to_int(c.get("aa_start")),
            "aa_end": to_int(c.get("aa_end")),
            "length": to_int(c.get("length")),
            "status": c.get("status", "exploratory_not_validated"),
            "overall_score": c.get("overall_score"),
        })
    return out


def _block_from_structure_exon(e: Dict[str, Any], transcript_id: str) -> Optional[Dict[str, Any]]:
    """A coding exon projected onto protein aa coordinates, with full provenance."""
    start = to_int(e.get("protein_start_aa"))
    end = to_int(e.get("protein_end_aa"))
    if start is None or end is None:
        return None
    exon_num = to_int(e.get("transcript_exon_number"))
    return {
        "feature_type": "coding_exon",
        "id": e.get("exon_id", ""),
        "label": e.get("exon_label") or (f"E{exon_num}" if exon_num else e.get("exon_id", "")),
        "start": start,
        "end": end,
        "exon_number": exon_num,
        "transcript_exon_number": exon_num,
        "transcript_id": transcript_id,
        "shared_exon_group_id": e.get("shared_exon_group_id", ""),
        "genomic_start": to_int(e.get("genomic_start")),
        "genomic_end": to_int(e.get("genomic_end")),
        "cds_start": to_int(e.get("cds_start")),
        "cds_end": to_int(e.get("cds_end")),
        "phase": e.get("phase", ""),
        "strand": e.get("strand", ""),
        "source": e.get("source", ""),
        "coding_status": e.get("coding_status", "coding"),
        "is_iiib_cassette": False,
        "is_iiic_cassette": False,
        "in_cassette": False,
    }


def _blocks_from_map(rows: List[Dict[str, Any]], pid: str, tid: str) -> List[Dict[str, Any]]:
    """Fallback: build blocks from exon_protein_map.tsv.

    ``exon_protein_map.tsv`` carries genomic CDS coordinates and the normalized strand,
    so this path derives the same genomic exon identity as the richer structure index
    rather than leaving it empty. Without it a species built through this path had no
    shared exon group at all, and "Compare transcript models" fell back to the
    protein-scoped exon id — which no two models can ever share, so every model was
    reported as differing from the primary.
    """
    blocks: List[Dict[str, Any]] = []
    for r in rows:
        if r.get("protein_id") != pid:
            continue
        start = to_int(r.get("protein_start_aa"))
        end = to_int(r.get("protein_end_aa"))
        exon_num = to_int(r.get("exon_number"))
        g_s, g_e = to_int(r.get("cds_start")), to_int(r.get("cds_end"))
        strand = str(r.get("normalized_strand") or r.get("strand") or "")
        blocks.append({
            "feature_type": "coding_exon",
            "id": r.get("exon_id", ""),
            "label": f"E{exon_num}" if exon_num else r.get("exon_id", ""),
            "start": start,
            "end": end,
            "exon_number": exon_num,
            "transcript_exon_number": exon_num,
            "transcript_id": r.get("transcript_id", tid),
            "shared_exon_group_id": shared_exon_group_id(g_s, g_e, strand),
            "genomic_start": g_s,
            "genomic_end": g_e,
            "cds_start": g_s,
            "cds_end": g_e,
            "phase": r.get("phase", ""),
            "strand": strand,
            "source": r.get("source", ""),
            "coding_status": "coding",
            "is_iiib_cassette": False,
            "is_iiic_cassette": False,
            "in_cassette": False,
        })
    return blocks


def _panel(pid: str, tid: str, sid: str, plen: Optional[int], blocks: List[Dict[str, Any]],
           *, is_primary: bool, curation: str, candidates: List[Dict[str, Any]],
           source_table: str) -> Dict[str, Any]:
    blocks = sorted(blocks, key=lambda b: (b["start"] is None, b["start"] or 0))
    if plen is None and blocks:
        plen = max((b["end"] or 0) for b in blocks)
    return {
        "species": sid,
        "isoform": "primary" if is_primary else pid,
        "display_species_name": display_species(sid),
        "final_isoform_label": "Primary" if is_primary else "Alternative",
        "protein_id": pid,
        "transcript_id": tid,
        "protein_length": plen,
        "curation_status": curation,
        "is_primary": is_primary,
        "role": "primary" if is_primary else "alternative",
        "cassette_start_aa": None,
        "cassette_end_aa": None,
        "cassette_available": False,
        "boundary_left_precision": "",
        "boundary_right_precision": "",
        "final_plot_status": "primary" if is_primary else "alternative",
        "final_claim_status_after_rescue": "primary" if is_primary else "alternative",
        "claim_class": "accepted" if is_primary else "neutral",
        "readiness_class": "neutral",
        "is_review": False,
        "blocks": blocks,
        "candidate_regions": candidates,
        "source_table": source_table,
    }


def _curation_from_accession(protein_id: str, transcript_id: str) -> str:
    probe = (protein_id or transcript_id or "").upper()
    if probe.startswith(("NM_", "NP_", "NR_")):
        return "curated"
    if probe.startswith(("XM_", "XP_", "XR_")):
        return "predicted"
    return ""


def _species_order(ctx: SharedRunContext, fallback: List[str]) -> List[str]:
    slist = ctx.run_dir / "species_list.txt"
    if slist.is_file():
        order = [ln.strip() for ln in slist.read_text(encoding="utf-8").splitlines() if ln.strip()]
        if order:
            return order
    return fallback


def build_coordinate_track_index(
    ctx: SharedRunContext,
    *,
    primary_protein_id: str = "",
    primary_transcript_id: str = "",
    protein_length: Optional[int] = None,
    species_id: str = "",
) -> Dict[str, Any]:
    """Per-species primary-protein exon map (FGFR2 CoordinateTrack contract).

    Emits one ``species[]`` entry per analysed species, each with its own primary
    protein + per-model tracks. Fully generic (no gene/species assumptions); a
    single-species run yields exactly one species entry.
    """
    exon_map = ctx.core_dir / "exon_protein_map.tsv"
    map_rows = read_tsv(exon_map)
    iso_rows = read_tsv(ctx.core_dir / "protein_isoform_index.tsv")
    wi = ctx.website_indices
    primary_sel = read_json(wi / "primary_selection_index.json", {})
    candidates_blob = read_json(wi / "event_candidate_evidence_index.json", {})
    structure = read_json(wi / "transcript_exon_structure_index.json", {})

    if not map_rows and not (structure.get("transcripts") if isinstance(structure, dict) else None):
        return {"available": False, "domain_layer": "pending_interproscan", "species": [], "models": []}

    all_candidates = _candidate_regions(candidates_blob)
    candidate_reference = str(candidates_blob.get("reference_protein") or "")
    source_table = rel(ctx.run_dir, exon_map)

    # Index transcript_exon_structure by protein_id (richer provenance).
    tx_by_protein: Dict[str, Dict[str, Any]] = {}
    for tx in (structure.get("transcripts") if isinstance(structure, dict) else []) or []:
        tx_by_protein.setdefault(str(tx.get("protein_id") or ""), tx)

    def blocks_for(protein_id: str, transcript_id: str) -> List[Dict[str, Any]]:
        tx = tx_by_protein.get(protein_id)
        if tx:
            out = []
            for e in tx.get("exons", []):
                blk = _block_from_structure_exon(e, tx.get("transcript_id", transcript_id))
                if blk:
                    out.append(blk)
            if out:
                return out
        return _blocks_from_map(map_rows, protein_id, transcript_id)

    # Group protein-coding models by species (from the isoform index, which carries
    # species_id + per-species primary_status + protein_length).
    species_models: Dict[str, List[Dict[str, str]]] = {}
    for r in iso_rows:
        sid_r = str(r.get("species_id") or "")
        species_models.setdefault(sid_r, []).append(r)

    # Reference species (first) drives the top-level fields.
    order = _species_order(ctx, list(species_models.keys()))
    ordered_sids = [s for s in order if s in species_models] or list(species_models.keys())

    species_entries: List[Dict[str, Any]] = []
    top_primary_pid = primary_protein_id or str(primary_sel.get("primary_protein_id") or "")
    top_primary_tid = primary_transcript_id or str(primary_sel.get("primary_transcript_id") or "")
    top_plen = protein_length if protein_length is not None else to_int(primary_sel.get("primary_length_aa"))
    top_models: List[Dict[str, Any]] = []

    for sid in ordered_sids:
        rows = species_models[sid]
        # per-species primary: the row flagged primary, else first
        prim = next((r for r in rows if str(r.get("primary_status", "")).lower() == "primary"), rows[0] if rows else {})
        pid = str(prim.get("protein_id") or "")
        tid = str(prim.get("transcript_id") or "")
        plen = to_int(prim.get("protein_length"))

        # protein order: primary first, then remaining models for this species
        pid_order: List[str] = [pid] if pid else []
        for r in rows:
            rp = str(r.get("protein_id") or "")
            if rp and rp not in pid_order:
                pid_order.append(rp)

        # candidates only overlay the species owning the candidate reference protein
        sp_pids = {str(r.get("protein_id") or "") for r in rows}
        sp_candidates = all_candidates if (candidate_reference in sp_pids or not candidate_reference) else []
        # if the candidate reference does not belong to this species, no overlay
        if candidate_reference and candidate_reference not in sp_pids:
            sp_candidates = []

        models: List[Dict[str, Any]] = []
        for mpid in pid_order:
            mrow = next((r for r in rows if str(r.get("protein_id")) == mpid), {})
            mtid = str(mrow.get("transcript_id") or "")
            if not mtid:
                tx = tx_by_protein.get(mpid)
                mtid = str((tx or {}).get("transcript_id") or "")
            mblocks = blocks_for(mpid, mtid)
            models.append(_panel(
                mpid, mtid, sid, to_int(mrow.get("protein_length")), mblocks,
                is_primary=(mpid == pid),
                curation=_curation_from_accession(mpid, mtid),
                candidates=sp_candidates, source_table=source_table))

        primary_blocks = blocks_for(pid, tid)
        primary_panel = _panel(pid, tid, sid, plen, primary_blocks, is_primary=True,
                               curation=_curation_from_accession(pid, tid),
                               candidates=sp_candidates, source_table=source_table)
        species_entries.append({
            "species": sid,
            "display_species_name": display_species(sid),
            "primary_protein_id": pid,
            "primary_transcript_id": tid,
            "protein_length": plen,
            "panels": {"primary": primary_panel},
            "models": models,
        })
        if sid == ordered_sids[0]:
            top_models = models
            if not top_primary_pid:
                top_primary_pid, top_primary_tid, top_plen = pid, tid, plen

    return {
        "available": bool(species_entries and species_entries[0]["models"]),
        "mode": "primary_protein",
        "domain_layer": "pending_interproscan",
        "gene_symbol": ctx.gene_symbol,
        "n_species": len(species_entries),
        "primary_protein_id": top_primary_pid,
        "primary_transcript_id": top_primary_tid,
        "protein_length": top_plen,
        "species": species_entries,
        "models": top_models,
    }
