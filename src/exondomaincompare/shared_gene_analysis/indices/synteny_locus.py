"""Build synteny_locus_index.json from the shared canonical synteny contract.

One contract and one renderer serve every gene: the per-species rows produced
here have the same shape whether the run holds one species or thirty, and
whether the gene is FGFR2 or any other. Gene-specific differences travel as
data, never as a second implementation.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from .. import species_order as so
from .. import synteny_contract as sc
from ..common import SharedRunContext, display_species, read_json, read_tsv, rel, to_int
from ..strand import is_reverse


def _neighbour(row: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    side = str(row.get("side") or "").strip()
    if side not in ("upstream", "downstream"):
        return None
    symbol = str(row.get("neighbor_symbol") or row.get("symbol") or "").strip()
    resolved = (row.get("resolved") is True
                or str(row.get("status", "")).lower() == "resolved")
    placeholder = sc.is_placeholder_locus(symbol)
    # A locus with a curated symbol is a resolved neighbour; an NCBI LOC
    # identifier is a real locus with an unavailable symbol, not a failure.
    cls = "placeholder" if placeholder else "exact" if resolved else "unresolved"
    return sc.neighbour_locus(
        side=side,
        rank=to_int(row.get("order") or row.get("rank")) or 0,
        source_symbol=symbol,
        resolved_symbol="" if placeholder else symbol,
        gene_id=str(row.get("gene_id") or ""),
        protein_id=str(row.get("protein_id") or ""),
        strand=row.get("orientation") or row.get("strand") or "",
        orthology_class=cls,
        identity_status="resolved" if resolved else "unresolved",
        distance=row.get("distance_to_target") or row.get("distance"),
        seqid=str(row.get("seqid") or ""),
        start=row.get("genomic_start"),
        end=row.get("genomic_end"),
        method=str(row.get("source") or "ncbi_gff"),
    )


def _derive_target_strand(neighbours: List[Dict[str, Any]]) -> str:
    """Recover the target's transcription direction from the recorded loci.

    The extractor labels a genomically later locus ``downstream`` only on the
    plus strand, so comparing the nearest neighbour on each side against their
    genomic coordinates recovers the strand exactly. Nothing is guessed: when
    the coordinates are missing the strand stays empty and no arrow is drawn.
    """
    up = next((n for n in sorted((x for x in neighbours if x["side"] == "upstream"),
                                 key=lambda x: x["rank"]) if n["genomic_start"] is not None),
              None)
    down = next((n for n in sorted((x for x in neighbours if x["side"] == "downstream"),
                                   key=lambda x: x["rank"]) if n["genomic_start"] is not None),
                None)
    if not up or not down:
        return ""
    return "+" if down["genomic_start"] > up["genomic_start"] else "-"


def _derive_target_span(neighbours: List[Dict[str, Any]],
                        plus: bool) -> Tuple[Optional[int], Optional[int]]:
    """Target start/end reconstructed from the recorded neighbour distances.

    ``distance_to_target`` is the gap between the target and that locus, so the
    nearest locus on each genomic side pins one target edge. This only fills in
    a legacy run that predates ``synteny_target.tsv``; a fresh run reads the
    real coordinates from that table instead.
    """
    # On the plus strand "upstream" is the lower coordinate; on the minus strand
    # transcription runs the other way, so the sides swap.
    lower_side, upper_side = ("upstream", "downstream") if plus else ("downstream", "upstream")
    lower = next((n for n in sorted((x for x in neighbours if x["side"] == lower_side),
                                    key=lambda x: x["rank"])
                  if n["genomic_end"] is not None and n["distance"]), None)
    upper = next((n for n in sorted((x for x in neighbours if x["side"] == upper_side),
                                    key=lambda x: x["rank"])
                  if n["genomic_start"] is not None and n["distance"]), None)
    if not lower or not upper:
        return None, None
    start = lower["genomic_end"] + lower["distance"]
    end = upper["genomic_start"] - upper["distance"]
    return (start, end) if end > start else (None, None)


def _target(sid: str, gene: str, neighbours: List[Dict[str, Any]],
            table: Dict[str, Dict[str, str]],
            gene_ids: Dict[str, str]) -> Dict[str, Any]:
    row = table.get(sid) or {}
    if row:
        return sc.target_locus(
            gene_symbol=row.get("target_symbol") or gene,
            gene_id=row.get("target_gene_id") or gene_ids.get(sid, ""),
            strand=row.get("strand", ""),
            seqid=row.get("seqid", ""),
            start=to_int(row.get("genomic_start")),
            end=to_int(row.get("genomic_end")),
        )
    strand = _derive_target_strand(neighbours)
    start, end = (_derive_target_span(neighbours, not is_reverse(strand))
                  if strand else (None, None))
    seqid = next((n["seqid"] for n in neighbours if n["seqid"]), "")
    return sc.target_locus(
        gene_symbol=gene, gene_id=gene_ids.get(sid, ""), strand=strand,
        seqid=seqid, start=start, end=end,
        coordinate_source=("derived_from_neighbour_offsets" if start is not None
                           else "unavailable"))


def _gene_ids(ctx: SharedRunContext) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for r in read_tsv(ctx.core_dir / "gene_model_index.tsv"):
        sid = str(r.get("species_id") or "")
        if sid and r.get("gene_id"):
            out.setdefault(sid, str(r["gene_id"]))
    return out


def build_synteny_locus_index(ctx: SharedRunContext) -> Dict[str, Any]:
    """Per-species local gene neighbourhood in the shared synteny contract.

    Groups ``synteny_neighbors.tsv`` by species so a multi-species run renders
    one comparative row per species; a single-species run yields one row. The
    target locus comes from ``synteny_target.tsv`` when the run wrote one, and
    is otherwise reconstructed from the recorded neighbour geometry.
    """
    wi = ctx.website_indices
    flat = read_json(wi / "synteny_index.json", {})
    gene = ctx.gene_symbol or str(flat.get("gene_symbol") or "")

    by_species: Dict[str, List[Dict[str, Any]]] = {}
    for r in read_tsv(ctx.core_dir / "synteny_neighbors.tsv"):
        sid = str(r.get("species_id") or "")
        node = _neighbour(r)
        if sid and node:
            by_species.setdefault(sid, []).append(node)

    # Fallback to the flat single-species index if no core table exists.
    if not by_species:
        sid = str(flat.get("species_id") or "unknown_species")
        for n in (flat.get("neighbours") or flat.get("neighbors") or []):
            node = _neighbour(n)
            if node:
                by_species.setdefault(sid, []).append(node)

    target_table = {str(r.get("species_id") or ""): r
                    for r in read_tsv(ctx.core_dir / "synteny_target.tsv")}
    gene_ids = _gene_ids(ctx)

    # The canonical taxonomic order, so the synteny rows line up with every other
    # comparative view of the same dataset.
    ordered_sids = so.order_species(by_species.keys())

    rows: List[Dict[str, Any]] = []
    for sid in ordered_sids:
        neighbours = by_species[sid]
        row = sc.species_row(
            sid, gene_symbol=gene,
            target=_target(sid, gene, neighbours, target_table, gene_ids),
            neighbours=neighbours,
            display_name=display_species(sid),
            taxon_group=so.taxon_group(sid),
            clade=so.clade_of(sid),
            requested=to_int((target_table.get(sid) or {}).get("requested_neighbour_count"))
            or sc.REQUESTED_NEIGHBOUR_COUNT,
            synteny_status="local_neighbourhood",
            comparison_available=False,
        )
        legacy = sc.legacy_nodes(row)
        row.update({
            "synteny_class": "local_neighbourhood",
            "synteny_status_class": "neutral",
            "has_resolved": any(n["resolved"] and not n["is_anchor"] for n in legacy),
            "n_resolved": sum(1 for n in legacy if n["resolved"] and not n["is_anchor"]),
            "n_neighbors": row["displayed_flanking_count"],
            "neighbors5": legacy,
            "neighbors10": legacy,
        })
        rows.append(row)

    # Human is only a comparison reference when it is one of the real species.
    human = next((r for r in rows if r["species_id"] == "homo_sapiens"), None)
    for r in rows:
        r["comparison_available"] = bool(human) and r["species_id"] != "homo_sapiens"

    totals = sc.summarise(rows)
    multi = len(rows) > 1
    source_table = flat.get("source_table") or rel(
        ctx.run_dir, ctx.core_dir / "synteny_neighbors.tsv")
    return {
        "schema_version": 3,
        "contract": "shared_synteny_v1",
        "available": bool(rows and totals["n_flanking_loci"]),
        "synteny_status": "computed" if totals["n_flanking_loci"] else "not_computed",
        "synteny_reason": (
            ("Comparative local genomic neighbourhood, one row per species; "
             "target-centred loci, not a whole-genome conservation claim." if multi else
             "Local genomic neighbourhood around the target gene; not a cross-species "
             "conservation claim.")
            if totals["n_flanking_loci"] else "No synteny neighbour table was produced."
        ),
        "n_resolved_neighbors": totals["n_resolved_neighbours"],
        "n_flanking_loci": totals["n_flanking_loci"],
        "orthology_classes_present": totals["classes_present"],
        "extraction_warning": "",
        "requested_neighbour_count": sc.REQUESTED_NEIGHBOUR_COUNT,
        "has_10neighbor": any(r["displayed_flanking_count"] >= 10 for r in rows),
        "human_reference": human,
        "human_reference_available": bool(human),
        "scope": "comparative_synteny" if multi else "single_species_local_neighbourhood",
        "gene_symbol": gene,
        "target_symbol": gene,
        "assembly_accession": flat.get("assembly_accession", ""),
        "chromosome_or_scaffold": flat.get("chromosome_or_scaffold", ""),
        "source_tables": {
            "neighbours": source_table,
            "target": rel(ctx.run_dir, ctx.core_dir / "synteny_target.tsv"),
        },
        "species_order": so.build_species_order([r["species_id"] for r in rows]),
        "species": rows,
    }
