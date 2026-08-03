#!/usr/bin/env python3
"""Single-species Gallery figures for transcript structure, genomic context and
exploratory candidates, drawn from the validated protein-coordinate model.

Figures produced per species model (nothing is drawn that the data do not support):

  A  transcript structure and translated protein product   two panels: nucleotides
                                                           and amino acids
  B1 transcript-model comparison, all protein models
  B2 transcript-model comparison, differences from the primary model
  D  local genomic neighbourhood
  E1 exploratory candidate ranking
  E2 exploratory candidates in their domain context

Each figure is exported as SVG, a true-vector PDF, a 300 dpi PNG and the source
table behind it, and registers exactly ONE Gallery card — the formats are formats
of that one card. The five integrated main figures are owned by
``generate_shared_main_figures``; this stage deliberately produces no competing
version of them and retires the cards they replaced.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from exondomaincompare.presentation import figure_captions as fc  # noqa: E402
from exondomaincompare.presentation import shared_gene_plots as sgp  # noqa: E402

FIGURE_DIR = "results/generic_gene_analysis/figures/exon_map"
ANALYSIS_DIR = "results/generic_gene_analysis"

# Cards and files this stage replaces. Only ids are matched, so no other generator's
# cards can be affected; the redesigned figures below take their place. The primary
# exon-to-protein projection is a shared main figure now, so this stage stops drawing
# its own version of it — retiring that card is the shared stage's business.
SUPERSEDED_FIGURE_IDS = (
    "exon_map_{sp}_primary_projection",
    "exon_map_{sp}_selected_candidate_detail",
    "generic_exploratory_event_candidates",
    "generic_synteny_neighbourhood",
)
RETIRED_STEMS = (
    "exon_map_{sp}_primary_projection",
    "exon_map_{sp}_model_comparison_diff",
    "exon_map_{sp}_selected_candidate_detail",
)

# The scientific figure types this stage owns. In a multi-species run each becomes
# one card per species, so the species-independent id has to be retired.
FIGURE_TYPES = (
    "transcript_exon_structure",
    "transcript_model_comparison",
    "transcript_model_comparison_differences",
    "local_gene_neighbourhood",
    "exploratory_candidate_ranking",
    "generic_candidate_domain_context",
)

EXON_SOURCE = f"{ANALYSIS_DIR}/exon_protein_architecture.tsv"


def _rel(p: Path) -> str:
    try:
        return str(p.resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


def _read_tsv(path: Path) -> List[Dict[str, Any]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _candidate_rows(run_dir: Path, model: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Exploratory candidates with their evidence, keyed by amino-acid interval.

    The ranking table carries the evidence scores, the candidate clusters carry the
    supporting comparison counts and the coordinate model carries the affected
    protein models and the display label. All three describe the same candidate, so
    they are merged on the amino-acid interval rather than shown separately.
    """
    ranking = _read_tsv(run_dir / ANALYSIS_DIR / "event_candidate_ranking.tsv")
    clusters = _read_tsv(run_dir / ANALYSIS_DIR / "event_region_candidate_clusters.tsv")
    by_interval = {(sgp._int(c.get("representative_start_aa")),
                    sgp._int(c.get("representative_end_aa"))): c for c in clusters}
    model_regions = {(sgp._int(c.get("start")), sgp._int(c.get("end"))): c
                     for c in (model.get("candidate_regions") or [])}
    protein_id = model.get("protein_id") or ""

    rows: List[Dict[str, Any]] = []
    source = ranking or [
        {"candidate_id": c.get("candidate_cluster_id"),
         "aa_start": c.get("representative_start_aa"),
         "aa_end": c.get("representative_end_aa"),
         "length": c.get("representative_length_aa"),
         "confidence_class": c.get("confidence")} for c in clusters]
    for r in source:
        if r.get("reference_protein") and protein_id \
                and r.get("reference_protein") != protein_id:
            continue
        key = (sgp._int(r.get("aa_start")), sgp._int(r.get("aa_end")))
        cluster = by_interval.get(key, {})
        region = model_regions.get(key, {})
        rows.append({
            "candidate_id": r.get("candidate_id") or cluster.get("candidate_cluster_id"),
            "candidate_label": region.get("id") or "",
            "aa_start": key[0], "aa_end": key[1],
            "length": r.get("length") or cluster.get("representative_length_aa"),
            "affected_proteins": region.get("affected_proteins")
                                 or cluster.get("proteins_involved") or "",
            "support_count": cluster.get("support_count"),
            "overall_score": r.get("overall_score"),
            "confidence_class": r.get("confidence_class") or cluster.get("confidence"),
        })
    return rows


def _footnote(doc: Dict[str, Any], sources: Sequence[str]) -> str:
    """Provenance line under a figure: what the numbers came from."""
    return ("Source: " + fc.join_sources([s for s in sources if s])
            + (f" · run {doc.get('run_id')}" if doc.get("run_id") else ""))


def generate(run_dir: Path, model_json: Path) -> dict:
    sgp.apply_style()
    run_dir = Path(run_dir)
    run_id = run_dir.name
    index = json.loads(model_json.read_text())
    models = index.get("models") or []
    if not models:
        raise SystemExit("no models in coordinate model")
    gene = index.get("gene_symbol") or "GENE"
    fig_dir = run_dir / FIGURE_DIR
    fig_dir.mkdir(parents=True, exist_ok=True)

    manifest: List[Dict[str, Any]] = []
    cards: List[Dict[str, Any]] = []
    drop: List[str] = []

    # In a multi-species run one card per species per figure type is registered, so
    # the card id has to carry the species. A single-species run keeps the bare id
    # the accepted Gallery was validated with.
    multi_species = len(models) > 1
    if multi_species:
        drop += list(FIGURE_TYPES)

    for model in models:
        sp = model.get("species_id") or "sp"
        species = model.get("scientific_name") or ""
        protein_id = model.get("protein_id") or ""
        transcript_id = model.get("transcript_id") or ""
        protein_length = sgp._int(model.get("protein_length"))
        tmodels = model.get("transcript_models") or []
        exons = model.get("exons") or []
        boundaries = model.get("exon_boundaries") or []
        domains = (model.get("representative_domains") or []) \
            if model.get("status") == "available" else []
        candidates = _candidate_rows(run_dir, model)
        classified = sgp.classify_transcript_models(tmodels)
        neighbours = _read_tsv(run_dir / ANALYSIS_DIR / "synteny_neighbourhood.tsv")
        neighbours = [n for n in neighbours
                      if not n.get("species_id") or n.get("species_id") == sp]
        strand = classified[0]["strand"] if classified else "+"

        drop += [t.format(sp=sp) for t in SUPERSEDED_FIGURE_IDS]
        fc.remove_retired_figure_files(fig_dir, [t.format(sp=sp)
                                                for t in RETIRED_STEMS])

        def emit(figure_id: str, stem: str, drawn: bool, *, title: str, category: str,
                 question: str, interpretation: str, caption: str,
                 table: Sequence[Dict[str, Any]] = (), columns: Sequence[str] = (),
                 sources: Sequence[str] = (), kind: str = "main") -> None:
            if not drawn:
                return
            has_table = bool(table) and bool(
                sgp.write_source_table(fig_dir, stem, columns, table))
            fc.write_caption_file(fig_dir, stem, caption)
            card_id = f"{figure_id}__{sp}" if multi_species else figure_id
            cards.append(fc.figure_card(
                figure_id=card_id, figure_type=figure_id,
                title=title, category=category, kind=kind,
                run_id=run_id, figure_dir=FIGURE_DIR, stem=stem,
                scientific_question=question, interpretation=interpretation,
                caption=caption, gene_symbol=gene, species=species, species_id=sp,
                protein_id=protein_id, transcript_id=transcript_id,
                has_table=has_table, source_files=list(sources)))
            for ext in ("svg", "pdf", "png", "tsv"):
                path = fig_dir / f"{stem}.{ext}"
                if path.exists():
                    manifest.append({"group": category, "kind": figure_id,
                                     "title": title, "format": ext,
                                     "protein_id": protein_id, "path": _rel(path)})

        # ---- A. transcript structure and translated protein product ---------
        stem = f"exon_map_{sp}_transcript_and_protein_structure"
        emit("transcript_exon_structure", stem,
             sgp.plot_transcript_exon_structure(
                 fig_dir, stem, gene_symbol=gene, transcripts=tmodels,
                 species_name=species,
                 footnote=_footnote(index, ["NCBI RefSeq GFF CDS features",
                                            "coordinate model"])),
             title=f"{gene} · Transcript structure and translated protein product",
             category="Exon structure",
             question="Which coding exons does each annotated transcript use, and "
                      "which part of the translated protein do they encode?",
             interpretation="Panel A shows the annotated transcripts on genomic "
                            "nucleotide coordinates; panel B shows the translated "
                            "protein models on amino-acid coordinates. Differences "
                            "between models are annotation differences, not "
                            "validated splicing events.",
             caption=fc.build_caption(
                 gene=gene, species=species, protein_id=protein_id,
                 description=f"transcript structure of {len(tmodels)} annotated "
                             f"transcript(s) on genomic coordinates and the "
                             f"corresponding translated protein models on "
                             f"amino-acid coordinates",
                 coordinate_system="panel A " + fc.GENOMIC_COORDINATE_SYSTEM
                                   + "; panel B " + fc.COORDINATE_SYSTEM,
                 annotation_source="NCBI RefSeq GFF CDS features",
                 status=fc.DESCRIPTIVE_STATUS),
             table=sgp.exon_identity_table(classified),
             columns=sgp.EXON_IDENTITY_TABLE_COLUMNS,
             sources=[EXON_SOURCE])

        # ---- B. transcript-model comparison ---------------------------------
        top_candidate = None
        for region in model.get("candidate_regions") or []:
            top_candidate = region
            break
        for figure_id, suffix, diff_only, mode_title, mode_note in (
            ("transcript_model_comparison", "model_comparison_all", False,
             "all protein models",
             "Every annotated protein model is shown, so shared exon structure is "
             "visible alongside the differences."),
            ("transcript_model_comparison_differences", "model_comparison_differences",
             True, "differences from the primary model",
             "Only models whose exon set really differs from the primary are kept."),
        ):
            stem = f"exon_map_{sp}_{suffix}"
            emit(figure_id, stem,
                 sgp.plot_transcript_model_comparison(
                     fig_dir, stem, gene_symbol=gene, models=tmodels,
                     candidate=top_candidate, diff_only=diff_only,
                     species_name=species,
                     footnote=_footnote(index, ["NCBI RefSeq GFF CDS features"])),
                 title=f"{gene} · Transcript-model comparison ({mode_title})",
                 category="Exon structure",
                 question="Which coding exons do the alternative protein models "
                          "share with the primary model, and where do they diverge?",
                 interpretation=f"{mode_note} Exon identity is compared on genomic "
                                f"CDS coordinates, so an upstream deletion does not "
                                f"mark genomically identical downstream exons as "
                                f"altered.",
                 caption=fc.build_caption(
                     gene=gene, species=species, protein_id=protein_id,
                     description=f"coding-exon identity of {len(tmodels)} protein "
                                 f"model(s) relative to the primary model "
                                 f"({mode_title})",
                     annotation_source="NCBI RefSeq GFF CDS features compared on "
                                       "genomic CDS intervals",
                     status=fc.DESCRIPTIVE_STATUS),
                 table=sgp.exon_identity_table(classified),
                 columns=sgp.EXON_IDENTITY_TABLE_COLUMNS,
                 sources=[EXON_SOURCE])

        # ---- D. local genomic neighbourhood --------------------------------
        stem = f"exon_map_{sp}_local_gene_neighbourhood"
        layout = sgp.neighbourhood_layout(gene, neighbours, strand)
        emit("local_gene_neighbourhood", stem,
             sgp.plot_synteny_neighbourhood(
                 fig_dir, stem, gene_symbol=gene, neighbours=neighbours,
                 species_name=species, target_orientation=strand,
                 footnote=_footnote(index, ["NCBI RefSeq gene annotation of this "
                                            "assembly"])),
             title=f"{gene} · Local genomic neighbourhood",
             category="Genomic context",
             question=f"Which annotated loci flank {gene} in this assembly, and in "
                      f"which direction is each transcribed?",
             interpretation="The target gene is centred between its annotated "
                            "upstream and downstream loci, each drawn in its "
                            "annotated transcription direction. Loci without an "
                            "approved symbol are shown as placeholder loci. Spacing "
                            "is ordinal, not to scale.",
             caption=fc.build_caption(
                 gene=gene, species=species, protein_id="",
                 description=f"{len(layout) - 1} annotated flanking loci around "
                             f"{gene} with their transcription directions",
                 coordinate_system="ordinal gene order along the assembly",
                 annotation_source="NCBI RefSeq gene annotation",
                 status=fc.DESCRIPTIVE_STATUS),
             table=[{
                 "position_in_figure": i + 1, "side": locus["side"],
                 "order_from_target": locus["order"],
                 "locus_symbol": locus["symbol"], "locus_kind": locus["kind"],
                 "transcription_direction": locus["orientation"],
                 "annotation_source": locus["source"] or "ncbi_gff",
             } for i, locus in enumerate(layout)],
             columns=sgp.NEIGHBOURHOOD_TABLE_COLUMNS,
             sources=[f"{ANALYSIS_DIR}/synteny_neighbourhood.tsv"])

        # ---- E1. exploratory candidate ranking ------------------------------
        stem = f"exon_map_{sp}_exploratory_candidate_ranking"
        emit("exploratory_candidate_ranking", stem,
             sgp.plot_evidence_regions_on_protein(
                 fig_dir, stem, gene_symbol=gene, clusters=candidates,
                 max_aa=protein_length, species_name=species,
                 protein_id=protein_id, n_models=len(tmodels),
                 footnote=_footnote(index, ["protein-isoform difference scan",
                                            "exon/protein boundary check"])),
             title=f"{gene} · Exploratory candidate ranking",
             category="Exploratory candidates",
             question="Which isoform-difference regions carry the most transparent "
                      "support?",
             interpretation=f"{sgp.EXPLORATORY_TAG} regions ranked by evidence "
                            f"score, with the amino-acid interval, the affected "
                            f"protein models and the supporting comparisons. "
                            f"{sgp.VALIDATION_TAG}: the score measures support, not "
                            f"biological validation.",
             caption=fc.build_caption(
                 gene=gene, species=species, protein_id=protein_id,
                 description=f"{len(candidates)} exploratory candidate region(s) "
                             f"ranked by evidence score, with affected isoform and "
                             f"supporting-comparison counts",
                 annotation_source="protein-isoform difference scan with "
                                   "exon-boundary and domain-context checks",
                 status=f"{sgp.EXPLORATORY_TAG} regions. {sgp.VALIDATION_TAG}."),
             table=sgp.candidate_ranking_rows(candidates),
             columns=sgp.CANDIDATE_RANKING_TABLE_COLUMNS,
             sources=[f"{ANALYSIS_DIR}/event_candidate_ranking.tsv",
                      f"{ANALYSIS_DIR}/event_region_candidate_clusters.tsv"])

        # ---- E2. exploratory candidates in their domain context -------------
        stem = f"exon_map_{sp}_candidate_domain_context"
        context_rows = [{
            "candidate_label": (r.get("candidate_label")
                                or f"C{i}"),
            "candidate_id": r.get("candidate_id") or "",
            "aa_start": r["aa_start"], "aa_end": r["aa_end"],
        } for i, r in enumerate(sgp.candidate_ranking_rows(candidates), start=1)]
        emit("generic_candidate_domain_context", stem,
             sgp.plot_candidate_domain_context(
                 fig_dir, stem, gene_symbol=gene, contexts=context_rows,
                 domains=domains, exon_blocks=exons, boundaries=boundaries,
                 protein_id=protein_id, protein_length=protein_length,
                 species_name=species,
                 footnote=_footnote(index, ["InterProScan representative domain "
                                            "instances", "GFF CDS coding exons"])),
             title=f"{gene} · Exploratory candidates in their domain context",
             category="Exploratory candidates",
             question="Do the exploratory candidate regions fall inside annotated "
                      "domains or between them?",
             interpretation=f"Each {sgp.EXPLORATORY_TAG.lower()} region is shown "
                            f"against the coding exons, the internal coding-exon "
                            f"boundaries and the representative domain instances it "
                            f"overlaps, with the distance to the nearest domain "
                            f"edge. Domains are resolved per instance, so repeated "
                            f"domains stay distinct. {sgp.VALIDATION_TAG}: a "
                            f"positional overlap is an observation, not a "
                            f"functional claim.",
             caption=fc.build_caption(
                 gene=gene, species=species, protein_id=protein_id,
                 description=f"{len(context_rows)} exploratory candidate region(s) "
                             f"against {len(domains)} representative domain "
                             f"instance(s) and the coding exons of {protein_id}",
                 annotation_source="InterProScan representative domain instances "
                                   "(resolved per instance) + GFF CDS exons",
                 status=f"{sgp.EXPLORATORY_TAG} regions. {sgp.VALIDATION_TAG}."),
             table=sgp.candidate_context_table(context_rows, domains, boundaries),
             columns=sgp.CANDIDATE_CONTEXT_TABLE_COLUMNS,
             sources=[f"{ANALYSIS_DIR}/event_candidate_ranking.tsv",
                      "results/core_gene_analysis/domain_features.tsv"])

    # The served coordinate model carries the export manifest for the Gene Explorer
    # export bar; cards other generators registered are preserved.
    own_kinds = {m["kind"] for m in manifest}
    other = [f for f in (index.get("publication_figures") or [])
             if f.get("group") != "exon_map" and f.get("kind") not in own_kinds]
    index["publication_figures"] = manifest + other
    model_json.write_text(json.dumps(index, indent=2))
    fc.register_gallery_cards(run_dir, cards, drop_figure_ids=drop)
    return {"figures": len(manifest), "cards": len(cards), "manifest": manifest}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--model", type=Path, default=None,
                    help="coordinate model JSON (default: "
                         "<run>/website_indices/generic/protein_coordinate_model.json)")
    args = ap.parse_args()
    run_dir = args.run_dir.resolve()
    model_json = args.model or (run_dir / "website_indices" / "generic"
                                / "protein_coordinate_model.json")
    if not model_json.exists():
        raise SystemExit(f"coordinate model not found: {model_json}")
    res = generate(run_dir, model_json)
    print(f"OK — wrote {res['figures']} figure file(s) in {res['cards']} Gallery "
          f"card(s) for {run_dir.name}")
    for f in res["manifest"]:
        if f["format"] == "svg":
            print(f"  · {f['kind']}: {f['path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
