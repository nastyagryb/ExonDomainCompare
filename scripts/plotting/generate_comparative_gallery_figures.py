#!/usr/bin/env python3
"""Generate and register Comparative Figure Gallery cards for multi-species runs.

Renders through the shared JS figure builders (``comparativeGalleryFigures.js`` /
``comparativeFigures.js``) so Gallery cards and Explorer exports stay one
implementation. Single-species runs produce nothing — empty comparative figures
are never registered.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from plotting.generate_shared_main_figures import (  # noqa: E402
    EXPORT_DPI, _api, _rasterise_png, _rel,
)
from exondomaincompare.runs.registry import RegistryError, resolve_run_record  # noqa: E402
from exondomaincompare.config import load_config  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=ROOT)

RENDERER = ROOT / "scripts" / "plotting" / "render_comparative_gallery_figures.mjs"
FIGURE_SUBDIR = Path("results") / "generic_gene_analysis" / "figures" / "comparative"
GROUP = "comparative_figures"
SCOPE = "comparative"

# The final comparative inventory (Part 5). Insertion order is the Gallery
# reading order. Every entry is a scientific visualisation: the three text-only
# cards that used to sit here — the standalone pairwise-identity page, the
# exon-boundary alignment summary page and the textual isoform-diversity page —
# were removed. Their numbers now live inside the MSA overview, in the card
# source tables and in the workbook.
RETIRED_FIGURE_IDS = (
    "cmp_pairwise_identity",
    "cmp_exon_boundary_alignment_summary",
)

#: Cards withdrawn from one gene's catalogue rather than from the shared figure set.
#: The per-neighbour synteny conservation matrix asks the same question as the main
#: synteny neighbourhood figure and its blank cells are annotation gaps in single
#: assemblies, so the FGFR2 catalogue no longer offers it as a reader's entry point.
#: The figure is still rendered and exported; only the card is withheld, and only for
#: the gene whose catalogue withdrew it — silently dropping it for every gene would
#: change galleries nobody asked to change.
WITHDRAWN_BY_GENE: Dict[str, tuple] = {
    "FGFR2": ("cmp_synteny_neighbour_conservation",),
}


def _withdrawn_comparative_stems(gene_symbol: str) -> set:
    return set(WITHDRAWN_BY_GENE.get(str(gene_symbol or "").upper(), ()))

FIGURE_META: Dict[str, Dict[str, str]] = {
    "cmp_msa_aligned_exon_architecture": {
        "title": "MSA-aligned exon architecture",
        "category": "Comparative exon structure",
        "kind": "main",
        "question": "Where do the coding exons of the selected primary proteins fall "
                    "on a common alignment coordinate, and which exon intervals are "
                    "matched across species?",
        "interpretation": "Coding exons are projected onto MSA columns. Exons whose "
                          "aligned interval is matched in every species carry the "
                          "shared-exon colour; the rest carry the alternative-exon "
                          "colour. A shared column means residues were aligned, not "
                          "that they are functionally equivalent.",
    },
    "cmp_native_exon_architecture": {
        "title": "Native-coordinate exon architecture comparison",
        "category": "Comparative exon structure",
        "kind": "supplement",
        "question": "How long is each primary protein, and how are its coding exons "
                    "and exon boundaries distributed along its own sequence?",
        "interpretation": "Native coordinates are not homologous alignment "
                          "coordinates: the same x position in two lanes is not the "
                          "same residue. Use the MSA-aligned exon architecture for "
                          "positional comparison.",
    },
    "cmp_primary_msa_overview": {
        "title": "Cross-species primary-protein MSA overview",
        "category": "Comparative sequence analysis",
        "kind": "main",
        "question": "Where do the selected primary proteins agree, where do residues "
                    "differ, and where does one species carry an insertion or "
                    "deletion?",
        "interpretation": "Column-resolved alignment summary: per-species coverage, "
                          "identical / mismatch / indel columns, a windowed "
                          "conservation curve and the variable blocks. Pairwise "
                          "identity is reported as a metric of this figure and is "
                          "computed on ungapped columns only.",
    },
    "cmp_domain_architecture_msa": {
        "title": "Comparative domain architecture (MSA-aligned coordinates)",
        "category": "Comparative domain architecture",
        "kind": "main",
        "question": "Which representative domain instances does each primary protein "
                    "carry, and do they occupy comparable aligned positions?",
        "interpretation": "Domain instances are projected onto MSA columns and "
                          "labelled individually; repeated domains of one accession "
                          "stay separate instances. MSA intervals are geometry from "
                          "the primary-protein alignment; they do not assert "
                          "functional equivalence.",
    },
    "cmp_domain_architecture_native": {
        "title": "Comparative domain architecture (native coordinates)",
        "category": "Comparative domain architecture",
        "kind": "main",
        "question": "Which representative domain instances and membrane segments does "
                    "each primary protein carry on its own amino-acid axis?",
        "interpretation": "The native-coordinate counterpart of the MSA-aligned "
                          "panel, with the transmembrane helix included. Native axes "
                          "are species-specific and not directly comparable position "
                          "by position.",
    },
    "cmp_exon_domain_architecture_native": {
        "title": "Exon and domain architecture per species (native coordinates)",
        "category": "Comparative exon–domain architecture",
        "kind": "main",
        "question": "For each species, where do the coding exon boundaries of the "
                    "primary protein fall relative to the edges of its annotated "
                    "domains?",
        "interpretation": "Each species owns a two-track group on its own amino-acid "
                          "axis: domain instances above, coding exons below, with a "
                          "connector carrying every exon boundary through the domain "
                          "track. Where the boundary analysis classified a boundary, "
                          "the connector carries that class colour. Native axes are "
                          "species-specific and not comparable position by position.",
    },
    "cmp_exon_domain_architecture_msa": {
        "title": "Exon and domain architecture per species (MSA-aligned coordinates)",
        "category": "Comparative exon–domain architecture",
        "kind": "main",
        "question": "On a common alignment axis, where do each species' exon "
                    "boundaries fall relative to its domain edges?",
        "interpretation": "The MSA-aligned counterpart of the native panel, built by "
                          "the same code from the same data. Both tracks of a species "
                          "share one alignment axis, so the connectors are directly "
                          "comparable between species. An aligned column means "
                          "residues were aligned, not that they are equivalent.",
    },
    "cmp_domain_annotation_matrix": {
        "title": "Domain annotation matrix",
        "category": "Comparative domain architecture",
        "kind": "supplement",
        "question": "Which comparable domain instance is detected in which species, "
                    "and where is the annotation still pending?",
        "interpretation": "States distinguish detected, not detected, pending, "
                          "unavailable and uncertain mapping; detected cells carry the "
                          "native interval. ‘Not detected’ is not biological absence.",
    },
    "cmp_boundary_matrix": {
        "title": "Comparative exon–domain boundary matrix",
        "category": "Comparative exon–domain boundaries",
        "kind": "main",
        "question": "Which exon–domain boundaries are comparable across species, and "
                    "at what signed distance to a domain edge?",
        "interpretation": "Cells mirror the canonical comparative index used by the "
                          "interactive Comparative Boundary Explorer. Tentative "
                          "mappings are marked and are not treated as confirmed pairs.",
    },
    "cmp_paired_signed_distance": {
        "title": "Paired signed-distance plot",
        "category": "Comparative exon–domain boundaries",
        "kind": "main",
        "question": "Do the species place the same comparable boundary at the same "
                    "distance, and on the same side, of the nearest domain edge?",
        "interpretation": "One marker per species per comparable group. A solid "
                          "connector marks a supported mapping; a dotted connector "
                          "marks a tentative one and asserts no equivalence.",
    },
    "cmp_boundary_position_consistency": {
        "title": "Boundary-position consistency",
        "category": "Comparative exon–domain boundaries",
        "kind": "main",
        "question": "Where does each species place a comparable boundary, how far "
                    "apart are those placements, and how well supported is the "
                    "comparison?",
        "interpretation": "Three aligned panels: the raw signed distances, the "
                          "cross-species difference, and class agreement with mapping "
                          "confidence and species coverage. A small difference is "
                          "evidence of a consistent boundary position, not of "
                          "evolutionary conservation.",
    },
    "cmp_local_boundary_architecture": {
        "title": "Comparative local boundary architecture",
        "category": "Comparative exon–domain boundaries",
        "kind": "supplement",
        "question": "How does the local exon and domain context around a comparable "
                    "boundary look in each species?",
        "interpretation": "Tracks are aligned on the selected boundary itself, so the "
                          "domain edges either line up or visibly do not. Tentative "
                          "mapping does not imply equivalence.",
    },
    "cmp_comparative_synteny": {
        "title": "Comparative local genomic context",
        "category": "Comparative genomic context",
        "kind": "main",
        "question": "Do the species share the same immediate gene neighbourhood "
                    "around the target locus, and in which orientation?",
        "interpretation": "The target locus occupies its own central slot and is "
                          "never counted as a neighbour. Loci are placed by "
                          "annotated neighbour rank, never by a fabricated "
                          "coordinate; where a species provides fewer real loci the "
                          "slot stays empty. A matching symbol is a nomenclature "
                          "match, not an orthology assignment.",
    },
    "cmp_synteny_neighbour_conservation": {
        "title": "Neighbour conservation across species",
        "category": "Comparative genomic context",
        "kind": "supplement",
        "question": "Which flanking loci recur across the analysed species, and "
                    "where is the neighbourhood incomplete?",
        "interpretation": "A filled cell means the annotation of that species "
                          "carries the same gene symbol somewhere in the flanking "
                          "window. A shared symbol supports locus context; it is "
                          "not an orthology assignment.",
    },
    "cmp_isoform_diversity": {
        "title": "Comparative isoform diversity",
        "category": "Comparative isoform diversity",
        "kind": "main",
        "question": "How many protein models does each species have, how are they "
                    "curated, how much do their lengths vary, and how many "
                    "exploratory candidates were derived from them?",
        "interpretation": "Model counts are annotation counts, not validated splice "
                          "products. Complements the within-species isoform "
                          "alignments; it does not replace them.",
    },
}


# Which comparative artefact holds the numbers behind each card.
SOURCE_ARTEFACT: Dict[str, str] = {
    "cmp_msa_aligned_exon_architecture": "msa_aligned_exons",
    "cmp_native_exon_architecture": "msa_aligned_exons",
    "cmp_primary_msa_overview": "pairwise_identity",
    "cmp_domain_architecture_msa": "msa_aligned_domains",
    "cmp_domain_architecture_native": "msa_aligned_domains",
    "cmp_exon_domain_architecture_native": "msa_aligned_domains",
    "cmp_exon_domain_architecture_msa": "msa_aligned_domains",
    "cmp_domain_annotation_matrix": "domain_annotation_matrix",
    "cmp_boundary_matrix": "exon_domain_boundaries_long",
    "cmp_paired_signed_distance": "exon_domain_boundaries_long",
    "cmp_boundary_position_consistency": "boundary_consistency_summary",
    "cmp_local_boundary_architecture": "exon_domain_boundaries_long",
    "cmp_comparative_synteny": "comparative_synteny",
    "cmp_synteny_neighbour_conservation": "comparative_synteny",
    "cmp_isoform_diversity": "isoform_diversity",
}


def _render(model_json: Path, comparative_json: Path, out_dir: Path) -> List[dict]:
    if shutil.which("node") is None:
        raise SystemExit("node is required to render comparative gallery figures")
    result = subprocess.run(
        ["node", str(RENDERER), str(model_json), str(comparative_json), str(out_dir)],
        capture_output=True, text=True, cwd=ROOT,
    )
    if result.returncode != 0:
        raise SystemExit(
            f"comparative gallery renderer failed:\n{result.stdout}\n{result.stderr}")
    print(result.stdout.rstrip())
    summary = out_dir / "render_summary.json"
    return json.loads(summary.read_text()) if summary.exists() else []


def _register_gallery(run_dir: Path, cards: List[Dict[str, Any]]) -> int:
    new_ids = {c["figure_id"] for c in cards} | set(RETIRED_FIGURE_IDS)
    written = 0
    for name in ("figures_index.json", "generic/figures_index.json"):
        fp = run_dir / "website_indices" / name
        if not fp.exists():
            continue
        try:
            doc = json.loads(fp.read_text())
        except (OSError, ValueError):
            continue
        kept = [f for f in (doc.get("figures") or [])
                if f.get("figure_id") not in new_ids
                and f.get("scope") != SCOPE]
        # Comparative cards lead the catalogue for multi-species reading order.
        doc["figures"] = cards + kept
        avail = doc.get("available")
        if isinstance(avail, list):
            def _id(a):
                return a.get("figure_id") if isinstance(a, dict) else a
            doc["available"] = (
                [c["figure_id"] for c in cards]
                + [a for a in avail if _id(a) not in new_ids]
            )
        fp.write_text(json.dumps(doc, indent=2))
        written += 1
    return written


def generate(run_dir: Path, model_json: Path) -> dict:
    from exondomaincompare.shared_gene_analysis.comparative_dataset import build_comparative_dataset

    run_dir = Path(run_dir)
    model_json = Path(model_json)
    index = json.loads(model_json.read_text())
    models = index.get("models") or []
    if len(models) < 2:
        # A one-species dataset has nothing to compare, so it registers no
        # comparative card at all rather than a cross-species card with one row.
        return {"figures": 0, "cards": 0, "skipped": "single_species"}

    comparative = build_comparative_dataset(run_dir, coordinate_index=index)
    cmp_path = (run_dir / "website_indices" / "generic"
                / "comparative_dataset_index.json")

    gene = index.get("gene_symbol") or "GENE"
    run_id = run_dir.name
    fig_dir = run_dir / FIGURE_SUBDIR
    fig_dir.mkdir(parents=True, exist_ok=True)
    # A retired card must not leave a downloadable file behind.
    from exondomaincompare.presentation.figure_captions import remove_retired_figure_files
    remove_retired_figure_files(fig_dir, RETIRED_FIGURE_IDS)
    rendered = _render(model_json, cmp_path, fig_dir)

    cards: List[Dict[str, Any]] = []
    manifest: List[Dict[str, Any]] = []
    warnings: List[str] = []
    proteins = ", ".join(
        f"{m.get('scientific_name') or m.get('species_id')} "
        f"({m.get('protein_id')})" for m in models)

    withdrawn = _withdrawn_comparative_stems(gene)

    for entry in rendered:
        stem = entry["stem"]
        if stem in withdrawn:
            continue
        meta = FIGURE_META.get(stem)
        if meta is None:
            warnings.append(f"unrecognised comparative stem: {stem}")
            continue
        svg = fig_dir / f"{stem}.svg"
        pdf = fig_dir / f"{stem}.pdf"
        png = fig_dir / f"{stem}.png"
        tsv = fig_dir / f"{stem}.tsv"
        if not _rasterise_png(svg, png):
            warnings.append(f"could not rasterise {stem}.png at {EXPORT_DPI} dpi")

        rel = {}
        for ext, path in (("svg", svg), ("pdf", pdf), ("png", png), ("tsv", tsv)):
            if path.exists():
                rel[ext] = _rel(path).split(f"runs/{run_id}/", 1)[-1]
                manifest.append({
                    "group": GROUP, "kind": stem, "title": meta["title"],
                    "format": ext, "scope": SCOPE, "path": _rel(path),
                })

        arts = comparative.get("artefacts") or {}
        source = SOURCE_ARTEFACT.get(stem, "")
        source = (arts.get(source) or "") if source else ""
        if stem == "cmp_primary_msa_overview" and not source:
            source = (comparative.get("msa") or {}).get("alignment_file") or ""

        cards.append({
            "figure_id": stem,
            "figure_type": stem,
            "title": f"{gene} · {meta['title']}",
            "category": meta["category"],
            "section": meta["category"],
            "kind": meta.get("kind", "main"),
            "scope": SCOPE,
            "species_id": "",
            "species": "comparative",
            "scientific_question": meta["question"],
            "interpretation": meta["interpretation"],
            "caption": (f"{gene}. Comparative scope across "
                        f"{comparative.get('n_species')} species. {proteins}. "
                        f"{meta['interpretation']}"),
            "gene_symbol": gene,
            "protein_id": "",
            "proteins": proteins,
            "stage": "post_cluster",
            "status": "available",
            "data_availability": "available",
            "analysis_status": "available",
            "renderer": "shared_comparative_figure_specification",
            "width_pt": entry.get("width"),
            "height_pt": entry.get("height"),
            "png_url": _api(run_id, rel["png"], inline=True) if "png" in rel else "",
            "svg_url": _api(run_id, rel["svg"]) if "svg" in rel else "",
            "pdf_url": _api(run_id, rel["pdf"]) if "pdf" in rel else "",
            "table_url": _api(run_id, rel["tsv"]) if "tsv" in rel else "",
            "source_table": source,
            "source_files": [source] if source else [],
            "error": "",
        })
        if entry.get("warnings"):
            warnings.extend(f"{stem}: {w}" for w in entry["warnings"])

    order = list(FIGURE_META)
    cards.sort(key=lambda c: order.index(c["figure_id"]) if c["figure_id"] in order else 999)

    existing = [f for f in (index.get("publication_figures") or [])
                if f.get("group") != GROUP]
    index["publication_figures"] = existing + manifest
    model_json.write_text(json.dumps(index, indent=2))
    n_idx = _register_gallery(run_dir, cards)

    n_files = sum(1 for c in cards for ext in ("svg", "pdf", "png")
                  if (fig_dir / f"{c['figure_id']}.{ext}").exists())
    return {
        "figures": n_files,
        "cards": len(cards),
        "indices_updated": n_idx,
        "warnings": warnings,
        "comparative_available": comparative.get("available"),
    }


def main(argv=None) -> int:
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-id", required=True)
    p.add_argument("--runs-root", type=Path, default=None)
    args = p.parse_args(argv)
    if args.runs_root is not None:
        run_dir = Path(args.runs_root).expanduser().resolve() / args.run_id
    else:
        try:
            record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
        except RegistryError as exc:
            raise SystemExit(str(exc)) from exc
        if record is None:
            raise SystemExit(f"run not found: {args.run_id}")
        if record.read_only:
            raise SystemExit("run is registered read-only; copy it before rendering")
        run_dir = record.path
    model = run_dir / "website_indices" / "generic" / "protein_coordinate_model.json"
    if not model.is_file():
        raise SystemExit(f"coordinate model missing: {model}")
    res = generate(run_dir, model)
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
