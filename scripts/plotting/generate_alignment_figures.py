#!/usr/bin/env python3
"""Generate the Figure Gallery's isoform-analysis figures.

The three figures are rendered by the *same* builders the Gene Explorer's alignment
view exports through: this stage shells out to
``scripts/plotting/render_alignment_figures.mjs``, which imports
``webapp/frontend/src/pages/viewers/alignmentFigure.js``. A Gallery figure and the
corresponding Gene Explorer export are therefore one figure with one
implementation, not two renderers kept in agreement by hand.

Figures produced per species with an available alignment:
  1. full isoform-alignment overview     — every column, every model
  2. wrapped residue-level alignment     — one multi-page vector PDF
  3. candidate-associated alignment detail (only when a candidate maps to columns)

Formats: standalone SVG, true vector PDF, 300 dpi PNG, and the alignment summary
table. The wrapped alignment additionally ships one SVG per page.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.plotting.figure_captions import (  # noqa: E402
    build_caption, figure_card, register_gallery_cards, write_caption_file,
)
from scripts.plotting.generate_shared_main_figures import (  # noqa: E402
    EXPORT_DPI, _rasterise_png,
)

RENDERER = ROOT / "scripts" / "plotting" / "render_alignment_figures.mjs"
FIGURE_SUBDIR = Path("results") / "generic_gene_analysis" / "figures" / "main"
CATEGORY = "Isoform analysis"
GROUP = "alignment_figures"

# The card that used to represent the alignment: a coverage bar per isoform in an
# arbitrary categorical colour, with no residues, no conservation and no identity.
# The three figures below replace it, so it is retired rather than left in place.
SUPERSEDED_FIGURE_IDS = {"isoform_alignment", "generic_isoform_alignment"}


def _rel(path: Path) -> str:
    """Repository-relative where possible; a run directory may live anywhere."""
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)

FIGURE_META = {
    "full_isoform_alignment": {
        "title": "Isoform alignment overview",
        "question": "Where along the primary protein do the alternative isoforms "
                    "differ from it, and how much of each isoform aligns?",
        "interpretation": "Colour encodes alignment state only: matches are neutral, "
                          "residues differing from the primary protein are "
                          "highlighted, and gaps are left light. Differences are "
                          "alignment observations, not annotated events.",
        "description": "complete multiple alignment of all protein isoforms against "
                       "the primary protein, with difference, gap-density and "
                       "conservation tracks",
    },
    "wrapped_alignment": {
        "title": "Residue-level isoform alignment (wrapped)",
        "question": "Which individual residues differ between the isoforms across "
                    "the whole alignment?",
        "interpretation": "The alignment wrapped into blocks so every residue stays "
                          "legible in print. The card previews the first page; the "
                          "PDF holds the complete alignment across all pages.",
        "description": "the complete alignment at residue resolution, wrapped into "
                       "blocks over several pages",
    },
    "candidate_alignment_detail": {
        "title": "Candidate-associated alignment detail",
        "question": "Which isoforms carry the exploratory candidate interval, and "
                    "what happens to the residues inside it?",
        "interpretation": "Rows are labelled by their role for this candidate and "
                          "identity is computed inside the candidate interval only. "
                          "The candidate is an exploratory interval, not a validated "
                          "splice event.",
        "description": "the residues of the exploratory candidate interval across all "
                       "isoforms, with their exon context",
    },
}

ORDER = list(FIGURE_META)


def _render(index_json: Path, out_dir: Path, species_id: str = "") -> dict:
    if shutil.which("node") is None:
        raise SystemExit("node is required to render the alignment figures")
    cmd = ["node", str(RENDERER), str(index_json), str(out_dir), "--gallery-stems"]
    if species_id:
        cmd.append(f"--species={species_id}")
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=ROOT,
    )
    if result.returncode != 0:
        raise SystemExit(f"alignment renderer failed:\n{result.stdout}\n{result.stderr}")
    summary = out_dir / "alignment_render_summary.json"
    if not summary.exists():
        raise SystemExit("the alignment renderer wrote no summary")
    return json.loads(summary.read_text())


def _caption(kind: str, summary: dict, threshold: str = "") -> str:
    """Compose the caption from the numbers the figure was actually drawn from."""
    meta = FIGURE_META[kind]
    n_rows = summary.get("n_rows") or 0
    n_cols = summary.get("n_columns") or 0
    tool = summary.get("tool") or "MAFFT"
    detail = f"{meta['description']} ({n_rows} isoforms, {n_cols} alignment columns)"

    if kind == "wrapped_alignment":
        w = summary.get("wrapped") or {}
        detail += (f", {w.get('nBlocks')} blocks of {w.get('colsPerBlock')} columns "
                   f"on {w.get('nPages')} pages")
    if kind == "candidate_alignment_detail":
        cols = summary.get("candidate_columns") or []
        aa = summary.get("candidate_aa") or []
        affected = summary.get("affected") or []
        if aa and cols:
            detail += (f"; candidate {summary.get('candidate')} spans aa {aa[0]}–{aa[1]} "
                       f"(columns {cols[0]}–{cols[1]}) and is carried by "
                       f"{len(affected)} of {max(n_rows - 1, 0)} alternative isoforms")

    return build_caption(
        gene=summary.get("gene") or "",
        species=summary.get("species") or "",
        protein_id=summary.get("primary") or "",
        description=detail,
        coordinate_system=f"alignment columns, mapped to residue positions on "
                          f"{summary.get('primary') or 'the primary protein'}",
        annotation_source=f"{tool} multiple sequence alignment",
        threshold=threshold,
    )


def alignment_index(run_dir: Path) -> Path:
    return run_dir / "website_indices" / "isoform_alignment_index.json"


def generate(run_dir: Path, model_json: Path | None = None, *,
             index_json: Path | None = None) -> dict:
    """Render the isoform-analysis figures and register their Gallery cards.

    ``model_json`` is accepted so this stage has the same signature as the other
    figure stages the pipeline calls; these figures are built from the run's
    isoform-alignment index, which is located here when not given explicitly.
    """
    del model_json  # the alignment index is this stage's input
    index_json = index_json or alignment_index(run_dir)
    if not index_json.exists():
        return {"cards": 0, "figures": 0, "indices": 0,
                "warnings": [f"no isoform alignment index at {index_json}"]}

    fig_dir = run_dir / FIGURE_SUBDIR
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig_rel = FIGURE_SUBDIR.as_posix()
    run_id = run_dir.name

    # Every species with its own within-species alignment gets its own three
    # figures, rendered from that species' own alignment. Nothing is shared or
    # copied between species; a species without an alignment is reported with the
    # exact reason instead of being given another species' figures.
    index_doc = json.loads(index_json.read_text())
    entries = index_doc.get("species") or []
    renderable = [e for e in entries if e.get("status") == "available"]

    cards: list[dict] = []
    manifest: list[dict] = []
    warnings: list[str] = []
    unavailable: list[dict] = []
    for entry in entries:
        if entry.get("status") != "available":
            unavailable.append({
                "species_id": entry.get("species_id"),
                "reason": entry.get("status_reason")
                          or f"no within-species isoform alignment "
                             f"({entry.get('status') or 'unavailable'})",
            })
    for entry in renderable:
        sid = entry.get("species_id") or ""
        _render_species(index_json, fig_dir, fig_rel, run_id, sid,
                        cards, manifest, warnings,
                        source_alignment=entry.get("alignment_file") or "")

    indices = register_gallery_cards(run_dir, cards,
                                     drop_figure_ids=SUPERSEDED_FIGURE_IDS)
    return {"cards": len(cards), "figures": len(manifest), "indices": indices,
            "warnings": warnings, "unavailable": unavailable}


def _render_species(index_json: Path, fig_dir: Path, fig_rel: str,
                    run_id: str, species_filter: str,
                    cards: list, manifest: list, warnings: list, *,
                    source_alignment: str = "") -> None:
    """Render and register the three isoform figures for exactly one species."""
    summary = _render(index_json, fig_dir, species_filter)
    species_id = summary.get("species_id") or species_filter
    gene = summary.get("gene") or "GENE"
    species = summary.get("species") or species_id

    produced = {s["name"] for s in summary.get("render") or []}
    # The wrapped export is reported per page, so it is keyed by its own stem.
    if any(n.startswith("wrapped_alignment_p") for n in produced):
        produced.add("wrapped_alignment")

    for kind in ORDER:
        if kind not in produced:
            continue
        meta = FIGURE_META[kind]
        stem = f"main_{species_id}_{kind}"
        svg, pdf = fig_dir / f"{stem}.svg", fig_dir / f"{stem}.pdf"
        png, tsv = fig_dir / f"{stem}.png", fig_dir / f"{stem}.tsv"
        if not svg.exists():
            warnings.append(f"{stem}: the renderer produced no SVG")
            continue
        if not _rasterise_png(svg, png):
            warnings.append(f"could not rasterise {stem}.png at {EXPORT_DPI} dpi")

        caption = _caption(kind, summary)
        write_caption_file(fig_dir, stem, caption)
        card = figure_card(
            figure_id=stem,
            figure_type=kind,
            title=f"{gene} · {meta['title']}",
            category=CATEGORY,
            run_id=run_id,
            figure_dir=fig_rel,
            stem=stem,
            scientific_question=meta["question"],
            interpretation=meta["interpretation"],
            caption=caption,
            gene_symbol=gene,
            species=species,
            species_id=species_id,
            protein_id=summary.get("primary") or "",
            transcript_id=summary.get("transcript") or "",
            has_table=tsv.exists(),
            source_files=[source_alignment] if source_alignment else [],
        )
        # These figures come from the shared figure specification, not matplotlib.
        card["renderer"] = "shared_figure_specification"
        entry = next((s for s in summary.get("render") or []
                      if s["name"] == kind), None)
        if entry:
            card["width_pt"], card["height_pt"] = entry.get("width"), entry.get("height")
        if kind == "wrapped_alignment":
            card["n_pages"] = (summary.get("wrapped") or {}).get("nPages")
        cards.append(card)

        for ext, path in (("svg", svg), ("pdf", pdf), ("png", png), ("tsv", tsv)):
            if path.exists():
                manifest.append({
                    "group": GROUP, "kind": kind, "title": meta["title"], "format": ext,
                    "species_id": species_id, "protein_id": summary.get("primary"),
                    "path": _rel(path),
                })

    for entry in summary.get("render") or []:
        warnings.extend(f"{species_id}/{entry['name']}: {w}"
                        for w in entry.get("warnings") or [])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--index", type=Path, default=None,
                    help="isoform alignment index (defaults to the run's own)")
    args = ap.parse_args()
    run_dir = args.run_dir.resolve()
    index_json = args.index or alignment_index(run_dir)
    if not index_json.exists():
        raise SystemExit(f"alignment index not found: {index_json}")
    res = generate(run_dir, index_json=index_json)
    print(f"OK — {res['cards']} isoform-analysis card(s), {res['figures']} file(s), "
          f"{res['indices']} index file(s) updated for {run_dir.name}")
    for w in res["warnings"]:
        print(f"  warning: {w}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
