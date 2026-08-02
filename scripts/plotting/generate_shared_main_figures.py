#!/usr/bin/env python3
"""Generate the Figure Gallery's main single-species figures.

These figures are rendered by the *same* code the Gene Explorer exports through:
this stage shells out to ``scripts/plotting/render_main_figures.mjs``, which imports
the shared figure builders (``webapp/frontend/src/pages/viewers/mainFigures.js``)
through the shared adapter (``figureData.js``) and is handed the same coordinate
model objects the React views receive. A Gallery figure and the corresponding Gene
Explorer figure are therefore one figure with one implementation, rather than two
renderers that have to be kept in agreement by hand.

Figures produced per species (only what the data supports — nothing is fabricated):
  1. primary exon-to-protein projection      (also available pre-cluster)
  2. integrated domain architecture          (needs the InterPro layer)
  3. boundary-on-architecture                (needs classified boundaries)
  4. signed boundary distances               (       "                    )
  5. boundary-class summary                  (       "                    )

Formats: standalone SVG, true vector PDF, 300 dpi PNG, and a source TSV where the
figure has an underlying table. The PDFs contain real path and text operators — no
embedded raster and no page sized from pixel counts.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

RENDERER = ROOT / "scripts" / "plotting" / "render_main_figures.mjs"
FIGURE_SUBDIR = Path("results") / "generic_gene_analysis" / "figures" / "main"
GROUP = "main_figures"
EXPORT_DPI = 300

# One card per figure, in reading order. Each states the single scientific question
# the figure answers and a deliberately cautious interpretation.
FIGURE_META = {
    "primary_exon_projection": {
        "title": "Primary exon-to-protein projection",
        "category": "Exon structure",
        "question": "Which coding exons produce which regions of the selected primary "
                    "protein?",
        "interpretation": "Coding exons are projected onto the primary protein "
                          "sequence. Exploratory candidate intervals are shown for "
                          "orientation only and are not validated events.",
        "kind": "main",
        "source_table": "results/core_gene_analysis/exon_protein_map.tsv",
    },
    "integrated_domain_architecture": {
        "title": "Integrated domain architecture",
        "category": "Domain architecture",
        "question": "How do annotated domain instances, family assignments, membrane "
                    "topology and the coding-exon structure align along the primary "
                    "protein?",
        "interpretation": "Representative InterPro domain instances are resolved by "
                          "coordinate, so repeated domains sharing one accession stay "
                          "distinct. Family and superfamily assignments are shown in a "
                          "separate neutral row and are not structural domains.",
        "kind": "main",
        "source_table": "results/core_gene_analysis/domain_features.tsv",
    },
    "boundary_on_architecture": {
        "title": "Exon boundaries on the domain architecture",
        "category": "Exon–domain boundaries",
        "question": "Where do internal coding-exon boundaries fall relative to "
                    "annotated domain edges?",
        "interpretation": "Boundary markers are coloured by their mutually exclusive "
                          "relation class. Co-location of a boundary with a domain edge "
                          "is a positional observation, not evidence of a functional "
                          "relationship.",
        "kind": "main",
        "source_table": "results/core_gene_analysis/exon_domain_boundary_distances.tsv",
    },
    "signed_boundary_distances": {
        "title": "Signed boundary distances to the nearest domain edge",
        "category": "Exon–domain boundaries",
        "question": "How far, and on which side, does each internal coding-exon "
                    "boundary sit relative to the nearest domain edge?",
        "interpretation": "Distances are signed against the specific domain instance "
                          "used in the calculation, so the direction of the offset and "
                          "the domain identity are both preserved.",
        "kind": "main",
        "source_table": "results/core_gene_analysis/exon_domain_boundary_distances.tsv",
    },
    "boundary_class_summary": {
        "title": "Boundary-class summary",
        "category": "Exon–domain boundaries",
        "question": "How are internal coding-exon boundaries distributed across the "
                    "domain-relation classes?",
        "interpretation": "Counts of the mutually exclusive classes for one protein. "
                          "Class membership depends on the near-edge threshold stated "
                          "in the figure.",
        "kind": "main",
        "source_table": "results/core_gene_analysis/exon_domain_boundary_distances.tsv",
    },
}

# Cards superseded by the figures above, or scientifically redundant with them.
# They are removed from the Gallery rather than left as duplicate entry points; the
# relationships they showed are all present in the integrated main figures.
SUPERSEDED_FIGURE_IDS = {
    # replaced by the shared main figures
    "generic_domain_architecture",
    "generic_exon_domain_boundary_distribution",
    "primary_protein_exon_projection",
    "domain_arch_{sp}_representative_architecture",
    "boundary_{sp}_boundary_on_architecture",
    "boundary_{sp}_signed_boundary_distances",
    "boundary_{sp}_boundary_class_summary",
    # pairwise overlays already contained in the integrated architecture figure
    "domain_arch_{sp}_domain_exon_projection",
    "domain_arch_{sp}_domain_boundary_overlay",
    "domain_arch_{sp}_domain_candidate_overlay",
    # a table rendered as a plot, and an on-demand export that was pinned as a card
    "boundary_{sp}_boundary_evidence_supplement",
    "boundary_{sp}_selected_boundary_detail",
}


CLASS_CAPTION_NAMES = {
    "near_domain_edge": "near a domain edge",
    "inside_domain": "inside a domain",
    "outside_annotated_domains": "outside annotated domains",
    "exact_domain_edge": "exactly at a domain edge",
    "uncertain": "unresolved",
}


def _class_counts(model: dict) -> dict[str, int]:
    counts: dict[str, int] = {}
    for b in model.get("exon_boundaries") or []:
        key = b.get("boundary_class") or b.get("classification") or "uncertain"
        counts[key] = counts.get(key, 0) + 1
    return counts


def _caption(kind: str, model: dict, gene: str, species: str) -> str:
    """Compose the figure caption from the figure's own data.

    Captions are derived rather than written per run, so the numbers a reader sees
    underneath a figure are the numbers the figure was drawn from.
    """
    protein = model.get("protein_id") or "the primary protein"
    transcript = model.get("transcript_id")
    length = model.get("protein_length")
    n_exons = len(model.get("exons") or [])
    n_cand = len(model.get("candidate_regions") or [])
    n_dom = len(model.get("representative_domains") or [])
    n_tm = len(model.get("tm_regions") or [])
    counts = _class_counts(model)
    n_bnd = sum(counts.values())
    threshold = model.get("near_edge_threshold_aa")
    lead = f"{gene}, {species}."

    def distribution() -> str:
        parts = [f"{n} {CLASS_CAPTION_NAMES.get(k, k)}"
                 for k, n in sorted(counts.items(), key=lambda kv: -kv[1])]
        return "; ".join(parts)

    if kind == "primary_exon_projection":
        src = f" of transcript {transcript}" if transcript else ""
        tail = (f" {n_cand} exploratory candidate interval"
                f"{'s' if n_cand != 1 else ''} are shown for orientation."
                if n_cand else "")
        return (f"{lead} The {n_exons} coding exons{src} projected onto primary protein "
                f"{protein} ({length} aa).{tail}")

    if kind == "integrated_domain_architecture":
        tm = (f", {n_tm} predicted transmembrane helix"
              f"{'es' if n_tm != 1 else ''}" if n_tm else "")
        return (f"{lead} {n_dom} representative InterPro domain instance"
                f"{'s' if n_dom != 1 else ''}{tm} and {n_exons} coding exons along "
                f"primary protein {protein} ({length} aa). Repeated domains sharing one "
                f"accession are kept apart by their coordinates.")

    if kind == "boundary_on_architecture":
        return (f"{lead} The {n_bnd} internal coding-exon boundaries of {protein} "
                f"({length} aa), coloured by their relation to the nearest "
                f"representative domain edge: {distribution()}.")

    if kind == "signed_boundary_distances":
        band = f" The shaded band marks ±{threshold} aa." if threshold else ""
        return (f"{lead} Signed distance from each of the {n_bnd} internal coding-exon "
                f"boundaries of {protein} to the nearest edge of the domain instance "
                f"used in the calculation. Negative values lie before the edge, positive "
                f"values after it.{band}")

    if kind == "boundary_class_summary":
        thr = (f" Classes are mutually exclusive at a near-edge threshold of "
               f"{threshold} aa." if threshold else "")
        return (f"{lead} Distribution of the {n_bnd} internal coding-exon boundaries of "
                f"{protein} across the domain-relation classes: {distribution()}.{thr}")

    return f"{lead} {protein}."


def _rel(p: Path) -> str:
    try:
        return str(p.resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


def _api(run_id: str, rel_to_run: str, inline: bool = False) -> str:
    q = f"path={rel_to_run}" + ("&inline=true" if inline else "")
    return f"/api/runs/{run_id}/files?{q}"


def _stamp_dpi(png: Path, dpi: int) -> None:
    """Record the physical resolution in the PNG's pHYs chunk.

    Neither rsvg-convert nor sips writes it, so a downstream tool would otherwise
    have to guess the figure's physical size from the pixel count alone.
    """
    try:
        from PIL import Image

        with Image.open(png) as im:
            im.save(png, dpi=(dpi, dpi))
    except Exception:
        pass


def _rasterise_png(svg: Path, png: Path, dpi: int = EXPORT_DPI) -> bool:
    """Rasterise a figure SVG at publication resolution.

    The figure is specified in points, so the pixel size follows from the physical
    geometry and the requested resolution — never from an arbitrary zoom factor.
    """
    zoom = dpi / 72.0
    if shutil.which("rsvg-convert"):
        result = subprocess.run(
            ["rsvg-convert", "--dpi-x", str(dpi), "--dpi-y", str(dpi),
             "--zoom", f"{zoom:.4f}", str(svg), "-o", str(png)],
            capture_output=True,
        )
        if result.returncode == 0 and png.exists():
            _stamp_dpi(png, dpi)
            return True
    try:  # optional dependency; used when the CLI tool is unavailable
        import cairosvg  # type: ignore

        cairosvg.svg2png(url=str(svg), write_to=str(png), scale=zoom)
        if png.exists():
            _stamp_dpi(png, dpi)
            return True
        return False
    except Exception:
        pass
    # Last resort on macOS: rasterise the vector PDF that accompanies the SVG.
    pdf = svg.with_suffix(".pdf")
    if pdf.exists() and shutil.which("sips"):
        result = subprocess.run(
            ["sips", "-s", "format", "png", str(pdf), "--out", str(png)],
            capture_output=True,
        )
        if result.returncode == 0 and png.exists():
            _stamp_dpi(png, dpi)
            return True
    return False


def _render(model_json: Path, out_dir: Path) -> list[dict]:
    if shutil.which("node") is None:
        raise SystemExit("node is required to render the shared main figures")
    result = subprocess.run(
        ["node", str(RENDERER), str(model_json), str(out_dir)],
        capture_output=True, text=True, cwd=ROOT,
    )
    if result.returncode != 0:
        raise SystemExit(f"figure renderer failed:\n{result.stdout}\n{result.stderr}")
    print(result.stdout.rstrip())
    summary = out_dir / "render_summary.json"
    return json.loads(summary.read_text()) if summary.exists() else []


def _register_gallery(run_dir: Path, cards: list[dict], species_ids: list[str]) -> int:
    """Insert the main cards and drop every superseded card."""
    drop = set()
    for template in SUPERSEDED_FIGURE_IDS:
        if "{sp}" in template:
            drop.update(template.format(sp=sp) for sp in species_ids)
        else:
            drop.add(template)
    new_ids = {c["figure_id"] for c in cards}

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
                if f.get("figure_id") not in drop and f.get("figure_id") not in new_ids]
        doc["figures"] = cards + kept
        avail = doc.get("available")
        if isinstance(avail, list):
            def _id(a):
                return a.get("figure_id") if isinstance(a, dict) else a
            doc["available"] = [a for a in avail if _id(a) not in drop and _id(a) not in new_ids] \
                + [c["figure_id"] for c in cards]
        fp.write_text(json.dumps(doc, indent=2))
        written += 1
    return written


def generate(run_dir: Path, model_json: Path) -> dict:
    index = json.loads(model_json.read_text())
    models = index.get("models") or []
    if not models:
        raise SystemExit("no models in the coordinate model")
    gene = index.get("gene_symbol") or "GENE"
    run_id = run_dir.name

    fig_dir = run_dir / FIGURE_SUBDIR
    fig_dir.mkdir(parents=True, exist_ok=True)
    rendered = _render(model_json, fig_dir)

    by_species = {m.get("species_id"): m for m in models}
    cards: list[dict] = []
    manifest: list[dict] = []
    warnings: list[str] = []

    for entry in rendered:
        stem = entry["stem"]
        # stems are main_<species_id>_<kind>
        body = stem[len("main_"):]
        kind = next((k for k in FIGURE_META if body.endswith(k)), None)
        if kind is None:
            warnings.append(f"unrecognised figure stem: {stem}")
            continue
        sp = body[: -(len(kind) + 1)]
        model = by_species.get(sp) or {}
        meta = FIGURE_META[kind]

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
                    "group": GROUP, "kind": kind, "title": meta["title"], "format": ext,
                    "species_id": sp, "protein_id": model.get("protein_id"),
                    "path": _rel(path),
                })

        sci = model.get("scientific_name") or sp
        cards.append({
            "figure_id": stem,
            "figure_type": kind,
            "title": f"{gene} · {meta['title']}",
            "category": meta["category"],
            "section": meta["category"],
            "kind": meta["kind"],
            "scientific_question": meta["question"],
            "interpretation": meta["interpretation"],
            "caption": _caption(kind, model, gene, sci),
            "gene_symbol": gene,
            "species": sci,
            "species_id": sp,
            "protein_id": model.get("protein_id"),
            "transcript_id": model.get("transcript_id"),
            "stage": "post_cluster" if model.get("status") == "available" else "pre_cluster",
            "status": "available",
            "renderer": "shared_figure_specification",
            "width_pt": entry.get("width"),
            "height_pt": entry.get("height"),
            "png_url": _api(run_id, rel["png"], inline=True) if "png" in rel else "",
            "svg_url": _api(run_id, rel["svg"]) if "svg" in rel else "",
            "pdf_url": _api(run_id, rel["pdf"]) if "pdf" in rel else "",
            "table_url": _api(run_id, rel["tsv"]) if "tsv" in rel else "",
            "source_table": meta["source_table"],
            "source_files": [meta["source_table"]],
            "error": "",
        })
        if entry.get("warnings"):
            warnings.extend(f"{stem}: {w}" for w in entry["warnings"])

    # Order the cards the way the Gallery reads them.
    order = list(FIGURE_META)
    cards.sort(key=lambda c: (c["species_id"] or "",
                              order.index(next(k for k in order
                                               if c["figure_id"].endswith(k)))))

    existing = [f for f in (index.get("publication_figures") or [])
                if f.get("group") != GROUP]
    index["publication_figures"] = existing + manifest
    model_json.write_text(json.dumps(index, indent=2))

    indices = _register_gallery(run_dir, cards, [m.get("species_id") for m in models])
    return {"figures": len(manifest), "cards": len(cards), "indices": indices,
            "warnings": warnings}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--model", type=Path, default=None)
    args = ap.parse_args()
    run_dir = args.run_dir.resolve()
    model_json = args.model or (
        run_dir / "website_indices" / "generic" / "protein_coordinate_model.json")
    if not model_json.exists():
        raise SystemExit(f"coordinate model not found: {model_json}")
    res = generate(run_dir, model_json)
    print(f"OK — {res['cards']} main figure card(s), {res['figures']} file(s), "
          f"{res['indices']} index file(s) updated for {run_dir.name}")
    for w in res["warnings"]:
        print(f"  warning: {w}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
