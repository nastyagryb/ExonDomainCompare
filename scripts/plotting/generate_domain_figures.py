#!/usr/bin/env python3
"""Domain-architecture supplement figure for the single-species Gallery.

The integrated domain architecture itself — representative domains, families,
transmembrane topology, coding exons, boundaries and candidate overlays on one
protein axis — is a main figure owned by ``generate_shared_main_figures``, which
the Gene Explorer exports through as well. This stage therefore produces only what
that integrated figure deliberately leaves out:

  C  member-database signature supplement   every signature behind the
                                            representative domains, grouped by its
                                            source database

The supplement is exported as SVG, a true-vector PDF, a 300 dpi PNG and its source
table, and registers exactly one Gallery card marked as a supplement. The cards of
the overlay figures it replaced are retired by id.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from plotting import figure_captions as fc  # noqa: E402
from plotting import shared_gene_plots as sgp  # noqa: E402

FIGURE_DIR = "results/generic_gene_analysis/figures/domain_architecture"
CATEGORY = "Domain architecture"
DOMAIN_SOURCE = "results/core_gene_analysis/domain_features.tsv"
SIGNATURE_SOURCE = "results/core_gene_analysis/interpro_annotations.tsv"

# Overlay cards this stage no longer produces: each relationship they showed is
# present in the integrated main domain-architecture figure. Matched by id only.
SUPERSEDED_FIGURE_IDS = (
    "domain_arch_{sp}_representative_architecture",
    "domain_arch_{sp}_domain_exon_projection",
    "domain_arch_{sp}_domain_boundary_overlay",
    "domain_arch_{sp}_domain_candidate_overlay",
    "generic_domain_architecture",
)


def _rel(p: Path) -> str:
    try:
        return str(p.resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


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

    manifest: List[Dict[str, Any]] = []
    cards: List[Dict[str, Any]] = []
    drop: List[str] = []

    for model in models:
        sp = model.get("species_id") or "sp"
        drop += [t.format(sp=sp) for t in SUPERSEDED_FIGURE_IDS]
        fc.remove_retired_figure_files(
            fig_dir, [t.format(sp=sp) for t in SUPERSEDED_FIGURE_IDS
                      if t.startswith("domain_arch_")])
        if model.get("status") != "available":
            continue  # pre-cluster: no fabricated domain figures
        species = model.get("scientific_name") or sp
        protein_id = model.get("protein_id") or ""
        protein_length = sgp._int(model.get("protein_length"))
        signatures = model.get("member_signatures") or []
        domains = model.get("representative_domains") or []
        if not signatures or not protein_length:
            continue

        fig_dir.mkdir(parents=True, exist_ok=True)
        stem = f"domain_arch_{sp}_member_signature_supplement"
        table = sgp.member_signature_table(signatures)
        n_db = len({r["member_database"] for r in table})
        n_integrated = sum(1 for r in table
                           if r["integration_status"] == "integrated")
        drawn = sgp.plot_member_signature_supplement(
            fig_dir, stem, gene_symbol=gene, species_name=species,
            protein_id=protein_id, protein_length=protein_length,
            signatures=signatures,
            footnote=f"Source: InterProScan member-database matches · "
                     f"{len(domains)} representative domain instance(s) are derived "
                     f"from these signatures · run {index.get('run_id', run_id)}")
        if not drawn:
            continue
        sgp.write_source_table(fig_dir, stem, sgp.MEMBER_SIGNATURE_TABLE_COLUMNS, table)
        caption = fc.build_caption(
            gene=gene, species=species, protein_id=protein_id,
            description=f"all {len(table)} member-database signature(s) from {n_db} "
                        f"database(s) behind the representative domain annotation, "
                        f"{n_integrated} of them integrated into an InterPro entry",
            annotation_source="InterProScan member-database matches",
            status=fc.DESCRIPTIVE_STATUS)
        fc.write_caption_file(fig_dir, stem, caption)
        cards.append(fc.figure_card(
            figure_id=stem, figure_type="member_signature_supplement",
            title=f"{gene} · Member-database signature supplement",
            category=CATEGORY, kind="supplement", run_id=run_id,
            figure_dir=FIGURE_DIR, stem=stem,
            scientific_question="Which member-database signatures underlie the "
                                "representative domain annotation of this protein?",
            interpretation="Every signature keeps its own row, grouped by its source "
                           "database. Signatures already integrated into an InterPro "
                           "entry are filled; unintegrated signatures are drawn open "
                           "because they are not annotated domains. Overlapping "
                           "signatures are not independent observations.",
            caption=caption, gene_symbol=gene, species=species, species_id=sp,
            protein_id=protein_id, transcript_id=model.get("transcript_id") or "",
            source_files=[SIGNATURE_SOURCE, DOMAIN_SOURCE]))
        for ext in ("svg", "pdf", "png", "tsv"):
            path = fig_dir / f"{stem}.{ext}"
            if path.exists():
                manifest.append({"group": CATEGORY, "kind": stem,
                                 "title": "Member-database signature supplement",
                                 "format": ext, "species_id": sp,
                                 "scientific_name": species,
                                 "protein_id": protein_id, "path": _rel(path)})

    own_kinds = {m["kind"] for m in manifest}
    other = [f for f in (index.get("publication_figures") or [])
             if f.get("group") != "domain_architecture" and f.get("kind") not in own_kinds]
    index["publication_figures"] = other + manifest
    model_json.write_text(json.dumps(index, indent=2))
    fc.register_gallery_cards(run_dir, cards, drop_figure_ids=drop)
    return {"figures": len(manifest), "cards": len(cards), "manifest": manifest}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--model", type=Path, default=None)
    args = ap.parse_args()
    run_dir = args.run_dir.resolve()
    model_json = args.model or (run_dir / "website_indices" / "generic"
                                / "protein_coordinate_model.json")
    if not model_json.exists():
        raise SystemExit(f"coordinate model not found: {model_json}")
    res = generate(run_dir, model_json)
    print(f"OK — wrote {res['figures']} domain figure file(s) in {res['cards']} "
          f"Gallery card(s) for {run_dir.name}")
    for f in res["manifest"]:
        if f["format"] == "svg":
            print(f"  · {f['kind']}: {f['path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
