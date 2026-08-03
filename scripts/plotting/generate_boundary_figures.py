#!/usr/bin/env python3
"""Retire the superseded exon–domain-boundary figure cards of a run.

All five boundary figures this stage used to draw are now main figures produced by
``generate_shared_main_figures`` — boundary-on-architecture, signed boundary
distances and the boundary-class summary — or are on-demand exports of the Boundary
page rather than Gallery cards: the selected-boundary detail and the per-boundary
evidence supplement. Keeping a second matplotlib version of them would give the
reader two entry points to the same analysis that could silently disagree.

This stage therefore draws nothing. It removes the cards it used to own, by
``figure_id`` only, so no other generator's cards can be affected, and drops the
matching entries from the export manifest in the coordinate model.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from exondomaincompare.presentation import figure_captions as fc  # noqa: E402

LEGACY_GROUP = "exon_domain_boundaries"

# Cards this stage retires. The first three live on as shared main figures; the last
# two are exports the Boundary page produces on request for a selected boundary.
SUPERSEDED_FIGURE_IDS = (
    "boundary_{sp}_boundary_on_architecture",
    "boundary_{sp}_signed_boundary_distances",
    "boundary_{sp}_boundary_class_summary",
    "boundary_{sp}_selected_boundary_detail",
    "boundary_{sp}_boundary_evidence_supplement",
    "generic_exon_domain_boundary_distribution",
)


def generate(run_dir: Path, model_json: Path) -> dict:
    run_dir = Path(run_dir)
    index = json.loads(model_json.read_text())
    models = index.get("models") or []
    species_ids = [m.get("species_id") or "sp" for m in models] or ["sp"]
    drop: List[str] = []
    for template in SUPERSEDED_FIGURE_IDS:
        if "{sp}" in template:
            drop += [template.format(sp=sp) for sp in species_ids]
        else:
            drop.append(template)

    fc.remove_retired_figure_files(
        run_dir / "results" / "generic_gene_analysis" / "figures"
        / "exon_domain_boundaries", drop)

    kept: List[Dict[str, Any]] = [f for f in (index.get("publication_figures") or [])
                                  if f.get("group") != LEGACY_GROUP
                                  and f.get("kind") not in set(drop)]
    if len(kept) != len(index.get("publication_figures") or []):
        index["publication_figures"] = kept
        model_json.write_text(json.dumps(index, indent=2))
    fc.register_gallery_cards(run_dir, [], drop_figure_ids=drop)
    return {"figures": 0, "cards": 0, "retired": len(drop), "manifest": []}


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
    print(f"OK — retired {res['retired']} superseded boundary card id(s) for "
          f"{run_dir.name}; the boundary main figures are owned by the shared "
          f"renderer")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
