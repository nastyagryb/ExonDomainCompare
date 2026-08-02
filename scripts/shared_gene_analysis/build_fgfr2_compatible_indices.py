"""Orchestrator: write FGFR2-compatible website indices for shared gene runs."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List

from .common import SharedRunContext, write_json
from .public_paths import sanitize_public_payload, write_public_download_projections
from .indices.coordinate_track import build_coordinate_track_index
from .indices.msa import build_msa_index
from .indices.synteny_locus import build_synteny_locus_index

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
from framework.local_registry import discover_run  # noqa: E402
from framework.portable_config import load_config  # noqa: E402

OUTPUT_FILES = {
    "coordinate_track_index": build_coordinate_track_index,
    "msa_index": build_msa_index,
    "synteny_locus_index": build_synteny_locus_index,
}


def build_fgfr2_compatible_indices(run_dir: Path) -> Dict[str, Any]:
    """Build coordinate_track, msa and synteny_locus indices for one run."""
    ctx = SharedRunContext.from_run_dir(run_dir)
    out_dir = ctx.website_indices
    out_dir.mkdir(parents=True, exist_ok=True)
    written: List[str] = []
    payloads: Dict[str, Any] = {}

    for stem, builder in OUTPUT_FILES.items():
        if stem == "coordinate_track_index":
            data = builder(ctx)
        else:
            data = builder(ctx)
        path = out_dir / f"{stem}.json"
        data = sanitize_public_payload(data)
        write_json(path, data)
        written.append(stem)
        payloads[stem] = data

    write_public_download_projections(ctx.run_dir)
    return {"written": written, "indices": payloads}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True, help="Run directory name under runs/")
    args = ap.parse_args()
    config = load_config(repository_root=ROOT)
    run_dir = discover_run(config, args.run_id) or config.runs_root / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"Run not found: {run_dir}")
    result = build_fgfr2_compatible_indices(run_dir)
    print(f"OK fgfr2_compatible_indices  written={','.join(result['written'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
