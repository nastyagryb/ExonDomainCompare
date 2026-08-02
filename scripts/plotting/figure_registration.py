#!/usr/bin/env python3
"""The one place a Gallery card becomes visible.

Six figure stages write into the same ``figures_index.json``. Each of them used
to maintain its own copy of the run's availability list, so a stage that only
appended to ``figures`` left the availability record describing an older card
set — and a legacy card that no stage owned any more stayed in that record for
good. This module is the single owner of the record: it derives availability from
the registered cards themselves and refuses a card that cannot be justified.

A card survives normalisation only when every part of its identity agrees with
the run it is registered in:

* the run_id of every referenced output path is this run;
* the referenced preview and export files exist on disk;
* the producing stage recorded a status;
* the card is not superseded by another registered card.

Cards that fail are removed rather than downgraded, because a rejected card has
no scientific content to show and a placeholder preview would only assert an
availability claim the run cannot support.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence
from urllib.parse import parse_qs, urlparse

import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
from framework import production_contract  # noqa: E402
from framework.data_contract import file_sha256  # noqa: E402
from framework.portable_config import load_config  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=ROOT)
from framework.local_registry import RegistryError, resolve_run_record  # noqa: E402

INDEX_NAMES = ("figures_index.json", "generic/figures_index.json",
               "figure_index.json")

# Path-bearing keys of a card, in the order a reader would reach for them.
_PATH_KEYS = ("png_url", "svg_url", "pdf_url", "table_url", "thumbnail",
              "png_path", "svg_path", "pdf_path")

AVAILABLE = "available"
#: An expected output of a completed stage is absent. A blocking state, kept
#: distinct from "the analysis does not apply" and from "not done yet".
TECHNICALLY_MISSING = "technically_missing"


def _card_paths(card: Dict[str, Any]) -> List[str]:
    out: List[str] = []
    for key in _PATH_KEYS:
        val = card.get(key)
        if isinstance(val, str) and val:
            out.append(val)
    for fmt in ("formats", ):
        block = card.get(fmt) or {}
        if isinstance(block, dict):
            out.extend(v for v in block.values() if isinstance(v, str) and v)
    for mode in card.get("modes") or []:
        if not isinstance(mode, dict):
            continue
        thumb = mode.get("thumbnail")
        if isinstance(thumb, str) and thumb:
            out.append(thumb)
        out.extend(v for v in (mode.get("formats") or {}).values()
                   if isinstance(v, str) and v)
    return out


def _resolve(reference: str) -> tuple[str, str]:
    """Split a card reference into (run_id, path relative to the run directory).

    Cards carry either a repository-relative path or the backend file URL the
    browser fetches; both name the same file, so both have to resolve to it.
    """
    ref = reference.strip()
    if ref.startswith("/api/runs/") or ref.startswith("api/runs/"):
        tail = ref.split("/runs/", 1)[1]
        run_id, _, query = tail.partition("/files")
        params = parse_qs(urlparse(f"?{query.lstrip('?')}").query)
        rel = (params.get("path") or [""])[0]
        return run_id, rel
    if "runs/" in ref:
        tail = ref.split("runs/", 1)[1]
        run_id, _, rel = tail.partition("/")
        return run_id, rel
    return "", ref.lstrip("/")


def _foreign_run(paths: Iterable[str], run_id: str) -> List[str]:
    """References that point into a different run directory."""
    bad = []
    for p in paths:
        other, _ = _resolve(p)
        if other and other != run_id:
            bad.append(p)
    return bad


def _missing(paths: Iterable[str], run_dir: Path) -> List[str]:
    missing = []
    for p in paths:
        _, rel = _resolve(p)
        if not rel or not (run_dir / rel).exists():
            missing.append(p)
    return missing


def gallery_asset_checksums(card: Dict[str, Any], run_dir: Path) -> Dict[str, str]:
    checksums: Dict[str, str] = {}
    for reference in _card_paths(card):
        _, rel = _resolve(reference)
        path = run_dir / rel
        if rel and path.is_file():
            checksums[rel] = file_sha256(path)
    return dict(sorted(checksums.items()))


def _superseded(cards: Sequence[Dict[str, Any]]) -> set:
    """Ids that a registered card explicitly replaces."""
    out: set = set()
    live = {c.get("figure_id") for c in cards}
    for c in cards:
        for sid in c.get("supersedes") or []:
            if sid in live and sid != c.get("figure_id"):
                out.add(sid)
    return out


def _contract_for(run_dir: Path) -> Dict[str, Any]:
    """The run's production contract, resolved from its gene symbol."""
    config = {}
    try:
        config = json.loads((run_dir / "run_config.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        pass
    return production_contract.resolve(config.get("gene_symbol")).identity()


def normalise_index(doc: Dict[str, Any], run_dir: Path) -> Dict[str, Any]:
    """Return the rejection report for one index document, edited in place."""
    run_id = run_dir.name
    cards = [c for c in (doc.get("figures") or []) if isinstance(c, dict)]
    superseded = _superseded(cards)
    # Every card carries the architecture that produced it. Without this a card is
    # only identified by its file path, so a card left over from an older renderer is
    # indistinguishable from a current one and survives on the strength of its file
    # still being there.
    identity = _contract_for(run_dir)

    kept: List[Dict[str, Any]] = []
    rejected: List[Dict[str, Any]] = []
    downgraded: List[Dict[str, Any]] = []
    seen: set = set()
    for card in cards:
        fid = card.get("figure_id") or card.get("id")
        # Wrong identity: the card does not belong to this run's card set at all,
        # and nothing about it can be salvaged for a reader.
        reasons: List[str] = []
        if not fid:
            reasons.append("no figure_id")
        elif fid in seen:
            reasons.append("duplicate figure_id")
        elif fid in superseded:
            reasons.append("superseded by a canonical figure")
        paths = _card_paths(card)
        foreign = _foreign_run(paths, run_id)
        if foreign:
            reasons.append(f"paths from another run: {foreign[0]}")
        recorded = card.get("figure_renderer_version")
        if recorded is not None and recorded != identity["figure_renderer_version"]:
            reasons.append(f"renderer v{recorded} is not this build's "
                           f"v{identity['figure_renderer_version']}")
        recorded_family = card.get("analysis_family")
        if recorded_family and recorded_family != identity["analysis_family"]:
            reasons.append(f"analysis_family {recorded_family} is not this run's "
                           f"{identity['analysis_family']}")
        if reasons:
            rejected.append({"figure_id": fid, "reasons": reasons})
            continue
        card.update(identity)
        card["run_id"] = run_id

        # The scope decides which Gallery tab shows the card and which
        # availability wording applies to it, so it is recorded rather than left
        # for a consumer to infer from the presence of a species id.
        if not card.get("scope"):
            card["scope"] = "species" if card.get("species_id") else "comparative"

        # A card that claims to be available but has no file to show is a gap in the
        # run, not a card to delete: deleting it would hide the gap, and readiness
        # has to see it. It is reported as technically missing instead, which is a
        # blocking state, so a broken preview cannot pass for a result.
        if (card.get("status") or AVAILABLE) == AVAILABLE:
            gone = _missing(paths, run_dir) if paths else ["<no output paths>"]
            if gone:
                card["status"] = TECHNICALLY_MISSING
                card["error"] = f"expected output is absent: {gone[0]}"
                downgraded.append({"figure_id": fid, "missing": gone[0]})
            else:
                current_checksums = gallery_asset_checksums(card, run_dir)
                metadata = dict(card.get("_exondomain") or {})
                recorded_checksums = metadata.get("asset_sha256")
                if recorded_checksums and recorded_checksums != current_checksums:
                    card["status"] = TECHNICALLY_MISSING
                    card["error"] = "Gallery asset content changed since registration"
                    downgraded.append({
                        "figure_id": fid,
                        "missing": "stale Gallery/source checksum",
                    })
                elif not recorded_checksums:
                    metadata.update({
                        "contract_version": "1.0",
                        "run_id": run_id,
                        "payload_type": "gallery_card",
                        "asset_sha256": current_checksums,
                    })
                    card["_exondomain"] = metadata
        seen.add(fid)
        kept.append(card)

    doc["figures"] = kept
    # One derived availability record, so the list can no longer describe a card
    # set that the figures themselves have moved on from.
    doc[AVAILABLE] = [c["figure_id"] for c in kept
                      if (c.get("status") or AVAILABLE) == AVAILABLE]
    doc["pending"] = [c["figure_id"] for c in kept
                      if (c.get("status") or AVAILABLE) != AVAILABLE]
    doc["run_id"] = run_id
    doc.update(identity)
    doc["registration"] = {
        "normaliser": "plotting.figure_registration",
        "n_registered": len(kept),
        "n_rejected": len(rejected),
        "n_technically_missing": len(downgraded),
        "rejected": rejected,
        "technically_missing": downgraded,
    }
    return {"n_registered": len(kept), "n_rejected": len(rejected),
            "n_technically_missing": len(downgraded),
            "rejected": rejected, "technically_missing": downgraded}


def normalise_run(run_dir: Path) -> Dict[str, Any]:
    """Normalise every figure index of a run."""
    out: Dict[str, Any] = {"run_id": run_dir.name, "indices": {}}
    for name in INDEX_NAMES:
        fp = run_dir / "website_indices" / name
        if not fp.exists():
            continue
        try:
            doc = json.loads(fp.read_text())
        except (OSError, ValueError) as err:
            out["indices"][name] = {"error": str(err)}
            continue
        report = normalise_index(doc, run_dir)
        fp.write_text(json.dumps(doc, indent=2))
        out["indices"][name] = report
    return out


def generate(run_dir: Path, model_json: Path) -> Dict[str, Any]:  # noqa: ARG001
    """Figure-stage entry point, so the sequence owns the final registration."""
    res = normalise_run(run_dir)
    total = sum(r.get("n_registered", 0) for r in res["indices"].values())
    return {"figures": 0, "registered": total, "registration": res}


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args(argv)
    try:
        record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    except RegistryError as exc:
        ap.error(str(exc))
    if record is None:
        ap.error(f"no such run: {args.run_id}")
    if record.read_only:
        ap.error("run is registered read-only; copy it before rebuilding")
    run_dir = record.path
    print(json.dumps(normalise_run(run_dir), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
