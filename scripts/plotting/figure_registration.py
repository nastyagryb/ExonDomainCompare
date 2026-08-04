#!/usr/bin/env python3
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
from exondomaincompare.framework import production_contract
from exondomaincompare.contracts import file_sha256
from exondomaincompare.config import load_config

RUNTIME_CONFIG = load_config(repository_root=ROOT)
from exondomaincompare.runs.registry import RegistryError, resolve_run_record

INDEX_NAMES = ("figures_index.json", "generic/figures_index.json",
               "figure_index.json")

_PATH_KEYS = ("png_url", "svg_url", "pdf_url", "table_url", "thumbnail",
              "png_path", "svg_path", "pdf_path")

AVAILABLE = "available"

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
    out: set = set()
    live = {c.get("figure_id") for c in cards}
    for c in cards:
        for sid in c.get("supersedes") or []:
            if sid in live and sid != c.get("figure_id"):
                out.add(sid)
    return out


def _contract_for(run_dir: Path) -> Dict[str, Any]:
    config = {}
    try:
        config = json.loads((run_dir / "run_config.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        pass
    return production_contract.resolve(config.get("gene_symbol")).identity()


def normalise_index(doc: Dict[str, Any], run_dir: Path) -> Dict[str, Any]:
    run_id = run_dir.name
    cards = [c for c in (doc.get("figures") or []) if isinstance(c, dict)]
    superseded = _superseded(cards)
    identity = _contract_for(run_dir)

    kept: List[Dict[str, Any]] = []
    rejected: List[Dict[str, Any]] = []
    downgraded: List[Dict[str, Any]] = []
    seen: set = set()
    for card in cards:
        fid = card.get("figure_id") or card.get("id")
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

        if not card.get("scope"):
            card["scope"] = "species" if card.get("species_id") else "comparative"

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


def generate(run_dir: Path, model_json: Path) -> Dict[str, Any]:
    res = normalise_run(run_dir)
    total = sum(r.get("n_registered", 0) for r in res["indices"].values())
    return {"figures": 0, "registered": total, "registration": res}


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description='The one place a Gallery card becomes visible.',
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
