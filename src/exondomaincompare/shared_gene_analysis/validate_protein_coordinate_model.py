from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import model_roles
from .protein_coordinate_model import build_models_for_run

BOUNDARY_MATCH_TOLERANCE_AA = 0  # boundary must land exactly on an exon-projection edge


class CoordinateModelError(ValueError):
    pass


_TRACKS = (
    "exons",
    "exon_boundaries",
    "representative_domains",
    "families_superfamilies",
    "member_signatures",
    "functional_sites",
    "disorder_regions",
    "tm_regions",
    "candidate_regions",
)


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _int(v: Any) -> Optional[int]:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def validate_model(model: Dict[str, Any], *, core_dir: Optional[Path] = None) -> List[str]:
    errors: List[str] = []
    sp = model.get("species_id", "?")
    pid = model.get("protein_id", "?")
    tag = f"[{sp}/{pid}]"
    length = _int(model.get("protein_length")) or 0

    if model.get("coordinate_system") != "protein_1_based_inclusive":
        errors.append(f"{tag} coordinate_system must be protein_1_based_inclusive")
    if length <= 0:
        errors.append(f"{tag} protein_length must be positive (got {model.get('protein_length')})")
    if not pid or pid == "?":
        errors.append(f"{tag} missing selected primary protein_id")

    # Rules 1 + 2: bounds and ordering for every feature on every track
    for track in _TRACKS:
        for f in model.get(track, []) or []:
            s, e = _int(f.get("start")), _int(f.get("end"))
            fid = f.get("id", "?")
            if s is None or e is None:
                errors.append(f"{tag} {track}:{fid} has non-numeric start/end")
                continue
            if s > e:
                errors.append(f"{tag} {track}:{fid} start {s} > end {e}")
            if length and (s < 1 or e > length):
                errors.append(f"{tag} {track}:{fid} coordinates [{s},{e}] outside [1,{length}]")
            for req in ("id", "label", "source", "source_file"):
                if not f.get(req):
                    errors.append(f"{tag} {track}:{fid} missing required field '{req}'")
            if "status" not in f and "confidence" not in f:
                errors.append(f"{tag} {track}:{fid} missing confidence/status")

    # Rule 4: exon-boundary positions agree with exon projections
    exon_edges = set()
    for ex in model.get("exons", []) or []:
        s, e = _int(ex.get("start")), _int(ex.get("end"))
        if s is not None:
            exon_edges.add(s - 1)  # boundary sits between exons
            exon_edges.add(s)
        if e is not None:
            exon_edges.add(e)
    for b in model.get("exon_boundaries", []) or []:
        pos = _int(b.get("start"))
        if pos is None:
            continue
        if not any(abs(pos - edge) <= BOUNDARY_MATCH_TOLERANCE_AA for edge in exon_edges):
            errors.append(
                f"{tag} boundary {b.get('id')} at aa {pos} does not coincide with any exon projection edge"
            )

    # Rule 8: domain-instance identity + boundary/instance consistency
    errors += domain_instance_errors(model)
    errors += boundary_instance_errors(model)

    # Rules 5/6/7 cross-check against the normalized source tables when available
    if core_dir is not None:
        errors += _cross_check_sources(model, core_dir, tag)

    return errors


def domain_instance_errors(model: Dict[str, Any]) -> List[str]:
    tag = f"[{model.get('species_id', '?')}/{model.get('protein_id', '?')}]"
    errs: List[str] = []
    seen: Dict[str, int] = {}
    per_acc: Dict[str, List[Dict[str, Any]]] = {}
    for d in model.get("representative_domains", []) or []:
        iid = d.get("domain_instance_id")
        if not iid:
            errs.append(f"{tag} representative domain {d.get('id')} has no domain_instance_id")
            continue
        expected = f"{d.get('interpro_accession') or 'NA'}:{d.get('start')}-{d.get('end')}"
        if iid != expected:
            errs.append(f"{tag} domain_instance_id {iid} does not match coordinates ({expected})")
        seen[iid] = seen.get(iid, 0) + 1
        per_acc.setdefault(d.get("interpro_accession") or "NA", []).append(d)
    for iid, n in seen.items():
        if n != 1:
            errs.append(f"{tag} domain_instance_id {iid} used by {n} domains (must be unique)")
    for acc, items in per_acc.items():
        by_start = sorted(items, key=lambda d: (_int(d.get("start")) or 0, _int(d.get("end")) or 0))
        for n, d in enumerate(by_start, start=1):
            if d.get("instance_number") != n:
                errs.append(f"{tag} {acc} instance at aa {d.get('start')}–{d.get('end')} has "
                            f"instance_number {d.get('instance_number')}, expected {n} "
                            "(numbering must follow sorted start coordinates)")
    return errs


def boundary_instance_errors(model: Dict[str, Any]) -> List[str]:
    tag = f"[{model.get('species_id', '?')}/{model.get('protein_id', '?')}]"
    thr = _int(model.get("near_edge_threshold_aa")) or 5
    by_instance = {d.get("domain_instance_id"): d
                   for d in model.get("representative_domains", []) or []}
    errs: List[str] = []
    for b in model.get("exon_boundaries", []) or []:
        iid = b.get("nearest_domain_instance_id")
        pos = _int(b.get("protein_position"))
        signed = b.get("signed_distance")
        if iid is None:
            if signed is not None:
                errs.append(f"{tag} boundary {b.get('id')} has a signed distance but no "
                            "nearest_domain_instance_id")
            continue
        dom = by_instance.get(iid)
        if dom is None:
            errs.append(f"{tag} boundary {b.get('id')} references unknown domain instance {iid}")
            continue
        ds, de = _int(dom.get("start")), _int(dom.get("end"))
        if (_int(b.get("nearest_domain_start")), _int(b.get("nearest_domain_end"))) != (ds, de):
            errs.append(f"{tag} boundary {b.get('id')} stores domain span "
                        f"{b.get('nearest_domain_start')}–{b.get('nearest_domain_end')} but "
                        f"instance {iid} spans {ds}–{de}")
        if b.get("nearest_domain_accession") != dom.get("interpro_accession"):
            errs.append(f"{tag} boundary {b.get('id')} accession "
                        f"{b.get('nearest_domain_accession')} disagrees with instance {iid}")
        edge_type = b.get("nearest_edge_type")
        expected_edge = ds if edge_type == "start" else de if edge_type == "end" else None
        if expected_edge is None:
            errs.append(f"{tag} boundary {b.get('id')} has no start/end nearest_edge_type")
            continue
        if _int(b.get("nearest_edge_position")) != expected_edge:
            errs.append(f"{tag} boundary {b.get('id')} nearest_edge_position "
                        f"{b.get('nearest_edge_position')} is not the {edge_type} of {iid} "
                        f"({expected_edge})")
        if pos is not None and signed is not None and signed != pos - expected_edge:
            errs.append(f"{tag} boundary {b.get('id')} signed_distance {signed} != "
                        f"{pos} - {expected_edge} ({pos - expected_edge})")
        if signed is not None and b.get("absolute_distance") != abs(signed):
            errs.append(f"{tag} boundary {b.get('id')} absolute_distance "
                        f"{b.get('absolute_distance')} != |{signed}|")
        cls = b.get("boundary_class") or b.get("class")
        if cls == "inside_domain" and pos is not None and not (ds <= pos <= de):
            errs.append(f"{tag} boundary {b.get('id')} is inside_domain but aa {pos} is outside "
                        f"instance {iid} ({ds}–{de})")
        if cls == "near_domain_edge" and signed is not None and abs(signed) > thr:
            errs.append(f"{tag} boundary {b.get('id')} is near_domain_edge but |{signed}| > {thr}")
    return errs


def _cross_check_sources(model: Dict[str, Any], core: Path, tag: str) -> List[str]:
    errs: List[str] = []
    sp, pid = model.get("species_id"), model.get("protein_id")

    # domains vs domain_features.tsv (representative layer)
    src_domains = {
        (_int(r.get("start_aa")), _int(r.get("end_aa")))
        for r in _read_tsv(core / "domain_features.tsv")
        if r.get("species_id") == sp and r.get("protein_id") == pid
        and str(r.get("layer", "domain")).lower() in ("domain", "representative_domain")
    }
    model_domains = {(_int(d.get("start")), _int(d.get("end"))) for d in model.get("representative_domains", [])}
    if model_domains - src_domains:
        errs.append(f"{tag} representative_domains not backed by domain_features.tsv: {model_domains - src_domains}")

    # tm vs tm_features.tsv
    src_tm = {
        (_int(r.get("start_aa")), _int(r.get("end_aa")))
        for r in _read_tsv(core / "tm_features.tsv")
        if r.get("species_id") == sp and r.get("protein_id") == pid
        and _int(r.get("start_aa")) is not None
    }
    model_tm = {(_int(t.get("start")), _int(t.get("end"))) for t in model.get("tm_regions", [])}
    if model_tm != src_tm:
        errs.append(f"{tag} tm_regions {model_tm} disagree with pyTMHMM table {src_tm}")

    return errs


def validate_index(index: Dict[str, Any], *, core_dir: Optional[Path] = None) -> List[str]:
    errors: List[str] = []
    if index.get("schema_version") != 1:
        errors.append(f"index schema_version must be 1 (got {index.get('schema_version')})")
    # Rule 3: one primary reference per species, and every model named.
    #
    # The rule is not "one model per species". A gene whose analysis is about
    # alternative proteins has several real models for one species, and FGFR2's IIIb
    # and IIIc are two proteins rather than two views of one. What must hold is the
    # hierarchy: each model says what it is (``model_role``), each is individually
    # addressable (``model_id``), and exactly one of a species' models is the
    # primary reference a comparative row speaks for. That is stricter than the
    # original single-model check, not looser: a duplicate is still caught, and now
    # a silently-chosen comparative row is too.
    models = index.get("models", [])
    for m in models:
        errors += validate_model(m, core_dir=core_dir)
    errors += model_roles.role_errors(models)

    proteins: Dict[str, List[str]] = {}
    for m in models:
        proteins.setdefault(str(m.get("species_id") or "?"), []).append(
            str(m.get("protein_id") or ""))
    for sp, pids in proteins.items():
        if len(set(pids)) != len(pids):
            errors.append(f"[{sp}] the same protein appears in {len(pids)} models: "
                          f"{sorted(pids)}")
    return errors


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description='validate_protein_coordinate_model.py — loud validator for the coordinate contract.')
    ap.add_argument("run_dir", nargs="?", help="run directory (builds the model)")
    ap.add_argument("--model-json", help="validate a pre-built model index JSON instead")
    ap.add_argument("--write", help="write the built model index JSON to this path")
    args = ap.parse_args(argv)

    core_dir: Optional[Path] = None
    if args.model_json:
        index = json.loads(Path(args.model_json).read_text(encoding="utf-8"))
    elif args.run_dir:
        run_dir = Path(args.run_dir)
        core_dir = run_dir / "results" / "core_gene_analysis"
        index = build_models_for_run(run_dir)
        if args.write:
            Path(args.write).write_text(json.dumps(index, indent=2), encoding="utf-8")
    else:
        ap.error("provide a run_dir or --model-json")
        return 2

    errors = validate_index(index, core_dir=core_dir)
    n = index.get("n_models", len(index.get("models", [])))
    if errors:
        print(f"FAIL — {len(errors)} coordinate-model violation(s) across {n} model(s):", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        raise CoordinateModelError(f"{len(errors)} coordinate-model violation(s)")
    print(f"OK — {n} coordinate model(s) validated; no coordinate inconsistencies.")
    for m in index.get("models", []):
        print(
            f"  · {m['species_id']}/{m['protein_id']} len={m['protein_length']}aa "
            f"exons={len(m['exons'])} domains={len(m['representative_domains'])} "
            f"tm={len(m['tm_regions'])} boundaries={len(m['exon_boundaries'])} status={m['status']}"
        )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CoordinateModelError:
        sys.exit(1)
