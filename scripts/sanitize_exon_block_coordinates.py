#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO = Path(__file__).resolve().parent.parent


def display_path(path) -> str:
    p = Path(path)
    try:
        return str(p.resolve().relative_to(REPO.resolve()))
    except Exception:
        return str(p)
import os as _os  # run-folder path override (RESULTS_DIR/BASE); legacy default preserved
BASE = Path(_os.environ.get("FGFR2_RESULTS_DIR") or _os.environ.get("RESULTS_DIR")
            or _os.environ.get("BASE") or (REPO / "results" / "final_30_until_interpro_prepare"))
CLOSURE = BASE / "13_final_pre_interpro_closure"
POST = BASE / "15_exon_domain_boundary_post_interpro"

MANIFEST = CLOSURE / "freeze" / "final_pre_interpro_sequence_manifest.tsv"
TRUTH = CLOSURE / "final_pre_interpro_truth_table.tsv"
FEATURES = POST / "tables" / "exon_domain_architecture_features.tsv"
CDS_FEATURES = BASE / "02_models" / "cds_features.tsv"
CASSETTE_MAP = (BASE / "09_paper_ready_qc_v2_9" / "figures_v2_22_final_qc_display"
                / "fgfr2_cassette_cds_block_map.tsv")

AUDIT_OUT = POST / "tables" / "exon_block_length_consistency_audit.tsv"
OVERRIDE_OUT = POST / "tables" / "exon_block_reconstruction_overrides.json"

# keys already handled by reconstruct_exon_blocks_post_interpro.py (the three
# major coordinate-artifact cases). Preserved verbatim; never re-processed here.
PRESERVE_KEYS = {
    ("canis_lupus_familiaris", "IIIb"),
    ("gorilla_gorilla_gorilla", "IIIb"),
    ("xenopus_tropicalis", "IIIb"),
}

MINOR_MAX = 2       # +1/+2 -> clamp
MODERATE_MAX = 15   # 3..15 moderate, >15 severe (both reconstruct/hide)


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def to_int(v, default=None):
    try:
        return int(float(str(v).strip()))
    except (ValueError, TypeError):
        return default


def _norm_tx(tx: str) -> str:
    return (tx or "").replace("rna-", "").split(".")[0]


def load_native(cds_rows: List[dict], transcript_id: str,
                protein_id: str) -> List[dict]:
    base = _norm_tx(transcript_id)
    pid_base = (protein_id or "").split(".")[0]
    blocks: List[dict] = []
    for r in cds_rows:
        tsrc = r.get("transcript_id_source", "")
        psrc = r.get("translation_id_source", "")
        if _norm_tx(tsrc) != base and (not pid_base or psrc.split(".")[0] != pid_base):
            continue
        pa, pe, rank = (to_int(r.get("protein_start_aa")), to_int(r.get("protein_end_aa")),
                        to_int(r.get("cds_rank")))
        if pa is None or pe is None or rank is None:
            continue
        blocks.append({"rank": rank, "exon_id": r.get("cds_id_source", ""),
                       "gstart": to_int(r.get("start")), "gend": to_int(r.get("end")),
                       "pa": pa, "pe": pe})
    blocks.sort(key=lambda b: b["rank"])
    return blocks


def native_cassette_rank(blocks: List[dict], g0: Optional[int],
                         g1: Optional[int]) -> Optional[dict]:
    if g0 is None or g1 is None:
        return None
    best, best_ov = None, 0
    for b in blocks:
        if b["gstart"] is None or b["gend"] is None:
            continue
        ov = min(b["gend"], g1) - max(b["gstart"], g0)
        if ov > best_ov:
            best, best_ov = b, ov
    return best


def main() -> int:
    manifest = {(r["species"], r["isoform"]): r for r in read_tsv(MANIFEST)}
    primary = {k for k, r in manifest.items()
               if str(r.get("included_in_primary_interpro", "")).lower() == "true"}
    truth = {(r["species"], r["isoform"]): r for r in read_tsv(TRUTH)}
    feat_rows = read_tsv(FEATURES)
    cds_rows = read_tsv(CDS_FEATURES)
    cassette_map = {(r["species"], r["isoform"]): r for r in read_tsv(CASSETTE_MAP)}

    # current coding-exon blocks + cassette per (species, isoform)
    exons: Dict[Tuple[str, str], List[dict]] = {}
    cassette_block: Dict[Tuple[str, str], dict] = {}
    meta: Dict[Tuple[str, str], dict] = {}
    for r in feat_rows:
        key = (r.get("species", ""), r.get("isoform", ""))
        ft = r.get("feature_type", "")
        s, e = to_int(r.get("start_aa")), to_int(r.get("end_aa"))
        if ft == "coding_exon":
            exons.setdefault(key, []).append(
                {"label": r.get("feature_label", ""), "start": s, "end": e,
                 "source": r.get("source", "")})
            meta.setdefault(key, {"tx": r.get("transcript_id", ""),
                                  "pid": r.get("protein_id", ""),
                                  "L": to_int(r.get("protein_length"))})
        elif ft in ("IIIb_slot", "IIIc_slot"):
            cassette_block[key] = {"label": r.get("feature_label", ""),
                                   "start": s, "end": e, "type": ft,
                                   "source": r.get("source", "")}
            meta.setdefault(key, {"tx": r.get("transcript_id", ""),
                                  "pid": r.get("protein_id", ""),
                                  "L": to_int(r.get("protein_length"))})

    # start from existing overrides so the three preserved cases survive
    overrides: Dict[str, dict] = {}
    if OVERRIDE_OUT.exists():
        try:
            overrides = json.loads(OVERRIDE_OUT.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            overrides = {}

    audit_rows: List[dict] = []
    counts = {"audited": 0, "no_issue": 0, "minor": 0, "reconstructed": 0,
              "hidden": 0, "cassette_only": 0, "preserved": 0}
    moderate_severe: List[dict] = []

    for key in sorted(primary):
        sp, iso = key
        m = meta.get(key, {})
        tx, pid, L = m.get("tx", ""), m.get("pid", ""), m.get("L")
        blks = sorted(exons.get(key, []), key=lambda b: (b["start"] is None, b["start"]))
        counts["audited"] += 1

        tr = truth.get(key, {})
        iso_label = tr.get("final_isoform_label", iso)

        # the three preserved coordinate-artifact cases -> already sanitary
        if key in PRESERVE_KEYS:
            counts["preserved"] += 1
            ov = overrides.get(f"{sp}|{iso}", {})
            status = ov.get("final_display_status", "preserved")
            for i, b in enumerate(blks, 1):
                audit_rows.append(_audit_row(sp, iso_label, tx, pid, L, i, b, 0,
                                             b["source"], "", "no_issue", "keep",
                                             b["start"], b["end"],
                                             f"preserved earlier reconstruction ({status})"))
            if not blks:
                audit_rows.append(_audit_row(sp, iso_label, tx, pid, L, 0, None, 0,
                                             "reconstruction_override", "", "no_issue",
                                             ov.get("final_display_status", "cassette_only_display"),
                                             "", "", "coding exons hidden; validated cassette shown"))
            continue

        max_end = max((b["end"] for b in blks if b["end"] is not None), default=0)
        overflow = (max_end - L) if (L and max_end) else 0

        if not blks:
            counts["no_issue"] += 1
            continue

        if overflow <= 0:
            counts["no_issue"] += 1
            for i, b in enumerate(blks, 1):
                audit_rows.append(_audit_row(sp, iso_label, tx, pid, L, i, b,
                                             max(0, (b["end"] or 0) - (L or 0)),
                                             b["source"], "", "no_issue", "keep",
                                             b["start"], b["end"], ""))
            continue

        if overflow <= MINOR_MAX:
            # clamp final block end to protein length (generator applies clamp)
            counts["minor"] += 1
            for i, b in enumerate(blks, 1):
                ovi = max(0, (b["end"] or 0) - (L or 0))
                if ovi > 0:
                    audit_rows.append(_audit_row(
                        sp, iso_label, tx, pid, L, i, b, ovi, b["source"],
                        "codon_boundary_rounding",
                        "minor_plus1_plus2_rounding",
                        "clamp_to_protein_length_with_minor_flag",
                        b["start"], L,
                        f"final block end clamped {b['end']}->{L} (+{ovi} aa rounding)"))
                else:
                    audit_rows.append(_audit_row(sp, iso_label, tx, pid, L, i, b, 0,
                                                 b["source"], "", "no_issue", "keep",
                                                 b["start"], b["end"], ""))
            continue

        # ---- moderate / severe -> reconstruct from native CDS, else hide ----
        issue_class = "moderate_overflow" if overflow <= MODERATE_MAX else "severe_overflow"
        native = load_native(cds_rows, tx, pid)
        cm = cassette_map.get(key, {})
        g0, g1 = to_int(cm.get("resolver_genomic_start")), to_int(cm.get("resolver_genomic_end"))
        rec = {"species": sp, "isoform": iso_label, "tx": tx, "pid": pid, "L": L,
               "overflow": overflow, "issue_class": issue_class}

        if native:
            cass = native_cassette_rank(native, g0, g1)
            cass_rank = cass["rank"] if cass else None
            # fall back: match existing cassette block by AA overlap
            if cass_rank is None and key in cassette_block:
                cb = cassette_block[key]
                best_ov, best_r = 0, None
                for b in native:
                    ov = min(b["pe"], cb["end"] or 0) - max(b["pa"], cb["start"] or 0)
                    if ov > best_ov:
                        best_ov, best_r = ov, b["rank"]
                cass_rank = best_r
            out_blocks = []
            _slot_type = "IIIb_slot" if iso.upper() == "IIIB" else "IIIc_slot"
            for i, b in enumerate(native, 1):
                s = max(1, b["pa"])
                e = min(L, b["pe"]) if L else b["pe"]
                if e < s:
                    e = s
                is_cass = (b["rank"] == cass_rank)
                out_blocks.append({"number": i, "label": f"CDS{b['rank']}",
                                   "exon_id": b["exon_id"], "start": s, "end": e,
                                   "is_cassette": is_cass})
            cass_coords = None
            if cass_rank is not None:
                cb = next(b for b in out_blocks if b["is_cassette"])
                cass_coords = {"start": cb["start"], "end": cb["end"]}
            max_native = max(b["end"] for b in out_blocks)
            note = (f"coding-exon blocks reconstructed from native local CDS features "
                    f"(figure3C blocks projected from a different transcript, overflow "
                    f"+{overflow} aa); native max end {max_native} <= length {L}")
            overrides[f"{sp}|{iso}"] = {
                "final_display_status": "native_exon_blocks_reconstructed",
                "exon_blocks": out_blocks,
                "cassette": cass_coords,
                "recon_note": "coding-exon blocks reconstructed from native CDS "
                              "(display-coordinate sanitation); biology unchanged",
                "exon_block_source": f"cds_features_native ({tx})",
            }
            counts["reconstructed"] += 1
            action = "reconstruct_from_native_cds"
            # per-exon audit against native (map original number -> native rank i)
            for i, b in enumerate(blks, 1):
                nb = out_blocks[i - 1] if i - 1 < len(out_blocks) else None
                fe = (nb["end"] if nb else "")
                fs = (nb["start"] if nb else "")
                ovi = max(0, (b["end"] or 0) - (L or 0))
                audit_rows.append(_audit_row(
                    sp, iso_label, tx, pid, L, i, b, ovi, b["source"],
                    f"native_cds_rank_{nb['label']}" if nb else "", issue_class,
                    action, fs, fe, note if i == 1 else "reconstructed native block"))
            rec.update({"action": action, "native": True,
                        "cassette_rank": cass_rank, "max_native": max_native})
        else:
            # hide untrusted coding blocks; keep validated cassette if present
            has_cass = key in cassette_block
            status = "cassette_only_display" if has_cass else "hide_untrusted_exon_block"
            overrides[f"{sp}|{iso}"] = {
                "final_display_status": ("cassette_only_high_confidence" if has_cass
                                         else "exon_blocks_hidden_untrusted"),
                "exon_blocks": [],
                "cassette": ({"start": cassette_block[key]["start"],
                              "end": cassette_block[key]["end"]} if has_cass else None),
                "recon_note": "coding-exon blocks hidden (no native CDS coordinates "
                              "for this transcript; display-coordinate sanitation)",
                "exon_block_source": "none_native_cds_unavailable",
            }
            counts["hidden"] += 1
            if has_cass:
                counts["cassette_only"] += 1
            action = status
            for i, b in enumerate(blks, 1):
                ovi = max(0, (b["end"] or 0) - (L or 0))
                audit_rows.append(_audit_row(
                    sp, iso_label, tx, pid, L, i, b, ovi, b["source"], "",
                    issue_class, action, "", "",
                    "native CDS unavailable; block hidden"))
            rec.update({"action": action, "native": False})
        moderate_severe.append(rec)

    _write_audit(audit_rows)
    OVERRIDE_OUT.write_text(json.dumps(overrides, indent=2, ensure_ascii=False),
                            encoding="utf-8")

    print(f"[ok] wrote {display_path(AUDIT_OUT)} ({len(audit_rows)} exon rows)")
    print(f"[ok] wrote {display_path(OVERRIDE_OUT)} ({len(overrides)} override keys)")
    print(f"[summary] audited={counts['audited']} no_issue={counts['no_issue']} "
          f"minor_clamp={counts['minor']} reconstructed={counts['reconstructed']} "
          f"hidden={counts['hidden']} (cassette_only={counts['cassette_only']}) "
          f"preserved={counts['preserved']}")
    print("[moderate/severe cases]")
    for r in moderate_severe:
        print(f"  {r['species']} {r['isoform']}: +{r['overflow']} ({r['issue_class']}) "
              f"-> {r['action']}"
              + (f" cassette_rank={r.get('cassette_rank')} max_native={r.get('max_native')}"
                 if r.get("native") else ""))
    return 0


def _audit_row(sp, iso, tx, pid, L, num, b, overflow, source, phase_info,
               issue_class, action, fstart, fend, notes) -> dict:
    return {
        "species": sp, "isoform": iso, "transcript_id": tx, "protein_id": pid,
        "protein_length": L if L is not None else "",
        "exon_number": num if num else "",
        "exon_start_aa": (b["start"] if b else ""),
        "exon_end_aa": (b["end"] if b else ""),
        "overflow_aa": overflow,
        "source": source,
        "phase_or_rounding_info_if_available": phase_info,
        "issue_class": issue_class,
        "action": action,
        "final_exon_start_aa": fstart,
        "final_exon_end_aa": fend,
        "notes": notes,
    }


def _write_audit(rows: List[dict]) -> None:
    cols = ["species", "isoform", "transcript_id", "protein_id", "protein_length",
            "exon_number", "exon_start_aa", "exon_end_aa", "overflow_aa", "source",
            "phase_or_rounding_info_if_available", "issue_class", "action",
            "final_exon_start_aa", "final_exon_end_aa", "notes"]
    AUDIT_OUT.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT_OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
