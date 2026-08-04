#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

MIN_CASSETTE_START = 150

# Sequence-calibrated decision thresholds (mirror the human-calibrated classifier).
MIN_IDENTITY_FLOOR = 0.45     # min identity*coverage for a usable assignment
DIRECTION_MARGIN = 0.10       # IIIb vs IIIc separation for a confident call
MIN_COVERAGE = 0.50           # min cassette-reference coverage for a plausible slot

ISO = ("IIIb", "IIIc")

AUDIT_COLS = [
    "species", "transcript_id", "protein_id", "upstream_isoform",
    "fallback_attempted", "fallback_validated_exon_type", "fallback_final_isoform_label",
    "human_IIIb_identity", "human_IIIc_identity",
    "IIIb_marker_support", "IIIc_marker_support", "MSA_discriminating_support",
    "fallback_decision", "fallback_confidence", "fallback_warning",
    "resolved_cassette_protein_start_aa", "resolved_cassette_protein_end_aa",
    "cassette_reference_coverage", "evidence_used",
]

ALLOWED_DECISIONS = {
    "assigned_by_sequence_calibrated_fallback",
    "assigned_by_validated_rescue_override",
    "manual_review_required",
    "unresolved_no_safe_mapping",
}

def smith_waterman(query: str, target: str, match: int = 2, mismatch: int = -1,
                   gap: int = -2) -> Tuple[int, int, float, float]:
    q = (query or "").upper()
    t = (target or "").upper()
    n, m = len(q), len(t)
    if n == 0 or m == 0:
        return (0, 0, 0.0, 0.0)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best = bi = bj = 0
    for i in range(1, n + 1):
        qi = q[i - 1]
        Hi, Hp = H[i], H[i - 1]
        for j in range(1, m + 1):
            v = Hp[j - 1] + (match if qi == t[j - 1] else mismatch)
            d = Hp[j] + gap
            if d > v:
                v = d
            e = Hi[j - 1] + gap
            if e > v:
                v = e
            if v < 0:
                v = 0
            Hi[j] = v
            if v > best:
                best, bi, bj = v, i, j
    if best == 0:
        return (0, 0, 0.0, 0.0)
    i, j = bi, bj
    t_end = j
    q_end = i
    n_match = aln_len = 0
    while i > 0 and j > 0 and H[i][j] > 0:
        cur = H[i][j]
        diag = H[i - 1][j - 1] + (match if q[i - 1] == t[j - 1] else mismatch)
        if cur == diag:
            aln_len += 1
            if q[i - 1] == t[j - 1]:
                n_match += 1
            i -= 1
            j -= 1
        elif cur == H[i - 1][j] + gap:
            i -= 1
            aln_len += 1
        else:
            j -= 1
            aln_len += 1
    t_start = j + 1
    q_start = i + 1
    identity = (n_match / aln_len) if aln_len else 0.0
    coverage = ((q_end - q_start + 1) / n) if n else 0.0
    return (t_start, t_end, round(min(1.0, identity), 4), round(min(1.0, coverage), 4))


def _score(identity: float, coverage: float) -> float:
    return round(identity * coverage, 6)


def read_single_fasta(path: Path) -> str:
    parts: List[str] = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if parts:
                    break
                continue
            parts.append("".join(c for c in line.strip() if c.isalpha()))
    return "".join(parts).upper()


def parse_header_kv(header: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for tok in header.split("|"):
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def load_proteins_by_species(path: Optional[Path]):
    by_sp: Dict[str, List[Dict[str, str]]] = {}
    if not path or not Path(path).exists():
        return by_sp
    cur_meta: Optional[Dict[str, str]] = None
    buf: List[str] = []

    def flush():
        if cur_meta is None:
            return
        seq = "".join(buf).upper().replace("*", "")
        sp = (cur_meta.get("species") or "").lower()
        if not sp or not seq:
            return
        by_sp.setdefault(sp, []).append({
            "protein_id": cur_meta.get("protein", ""),
            "transcript": cur_meta.get("transcript", ""),
            "isoform": cur_meta.get("isoform", ""),
            "role": cur_meta.get("role", ""),
            "seq": seq,
        })

    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                cur_meta = parse_header_kv(line[1:].strip())
                buf = []
            else:
                buf.append("".join(c for c in line.strip() if c.isalpha()))
    flush()
    return by_sp


def read_tsv(path: Optional[Path]) -> List[Dict[str, str]]:
    if not path or not Path(path).exists():
        return []
    with Path(path).open(encoding="utf-8", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return [dict(x) for x in r] if r.fieldnames else []


def write_tsv(path: Path, rows: List[Dict[str, object]], cols: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in cols})


def to_int(v) -> Optional[int]:
    try:
        if v in (None, "") or str(v).lower() == "nan":
            return None
        return int(float(str(v)))
    except Exception:
        return None


def load_rescue_overrides(path: Optional[Path]):
    out: Dict[Tuple[str, str], str] = {}
    for r in read_tsv(path):
        sp = (r.get("species") or r.get("species_canonical") or "").lower()
        iso = r.get("isoform") or r.get("final_isoform_label") or r.get("validated_exon_type") or ""
        decision = (r.get("rescue_decision") or r.get("maximal_rescue_decision")
                    or r.get("decision") or "").lower()
        status = (r.get("validation_status") or r.get("rescue_status") or "").lower()
        validated = ("validated" in decision or "validated" in status
                     or "rescued_ok" in decision or "accept" in decision)
        et = r.get("validated_exon_type") or iso
        if sp and et in ISO and validated:
            out[(sp, et)] = et
    return out


def load_markers(path: Optional[Path]):
    out: Dict[str, Dict[str, str]] = {}
    for r in read_tsv(path):
        sp = (r.get("species") or r.get("species_canonical") or "").lower()
        if not sp:
            continue
        prefers = (r.get("protein_prefers") or "").strip()
        idb = r.get("human_IIIb_identity", "")
        idc = r.get("human_IIIc_identity", "")
        if sp not in out or r.get("protein_validation_is_decisive", "").lower() in ("true", "1", "yes"):
            out[sp] = {
                "IIIb_marker_support": f"prefers={prefers};human_IIIb_identity={idb}" if prefers else "",
                "IIIc_marker_support": f"prefers={prefers};human_IIIc_identity={idc}" if prefers else "",
                "prefers": prefers,
            }
    return out


def load_msa_discriminating(path: Optional[Path]) -> bool:
    return bool(read_tsv(path))


def needs_fallback(row: Dict[str, str]) -> bool:
    iso = (row.get("inferred_isoform") or "").strip()
    if iso not in ISO:
        return False
    pid = (row.get("protein_id") or "").strip()
    start = to_int(row.get("native_protein_start_aa"))
    plen = to_int(row.get("protein_length_aa"))
    return (not pid) or (start is None) or (start <= 1) or (plen is None)


def best_protein_for_iso(prots: List[Dict[str, str]], ref: str):
    best = None
    for p in prots:
        ts, te, ident, cov = smith_waterman(ref, p["seq"])
        sc = _score(ident, cov)
        if best is None or sc > best[5]:
            best = (p, ts, te, ident, cov, sc)
    return best


def assign_row(row, prots, iiib_ref, iiic_ref, rescue, markers, msa_avail):
    sp = (row.get("species_canonical") or row.get("species") or "").lower()
    iso = (row.get("inferred_isoform") or "").strip()
    upstream = (row.get("upstream_label") or row.get("isoform")
                or (prots[0]["isoform"] if prots else "") or "unclassified")
    audit = {
        "species": row.get("species_canonical", row.get("species", "")),
        "transcript_id": row.get("transcript_id_source", ""),
        "protein_id": row.get("protein_id", ""),
        "upstream_isoform": upstream or "unclassified",
        "fallback_attempted": "true",
        "fallback_validated_exon_type": "",
        "fallback_final_isoform_label": "",
        "human_IIIb_identity": "", "human_IIIc_identity": "",
        "IIIb_marker_support": "", "IIIc_marker_support": "",
        "MSA_discriminating_support": "available" if msa_avail else "not_available_pre_msa",
        "fallback_decision": "unresolved_no_safe_mapping",
        "fallback_confidence": "none",
        "fallback_warning": "",
        "resolved_cassette_protein_start_aa": "", "resolved_cassette_protein_end_aa": "",
        "cassette_reference_coverage": "", "evidence_used": "",
    }

    mk = markers.get(sp, {})
    audit["IIIb_marker_support"] = mk.get("IIIb_marker_support", "") or "not_available"
    audit["IIIc_marker_support"] = mk.get("IIIc_marker_support", "") or "not_available"

    if not prots:
        audit["fallback_warning"] = "no_species_protein_available_for_sequence_fallback"
        return audit, None

    best_b = best_protein_for_iso(prots, iiib_ref)
    best_c = best_protein_for_iso(prots, iiic_ref)
    sb = best_b[5] if best_b else 0.0
    sc = best_c[5] if best_c else 0.0
    audit["human_IIIb_identity"] = best_b[3] if best_b else 0.0
    audit["human_IIIc_identity"] = best_c[3] if best_c else 0.0

    seq_prefers = "ambiguous"
    if max(sb, sc) >= MIN_IDENTITY_FLOOR:
        if sb - sc >= DIRECTION_MARGIN:
            seq_prefers = "IIIb"
        elif sc - sb >= DIRECTION_MARGIN:
            seq_prefers = "IIIc"

    evidence = []

    decision = None
    if (sp, iso) in rescue:
        decision = "assigned_by_validated_rescue_override"
        evidence.append("validated_rescue_override")

    ref = iiib_ref if iso == "IIIb" else iiic_ref
    p, ts, te, _ident, cov, score = best_protein_for_iso(prots, ref)
    best_overall = max(sb, sc)

    if mk.get("prefers"):
        evidence.append("protein_marker_support")
    if seq_prefers in ISO and seq_prefers != iso:
        audit["fallback_warning"] = (
            f"sequence_prefers_{seq_prefers}_(single_available_protein);"
            f"exon_structure_call_{iso}_preserved")

    audit["fallback_validated_exon_type"] = seq_prefers if seq_prefers in ISO else iso
    audit["fallback_final_isoform_label"] = iso
    audit["cassette_reference_coverage"] = cov
    audit["evidence_used"] = ";".join(evidence) or "none"

    plausible_slot = (ts >= MIN_CASSETTE_START)
    real_cassette_region = (best_overall >= MIN_IDENTITY_FLOOR)

    if decision == "assigned_by_validated_rescue_override" and plausible_slot:
        conf = "high"
    elif plausible_slot and real_cassette_region:
        decision = "assigned_by_sequence_calibrated_fallback"
        evidence.append("human_cassette_reference_similarity")
        audit["evidence_used"] = ";".join(evidence)
        if score >= MIN_IDENTITY_FLOOR and cov >= 0.8:
            conf = "high"
        elif cov >= MIN_COVERAGE or best_overall >= 0.8:
            conf = "medium"
        else:
            conf = "low"
    elif plausible_slot and not real_cassette_region:
        decision = "manual_review_required"
        conf = "low"
    else:
        decision = "unresolved_no_safe_mapping"
        conf = "none"

    audit["fallback_decision"] = decision
    audit["fallback_confidence"] = conf

    if decision in ("assigned_by_sequence_calibrated_fallback",
                    "assigned_by_validated_rescue_override"):
        audit["resolved_cassette_protein_start_aa"] = ts
        audit["resolved_cassette_protein_end_aa"] = te
        patch = {
            "protein_id": p["protein_id"],
            "transcript_id_source": row.get("transcript_id_source") or p["transcript"],
            "protein_length_aa": len(p["seq"]),
            "native_protein_start_aa": ts, "native_protein_end_aa": te,
            "native_protein_center_aa": (ts + te) // 2,
            "native_protein_length_aa": len(p["seq"]),
        }
        return audit, patch

    if decision == "manual_review_required":
        audit["resolved_cassette_protein_start_aa"] = ts
        audit["resolved_cassette_protein_end_aa"] = te
        patch = {
            "protein_id": p["protein_id"],
            "protein_length_aa": len(p["seq"]),
            "native_protein_start_aa": ts, "native_protein_end_aa": te,
            "native_protein_center_aa": (ts + te) // 2,
            "native_protein_length_aa": len(p["seq"]),
            "__review__": "true",
        }
        return audit, patch

    if not audit["fallback_warning"]:
        audit["fallback_warning"] = "no_plausible_mid_protein_cassette_slot_found"
    patch = {
        "native_protein_start_aa": "", "native_protein_end_aa": "",
        "native_protein_center_aa": "", "__review__": "true",
    }
    return audit, patch


def summarise(audit_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    from collections import Counter
    dec = Counter(str(r["fallback_decision"]) for r in audit_rows)
    out = []
    for k in sorted(ALLOWED_DECISIONS):
        out.append({"fallback_decision": k, "n_rows": dec.get(k, 0)})
    out.append({"fallback_decision": "total_fallback_rows", "n_rows": len(audit_rows)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Sequence-calibrated fallback for unclassified FGFR2 IIIb/IIIc "
                    "cassettes (no aa-1 coordinate defaults; gate stays strict).")
    ap.add_argument("--coordinate_audit", type=Path, required=True)
    ap.add_argument("--proteins", type=Path, required=True)
    ap.add_argument("--human_iiib_segment", type=Path,
                    default=Path("references/fgfr2_iii_segments/human_FGFR2_IIIb_segment.fasta"))
    ap.add_argument("--human_iiic_segment", type=Path,
                    default=Path("references/fgfr2_iii_segments/human_FGFR2_IIIc_segment.fasta"))
    ap.add_argument("--protein_validation_summary", type=Path, default=None)
    ap.add_argument("--rescue_overrides", type=Path, default=None)
    ap.add_argument("--msa_discriminating", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    iiib_ref = read_single_fasta(args.human_iiib_segment)
    iiic_ref = read_single_fasta(args.human_iiic_segment)
    if not iiib_ref or not iiic_ref:
        print("[FAIL] human IIIb/IIIc cassette reference segments not found/empty.", file=sys.stderr)
        return 2

    prots_by_sp = load_proteins_by_species(args.proteins)
    rescue = load_rescue_overrides(args.rescue_overrides)
    markers = load_markers(args.protein_validation_summary)
    msa_avail = load_msa_discriminating(args.msa_discriminating)

    with args.coordinate_audit.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]
    if not fieldnames:
        print("[FAIL] coordinate audit is empty.", file=sys.stderr)
        return 2

    audit_rows: List[Dict[str, object]] = []
    n_patched = n_review = n_unresolved = 0
    for row in rows:
        if not needs_fallback(row):
            continue
        sp = (row.get("species_canonical") or row.get("species") or "").lower()
        prots = prots_by_sp.get(sp, [])
        audit, patch = assign_row(row, prots, iiib_ref, iiic_ref, rescue, markers, msa_avail)
        assert audit["fallback_decision"] in ALLOWED_DECISIONS, audit["fallback_decision"]
        audit_rows.append(audit)
        if patch is None:
            n_unresolved += 1
            continue
        review = patch.pop("__review__", "")
        for k, v in patch.items():
            if k in row:
                row[k] = v
            else:
                pass
        note = f"unclassified_isoform_fallback:{audit['fallback_decision']}"
        for col in ("iii_slot_coordinate_note", "coordinate_warning", "resolver_warning"):
            if col in row:
                prev = (row.get(col) or "").strip()
                row[col] = (prev + ";" if prev else "") + note
                break
        if review == "true":
            n_review += 1
            for col in ("display_uncertainty_class", "annotation_review_state"):
                if col in row:
                    row[col] = "review_unclassified_isoform_fallback"
        else:
            n_patched += 1

    write_tsv(args.coordinate_audit, rows, fieldnames)

    write_tsv(args.outdir / "unclassified_isoform_fallback_audit.tsv", audit_rows, AUDIT_COLS)
    write_tsv(args.outdir / "unclassified_isoform_fallback_summary.tsv",
              summarise(audit_rows), ["fallback_decision", "n_rows"])

    aa1 = [a["species"] for a in audit_rows
           if str(a.get("resolved_cassette_protein_start_aa")) == "1"]
    print(f"[OK] unclassified-isoform fallback: {len(audit_rows)} cassette rows examined "
          f"| resolved={n_patched} review={n_review} unresolved={n_unresolved}")
    if aa1:
        print(f"[FAIL] fallback produced aa-1 coordinate for {aa1}", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
