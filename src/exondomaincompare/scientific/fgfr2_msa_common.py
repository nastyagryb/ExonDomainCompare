#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

MODULE_REL = "12_msa_boundary_robustness_pre_interpro"
CLOSURE_REL = "13_final_pre_interpro_closure"
SUBDIRS = ["inputs", "alignments", "maps", "conservation", "robustness", "splice_qc",
           "protein_integrity", "review_diagnostics", "figures", "tables", "captions",
           "metadata", "synteny"]
CLOSURE_SUBDIRS = ["MSA", "figures", "tables", "reports", "gates", "freeze", "archive", "metadata"]

AA_VALID = set("ACDEFGHIKLMNPQRSTVWY")
AA_AMBIG = set("XBZJUO")  # accepted ambiguous / rare residues
GAP_CHARS = set("-.")

# stable IIIb/IIIc colors (identical across all thesis figures)
C_IIIB = "#0072B2"
C_IIIC = "#E69F00"
C_REVIEW = "#D55E00"


def now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def module_dir(base: Path) -> Path:
    return Path(base) / MODULE_REL


def closure_dir(base: Path) -> Path:
    return Path(base) / CLOSURE_REL


def ensure_closure_dirs(base: Path) -> Dict[str, Path]:
    cd = closure_dir(base)
    out = {}
    for sd in CLOSURE_SUBDIRS:
        p = cd / sd
        p.mkdir(parents=True, exist_ok=True)
        out[sd] = p
    return out


def ensure_module_dirs(base: Path) -> Dict[str, Path]:
    md = module_dir(base)
    out = {}
    for sd in SUBDIRS:
        p = md / sd
        p.mkdir(parents=True, exist_ok=True)
        out[sd] = p
    return out


# ---------------------------------------------------------------------------
# TSV / FASTA IO
# ---------------------------------------------------------------------------
def read_tsv(path: Path) -> List[Dict[str, str]]:
    path = Path(path)
    if not path.exists():
        return []
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, rows: List[Dict[str, object]], fields: List[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def read_fasta(path: Path) -> "List[Tuple[str, str]]":
    path = Path(path)
    items: List[Tuple[str, str]] = []
    if not path.exists():
        return items
    cur_id, cur_seq = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if cur_id is not None:
                    items.append((cur_id, "".join(cur_seq)))
                cur_id = line[1:].split()[0] if len(line) > 1 else ""
                cur_seq = []
            elif line:
                cur_seq.append(line.strip())
    if cur_id is not None:
        items.append((cur_id, "".join(cur_seq)))
    return items


def write_fasta(path: Path, items: Iterable[Tuple[str, str]], width: int = 60) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for sid, seq in items:
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
    except OSError:
        return ""
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Locating final-pipeline inputs
# ---------------------------------------------------------------------------
def locate(base: Path, name: str, hint: str = "") -> Optional[Path]:
    matches = sorted(Path(base).rglob(name))
    # never descend into our own module outputs when locating upstream inputs
    matches = [m for m in matches if MODULE_REL not in str(m)]
    if not matches:
        return None
    if hint:
        for m in matches:
            if hint in str(m):
                return m
    return sorted(matches, key=lambda p: len(p.parts))[0]


def require(base: Path, name: str, hint: str = "") -> Path:
    p = locate(base, name, hint)
    if p is None:
        raise RuntimeError(f"Required input not found under {base}: {name}")
    return p


# ---------------------------------------------------------------------------
# Small shared vocab / utilities
# ---------------------------------------------------------------------------
def to_int(v: object, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(float(str(v).strip()))
    except Exception:
        return default


def to_float(v: object, default: Optional[float] = None) -> Optional[float]:
    try:
        return float(str(v).strip())
    except Exception:
        return default


RECOMMENDED_USE_TOKEN = {
    "main_text_primary_claim": "main_figure",
    "main_text_with_footnote": "main_with_footnote",
    "supplementary_only": "supplement",
}


def recommended_use_token(recommended_use: str) -> str:
    return RECOMMENDED_USE_TOKEN.get((recommended_use or "").strip(), "review")


def is_main_use(recommended_use: str) -> bool:
    return (recommended_use or "").strip() in (
        "main_text_primary_claim", "main_text_with_footnote")


def load_label_reconciliation(base: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    p = module_dir(base) / "maps" / "fgfr2_exon_type_label_reconciliation.tsv"
    out: Dict[Tuple[str, str], Dict[str, str]] = {}
    for r in read_tsv(p):
        out[((r.get("species") or "").lower(), r.get("upstream_label") or "")] = r
    return out


def final_label(recon: Dict[Tuple[str, str], Dict[str, str]], species: str,
                upstream_label: str) -> str:
    r = recon.get(((species or "").lower(), upstream_label or ""))
    fl = (r or {}).get("final_isoform_label") or ""
    return fl or upstream_label


def label_gate(base: Path, controls: Iterable[str] = ("homo_sapiens", "mus_musculus")
               ) -> Tuple[bool, List[str]]:
    p = module_dir(base) / "maps" / "fgfr2_exon_type_label_reconciliation.tsv"
    rec = read_tsv(p)
    msgs: List[str] = []
    if not rec:
        return False, ["label reconciliation table missing; cannot validate final labels"]
    for r in rec:
        action = r.get("label_reconciliation_action", "")
        final = r.get("final_isoform_label", "")
        validated = r.get("validated_exon_type", "")
        if action in ("keep_upstream_label", "correct_final_label_from_sequence"):
            if final != validated:
                msgs.append(f"{r.get('species')}/{r.get('upstream_label')}: "
                            f"final={final} != validated={validated} (action={action})")
    for ctrl in controls:
        crows = [r for r in rec if (r.get("species") or "").lower() == ctrl]
        vals = {r.get("validated_exon_type") for r in crows}
        if not crows or vals != {"IIIb", "IIIc"}:
            msgs.append(f"control {ctrl}: validated types={sorted(vals)} (expected IIIb & IIIc)")
        for r in crows:
            if r.get("final_isoform_label") != r.get("validated_exon_type"):
                msgs.append(f"control {ctrl}/{r.get('upstream_label')}: final != validated")
    return (len(msgs) == 0), msgs


def load_claim_status(base: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    p = module_dir(base) / "maps" / "fgfr2_exon_type_label_reconciliation.tsv"
    out: Dict[Tuple[str, str], Dict[str, str]] = {}
    for r in read_tsv(p):
        sp = (r.get("species") or "").lower()
        fl = r.get("final_isoform_label") or r.get("upstream_label") or ""
        out[(sp, fl)] = r
    return out


def claim_is_primary(claim: str) -> bool:
    return (claim or "").startswith("primary_claim")


def claim_value(r: Dict[str, str]) -> str:
    return r.get("final_claim_status_after_rescue") or r.get("final_claim_status") or ""


def species_claim(claims: Dict[Tuple[str, str], Dict[str, str]], species: str) -> str:
    rank = {"excluded_from_primary_claim": 3, "supplement_review": 2,
            "primary_claim_supported_with_minor_flags": 1, "primary_claim_supported": 0}
    vals = [claim_value(r) for (sp, _), r in claims.items()
            if sp == (species or "").lower()]
    return max(vals, key=lambda v: rank.get(v, -1), default="")


def maximal_rescue_gate(base: Path) -> Tuple[bool, List[str]]:
    p = module_dir(base) / "maps" / "fgfr2_maximal_rescue_validation_gate.tsv"
    rows = read_tsv(p)
    if not rows:
        return True, ["maximal rescue gate table absent (maximal rescue not run yet)"]
    fails = [f"{r.get('check')} [{r.get('scope')}]: {r.get('detail')}"
             for r in rows if r.get("status") != "pass"]
    return (len(fails) == 0), fails


def load_post_rescue_truth(base: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    p = module_dir(base) / "maps" / "fgfr2_post_rescue_final_truth_table.tsv"
    out: Dict[Tuple[str, str], Dict[str, str]] = {}
    for r in read_tsv(p):
        out[((r.get("species") or "").lower(), r.get("isoform")
             or r.get("final_isoform_label") or "")] = r
    return out




def synteny_gate(base: Path) -> Tuple[bool, List[str]]:
    p = module_dir(base) / "synteny" / "fgfr2_5neighbor_synteny_validation_gate.tsv"
    rows = read_tsv(p)
    if not rows:
        return True, ["synteny gate table absent (synteny not run yet)"]
    fails = [f"{r.get('check')} [{r.get('scope')}]: {r.get('detail')}"
             for r in rows if r.get("status") != "pass"]
    return (len(fails) == 0), fails


def post_rescue_consistency_gate(base: Path, write: bool = True) -> Tuple[bool, List[str]]:
    md = module_dir(base)
    maps, tabd, robd = md / "maps", md / "tables", md / "robustness"
    truth = load_post_rescue_truth(base)
    checks: List[Dict[str, str]] = []

    def add(check, scope, ok, detail=""):
        checks.append({"check": check, "scope": scope,
                       "status": "pass" if ok else "FAIL", "detail": detail})

    if not truth:
        add("post_rescue_truth_table_present", "global", False, "truth table absent")
        if write:
            _write_consistency_gate(maps, checks)
        return False, ["post-rescue truth table absent"]

    def keyed(rows):
        return {((r.get("species") or "").lower(),
                 r.get("isoform") or r.get("final_isoform_label") or ""): r for r in rows}

    recon = {((r.get("species") or "").lower(), r.get("final_isoform_label")
              or r.get("upstream_label") or ""): r
             for r in read_tsv(maps / "fgfr2_exon_type_label_reconciliation.tsv")}
    rob = keyed(read_tsv(robd / "fgfr2_boundary_robustness_scores.tsv"))
    f6 = keyed(read_tsv(tabd / "figure6_msa_projected_boundary_map.tsv"))
    f8 = keyed(read_tsv(tabd / "figure8_boundary_robustness_evidence_stack.tsv"))
    f6c = keyed(read_tsv(tabd / "figure6C_human_referenced_residue_agreement_map.tsv"))
    _master = {(r.get("species") or "").lower(): r
              for r in read_tsv(require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}

    def tclaim(k):
        # Safe lookup: species-specific control checks below reference fixed reference
        # species (canis/pongo/gorilla/human/mouse) that may be ABSENT from a custom
        # run panel. Return "" for missing keys instead of raising KeyError.
        return truth[k].get("final_claim_status_after_rescue", "") if k in truth else ""

    # 1 final_isoform_label consistent (truth vs reconciliation vs robustness)
    bad = [f"{k[0]}/{k[1]}" for k in truth
           if (recon.get(k, {}).get("final_isoform_label", k[1]) != k[1]
               or (k in rob and rob[k].get("final_isoform_label", k[1]) != k[1]))]
    add("final_isoform_label_consistent", "all_tables", not bad, "; ".join(bad) or "ok")
    # 2 transcript/protein consistent (truth vs reconciliation vs robustness)
    badp = [f"{k[0]}/{k[1]}" for k in truth
            if (recon.get(k, {}).get("protein_id", truth[k]["protein_id"]) != truth[k]["protein_id"]
                or (k in rob and rob[k].get("protein_id", truth[k]["protein_id"])
                    != truth[k]["protein_id"]))]
    add("transcript_protein_consistent", "final_major_tables", not badp, "; ".join(badp) or "ok")
    # 3 claim consistent across all final figure tables
    badc = []
    for name, tbl in (("figure6", f6), ("figure8", f8), ("figure6C", f6c)):
        for k, r in tbl.items():
            if k in truth and r.get("final_claim_status_after_rescue", tclaim(k)) != tclaim(k):
                badc.append(f"{name}:{k[0]}/{k[1]}")
    add("claim_consistent_across_figure_tables", "figure_tables", not badc, "; ".join(badc) or "ok")
    # 4 no table uses upstream_label as final biology
    badu = [f"{k[0]}/{k[1]}" for k in truth
            if truth[k].get("final_label_source", "") in ("upstream_label", "upstream")]
    add("no_upstream_label_as_final_biology", "all_tables", not badu, "; ".join(badu) or "ok")
    # 5 no primary (6C) row is unresolved/excluded
    bad6c = [f"{k[0]}/{k[1]}" for k in f6c if not (tclaim(k) or "").startswith("primary_claim")]
    add("primary_figure_only_primary_rows", "figure6C", not bad6c, "; ".join(bad6c) or "ok")
    # 6 no rescued row remains supplementary_only due to stale recommended_use
    badr = [f"{k[0]}/{k[1]}" for k in truth
            if truth[k].get("rescue_decision", "").startswith("rescued")
            and (not tclaim(k).startswith("primary_claim")
                 or truth[k].get("recommended_use_post_rescue") != "main_analysis")]
    add("no_stale_supplementary_for_rescued", "rescued", not badr, "; ".join(badr) or "ok")
    # 7 figure6C/6D/8C agree on primary inclusion: 6C set == primary set in truth that have a map
    prim = {k for k in truth if tclaim(k).startswith("primary_claim")}
    badset = [f"{k[0]}/{k[1]}" for k in f6c if k not in prim]
    add("figure_tables_agree_primary_inclusion", "6C_6D_8C", not badset, "; ".join(badset) or "ok")
    # 8 human/mouse controls
    badctrl = []
    for ctrl in ("homo_sapiens", "mus_musculus"):
        ks = [k for k in truth if k[0] == ctrl]
        for k in ks:
            if not tclaim(k).startswith("primary_claim") or \
               truth[k].get("final_isoform_label") != truth[k].get("validated_exon_type"):
                badctrl.append(f"{k[0]}/{k[1]}")
    add("controls_pass", "human_mouse", not badctrl, "; ".join(badctrl) or "ok")

    # 9/10/11 species-specific consistency (across truth + robustness + figure tables)
    def claim_everywhere(sp, iso):
        vals = {tclaim((sp, iso))}
        for tbl in (rob, f6, f8, f6c):
            if (sp, iso) in tbl:
                vals.add(tbl[(sp, iso)].get("final_claim_status_after_rescue", tclaim((sp, iso))))
        return vals

    def is_primary(sp, iso):
        v = claim_everywhere(sp, iso)
        return all(x.startswith("primary_claim") for x in v) and len(v) >= 1

    # Missing reference species make species-specific controls not applicable.
    if ("gorilla_gorilla_gorilla", "IIIb") in truth or ("gorilla_gorilla_gorilla", "IIIc") in truth:
        gor_ok = (is_primary("gorilla_gorilla_gorilla", "IIIb")
                  and is_primary("gorilla_gorilla_gorilla", "IIIc"))
        add("gorilla_consistent_primary_pair", "gorilla_gorilla_gorilla", gor_ok,
            "both IIIb+IIIc primary across tables" if gor_ok else "gorilla pair not consistently primary")
    else:
        add("gorilla_consistent_primary_pair", "gorilla_gorilla_gorilla", True,
            "not_applicable: gorilla_gorilla_gorilla absent from this run panel (custom run)")
    if ("canis_lupus_familiaris", "IIIb") in truth or ("canis_lupus_familiaris", "IIIc") in truth:
        can_ok = (is_primary("canis_lupus_familiaris", "IIIb")
                  and not tclaim(("canis_lupus_familiaris", "IIIc")).startswith("primary_claim"))
        add("canis_partial_rescue_consistent", "canis_lupus_familiaris", can_ok,
            "IIIb primary, IIIc supplement/review" if can_ok else "canis partial state inconsistent")
    else:
        add("canis_partial_rescue_consistent", "canis_lupus_familiaris", True,
            "not_applicable: canis_lupus_familiaris absent from this run panel (custom run)")
    if ("pongo_abelii", "IIIb") in truth or ("pongo_abelii", "IIIc") in truth:
        pon_ok = (is_primary("pongo_abelii", "IIIc")
                  and not tclaim(("pongo_abelii", "IIIb")).startswith("primary_claim"))
        add("pongo_partial_rescue_consistent", "pongo_abelii", pon_ok,
            "IIIc primary, IIIb supplement/review" if pon_ok else "pongo partial state inconsistent")
    else:
        add("pongo_partial_rescue_consistent", "pongo_abelii", True,
            "not_applicable: pongo_abelii absent from this run panel (custom run)")

    if write:
        _write_consistency_gate(maps, checks)
    fails = [f"{c['check']} [{c['scope']}]: {c['detail']}" for c in checks if c["status"] != "pass"]
    return (len(fails) == 0), fails


def _write_consistency_gate(maps: Path, checks: List[Dict[str, str]]) -> None:
    import json as _json
    write_tsv(maps / "fgfr2_post_rescue_cross_table_consistency_gate.tsv", checks,
              ["check", "scope", "status", "detail"])
    hard = any(c["status"] != "pass" for c in checks)
    (maps / "fgfr2_post_rescue_cross_table_consistency_gate.json").write_text(
        _json.dumps({"checks": checks, "hard_fail": hard, "timestamp": now_iso()}, indent=2),
        encoding="utf-8")


def general_rescue_gate(base: Path) -> Tuple[bool, List[str]]:
    p = module_dir(base) / "maps" / "fgfr2_general_rescue_validation_gate.tsv"
    rows = read_tsv(p)
    if not rows:
        return True, ["general rescue gate table absent (validation/rescue not run yet)"]
    fails = [f"{r.get('check')} [{r.get('scope')}]: {r.get('detail')}"
             for r in rows if r.get("status") != "pass"]
    return (len(fails) == 0), fails


def clean_alignment_seq(seq: str) -> str:
    return "".join(seq.split()).upper()


def ungapped(seq: str) -> str:
    return "".join(c for c in seq if c not in GAP_CHARS)


def invalid_residues(seq: str) -> str:
    bad = set()
    for c in seq.upper():
        if c in AA_VALID or c in AA_AMBIG or c in GAP_CHARS or c == "*":
            continue
        bad.add(c)
    return "".join(sorted(bad))
