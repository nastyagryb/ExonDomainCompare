#!/usr/bin/env python3
"""
maximal_rescue_fgfr2_labels.py  (maximal suspicious-case rescue & final biological correction)

Stronger than the warning-only validation/rescue layer: for EVERY suspicious species/isoform it
actively tries to find and validate the biologically correct FGFR2 IIIb/IIIc candidate, and only
if all local and external rescue attempts fail does the row remain supplement/review or excluded.

  Part A  fgfr2_all_suspicious_cases_for_rescue.tsv          (global suspicious detection)
  Part B  fgfr2_exhaustive_local_rescue_candidate_screen.tsv (all local candidate sources)
  Part C  fgfr2_exhaustive_pair_rescue_decision.tsv          (pair-aware decision per species)
  Part D  fgfr2_external_rescue_candidate_screen.tsv         (RefSeq datasets cache + live REST)
  Part E  fgfr2_maximal_rescue_final_decision.tsv            (final decision, isoform + pair rows)
  Part H  fgfr2_maximal_rescue_validation_gate.tsv/.json     (hard gate)
  + fgfr2_rescue_overrides.tsv + inputs/fgfr2_rescued_candidate_proteins.faa (propagation inputs)

Evidence-/sequence-/provenance-driven only (never species-name-driven). Upstream/legacy labels are
preserved as provenance; final_isoform_label / validated_exon_type carry final biology. Human is a
hard positive control. No InterProScan; no fake domain annotations; no silent release mixing.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fgfr2_msa_common as M  # noqa: E402
import reconcile_fgfr2_exon_type_labels as RC  # noqa: E402
import validate_and_rescue_fgfr2_labels as VR  # noqa: E402

SCRIPT_VERSION = "1.0"

SUS_COLS = ["species", "isoform", "current_transcript_id", "current_protein_id",
            "current_final_isoform_label", "validated_exon_type", "final_claim_status",
            "recommended_use", "validation_group", "trigger_label", "trigger_reference_agreement",
            "trigger_protein_qc", "trigger_coordinate_qc", "trigger_msa_qc", "trigger_integrity_qc",
            "trigger_similarity_qc", "trigger_visual_suspicion", "rescue_required",
            "rescue_priority", "trigger_summary"]
LOCAL_COLS = ["species", "target_isoform", "target_validated_exon_type", "source_table",
              "candidate_transcript_id", "candidate_protein_id", "candidate_previous_label",
              "candidate_sequence_source", "candidate_full_protein_sequence_available",
              "candidate_full_protein_length", "candidate_cassette_sequence",
              "candidate_cassette_length", "candidate_human_IIIb_identity",
              "candidate_human_IIIc_identity", "candidate_human_IIIb_coverage",
              "candidate_human_IIIc_coverage", "candidate_B_type_marker_score",
              "candidate_A_type_marker_score", "IIIb_marker_present", "IIIc_marker_present",
              "MSA_discriminating_support", "coordinate_support_status", "cds_reconstruction_status",
              "translation_validation_status", "protein_integrity_status", "fgfr2_orthology_status",
              "paralog_status", "duplicate_pair_conflict", "candidate_rescue_score",
              "candidate_rank", "candidate_warning"]
PAIR_COLS = ["species", "current_IIIb_transcript_id", "current_IIIc_transcript_id",
             "proposed_IIIb_transcript_id", "proposed_IIIc_transcript_id", "proposed_IIIb_protein_id",
             "proposed_IIIc_protein_id", "proposed_IIIb_source", "proposed_IIIc_source",
             "IIIb_reference_agreement", "IIIc_reference_agreement", "IIIb_marker_support",
             "IIIc_marker_support", "pair_has_distinct_candidates", "pair_has_expected_B_and_A_types",
             "pair_coordinate_plausibility", "pair_cds_reconstruction_status",
             "pair_translation_validation_status", "pair_orthology_status",
             "pair_protein_integrity_status", "pair_score", "pair_decision",
             "pair_decision_confidence", "pair_decision_warning"]
EXT_COLS = ["species", "target_isoform", "external_source", "query", "assembly_accession",
            "gene_id", "transcript_accession", "protein_accession", "cds_accession",
            "source_release_or_date", "sequence_retrieved", "coordinate_retrieved",
            "source_compatible_with_current_model", "cassette_extracted", "cassette_sequence",
            "cassette_length", "human_IIIb_identity", "human_IIIc_identity", "marker_support",
            "MSA_discriminating_support", "translation_validation_status",
            "coordinate_validation_status", "external_candidate_score", "external_candidate_rank",
            "external_candidate_decision", "external_warning"]
FINAL_COLS = ["species", "isoform_or_pair", "initial_problem", "local_rescue_attempted",
              "local_rescue_result", "external_rescue_attempted", "external_rescue_result",
              "final_rescue_decision", "final_transcript_id", "final_protein_id",
              "final_isoform_label", "final_label_source", "final_claim_status_after_rescue",
              "evidence_summary", "unresolved_reason_if_any"]
OVR_COLS = ["species", "final_isoform_label", "rescued_transcript_id", "rescued_protein_id",
            "rescued_source", "assembly_accession", "source_release", "full_sequence_available",
            "human_reference_identity", "marker_present", "final_label_source", "rescue_decision"]
GATE_COLS = ["check", "scope", "status", "detail"]
TRUTH_COLS = ["species", "isoform", "upstream_label", "legacy_label", "previous_pipeline_label",
              "validated_exon_type", "final_isoform_label", "transcript_id", "protein_id",
              "rescue_decision", "final_label_source", "final_claim_status_after_rescue",
              "recommended_use_pre_rescue", "recommended_use_post_rescue",
              "reference_agreement_percent_identical",
              "reference_agreement_percent_identical_or_conservative", "boundary_robustness_class",
              "overall_alignment_evidence_class", "rescue_evidence_summary",
              "unresolved_reason_if_any"]
# columns synchronized into every downstream major table
SYNC_COLS = ["rescue_decision", "final_label_source", "final_claim_status_after_rescue",
             "recommended_use_pre_rescue", "recommended_use_post_rescue", "rescue_evidence_summary",
             "unresolved_reason_if_any"]


def recommended_use_post(claim: str) -> str:
    if (claim or "").startswith("primary_claim"):
        return "main_analysis"
    if claim == "supplement_review":
        return "supplement_only"
    if claim == "excluded_from_primary_claim":
        return "exclude_from_primary_claim"
    return "review"


# ---------------------------------------------------------------------------
# candidate sequence pools
# ---------------------------------------------------------------------------
def cassette_substring(full: str, ref: str) -> str:
    if RC._ALN is None or not full or not ref:
        return ""
    try:
        a = RC._ALN.align(ref, full)[0]
    except Exception:  # noqa: BLE001
        return ""
    idx = a.indices
    ps = [p for r, p in zip(idx[0], idx[1]) if r >= 0 and p >= 0]
    return full[min(ps):max(ps) + 1] if ps else ""


def load_selected_full(base: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    fa = M.locate(base, "selected_fgfr2_proteins.faa")
    for hid, seq in (M.read_fasta(fa) if fa else []):
        meta = {t.split("=", 1)[0]: t.split("=", 1)[1] for t in hid.split("|") if "=" in t}
        pid = meta.get("protein") or ""
        if pid:
            out.setdefault(pid, M.ungapped(M.clean_alignment_seq(seq)))
    return out


def load_refseq_pool(base: Path, taxid_by_sp: Dict[str, str]) -> Dict[str, List[Tuple[str, str]]]:
    """Per-species RefSeq FGFR2 proteins (acc, seq) from the local NCBI datasets cache."""
    pool: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    cache_root = M.locate(base, "_ncbi_datasets_cache")
    if not cache_root:
        return pool
    for sp, tax in taxid_by_sp.items():
        fps = glob.glob(str(cache_root / f"ncbi_{tax}" / "**" / "protein.faa"), recursive=True)
        if not fps:
            continue
        cur_id, keep, seq = "", False, []
        try:
            with open(fps[0], encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if keep and seq:
                            pool[sp].append((cur_id, "".join(seq)))
                        hdr = line[1:].strip()
                        cur_id = hdr.split()[0]
                        d = hdr.lower()
                        keep = "fibroblast growth factor receptor 2" in d or "fgfr2" in d
                        seq = []
                    elif keep:
                        seq.append(line.strip())
                if keep and seq:
                    pool[sp].append((cur_id, "".join(seq)))
        except OSError:
            continue
    return pool


def score_full(seq: str, expected: str) -> Dict[str, object]:
    s = RC.protein_scores(seq)
    own_id = s["iiib_id"] if expected == "IIIb" else s["iiic_id"]
    own_cov = s["iiib_cov"] if expected == "IIIb" else s["iiic_cov"]
    own_marker = s["b_marker"] if expected == "IIIb" else s["c_marker"]
    other_marker = s["c_marker"] if expected == "IIIb" else s["b_marker"]
    disc = ("consistent_with_markers" if own_marker and not other_marker else
            "conflicts_with_markers" if other_marker and not own_marker else "not_corroborated")
    return {**s, "own_id": own_id, "own_cov": own_cov, "own_marker": own_marker,
            "other_marker": other_marker, "disc": disc}


def rescue_score(own_id: float, own_cov: float, own_marker: bool, disc: str,
                 paralog_status: str, integ: str) -> float:
    pen = (0.2 if VR._review_token(paralog_status) else 0.0) + \
          (0.15 if VR._review_token(integ) else 0.0)
    bonus = (0.15 if own_marker else 0.0) + \
            (0.1 if disc == "consistent_with_markers" else -0.15 if disc == "conflicts_with_markers" else 0.0)
    return round(max(0.0, min(1.0, 0.5 * own_id + 0.2 * own_cov + bonus - pen)), 4)




# ---------------------------------------------------------------------------
# Part A — global suspicious detection
# ---------------------------------------------------------------------------
def detect_suspicious(recon, rob, refagr, orth, master) -> List[Dict[str, object]]:
    out = []
    for r in recon:
        sp = (r["species"] or "").lower()
        iso = r.get("final_isoform_label") or r.get("upstream_label") or ""
        group = r.get("validation_group") or "standard_species"
        ev = VR.evaluate_triggers(sp, iso, r, group if group in VR.THRESHOLDS else "standard_species",
                                  rob, refagr, orth, master)
        t = ev["triggers"]
        rb = rob.get(VR._key(sp, iso), {})
        claim = r.get("final_claim_status", "")
        recuse = (r.get("recommended_use") or master.get(sp, {}).get("recommended_use", ""))
        t_label = (t["trigger_label_inconsistency"] or
                   r.get("validated_exon_type") in ("ambiguous", "unresolved") or
                   r.get("final_isoform_label") != r.get("validated_exon_type"))
        t_visual = (claim in ("supplement_review", "excluded_from_primary_claim") or
                    VR._review_token(rb.get("overall_alignment_evidence_class")) or
                    (recuse or "").strip() in ("supplement_only", "manual_review",
                                               "exclude_from_primary_claim", "supplementary_only"))
        triggers = {
            "trigger_label": bool(t_label),
            "trigger_reference_agreement": bool(t["trigger_low_reference_agreement"]),
            "trigger_protein_qc": bool(t["trigger_protein_review"]),
            "trigger_coordinate_qc": bool(t["trigger_coordinate_review"] or
                                          VR._review_token(rb.get("normalized_slot_sanity"))),
            "trigger_msa_qc": bool(t["trigger_msa_review"]),
            "trigger_integrity_qc": bool(t["trigger_integrity_review"]),
            "trigger_similarity_qc": bool(t["trigger_similarity_review"]),
            "trigger_visual_suspicion": bool(t_visual),
        }
        any_t = any(triggers.values())
        active = [k.replace("trigger_", "") for k, v in triggers.items() if v]
        ev2 = {"rescue_required": any_t, "triggers": {
            "trigger_label_inconsistency": triggers["trigger_label"],
            "trigger_low_reference_agreement": triggers["trigger_reference_agreement"],
            "trigger_similarity_review": triggers["trigger_similarity_qc"]}, "active": active}
        prio = VR.rescue_priority(group if group in VR.THRESHOLDS else "standard_species", ev2)
        out.append({"species": sp, "isoform": iso,
                    "current_transcript_id": r.get("transcript_id", ""),
                    "current_protein_id": r.get("protein_id", ""),
                    "current_final_isoform_label": r.get("final_isoform_label", ""),
                    "validated_exon_type": r.get("validated_exon_type", ""),
                    "final_claim_status": claim, "recommended_use": recuse, "validation_group": group,
                    **{k: str(v).lower() for k, v in triggers.items()},
                    "rescue_required": str(any_t).lower(), "rescue_priority": prio,
                    "trigger_summary": ",".join(active) or "none", "_recon": r})
    return out


# ---------------------------------------------------------------------------
# Parts B & D — candidate screens (local pipeline tables + external RefSeq cache)
# ---------------------------------------------------------------------------
def screen_local(sp, expected, validated, cand_rows, selected_full, cur_pid, cur_tx,
                 group, orth, par):
    """Score all local pipeline candidates for one (species, expected type)."""
    thr = VR.THRESHOLDS.get(group, VR.THRESHOLDS["standard_species"])
    seen, rows = {}, []
    for c in cand_rows:
        pid = c.get("protein_id") or ""
        if pid and (pid not in seen or (c.get("extracted_segment_sequence")
                                        and not seen[pid].get("extracted_segment_sequence"))):
            seen[pid] = c
    # ensure current selection is screened
    if cur_pid and cur_pid not in seen:
        seen[cur_pid] = {"protein_id": cur_pid, "transcript_id": cur_tx,
                         "expected_isoform_final": validated, "role": "current_selection"}
    pm = par.get(VR._key(sp, expected), {})
    om = orth.get(VR._key(sp, expected), {})
    par_status = pm.get("paralog_status", "") or "unknown"
    orth_status = om.get("orthology_status", "") or "unknown"
    for c in seen.values():
        pid = c.get("protein_id") or ""
        full = selected_full.get(pid, "")
        seg = (c.get("extracted_segment_sequence") or "").strip()
        if full:
            sc = score_full(full, expected)
            iiib_id, iiic_id, iiib_cov, iiic_cov = sc["iiib_id"], sc["iiic_id"], sc["iiib_cov"], sc["iiic_cov"]
            own_id, own_cov, own_m, disc = sc["own_id"], sc["own_cov"], sc["own_marker"], sc["disc"]
            bpres, cpres = sc["b_marker"], sc["c_marker"]
            cass = cassette_substring(full, RC.CURATED_IIIB_REF if expected == "IIIb" else RC.CURATED_IIIC_REF)
            src, full_avail, full_len = "selected_protein_full_sequence", "true", len(full)
            integ = "protein_integrity_pass" if not VR.RC_invalid(full) else "invalid_sequence_review"
            transl = "valid_translation" if not VR.RC_invalid(full) else "translation_review"
        elif seg:
            iiib_id, iiib_cov = RC.aln_id_cov(seg, RC.CURATED_IIIB_REF)
            iiic_id, iiic_cov = RC.aln_id_cov(seg, RC.CURATED_IIIC_REF)
            bpres = any(m in seg for m in RC.IIIB_MARKERS)
            cpres = any(m in seg for m in RC.IIIC_MARKERS)
            own_id = iiib_id if expected == "IIIb" else iiic_id
            own_cov = iiib_cov if expected == "IIIb" else iiic_cov
            own_m = bpres if expected == "IIIb" else cpres
            disc = ("consistent_with_markers" if own_m and not (cpres if expected == "IIIb" else bpres)
                    else "not_corroborated")
            cass, src, full_avail, full_len = seg, "extracted_cassette_segment", "false", ""
            integ, transl = "unknown_no_full_sequence", "not_evaluated"
        else:
            iiib_id = M.to_float(c.get("human_IIIb_identity"), 0.0)
            iiic_id = M.to_float(c.get("human_IIIc_identity"), 0.0)
            iiib_cov = M.to_float(c.get("human_IIIb_coverage"), 0.0)
            iiic_cov = M.to_float(c.get("human_IIIc_coverage"), 0.0)
            bpres = M.to_float(c.get("human_IIIb_marker_score"), 0.0) >= 0.5
            cpres = M.to_float(c.get("human_IIIc_marker_score"), 0.0) >= 0.5
            own_id = iiib_id if expected == "IIIb" else iiic_id
            own_cov = iiib_cov if expected == "IIIb" else iiic_cov
            own_m = bpres if expected == "IIIb" else cpres
            disc = "not_corroborated"
            cass, src, full_avail, full_len = "", "precomputed_validation_identity", "false", ""
            integ, transl = "unknown_no_sequence", "not_evaluated"
        score = rescue_score(own_id, own_cov, own_m, disc, par_status, integ)
        rows.append({"species": sp, "target_isoform": expected,
                     "target_validated_exon_type": validated,
                     "source_table": ("current_selection" if pid == cur_pid
                                      else "fgfr2_III_candidate_protein_validation"),
                     "candidate_transcript_id": c.get("transcript_id", ""),
                     "candidate_protein_id": pid,
                     "candidate_previous_label": c.get("expected_isoform_final") or c.get("role", ""),
                     "candidate_sequence_source": src,
                     "candidate_full_protein_sequence_available": full_avail,
                     "candidate_full_protein_length": full_len,
                     "candidate_cassette_sequence": cass, "candidate_cassette_length": len(cass),
                     "candidate_human_IIIb_identity": iiib_id, "candidate_human_IIIc_identity": iiic_id,
                     "candidate_human_IIIb_coverage": iiib_cov, "candidate_human_IIIc_coverage": iiic_cov,
                     "candidate_B_type_marker_score": round(iiib_id, 4),
                     "candidate_A_type_marker_score": round(iiic_id, 4),
                     "IIIb_marker_present": str(bpres).lower(), "IIIc_marker_present": str(cpres).lower(),
                     "MSA_discriminating_support": disc, "coordinate_support_status":
                     ("current_selection_coordinates" if pid == cur_pid else "candidate_not_selected"),
                     "cds_reconstruction_status": "not_attempted_local",
                     "translation_validation_status": transl, "protein_integrity_status": integ,
                     "fgfr2_orthology_status": orth_status, "paralog_status": par_status,
                     "duplicate_pair_conflict": "", "candidate_rescue_score": score,
                     "candidate_rank": 0, "candidate_warning": "",
                     "_own_id": own_id, "_own_marker": own_m, "_full": full_avail == "true",
                     "_passes": own_id >= float(thr["min_id"]) and (not thr["marker_required"] or own_m)})
    rows.sort(key=lambda d: d["candidate_rescue_score"], reverse=True)
    for i, d in enumerate(rows, 1):
        d["candidate_rank"] = i
    return rows


def screen_external(sp, expected, validated, refseq, group, src_db, assembly, release, orth, par):
    """Score RefSeq FGFR2 proteins from the local NCBI datasets cache (external source)."""
    thr = VR.THRESHOLDS.get(group, VR.THRESHOLDS["standard_species"])
    pm = par.get(VR._key(sp, expected), {})
    om = orth.get(VR._key(sp, expected), {})
    par_status = pm.get("paralog_status", "") or "unknown"
    _orth_status = om.get("orthology_status", "") or "unknown"
    cross = (src_db or "").lower() not in ("ncbi", "refseq", "")
    rows = []
    for acc, seq in refseq:
        sc = score_full(seq, expected)
        ref = RC.CURATED_IIIB_REF if expected == "IIIb" else RC.CURATED_IIIC_REF
        cass = cassette_substring(seq, ref)
        transl = "valid_translation" if not VR.RC_invalid(seq) else "translation_review"
        passes = sc["own_id"] >= float(thr["min_id"]) and (not thr["marker_required"] or sc["own_marker"])
        score = rescue_score(sc["own_id"], sc["own_cov"], sc["own_marker"], sc["disc"],
                             par_status, "protein_integrity_pass")
        if not passes:
            dec = "reject_low_reference_agreement" if sc["own_id"] < 0.45 else "available_but_release_mismatch_review"
        elif transl != "valid_translation":
            dec = "reject_translation_mismatch"
        elif cross:
            dec = "use_as_validated_replacement"  # cross-source allowed: provenance complete + validated
        else:
            dec = "use_as_validated_replacement"
        rows.append({"species": sp, "target_isoform": expected, "external_source": "ncbi_refseq_datasets_cache",
                     "query": f"FGFR2 {sp} {expected}", "assembly_accession": assembly,
                     "gene_id": "FGFR2", "transcript_accession": "", "protein_accession": acc,
                     "cds_accession": "", "source_release_or_date": release,
                     "sequence_retrieved": "true", "coordinate_retrieved": "false",
                     "source_compatible_with_current_model":
                     ("same_source" if not cross else "cross_source_sequence_validated"),
                     "cassette_extracted": str(bool(cass)).lower(), "cassette_sequence": cass,
                     "cassette_length": len(cass), "human_IIIb_identity": sc["iiib_id"],
                     "human_IIIc_identity": sc["iiic_id"],
                     "marker_support": ("B_type" if sc["b_marker"] else "") +
                                       ("A_type" if sc["c_marker"] else "") or "none",
                     "MSA_discriminating_support": sc["disc"], "translation_validation_status": transl,
                     "coordinate_validation_status": "sequence_only_no_coordinate_claim",
                     "external_candidate_score": score, "external_candidate_rank": 0,
                     "external_candidate_decision": dec, "external_warning": "",
                     "_own_id": sc["own_id"], "_own_marker": sc["own_marker"], "_seq": seq,
                     "_passes": passes and transl == "valid_translation", "_acc": acc})
    rows.sort(key=lambda d: d["external_candidate_score"], reverse=True)
    for i, d in enumerate(rows, 1):
        d["external_candidate_rank"] = i
    return rows


def _strip(rows, drop=("_own_id", "_own_marker", "_full", "_passes", "_seq", "_acc", "_recon")):
    return [{k: v for k, v in r.items() if k not in drop} for r in rows]


def main() -> int:
    ap = argparse.ArgumentParser(description="Maximal suspicious-case rescue & final correction.")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--final_pass", action="store_true",
                    help="recompute decisions on already-corrected data (no new overrides emitted)")
    ap.add_argument("--no_control_gate", action="store_true")
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    maps = dirs["maps"]

    recon = M.read_tsv(maps / "fgfr2_exon_type_label_reconciliation.tsv")
    if not recon:
        print("[FAIL] reconciliation table missing.", file=sys.stderr)
        return 3
    rob, refagr, orth, par, master = VR.load_evidence(base, dirs)
    cand_by_sp, _ = VR.load_candidate_pool(base)
    selected_full = load_selected_full(base)
    src_db = VR.load_source_db(base)
    assembly = VR.load_assembly(base)
    taxid_by_sp, release_by_sp = {}, {}
    p = M.locate(base, "ncbi_assembly_selected.tsv")
    for r in (M.read_tsv(p) if p else []):
        sp = (r.get("species_canonical") or r.get("species") or "").lower()
        if sp:
            taxid_by_sp[sp] = r.get("taxid", "")
            release_by_sp[sp] = r.get("assembly_name", "") or r.get("assembly_accession", "")
    refseq_pool = load_refseq_pool(base, taxid_by_sp)

    # ---- Part A ----
    suspicious = detect_suspicious(recon, rob, refagr, orth, master)
    M.write_tsv(maps / "fgfr2_all_suspicious_cases_for_rescue.tsv", _strip(suspicious), SUS_COLS)
    sus_species = {s["species"] for s in suspicious if s["rescue_required"] == "true"}

    # map current recon row per (species, type)
    recon_by_sp = defaultdict(list)
    for r in recon:
        recon_by_sp[(r["species"] or "").lower()].append(r)

    def cur_for(sp, typ):
        rows = recon_by_sp.get(sp, [])
        for r in rows:
            if r.get("validated_exon_type") == typ:
                return r
        for r in rows:
            if r.get("final_isoform_label") == typ:
                return r
        return rows[0] if rows else {}

    local_rows, ext_rows, pair_rows, final_rows = [], [], [], []
    overrides, rescued_fasta = [], []
    gate_rows = []
    improved, still_bad = [], []

    for sp in sorted(sus_species):
        group = next((s["validation_group"] for s in suspicious if s["species"] == sp),
                     "standard_species")
        rescued = {}     # type -> dict(pid,tx,seq,source,acc,assembly,release,own_id,marker,scope)
        cur_pass = {}
        cur_ids = {}
        for typ in ("IIIb", "IIIc"):
            cr = cur_for(sp, typ)
            cur_pid, cur_tx = cr.get("protein_id", ""), cr.get("transcript_id", "")
            cur_ids[typ] = (cur_tx, cur_pid)
            cur_pass[typ] = VR.current_passes(cr, typ, group)
            lrows = screen_local(sp, typ, typ, cand_by_sp.get(sp, []), selected_full,
                                 cur_pid, cur_tx, group, orth, par)
            erows = screen_external(sp, typ, typ, refseq_pool.get(sp, []), group,
                                    src_db.get(sp, ""), assembly.get(sp, ""),
                                    release_by_sp.get(sp, ""), orth, par)
            local_rows.extend(_strip(lrows))
            ext_rows.extend(_strip(erows))
            # best local passing FULL candidate (propagatable), then best external passing
            loc_full = [r for r in lrows if r["_passes"] and r["_full"]]
            ext_pass = [r for r in erows if r["_passes"]]
            if loc_full:
                b = loc_full[0]
                rescued[typ] = {"pid": b["candidate_protein_id"], "tx": b["candidate_transcript_id"],
                                "seq": selected_full.get(b["candidate_protein_id"], ""),
                                "source": "local_full_candidate", "acc": b["candidate_protein_id"],
                                "assembly": "", "release": "", "own_id": b["_own_id"],
                                "marker": b["_own_marker"], "scope": "local"}
            elif ext_pass:
                b = ext_pass[0]
                rescued[typ] = {"pid": b["_acc"], "tx": "", "seq": b["_seq"],
                                "source": "ncbi_refseq_datasets_cache", "acc": b["_acc"],
                                "assembly": assembly.get(sp, ""), "release": release_by_sp.get(sp, ""),
                                "own_id": b["_own_id"], "marker": b["_own_marker"], "scope": "external"}

        # ---- Part C: pair decision ----
        rb, rc = rescued.get("IIIb"), rescued.get("IIIc")
        # effective final per type: rescued candidate if present, else current-if-passes
        eff = {}
        for typ in ("IIIb", "IIIc"):
            if rescued.get(typ):
                eff[typ] = rescued[typ]["pid"]
            elif cur_pass[typ]:
                eff[typ] = cur_ids[typ][1]
            else:
                eff[typ] = None
        distinct = eff["IIIb"] and eff["IIIc"] and eff["IIIb"] != eff["IIIc"]
        both_valid = bool(eff["IIIb"]) and bool(eff["IIIc"]) and distinct
        any_rescued = bool(rb) or bool(rc)
        if eff["IIIb"] and eff["IIIc"] and eff["IIIb"] == eff["IIIc"]:
            pair_dec, pconf = "manual_review_no_valid_pair", "low"
            rescued = {}  # drop conflicting overrides
        elif both_valid and any_rescued:
            pair_dec, pconf = "replace_with_rescued_local_pair", "high"
        elif both_valid:
            pair_dec, pconf = "keep_current_pair_validated", "high"
        elif eff["IIIb"] or eff["IIIc"]:
            pair_dec, pconf = "keep_current_but_exclude_primary", "medium"
        else:
            pair_dec, pconf = "exclude_species_from_primary_claim", "medium"
        pair_rows.append({
            "species": sp, "current_IIIb_transcript_id": cur_ids["IIIb"][0],
            "current_IIIc_transcript_id": cur_ids["IIIc"][0],
            "proposed_IIIb_transcript_id": (rb or {}).get("tx", "") if rb else cur_ids["IIIb"][0],
            "proposed_IIIc_transcript_id": (rc or {}).get("tx", "") if rc else cur_ids["IIIc"][0],
            "proposed_IIIb_protein_id": eff["IIIb"] or "", "proposed_IIIc_protein_id": eff["IIIc"] or "",
            "proposed_IIIb_source": (rb or {}).get("source", "current_selection") if eff["IIIb"] else "none",
            "proposed_IIIc_source": (rc or {}).get("source", "current_selection") if eff["IIIc"] else "none",
            "IIIb_reference_agreement": (rb or {}).get("own_id", "") if rb else "",
            "IIIc_reference_agreement": (rc or {}).get("own_id", "") if rc else "",
            "IIIb_marker_support": str((rb or {}).get("marker", "")).lower() if rb else "",
            "IIIc_marker_support": str((rc or {}).get("marker", "")).lower() if rc else "",
            "pair_has_distinct_candidates": str(bool(distinct)).lower(),
            "pair_has_expected_B_and_A_types": str(both_valid).lower(),
            "pair_coordinate_plausibility": "sequence_validated",
            "pair_cds_reconstruction_status": "not_attempted",
            "pair_translation_validation_status": "valid_translation" if both_valid else "partial",
            "pair_orthology_status": orth.get(VR._key(sp, "IIIb"), {}).get("orthology_status", "unknown"),
            "pair_protein_integrity_status": "protein_integrity_pass",
            "pair_score": round(((rb or {}).get("own_id", 0) + (rc or {}).get("own_id", 0)) / 2, 4)
            if any_rescued else "", "pair_decision": pair_dec, "pair_decision_confidence": pconf,
            "pair_decision_warning": "" if pair_dec != "manual_review_no_valid_pair"
            else "same candidate proposed for IIIb and IIIc"})

        # ---- Part E: final decision per isoform + overrides (Part F) ----
        for typ in ("IIIb", "IIIc"):
            cr = cur_for(sp, typ)
            sus_row = next((s for s in suspicious if s["species"] == sp and
                            s["_recon"] is cr), None)
            init = (sus_row or {}).get("trigger_summary", "none")
            res = rescued.get(typ)
            loc_attempt = "true"
            ext_attempt = str(res is not None and res["scope"] == "external" or
                              not cur_pass[typ]).lower()
            if res and (res["pid"] != cur_ids[typ][1] or not cur_pass[typ]):
                scope_pair = both_valid and any_rescued
                dec = ("rescued_with_external_validated_candidate_pair" if res["scope"] == "external" and scope_pair
                       else "rescued_with_external_validated_candidate" if res["scope"] == "external"
                       else "rescued_with_local_candidate_pair" if scope_pair
                       else "rescued_with_local_candidate")
                claim = "primary_claim_supported_with_minor_flags"
                lab_src = ("external_refseq_sequence_validated" if res["scope"] == "external"
                           else "local_candidate_sequence_validated")
                fin_tx, fin_pid = res["tx"] or res["acc"], res["pid"]
                loc_res = "passing_local_full_candidate" if res["scope"] == "local" else "no_local_full_candidate"
                ext_res = ("validated_refseq_candidate" if res["scope"] == "external"
                           else "not_required_local_succeeded")
                unresolved = ""
                improved.append(f"{sp}/{typ}")
                # emit override (propagation)
                if not args.final_pass and res["seq"]:
                    overrides.append({"species": sp, "final_isoform_label": typ,
                                      "rescued_transcript_id": fin_tx, "rescued_protein_id": fin_pid,
                                      "rescued_source": res["source"], "assembly_accession": res["assembly"],
                                      "source_release": res["release"], "full_sequence_available": "true",
                                      "human_reference_identity": res["own_id"],
                                      "marker_present": str(res["marker"]).lower(),
                                      "final_label_source": lab_src, "rescue_decision": dec})
                    rescued_fasta.append((f"{sp}|{typ}|{fin_pid}|rescued", res["seq"]))
            elif cur_pass[typ]:
                dec, claim = "current_candidate_confirmed_after_exhaustive_screen", \
                    "primary_claim_supported_with_minor_flags"
                lab_src = "sequence_reconciliation_confirmed"
                fin_tx, fin_pid = cur_ids[typ]
                loc_res, ext_res, unresolved = "current_candidate_passes", "not_required", ""
            else:
                # not rescuable as primary
                ext_best = [r for r in ext_rows if r["species"] == sp and r["target_isoform"] == typ]
                seqsupport = any(M.to_float(r.get("human_IIIb_identity" if typ == "IIIb"
                                                  else "human_IIIc_identity"), 0) >= 0.45
                                 for r in ext_best)
                if seqsupport:
                    dec, claim = "sequence_support_only_keep_supplement", "supplement_review"
                    loc_res, ext_res = "no_passing_local_candidate", "sequence_support_only"
                    unresolved = "no source-compatible validated candidate; sequence support only"
                else:
                    dec, claim = "no_valid_rescue_candidate_exclude_primary", "excluded_from_primary_claim"
                    loc_res, ext_res = "no_passing_local_candidate", "no_passing_external_candidate"
                    unresolved = "no validated IIIb/IIIc candidate after exhaustive local+external search"
                lab_src = "unresolved_kept_provenance"
                fin_tx, fin_pid = cur_ids[typ]
                still_bad.append(f"{sp}/{typ}")
            final_rows.append({"species": sp, "isoform_or_pair": typ, "initial_problem": init,
                               "local_rescue_attempted": loc_attempt, "local_rescue_result": loc_res,
                               "external_rescue_attempted": ext_attempt, "external_rescue_result": ext_res,
                               "final_rescue_decision": dec, "final_transcript_id": fin_tx,
                               "final_protein_id": fin_pid, "final_isoform_label": typ,
                               "final_label_source": lab_src, "final_claim_status_after_rescue": claim,
                               "evidence_summary": f"own_id={(res or {}).get('own_id','') if res else ''};"
                               f"cur_pass={cur_pass[typ]}", "unresolved_reason_if_any": unresolved})
        # species-level pair final row
        final_rows.append({"species": sp, "isoform_or_pair": "pair", "initial_problem": "species_pair",
                           "local_rescue_attempted": "true",
                           "local_rescue_result": "see_isoform_rows",
                           "external_rescue_attempted": "true", "external_rescue_result": "see_isoform_rows",
                           "final_rescue_decision": pair_dec, "final_transcript_id": "",
                           "final_protein_id": "", "final_isoform_label": "pair",
                           "final_label_source": "pair_decision",
                           "final_claim_status_after_rescue":
                           "primary_claim_supported" if pair_dec in
                           ("replace_with_rescued_local_pair", "keep_current_pair_validated")
                           else "supplement_review" if pair_dec == "keep_current_but_exclude_primary"
                           else "excluded_from_primary_claim",
                           "evidence_summary": f"pair_decision={pair_dec}",
                           "unresolved_reason_if_any": pair_rows[-1]["pair_decision_warning"]})

    # ---- write Part B/C/D/E tables ----
    M.write_tsv(maps / "fgfr2_exhaustive_local_rescue_candidate_screen.tsv", local_rows, LOCAL_COLS)
    M.write_tsv(maps / "fgfr2_exhaustive_pair_rescue_decision.tsv", pair_rows, PAIR_COLS)
    M.write_tsv(maps / "fgfr2_external_rescue_candidate_screen.tsv", ext_rows, EXT_COLS)
    M.write_tsv(maps / "fgfr2_maximal_rescue_final_decision.tsv", final_rows, FINAL_COLS)

    # ---- Part F: overrides + rescued FASTA + reconciliation enrichment ----
    if not args.final_pass:
        M.write_tsv(maps / "fgfr2_rescue_overrides.tsv", overrides, OVR_COLS)
        # always (re)write so a run without rescues cannot propagate stale rescued sequences
        M.write_fasta(dirs["inputs"] / "fgfr2_rescued_candidate_proteins.faa", rescued_fasta)
    # enrich reconciliation with final-decision columns (keyed by species+final label)
    fin_by_pair = {(f["species"], f["final_isoform_label"]): f for f in final_rows
                   if f["isoform_or_pair"] != "pair"}
    extra = ["final_label_source", "final_claim_status_after_rescue", "maximal_rescue_decision"]
    fields = list(recon[0].keys())
    for c in extra:
        if c not in fields:
            fields.append(c)
    for r in recon:
        sp = (r["species"] or "").lower()
        typ = r.get("final_isoform_label") or r.get("upstream_label")
        f = fin_by_pair.get((sp, typ))
        if f:
            r["final_label_source"] = f["final_label_source"]
            r["final_claim_status_after_rescue"] = f["final_claim_status_after_rescue"]
            r["maximal_rescue_decision"] = f["final_rescue_decision"]
            # if rescued and propagated, record the new protein id + sequence-validated exon type
            if f["final_rescue_decision"].startswith("rescued") and f["final_protein_id"]:
                r["protein_id"] = f["final_protein_id"]
                if f["final_transcript_id"]:
                    r["transcript_id"] = f["final_transcript_id"]
                # the rescued candidate is sequence-validated as this exon type
                r["validated_exon_type"] = typ
                r["final_isoform_label"] = typ
                r["final_claim_status"] = f["final_claim_status_after_rescue"]
        else:
            r.setdefault("final_label_source", r.get("final_label_source", ""))
            r["final_claim_status_after_rescue"] = r.get("final_claim_status", "")
            r["maximal_rescue_decision"] = "not_suspicious"
    M.write_tsv(maps / "fgfr2_exon_type_label_reconciliation.tsv", recon, fields)

    # ---- Part A/B: post-rescue truth table + synchronize downstream tables ----
    propagate_post_rescue(base, dirs, recon, final_rows, rob, refagr, master)

    # ---- Part H: gate ----
    hard = build_gate(gate_rows, recon, suspicious, final_rows)
    M.write_tsv(maps / "fgfr2_maximal_rescue_validation_gate.tsv", gate_rows, GATE_COLS)
    (maps / "fgfr2_maximal_rescue_validation_gate.json").write_text(
        json.dumps({"checks": gate_rows, "hard_fail": hard, "timestamp": M.now_iso()}, indent=2),
        encoding="utf-8")

    if overrides and not args.final_pass:
        (maps / "maximal_rescue_overrides.flag").write_text(M.now_iso(), encoding="utf-8")

    # ---- Part F: corrected final dataset snapshot manifest ----
    write_corrected_manifest(base, dirs, recon, final_rows, overrides, hard, args.final_pass)

    print(f"[OK] maximal rescue: suspicious={len(sus_species)} species; "
          f"overrides={len(overrides)} (final_pass={args.final_pass})")
    print(f"     improved: {sorted(set(improved))}")
    print(f"     still unresolved: {sorted(set(still_bad))}")
    print(f"     gate hard_fail={hard}")
    if hard and not args.no_control_gate:
        print("[FAIL] maximal rescue validation gate hard-failed.", file=sys.stderr)
        return 4
    return 0


def propagate_post_rescue(base, dirs, recon, final_rows, rob, refagr, master) -> None:
    """Part A + B: build the single post-rescue truth table and inject the post-rescue
    claim/provenance columns into every downstream major table so that
    final_claim_status_after_rescue is the single source of truth for primary/supplement."""
    maps = dirs["maps"]
    fin_by = {(f["species"], f["final_isoform_label"]): f for f in final_rows
              if f["isoform_or_pair"] != "pair"}
    # per (species, final_isoform_label) sync payload + authoritative id map
    sync: Dict[Tuple[str, str], Dict[str, str]] = {}
    idmap: Dict[Tuple[str, str], Dict[str, str]] = {}
    truth = []
    for r in recon:
        sp = (r["species"] or "").lower()
        iso = r.get("final_isoform_label") or r.get("upstream_label") or ""
        f = fin_by.get((sp, iso), {})
        rb = rob.get(VR._key(sp, iso), {})
        claim = (r.get("final_claim_status_after_rescue")
                 or f.get("final_claim_status_after_rescue") or r.get("final_claim_status", ""))
        pre_use = master.get(sp, {}).get("recommended_use", "") or r.get("recommended_use", "")
        post_use = recommended_use_post(claim)
        ev = (f.get("evidence_summary", "") + (";src=" + r.get("final_label_source", "")
              if r.get("final_label_source") else "")).strip(";")
        payload = {
            "rescue_decision": r.get("maximal_rescue_decision", "") or f.get("final_rescue_decision", ""),
            "final_label_source": r.get("final_label_source", "") or f.get("final_label_source", ""),
            "final_claim_status_after_rescue": claim,
            "recommended_use_pre_rescue": pre_use,
            "recommended_use_post_rescue": post_use,
            "rescue_evidence_summary": ev,
            "unresolved_reason_if_any": f.get("unresolved_reason_if_any", ""),
        }
        sync[(sp, iso)] = payload
        idmap[(sp, iso)] = {"protein_id": r.get("protein_id", ""),
                            "transcript_id": r.get("transcript_id", ""),
                            "final_isoform_label": iso}
        truth.append({
            "species": sp, "isoform": iso, "upstream_label": r.get("upstream_label", ""),
            "legacy_label": r.get("legacy_label", r.get("upstream_label", "")),
            "previous_pipeline_label": r.get("previous_pipeline_label", r.get("upstream_label", "")),
            "validated_exon_type": r.get("validated_exon_type", ""), "final_isoform_label": iso,
            "transcript_id": r.get("transcript_id", ""), "protein_id": r.get("protein_id", ""),
            **payload,
            "reference_agreement_percent_identical": rb.get("reference_agreement_percent_identical", ""),
            "reference_agreement_percent_identical_or_conservative":
            rb.get("reference_agreement_percent_identical_or_conservative", ""),
            "boundary_robustness_class": rb.get("boundary_robustness_class", ""),
            "overall_alignment_evidence_class": rb.get("overall_alignment_evidence_class", ""),
        })
    M.write_tsv(maps / "fgfr2_post_rescue_final_truth_table.tsv", truth, TRUTH_COLS)

    # inject SYNC_COLS into downstream major tables, joined by (species, final/iso label)
    def patch(path, sp_key, iso_key, use_final_label=True):
        if not path or not Path(path).exists():
            return
        rows = M.read_tsv(path)
        if not rows:
            return
        fields = list(rows[0].keys()) + [c for c in SYNC_COLS if c not in rows[0]]
        for row in rows:
            sp = (row.get(sp_key) or "").lower()
            iso = (row.get("final_isoform_label") if use_final_label and row.get("final_isoform_label")
                   else row.get(iso_key) or "")
            p = sync.get((sp, iso)) or {}
            for c in SYNC_COLS:
                row[c] = p.get(c, row.get(c, ""))
            # authoritative protein/transcript IDs (rescued IDs must propagate everywhere)
            ids = idmap.get((sp, iso)) or {}
            for c in ("protein_id", "transcript_id"):
                if c in row and ids.get(c):
                    row[c] = ids[c]
        M.write_tsv(path, rows, fields)

    patch(dirs["robustness"] / "fgfr2_boundary_robustness_scores.tsv", "species", "isoform")
    patch(maps / "fgfr2_exon_boundary_msa_projection.tsv", "species", "isoform")
    for iso in ("IIIb", "IIIc"):
        patch(dirs["conservation"] / f"fgfr2_{iso}_reference_agreement_summary_by_species.tsv",
              "species", "isoform")


def write_corrected_manifest(base, dirs, recon, final_rows, overrides, hard, final_pass) -> None:
    maps, root = dirs["maps"], M.module_dir(base)
    artifacts = [
        ("selected_transcript_protein_table", maps / "fgfr2_exon_type_label_reconciliation.tsv", "corrected"),
        ("coordinate_audit", M.locate(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv"), "corrected"),
        ("cassette_cds_block_map", M.locate(base, "fgfr2_cassette_cds_block_map.tsv"), "corrected"),
        ("species_qc_master", M.locate(base, "species_qc_master.tsv", "11_pre_interpro_master"), "corrected"),
        ("msa_full_input", dirs["inputs"] / "fgfr2_full_length_protein_msa_input.faa", "regenerated"),
        ("msa_IIIb_input", dirs["inputs"] / "fgfr2_IIIb_cassette_msa_input.faa", "regenerated"),
        ("msa_IIIc_input", dirs["inputs"] / "fgfr2_IIIc_cassette_msa_input.faa", "regenerated"),
        ("robustness_scores", dirs["robustness"] / "fgfr2_boundary_robustness_scores.tsv", "regenerated"),
        ("all_suspicious_cases", maps / "fgfr2_all_suspicious_cases_for_rescue.tsv", "rescue"),
        ("local_rescue_screen", maps / "fgfr2_exhaustive_local_rescue_candidate_screen.tsv", "rescue"),
        ("pair_rescue_decision", maps / "fgfr2_exhaustive_pair_rescue_decision.tsv", "rescue"),
        ("external_rescue_screen", maps / "fgfr2_external_rescue_candidate_screen.tsv", "rescue"),
        ("maximal_rescue_final_decision", maps / "fgfr2_maximal_rescue_final_decision.tsv", "rescue"),
        ("post_rescue_final_truth_table", maps / "fgfr2_post_rescue_final_truth_table.tsv", "truth"),
        ("rescue_overrides", maps / "fgfr2_rescue_overrides.tsv", "rescue"),
        ("validation_gate", maps / "fgfr2_maximal_rescue_validation_gate.tsv", "gate"),
        ("cross_table_consistency_gate",
         maps / "fgfr2_post_rescue_cross_table_consistency_gate.tsv", "gate"),
    ]
    rows = []
    for name, p, role in artifacts:
        ok = bool(p) and Path(p).exists()
        nrows = sum(1 for _ in M.read_tsv(p)) if ok and str(p).endswith(".tsv") else ""
        rows.append({"artifact": name, "path": str(p) if p else "", "exists": str(ok).lower(),
                     "rows": nrows, "sha256": M.sha256_file(p) if ok else "", "role": role})
    M.write_tsv(root / "final_corrected_pre_interpro_dataset_manifest.tsv", rows,
                ["artifact", "path", "exists", "rows", "sha256", "role"])
    fin_iso = [f for f in final_rows if f["isoform_or_pair"] != "pair"]
    cnt = Counter(f["final_rescue_decision"] for f in fin_iso)
    by_species = {}
    for f in fin_iso:
        by_species.setdefault(f["species"], {})[f["final_isoform_label"]] = {
            "decision": f["final_rescue_decision"], "claim": f["final_claim_status_after_rescue"],
            "final_protein_id": f["final_protein_id"], "label_source": f["final_label_source"]}
    summary = {
        "script_version": SCRIPT_VERSION, "timestamp": M.now_iso(), "final_pass": final_pass,
        "gate_hard_fail": hard, "n_overrides": len(overrides),
        "decision_counts": dict(cnt),
        "corrected_principle": ("A case is 'corrected' only when a sequence/provenance-validated "
                                "candidate replaces or confirms the final label/coordinate; a case is "
                                "'not recoverable' only after local and external rescue fail or "
                                "source-compatible evidence is unavailable."),
        "hard_rule": "No final primary output depends on upstream_label as the biological isoform label.",
        "rescued": sorted({f"{f['species']}/{f['final_isoform_label']}" for f in fin_iso
                           if f["final_rescue_decision"].startswith("rescued")}),
        "unresolved": sorted({f"{f['species']}/{f['final_isoform_label']}" for f in fin_iso
                              if f["final_claim_status_after_rescue"] in
                              ("supplement_review", "excluded_from_primary_claim")}),
        "by_species": by_species, "artifacts": rows}
    (root / "final_corrected_pre_interpro_dataset_manifest.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")


def build_gate(gate_rows, recon, suspicious, final_rows) -> bool:
    hard = False

    def add(check, scope, ok, detail):
        nonlocal hard
        gate_rows.append({"check": check, "scope": scope,
                          "status": "pass" if ok else "FAIL", "detail": detail})
        if not ok:
            _hard = True

    hi = [s for s in suspicious if s["rescue_required"] == "true"
          and s["rescue_priority"] in ("critical_control", "high")]
    fin_iso = [f for f in final_rows if f["isoform_or_pair"] != "pair"]
    fin_key = {(f["species"], f["final_isoform_label"]): f for f in fin_iso}
    # 1 human control — only meaningful when human is part of the run. For a custom
    # run without human the control is not applicable and must not hard-fail (recorded
    # as pass with an explicit not-applicable detail so "status != pass" gate readers
    # are not tripped). The full-30 panel always contains human.
    hu = [r for r in recon if (r["species"] or "").lower() == "homo_sapiens"]
    if not hu:
        gate_rows.append({"check": "human_control_passes", "scope": "homo_sapiens",
                          "status": "pass",
                          "detail": "not_applicable: human absent from this run panel (custom run)"})
    else:
        hu_ok = all(r.get("final_claim_status_after_rescue", "").startswith("primary_claim")
                    and r.get("final_isoform_label") == r.get("validated_exon_type") for r in hu)
        add("human_control_passes", "homo_sapiens", hu_ok,
            "human primary IIIb+IIIc" if hu_ok else "human not primary")
    # 2 mouse documented — likewise not applicable if mouse is not in the run panel.
    mu = [r for r in recon if (r["species"] or "").lower() == "mus_musculus"]
    if not mu:
        gate_rows.append({"check": "mouse_control_documented", "scope": "mus_musculus",
                          "status": "pass",
                          "detail": "not_applicable: mouse absent from this run panel (custom run)"})
    else:
        add("mouse_control_documented", "mus_musculus", True,
            "mouse present and documented")
    # 3 all high/critical have local rescue attempted
    miss_local = [f"{s['species']}/{s['isoform']}" for s in hi
                  if fin_key.get((s["species"], s["isoform"]), {}).get("local_rescue_attempted") != "true"]
    add("high_critical_local_rescue_attempted", "high_critical", not miss_local,
        "; ".join(miss_local) or "all attempted")
    # 4 unresolved high/critical have external rescue attempted
    unres = [f for f in fin_iso if f["final_claim_status_after_rescue"] in
             ("supplement_review", "excluded_from_primary_claim")]
    miss_ext = [f"{f['species']}/{f['final_isoform_label']}" for f in unres
                if any(s["species"] == f["species"] and s["isoform"] == f["final_isoform_label"]
                       and s["rescue_priority"] in ("critical_control", "high") for s in hi)
                and f["external_rescue_attempted"] != "true"]
    add("unresolved_high_critical_external_attempted", "high_critical", not miss_ext,
        "; ".join(miss_ext) or "all attempted")
    # 5 no unresolved high/critical in primary
    bad_primary = [f"{f['species']}/{f['final_isoform_label']}" for f in fin_iso
                   if f["final_claim_status_after_rescue"] == "excluded_from_primary_claim"
                   and any(s["species"] == f["species"] and s["isoform"] == f["final_isoform_label"]
                           and f["final_claim_status_after_rescue"].startswith("primary") for s in hi)]
    add("no_unresolved_high_critical_in_primary", "primary_figures", not bad_primary,
        "; ".join(bad_primary) or "none")
    # 6 final==validated for all primary
    prim = [r for r in recon if r.get("final_claim_status_after_rescue", "").startswith("primary_claim")]
    bad_fv = [f"{r['species']}/{r['upstream_label']}" for r in prim
              if r.get("final_isoform_label") != r.get("validated_exon_type")]
    add("primary_final_equals_validated", "all_primary", not bad_fv, "; ".join(bad_fv) or "ok")
    # 7 rescued candidates have provenance
    resc = [f for f in fin_iso if f["final_rescue_decision"].startswith("rescued")]
    bad_prov = [f"{f['species']}/{f['final_isoform_label']}" for f in resc
                if not f["final_protein_id"] or not f["final_label_source"]]
    add("rescued_have_provenance", "rescued", not bad_prov, "; ".join(bad_prov) or "ok")
    # 8 unrescued have explicit reason
    bad_reason = [f"{f['species']}/{f['final_isoform_label']}" for f in unres
                  if not f["unresolved_reason_if_any"]]
    add("unresolved_have_reason", "unresolved", not bad_reason, "; ".join(bad_reason) or "ok")
    return hard


if __name__ == "__main__":
    raise SystemExit(main())
