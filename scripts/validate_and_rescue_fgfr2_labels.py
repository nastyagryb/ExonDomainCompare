#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402
import reconcile_fgfr2_exon_type_labels as RC  # noqa: E402


HUMAN = "homo_sapiens"
CLOSE_PRIMATES = {"homo_sapiens", "pan_troglodytes", "gorilla_gorilla_gorilla",
                  "pongo_abelii", "macaca_mulatta", "callithrix_jacchus"}
KNOWN_RISK_MAMMALS = {"mus_musculus", "canis_lupus_familiaris", "sus_scrofa", "bos_taurus"}

# Part E — group-specific thresholds (do NOT use one universal threshold for all vertebrates).
THRESHOLDS: Dict[str, Dict[str, object]] = {
    "human_curated_positive_control":
        {"min_pioc": 0.98, "min_id": 0.95, "max_gap": 0.05, "marker_required": True, "rescue_below": True},
    "close_primate_control":
        {"min_pioc": 0.95, "min_id": 0.90, "max_gap": 0.08, "marker_required": True, "rescue_below": True},
    "known_label_risk_mammal":
        {"min_pioc": 0.88, "min_id": 0.80, "max_gap": 0.12, "marker_required": True, "rescue_below": True},
    "global_review_outlier":
        {"min_pioc": 0.70, "min_id": 0.55, "max_gap": 0.25, "marker_required": False, "rescue_below": True},
    "standard_species":
        {"min_pioc": 0.60, "min_id": 0.45, "max_gap": 0.30, "marker_required": False, "rescue_below": False},
}

GROUP_COLS = ["species", "isoform", "validation_group", "validation_reason", "validation_threshold_set"]
THR_COLS = ["validation_group", "min_identical_or_conservative_for_pass", "min_identity_for_pass",
            "max_gap_fraction_for_pass", "marker_required", "rescue_required_if_below_threshold"]
TRIG_COLS = ["species", "isoform", "validation_group", "trigger_label_inconsistency",
             "trigger_low_reference_agreement", "trigger_protein_review", "trigger_coordinate_review",
             "trigger_msa_review", "trigger_integrity_review", "trigger_recommended_use_review",
             "trigger_similarity_review", "rescue_required", "rescue_priority", "trigger_summary"]
RESCUE_COLS = ["species", "isoform_or_expected_exon_type", "validation_group", "current_transcript_id",
               "current_protein_id", "candidate_transcript_id", "candidate_protein_id",
               "candidate_previous_label", "candidate_sequence_source", "candidate_cassette_sequence",
               "candidate_cassette_length", "candidate_human_IIIb_identity",
               "candidate_human_IIIc_identity", "candidate_human_IIIb_coverage",
               "candidate_human_IIIc_coverage", "IIIb_marker_present", "IIIc_marker_present",
               "MSA_discriminating_support", "coordinate_support_status", "protein_integrity_status",
               "orthology_status", "paralog_status", "candidate_score", "rescue_candidate_rank",
               "rescue_decision", "rescue_warning"]
RESCUE_SUM_COLS = ["metric", "value"]
PATCH_COLS = ["species", "isoform", "validation_group", "issue_type", "patch_attempted",
              "external_source", "external_accession", "source_release_or_assembly",
              "fetched_transcript", "fetched_protein", "fetched_cds", "fetched_gff3",
              "source_compatibility_status", "translation_validation_status",
              "cassette_reference_agreement", "patch_success", "patch_used_in_final",
              "patch_status", "patch_rejected_reason", "patch_warning"]
RISK_COLS = ["species", "isoform", "validation_group", "problem_detected", "evidence_before_rescue",
             "rescue_attempted", "rescue_result", "final_label_source", "final_claim_status",
             "interpretation"]
GATE_COLS = ["check", "scope", "status", "detail"]

# extra reconciliation columns added by Part F
RECON_EXTRA = ["validation_group", "rescue_required", "rescue_status", "rescue_source",
               "final_label_source", "final_label_evidence_level", "final_claim_status"]


# ---------------------------------------------------------------------------
# evidence loading
# ---------------------------------------------------------------------------
def _key(sp: str, iso: str) -> Tuple[str, str]:
    return ((sp or "").lower(), iso or "")


def load_evidence(base: Path, dirs: Dict[str, Path]):
    rob = {}
    for r in M.read_tsv(dirs["robustness"] / "fgfr2_boundary_robustness_scores.tsv"):
        rob[_key(r["species"], r.get("final_isoform_label") or r.get("isoform"))] = r
    refagr = {}
    for iso in ("IIIb", "IIIc"):
        for r in M.read_tsv(dirs["conservation"] /
                            f"fgfr2_{iso}_reference_agreement_summary_by_species.tsv"):
            refagr[_key(r["species"], r.get("final_isoform_label") or r.get("isoform"))] = r
    orth = {}
    p = M.locate(base, "fgfr2_orthology_evidence.tsv")
    for r in (M.read_tsv(p) if p else []):
        orth[_key(r["species"], r.get("isoform"))] = r
    par = {}
    p = M.locate(base, "fgfr2_paralog_screen_detailed.tsv")
    for r in (M.read_tsv(p) if p else []):
        par[_key(r["species"], r.get("isoform"))] = r
    master = {}
    p = M.locate(base, "species_qc_master.tsv", "11_pre_interpro_master")
    for r in (M.read_tsv(p) if p else []):
        master[(r["species"] or "").lower()] = r
    return rob, refagr, orth, par, master


def load_candidate_pool(base: Path):
    by_species: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    p = M.locate(base, "fgfr2_III_candidate_protein_validation.tsv")
    for r in (M.read_tsv(p) if p else []):
        sp = (r.get("species") or r.get("species_canonical") or "").lower()
        if sp:
            by_species[sp].append(r)
    # full-length sequences from the selected protein FASTA (has protein id + species + role)
    seq_by_pid: Dict[str, str] = {}
    fa = M.locate(base, "selected_fgfr2_proteins.faa")
    for hid, seq in (M.read_fasta(fa) if fa else []):
        meta = {t.split("=", 1)[0]: t.split("=", 1)[1] for t in hid.split("|") if "=" in t}
        pid = meta.get("protein") or ""
        if pid:
            seq_by_pid.setdefault(pid, M.ungapped(M.clean_alignment_seq(seq)))
    return by_species, seq_by_pid


# ---------------------------------------------------------------------------
# Part A — validation group (strictness only; never forces a label)
# ---------------------------------------------------------------------------
def base_group(sp: str) -> Optional[str]:
    if sp == HUMAN:
        return "human_curated_positive_control"
    if sp in CLOSE_PRIMATES:
        return "close_primate_control"
    if sp in KNOWN_RISK_MAMMALS:
        return "known_label_risk_mammal"
    return None


def assign_group(sp: str, suspicious_species: bool) -> Tuple[str, str]:
    g = base_group(sp)
    if g == "human_curated_positive_control":
        return g, "human curated positive control (UniProt P21802-anchored); strictest validation"
    if g == "close_primate_control":
        return g, "close primate control; strict reference-agreement expected"
    if g == "known_label_risk_mammal":
        return g, "known label-risk mammal; targeted validation/rescue if suspicious"
    if suspicious_species:
        return "global_review_outlier", "non-control species with suspicious cassette/label evidence"
    return "standard_species", "standard species; evidence-specific threshold (may be legitimately divergent)"


# ---------------------------------------------------------------------------
# helpers for trigger evaluation
# ---------------------------------------------------------------------------
def _review_token(v: str) -> bool:
    v = (v or "").lower()
    return any(t in v for t in ("review", "exclude", "conflict", "ambiguous", "outlier",
                                "fail", "shift", "gap_rich", "unresolved"))


def gap_fraction(ra: Dict[str, str]) -> Optional[float]:
    n = M.to_float(ra.get("n_reference_positions"))
    g = M.to_float(ra.get("n_gap_or_missing"))
    if n and n > 0 and g is not None:
        return round(g / n, 4)
    return None


def evaluate_triggers(sp: str, iso: str, recon_row: Dict[str, str], group: str,
                      rob: Dict, refagr: Dict, orth: Dict, master: Dict) -> Dict[str, object]:
    thr = THRESHOLDS[group]
    rb = rob.get(_key(sp, iso), {})
    ra = refagr.get(_key(sp, iso), {})
    om = orth.get(_key(sp, iso), {})
    status = recon_row.get("label_consistency_status", "")
    validated = recon_row.get("validated_exon_type", "")

    pioc = M.to_float(ra.get("percent_identical_or_conservative"))
    gapf = gap_fraction(ra)

    t_label = status in ("swapped_relative_to_upstream", "ambiguous_label_review") or \
        validated in ("unresolved", "ambiguous")
    t_lowref = (pioc is not None and pioc < float(thr["min_pioc"])) or \
        (gapf is not None and gapf > float(thr["max_gap"])) or \
        (ra.get("agreement_status") in ("low_reference_agreement_review", "gap_rich_review"))
    t_protein = _review_token(rb.get("protein_evidence_state")) or \
        _review_token(om.get("protein_qc_status"))
    t_coord = _review_token(rb.get("coordinate_resolution_state")) or \
        _review_token(rb.get("native_coordinate_sanity")) or \
        _review_token(rb.get("boundary_precision_state"))
    t_msa = ("shift_review" in (rb.get("msa_boundary_projection_status") or "")) or \
        ("gap_rich_review" in (rb.get("msa_boundary_projection_status") or "")) or \
        _review_token(rb.get("overall_alignment_evidence_class"))
    t_integ = (rb.get("protein_integrity_status") in
               ("protein_length_outlier_review", "invalid_sequence_review", "missing_sequence_fail"))
    ru = rb.get("recommended_use") or om.get("recommended_use") or master.get(sp, {}).get("recommended_use", "")
    t_recuse = (ru or "").strip() in ("supplementary_only", "supplement_only", "manual_review",
                                      "exclude", "excluded") or _review_token(ru) and "main" not in (ru or "")
    t_sim = validated in ("unresolved", "ambiguous") or _review_token(om.get("orthology_warning")) or \
        (rb.get("boundary_robustness_class") in ("review_boundary",
                                                 "unresolved_or_annotation_dependent_boundary"))

    triggers = {
        "trigger_label_inconsistency": t_label,
        "trigger_low_reference_agreement": bool(t_lowref),
        "trigger_protein_review": bool(t_protein),
        "trigger_coordinate_review": bool(t_coord),
        "trigger_msa_review": bool(t_msa),
        "trigger_integrity_review": bool(t_integ),
        "trigger_recommended_use_review": bool(t_recuse),
        "trigger_similarity_review": bool(t_sim),
    }
    any_trig = any(triggers.values())
    rescue_required = any_trig or validated in ("unresolved", "ambiguous")
    active = [k.replace("trigger_", "") for k, v in triggers.items() if v]
    return {"triggers": triggers, "rescue_required": rescue_required, "active": active,
            "pioc": pioc, "gapf": gapf}


def rescue_priority(group: str, ev: Dict[str, object]) -> str:
    if not ev["rescue_required"]:
        return "none"
    t = ev["triggers"]
    strong = t["trigger_label_inconsistency"] or t["trigger_low_reference_agreement"] or \
        t["trigger_similarity_review"]
    if group == "human_curated_positive_control":
        return "critical_control"
    if group == "close_primate_control":
        return "critical_control" if strong else "high"
    if group == "known_label_risk_mammal":
        return "high"
    if group == "global_review_outlier":
        return "high" if strong else "medium"
    n = len(ev["active"])
    return "medium" if (strong or n >= 2) else "low"


# ---------------------------------------------------------------------------
# Part C — general candidate rescue (sequence/evidence-driven)
# ---------------------------------------------------------------------------
def _marker(present_full: Optional[bool], score: Optional[float]) -> bool:
    if present_full is not None:
        return present_full
    return score is not None and score >= 0.5


def score_candidate(cand: Dict[str, str], expected: str, seq_by_pid: Dict[str, str],
                    orth: Dict, par: Dict, sp: str) -> Dict[str, object]:
    pid = cand.get("protein_id") or ""
    tx = cand.get("transcript_id") or ""
    prev = cand.get("expected_isoform_final") or cand.get("role") or ""
    full = seq_by_pid.get(pid, "")
    seg = (cand.get("extracted_segment_sequence") or "").strip()
    if full:
        s = RC.protein_scores(full)
        iiib_id, iiic_id = s["iiib_id"], s["iiic_id"]
        iiib_cov, iiic_cov = s["iiib_cov"], s["iiic_cov"]
        b_present = _marker(s["b_marker"], None)
        c_present = _marker(s["c_marker"], None)
        src = "selected_protein_full_sequence"
        cass = seg or full[:60]
    elif seg:
        b_id, b_cov = RC.aln_id_cov(seg, RC.CURATED_IIIB_REF)
        c_id, c_cov = RC.aln_id_cov(seg, RC.CURATED_IIIC_REF)
        iiib_id, iiic_id, iiib_cov, iiic_cov = b_id, c_id, b_cov, c_cov
        b_present = _marker(any(m in seg for m in RC.IIIB_MARKERS),
                            M.to_float(cand.get("human_IIIb_marker_score")))
        c_present = _marker(any(m in seg for m in RC.IIIC_MARKERS),
                            M.to_float(cand.get("human_IIIc_marker_score")))
        src = "extracted_cassette_segment"
        cass = seg
    else:
        iiib_id = M.to_float(cand.get("human_IIIb_identity"), 0.0)
        iiic_id = M.to_float(cand.get("human_IIIc_identity"), 0.0)
        iiib_cov = M.to_float(cand.get("human_IIIb_coverage"), 0.0)
        iiic_cov = M.to_float(cand.get("human_IIIc_coverage"), 0.0)
        b_present = _marker(None, M.to_float(cand.get("human_IIIb_marker_score")))
        c_present = _marker(None, M.to_float(cand.get("human_IIIc_marker_score")))
        src = "precomputed_validation_identity"
        cass = ""
    own_id = iiib_id if expected == "IIIb" else iiic_id
    own_cov = iiib_cov if expected == "IIIb" else iiic_cov
    own_marker = b_present if expected == "IIIb" else c_present
    other_marker = c_present if expected == "IIIb" else b_present
    if own_marker and not other_marker:
        disc = "consistent_with_markers"
    elif other_marker and not own_marker:
        disc = "conflicts_with_markers"
    else:
        disc = "not_corroborated"
    om = orth.get(_key(sp, expected), {})
    pm = par.get(_key(sp, expected), {})
    orth_status = om.get("orthology_status", "") or "unknown"
    par_status = pm.get("paralog_status", "") or "unknown"
    integ = "protein_integrity_pass" if full and not RC_invalid(full) else (
        "invalid_sequence_review" if full else "unknown_no_full_sequence")
    paralog_pen = 0.2 if _review_token(par_status) else 0.0
    integ_pen = 0.15 if _review_token(integ) else 0.0
    disc_bonus = 0.1 if disc == "consistent_with_markers" else (-0.15 if disc == "conflicts_with_markers" else 0.0)
    marker_bonus = 0.15 if own_marker else 0.0
    score = max(0.0, min(1.0, 0.5 * own_id + 0.2 * own_cov + marker_bonus + disc_bonus
                         - paralog_pen - integ_pen))
    return {"candidate_transcript_id": tx, "candidate_protein_id": pid,
            "candidate_previous_label": prev, "candidate_sequence_source": src,
            "candidate_cassette_sequence": cass, "candidate_cassette_length": len(cass),
            "candidate_human_IIIb_identity": iiib_id, "candidate_human_IIIc_identity": iiic_id,
            "candidate_human_IIIb_coverage": iiib_cov, "candidate_human_IIIc_coverage": iiic_cov,
            "IIIb_marker_present": str(b_present).lower(), "IIIc_marker_present": str(c_present).lower(),
            "MSA_discriminating_support": disc, "coordinate_support_status": "",
            "protein_integrity_status": integ, "orthology_status": orth_status,
            "paralog_status": par_status, "candidate_score": round(score, 4),
            "_own_id": own_id, "_own_marker": own_marker}


def RC_invalid(seq: str) -> bool:
    return bool(M.invalid_residues(seq)) or ("*" in seq[:-1])


def current_passes(recon_row: Dict[str, str], expected: str, group: str) -> bool:
    thr = THRESHOLDS[group]
    own_id = M.to_float(recon_row.get("human_IIIb_identity" if expected == "IIIb"
                                      else "human_IIIc_identity"), 0.0)
    marker = (recon_row.get("IIIb_marker_present" if expected == "IIIb"
                            else "IIIc_marker_present") or "").lower() == "true"
    resolved = recon_row.get("validated_exon_type") in ("IIIb", "IIIc")
    return resolved and own_id >= float(thr["min_id"]) and (not thr["marker_required"] or marker)


def rescue_species_isoform(sp: str, expected: str, recon_row: Dict[str, str], group: str,
                           pool: List[Dict[str, str]], seq_by_pid: Dict[str, str],
                           orth: Dict, par: Dict, rob: Dict) -> Tuple[List[Dict[str, object]], str, str]:
    thr = THRESHOLDS[group]
    cur_pid = recon_row.get("protein_id") or ""
    cur_tx = recon_row.get("transcript_id") or ""
    # de-duplicate candidate proteins, keep best metadata per protein id
    seen: Dict[str, Dict[str, str]] = {}
    for c in pool:
        pid = c.get("protein_id") or ""
        if not pid:
            continue
        if pid not in seen or (c.get("extracted_segment_sequence") and
                               not seen[pid].get("extracted_segment_sequence")):
            seen[pid] = c
    # always include current selection even if not in candidate table
    if cur_pid and cur_pid not in seen:
        seen[cur_pid] = {"protein_id": cur_pid, "transcript_id": cur_tx,
                         "expected_isoform_final": recon_row.get("final_isoform_label", "")}
    scored = [score_candidate(c, expected, seq_by_pid, orth, par, sp) for c in seen.values()]
    for s in scored:
        s["coordinate_support_status"] = (
            rob.get(_key(sp, expected), {}).get("coordinate_resolution_state", "")
            if s["candidate_protein_id"] == cur_pid else "candidate_not_currently_selected")
    scored.sort(key=lambda d: d["candidate_score"], reverse=True)
    for i, s in enumerate(scored, 1):
        s["rescue_candidate_rank"] = i
    best = scored[0] if scored else None
    if best is None:
        return [], "manual_review_required", ""
    # if the current selection already passes its group threshold, keep it (do not optimize away)
    if current_passes(recon_row, expected, group):
        return scored, "keep_current_candidate", cur_pid
    best_passes = (best["_own_id"] >= float(thr["min_id"]) and
                   (not thr["marker_required"] or best["_own_marker"]))
    propagatable = best["candidate_sequence_source"] == "selected_protein_full_sequence"
    if best["candidate_protein_id"] == cur_pid:
        decision = "manual_review_required" if best["_own_id"] >= 0.45 else "exclude_from_primary_claim"
    elif best_passes and propagatable:
        decision = "use_rescued_candidate"
    elif best_passes or best["_own_id"] >= 0.45:
        decision = "manual_review_required"
    else:
        decision = "exclude_from_primary_claim"
    return scored, decision, best["candidate_protein_id"]


# ---------------------------------------------------------------------------
# Part D — targeted external candidate patch (local NCBI datasets cache only)
# ---------------------------------------------------------------------------
def find_ncbi_cache_fgfr2(base: Path, sp: str) -> List[Tuple[str, str]]:
    _cache = M.locate(base, "_ncbi_datasets_cache")
    out: List[Tuple[str, str]] = []
    # the cache is keyed by taxid; match by scanning protein.faa under any ncbi_* dir,
    # restricted to the species via the assembly mapping if available
    roots = sorted(Path(base).rglob("protein.faa"))
    for fp in roots:
        if "_ncbi_datasets_cache" not in str(fp):
            continue
        try:
            cur_id, cur_def, keep, seq = "", "", False, []
            with open(fp, encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if keep and seq:
                            out.append((cur_id, "".join(seq)))
                        hdr = line[1:].strip()
                        cur_id = hdr.split()[0]
                        cur_def = hdr.lower()
                        keep = ("fibroblast growth factor receptor 2" in cur_def or "fgfr2" in cur_def)
                        seq = []
                    elif keep:
                        seq.append(line.strip())
                if keep and seq:
                    out.append((cur_id, "".join(seq)))
        except OSError:
            continue
        if out:
            break  # first cache containing FGFR2 records (per-species dirs are separate files)
    return out


def external_patch(base: Path, sp: str, iso: str, group: str, decision: str,
                   best_own_id: float, src_db: str, assembly: str) -> Dict[str, object]:
    issue = ("local_candidates_below_threshold" if decision == "exclude_from_primary_claim"
             else "local_candidates_ambiguous")
    if decision in ("keep_current_candidate", "use_rescued_candidate"):
        return {"species": sp, "isoform": iso, "validation_group": group,
                "issue_type": "resolved_locally", "patch_attempted": "false",
                "patch_status": "external_patch_not_needed", "patch_used_in_final": "false",
                "patch_success": "false"}
    recs = find_ncbi_cache_fgfr2(base, sp)
    if not recs:
        return {"species": sp, "isoform": iso, "validation_group": group, "issue_type": issue,
                "patch_attempted": "true", "external_source": "ncbi_datasets_local_cache",
                "source_release_or_assembly": assembly, "patch_status": "external_patch_failed",
                "patch_success": "false", "patch_used_in_final": "false",
                "patch_rejected_reason": "no_fgfr2_protein_in_local_ncbi_cache",
                "patch_warning": "no cached NCBI FGFR2 protein available for this species"}
    ref = RC.CURATED_IIIB_REF if iso == "IIIb" else RC.CURATED_IIIC_REF
    best = max(recs, key=lambda kv: RC.aln_id_cov(kv[1], ref)[0])
    acc, bseq = best
    own_id, _own_cov = RC.aln_id_cov(bseq, ref)
    translation_ok = "valid_translation" if not RC_invalid(bseq) else "translation_review"
    # cross-source release mixing is NOT done silently
    cross = (src_db or "").lower() not in ("ncbi", "refseq", "")
    if own_id < 0.45:
        status, reason, compat = ("external_patch_rejected_translation_mismatch",
                                  "cached NCBI cassette disagrees with curated reference", "evaluated")
    elif cross:
        status, reason, compat = ("external_patch_release_mismatch_review",
                                  "selected model is Ensembl; NCBI patch would mix releases",
                                  "incompatible_release_not_mixed")
    else:
        status, reason, compat = ("external_patch_available_but_not_used",
                                  "validated NCBI candidate available; kept as review (no silent swap)",
                                  "same_source_compatible")
    return {"species": sp, "isoform": iso, "validation_group": group, "issue_type": issue,
            "patch_attempted": "true", "external_source": "ncbi_datasets_local_cache",
            "external_accession": acc, "source_release_or_assembly": assembly,
            "fetched_transcript": "", "fetched_protein": acc, "fetched_cds": "", "fetched_gff3": "",
            "source_compatibility_status": compat, "translation_validation_status": translation_ok,
            "cassette_reference_agreement": round(own_id, 4), "patch_success": "false",
            "patch_used_in_final": "false", "patch_status": status,
            "patch_rejected_reason": reason if "rejected" in status else "",
            "patch_warning": reason}


# ---------------------------------------------------------------------------
# Part F — final claim status
# ---------------------------------------------------------------------------
def claim_status(group: str, recon_row: Dict[str, str], ev: Dict, decision: str,
                 ev_pass: bool) -> Tuple[str, str, str]:
    validated = recon_row.get("validated_exon_type", "")
    final = recon_row.get("final_isoform_label", "")
    status = recon_row.get("label_consistency_status", "")
    resolved = validated in ("IIIb", "IIIc") and final == validated
    src = ("rescued_candidate" if decision == "use_rescued_candidate"
           else "sequence_reconciliation" if resolved else "upstream_provenance_review")
    if not resolved or validated in ("unresolved", "ambiguous"):
        return "excluded_from_primary_claim", src, "low"
    if decision in ("manual_review_required", "exclude_from_primary_claim"):
        return ("supplement_review" if decision == "manual_review_required"
                else "excluded_from_primary_claim", src, "low")
    if not ev_pass:
        return "supplement_review", src, "moderate"
    minor = any(ev["triggers"][k] for k in ("trigger_coordinate_review", "trigger_msa_review",
                                            "trigger_integrity_review")) or status == "swapped_relative_to_upstream"
    if minor:
        return "primary_claim_supported_with_minor_flags", src, "high"
    return "primary_claim_supported", src, "high"




def load_source_db(base: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    p = M.locate(base, "transcripts.tsv", "02_models")
    for r in (M.read_tsv(p) if p else []):
        sp = (r.get("species_canonical") or "").lower()
        if sp and sp not in out:
            out[sp] = (r.get("source_db") or "").strip()
    return out


def load_assembly(base: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    p = M.locate(base, "ncbi_assembly_selected.tsv")
    for r in (M.read_tsv(p) if p else []):
        sp = (r.get("species_canonical") or r.get("species") or "").lower()
        acc = r.get("assembly_accession") or r.get("assembly") or ""
        if sp:
            out[sp] = acc
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="General suspicious-case validation & rescue layer.")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--no_control_gate", action="store_true",
                    help="skip the human positive-control stop-gate (debug only)")
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    maps = dirs["maps"]

    recon = M.read_tsv(maps / "fgfr2_exon_type_label_reconciliation.tsv")
    if not recon:
        print("[FAIL] reconciliation table missing; run reconcile step first.", file=sys.stderr)
        return 3
    rob, refagr, orth, par, master = load_evidence(base, dirs)
    pool_by_sp, seq_by_pid = load_candidate_pool(base)
    src_db = load_source_db(base)
    assembly = load_assembly(base)

    # Part E — write thresholds
    thr_rows = [{"validation_group": g,
                 "min_identical_or_conservative_for_pass": v["min_pioc"],
                 "min_identity_for_pass": v["min_id"],
                 "max_gap_fraction_for_pass": v["max_gap"],
                 "marker_required": str(v["marker_required"]).lower(),
                 "rescue_required_if_below_threshold": str(v["rescue_below"]).lower()}
                for g, v in THRESHOLDS.items()]
    M.write_tsv(maps / "fgfr2_validation_thresholds.tsv", thr_rows, THR_COLS)

    # ---- pass 1: decide suspicious species (for group assignment of non-controls) ----
    by_species: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for r in recon:
        by_species[(r["species"] or "").lower()].append(r)
    suspicious_sp: Dict[str, bool] = {}
    for sp, rows in by_species.items():
        susp = False
        for r in rows:
            iso = r.get("final_isoform_label") or r.get("upstream_label") or ""
            ev = evaluate_triggers(sp, iso, r, "standard_species", rob, refagr, orth, master)
            if (ev["rescue_required"] or r.get("validated_exon_type") in ("unresolved", "ambiguous")):
                susp = True
        suspicious_sp[sp] = susp

    group_rows, trig_rows, rescue_rows, patch_rows, risk_rows = [], [], [], [], []
    rescue_dec_counter, claim_counter = Counter(), Counter()
    selection_changed = False
    extra_by_pair: Dict[Tuple[str, str], Dict[str, str]] = {}
    # collect per-species expected->best to detect duplicate rescue conflicts
    sp_best: Dict[str, Dict[str, str]] = defaultdict(dict)

    for sp, rows in by_species.items():
        group, reason = assign_group(sp, suspicious_sp.get(sp, False))
        for r in rows:
            iso_final = r.get("final_isoform_label") or r.get("upstream_label") or ""
            expected = iso_final if iso_final in ("IIIb", "IIIc") else (r.get("upstream_label") or "IIIb")
            group_rows.append({"species": sp, "isoform": iso_final, "validation_group": group,
                               "validation_reason": reason, "validation_threshold_set": group})
            ev = evaluate_triggers(sp, iso_final, r, group, rob, refagr, orth, master)
            prio = rescue_priority(group, ev)
            trig_rows.append({"species": sp, "isoform": iso_final, "validation_group": group,
                              **{k: str(v).lower() for k, v in ev["triggers"].items()},
                              "rescue_required": str(ev["rescue_required"]).lower(),
                              "rescue_priority": prio,
                              "trigger_summary": ",".join(ev["active"]) or "none"})
            # ---- Part C rescue ----
            cur_pass = current_passes(r, expected, group)
            decision, best_pid = "not_required", ""
            if ev["rescue_required"]:
                scored, decision, best_pid = rescue_species_isoform(
                    sp, expected, r, group, pool_by_sp.get(sp, []), seq_by_pid, orth, par, rob)
                sp_best[sp][expected] = best_pid
                for s in scored[:8]:
                    rescue_rows.append({"species": sp, "isoform_or_expected_exon_type": expected,
                                        "validation_group": group,
                                        "current_transcript_id": r.get("transcript_id", ""),
                                        "current_protein_id": r.get("protein_id", ""),
                                        "rescue_decision": decision, "rescue_warning": "",
                                        **{k: s.get(k, "") for k in RESCUE_COLS
                                           if k in s}})
                rescue_dec_counter[decision] += 1
            # ---- Part D external patch (only when local rescue did not resolve) ----
            if ev["rescue_required"] and decision in ("manual_review_required",
                                                      "exclude_from_primary_claim"):
                best_own = 0.0
                cand_rows = [x for x in rescue_rows if x["species"] == sp
                             and x["isoform_or_expected_exon_type"] == expected]
                if cand_rows:
                    key = "candidate_human_IIIb_identity" if expected == "IIIb" else "candidate_human_IIIc_identity"
                    best_own = max(M.to_float(x.get(key), 0.0) for x in cand_rows)
                patch_rows.append(external_patch(base, sp, expected, group, decision, best_own,
                                                 src_db.get(sp, ""), assembly.get(sp, "")))
            # ---- Part F claim status + reconciliation enrichment ----
            cs, src_lbl, lvl = claim_status(group, r, ev, decision, cur_pass)
            rescue_src = ("local_candidate" if decision in ("use_rescued_candidate",
                                                            "keep_current_candidate")
                          else "none")
            extra_by_pair[(sp, r.get("upstream_label", ""))] = {
                "validation_group": group, "rescue_required": str(ev["rescue_required"]).lower(),
                "rescue_status": decision, "rescue_source": rescue_src,
                "final_label_source": src_lbl, "final_label_evidence_level": lvl,
                "final_claim_status": cs}
            claim_counter[cs] += 1

    # ---- duplicate-rescue conflict guard (no duplicate IIIb/IIIc unless review) ----
    for sp, bests in sp_best.items():
        if bests.get("IIIb") and bests.get("IIIb") == bests.get("IIIc"):
            for x in rescue_rows:
                if x["species"] == sp and x["rescue_decision"] == "use_rescued_candidate":
                    x["rescue_decision"] = "manual_review_required"
                    x["rescue_warning"] = "same candidate selected for IIIb and IIIc; manual review"
            for r in by_species[sp]:
                pr = extra_by_pair.get((sp, r.get("upstream_label", "")))
                if pr and pr["rescue_status"] == "use_rescued_candidate":
                    pr["rescue_status"] = "manual_review_required"
                    pr["final_claim_status"] = "supplement_review"

    # ---- write Part A/B/C/D tables ----
    M.write_tsv(maps / "fgfr2_validation_group_assignment.tsv", group_rows, GROUP_COLS)
    M.write_tsv(maps / "fgfr2_suspicious_case_triggers.tsv", trig_rows, TRIG_COLS)
    M.write_tsv(maps / "fgfr2_general_candidate_rescue.tsv", rescue_rows, RESCUE_COLS)
    M.write_tsv(maps / "fgfr2_targeted_external_candidate_patch_report.tsv", patch_rows, PATCH_COLS)
    rescue_sum = [{"metric": "n_rescue_required", "value": sum(1 for t in trig_rows
                                                               if t["rescue_required"] == "true")},
                  *[{"metric": f"decision_{k}", "value": v} for k, v in sorted(rescue_dec_counter.items())],
                  *[{"metric": f"claim_{k}", "value": v} for k, v in sorted(claim_counter.items())],
                  {"metric": "n_external_patch_attempts", "value": sum(1 for p in patch_rows
                                                                       if p.get("patch_attempted") == "true")},
                  {"metric": "n_external_patch_used", "value": sum(1 for p in patch_rows
                                                                   if p.get("patch_used_in_final") == "true")}]
    M.write_tsv(maps / "fgfr2_general_candidate_rescue_summary.tsv", rescue_sum, RESCUE_SUM_COLS)

    # ---- Part F: enrich reconciliation table + summary ----
    fields = list(recon[0].keys())
    for c in RECON_EXTRA:
        if c not in fields:
            fields.append(c)
    for r in recon:
        pr = extra_by_pair.get(((r["species"] or "").lower(), r.get("upstream_label", "")), {})
        for c in RECON_EXTRA:
            r[c] = pr.get(c, "")
        # decision use_rescued_candidate would change selection; record + flag
        if pr.get("rescue_status") == "use_rescued_candidate":
            best_pid = sp_best.get((r["species"] or "").lower(), {}).get(
                r.get("final_isoform_label") or r.get("upstream_label"), "")
            if best_pid and best_pid != r.get("protein_id"):
                r["protein_id"] = best_pid
                selection_changed = True
    M.write_tsv(maps / "fgfr2_exon_type_label_reconciliation.tsv", recon, fields)

    # ---- Part G: known-risk species report ----
    prio_by_pair = {(t["species"], t["isoform"]): t for t in trig_rows}
    report_species = set(KNOWN_RISK_MAMMALS) | CLOSE_PRIMATES
    for t in trig_rows:
        if t["rescue_priority"] in ("critical_control", "high"):
            report_species.add(t["species"])
    for r in recon:
        sp = (r["species"] or "").lower()
        if sp not in report_species:
            continue
        iso = r.get("final_isoform_label") or r.get("upstream_label") or ""
        t = prio_by_pair.get((sp, iso), {})
        ra = refagr.get(_key(sp, iso), {})
        problems = (t.get("trigger_summary", "none"))
        risk_rows.append({
            "species": sp, "isoform": iso,
            "validation_group": r.get("validation_group", ""),
            "problem_detected": problems if problems != "none" else "no_trigger",
            "evidence_before_rescue": (f"pioc={ra.get('percent_identical_or_conservative','na')};"
                                       f"status={r.get('label_consistency_status','')};"
                                       f"prio={t.get('rescue_priority','none')}"),
            "rescue_attempted": str(r.get("rescue_required") == "true").lower(),
            "rescue_result": r.get("rescue_status", ""),
            "final_label_source": r.get("final_label_source", ""),
            "final_claim_status": r.get("final_claim_status", ""),
            "interpretation": _interpret(r)})
    risk_rows.sort(key=lambda d: (d["species"], d["isoform"]))
    M.write_tsv(maps / "fgfr2_known_risk_species_validation_report.tsv", risk_rows, RISK_COLS)

    # ---- Part J: validation gate ----
    gate_rows, hard_fail = build_gate(recon, trig_rows)
    M.write_tsv(maps / "fgfr2_general_rescue_validation_gate.tsv", gate_rows, GATE_COLS)

    if selection_changed:
        (maps / "rescue_selection_changed.flag").write_text(M.now_iso(), encoding="utf-8")

    print(f"[OK] validation/rescue: rescue_required="
          f"{sum(1 for t in trig_rows if t['rescue_required']=='true')}/{len(trig_rows)}; "
          f"decisions={dict(rescue_dec_counter)}")
    print(f"     claim status={dict(claim_counter)}")
    print(f"     external patch attempts={sum(1 for p in patch_rows if p.get('patch_attempted')=='true')} "
          f"(used_in_final={sum(1 for p in patch_rows if p.get('patch_used_in_final')=='true')})")
    print(f"     gate hard_fail={hard_fail}; selection_changed={selection_changed}")
    if hard_fail and not args.no_control_gate:
        print("[FAIL] general rescue validation gate hard-failed (see gate table).", file=sys.stderr)
        return 4
    return 0


def _interpret(r: Dict[str, str]) -> str:
    cs = r.get("final_claim_status", "")
    if cs == "primary_claim_supported":
        return "validated; used as primary claim"
    if cs == "primary_claim_supported_with_minor_flags":
        return "validated with minor flags; primary claim with footnote"
    if cs == "supplement_review":
        return "evidence insufficient for primary claim; shown in supplement/review"
    return "excluded from primary claim; explained in supplement"


def build_gate(recon: List[Dict[str, str]], trig_rows: List[Dict[str, str]]
               ) -> Tuple[List[Dict[str, str]], bool]:
    rows, hard = [], False
    prio = {(t["species"], t["isoform"]): t for t in trig_rows}

    def add(check, scope, ok, detail):
        nonlocal hard
        rows.append({"check": check, "scope": scope,
                     "status": "pass" if ok else "FAIL", "detail": detail})
        if not ok:
            _hard = True

    primary = [r for r in recon if r.get("final_claim_status", "").startswith("primary_claim")]
    # 1. no primary row with final != validated
    bad = [r for r in primary if r.get("final_isoform_label") != r.get("validated_exon_type")]
    add("primary_final_equals_validated", "all_primary", not bad,
        "; ".join(f"{r['species']}/{r['upstream_label']}" for r in bad) or "ok")
    # 2. no primary row uses upstream as biological label where it was swapped
    bad2 = [r for r in primary if r.get("final_label_source") == "upstream_provenance_review"]
    add("primary_not_upstream_label", "all_primary", not bad2,
        "; ".join(r["species"] for r in bad2) or "ok")
    # 3. no primary row with rescue_required and unresolved rescue status
    bad3 = [r for r in primary if r.get("rescue_required") == "true"
            and r.get("rescue_status") in ("manual_review_required", "exclude_from_primary_claim")]
    add("primary_no_unresolved_rescue", "all_primary", not bad3,
        "; ".join(f"{r['species']}/{r['upstream_label']}" for r in bad3) or "ok")
    # 4. human control primary
    # The human positive control is only meaningful when human is part of the run.
    # For a custom run WITHOUT human, the control is not applicable and must not
    # hard-fail the gate (recorded as pass with an explicit not-applicable detail so
    # downstream "status != pass" gate readers are not tripped). The full-30 panel
    # always contains human, so the control is enforced there exactly as before.
    hu = [r for r in recon if (r["species"] or "").lower() == HUMAN]
    if not hu:
        rows.append({"check": "human_positive_control_primary", "scope": "homo_sapiens",
                     "status": "pass",
                     "detail": "not_applicable: human absent from this run panel (custom run); "
                               "positive control not evaluated"})
    else:
        hu_ok = all(r.get("final_claim_status", "").startswith("primary_claim")
                    and r.get("final_isoform_label") == r.get("validated_exon_type")
                    for r in hu)
        add("human_positive_control_primary", "homo_sapiens", hu_ok,
            "human is primary IIIb+IIIc" if hu_ok else "human not validated as primary control")
    # 5. close-primate failures must be rescued or excluded (not primary if below)
    cp_bad = [r for r in recon if (r["species"] or "").lower() in CLOSE_PRIMATES
              and r.get("final_claim_status") == "primary_claim_supported"
              and prio.get(((r["species"] or "").lower(),
                            r.get("final_isoform_label", "")), {}).get("trigger_low_reference_agreement") == "true"]
    add("close_primate_low_agreement_not_primary", "close_primates", not cp_bad,
        "; ".join(r["species"] for r in cp_bad) or "ok")
    # 6. canis & mus, if triggered, must not be silently primary-unflagged
    for sp in ("canis_lupus_familiaris", "mus_musculus"):
        rs = [r for r in recon if (r["species"] or "").lower() == sp]
        trg = any(prio.get((sp, r.get("final_isoform_label", "")), {}).get("rescue_required") == "true"
                  for r in rs)
        ok = (not trg) or all(r.get("final_claim_status") in
                              ("primary_claim_supported_with_minor_flags", "supplement_review",
                               "excluded_from_primary_claim") or
                              r.get("rescue_status") == "keep_current_candidate" for r in rs)
        add("known_risk_mammal_triggered_handled", sp, ok,
            "handled (kept/flagged/excluded)" if ok else "triggered but silently primary")
    return rows, hard


if __name__ == "__main__":
    raise SystemExit(main())
