#!/usr/bin/env python3
"""
score_fgfr2_boundary_robustness.py  (MSA boundary-robustness sprint, Parts 7 + 11)

Transparent, component-based boundary robustness score per species/isoform, combining
annotation/coordinate resolution, codon-phase/boundary precision, protein QC, MSA boundary
projection, conservation/gap evidence and protein integrity. All component values and the
exact weights are written out; nothing is hidden and uncertain cases are not forced to exact.

Also emits review-case MSA diagnostics (Part 11): difficult cases are explained, not hidden.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fgfr2_msa_common as M  # noqa: E402


WEIGHTS = {
    "component_annotation_score": 0.20,
    "component_codon_phase_score": 0.15,
    "component_protein_qc_score": 0.15,
    "component_msa_projection_score": 0.25,
    "component_conservation_score": 0.15,
    "component_integrity_score": 0.10,
}
SCORE_COLS = ["species", "isoform", "final_isoform_label", "upstream_label",
              "label_consistency_status", "protein_id", "transcript_id", "recommended_use",
              "coordinate_resolution_state", "boundary_precision_state", "cds_boundary_confidence",
              "protein_evidence_state", "native_coordinate_sanity", "normalized_slot_sanity",
              "msa_boundary_projection_status", "boundary_window_gap_fraction",
              "boundary_window_conservation_score", "cassette_conservation_score",
              "protein_integrity_status", "splice_site_qc_status_if_available",
              "component_annotation_score", "component_coordinate_score",
              "component_codon_phase_score", "component_protein_qc_score",
              "component_msa_projection_score", "component_conservation_score",
              "component_integrity_score", "boundary_robustness_score",
              "boundary_robustness_class",
              # Part F — explicit reference-agreement components
              "reference_agreement_percent_identical",
              "reference_agreement_percent_identical_or_conservative",
              "left_boundary_reference_agreement", "right_boundary_reference_agreement",
              "cassette_core_reference_agreement", "discriminating_residue_support",
              "gap_rich_penalty", "reference_guided_boundary_score",
              "overall_alignment_evidence_class", "robustness_warning"]
DIAG_COLS = ["species", "isoform", "recommended_use", "review_reason_short",
             "native_coordinate_sanity", "normalized_slot_sanity", "protein_evidence_state",
             "boundary_precision_state", "msa_boundary_projection_status",
             "boundary_window_gap_fraction", "cassette_conservation_score",
             "boundary_robustness_score", "splice_site_qc_status_if_available",
             "protein_integrity_status", "final_interpretation", "suggested_display_location"]

ANNOT = {"cds_coordinate_resolved": 1.0, "protein_overlay_no_cds_model": 0.4,
         "coordinate_conflict_review": 0.3, "coordinate_unresolved": 0.0}
PHASE = {"codon_boundary_exact": 1.0, "known_split_codon_boundary": 0.9,
         "phase_not_available_but_coordinate_resolved": 0.7,
         "nucleotide_sequence_unavailable": 0.3, "boundary_unresolved": 0.0}
PROT = {"protein_validated": 1.0, "protein_transcript_level_only": 0.8,
        "protein_ambiguous_review": 0.5, "protein_conflict_review": 0.3, "protein_unavailable": 0.0}
MSA = {"msa_boundary_projected_high_confidence": 1.0,
       "msa_boundary_projected_with_minor_gaps": 0.8,
       "msa_boundary_projected_gap_rich_review": 0.4, "msa_boundary_shift_review": 0.3,
       "msa_boundary_unresolved": 0.0, "msa_boundary_not_applicable": 0.5}
INTEG = {"protein_integrity_pass": 1.0, "protein_integrity_pass_with_minor_warning": 0.85,
         "protein_length_outlier_review": 0.5, "invalid_sequence_review": 0.2,
         "missing_sequence_fail": 0.0}


def cls(score: float) -> str:
    if score >= 0.85:
        return "robust_boundary"
    if score >= 0.70:
        return "supported_boundary_with_minor_flags"
    if score >= 0.50:
        return "review_boundary"
    return "unresolved_or_annotation_dependent_boundary"


def main() -> int:
    ap = argparse.ArgumentParser(description="Boundary robustness scoring (Parts 7,11).")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)

    recon = M.load_label_reconciliation(base)

    def fin(sp, up):
        return M.final_label(recon, sp, up)

    refined = {(r["species"].lower(), r["isoform"]): r for r in
               M.read_tsv(M.require(base, "fgfr2_refined_uncertainty_classes.tsv"))}
    cds_audit = {(r["species"].lower(), r["isoform"]): r for r in
                 M.read_tsv(M.require(base, "cds_phase_boundary_audit.tsv"))}
    master = {r["species"].lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}
    # MSA-derived tables already use the FINAL biological label (keyed by final isoform)
    proj = {(r["species"].lower(), r["isoform"]): r for r in
            M.read_tsv(dirs["maps"] / "fgfr2_exon_boundary_msa_projection.tsv")}
    integ = {(r["species"].lower(), r["isoform"]): r for r in
             M.read_tsv(dirs["protein_integrity"] / "fgfr2_pre_interpro_protein_integrity_qc.tsv")}
    splice = {(r["species"].lower(), r["isoform"]): r for r in
              M.read_tsv(dirs["splice_qc"] / "fgfr2_splice_site_boundary_qc.tsv")}
    region = M.read_tsv(dirs["conservation"] / "fgfr2_msa_region_conservation_summary.tsv")
    reg = {}
    for r in region:
        reg.setdefault((r["species"].lower(), r["isoform"]), {})[r["region_type"]] = r

    # Part F — reference-agreement evidence (keyed by final isoform)
    cons = dirs["conservation"]
    ref_sp = {}
    for fn in ("fgfr2_IIIb_reference_agreement_summary_by_species.tsv",
               "fgfr2_IIIc_reference_agreement_summary_by_species.tsv"):
        for r in M.read_tsv(cons / fn):
            ref_sp[(r["species"].lower(), r["isoform"])] = r
    seg = {}
    for r in M.read_tsv(cons / "fgfr2_cassette_segment_agreement_summary.tsv"):
        seg[(r["species"].lower(), r["isoform"], r["segment_type"])] = r
    disc_cols = {M.to_int(r["combined_alignment_col"]) for r in
                 M.read_tsv(cons / "fgfr2_IIIb_IIIc_discriminating_positions_informative.tsv")
                 if r.get("position_class") == "isoform_discriminating_conserved"}
    comb_agree = {}
    for r in M.read_tsv(cons / "fgfr2_combined_human_reference_residue_agreement.tsv"):
        comb_agree.setdefault((r["species"].lower(), r["isoform"]), []).append(r)

    rows: List[Dict[str, object]] = []
    for (sp, up_iso), rf in refined.items():
        iso = fin(sp, up_iso)  # FINAL biological label used for all MSA joins/output
        mr = master.get(sp, {})
        pr = proj.get((sp, iso), {})
        ig = integ.get((sp, iso), {})
        sq = splice.get((sp, iso), {})
        rg = reg.get((sp, iso), {})

        coord_state = rf.get("coordinate_resolution_state", "")
        bprec = rf.get("boundary_precision_state", "")
        pevid = rf.get("protein_evidence_state", "")
        cds_conf = cds_audit.get((sp, iso), {}).get("cds_boundary_confidence", "")
        msa_status = pr.get("boundary_projection_status", "msa_boundary_unresolved")
        lw = M.to_float(pr.get("left_boundary_gap_fraction_w5"))
        rw = M.to_float(pr.get("right_boundary_gap_fraction_w5"))
        bw_gap = max([v for v in (lw, rw) if v is not None], default=None)
        cass_cons = M.to_float((rg.get("cassette", {}) or {}).get("mean_conservation_score"))
        lcons = M.to_float((rg.get("left_boundary_window", {}) or {}).get("mean_conservation_score"))
        rcons = M.to_float((rg.get("right_boundary_window", {}) or {}).get("mean_conservation_score"))
        bw_cons = None
        bwl = [v for v in (lcons, rcons) if v is not None]
        if bwl:
            bw_cons = round(sum(bwl) / len(bwl), 4)
        integ_status = ig.get("protein_integrity_status", "missing_sequence_fail")
        splice_status = sq.get("splice_site_qc_status", "splice_site_sequence_unavailable")

        c_annot = ANNOT.get(coord_state, 0.5)
        c_phase = PHASE.get(bprec, 0.5)
        c_prot = PROT.get(pevid, 0.5)
        c_msa = MSA.get(msa_status, 0.0)
        if cass_cons is not None or bw_gap is not None:
            cc = cass_cons if cass_cons is not None else 0.5
            gg = bw_gap if bw_gap is not None else 0.5
            c_cons = round(max(0.0, min(1.0, 0.5 * cc + 0.5 * (1.0 - gg))), 4)
        else:
            c_cons = 0.5
        c_integ = INTEG.get(integ_status, 0.5)

        score = round(
            WEIGHTS["component_annotation_score"] * c_annot
            + WEIGHTS["component_codon_phase_score"] * c_phase
            + WEIGHTS["component_protein_qc_score"] * c_prot
            + WEIGHTS["component_msa_projection_score"] * c_msa
            + WEIGHTS["component_conservation_score"] * c_cons
            + WEIGHTS["component_integrity_score"] * c_integ, 4)
        rclass = cls(score)
        warns = []
        for label, val in (("annotation", c_annot), ("codon_phase", c_phase),
                           ("protein_qc", c_prot), ("msa_projection", c_msa),
                           ("conservation_gap", c_cons), ("integrity", c_integ)):
            if val < 0.5:
                warns.append(f"low_{label}({round(val,2)})")

        # ---- Part F: explicit reference-agreement components ----
        rsp = ref_sp.get((sp, iso), {})
        ref_pid = M.to_float(rsp.get("percent_identical"))
        ref_pic = M.to_float(rsp.get("percent_identical_or_conservative"))

        def seg_agr(stype):
            return M.to_float((seg.get((sp, iso, stype), {}) or {}).get(
                "percent_identical_or_conservative"))
        left_agr = seg_agr("left_boundary_window")
        right_agr = seg_agr("right_boundary_window")
        core_agr = seg_agr("cassette_core")
        full_gap = M.to_float((seg.get((sp, iso, "full_cassette"), {}) or {}).get("gap_fraction"))
        gap_pen = round(full_gap, 4) if full_gap is not None else (
            round(bw_gap, 4) if bw_gap is not None else "")
        # discriminating-residue support: fraction of isoform-discriminating columns where this
        # species matches its own isoform's human reference (identical)
        cagr = comb_agree.get((sp, iso), [])
        disc_hits = [r for r in cagr if M.to_int(r.get("alignment_col")) in disc_cols]
        disc_support = (round(sum(1 for r in disc_hits
                                  if r.get("agreement_class") == "identical_to_human") / len(disc_hits), 4)
                        if disc_hits else "")
        comps = [v for v in (left_agr, right_agr, core_agr,
                             (disc_support if isinstance(disc_support, float) else None)) if v is not None]
        if comps:
            base_ref = sum(comps) / len(comps)
            gp = gap_pen if isinstance(gap_pen, float) else 0.0
            ref_guided = round(max(0.0, min(1.0, base_ref * (1.0 - 0.5 * gp))), 4)
        else:
            ref_guided = ""
        # overall alignment-evidence class
        if ref_pic is None and not comps:
            evid = "alignment_unresolved"
        elif msa_status == "msa_boundary_shift_review":
            evid = "alignment_shift_review"
        elif isinstance(gap_pen, float) and gap_pen >= 0.5:
            evid = "alignment_gap_rich_review"
        elif isinstance(ref_guided, float) and ref_guided >= 0.8 and c_msa >= 0.8:
            evid = "alignment_supports_boundary"
        elif isinstance(ref_guided, float) and ref_guided >= 0.6:
            evid = "alignment_supports_boundary_with_minor_variation"
        else:
            evid = "alignment_unresolved"

        rows.append({
            "species": rf.get("species", sp), "isoform": iso,
            "final_isoform_label": iso, "upstream_label": up_iso,
            "label_consistency_status": recon.get((sp, up_iso), {}).get(
                "label_consistency_status", "no_reconciliation"),
            "protein_id": rf.get("protein_id", ""), "transcript_id": rf.get("transcript_id", ""),
            "recommended_use": mr.get("recommended_use", ""),
            "coordinate_resolution_state": coord_state, "boundary_precision_state": bprec,
            "cds_boundary_confidence": cds_conf, "protein_evidence_state": pevid,
            "native_coordinate_sanity": mr.get("native_coordinate_sanity", ""),
            "normalized_slot_sanity": mr.get("normalized_slot_sanity", ""),
            "msa_boundary_projection_status": msa_status,
            "boundary_window_gap_fraction": bw_gap if bw_gap is not None else "",
            "boundary_window_conservation_score": bw_cons if bw_cons is not None else "",
            "cassette_conservation_score": cass_cons if cass_cons is not None else "",
            "protein_integrity_status": integ_status,
            "splice_site_qc_status_if_available": splice_status,
            "component_annotation_score": c_annot, "component_coordinate_score": c_annot,
            "component_codon_phase_score": c_phase, "component_protein_qc_score": c_prot,
            "component_msa_projection_score": c_msa, "component_conservation_score": c_cons,
            "component_integrity_score": c_integ, "boundary_robustness_score": score,
            "boundary_robustness_class": rclass,
            "reference_agreement_percent_identical": ref_pid if ref_pid is not None else "",
            "reference_agreement_percent_identical_or_conservative": ref_pic if ref_pic is not None else "",
            "left_boundary_reference_agreement": left_agr if left_agr is not None else "",
            "right_boundary_reference_agreement": right_agr if right_agr is not None else "",
            "cassette_core_reference_agreement": core_agr if core_agr is not None else "",
            "discriminating_residue_support": disc_support,
            "gap_rich_penalty": gap_pen, "reference_guided_boundary_score": ref_guided,
            "overall_alignment_evidence_class": evid,
            "robustness_warning": ";".join(warns)})

    M.write_tsv(dirs["robustness"] / "fgfr2_boundary_robustness_scores.tsv", rows, SCORE_COLS)
    M.write_tsv(dirs["robustness"] / "fgfr2_boundary_robustness_component_weights.tsv",
                [{"component": k, "weight": v,
                  "description": d} for (k, v), d in zip(WEIGHTS.items(), [
                      "annotation / coordinate resolution state",
                      "codon phase / boundary precision (split codon is not an error)",
                      "protein marker QC evidence", "MSA boundary projection quality",
                      "cassette conservation and boundary-window gap evidence",
                      "pre-InterPro protein integrity"])],
                ["component", "weight", "description"])

    # ---- Part 11: review-case diagnostics ----
    score_by = {(r["species"].lower(), r["isoform"]): r for r in rows}
    diag: List[Dict[str, object]] = []
    for (sp, iso), r in score_by.items():
        ruse = r["recommended_use"]
        proj_review = "review" in str(r["msa_boundary_projection_status"])
        if M.recommended_use_token(ruse) == "main_figure" and not proj_review:
            continue  # robust main species without MSA review -> not a review case
        mr = master.get(sp, {})
        rclass = r["boundary_robustness_class"]
        if rclass == "unresolved_or_annotation_dependent_boundary":
            interp, loc = "annotation-dependent / unresolved boundary", "exclude_from_primary_claim"
        elif M.recommended_use_token(ruse) == "supplement":
            interp, loc = "supplement/review species, interpreted separately", "supplement_review_panel"
        elif proj_review:
            interp, loc = "main-eligible but MSA boundary projection flagged for review", "main_figure_with_warning"
        else:
            interp, loc = "main with minor flags", "main_figure_with_warning"
        diag.append({
            "species": r["species"], "isoform": iso, "recommended_use": ruse,
            "review_reason_short": mr.get("review_reason_short", ""),
            "native_coordinate_sanity": r["native_coordinate_sanity"],
            "normalized_slot_sanity": r["normalized_slot_sanity"],
            "protein_evidence_state": r["protein_evidence_state"],
            "boundary_precision_state": r["boundary_precision_state"],
            "msa_boundary_projection_status": r["msa_boundary_projection_status"],
            "boundary_window_gap_fraction": r["boundary_window_gap_fraction"],
            "cassette_conservation_score": r["cassette_conservation_score"],
            "boundary_robustness_score": r["boundary_robustness_score"],
            "splice_site_qc_status_if_available": r["splice_site_qc_status_if_available"],
            "protein_integrity_status": r["protein_integrity_status"],
            "final_interpretation": interp, "suggested_display_location": loc})
    M.write_tsv(dirs["review_diagnostics"] / "fgfr2_msa_review_case_diagnostics.tsv", diag, DIAG_COLS)

    cc = Counter(r["boundary_robustness_class"] for r in rows)
    print(f"[OK] boundary robustness: {len(rows)} rows; class={dict(cc)}")
    print(f"     review-case diagnostics rows={len(diag)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
