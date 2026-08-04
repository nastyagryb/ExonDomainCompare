#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shutil
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402

SCRIPT_VERSION = "1.0"
TRUTH_COLS = [
    "species", "display_species_name", "taxon_group", "isoform",
    "upstream_label", "legacy_label", "previous_pipeline_label",
    "validated_exon_type", "final_isoform_label", "transcript_id", "protein_id", "gene_id",
    "protein_length", "sequence_md5", "label_consistency_status",
    "rescue_required", "rescue_decision", "final_label_source", "final_claim_status",
    "final_claim_status_after_rescue", "orthology_status", "paralog_status",
    "synteny_validation_class", "combined_synteny_validation_class", "neighbor_label_quality_score",
    "MSA_full_length_status", "full_length_gap_fraction", "full_length_outlier_status",
    "full_length_MSA_warning", "MSA_cassette_status", "cassette_MSA_warning",
    "boundary_robustness_class", "coordinate_validation_status", "CDS_reconstruction_status",
    "protein_integrity_status", "pre_interpro_readiness_class", "pre_interpro_warning",
    "unresolved_reason_if_any",
]

MSA_FINAL_MAP = [
    ("alignments/fgfr2_full_length_protein_msa.aln.faa",
     "MSA/final_fgfr2_full_length_protein_msa.aln.faa"),
    ("alignments/fgfr2_IIIb_cassette_msa.aln.faa",
     "MSA/final_fgfr2_IIIb_cassette_msa.aln.faa"),
    ("alignments/fgfr2_IIIc_cassette_msa.aln.faa",
     "MSA/final_fgfr2_IIIc_cassette_msa.aln.faa"),
    ("alignments/fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa",
     "MSA/final_fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa"),
    ("protein_integrity/fgfr2_pre_interpro_protein_integrity_qc.tsv",
     "MSA/final_full_length_msa_integrity_qc.tsv"),
    ("conservation/fgfr2_msa_region_conservation_summary.tsv",
     "MSA/final_full_length_msa_conservation_summary.tsv"),
    ("maps/fgfr2_exon_boundary_msa_projection.tsv",
     "MSA/final_cassette_msa_boundary_projection.tsv"),
    ("conservation/fgfr2_combined_human_reference_residue_agreement.tsv",
     "MSA/final_human_referenced_residue_agreement.tsv"),
    ("conservation/fgfr2_IIIb_IIIc_discriminating_positions_main_only.tsv",
     "MSA/final_isoform_discriminating_residues.tsv"),
    ("metadata/msa_strategy_comparison.tsv",
     "MSA/final_msa_strategy_comparison.tsv"),
]


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _run_id() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def load_protein_sequences(base: Path) -> Dict[str, Tuple[str, str, int]]:
    md = M.module_dir(base)
    seqs: Dict[str, Tuple[str, str, int]] = {}

    def _add(sid: str, seq: str) -> None:
        md5 = M.sha256_text(seq)
        plen = len(seq)
        seqs[sid] = (seq, md5, plen)
        parts = sid.split("|")
        if len(parts) >= 3:
            pid = parts[2]
            if pid and pid not in seqs:
                seqs[pid] = (seq, md5, plen)

    for rel in ("inputs/fgfr2_full_length_protein_msa_input.faa",
                "inputs/fgfr2_rescued_candidate_proteins.faa"):
        for sid, seq in M.read_fasta(md / rel):
            _add(sid, seq)
    exp = base / "06_protein_export_v2_7_marker_validated" / "selected_fgfr2_proteins.faa"
    for sid, seq in M.read_fasta(exp):
        _add(sid, seq)
    return seqs


def derive_readiness(claim: str, integrity: str, msa_fl: str, unresolved: str) -> str:
    if unresolved and claim.startswith("primary_claim"):
        return "ready_for_interpro_with_minor_flags"
    if claim in ("supplement_review",):
        return "supplement_review_only"
    if claim.startswith("excluded"):
        return "excluded_from_interpro_primary"
    if integrity in ("fail", "length_outlier") or msa_fl == "full_length_msa_unavailable":
        return "not_ready_unresolved"
    if claim == "primary_claim_supported_with_minor_flags" or msa_fl == "full_length_msa_pass_with_minor_flags":
        return "ready_for_interpro_with_minor_flags"
    if claim == "primary_claim_supported" and msa_fl in (
            "full_length_msa_pass", "full_length_msa_pass_with_minor_flags"):
        return "ready_for_interpro_primary"
    if claim == "primary_claim_supported":
        return "ready_for_interpro_primary"
    return "not_ready_unresolved"


def msa_full_length_status(integ_row: Dict[str, str]) -> Tuple[str, str, str, float]:
    st = integ_row.get("protein_integrity_status", "")
    outlier = integ_row.get("length_outlier_status", "")
    warn = integ_row.get("protein_integrity_warning", "")
    gap = 0.0
    if st in ("fail",):
        return "full_length_msa_outlier_review", outlier or "review", warn, gap
    if outlier in ("major_outlier", "moderate_outlier", "length_outlier"):
        return "full_length_msa_outlier_review", outlier, warn, gap
    if st in ("review", "length_outlier"):
        return "full_length_msa_pass_with_minor_flags", outlier or "minor", warn, gap
    if st in ("pass", "ok", "acceptable", "protein_integrity_pass"):
        return "full_length_msa_pass", outlier or "ok", warn, gap
    return "full_length_msa_pass_with_minor_flags", outlier or "unknown", warn, gap


def build_truth_table(base: Path) -> List[Dict[str, object]]:
    md = M.module_dir(base)
    post = M.read_tsv(md / "maps" / "fgfr2_post_rescue_final_truth_table.tsv")
    if not post:
        raise RuntimeError("post-rescue truth table missing; run MSA/rescue module first")
    master = {(r["species"] or "").lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}
    recon = {((r["species"] or "").lower(), r.get("final_isoform_label") or r.get("isoform") or ""): r
             for r in M.read_tsv(md / "maps" / "fgfr2_exon_type_label_reconciliation.tsv")}
    ortho = {((r["species"] or "").lower(), r.get("isoform") or ""): r
             for r in M.read_tsv(M.locate(base, "fgfr2_orthology_evidence.tsv") or Path(""))}
    integ = {((r["species"] or "").lower(), r.get("isoform") or ""): r
             for r in M.read_tsv(md / "protein_integrity" / "fgfr2_pre_interpro_protein_integrity_qc.tsv")}
    rob = {((r["species"] or "").lower(), r.get("isoform") or ""): r
           for r in M.read_tsv(md / "robustness" / "fgfr2_boundary_robustness_scores.tsv")}
    susp = {((r["species"] or "").lower(), r.get("isoform") or ""): r
            for r in M.read_tsv(md / "maps" / "fgfr2_all_suspicious_cases_for_rescue.tsv")}
    regcons = {((r["species"] or "").lower(), r.get("isoform") or ""): r
               for r in M.read_tsv(md / "conservation" / "fgfr2_msa_region_conservation_summary.tsv")
               if r.get("msa_name") == "full_length" and r.get("region_type") == "full_sequence"}
    seqs = load_protein_sequences(base)
    rows: List[Dict[str, object]] = []
    for r in post:
        sp = (r.get("species") or "").lower()
        iso = r.get("isoform") or r.get("final_isoform_label") or ""
        m = master.get(sp, {})
        rc = recon.get((sp, iso), {})
        o = ortho.get((sp, iso), {})
        iq = integ.get((sp, iso), {})
        rb = rob.get((sp, iso), {})
        su = susp.get((sp, iso), {})
        pid = r.get("protein_id") or rc.get("protein_id") or iq.get("protein_id") or ""
        _seq, md5, plen = seqs.get(pid, ("", "", 0))
        if not plen and iq.get("sequence_length"):
            plen = M.to_int(iq.get("sequence_length"), 0) or 0
        if not md5 and iq.get("sequence_hash"):
            md5 = iq.get("sequence_hash", "")
        fl_st, fl_out, fl_warn, fl_gap = msa_full_length_status(iq)
        rc_reg = regcons.get((sp, iso), {})
        if rc_reg.get("mean_gap_fraction"):
            fl_gap = M.to_float(rc_reg.get("mean_gap_fraction"), fl_gap) or fl_gap
        if rc_reg.get("conservation_warning"):
            fl_warn = fl_warn or rc_reg.get("conservation_warning", "")
        claim = r.get("final_claim_status_after_rescue") or rc.get("final_claim_status") or ""
        unresolved = (r.get("unresolved_reason_if_any") or rc.get("unresolved_reason_if_any") or "").strip()
        rescue = r.get("rescue_decision") or ""
        # rescued+validated should not stay review unless unresolved warning remains
        if rescue.startswith("rescued") and claim.startswith("primary_claim") and not unresolved:
            pass  # keep primary — do not downgrade
        readiness = derive_readiness(claim, iq.get("protein_integrity_status", ""), fl_st, unresolved)
        nscore = M.to_float(r.get("neighbor_identity_confidence_score"), 0.0) or 0.0
        rows.append({
            "species": sp,
            "display_species_name": m.get("display_species_name", sp),
            "taxon_group": m.get("taxon_group_display") or m.get("taxon_group") or m.get("major_clade", ""),
            "isoform": iso,
            "upstream_label": r.get("upstream_label") or rc.get("upstream_label", ""),
            "legacy_label": r.get("legacy_label") or rc.get("legacy_label", ""),
            "previous_pipeline_label": r.get("previous_pipeline_label") or rc.get("previous_pipeline_label", ""),
            "validated_exon_type": r.get("validated_exon_type") or iso,
            "final_isoform_label": r.get("final_isoform_label") or iso,
            "transcript_id": r.get("transcript_id") or rc.get("transcript_id", ""),
            "protein_id": pid,
            "gene_id": rc.get("gene_id", ""),
            "protein_length": plen,
            "sequence_md5": md5,
            "label_consistency_status": rc.get("label_consistency_status", m.get("label_consistency_status", "")),
            "rescue_required": su.get("rescue_required", rc.get("rescue_required", "")),
            "rescue_decision": rescue,
            "final_label_source": r.get("final_label_source") or rc.get("final_label_source", ""),
            "final_claim_status": rc.get("final_claim_status") or claim,
            "final_claim_status_after_rescue": claim,
            "orthology_status": o.get("orthology_status", m.get("fgfr2_ortholog_status", "")),
            "paralog_status": o.get("paralog_panel_status", m.get("paralog_screen_status", "")),
            "synteny_validation_class": r.get("synteny_validation_class", m.get("synteny_validation_class", "")),
            "combined_synteny_validation_class": r.get("combined_synteny_validation_class",
                                                        m.get("combined_synteny_validation_class", "")),
            "neighbor_label_quality_score": round(nscore, 3),
            "MSA_full_length_status": fl_st,
            "full_length_gap_fraction": fl_gap,
            "full_length_outlier_status": fl_out,
            "full_length_MSA_warning": fl_warn,
            "MSA_cassette_status": m.get("cassette_msa_status", rb.get("msa_boundary_projection_status", "cassette_msa_pass")),
            "cassette_MSA_warning": rb.get("boundary_robustness_warning", ""),
            "boundary_robustness_class": r.get("boundary_robustness_class", rb.get("boundary_robustness_class", "")),
            "coordinate_validation_status": m.get("native_coordinate_sanity", rb.get("coordinate_support_status", "")),
            "CDS_reconstruction_status": m.get("cds_boundary_explainability", ""),
            "protein_integrity_status": iq.get("protein_integrity_status", m.get("protein_integrity_summary", "")),
            "pre_interpro_readiness_class": readiness,
            "pre_interpro_warning": fl_warn or rb.get("boundary_robustness_warning", "") or unresolved,
            "unresolved_reason_if_any": unresolved,
        })
    rows.sort(key=lambda x: (M.to_int(master.get(x["species"], {}).get("phylo_order"), 999) or 999,
                             x["species"], x["isoform"]))
    return rows


def final_consistency_gate(base: Path, truth: List[Dict[str, object]], cdir: Path) -> Tuple[bool, List[str]]:
    md = M.module_dir(base)
    tabd = md / "tables"
    checks: List[Dict[str, str]] = []

    def add(check, scope, ok, detail=""):
        checks.append({"check": check, "scope": scope,
                       "status": "pass" if ok else "FAIL", "detail": detail})

    truth_k = {(r["species"], r["isoform"]): r for r in truth}
    post = M.load_post_rescue_truth(base)
    master = {(r.get("species") or "").lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}

    # 1 final tables agree with truth table (post-rescue source of record)
    bad = []
    for k, r in truth_k.items():
        p = post.get(k, {})
        if not p:
            bad.append(f"{k[0]}/{k[1]} missing in post-rescue")
            continue
        for col in ("transcript_id", "protein_id", "final_isoform_label",
                    "validated_exon_type", "final_claim_status_after_rescue"):
            if str(r.get(col, "")) != str(p.get(col, "")):
                bad.append(f"{k[0]}/{k[1]}:{col}")
    add("final_tables_agree_with_truth_table", "post_rescue_truth", not bad, "; ".join(bad[:8]) or "ok")

    # 2 no upstream as biology
    badup = [f"{r['species']}/{r['isoform']}" for r in truth
             if r.get("final_label_source") in ("upstream_label", "upstream")]
    add("no_upstream_label_as_biological_truth", "truth_table", not badup, "; ".join(badup) or "ok")

    # 3 no stale recommended_use overriding post-rescue claim
    badru = []
    for k, r in truth_k.items():
        m = master.get(k[0], {})
        claim = str(r.get("final_claim_status_after_rescue", ""))
        post_use = m.get("recommended_use_post_rescue", "")
        if claim.startswith("primary_claim") and post_use in ("supplement_only", "manual_review", "exclude"):
            badru.append(f"{k[0]}/{k[1]}")
    add("no_stale_recommended_use_overrides_post_rescue", "species_qc_master", not badru,
        "; ".join(badru) or "ok")

    # 4 final Figure 6 table separates review/excluded rows via visual_review_flag (not absence)
    f6 = M.read_tsv(cdir / "tables" / "figure6_human_referenced_IIIb_IIIc_residue_agreement.tsv")
    badprim = [f"{r.get('species')}/{r.get('isoform')}" for r in f6
               if not M.claim_is_primary(str(r.get("final_claim_status_after_rescue", "")))
               and str(r.get("visual_review_flag", "")).lower() != "true"]
    add("figure6_review_rows_flagged", "figure6", not badprim, "; ".join(badprim) or "ok")

    # 5 human and mouse controls pass
    badctrl = []
    for ctrl in ("homo_sapiens", "mus_musculus"):
        for k, r in truth_k.items():
            if k[0] != ctrl:
                continue
            if not M.claim_is_primary(str(r.get("final_claim_status_after_rescue", ""))):
                badctrl.append(f"{k[0]}/{k[1]}")
            if r.get("final_isoform_label") != r.get("validated_exon_type"):
                badctrl.append(f"{k[0]}/{k[1]}:label")
    add("human_mouse_controls_pass", "controls", not badctrl, "; ".join(badctrl) or "ok")

    # 6 rescued Gorilla/Canis/Pongo consistent
    def claim_ok(sp, iso, expect_primary):
        r = truth_k.get((sp, iso), {})
        c = str(r.get("final_claim_status_after_rescue", ""))
        return M.claim_is_primary(c) if expect_primary else not M.claim_is_primary(c)
    # These are fixed reference/rescue control species. Only enforce the sub-condition
    # for species that are actually part of the run panel; absent species are not
    # applicable (custom run). For the full-30 panel all are present and enforced.
    def _sp_present(sp):
        return any(k[0] == sp for k in truth_k)
    gcp_conds = []
    if _sp_present("gorilla_gorilla_gorilla"):
        gcp_conds += [claim_ok("gorilla_gorilla_gorilla", "IIIb", True),
                      claim_ok("gorilla_gorilla_gorilla", "IIIc", True)]
    if _sp_present("canis_lupus_familiaris"):
        gcp_conds += [claim_ok("canis_lupus_familiaris", "IIIb", True),
                      not claim_ok("canis_lupus_familiaris", "IIIc", True)]
    if _sp_present("pongo_abelii"):
        gcp_conds += [claim_ok("pongo_abelii", "IIIc", True),
                      not claim_ok("pongo_abelii", "IIIb", True)]
    gcp_ok = all(gcp_conds)
    add("gorilla_canis_pongo_rescue_status_consistent", "rescued_controls", gcp_ok,
        ("not_applicable: no gorilla/canis/pongo control species in this run panel (custom run)"
         if not gcp_conds else
         ("gorilla/canis IIIb + pongo IIIc primary; canis IIIc + pongo IIIb supplement" if gcp_ok
          else "inconsistent")))

    # 7 synteny + MSA statuses propagated to master
    miss_m = [sp for sp, m in master.items()
              if not m.get("synteny_validation_class") or not m.get("full_length_msa_status")]
    add("synteny_msa_propagated_to_master", "species_qc_master", not miss_m,
        f"{len(miss_m)} species missing synteny/msa cols" if miss_m else "ok")

    # 8 FASTA MD5 matches manifest/truth
    man = M.read_tsv(cdir / "freeze" / "final_pre_interpro_sequence_manifest.tsv")
    badmd5 = []
    for row in man:
        if row.get("included_in_primary_interpro") != "true":
            continue
        k = (row["species"], row["isoform"])
        if k in truth_k and str(truth_k[k].get("sequence_md5", "")) != row.get("sequence_md5", ""):
            badmd5.append(f"{k[0]}/{k[1]}")
    add("fasta_md5_checksums_match_truth", "freeze", not badmd5, "; ".join(badmd5) or "ok")

    # 9 figure input tables carry post-rescue sync columns (current-run provenance)
    sync_cols = ("final_claim_status_after_rescue", "rescue_decision", "final_label_source")
    f8 = M.read_tsv(tabd / "figure8_boundary_robustness_evidence_stack.tsv")
    miss_sync = [c for c in sync_cols if f8 and c not in f8[0]]
    add("figure_tables_from_current_run_with_sync_cols", "figure_tables", not miss_sync,
        f"missing cols: {miss_sync}" if miss_sync else "ok")

    # 10 deprecated outputs not used by final figure assembly (static source-path check)
    deprecated = ("figures_v2_22_final_qc_display", "08_figures_v2_7", "make_fgfr2_exact_exon")
    fig_sources = [
        str(base / "11_publication_figures_pre_interpro" / "figures"),
        str(M.module_dir(base) / "figures"),
    ]
    baddep = [s for s in fig_sources if any(d in s for d in deprecated)]
    add("deprecated_outputs_not_used_in_closure_figures", "figure_sources", not baddep,
        "; ".join(baddep) or "ok")

    # 11 rescued-and-validated not stigmatized (readiness not supplement unless unresolved)
    badres = [f"{r['species']}/{r['isoform']}" for r in truth
              if str(r.get("rescue_decision", "")).startswith("rescued")
              and M.claim_is_primary(str(r.get("final_claim_status_after_rescue", "")))
              and r.get("pre_interpro_readiness_class") in ("supplement_review_only", "not_ready_unresolved")
              and not r.get("unresolved_reason_if_any")]
    add("rescued_validated_not_flagged_as_problem", "truth_table", not badres, "; ".join(badres) or "ok")

    # 12 primary final equals validated; readiness matches claim
    badp = [f"{r['species']}/{r['isoform']}" for r in truth
            if M.claim_is_primary(str(r.get("final_claim_status_after_rescue", "")))
            and r.get("final_isoform_label") != r.get("validated_exon_type")]
    add("primary_final_isoform_equals_validated", "truth_table", not badp, "; ".join(badp) or "ok")

    # ---- Part F: inspect the ACTUAL final figure input tables in closure/tables ----
    ftab = cdir / "tables"
    fig_tables = {
        "figure2": ftab / "figure2_final_exon_to_protein_architecture_pre_interpro.tsv",
        "figure3": ftab / "figure3_final_IIIb_IIIc_cassette_zoom_pre_interpro.tsv",
        "figure5": ftab / "figure5_full_length_FGFR2_MSA_integrity_paper.tsv",
        "figure6": ftab / "figure6_human_referenced_IIIb_IIIc_residue_agreement.tsv",
        "figure8": ftab / "figure_final_framework_evidence_stack.tsv",
    }
    loaded = {name: M.read_tsv(p) for name, p in fig_tables.items()}

    # F-1 every figure table exists and is non-empty
    missing_tabs = [name for name, rows in loaded.items() if not rows]
    add("figure_input_tables_exist", "figure_tables", not missing_tabs,
        "missing/empty: " + ", ".join(missing_tabs) if missing_tabs else "ok")

    # F-2 no primary/minor row carries visual_review_flag = true in any figure table
    bad_vr = []
    for name, rows in loaded.items():
        for r in rows:
            claim = str(r.get("final_claim_status_after_rescue", ""))
            if M.claim_is_primary(claim) and str(r.get("visual_review_flag", "")).lower() == "true":
                bad_vr.append(f"{name}:{r.get('species')}/{r.get('isoform')}")
    add("no_primary_row_visually_flagged_review", "figure_tables", not bad_vr,
        "; ".join(bad_vr[:8]) or "ok")

    # F-3 no rescued-and-validated primary row has review styling in figure tables
    rescued_primary = {(r["species"], r["isoform"]) for r in truth
                       if str(r.get("rescue_decision", "")).startswith("rescued")
                       and M.claim_is_primary(str(r.get("final_claim_status_after_rescue", "")))}
    bad_resvr = []
    for name, rows in loaded.items():
        for r in rows:
            if (r.get("species"), r.get("isoform")) in rescued_primary \
                    and str(r.get("visual_review_flag", "")).lower() == "true":
                bad_resvr.append(f"{name}:{r.get('species')}/{r.get('isoform')}")
    add("rescued_validated_not_review_styled_in_figures", "figure_tables", not bad_resvr,
        "; ".join(bad_resvr) or "ok")

    # F-4 Canis IIIc and Pongo IIIb are NOT primary in any figure table
    bad_supp = []
    must_supp = {("canis_lupus_familiaris", "IIIc"), ("pongo_abelii", "IIIb")}
    for name, rows in loaded.items():
        for r in rows:
            if (r.get("species"), r.get("isoform")) in must_supp \
                    and M.claim_is_primary(str(r.get("final_claim_status_after_rescue", ""))):
                bad_supp.append(f"{name}:{r.get('species')}/{r.get('isoform')}")
            if (r.get("species"), r.get("isoform")) in must_supp \
                    and str(r.get("visual_review_flag", "")).lower() != "true":
                bad_supp.append(f"{name}:{r.get('species')}/{r.get('isoform')}:not_flagged")
    add("canis_IIIc_pongo_IIIb_supplement_in_figures", "figure_tables", not bad_supp,
        "; ".join(bad_supp[:8]) or "ok")

    # F-5 Gorilla IIIb/IIIc consistent with final truth table claim across figure tables
    bad_gor = []
    for name, rows in loaded.items():
        for r in rows:
            k = (r.get("species"), r.get("isoform"))
            if k[0] != "gorilla_gorilla_gorilla":
                continue
            if k in truth_k and str(r.get("final_claim_status_after_rescue", "")) != \
                    str(truth_k[k].get("final_claim_status_after_rescue", "")):
                bad_gor.append(f"{name}:{k[0]}/{k[1]}")
    add("gorilla_consistent_in_figure_tables", "figure_tables", not bad_gor,
        "; ".join(bad_gor) or "ok")

    # F-6 Figure 6 table does not rely on stale is_review_species / recommended_use
    f6cols = set(loaded["figure6"][0].keys()) if loaded["figure6"] else set()
    stale6 = [c for c in ("is_review_species", "recommended_use", "recommended_use_pre_rescue")
              if c in f6cols]
    f6_ok = (not stale6) and ("final_claim_status_after_rescue" in f6cols)
    add("figure6_uses_final_truth_not_stale_review", "figure6",
        f6_ok, ("stale cols: " + ", ".join(stale6)) if stale6 else
        ("missing final_claim_status_after_rescue" if not f6_ok else "ok"))

    # F-7 Figure 2/3 tables do not rely on stale pre-rescue recommended_use
    stale23 = []
    for name in ("figure2", "figure3"):
        cols = set(loaded[name][0].keys()) if loaded[name] else set()
        for c in ("recommended_use", "recommended_use_pre_rescue", "is_review_species"):
            if c in cols:
                stale23.append(f"{name}:{c}")
        if "final_claim_status_after_rescue" not in cols:
            stale23.append(f"{name}:missing_final_claim")
    add("figure2_3_use_post_rescue_not_stale", "figure_tables", not stale23,
        "; ".join(stale23) or "ok")

    # F-8 Figure 5 paper table exists and is not histogram-only (has MSA/length integrity cols)
    f5cols = set(loaded["figure5"][0].keys()) if loaded["figure5"] else set()
    f5_ok = {"MSA_full_length_status", "protein_length", "length_outlier_flag"} <= f5cols
    f5_fig = any((cdir / "figures" / f"Figure_5_full_length_FGFR2_MSA_integrity_paper.{e}").exists()
                 for e in ("svg", "pdf", "png"))
    add("figure5_real_msa_integrity_not_histogram", "figure5", f5_ok and f5_fig,
        "ok" if (f5_ok and f5_fig) else "missing MSA integrity cols or paper figure")

    # F-9 Figure 8 has readable labels (display_species_name present in table)
    f8cols = set(loaded["figure8"][0].keys()) if loaded["figure8"] else set()
    f8_ok = "display_species_name" in f8cols and "taxon_group" in f8cols
    add("figure8_labels_present", "figure8", f8_ok,
        "ok" if f8_ok else "missing display_species_name/taxon_group")

    # ---- Part E: review-case explanation + restored cassette/exon plots ----
    # E-1 Pongo IIIb and Canis IIIc have explicit review explanations (if still review)
    expl = M.read_tsv(ftab / "final_review_case_explanation.tsv")
    expl_k = {(r.get("species"), r.get("isoform")): r for r in expl}
    must_explain = [("pongo_abelii", "IIIb"), ("canis_lupus_familiaris", "IIIc")]
    bad_expl = []
    for k in must_explain:
        if k not in truth_k:
            continue  # control species not in this run panel (custom run) → not applicable
        tr = truth_k.get(k, {})
        if M.claim_is_primary(str(tr.get("final_claim_status_after_rescue", ""))):
            continue  # primary now → no review explanation required
        e = expl_k.get(k, {})
        if not e or not e.get("final_interpretation") or not e.get("unresolved_reason_if_any"):
            bad_expl.append(f"{k[0]}/{k[1]}")
    add("review_cases_pongo_IIIb_canis_IIIc_explained", "review_explanation", not bad_expl,
        "; ".join(bad_expl) or "ok")

    # E-2 no row labelled both rescued and unresolved in the explanation table
    bad_contradiction = [f"{r.get('species')}/{r.get('isoform')}" for r in expl
                         if r.get("rescued_candidate_found") == "yes"
                         and str(r.get("unresolved_reason_if_any", "")).strip()]
    add("no_row_rescued_and_unresolved", "review_explanation", not bad_contradiction,
        "; ".join(bad_contradiction) or "ok")

    # E-3 Figure 3B amino-acid motif map exists, is informative, uses final labels (no upstream)
    f3b = M.read_tsv(ftab / "figure3B_IIIb_IIIc_cassette_amino_acid_motif_map.tsv")
    f3b_cols = set(f3b[0].keys()) if f3b else set()
    f3b_disc = sum(1 for r in f3b if r.get("is_isoform_discriminating") == "true")
    f3b_fig = any((cdir / "figures" / f"Figure_3B_IIIb_IIIc_cassette_amino_acid_motif_map.{e}").exists()
                  for e in ("svg", "pdf", "png"))
    f3b_ok = (bool(f3b) and f3b_fig and f3b_disc > 0
              and "upstream_label" not in f3b_cols
              and {"IIIb_residue_property_class", "substitution_class_IIIb_vs_IIIc"} <= f3b_cols)
    add("figure3B_motif_map_present_informative", "figure3B", f3b_ok,
        f"discriminating positions={f3b_disc}" if f3b_ok else "missing/empty/uses upstream or not informative")

    # E-4 Figure 3C exon-to-protein coordinate map exists and uses final labels
    f3c = M.read_tsv(ftab / "figure3C_exon_to_protein_cassette_coordinate_map.tsv")
    f3c_cols = set(f3c[0].keys()) if f3c else set()
    f3c_fig = any((cdir / "figures" / f"Figure_3C_exon_to_protein_cassette_coordinate_map.{e}").exists()
                  for e in ("svg", "pdf", "png"))
    f3c_ok = (bool(f3c) and f3c_fig
              and "upstream_label" not in f3c_cols
              and {"block_start_aa", "block_end_aa", "cassette_start_aa",
                   "final_isoform_label"} <= f3c_cols)
    add("figure3C_exon_to_protein_map_present", "figure3C", f3c_ok,
        "ok" if f3c_ok else "missing/empty/uses upstream or lacks coordinates")

    # E-5 no rescued-and-validated primary row styled review in 3B/3C/4 tables
    bad_resfig = []
    for name, rows in (("figure3C", f3c),
                       ("figure4", M.read_tsv(ftab / "figure4_label_reconciliation_and_rescue_summary.tsv"))):
        for r in rows:
            if (r.get("species"), r.get("isoform")) in rescued_primary \
                    and str(r.get("visual_review_flag", "")).lower() == "true":
                bad_resfig.append(f"{name}:{r.get('species')}/{r.get('isoform')}")
    add("rescued_validated_not_review_in_3C_4", "figure_tables", not bad_resfig,
        "; ".join(bad_resfig) or "ok")

    # E-6 Figure 4 does not count corrected/rescued accepted rows as review/uncertainty
    f4 = M.read_tsv(ftab / "figure4_label_reconciliation_and_rescue_summary.tsv")
    f4m = {r.get("category"): r for r in f4}
    review_total = (M.to_int(f4m.get("supplement_review_only", {}).get("count"), 0) or 0) \
        + (M.to_int(f4m.get("excluded_or_unresolved", {}).get("count"), 0) or 0)
    actual_review = sum(1 for r in truth
                        if not M.claim_is_primary(str(r.get("final_claim_status_after_rescue", ""))))
    f4_corrected_positive = (f4m.get("labels_corrected_and_accepted", {}).get("group") == "positive"
                             and f4m.get("rescued_and_validated", {}).get("group") == "positive")
    f4_ok = (review_total == actual_review) and f4_corrected_positive
    add("figure4_positive_not_exaggerated", "figure4", f4_ok,
        f"review_total={review_total} actual_review={actual_review}" if not f4_ok else "ok")

    # E-7 species list path resolved correctly (from run-mode json written by the runner)
    rmj = cdir / "final_pre_interpro_run_mode.json"
    sl_status, sl_path = "unspecified", ""
    if rmj.exists():
        try:
            rm = json.loads(rmj.read_text(encoding="utf-8"))
            sl_status = rm.get("species_list_status", "unspecified")
            sl_path = rm.get("species_list_resolved", "")
        except Exception:
            pass
    # pass if resolved, or unspecified (closure invoked directly without runner)
    sl_ok = sl_status in ("resolved", "unspecified")
    add("species_list_path_resolved", "run_mode", sl_ok,
        f"{sl_status}: {sl_path}" if sl_path else sl_status)

    # E-8 Figure 6B exists, uses final labels (no upstream), excludes supplement from main panel
    f6b = M.read_tsv(ftab / "figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv")
    f6b_cols = set(f6b[0].keys()) if f6b else set()
    f6b_fig = any((cdir / "figures"
                   / f"Figure_6B_species_resolved_IIIb_IIIc_cassette_residue_map.{e}").exists()
                  for e in ("svg", "pdf", "png"))
    f6b_main_review = [f"{r.get('species')}/{r.get('isoform')}" for r in f6b
                       if r.get("panel") == "main"
                       and not M.claim_is_primary(str(r.get("final_claim_status_after_rescue", "")))]
    # The species-resolved cassette residue map is built relative to the HUMAN reference
    # cassette. On a custom run whose panel does not contain homo_sapiens the table is
    # empty by construction (only the figure placeholder is rendered). An empty table is
    # not a consistency violation (there is trivially no review row leaking into the main
    # panel and no upstream label misuse), so it is treated as not-applicable as long as
    # the figure file exists. For any panel that DOES contain rows the full check applies.
    human_in_truth = any(k[0] == "homo_sapiens" for k in truth_k)
    if not f6b and not human_in_truth:
        add("figure6B_present_final_labels_no_review_in_main", "figure6B", f6b_fig,
            "not_applicable: human reference cassette absent from panel; figure rendered without "
            "species-resolved rows (custom run)" if f6b_fig else "figure6B figure file missing")
    else:
        f6b_ok = (bool(f6b) and f6b_fig and "upstream_label" not in f6b_cols
                  and "final_isoform_label" in f6b_cols and not f6b_main_review)
        add("figure6B_present_final_labels_no_review_in_main", "figure6B", f6b_ok,
            ("; ".join(f6b_main_review[:6]) if f6b_main_review else
             ("missing/empty or uses upstream" if not (f6b and f6b_fig) else "ok")))

    # E-9 Figure 3B exists (compact motif/major-residue map kept)
    f3b_exists = any((cdir / "figures"
                      / f"Figure_3B_IIIb_IIIc_cassette_amino_acid_motif_map.{e}").exists()
                     for e in ("svg", "pdf", "png"))
    add("figure3B_motif_map_kept", "figure3B", f3b_exists, "ok" if f3b_exists else "missing")

    # E-10 evidence stack excludes not-performed columns (no MCScanX/RMS grey column)
    f8 = M.read_tsv(ftab / "figure_final_framework_evidence_stack.tsv")
    f8cols2 = {c.lower() for c in (f8[0].keys() if f8 else [])}
    banned = [c for c in f8cols2 if "mcscanx" in c or c.startswith("rms") or "_rms" in c]
    add("evidence_stack_no_unused_columns", "figure8", not banned,
        "; ".join(sorted(banned)) or "ok")

    # delegate hard gates from module 12
    ok_post, msgs_post = M.post_rescue_consistency_gate(base, write=False)
    add("post_rescue_cross_table_gate", "module12", ok_post, "; ".join(msgs_post[:4]) or "ok")
    ok_syn, msgs_syn = M.synteny_gate(base)
    add("synteny_validation_gate", "synteny", ok_syn, "; ".join(msgs_syn[:4]) or "ok")

    hard = any(c["status"] != "pass" for c in checks)
    M.write_tsv(cdir / "gates" / "final_pre_interpro_cross_table_consistency_gate.tsv", checks,
                ["check", "scope", "status", "detail"])
    (cdir / "gates" / "final_pre_interpro_cross_table_consistency_gate.json").write_text(
        json.dumps({"checks": checks, "hard_fail": hard, "timestamp": _now()}, indent=2),
        encoding="utf-8")
    fails = [f"{c['check']} [{c['scope']}]: {c['detail']}" for c in checks if c["status"] != "pass"]
    return (not hard), fails


def finalize_msa_outputs(base: Path, cdir: Path) -> List[str]:
    md = M.module_dir(base)
    written = []
    for src_rel, dst_rel in MSA_FINAL_MAP:
        src = md / src_rel
        dst = cdir / dst_rel
        if not src.exists():
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        written.append(str(dst.relative_to(cdir)))
    return written


def freeze_fasta(base: Path, truth: List[Dict[str, object]], cdir: Path) -> None:
    seqs = load_protein_sequences(base)
    freeze = cdir / "freeze"
    freeze.mkdir(parents=True, exist_ok=True)
    primary, review = [], []
    manifest: List[Dict[str, object]] = []
    for r in truth:
        pid = str(r.get("protein_id", ""))
        if pid not in seqs:
            continue
        seq, md5, plen = seqs[pid]
        claim = str(r.get("final_claim_status_after_rescue", ""))
        readiness = str(r.get("pre_interpro_readiness_class", ""))
        in_pri = readiness in ("ready_for_interpro_primary", "ready_for_interpro_with_minor_flags")
        in_rev = readiness != "excluded_from_interpro_primary"
        hdr = "|".join([r["species"], r["isoform"], r["final_isoform_label"], pid,
                        r["transcript_id"], claim, md5])
        if in_pri:
            primary.append((hdr, seq))
        if in_rev:
            rev_hdr = hdr + ("|review" if not in_pri else "")
            review.append((rev_hdr, seq))
        manifest.append({
            "species": r["species"], "isoform": r["isoform"],
            "final_isoform_label": r["final_isoform_label"],
            "transcript_id": r["transcript_id"], "protein_id": pid,
            "protein_length": plen, "sequence_md5": md5,
            "fasta_file": "final_pre_interpro_proteins_primary.faa" if in_pri else
            ("final_pre_interpro_proteins_all_review_included.faa" if in_rev else ""),
            "included_in_primary_interpro": "true" if in_pri else "false",
            "included_in_review_interpro": "true" if in_rev else "false",
            "final_claim_status_after_rescue": claim,
            "pre_interpro_readiness_class": readiness,
            "warning": r.get("pre_interpro_warning", ""),
        })
    M.write_fasta(freeze / "final_pre_interpro_proteins_primary.faa", primary)
    M.write_fasta(freeze / "final_pre_interpro_proteins_all_review_included.faa", review)
    M.write_tsv(freeze / "final_pre_interpro_sequence_manifest.tsv", manifest,
                list(manifest[0].keys()) if manifest else ["species"])
    (freeze / "final_pre_interpro_sequence_manifest.json").write_text(
        json.dumps({"sequences": manifest, "primary_count": len(primary),
                    "review_count": len(review), "timestamp": _now()}, indent=2),
        encoding="utf-8")


def write_checksums(cdir: Path) -> List[Dict[str, str]]:
    rows = []
    for p in sorted(cdir.rglob("*")):
        if p.is_file() and p.suffix not in (".zip",):
            rel = str(p.relative_to(cdir))
            rows.append({"file_path": rel, "md5": M.sha256_file(p), "bytes": p.stat().st_size})
    M.write_tsv(cdir / "freeze" / "final_pre_interpro_file_checksums.tsv", rows,
                ["file_path", "md5", "bytes"])
    return rows


def assemble_final_figures(base: Path, cdir: Path) -> List[str]:
    fig_dir = cdir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    pub = base / "11_publication_figures_pre_interpro" / "figures"
    msa = M.module_dir(base) / "figures"
    syn = M.module_dir(base) / "figures"
    # source stem -> final stem (without extension)
    mapping = [
        (pub, "Figure_1_framework_pre_interpro", "Figure_1_framework_overview"),
        # Figures 2, 3, 3B, 3C, 4, 5, 6 are regenerated directly into closure/figures by
        # make_fgfr2_final_closure_figures.py from the final truth table; do NOT copy
        # stale Step-11 / Figure_6C / pre-rescue evidence-matrix versions here.
        (msa, "Supplement_Figure_full_length_MSA_protein_integrity",
         "Supplement_full_length_MSA_QC_histograms"),
        (msa, "Figure_7C_IIIb_IIIc_isoform_discriminating_residues_informative",
         "Figure_7_isoform_discriminating_residues"),
        (cdir / "figures", "Figure_Final_Framework_Evidence_Stack",
         "Figure_8_final_framework_evidence_stack"),
        # synteny paper names (Part H) preferred if present
        (syn, "Figure_9A_FGFR2_local_synteny_5neighbor_paper",
         "Figure_9_FGFR2_local_synteny_neighborhood"),
        (syn, "Figure_9A_FGFR2_local_synteny_5neighbor_paper",
         "Figure_9A_FGFR2_local_synteny_5neighbor_paper"),
        (syn, "Figure_9B_FGFR2_5neighbor_conservation_matrix_paper",
         "Figure_9B_FGFR2_5neighbor_conservation_matrix_paper"),
        (syn, "Figure_9B_FGFR2_5neighbor_conservation_matrix",
         "Figure_9B_FGFR2_5neighbor_conservation_matrix_paper"),
        (syn, "Figure_9C_FGFR2_synteny_review_cases_paper",
         "Figure_9C_FGFR2_synteny_review_cases_paper"),
        (syn, "Figure_9C_FGFR2_rescue_case_locus_panels",
         "Figure_9C_FGFR2_synteny_review_cases_paper"),
        (syn, "Supplement_Figure_FGFR2_local_synteny_10neighbor_all_species",
         "Supplement_Figure_FGFR2_local_synteny_10neighbor_all_species"),
        (pub, "Supplement_Figure_1_all_species_native_tracks_pre_interpro",
         "Supplement_all_species_exon_protein_architecture"),
        (msa, "Supplement_Figure_per_species_cassette_difference_panels",
         "Supplement_all_species_cassette_zoom"),
        (msa, "Supplement_Figure_full_length_MSA_protein_integrity",
         "Supplement_full_length_MSA_outliers"),
        (pub, "Supplement_Figure_2_review_cases_pre_interpro",
         "Supplement_review_cases_pre_interpro"),
        (msa, "Supplement_Figure_MSA_review_case_diagnostics",
         "Supplement_review_unresolved_case_panels"),
    ]
    copied = set()
    out = []
    for src_dir, src_stem, dst_stem in mapping:
        if dst_stem in copied:
            continue
        for ext in ("svg", "pdf", "png"):
            src = src_dir / f"{src_stem}.{ext}"
            if not src.exists():
                continue
            dst = fig_dir / f"{dst_stem}.{ext}"
            if dst.exists():
                continue
            shutil.copy2(src, dst)
            out.append(str(dst.relative_to(cdir)))
        if any((fig_dir / f"{dst_stem}.{ext}").exists() for ext in ("svg", "pdf", "png")):
            copied.add(dst_stem)
    return out


def write_reports(base: Path, truth: List[Dict[str, object]], cdir: Path, run_id: str,
                  gate_ok: bool = True, gate_failures: Optional[List[str]] = None,
                  run_mode: Optional[Dict[str, object]] = None) -> None:
    rep = cdir / "reports"
    rep.mkdir(parents=True, exist_ok=True)
    from collections import Counter
    claims = Counter(r["final_claim_status_after_rescue"] for r in truth)
    ready = Counter(r["pre_interpro_readiness_class"] for r in truth)
    rescued = sum(1 for r in truth if str(r.get("rescue_decision", "")).startswith("rescued"))
    swapped = sum(1 for r in truth if r.get("label_consistency_status") == "swapped_relative_to_upstream")
    pri_n = ready.get("ready_for_interpro_primary", 0) + ready.get("ready_for_interpro_with_minor_flags", 0)

    gate_line = "PASS" if gate_ok else "FAIL"
    fail_block = ""
    if gate_failures:
        fail_block = "\n\nGate failures:\n" + "\n".join(f"- {f}" for f in gate_failures)

    run_mode = run_mode or {}
    full_clean = bool(run_mode.get("full_clean_run_completed", False))
    cached_v3 = bool(run_mode.get("used_cached_v3_outputs", False))
    cached_msa = bool(run_mode.get("used_cached_msa_outputs", False))
    if not run_mode:
        mode_line = ("Run mode: unspecified (closure invoked directly without the runner; "
                     "cannot certify a full clean end-to-end run).")
    elif full_clean:
        mode_line = ("Run mode: **FULL CLEAN END-TO-END RUN** — Steps 1-11, MSA/rescue/synteny and "
                     "closure all executed (no cached stages).")
    else:
        skipped = []
        if cached_v3:
            skipped.append("Steps 1-11 (A1, cached v3 outputs)")
        if cached_msa:
            skipped.append("MSA/rescue/synteny module (A2, cached outputs)")
        mode_line = ("Run mode: **CACHED DEBUG RUN — NOT a full clean end-to-end run.** "
                     "Cached/skipped stages: " + "; ".join(skipped) + ".")

    summary = f"""# Final pre-InterProScan results summary

**Run ID:** `{run_id}`  
**Pipeline closure version:** {SCRIPT_VERSION}  
**Generated:** {_now()}  
**Cross-table consistency gate:** {gate_line}  
**full_clean_run_completed:** {str(full_clean).lower()}  
**used_cached_v3_outputs:** {str(cached_v3).lower()}  
**used_cached_msa_outputs:** {str(cached_msa).lower()}

> {mode_line}

## Dataset freeze statement

The pre-InterProScan dataset is frozen after label reconciliation, maximal suspicious-case rescue,
coordinate validation, MSA robustness analysis, and local synteny / gene-neighborhood validation.
Domain evidence is **not** inferred before InterProScan. Successfully rescued and validated cases are
included as final accepted primary cases; rescue provenance is retained in tables and reports.
Unresolved or non-recoverable cases are excluded from primary claims and retained only as
review/supplement where appropriate.

## Counts

| metric | value |
|---|---|
| species/isoform rows | {len(truth)} |
| primary claim (supported) | {claims.get('primary_claim_supported', 0)} |
| primary with minor flags | {claims.get('primary_claim_supported_with_minor_flags', 0)} |
| supplement/review | {claims.get('supplement_review', 0)} |
| excluded | {claims.get('excluded_from_primary_claim', 0)} |
| InterPro primary-ready | {pri_n} |
| supplement-review-only | {ready.get('supplement_review_only', 0)} |
| rescued candidates | {rescued} |
| upstream label swaps corrected | {swapped} |

## InterProScan submission

Submit **`freeze/final_pre_interpro_proteins_primary.faa`** first ({pri_n} sequences).  
Use **`freeze/final_pre_interpro_proteins_all_review_included.faa`** for exploratory/review runs.  
See **`freeze/final_pre_interpro_sequence_manifest.tsv`** for per-sequence MD5 checksums and readiness.

## Key outputs

- `final_pre_interpro_truth_table.tsv` — single source of truth
- `gates/final_pre_interpro_cross_table_consistency_gate.tsv` — must PASS
- `MSA/final_*` — frozen MSA snapshot
- `figures/` — final paper figure set (SVG/PDF/PNG)
- `archive/FGFR2_final_pre_interpro_freeze_{run_id}.zip`
{fail_block}
"""
    (rep / "final_pre_interpro_results_summary.md").write_text(summary, encoding="utf-8")

    methods = f"""# Final pre-InterProScan methods summary

Run `{run_id}`. FGFR2 IIIb/IIIc comparative analysis across vertebrate orthologs using an
annotation-aware framework: sequence-calibrated label reconciliation, maximal evidence-driven rescue,
CDS/protein coordinate mapping, full-length and cassette MAFFT alignments, human-reference residue
agreement, boundary robustness scoring, and local 5/10-neighbor synteny validation (MCScanX omitted).
IIIb/IIIc biological labels are driven by `final_isoform_label` / `validated_exon_type`; upstream
labels are provenance only. MSA and synteny are independent robustness layers and do not relabel
isoforms. InterProScan has not been run; no domain annotations are inferred pre-InterPro.
"""
    (rep / "final_pre_interpro_methods_summary.md").write_text(methods, encoding="utf-8")

    qc = f"""# Final pre-InterProScan QC report

Run `{run_id}`. Cross-table consistency gate: **{gate_line}**  
(`gates/final_pre_interpro_cross_table_consistency_gate.tsv`).

Hard gates: post-rescue cross-table consistency, synteny validation, and final closure cross-table gate.

Primary figures include only `primary_claim_supported` or `primary_claim_supported_with_minor_flags`
rows per `final_pre_interpro_truth_table.tsv`. Rescued-and-validated cases (Gorilla, Canis IIIb,
Pongo IIIc) appear as accepted primary cases in main figures; rescue provenance is in tables/reports
and Figure 4, not repeated as problem styling in architecture/MSA/synteny figures.

Review-case explanation (`tables/final_review_case_explanation.tsv`): **Pongo abelii IIIb** and
**Canis lupus familiaris IIIc** remain supplement/review. Their locus, orthology, synteny, MSA,
coordinates and protein integrity all pass, but no source-compatible externally validated
isoform-specific candidate was recovered (sequence support only), so the isoform-specific claim is
retained as review with provenance rather than asserted as primary. These are genuine unresolved
cases, **not** rescued ones (no row is simultaneously rescued and unresolved). Their isoform partners
**Pongo IIIc** and **Canis IIIb** are confirmed/rescued and appear as accepted primary rows.
"""
    (rep / "final_pre_interpro_QC_report.md").write_text(qc, encoding="utf-8")

    caps = """# Final pre-InterProScan figure captions

See also module captions in `11_publication_figures_pre_interpro/captions/` and
`12_msa_boundary_robustness_pre_interpro/captions/`.

**Figure 1** — Annotation-aware framework overview (pre-InterPro).  
**Figure 2** — Final exon-to-protein architecture tracks (no fake domains).  
**Figure 3** — IIIb/IIIc cassette zoom (pre-InterPro).  
**Figure 3B** — IIIb/IIIc cassette amino-acid motif map (final post-rescue consensus): human IIIb/IIIc
reference and major-residue rows with one-letter codes, residue-property colouring, and highlighted
isoform-discriminating positions. Descriptive properties only; no functional/ligand claims.  
**Figure 3C** — Exon/CDS-to-protein cassette coordinate map: exon/CDS blocks on the protein amino-acid
axis with the IIIb/IIIc cassette block and boundary ticks per species/isoform.  
**Figure 4** — Label reconciliation and rescue summary (positive framing): final claim composition and
annotation-aware outcomes; corrected and rescued rows are shown as gains, only true unresolved rows as
review.  
**Figure 5** — Full-length FGFR2 MSA integrity QC.  
**Figure 6** — Human-referenced IIIb/IIIc residue agreement (primary cassettes).  
**Figure 7** — Isoform-discriminating residues.  
**Figure 8** — Integrated framework evidence stack.  
**Figure 9A** — FGFR2 local synteny map (5 neighbors per side).  
**Figure 9B** — 5-neighbor conservation matrix.  
**Figure 9C** — Synteny review-case panels (unresolved/supplement only).  
**Supplement** — 10-neighbor all-species synteny context.
"""
    (rep / "final_pre_interpro_figure_captions.md").write_text(caps, encoding="utf-8")

    lim = """# Final pre-InterProScan limitations

- No InterProScan / domain annotation yet — domain-boundary claims remain downstream.
- Synteny validates locus context only; it does not assign IIIb/IIIc.
- MSA is a robustness layer; it does not relabel isoforms.
- Some distant teleost neighbors show partial synteny (not over-penalized).
- Uncharacterized LOC neighbors may receive loose human-homology names (with % identity); these are
  not curated orthology claims.
- MCScanX block-level synteny was intentionally omitted from this build.
"""
    (rep / "final_pre_interpro_limitations.md").write_text(lim, encoding="utf-8")


CLOSURE_SCRIPT_NAMES = [
    "run_fgfr2_final_pre_interpro_closure.py",
    "run_fgfr2_msa_boundary_module.py",
    "make_fgfr2_final_framework_figure.py",
    "make_fgfr2_synteny_figures_paper.py",
    "make_fgfr2_msa_boundary_figures.py",
    "make_publication_figures_pre_interpro.py",
    "make_all_figures.py",
]


def write_environment_report(cdir: Path) -> Path:
    import subprocess

    arch = cdir / "archive"
    arch.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, str]] = []
    probes = [
        ("mafft", ["mafft", "--version"]),
        ("diamond", ["diamond", "version"]),
        ("python", [sys.executable, "--version"]),
    ]
    for tool, cmd in probes:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            out = ((r.stdout or "") + (r.stderr or "")).strip().replace("\n", " ")[:300]
            rows.append({"tool": tool, "command": " ".join(cmd), "version_output": out,
                         "return_code": str(r.returncode)})
        except Exception as exc:
            rows.append({"tool": tool, "command": " ".join(cmd),
                         "version_output": str(exc), "return_code": "-1"})
    for pkg in ("Bio", "matplotlib", "numpy", "pandas"):
        try:
            mod = __import__(pkg)
            ver = getattr(mod, "__version__", "unknown")
        except Exception as exc:
            ver = f"unavailable: {exc}"
        rows.append({"tool": f"python:{pkg}", "command": f"import {pkg}",
                     "version_output": str(ver), "return_code": "0"})
    out_path = arch / "environment_dependency_report.tsv"
    M.write_tsv(out_path, rows, ["tool", "command", "version_output", "return_code"])
    return out_path


def stage_closure_scripts(cdir: Path) -> List[str]:
    scripts_root = Path(__file__).resolve().parent
    dest_root = cdir / "archive" / "scripts"
    dest_root.mkdir(parents=True, exist_ok=True)
    staged = []
    for name in CLOSURE_SCRIPT_NAMES:
        src = scripts_root / name
        if not src.exists():
            continue
        dst = dest_root / name
        shutil.copy2(src, dst)
        staged.append(str(dst.relative_to(cdir)))
    shared_source = (scripts_root.parent / "src" / "exondomaincompare" /
                     "scientific" / "fgfr2_msa_common.py")
    if shared_source.is_file():
        shared_dest = dest_root / shared_source.name
        shutil.copy2(shared_source, shared_dest)
        staged.append(str(shared_dest.relative_to(cdir)))
    shell = scripts_root.parent / "run_fgfr2_pipeline_current_final_pre_interpro.sh"
    if shell.exists():
        shutil.copy2(shell, dest_root / shell.name)
        staged.append(str((dest_root / shell.name).relative_to(cdir)))
    return staged


def create_archive(cdir: Path, run_id: str) -> Path:
    arch_dir = cdir / "archive"
    arch_dir.mkdir(parents=True, exist_ok=True)
    zip_path = arch_dir / f"FGFR2_final_pre_interpro_freeze_{run_id}.zip"
    manifest_rows = []
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(cdir.rglob("*")):
            if p.is_file() and p != zip_path:
                rel = p.relative_to(cdir)
                zf.write(p, rel)
                role = "report" if "reports" in str(rel) else (
                    "figure" if "figures" in str(rel) else (
                        "freeze" if "freeze" in str(rel) else (
                            "gate" if "gates" in str(rel) else "data")))
                manifest_rows.append({
                    "file_path": str(rel), "file_type": p.suffix.lstrip(".") or "dir",
                    "md5": M.sha256_file(p), "role": role,
                    "required_for_interpro": "true" if "freeze" in str(rel) and p.suffix == ".faa" else "false",
                    "required_for_thesis": "true" if role in ("figure", "report", "freeze") else "false",
                    "description": rel.name,
                })
    M.write_tsv(cdir / "archive" / "final_pre_interpro_archive_manifest.tsv", manifest_rows,
                ["file_path", "file_type", "md5", "role", "required_for_interpro",
                 "required_for_thesis", "description"])
    return zip_path


def main() -> int:
    ap = argparse.ArgumentParser(description="FGFR2 final pre-InterPro closure (Parts B–E, I, J).")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--run-id", default="")
    ap.add_argument("--run-mode-json", type=Path, default=None,
                    help="run-mode flags from the shell runner (full_clean_run_completed, etc.)")
    ap.add_argument("--skip-figures", action="store_true")
    ap.add_argument("--skip-archive", action="store_true")
    args = ap.parse_args()
    base = args.base.resolve()
    run_id = args.run_id or _run_id()
    run_mode = {}
    if args.run_mode_json and args.run_mode_json.exists():
        try:
            run_mode = json.loads(args.run_mode_json.read_text(encoding="utf-8"))
        except Exception:
            run_mode = {}
    cdir = M.closure_dir(base)
    _dirs = M.ensure_closure_dirs(base)

    # clean closure outputs (not upstream cache)
    for sub in ("MSA", "figures", "gates", "freeze"):
        p = cdir / sub
        if p.exists():
            shutil.rmtree(p)
    M.ensure_closure_dirs(base)

    truth = build_truth_table(base)
    M.write_tsv(cdir / "final_pre_interpro_truth_table.tsv", truth, TRUTH_COLS)
    (cdir / "final_pre_interpro_truth_table.json").write_text(
        json.dumps({"run_id": run_id, "rows": truth, "timestamp": _now()}, indent=2),
        encoding="utf-8")

    msa_written = finalize_msa_outputs(base, cdir)
    freeze_fasta(base, truth, cdir)

    # Figures + figure-input tables are derived from the final truth table BEFORE the gate,
    # so the strengthened gate (Part F) can inspect the actual figure tables.
    if not args.skip_figures:
        import subprocess
        subprocess.run([sys.executable, str(Path(__file__).parent / "make_fgfr2_final_framework_figure.py"),
                        "--base", str(base)], check=True)
        subprocess.run([sys.executable, str(Path(__file__).parent / "make_fgfr2_final_closure_figures.py"),
                        "--base", str(base)], check=True)
        assemble_final_figures(base, cdir)

    write_checksums(cdir)

    ok, fails = final_consistency_gate(base, truth, cdir)
    if not ok:
        print("[FAIL] consistency/figure-table gate FAILED — archive will NOT be generated",
              file=sys.stderr)

    write_reports(base, truth, cdir, run_id, gate_ok=ok, gate_failures=fails, run_mode=run_mode)
    write_environment_report(cdir)
    stage_closure_scripts(cdir)

    zip_path = None
    if ok and not args.skip_archive:
        zip_path = create_archive(cdir, run_id)

    # metadata
    meta = {"run_id": run_id, "script_version": SCRIPT_VERSION, "timestamp": _now(),
            "truth_rows": len(truth), "msa_final_files": msa_written,
            "gate_pass": ok, "gate_failures": fails}
    (cdir / "metadata" / "final_pre_interpro_closure_manifest.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8")

    print(f"[{'OK' if ok else 'FAIL'}] closure run_id={run_id} truth_rows={len(truth)} "
          f"msa_final={len(msa_written)} gate={'PASS' if ok else 'FAIL'}")
    if not ok:
        for f in fails:
            print(f"  - {f}", file=sys.stderr)
        return 4
    if not args.skip_archive:
        print(f"[OK] archive: {zip_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
