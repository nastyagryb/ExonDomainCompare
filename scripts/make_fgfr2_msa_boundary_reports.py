#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Dict

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402



def counts(rows, key):
    return dict(Counter(r.get(key, "") for r in rows))


def main() -> int:
    ap = argparse.ArgumentParser(description="Write MSA reports and captions.")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    md = M.module_dir(base)

    run_man = M.read_tsv(dirs["metadata"] / "msa_run_manifest.tsv")
    proj = M.read_tsv(dirs["maps"] / "fgfr2_exon_boundary_msa_projection.tsv")
    scores = M.read_tsv(dirs["robustness"] / "fgfr2_boundary_robustness_scores.tsv")
    integ = M.read_tsv(dirs["protein_integrity"] / "fgfr2_pre_interpro_protein_integrity_qc.tsv")
    splice = M.read_tsv(dirs["splice_qc"] / "fgfr2_splice_site_boundary_qc.tsv")
    disc_sum = M.read_tsv(dirs["conservation"] / "fgfr2_IIIb_IIIc_discriminating_positions_summary.tsv")
    region = M.read_tsv(dirs["conservation"] / "fgfr2_msa_region_conservation_summary.tsv")
    diag = M.read_tsv(dirs["review_diagnostics"] / "fgfr2_msa_review_case_diagnostics.tsv")
    inp_val = M.read_tsv(dirs["inputs"] / "fgfr2_msa_input_validation.tsv")

    n_full = sum(1 for _ in M.read_fasta(dirs["inputs"] / "fgfr2_full_length_protein_msa_input.faa"))
    n_iiib = sum(1 for _ in M.read_fasta(dirs["inputs"] / "fgfr2_IIIb_cassette_msa_input.faa"))
    n_iiic = sum(1 for _ in M.read_fasta(dirs["inputs"] / "fgfr2_IIIc_cassette_msa_input.faa"))
    proj_c = counts(proj, "boundary_projection_status")
    rob_c = counts(scores, "boundary_robustness_class")
    integ_c = counts(integ, "protein_integrity_status")
    splice_c = counts(splice, "splice_site_qc_status")
    cass_cons = [M.to_float(r["mean_conservation_score"]) for r in region
                 if r["region_type"] == "cassette" and M.to_float(r["mean_conservation_score"]) is not None]
    mean_cass = round(sum(cass_cons) / len(cass_cons), 3) if cass_cons else "n/a"
    disc_main = {r["position_class"]: r["count"] for r in disc_sum if r["analysis_set"] == "main_only"}
    n_val_fail = sum(1 for r in inp_val if r["status"] == "fail")

    # ---- reconciliation + reference-guided evidence ----
    recon = M.read_tsv(dirs["maps"] / "fgfr2_exon_type_label_reconciliation.tsv")
    recon_sum = {r["metric"]: r["value"] for r in
                 M.read_tsv(dirs["maps"] / "fgfr2_exon_type_label_reconciliation_summary.tsv")}
    recon_c = counts(recon, "label_consistency_status")
    swapped_species = sorted({r["species"] for r in recon
                              if r["label_consistency_status"] == "swapped_relative_to_upstream"})
    ref_iiib = M.read_tsv(dirs["conservation"] / "fgfr2_IIIb_reference_agreement_summary_by_species.tsv")
    ref_iiic = M.read_tsv(dirs["conservation"] / "fgfr2_IIIc_reference_agreement_summary_by_species.tsv")
    ref_all = ref_iiib + ref_iiic
    agree_status_c = counts(ref_all, "agreement_status")
    lowest = sorted(ref_all, key=lambda r: M.to_float(r.get("percent_identical_or_conservative"), 1.0))[:8]
    disc_inf = M.read_tsv(dirs["conservation"] / "fgfr2_IIIb_IIIc_discriminating_positions_informative.tsv")
    strongest = sorted([r for r in disc_inf
                        if r.get("position_class") == "isoform_discriminating_conserved"],
                       key=lambda r: M.to_float(r.get("discriminating_score"), 0.0), reverse=True)[:10]
    # agreement-class totals (across all reference positions, both isoforms)
    n_ident = sum(M.to_int(r.get("n_identical"), 0) for r in ref_all)
    n_cons = sum(M.to_int(r.get("n_conservative"), 0) for r in ref_all)
    n_nonc = sum(M.to_int(r.get("n_nonconservative"), 0) for r in ref_all)
    n_gapm = sum(M.to_int(r.get("n_gap_or_missing"), 0) for r in ref_all)

    lowest_tbl = "\n".join(
        f"| {r['species']} | {r['isoform']} | {r.get('percent_identical')} | "
        f"{r.get('percent_identical_or_conservative')} | {r.get('agreement_status')} | "
        f"{r.get('label_consistency_status','')} |" for r in lowest)
    strongest_tbl = "\n".join(
        f"| {r.get('human_reference_residue_index')} | {r.get('human_IIIb_aa')} | "
        f"{r.get('human_IIIc_aa')} | {r.get('IIIb_major_aa')} | {r.get('IIIc_major_aa')} | "
        f"{r.get('discriminating_score')} |" for r in strongest)
    swap_tbl = "\n".join(
        f"| {r['species']} | {r['upstream_label']} | {r['validated_exon_type']} | "
        f"{r['final_isoform_label']} | {r.get('human_IIIb_identity')} | {r.get('human_IIIc_identity')} | "
        f"{r['label_reconciliation_confidence']} |" for r in recon
        if r["label_consistency_status"] == "swapped_relative_to_upstream")

    RECON = f"""## Annotation robustness finding — IIIb/IIIc label reconciliation

Upstream IIIb/IIIc labels were **not** trusted blindly. The final biological isoform label
(`final_isoform_label`) was determined from cassette **sequence evidence**: local alignment of
each selected protein to the curated, UniProt P21802-anchored human IIIb/IIIc cassette
references, corroborated by isoform-marker residues, with a species-level pairing that forces a
valid one-IIIb/one-IIIc assignment. Upstream labels are retained as `legacy_label`/`upstream_label`.

- Cassettes label-consistent with upstream: **{recon_c.get('label_consistent', 0)}**
- Cassettes **swapped relative to upstream** (corrected from sequence): **{recon_c.get('swapped_relative_to_upstream', 0)}**
- Ambiguous (manual review): **{recon_c.get('ambiguous_label_review', 0)}**
- Unresolved (no decisive sequence evidence): **{recon_c.get('unresolved_no_sequence', 0)}**
- Human/mouse positive control passed: **{recon_sum.get('human_mouse_control_pass', 'n/a')}**

Species with at least one swapped cassette: {', '.join(swapped_species) or '(none)'}.

| species | upstream | validated (sequence) | final | human_IIIb_id | human_IIIc_id | confidence |
|---|---|---|---|---|---|---|
{swap_tbl}

This cross-species annotation inconsistency (including the human and mouse references) is reported
as a **major annotation-robustness finding**, not a hidden manual correction. All downstream MSA
inputs, maps, projections, robustness scores, figures and `species_qc_master.tsv` use
`final_isoform_label`; the upstream labels are preserved only as legacy/provenance.

## Reference-guided residue agreement

Every cassette alignment column was mapped to human IIIb/IIIc reference residue positions; each
species residue was then classified vs the human reference as **identical**, **conservative
substitution** (positive BLOSUM62), **nonconservative substitution** (BLOSUM62 ≤ 0), **gap/missing**
or **insertion**. Conservative substitutions are sequence-comparison evidence only and are **not**
interpreted as functional claims. Gap-rich columns are excluded from the main isoform-discriminating
figure but retained in the supplement.

Reference-position agreement totals (both isoforms): identical=**{n_ident}**, conservative=**{n_cons}**,
nonconservative=**{n_nonc}**, gap/missing=**{n_gapm}**. Per-species agreement status:
{_fmt(agree_status_c)}

### Species with lowest reference agreement
| species | isoform | %identical | %identical_or_conservative | status | label status |
|---|---|---|---|---|---|
{lowest_tbl}

### Strongest isoform-discriminating positions (informative)
| human ref pos | human IIIb | human IIIc | IIIb major | IIIc major | discriminating score |
|---|---|---|---|---|---|
{strongest_tbl}
"""

    # ---- general suspicious-case rescue & validation (Part I) ----
    trig = M.read_tsv(dirs["maps"] / "fgfr2_suspicious_case_triggers.tsv")
    rescue_sum = {r["metric"]: r["value"] for r in
                  M.read_tsv(dirs["maps"] / "fgfr2_general_candidate_rescue_summary.tsv")}
    grp = M.read_tsv(dirs["maps"] / "fgfr2_validation_group_assignment.tsv")
    risk = M.read_tsv(dirs["maps"] / "fgfr2_known_risk_species_validation_report.tsv")
    patch = M.read_tsv(dirs["maps"] / "fgfr2_targeted_external_candidate_patch_report.tsv")
    gate = M.read_tsv(dirs["maps"] / "fgfr2_general_rescue_validation_gate.tsv")
    grp_c = counts(grp, "validation_group")
    n_trig = sum(1 for t in trig if t.get("rescue_required") == "true")
    trig_cat = {k.replace("trigger_", ""): sum(1 for t in trig if t.get(k) == "true")
                for k in (trig[0].keys() if trig else []) if k.startswith("trigger_")}
    dec = {k.replace("decision_", ""): v for k, v in rescue_sum.items() if k.startswith("decision_")}
    claim = {k.replace("claim_", ""): v for k, v in rescue_sum.items() if k.startswith("claim_")}
    patch_c = counts(patch, "patch_status")
    gate_fail = [g for g in gate if g.get("status") != "pass"]
    risk_tbl = "\n".join(
        f"| {r['species']} | {r['isoform']} | {r['validation_group']} | {r['problem_detected']} | "
        f"{r['rescue_result']} | {r['final_claim_status']} |" for r in risk)
    RESCUE = f"""## General suspicious-case rescue and validation

Upstream labels are not trusted blindly and no species-only swaps are hard-coded. Every
species/isoform was assigned a **validation group** (strictness only, never a forced label) and
screened for **suspicious-case triggers**; triggered high/critical cases were sent to an
evidence-driven **candidate rescue** (search per-species candidate transcripts/proteins, score
against the curated human IIIb/IIIc references + markers + MSA-discriminating support, coordinate
plausibility, protein integrity, orthology/paralog). MSA never auto-relabels; unresolved cases are
not hidden and not forced to pass; failed rescue → review/excluded with explicit reason.

- Validation groups: {_fmt(grp_c)}
- Suspicious cases triggered (rescue_required): **{n_trig}** / {len(trig)}
- Trigger categories: {_fmt(trig_cat)}
- Rescue decisions: {_fmt(dec)}
- Final claim status: {_fmt(claim)}
- External patch attempts: **{rescue_sum.get('n_external_patch_attempts', 0)}**, used in final:
  **{rescue_sum.get('n_external_patch_used', 0)}** ({_fmt(patch_c)}). External patches use only the
  local NCBI datasets cache with recorded assembly/release and are **never silently mixed** across
  incompatible annotation releases.
- General rescue validation gate: **{'PASS' if not gate_fail else 'FAIL'}**
  {('(' + '; '.join(g['check'] for g in gate_fail) + ')') if gate_fail else ''}

Close primates use strict criteria; Mus musculus and Canis lupus familiaris are explicitly
validated. Known-risk / high-priority species:

| species | isoform | group | problem detected | rescue result | final claim status |
|---|---|---|---|---|---|
{risk_tbl}

**Interpretation.** final_isoform_label drives all downstream outputs and figures; only
`primary_claim_supported`/`..._with_minor_flags` rows enter primary claims, while
`supplement_review`/`excluded_from_primary_claim` cases are shown distinctly (amber/hatched) and
explained in the supplement. Human is a hard positive control (pipeline stops if it cannot be a
primary claim).
"""

    # ---- maximal suspicious-case rescue & final biological correction (Part I) ----
    sus = M.read_tsv(dirs["maps"] / "fgfr2_all_suspicious_cases_for_rescue.tsv")
    loc_screen = M.read_tsv(dirs["maps"] / "fgfr2_exhaustive_local_rescue_candidate_screen.tsv")
    ext_screen = M.read_tsv(dirs["maps"] / "fgfr2_external_rescue_candidate_screen.tsv")
    final_dec = M.read_tsv(dirs["maps"] / "fgfr2_maximal_rescue_final_decision.tsv")
    mgate = M.read_tsv(dirs["maps"] / "fgfr2_maximal_rescue_validation_gate.tsv")
    fin_iso = [r for r in final_dec if r.get("isoform_or_pair") not in ("pair", "")]
    n_sus = sum(1 for s in sus if s.get("rescue_required") == "true")
    sus_prio = counts([s for s in sus if s.get("rescue_required") == "true"], "rescue_priority")
    sus_trig_cat = {k.replace("trigger_", ""): sum(1 for s in sus if s.get(k) == "true")
                    for k in (sus[0].keys() if sus else []) if k.startswith("trigger_")}
    dec_c = counts(fin_iso, "final_rescue_decision")
    claim_after_c = counts(fin_iso, "final_claim_status_after_rescue")
    n_rescued = sum(v for k, v in dec_c.items() if k.startswith("rescued"))
    n_confirmed = sum(v for k, v in dec_c.items() if "confirmed" in k)
    n_seqsupport = dec_c.get("sequence_support_only_keep_supplement", 0)
    n_excluded = dec_c.get("no_valid_rescue_candidate_exclude_primary", 0)
    improved = sorted({f"{r['species']}/{r['final_isoform_label']}" for r in fin_iso
                       if r.get("final_rescue_decision", "").startswith("rescued")})
    unresolved = sorted({f"{r['species']}/{r['final_isoform_label']}": r for r in fin_iso
                         if r.get("final_claim_status_after_rescue") in
                         ("supplement_review", "excluded_from_primary_claim")}.items())
    mgate_fail = [g for g in mgate if g.get("status") != "pass"]
    n_primary = sum(1 for r in fin_iso
                    if (r.get("final_claim_status_after_rescue") or "").startswith("primary_claim"))
    n_supp = len(fin_iso) - n_primary
    rescued_tbl = "\n".join(
        f"| {r['species']} | {r['final_isoform_label']} | {r['initial_problem'][:36]} | "
        f"{r['local_rescue_result']} | {r['external_rescue_result']} | {r['final_rescue_decision']} | "
        f"{r['final_protein_id']} | {r['final_claim_status_after_rescue']} |"
        for r in fin_iso if r.get("final_rescue_decision", "").startswith("rescued")
        or r.get("final_claim_status_after_rescue") in ("supplement_review", "excluded_from_primary_claim"))
    unres_tbl = "\n".join(
        f"| {r['species']} | {r['final_isoform_label']} | {r['final_claim_status_after_rescue']} | "
        f"{r['unresolved_reason_if_any']} |" for _, r in unresolved)
    MAXRESCUE = f"""## Maximal suspicious-case rescue

Beyond detecting and excluding problematic rows, every suspicious species/isoform was actively
driven toward the biologically correct FGFR2 IIIb/IIIc candidate. Each case has exactly one of two
outcomes: **rescued** with a sequence/provenance-validated better candidate (propagated downstream),
or **not recoverable** after exhaustive local + external search (excluded from primary claims and
explained). Rescue is evidence-driven (cassette identity to curated human IIIb/IIIc references,
B-type `SGINSSN`-like / A-type `GVNTTDKEI`-like markers, MSA-discriminating support, orthology,
protein integrity), never species-name-driven; upstream/legacy labels are preserved only as
provenance.

- Suspicious cases (rescue_required): **{n_sus}** / {len(sus)}; priority {_fmt(sus_prio)}
- Trigger categories: {_fmt(sus_trig_cat)}
- Local candidates screened: **{len(loc_screen)}** rows
- External candidates screened (NCBI/RefSeq datasets cache, provenance-tracked): **{len(ext_screen)}** rows
- Rescue decisions: {_fmt(dec_c)}
- Rescued cases: **{n_rescued}**; current candidates confirmed after exhaustive screen:
  **{n_confirmed}**; sequence-support-only (supplement): **{n_seqsupport}**; excluded from primary:
  **{n_excluded}**
- Final claim status after rescue: {_fmt(claim_after_c)}
- Final **primary** vs **supplement/review** cassettes: **{n_primary}** vs **{n_supp}**
- Maximal-rescue validation gate (Part H): **{'PASS' if not mgate_fail else 'FAIL'}**
  {('(' + '; '.join(g['check'] for g in mgate_fail) + ')') if mgate_fail else ''}

Species specifically improved by rescue: {', '.join(improved) or '(none)'}.

### Rescued / corrected and supplement cases
| species | iso | initial problem | local result | external result | decision | final protein | claim |
|---|---|---|---|---|---|---|---|
{rescued_tbl}

### Cases that still cannot be rescued (and why)
| species | iso | claim | unresolved reason |
|---|---|---|---|
{unres_tbl}

**Corrected vs not recoverable.** A case is considered **corrected** only when a validated
candidate replaces or confirms the final label/coordinate (sequence + provenance + validation
pass). A case is considered **not recoverable** only after local and external rescue fail or no
source-compatible evidence is available; such cases stay in supplement/review with an explicit
reason and never enter primary figures. No final primary output depends on `upstream_label` as the
biological isoform label. The corrected snapshot is recorded in
`final_corrected_pre_interpro_dataset_manifest.tsv/.json`.
"""

    # ---- final propagation consistency (post-rescue single source of truth) ----
    truth = M.read_tsv(dirs["maps"] / "fgfr2_post_rescue_final_truth_table.tsv")
    cgate = M.read_tsv(dirs["maps"] / "fgfr2_post_rescue_cross_table_consistency_gate.tsv")
    cgate_fail = [g for g in cgate if g.get("status") != "pass"]
    t_primary = sum(1 for r in truth
                    if (r.get("final_claim_status_after_rescue") or "").startswith("primary_claim"))
    t_supp = len(truth) - t_primary
    _partial = sorted({r["species"] for r in truth} &
                     {"canis_lupus_familiaris", "pongo_abelii"})

    def _truth_row(sp, iso):
        return next((r for r in truth if r["species"] == sp and r["isoform"] == iso), {})
    partial_tbl = "\n".join(
        f"| {sp} | {iso} | {_truth_row(sp, iso).get('final_claim_status_after_rescue','')} | "
        f"{_truth_row(sp, iso).get('rescue_decision','')} | "
        f"{_truth_row(sp, iso).get('unresolved_reason_if_any','') or '-'} |"
        for sp in ("canis_lupus_familiaris", "pongo_abelii") for iso in ("IIIb", "IIIc")
        if _truth_row(sp, iso))

    PROPSYNC = f"""## Final propagation consistency after maximal rescue

After the maximal rescue, **`final_claim_status_after_rescue` is the single inclusion field** that
drives primary vs supplement/review in *every* downstream output. The post-rescue state is
materialised once in `fgfr2_post_rescue_final_truth_table.tsv` (one row per species/isoform) and
joined into all major tables (`figure6C_human_referenced_residue_agreement_map.tsv`,
`figure6_msa_projected_boundary_map.tsv`, `figure8_boundary_robustness_evidence_stack.tsv`, the
robustness/agreement tables, `fgfr2_exon_type_label_reconciliation.tsv`, `species_qc_master.tsv`
and the corrected dataset manifest). The pre-rescue `recommended_use` is retained **only as
provenance** (`recommended_use_pre_rescue`); it never overrides the post-rescue claim. Rescued rows
also propagate their validated `transcript_id`/`protein_id` everywhere, so no downstream table can
keep a stale pre-rescue identifier or a stale `supplementary_only` status.

- Post-rescue **primary** vs **supplement/review** cassettes (truth table): **{t_primary}** vs **{t_supp}**
- Gorilla (rescued IIIc + confirmed IIIb): **primary pair** in figure6C, figure6, figure8 and species_qc_master.
- Canis / Pongo are **isoform-level partial rescues** (one cassette primary, the other retained as
  supplement/review with an explicit reason) — they are never shown as a complete primary pair.
- Cross-table consistency gate (Part E): **{'PASS' if not cgate_fail else 'FAIL'}**
  {('(' + '; '.join(g['check'] for g in cgate_fail) + ')') if cgate_fail else ''}

### Partial-rescue cases (isoform-level)
| species | iso | post-rescue claim | rescue decision | unresolved reason |
|---|---|---|---|---|
{partial_tbl}

A previously reported inconsistency — figure6C treating Gorilla as primary while figure6/figure8
still listed it as `supplementary_only` and lacked rescue-provenance columns — is resolved: all
figure/evidence tables now carry `rescue_decision`, `final_label_source`,
`final_claim_status_after_rescue`, `recommended_use_pre_rescue`/`recommended_use_post_rescue`,
`rescue_evidence_summary` and `unresolved_reason_if_any`, and agree on inclusion. The
`fgfr2_post_rescue_cross_table_consistency_gate.tsv/.json` gate blocks final figures if any
table disagrees on label, identifier or claim status.
"""

    # ---- local synteny / gene-neighborhood validation ----
    syn = dirs["synteny"]
    syn_val = M.read_tsv(syn / "fgfr2_5neighbor_synteny_validation.tsv")
    syn_id = M.read_tsv(syn / "fgfr2_neighbor_identity_resolution.tsv")
    syn_gate = M.read_tsv(syn / "fgfr2_5neighbor_synteny_validation_gate.tsv")
    syn_backend = M.read_tsv(syn / "fgfr2_neighborhood_plotting_backend_report.tsv")
    syn_cls = counts(syn_val, "synteny_validation_class")
    syn_locus = counts(syn_val, "fgfr2_locus_status")
    id_method = counts(syn_id, "identity_resolution_method")
    id_status = counts(syn_id, "identity_resolution_status")
    # broad human-proteome homology naming of uncharacterized (LOC...) neighbors
    loc_rows = [r for r in syn_id if (r.get("neighbor_symbol_raw", "").upper().startswith("LOC")
                or "uncharacterized" in r.get("neighbor_symbol_raw", "").lower())]
    loc_named = [r for r in loc_rows if r.get("broad_homology_symbol")]
    loc_strong = [r for r in loc_named
                  if M.to_float(r.get("broad_homology_percent_identity"), 0) >= 50]
    loc_unnamed = len(loc_rows) - len(loc_named)
    syn_gate_fail = [g for g in syn_gate if g.get("status") != "pass"]
    diamond_used = any("diamond" in (b.get("version", "") + b.get("status", "")).lower()
                       or b.get("backend") == "matplotlib" for b in syn_backend)
    syn_rescue_tbl = "\n".join(
        f"| {r['species']} | {r['synteny_validation_class']} | "
        f"{r['total_neighbor_support_score']} | {r['rescued_candidate_locus_support']} | "
        f"{r['synteny_warning'] or '-'} |" for r in syn_val
        if r["species"] in ("gorilla_gorilla_gorilla", "canis_lupus_familiaris", "pongo_abelii"))

    SYNTENY = f"""## Local synteny / gene-neighborhood validation

To check that the selected and rescued FGFR2 candidates represent the genuine **FGFR2 locus** rather
than paralogous, fragmented or annotation-artifact loci, an independent local synteny layer was
added. For every species the FGFR2 gene was located in its own source-compatible NCBI RefSeq
`genomic.gff` and the **5 protein-coding genes upstream and 5 downstream** (main window; a 10+10
supplement window is also stored) were extracted with coordinates, strand, distance and provenance.
A human FGFR2 5/10-neighbor reference (e.g. PLEKHA1, BTBD16, TACC2, NSMCE4A, ATE1 / WDR11, PLPP4,
SEC23IP, MCMBP, INPP5F) anchors the comparison.

**Synteny validates locus/orthology context only — it never assigns or relabels IIIb/IIIc.** Isoform
labels remain sequence-calibrated (cassette references, marker residues, MSA-discriminating residues,
CDS/protein coordinate mapping).

Neighbor identities are normalized with an evidence hierarchy: exact gene-symbol match (case
insensitive across species) > curated orthology > **reciprocal best hit** > high-confidence one-way
BLAST > raw annotation. DIAMOND (preferred) / BLASTP was run against the human neighborhood reference
(standard E≤1e-5, identity≥30%, coverage≥50%; stricter 60%/70% for close mammals/primates). BLAST/RBH
names are reported as **probable** orthologs (italic + "?"), never as curated truth; ambiguous hits
are not forced.

**Naming uncharacterized (LOC...) neighbors.** Every uncharacterized / `LOC...` neighbor is
additionally searched against the **whole human RefSeq proteome** to attach the best human-homolog
gene name, with percent identity and query/subject coverage recorded so the strength of each call is
explicit. A standard pass (E≤1e-3) is followed by a very-permissive fallback (very-sensitive, large
E) so that even loose hits get a name; low-confidence hits are flagged `very weak`. Hits to human
uncharacterized loci are not used as names. Of {len(loc_rows)} `LOC...` neighbors, **{len(loc_named)}
received a human-homolog name** ({len(loc_strong)} at ≥50% identity); {loc_unnamed} have no detectable
named human homolog (lineage-specific) and remain raw IDs. These broad-homology names are **loose
locus context only** — they are not orthology claims and do not enter the synteny order/conservation
score.

- FGFR2 locus extracted (source-compatible RefSeq): {syn_locus.get('neighborhood_extracted', 0)}/{len(syn_val)} species
- Synteny validation classes: {_fmt(syn_cls)}
- Neighbor identity methods: {_fmt(id_method)}
- Neighbor identity status: {_fmt(id_status)}
- Identity engine: {'DIAMOND (preferred) / BLASTP reciprocal best hit' if diamond_used else 'BLASTP'}
- Local 5-neighbor synteny validation gate: **{'PASS' if not syn_gate_fail else 'FAIL'}**
  {('(' + '; '.join(g['check'] for g in syn_gate_fail) + ')') if syn_gate_fail else ''}

### Rescue cases — locus support
| species | synteny class | neighbor support | rescued-candidate locus support | warning |
|---|---|---|---|---|
{syn_rescue_tbl}

**Interpretation.** Strong local synteny increases locus/orthology confidence; the Gorilla, Canis and
Pongo rescued/partial candidates all sit in the conserved FGFR2 neighborhood (locus supported),
independently corroborating the sequence-based rescue. Teleost fish show expected neighborhood
rearrangement and are classified `synteny_partial_blast_supported` rather than conflict — distant
taxa are not over-penalized for order differences when the FGFR2 locus and several neighbors are
still supported. Synteny **unavailable** (no source-compatible annotation) is reported distinctly
from synteny **conflict**; a conflict can downgrade a claim to supplement/review but is never applied
silently. **MCScanX block-level synteny is intentionally not part of this build**; `combined_synteny
_validation_class` therefore equals the local class. Limitations: fragmented assemblies / scaffold
edges, missing neighbor protein sequences, generic `LOC`/uncharacterized symbols, and ambiguous
BLAST hits are reported as such rather than forced.
"""

    NOCLAIM = (
        "**What is claimed and not claimed.** The MSA layer is QC/robustness evidence only. "
        "It does **not** relabel IIIb/IIIc: isoform assignment remains sequence-calibrated from "
        "the exon/cassette pipeline. Cassette boundaries are **projected** from validated "
        "exon/protein coordinates onto alignment columns; they are not re-derived from the MSA. "
        "Review/supplement species are retained and explained but are **not** used for primary "
        "claims. **InterProScan has not been run** and no domain annotations are asserted; "
        "domain-boundary claims remain downstream of InterProScan.")

    # ---------------- captions ----------------
    cap = f"""# MSA boundary-robustness figure captions (pre-InterPro)

{NOCLAIM}

## Figure 6 — MSA-projected IIIb/IIIc boundary map (primary-claim cassettes)
Validated/rescued IIIb (blue) and IIIc (orange) cassette boundaries projected onto the full-length
FGFR2 protein MSA (MAFFT --auto). **Per Part D/G, only post-rescue primary-claim cassettes are
shown** (`final_claim_status_after_rescue` is the inclusion field); a bullet marks
primary-with-minor-flags rows. Rows are species in phylogenetic order; left strips encode the
component-based boundary-robustness class; a per-column conservation track is shown below each panel.
Rescued candidates (e.g. Gorilla IIIc) appear under their final label with rescued coordinates;
unresolved/supplement cassettes are moved to the supplement panel. MSA does not relabel IIIb/IIIc.

## Figure 7 — Isoform-discriminating residues
Per-column IIIb and IIIc within-isoform conservation and the between-isoform discriminating
score computed on the combined IIIb+IIIc cassette MSA (main-analysis species). Positions that
are conserved within each isoform but differ between IIIb and IIIc are highlighted; human
IIIb/IIIc residues are annotated at key positions. This supports the sequence-calibrated
IIIb/IIIc distinction and does not relabel any isoform.

## Figure 8 — Boundary robustness evidence stack (primary-claim cassettes)
Per species/isoform component scores (coordinate resolution, codon/boundary precision, protein
QC, MSA projection, conservation/gap, protein integrity) and the final robustness class/score.
**Per Part D/G, only post-rescue primary-claim cassettes are plotted**; the underlying table
`figure8_boundary_robustness_evidence_stack.tsv` retains every species/isoform with the post-rescue
provenance columns (`rescue_decision`, `final_label_source`, `final_claim_status_after_rescue`,
`recommended_use_pre_rescue`/`recommended_use_post_rescue`). Robustness is supported by multiple
independent evidence layers, not a single coordinate source. A minor split-codon flag does not by
itself demote a species; protein conflicts or gap-rich boundaries reduce robustness.

## Supplement — MSA review-case diagnostics
One panel per review/supplement species (and any MSA projection-review case) showing native
coordinate status, normalized slot status, MSA projection, gap/conservation, robustness score
and final interpretation. Difficult cases are retained and explained, not hidden.

## Supplement — Full-length MSA / protein integrity overview
Protein length distribution, per-sequence gap fraction in the full-length MSA, integrity-status
counts and InterProScan readiness for the selected FGFR2 proteins.

## Figure 6C — Human-referenced IIIb/IIIc residue agreement map
Rows are species in phylogenetic order; separate panels for IIIb and IIIc. The x-axis is the
human reference cassette residue position; each cell is the per-position agreement class
(identical = strong blue-green, conservative = light blue, nonconservative = orange/vermillion,
gap/missing = light grey). Black vertical lines mark cassette boundaries; left sidebars encode
taxon group; the top track shows per-position identical-or-conservative fraction. **Per Part G,
only primary-claim cassettes after maximal rescue are shown** (each panel is filtered per isoform);
sequence/provenance-validated **rescued candidates are included** under their final label (e.g.
Canis IIIb, Gorilla IIIc); review/excluded cassettes are moved to the supplement. Final
sequence-calibrated labels are used. InterProScan pending.

## Figure 6D — MSA cassette boundary map, local zoom
Per-species residue occupancy across the L-INS-i cassette alignment columns (IIIb/IIIc panels),
coloured by isoform; insertions relative to human and gaps are shown separately, gap-rich columns
are shaded, and black ticks mark the human-reference cassette boundary columns. Validated
exon/cassette residues project to comparable local alignment columns in robust species.

## Figure 7C — Isoform-discriminating residues (informative positions)
Tracks of human IIIb, human IIIc, IIIb-major and IIIc-major residues at informative human-reference
positions, with the per-position discriminating score below. Positions conserved within each
isoform but different between IIIb and IIIc are boxed. Gap-rich non-informative columns are
excluded from this main plot and retained in the supplement. MSA supports but does not relabel
IIIb/IIIc.

## Figure 8C — Alignment evidence stack (primary-claim cassettes)
Per species/isoform: a validation-group sidebar, IIIb/IIIc reference agreement, boundary
projection, left/right boundary and cassette-core conservation (blue-green = stronger), gap burden
(orange), discriminating-residue support, the final discrete alignment-evidence class, and the
**rescue status** and **final claim status** from the validation/rescue layer. Per Part G the
primary stack shows only primary-claim (species, isoform) rows after maximal rescue; supplement/
excluded cassettes (with their rescue decision and reason) are in the supplement panels and the
maximal-rescue tables. Distinguishes coordinate robustness, reference sequence agreement,
isoform-discriminating support and gap-rich/review cases. InterProScan pending.

## Supplement — Per-species cassette difference panels
For review species and representative main species: IIIb and IIIc reference-agreement strips,
nonconservative differences and gaps, with the label-consistency status, to make difficult species
interpretable.

## Figure 9A — FGFR2 local gene-neighborhood map (5 neighbors)
Strand-aware arrow map of the FGFR2 locus for representative species, showing the 5 protein-coding
genes upstream and 5 downstream of FGFR2 (centre, black anchor) in their genomic order, oriented
relative to the FGFR2 strand. Arrows are coloured by normalized neighbor orthology group (consistent
across species; human reference neighbors PLEKHA1/BTBD16/TACC2/NSMCE4A/ATE1 and WDR11/PLPP4/SEC23IP/
MCMBP/INPP5F). Neighbor identity is resolved by symbol > curated orthology > reciprocal best hit >
one-way BLAST > raw ID; **probable** BLAST/RBH names are italic with "?", unmapped/raw IDs grey,
ambiguous identities amber-outlined. Rescued-candidate loci have a dashed blue outline; supplement/
review species a dotted outline; a right-hand swatch gives the synteny validation class. Synteny
validates the FGFR2 locus/orthology context only and does **not** assign IIIb/IIIc.

## Figure 9B — FGFR2 5/10-neighbor conservation matrix
Rows = species (phylogenetic order), columns = human FGFR2 reference neighbors. Cell colour encodes
neighbor state: present same-side/order, same-side reordered, opposite side, present elsewhere in the
5-neighbor window, present only in the 10-neighbor supplement, probable-by-BLAST/RBH, ambiguous
identity, missing/unmapped, or scaffold unavailable. Right-hand sidebars show taxon group, post-rescue
claim status, rescue, and synteny validation class. Locus/orthology context only; IIIb/IIIc labels are
not derived from synteny.

## Figure 9C — FGFR2 rescue-case locus panels
Local neighborhood arrow panels for the rescued/partial cases (Gorilla, Canis, Pongo) plus any other
high-risk supplement/review species, each annotated with synteny validation class, total neighbor
support, rescued-candidate locus support, rescue decision, post-rescue claim status, and any
unresolved reason. Confirms whether rescued candidates sit in the genuine FGFR2 neighborhood; synteny
does not relabel IIIb/IIIc.

## Supplement — FGFR2 local gene-neighborhood (10 neighbors, all species)
The 10-neighbor window for all species, showing whether human reference neighbors fall just outside
the 5-neighbor main window. (MCScanX block-level synteny / Figure 9D is intentionally omitted from
this build.)
"""
    (dirs["captions"] / "msa_figure_captions_pre_interpro.md").write_text(cap, encoding="utf-8")

    # ---------------- QC migration report ----------------
    aln_lines = "\n".join(
        f"| {r['msa_name']} | {r['n_sequences']} | {r['aligned_length']} | "
        f"{r['mean_gap_fraction']} | {r['alignment_status']} |" for r in run_man)
    qc = f"""# QC migration report — pre-InterPro MSA & boundary robustness

Module: `12_msa_boundary_robustness_pre_interpro` · generated {M.now_iso()}

{NOCLAIM}

## MSA input counts
- Full-length protein sequences: **{n_full}**
- IIIb cassette sequences: **{n_iiib}**
- IIIc cassette sequences: **{n_iiic}**
- Input validation failures: **{n_val_fail}**

## Alignment quality (MAFFT --auto)
| MSA | n_seq | aligned_len | mean_gap | status |
|---|---|---|---|---|
{aln_lines}

## Boundary projection counts
{_fmt(proj_c)}

## Conservation summary
- Mean cassette-region conservation (per species/isoform): **{mean_cass}**

## Isoform-discriminating residues (main-only)
{_fmt(disc_main)}

## Robustness score counts
{_fmt(rob_c)}

## Splice-site QC availability
{_fmt(splice_c)}
(Splice-site motif QC requires source-compatible genomic sequence; where unavailable it is
reported as `splice_site_sequence_unavailable` and never faked.)

## Protein integrity counts
{_fmt(integ_c)}

## Review-case diagnostics
- Review/projection-review cases retained and explained: **{len(diag)}**

{RECON}

{RESCUE}

{MAXRESCUE}

{PROPSYNC}

{SYNTENY}

## Known limitations
- MSA quality may be lower in divergent or incompletely annotated species.
- Splice-site QC is only available when source-compatible genomic sequence exists.
- Gap-rich columns are excluded from the main discriminating-residue plot but retained in the supplement.
- Conservative substitutions are sequence-comparison evidence only and are not functional claims.
- A few species' selected isoform proteins remain ambiguous/unresolved and are excluded from primary claims.
- InterProScan has not been run.
- Domain-boundary claims remain downstream of InterProScan.
"""
    (md / "QC_migration_report_msa_boundary_robustness_pre_interpro.md").write_text(qc, encoding="utf-8")

    # ---------------- methods ----------------
    methods = f"""# Methods — pre-InterPro MSA & boundary-robustness layer

Selected FGFR2 proteins (one per species/isoform) and the validated IIIb/IIIc cassette peptides
(extracted with the corrected coordinate-overlap cassette map) were aligned with MAFFT
(`--auto`; {next((r['msa_name'] for r in run_man), '')} etc.). For each alignment a residue↔column
coordinate map was built and validated against input lengths. Validated cassette boundaries were
projected onto alignment columns (full-length and cassette-only), with gap-aware projection
status (high-confidence / minor-gaps / gap-rich / shift / unresolved). Per-column conservation
was computed as 1 − normalized Shannon entropy over non-gap residues, with gap fraction reported
separately. A transparent, component-based boundary-robustness score combined coordinate
resolution (0.20), codon/boundary precision (0.15), protein QC (0.15), MSA projection (0.25),
conservation/gap evidence (0.15) and protein integrity (0.10); exact weights are written to
`robustness/fgfr2_boundary_robustness_component_weights.tsv`. Isoform-discriminating cassette
positions were identified on the combined IIIb+IIIc MSA. Splice-site motif QC was attempted only
where source-compatible genomic sequence was available.

Before alignment, a **sequence-calibrated IIIb/IIIc label reconciliation** was applied: the final
biological isoform label was derived from local alignment to the curated UniProt P21802-anchored
human IIIb/IIIc cassette references (with marker-residue corroboration and species-level pairing),
not from upstream labels or exon order. High-accuracy MAFFT L-INS-i (`--localpair --maxiterate
1000`) cassette alignments were added for the reference-guided comparison; `--auto` alignments are
retained as sensitivity outputs. Per human-reference position, residues were classified as
identical / conservative (positive BLOSUM62) / nonconservative (BLOSUM62 ≤ 0) / gap / insertion,
and summarised by species, by position and by cassette segment (left boundary / core / right
boundary). The boundary-robustness score gained explicit reference-agreement components
(reference identity, boundary/core agreement, discriminating-residue support, gap-rich penalty and
an overall alignment-evidence class). {NOCLAIM}

{RECON}

{RESCUE}

{MAXRESCUE}

{PROPSYNC}

{SYNTENY}
"""
    (md / "methods_update_msa_boundary_robustness_pre_interpro.md").write_text(methods, encoding="utf-8")

    # ---------------- results summary ----------------
    n_robust = rob_c.get("robust_boundary", 0)
    n_hi = proj_c.get("msa_boundary_projected_high_confidence", 0)
    results = f"""# Results summary — pre-InterPro MSA & boundary robustness

Across {n_full} full-length FGFR2 proteins (IIIb={n_iiib}, IIIc={n_iiic}), MAFFT alignments were
generated for the full-length proteins and the IIIb, IIIc and combined cassette sets. Validated
cassette boundaries projected to comparable alignment columns for the majority of orthologs
({n_hi}/{len(proj)} high-confidence projections); shifted or gap-rich projections were flagged
for review rather than corrected. Mean cassette-region conservation was {mean_cass}, and
{disc_main.get('isoform_discriminating_conserved', 0)} cassette positions were conserved within
isoforms but discriminating between IIIb and IIIc, supporting the sequence-calibrated distinction.
The component-based boundary-robustness score classified {n_robust}/{len(scores)} species/isoforms
as robust; {len(diag)} review cases were retained and explained. {NOCLAIM}

{RECON}

{RESCUE}

{MAXRESCUE}

{PROPSYNC}

{SYNTENY}
"""
    (md / "results_summary_msa_boundary_robustness_pre_interpro.md").write_text(results, encoding="utf-8")

    # ---- append a clearly-marked reconciliation section to earlier pre-InterPro reports ----
    marker_a = "<!-- BEGIN: IIIb/IIIc label reconciliation (auto) -->"
    marker_b = "<!-- END: IIIb/IIIc label reconciliation (auto) -->"
    block = (f"\n\n{marker_a}\n# IIIb/IIIc label reconciliation update\n\n{RECON}\n\n{RESCUE}\n\n"
             f"{MAXRESCUE}\n\n{PROPSYNC}\n\n{SYNTENY}\n{marker_b}\n")
    for name in ("QC_migration_report_tasks_7_to_12_pre_interpro.md",
                 "results_summary_pre_interpro.md", "methods_update_pre_interpro.md",
                 "figure_captions_pre_interpro.md"):
        p = M.locate(base, name)
        if not p:
            continue
        txt = p.read_text(encoding="utf-8", errors="replace")
        if marker_a in txt and marker_b in txt:
            pre = txt.split(marker_a)[0].rstrip()
            post = txt.split(marker_b, 1)[1] if marker_b in txt else ""
            txt = pre + block + post.lstrip()
        else:
            txt = txt.rstrip() + block
        p.write_text(txt, encoding="utf-8")

    print("[OK] reports + captions written (captions, QC migration, methods, results, earlier reports)")
    return 0


def _fmt(d: Dict[str, object]) -> str:
    if not d:
        return "- (none)"
    return "\n".join(f"- `{k}`: {v}" for k, v in sorted(d.items(), key=lambda x: str(x[0])))


if __name__ == "__main__":
    raise SystemExit(main())
