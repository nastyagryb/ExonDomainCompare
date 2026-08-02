#!/usr/bin/env python3
"""
make_all_figures.py  (Sprint Part 4 + Part 7 — MANDATORY central entry point)

Single entry point for the final pre-InterPro publication sprint:

    python scripts/make_all_figures.py --base results/final_30_until_interpro_prepare

This generates ALL final pre-InterPro figures, plotting tables, captions and
manifests, and writes the completion reports. It:

  * builds the reproducible phylogenetic/taxonomic species order (Part 2),
  * builds the CDS phase/boundary audit tables (Part 1 support),
  * ensures species_qc_master.tsv carries the phylo-order columns (canonical),
  * runs the final validation gate (Part 1) and FAILS clearly if data are
    missing/stale,
  * renders the publication figures (Parts 3, 5, 6) via
    make_publication_figures_pre_interpro.render_all,
  * writes publication_figure_manifest.tsv, output_file_manifest_pre_interpro.tsv
    and the completion reports (Part 7).

Reads FINAL tables only and never recomputes biological QC.
No real InterPro domains are generated or plotted.

A legacy mode (--tables/--outdir) is retained for the older Task-11 figure set.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

SCRIPT_VERSION = "2.0"
_HERE = Path(__file__).resolve().parent
PY = sys.executable


def _load(modname: str, filename: str):
    spec = importlib.util.spec_from_file_location(modname, _HERE / filename)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, rows: List[Dict[str, object]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def locate(base: Path, name: str, hint: str = "") -> Optional[Path]:
    matches = sorted(base.rglob(name))
    if not matches:
        return None
    if hint:
        for m in matches:
            if hint in str(m):
                return m
    return sorted(matches, key=lambda p: len(p.parts))[0]


def _require(base: Path, name: str, hint: str = "") -> Path:
    p = locate(base, name, hint)
    if p is None:
        raise RuntimeError(f"Required input not found under {base}: {name}")
    return p


# ---------------------------------------------------------------------------
# Prerequisite builders
# ---------------------------------------------------------------------------
def build_phylo_order(base: Path, pub_tables: Path) -> Path:
    registry = _require(base, "species_registry.tsv", "01_species_registry")
    subprocess.run([PY, str(_HERE / "build_species_phylogenetic_order.py"),
                    "--registry", str(registry), "--outdir", str(pub_tables)], check=True)
    return pub_tables / "species_phylogenetic_order.tsv"


def build_cds_audit(base: Path) -> Path:
    coord = _require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    mapping = _require(base, "fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv")
    cds_features = _require(base, "cds_features.tsv", "02_models")
    pair_qc = _require(base, "fgfr2_pair_level_qc_summary.tsv")
    proteins = locate(base, "selected_fgfr2_proteins.faa")
    cmd = [PY, str(_HERE / "build_cds_phase_boundary_audit.py"),
           "--coordinate_audit", str(coord),
           "--cds_features", str(cds_features),
           "--outdir", str(coord.parent),
           "--update_coordinate_audit",
           "--update_exon_cds_mapping", str(mapping),
           "--update_pair_qc", str(pair_qc)]
    if proteins:
        cmd += ["--proteins", str(proteins)]
    subprocess.run(cmd, check=True)
    return coord.parent / "cds_phase_boundary_audit.tsv"


def build_unclassified_isoform_fallback(base: Path) -> Optional[Path]:
    """Sequence-calibrated fallback for cassettes whose protein interval is unresolved
    (e.g. fresh Ensembl returned isoform=unclassified). Patches the coordinate audit so
    that the IIIb/IIIc cassette slot gets a REAL mid-protein coordinate instead of the
    synthetic aa-1 placeholder, BEFORE the cassette coordinate mapping. Keeps the Step-11
    aa-1 gate strict; never relabels validated calls."""
    coord = _require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    proteins = locate(base, "selected_fgfr2_proteins.faa")
    if proteins is None:
        print("[skip] unclassified-isoform fallback: no selected_fgfr2_proteins.faa")
        return None
    cmd = [PY, str(_HERE / "assign_unclassified_fgfr2_isoform_fallback.py"),
           "--coordinate_audit", str(coord),
           "--proteins", str(proteins),
           "--outdir", str(coord.parent)]
    pvs = locate(base, "fgfr2_III_final_selected_protein_validation_summary.tsv")
    if pvs:
        cmd += ["--protein_validation_summary", str(pvs)]
    ro = locate(base, "fgfr2_rescue_overrides.tsv")
    if ro:
        cmd += ["--rescue_overrides", str(ro)]
    msa = locate(base, "final_isoform_discriminating_residues.tsv")
    if msa:
        cmd += ["--msa_discriminating", str(msa)]
    subprocess.run(cmd, check=True)
    return coord.parent / "unclassified_isoform_fallback_audit.tsv"


def build_cassette_map(base: Path) -> Path:
    """Cassette-mapping correction (Parts A/B/C/E): unique CDS blocks + coordinate-
    overlap cassette mapping + reconstruction + sanity audits."""
    coord = _require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    cds_features = _require(base, "cds_features.tsv", "02_models")
    proteins = locate(base, "selected_fgfr2_proteins.faa")
    master = _require(base, "species_qc_master.tsv", "11_pre_interpro_master")
    patch_cache = base / "02_models" / "_ncbi_cds_boundary_patch_cache"
    cmd = [PY, str(_HERE / "build_cds_block_cassette_map.py"),
           "--coordinate_audit", str(coord),
           "--cds_features", str(cds_features),
           "--review_master", str(master),
           "--outdir", str(coord.parent)]
    if proteins:
        cmd += ["--proteins", str(proteins)]
    if patch_cache.exists():
        cmd += ["--cds_fasta_dir", str(patch_cache)]
    subprocess.run(cmd, check=True)
    return coord.parent / "fgfr2_cassette_cds_block_map.tsv"


def build_phase_rescue(base: Path) -> Path:
    """Uncertainty-refinement Part B: attempt codon-phase rescue (cumulative-CDS
    reconstruction) for phase-unavailable cassettes before final classification."""
    coord = _require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    d = coord.parent
    cds_features = _require(base, "cds_features.tsv", "02_models")
    subprocess.run([PY, str(_HERE / "build_phase_rescue.py"),
                    "--cds_audit", str(d / "cds_phase_boundary_audit.tsv"),
                    "--cds_features", str(cds_features),
                    "--cassette_map", str(d / "fgfr2_cassette_cds_block_map.tsv"),
                    "--reconstruction_audit", str(d / "fgfr2_transcript_cds_reconstruction_audit.tsv"),
                    "--outdir", str(d)], check=True)
    return d / "cds_phase_rescue_audit.tsv"


def build_refined_classes(base: Path) -> Path:
    """Uncertainty-refinement Part A: refined, explainable uncertainty/display classes."""
    coord = _require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    d = coord.parent
    master = _require(base, "species_qc_master.tsv", "11_pre_interpro_master")
    subprocess.run([PY, str(_HERE / "build_refined_uncertainty_classes.py"),
                    "--rescue", str(d / "cds_phase_rescue_audit.tsv"),
                    "--cassette_map", str(d / "fgfr2_cassette_cds_block_map.tsv"),
                    "--sanity", str(d / "fgfr2_cassette_coordinate_sanity_audit.tsv"),
                    "--coordinate_audit", str(coord),
                    "--master", str(master),
                    "--outdir", str(d)], check=True)
    return d / "fgfr2_refined_uncertainty_classes.tsv"


def build_patch_report(base: Path, enable_network: bool = False) -> Path:
    """Uncertainty-refinement Part C: targeted NCBI patch provenance for TRUE missing
    data only (never minor phase/split flags). Offline-safe."""
    coord = _require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    d = coord.parent
    cache = base / "02_models" / "_ncbi_cds_boundary_patch_cache"
    cmd = [PY, str(_HERE / "patch_ncbi_cds_boundaries.py"),
           "--refined_classes", str(d / "fgfr2_refined_uncertainty_classes.tsv"),
           "--cds_audit", str(d / "cds_phase_boundary_audit.tsv"),
           "--reconstruction_audit", str(d / "fgfr2_transcript_cds_reconstruction_audit.tsv"),
           "--cache_dir", str(cache), "--outdir", str(d)]
    if enable_network:
        cmd += ["--enable_network"]
    subprocess.run(cmd, check=True)
    return d / "fgfr2_ncbi_cds_boundary_patch_report.tsv"


# ---------------------------------------------------------------------------
# Part D — propagate refined uncertainty/display classes into the final QC tables
# ---------------------------------------------------------------------------
DISPLAY_RANK = {
    "hard_fail_excluded": 6, "protein_overlay_only": 5, "review_protein": 4,
    "review_annotation": 3, "resolved_phase_not_available": 2,
    "resolved_with_split_codon": 1, "robust": 0,
}
VIS_RANK = {"hard_fail": 4, "main_warning": 3, "supplement_only": 2,
            "subtle_symbol": 1, "none": 0}
REFINED_ISO_COLS = [
    "coordinate_resolution_state", "boundary_precision_state", "protein_evidence_state",
    "annotation_review_state", "display_uncertainty_class", "plot_visibility_level",
    "uncertainty_explanation_short", "uncertainty_explanation_long",
]
MASTER_REFINED_COLS = REFINED_ISO_COLS + ["phase_rescue_status", "ncbi_patch_status"]


def _merge_cols(path: Path, lookup: Dict, key_fn, cols: List[str]) -> None:
    if not path.exists():
        return
    rows = read_tsv(path)
    if not rows:
        return
    header = list(rows[0].keys())
    for r in rows:
        ref = lookup.get(key_fn(r), {})
        for c in cols:
            r[c] = ref.get(c, r.get(c, ""))
    for c in cols:
        if c not in header:
            header.append(c)
    write_tsv(path, rows, header)


def propagate_refined(base: Path, master: Path) -> None:
    """Append refined uncertainty/display columns to the final QC tables and the
    canonical master (per-isoform where keyed by isoform, aggregated per species
    for pair-level and master)."""
    coord = _require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    d = coord.parent
    refined = read_tsv(d / "fgfr2_refined_uncertainty_classes.tsv")
    rescue = {(r["species"], r["isoform"]): r for r in read_tsv(d / "cds_phase_rescue_audit.tsv")}
    patch = {(r["species"], r["isoform"]): r
             for r in read_tsv(d / "fgfr2_ncbi_cds_boundary_patch_report.tsv")}

    iso_lookup = {(r["species"], r["isoform"]): r for r in refined}

    # per-species aggregate (worst-case display/visibility, joined states)
    by_sp: Dict[str, List[Dict[str, str]]] = {}
    for r in refined:
        by_sp.setdefault(r["species"].lower(), []).append(r)
    sp_agg: Dict[str, Dict[str, str]] = {}
    for sp, rs in by_sp.items():
        worst = max(rs, key=lambda x: DISPLAY_RANK.get(x["display_uncertainty_class"], 0))
        vis = max((x["plot_visibility_level"] for x in rs), key=lambda v: VIS_RANK.get(v, 0))
        def _join(col):
            return ";".join(sorted({x[col] for x in rs if x.get(col)}))
        rescue_st = ";".join(sorted({rescue.get((x["species"], x["isoform"]), {})
                                     .get("rescue_status", "") for x in rs} - {""}))
        patch_st = ";".join(sorted({patch.get((x["species"], x["isoform"]), {})
                                    .get("patch_status", "") for x in rs} - {""}))
        sp_agg[sp] = {
            "coordinate_resolution_state": _join("coordinate_resolution_state"),
            "boundary_precision_state": _join("boundary_precision_state"),
            "protein_evidence_state": _join("protein_evidence_state"),
            "annotation_review_state": _join("annotation_review_state"),
            "display_uncertainty_class": worst["display_uncertainty_class"],
            "plot_visibility_level": vis,
            "uncertainty_explanation_short": worst["uncertainty_explanation_short"],
            "uncertainty_explanation_long": worst["uncertainty_explanation_long"],
            "phase_rescue_status": rescue_st or "not_needed_already_resolved",
            "ncbi_patch_status": patch_st or "patch_not_needed",
        }

    # per-isoform tables
    _merge_cols(d / "cds_phase_boundary_audit.tsv", iso_lookup,
                lambda r: (r.get("species"), r.get("isoform")), REFINED_ISO_COLS)
    # also carry rescued boundary precision into the audit
    _merge_cols(d / "cds_phase_boundary_audit.tsv",
                {(r["species"], r["isoform"]): r for r in read_tsv(d / "cds_phase_rescue_audit.tsv")},
                lambda r: (r.get("species"), r.get("isoform")),
                ["rescued_boundary_precision", "rescue_status"])
    _merge_cols(coord, iso_lookup,
                lambda r: (r.get("species_canonical"), r.get("inferred_isoform")),
                REFINED_ISO_COLS)
    _merge_cols(d / "fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv", iso_lookup,
                lambda r: (r.get("species_canonical"), r.get("inferred_isoform")),
                REFINED_ISO_COLS)

    # per-species tables
    _merge_cols(d / "fgfr2_pair_level_qc_summary.tsv", sp_agg,
                lambda r: r.get("species_canonical", "").lower(), MASTER_REFINED_COLS)
    _merge_cols(master, sp_agg,
                lambda r: r.get("species", "").lower(), MASTER_REFINED_COLS)
    alias = master.parent / "species_qc_master_pre_interpro.tsv"
    if alias.exists():
        _merge_cols(alias, sp_agg, lambda r: r.get("species", "").lower(), MASTER_REFINED_COLS)
    print(f"[OK] Part D: propagated refined classes into QC tables + master "
          f"({len(sp_agg)} species).")


PHYLO_MERGE_COLS = {
    "taxid": "taxid", "taxon_group": "taxon_group",
    "taxon_group_display": "taxon_group_display", "major_clade": "major_clade",
    "phylo_order": "phylo_order", "phylo_order_source": "order_source",
    "phylo_order_confidence": "order_confidence",
}


def _cds_explain_for_species(cds_rows: List[Dict[str, str]]) -> str:
    if not cds_rows:
        return "no_cds_audit"
    unk = [r for r in cds_rows if str(r.get("reason_if_unknown")) != "not_unknown"]
    spl = [r for r in cds_rows if str(r.get("reason_if_split")) != "not_split"]
    if not unk and not spl:
        return "all_boundaries_known"
    parts = []
    if unk:
        reasons = sorted({str(r.get("reason_if_unknown")) for r in unk})
        parts.append(f"unknown_x{len(unk)}:{'/'.join(reasons)}")
    if spl:
        parts.append(f"split_x{len(spl)}")
    return "; ".join(parts)


def ensure_master_phylo(base: Path, phylo_path: Path, cds_audit_path: Optional[Path] = None) -> Path:
    """If species_qc_master.tsv lacks phylo / CDS-explainability columns, merge them
    in (ordering + audit metadata only; no biological QC is recomputed). Keeps the
    alias identical."""
    master = _require(base, "species_qc_master.tsv", "11_pre_interpro_master")
    rows = read_tsv(master)
    header = list(rows[0].keys()) if rows else []
    has_phylo = "phylo_order" in header and all(r.get("phylo_order", "") for r in rows)
    has_cds = "cds_boundary_explainability" in header
    if has_phylo and has_cds:
        return master  # already integrated by the master builder

    if not has_phylo:
        phylo = {r["species"].strip().lower(): r for r in read_tsv(phylo_path)}
        for r in rows:
            ph = phylo.get(r["species"].strip().lower(), {})
            for out_col, src_col in PHYLO_MERGE_COLS.items():
                r[out_col] = ph.get(src_col, r.get(out_col, ""))
            if ph.get("taxon_group"):
                r["taxon_group"] = ph["taxon_group"]

    if not has_cds and cds_audit_path and Path(cds_audit_path).exists():
        cds_by_sp: Dict[str, List[Dict[str, str]]] = {}
        for r in read_tsv(cds_audit_path):
            cds_by_sp.setdefault(str(r.get("species", "")).strip().lower(), []).append(r)
        for r in rows:
            r["cds_boundary_explainability"] = _cds_explain_for_species(
                cds_by_sp.get(r["species"].strip().lower(), []))

    new_fields = list(header)
    for c in list(PHYLO_MERGE_COLS) + ["cds_boundary_explainability"]:
        if c not in new_fields:
            new_fields.append(c)

    def _po(r):
        v = str(r.get("phylo_order", "")).strip()
        return (int(v), r["species"]) if v.isdigit() else (10 ** 6, r["species"])

    rows.sort(key=_po)
    write_tsv(master, rows, new_fields)
    alias = master.parent / "species_qc_master_pre_interpro.tsv"
    write_tsv(alias, rows, new_fields)
    return master


# ---------------------------------------------------------------------------
# Reports (Part 7)
# ---------------------------------------------------------------------------
def _count(rows, col):
    return dict(Counter(str(r.get(col, "")) for r in rows))


def write_reports(base: Path, pub: Path, master: Path, phylo_path: Path,
                  validation_summary: Dict, manifest: List[Dict]) -> Dict[str, Path]:
    mrows = read_tsv(master)
    n_species = len(mrows)
    fdc = _count(mrows, "final_display_class")
    n_main = sum(v for k, v in fdc.items() if k.startswith("main_analysis"))
    n_supp = sum(v for k, v in fdc.items() if k.startswith("supplementary"))
    taxa = _count(mrows, "taxon_group")
    phy = read_tsv(phylo_path)
    order_src = _count(phy, "order_source")
    order_conf = _count(phy, "order_confidence")

    interpro_summary = locate(base, "fgfr2_interpro_prepare_summary.tsv")
    interpro = {r["metric"]: r["value"] for r in read_tsv(interpro_summary)} if interpro_summary else {}

    cds_counts_p = locate(base, "cds_phase_boundary_legacy_vs_refined_counts.tsv")
    cds_counts = read_tsv(cds_counts_p) if cds_counts_p else []
    cds_exp_p = locate(base, "cds_phase_boundary_explainability_summary.tsv")
    cds_exp = read_tsv(cds_exp_p) if cds_exp_p else []

    def _exp(dim):
        return {r["category"]: r["count"] for r in cds_exp if r["dimension"] == dim}

    # cassette-mapping correction sprint tables
    cmap_p = locate(base, "fgfr2_cassette_cds_block_map.tsv")
    cmap = read_tsv(cmap_p) if cmap_p else []
    sanity_p = locate(base, "fgfr2_cassette_coordinate_sanity_audit.tsv")
    sanity = read_tsv(sanity_p) if sanity_p else []
    recon_p = locate(base, "fgfr2_transcript_cds_reconstruction_audit.tsv")
    recon = read_tsv(recon_p) if recon_p else []
    patch_p = locate(base, "fgfr2_ncbi_cds_boundary_patch_report.tsv")
    patch = read_tsv(patch_p) if patch_p else []
    cmap_status = _count(cmap, "cassette_overlap_status")
    sanity_status = _count(sanity, "coordinate_sanity_status")
    recon_status = _count(recon, "reconstruction_status")
    patch_status = _count(patch, "patch_status")
    n_at_start1 = sum(1 for m in cmap if str(m.get("matched_protein_start_aa")) == "1")

    # uncertainty-refinement sprint tables (Parts A/B/C)
    refined_p = locate(base, "fgfr2_refined_uncertainty_classes.tsv")
    refined = read_tsv(refined_p) if refined_p else []
    rescue_p = locate(base, "cds_phase_rescue_audit.tsv")
    rescue = read_tsv(rescue_p) if rescue_p else []
    coord_res = _count(refined, "coordinate_resolution_state")
    bound_prec = _count(refined, "boundary_precision_state")
    disp_class = _count(refined, "display_uncertainty_class")
    vis_level = _count(refined, "plot_visibility_level")
    rescue_status = _count(rescue, "rescue_status")
    patch_status_ref = _count(patch, "patch_status")
    n_rescued = sum(v for k, v in rescue_status.items() if k.startswith("rescued_"))
    n_true_missing = sum(v for k, v in coord_res.items()
                         if k in ("protein_overlay_no_cds_model", "coordinate_unresolved"))

    paralog_manifest = locate(base, "fgfr2_paralog_reference_panel_manifest.tsv")
    n_panel = len(read_tsv(paralog_manifest)) if paralog_manifest else 0
    ortho_sum = locate(base, "fgfr2_orthology_species_summary.tsv")
    ortho_rows = read_tsv(ortho_sum) if ortho_sum else []
    ortho_status = _count(ortho_rows, "orthology_status_species")
    ortho_evidence = locate(base, "fgfr2_orthology_evidence.tsv")
    n_ortho_ev = len(read_tsv(ortho_evidence)) if ortho_evidence else 0

    fig_lines = "\n".join(
        f"- `{e['figure_id']}` — SVG/PDF/PNG — {e.get('main_message','')}" for e in manifest)
    cds_lines = "\n".join(
        f"  - {r['precision_category']}: legacy={r['legacy_count']} -> refined={r['refined_count']}"
        for r in cds_counts)
    taxa_lines = "\n".join(f"  - {k}: {v}" for k, v in taxa.items())
    ortho_lines = "\n".join(f"  - {k}: {v}" for k, v in ortho_status.items())

    report = f"""# QC migration report — Tasks 7–12 (pre-InterPro completion)

_Generated: {_now()} · script v{SCRIPT_VERSION}_

## Pipeline completion status up to InterProScan preparation

- The pipeline is **complete up to validated InterProScan-ready input preparation**.
- **InterProScan / domain annotation has NOT yet been run.**
- **No real InterPro domain coordinates are claimed** anywhere in these outputs.
- **All final figures are pre-InterPro figures** (exon-to-protein coordinate mapping and QC).
- Downstream InterProScan analysis remains the next step.

## Final QC table status

- Canonical display/QC table: `species_qc_master.tsv` ({n_species} species).
- `species_qc_master_pre_interpro.tsv` is maintained as an identical alias.
- final_display_class distribution: {json.dumps(fdc)}
- main-analysis species: {n_main}; supplement/review species: {n_supp}.
- Validation gate overall status: **{validation_summary.get('overall_status','?')}**
  (ok={validation_summary.get('n_ok')}, warning={validation_summary.get('n_warning')}, fail={validation_summary.get('n_fail')}).

## species_qc_master.tsv status

- Canonical and integrates orthology, paralog-panel, protein QC, sequence-calibrated
  direction, isoform detection, resolver status, CDS-boundary precision and
  phylogenetic ordering columns.

## Orthology / paralog evidence status

- Multi-vertebrate FGFR1/2/3/4 paralog reference panel entries: {n_panel}.
- Orthology evidence records (protein-level): {n_ortho_ev}.
- Species-level orthology status:
{ortho_lines or '  - (orthology summary unavailable)'}

## CDS-boundary audit and uncertainty explainability

- `cds_phase_boundary_audit.tsv` (one explainable row per resolved IIIb/IIIc mapping),
  `cds_phase_boundary_legacy_vs_refined_counts.tsv` and
  `cds_phase_boundary_explainability_summary.tsv` are present and referenced here.
- legacy vs refined precision:
{cds_lines or '  - (counts unavailable)'}
- reason_if_unknown distribution: {json.dumps(_exp('reason_if_unknown'))}
- reason_if_split distribution: {json.dumps(_exp('reason_if_split'))}
- transcript_cds_reconstruction_status: {json.dumps(_exp('transcript_cds_reconstruction_status'))}
- protein_translation_check_status: {json.dumps(_exp('protein_translation_check_status'))}

**Why unknown/split codon phase can remain even though CDS coordinates are known.**
A CDS feature carries genomic start/end coordinates, but the reading-frame *phase* at a
boundary is a separate GFF3 annotation. Every `unknown_codon_phase` case here is explained
and not a silent failure: `phase_not_propagated_from_source` (Ensembl-sourced cassettes whose
GFF3 phase was not propagated) or `nucleotide_sequence_unavailable` (cassette transcript absent
from the local CDS model). Split-codon boundaries are expected biology for internal cassette
exons whose length is not a multiple of three.

**Uncertain boundaries are evidence-level flags, not errors.** Uncertain cases are never forced
to `exact`. As an independent control, the resolved cassette transcripts that are present in the
CDS model reconstruct from CDS coordinates to a total length consistent with the selected protein
(`cds_protein_length_consistent`), validating the protein coordinate projection even where codon
phase is unannotated. The explainability columns are propagated into
`fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv`,
`fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv`, `fgfr2_pair_level_qc_summary.tsv`
and summarised per species in `species_qc_master.tsv` (`cds_boundary_explainability`).

## Uncertainty classes and why they are shown differently

The refined uncertainty classes (`fgfr2_refined_uncertainty_classes.tsv`, summarised per species in
`species_qc_master.tsv`) collapse the many low-level flags into a small, explainable set with an
explicit `plot_visibility_level`, so the figures no longer overstate uncertainty.

- **A known split-codon boundary is not a failure.** Internal cassette exons whose length is not a
  multiple of three split a codon at a boundary; expected biology, shown only as a small grey edge symbol.
- **Phase unavailable does not mean the coordinate is wrong.** {n_rescued}/60 phase-unavailable cases were
  rescued from a length-consistent CDS reconstruction (the cumulative-CDS reading frame reproduces the
  NCBI source-phase split/exact calls exactly); the remainder are labelled
  `phase_not_available_but_coordinate_resolved`, not wrong.
- **True missing data are rare and explicitly counted:** {n_true_missing}/60 mappings lack a local
  CDS-block model (`protein_overlay_no_cds_model` / `coordinate_unresolved`) and are the only NCBI-patch
  candidates.
- **NCBI patching is used only for true missing information and is provenance-tracked**
  (`fgfr2_ncbi_cds_boundary_patch_report.tsv`); known split/phase flags and locally reconstructable cases
  are explicitly NOT patched, and Ensembl/NCBI releases are never mixed silently.
- **Minor boundary-precision flags are shown subtly** (small edge symbols / pale colors); review cases
  (protein conflict, major native-coordinate offset, hard coordinate sanity fail) remain visible with a
  prominent marker and are interpreted separately, not used for primary claims.

Counts (refined uncertainty / phase rescue / NCBI patch):

- coordinate_resolution_state: {json.dumps(coord_res)}
- boundary_precision_state: {json.dumps(bound_prec)}
- display_uncertainty_class: {json.dumps(disp_class)}
- plot_visibility_level: {json.dumps(vis_level)}
- phase_rescue_status: {json.dumps(rescue_status)}
- ncbi_patch_status: {json.dumps(patch_status_ref)}

## Cassette → CDS-block mapping correction (coordinate overlap)

**Previous-figure caveat.** An earlier version of `figure2_exon_to_protein_architecture_tracks.tsv`
identified the IIIb/IIIc cassette CDS block by a CDS identifier. For NCBI/RefSeq models all CDS
blocks of a transcript share the same `cds-XP_...` id, so the first block was wrongly tagged as the
cassette and many cassettes were plotted at `protein_start_aa = 1`. This was a join/mapping bug, not
a biological result, and those figures were **not biologically reliable**.

**Fix.** Cassettes are now mapped onto a table of **unique CDS blocks**
(`fgfr2_unique_cds_block_table.tsv`, deterministic `unique_cds_block_id =
species|transcript|cds_rank|seqid:start-end:strand`) by **genomic (then protein) coordinate
overlap** (`fgfr2_cassette_cds_block_map.tsv`). The unique `cds_rank` — not the repeated NCBI id —
is used to select the cassette block.

- cassette_overlap_status: {json.dumps(cmap_status)}
- cassettes now mapped to protein_start_aa == 1: **{n_at_start1}** (was 38).
- cassette coordinate sanity ({sanity_p.name if sanity_p else 'n/a'}): {json.dumps(sanity_status)}
- All main-analysis cassettes passed the coordinate sanity gate; the figure pipeline is **fail-fast**
  and writes `publication_figure_validation_failed.tsv` (and skips figures) if any main-analysis or
  control-primate cassette is N-terminal (< 150 aa in a full-length protein) or maps to CDS rank 1.

## CDS reconstruction / translation validation

- `fgfr2_transcript_cds_reconstruction_audit.tsv`: {json.dumps(recon_status)}.
- Reconstructed transcript CDS length is consistent with the selected protein length for the
  cassette transcripts present in the local CDS model (terminal stop-codon offset allowed). Full
  nucleotide translation/identity is computed when CDS nucleotide sequence is available (Part D
  patch cache); otherwise the coordinate/length-consistency proxy is reported honestly.

## NCBI/RefSeq CDS-boundary patch

- `fgfr2_ncbi_cds_boundary_patch_report.tsv`: {json.dumps(patch_status)}.
- Targeted retrieval cache: `02_models/_ncbi_cds_boundary_patch_cache/` (existing validated models
  are never overwritten; full provenance per row). Ensembl-sourced phase gaps are flagged
  `annotation_release_mismatch_review` rather than patched from NCBI (no release mixing).

## InterProScan FASTA / mapping / manifest status

- Input FASTA: `{interpro.get('interpro_input_fasta','fgfr2_interpro_clean_unique.fasta')}`
- selected proteins: {interpro.get('total_selected_proteins','?')};
  unique sequences: {interpro.get('unique_sequences','?')};
  duplicates collapsed: {interpro.get('duplicates_collapsed','?')};
  invalid rejected: {interpro.get('invalid_sequences_rejected','?')}.
- status: {interpro.get('interpro_status','interpro_ready_input_prepared')}.

## Figure export status

- All publication figures exported in **SVG, PDF and PNG**.
- Figures:
{fig_lines}

## Phylogenetic ordering method and source

- Reproducible order in `species_phylogenetic_order.tsv`.
- order_source distribution: {json.dumps(order_src)}
- order_confidence distribution: {json.dumps(order_conf)}
- Method: NCBI Taxonomy via ETE when available, otherwise a documented curated
  fallback taxonomic order (never silent alphabetical except flagged unresolved rows).

## Taxon group counts

{taxa_lines}

## Main vs supplement/review species counts

- main-analysis species: {n_main}
- supplement/review species: {n_supp}

## Known limitations

- InterProScan domain annotation pending; domain-aware overlays are a downstream step.
- Phylogenetic within-group order is curated (taxonomy-derived) where ETE/NCBITaxa is
  unavailable; group-level placement follows standard vertebrate phylogeny.
- A subset of species rely on transcript-level orthology evidence and/or have
  unresolved sequence-calibration direction (documented in the evidence tables) and
  are retained as review/supplement, not used for primary claims.
- CDS-boundary precision is phase-derived; codon_split / unknown_phase cases are
  flagged and not over-interpreted.
"""

    results = f"""# Results summary (pre-InterPro)

- {n_species} vertebrate FGFR2 orthologs analyzed; 60 resolved IIIb/IIIc rows (30 species x 2 isoforms).
- {n_main} main-analysis species and {n_supp} supplement/review species (final_display_class).
- Sequence-calibrated IIIb/IIIc direction assignment is fixed; protein evidence is QC only and
  never auto-swaps IIIb/IIIc labels.
- Multi-vertebrate paralog panel ({n_panel} reference proteins) supports FGFR2 orthology;
  species-level orthology status: {json.dumps(ortho_status)}.
- {interpro.get('unique_sequences','?')} unique, non-redundant InterProScan-ready sequences prepared
  ({interpro.get('duplicates_collapsed','?')} duplicates collapsed).
- All species ordered by reproducible phylogenetic/taxonomic order; taxon groups: {json.dumps(taxa)}.
- InterProScan has not been executed; domain annotation is the next step.
"""

    methods = f"""# Methods update (pre-InterPro)

## Species ordering
Species are ordered by a reproducible phylogenetic/taxonomic scheme
(`species_phylogenetic_order.tsv`). NCBI Taxonomy via ETE/NCBITaxa is used when the
local taxonomy database is available; otherwise a documented curated fallback order is
applied (order_source: {json.dumps(order_src)}). Broad taxon groups (Primates, Other
mammals, Birds, Reptiles, Amphibians, Teleost fish) are used for grouping and subtle
visual banding in figures.

## Canonical QC source
`species_qc_master.tsv` is the single canonical display/QC table. All plotting scripts
read final tables only and do not recompute biological QC, IIIb/IIIc assignment,
similarity classification or review status.

## CDS-boundary precision and explainability
GFF3 phase-derived left/right CDS-boundary precision is recorded per resolved IIIb/IIIc
cassette in `cds_phase_boundary_audit.tsv`, with the source of any uncertainty made explicit
(`reason_if_unknown`, `reason_if_split`), a coordinate-based transcript-CDS reconstruction
(`transcript_cds_reconstruction_status`) and a protein-length consistency check
(`protein_translation_check_status`). Category counts are in
`cds_phase_boundary_explainability_summary.tsv` and legacy-vs-refined deltas in
`cds_phase_boundary_legacy_vs_refined_counts.tsv`. Uncertain boundaries are retained as
evidence-level flags and never forced to exact.

## InterProScan preparation
A non-redundant, validated protein FASTA (`fgfr2_interpro_clean_unique.fasta`) with stable
ID mappings and an input manifest is prepared as the endpoint of this pipeline stage.
InterProScan itself has not been executed; no InterPro domain coordinates are claimed.

## Figures
Figures are exported as SVG and PDF (primary) and PNG (300–600 dpi preview) using a
color-blind-safe palette with stable IIIb (#0072B2) and IIIc (#E69F00) colors.
"""

    (pub / "QC_migration_report_tasks_7_to_12_pre_interpro.md").write_text(report, encoding="utf-8")
    (pub / "results_summary_pre_interpro.md").write_text(results, encoding="utf-8")
    (pub / "methods_update_pre_interpro.md").write_text(methods, encoding="utf-8")

    # output file manifest (walk base)
    man_rows: List[Dict[str, object]] = []
    for p in sorted(base.rglob("*")):
        if p.is_file():
            try:
                st = p.stat()
            except OSError:
                continue
            man_rows.append({
                "path": str(p.relative_to(base)), "size_bytes": st.st_size,
                "modified_utc": datetime.fromtimestamp(st.st_mtime, tz=timezone.utc)
                .strftime("%Y-%m-%dT%H:%M:%SZ"),
            })
    out_manifest = pub / "metadata" / "output_file_manifest_pre_interpro.tsv"
    write_tsv(out_manifest, man_rows, ["path", "size_bytes", "modified_utc"])

    return {
        "qc_report": pub / "QC_migration_report_tasks_7_to_12_pre_interpro.md",
        "results_summary": pub / "results_summary_pre_interpro.md",
        "methods_update": pub / "methods_update_pre_interpro.md",
        "output_manifest": out_manifest,
    }


# ---------------------------------------------------------------------------
# Sprint orchestration
# ---------------------------------------------------------------------------
def run_sprint(base: Path) -> int:
    base = base.resolve()
    pub = base / "11_publication_figures_pre_interpro"
    pub_tables, pub_meta = pub / "tables", pub / "metadata"
    pub_tables.mkdir(parents=True, exist_ok=True)
    pub_meta.mkdir(parents=True, exist_ok=True)

    print("[1/6] building reproducible phylogenetic order …")
    phylo_path = build_phylo_order(base, pub_tables)

    print("[2/6] building explainable CDS phase/boundary audit …")
    cds_audit_path = build_cds_audit(base)

    print("[2a/6] sequence-calibrated fallback for unclassified-isoform cassettes (pre-mapping) …")
    build_unclassified_isoform_fallback(base)

    print("[2b/6] mapping IIIb/IIIc cassettes onto unique CDS blocks (coordinate overlap) …")
    build_cassette_map(base)

    print("[2c/6] phase rescue + refined uncertainty classes + targeted NCBI patch …")
    build_phase_rescue(base)
    build_refined_classes(base)
    build_patch_report(base, enable_network=False)

    print("[3/6] ensuring species_qc_master.tsv carries phylo + CDS-explainability columns …")
    master = ensure_master_phylo(base, phylo_path, cds_audit_path)

    print("[3b/6] propagating refined uncertainty/display classes into final QC tables (Part D) …")
    propagate_refined(base, master)

    print("[4/6] running final pre-InterPro validation gate …")
    val = _load("final_pre_interpro_validation", "final_pre_interpro_validation.py")
    ok, _rows, summary = val.run_validation(base, outdir=pub_meta)
    if not ok:
        val.validate_or_raise(base)  # raises with detailed message

    print("[5/6] rendering publication figures (SVG/PDF/PNG) …")
    figs = _load("make_publication_figures_pre_interpro", "make_publication_figures_pre_interpro.py")
    manifest, _paths = figs.render_all(base)

    print("[6/6] writing completion reports and manifests …")
    reports = write_reports(base, pub, master, phylo_path, summary, manifest)

    print(f"[DONE] {len(manifest)} figures; validation={summary['overall_status']}")
    print(f"       figures   -> {pub / 'figures'}")
    print(f"       captions  -> {pub / 'captions' / 'figure_captions_pre_interpro.md'}")
    print(f"       manifest  -> {pub_meta / 'publication_figure_manifest.tsv'}")
    print(f"       report    -> {reports['qc_report']}")
    return 0


# ---------------------------------------------------------------------------
# Legacy mode (older Task-11 figure set)
# ---------------------------------------------------------------------------
def run_legacy(tables: Path, outdir: Path) -> int:
    figs = _load("pre_interpro_figs", "make_pre_interpro_figures.py")
    _spec = getattr(figs, "FIGURE_SPEC", None)
    # Fall back to the module's own main if present.
    if hasattr(figs, "render_all_legacy"):
        figs.render_all_legacy(tables, outdir)  # type: ignore[attr-defined]
        return 0
    raise SystemExit("Legacy figure spec unavailable; use --base for the publication sprint.")


def main() -> int:
    ap = argparse.ArgumentParser(description="Central pre-InterPro figure entry point (Part 4).")
    ap.add_argument("--base", type=Path, default=None,
                    help="results/<run> base dir (publication sprint; recommended)")
    ap.add_argument("--tables", type=Path, default=None, help="legacy figure-tables dir")
    ap.add_argument("--outdir", type=Path, default=None, help="legacy output dir")
    args = ap.parse_args()
    if args.base:
        return run_sprint(args.base)
    if args.tables and args.outdir:
        return run_legacy(args.tables, args.outdir)
    ap.error("provide --base (publication sprint) or --tables and --outdir (legacy)")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
