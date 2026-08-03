#!/usr/bin/env python3
"""
Build refined uncertainty classes.

Collapses the many low-level QC flags into a small set of explainable uncertainty
classes and a plot-visibility level, so figures stop overstating uncertainty:
split codons and phase-unavailable-but-coordinate-resolved cases are NOT errors and
must be shown subtly, while true missing data / conflicts / hard fails stay prominent.

Reads final tables only (no biological QC recomputed):
  cds_phase_rescue_audit.tsv, fgfr2_cassette_cds_block_map.tsv,
  fgfr2_cassette_coordinate_sanity_audit.tsv, fgfr2_current_stage_..._coordinate_audit.tsv,
  species_qc_master.tsv

Output:
  fgfr2_refined_uncertainty_classes.tsv  (one row per resolved IIIb/IIIc mapping)
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List


COLS = [
    "species", "isoform", "transcript_id", "protein_id",
    "coordinate_resolution_state", "boundary_precision_state", "protein_evidence_state",
    "annotation_review_state", "display_uncertainty_class", "plot_visibility_level",
    "uncertainty_explanation_short", "uncertainty_explanation_long",
]


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path or not Path(path).exists():
        return []
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, rows, fields) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def coord_state(overlap_status: str) -> str:
    s = overlap_status or ""
    if "exact_overlap" in s or "partial_overlap" in s:
        return "cds_coordinate_resolved"
    if "overlay_from_resolved_protein_interval" in s:
        return "protein_overlay_no_cds_model"
    if "conflict" in s:
        return "coordinate_conflict_review"
    return "coordinate_unresolved"


def boundary_state(rescued: str) -> str:
    r = (rescued or "").lower()
    if "exact" in r:
        return "codon_boundary_exact"
    if "split" in r:
        return "known_split_codon_boundary"
    if "phase_not_available" in r:
        return "phase_not_available_but_coordinate_resolved"
    if "nucleotide_sequence_unavailable" in r:
        return "nucleotide_sequence_unavailable"
    return "boundary_unresolved"


def protein_state(summary: str) -> str:
    s = (summary or "").lower()
    if "conflict" in s:
        return "protein_conflict_review"
    if "ambiguous" in s or "below_threshold" in s:
        return "protein_ambiguous_review"
    if "transcript_only" in s and "validated" not in s:
        return "protein_transcript_level_only"
    if "validated" in s:
        return "protein_validated"
    if not s:
        return "protein_unavailable"
    return "protein_validated"


def annotation_state(sanity_status: str, native_sanity: str, iii_sim: str):
    """Return (state, severity) where severity in {hard, major, minor, none}."""
    if ("implausible" in (sanity_status or "") or "first_cds_block" in (sanity_status or "")):
        return "hard_coordinate_fail_review", "hard"
    ns = native_sanity or ""
    if "major_native_offset" in ns:
        return "native_coordinate_offset_review", "major"
    if "ambiguous_similarity_review" in (iii_sim or ""):
        return "iii_region_similarity_review", "major"
    if "moderate_native_offset" in ns:
        return "native_coordinate_offset_review", "minor"
    return "no_annotation_review", "none"


EXPL = {
    "robust": ("coordinates resolved; codon boundaries exact",
               "CDS coordinates are resolved and both cassette boundaries fall exactly on codon "
               "boundaries; no uncertainty to display."),
    "resolved_with_split_codon": ("coordinate-resolved; cassette boundary splits a codon (normal, not an error)",
               "The cassette is coordinate-resolved; one or both boundaries split a codon. Split codons "
               "are expected for internal cassette exons and are NOT errors; shown as a small edge symbol."),
    "resolved_phase_not_available": ("coordinate-resolved; codon phase not inferable (not a wrong coordinate)",
               "The cassette is coordinate-resolved but codon phase could not be inferred from available "
               "annotation/sequence. Phase-unavailable does NOT mean the coordinate is wrong; shown subtly."),
    "protein_overlay_only": ("true missing data: no local CDS-block model; protein-coordinate overlay",
               "No per-exon CDS-block model is available for this transcript in the local model, so the "
               "cassette is drawn as a protein-coordinate overlay. This is genuine missing data and is "
               "flagged prominently / handled by targeted NCBI patching with provenance."),
    "review_protein": ("protein-evidence review (conflict/ambiguous); interpreted separately",
               "Protein marker evidence is in conflict or ambiguous; the IIIb/IIIc assignment is retained "
               "(protein QC never auto-swaps labels) but this species is reviewed separately and not used "
               "for primary claims."),
    "review_annotation": ("annotation review (coordinate offset / III-region similarity)",
               "An annotation-level review flag applies (native-coordinate offset or III-region similarity). "
               "The species remains visible and is interpreted separately from primary claims."),
    "hard_fail_excluded": ("hard cassette-coordinate sanity failure; excluded from primary claims",
               "The mapped cassette failed a hard coordinate sanity check (e.g. N-terminal in a full-length "
               "protein). It is kept visible for transparency but excluded from primary claims."),
}


def main() -> int:
    ap = argparse.ArgumentParser(description="Refined uncertainty / display classes (Part A).")
    ap.add_argument("--rescue", type=Path, required=True)
    ap.add_argument("--cassette_map", type=Path, required=True)
    ap.add_argument("--sanity", type=Path, required=True)
    ap.add_argument("--coordinate_audit", type=Path, required=True)
    ap.add_argument("--master", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    rescue = {(r["species"], r["isoform"]): r for r in read_tsv(args.rescue)}
    cmap = {(r["species"], r["isoform"]): r for r in read_tsv(args.cassette_map)}
    sanity = {(r["species"], r["isoform"]): r for r in read_tsv(args.sanity)}
    master = {r["species"].lower(): r for r in read_tsv(args.master)}

    out = []
    for c in read_tsv(args.coordinate_audit):
        sp = (c.get("species_canonical") or "").lower()
        iso = c.get("inferred_isoform", "")
        key = (sp, iso)
        rc = rescue.get(key, {})
        m = cmap.get(key, {})
        sa = sanity.get(key, {})
        mr = master.get(sp, {})

        fdc = mr.get("final_display_class", "")
        is_supp = fdc.startswith("supplementary")
        _is_main = fdc.startswith("main_analysis")

        cstate = coord_state(m.get("cassette_overlap_status", ""))
        bstate = boundary_state(rc.get("rescued_boundary_precision",
                                        c.get("cds_boundary_precision_refined", "")))
        pstate = protein_state(mr.get("protein_validation_summary", ""))
        astate, severity = annotation_state(sa.get("coordinate_sanity_status", ""),
                                            mr.get("native_coordinate_sanity", ""),
                                            mr.get("iii_region_similarity_class", ""))

        # display class — priority of biological importance
        if severity == "hard":
            disp = "hard_fail_excluded"
        elif cstate in ("coordinate_conflict_review", "coordinate_unresolved"):
            disp = "review_annotation"
        elif cstate == "protein_overlay_no_cds_model" or bstate == "nucleotide_sequence_unavailable":
            disp = "protein_overlay_only"
        elif pstate == "protein_conflict_review":
            disp = "review_protein"
        elif severity == "major":
            disp = "review_annotation"
        elif pstate == "protein_ambiguous_review":
            disp = "review_protein"
        elif severity == "minor":
            disp = "review_annotation"
        elif bstate == "codon_boundary_exact":
            disp = "robust"
        elif bstate == "known_split_codon_boundary":
            disp = "resolved_with_split_codon"
        elif bstate == "phase_not_available_but_coordinate_resolved":
            disp = "resolved_phase_not_available"
        else:
            disp = "review_annotation"

        # plot visibility — keep minor technical flags subtle; reserve prominence for
        # true missing data / conflicts / major offsets / hard fails
        if disp == "hard_fail_excluded":
            vis = "hard_fail"
        elif disp == "protein_overlay_only":
            vis = "supplement_only" if is_supp else "main_warning"
        elif disp == "review_protein":
            vis = ("supplement_only" if is_supp else
                   "main_warning" if pstate == "protein_conflict_review" else "subtle_symbol")
        elif disp == "review_annotation":
            if severity == "major":
                vis = "supplement_only" if is_supp else "main_warning"
            elif cstate in ("coordinate_conflict_review", "coordinate_unresolved"):
                vis = "supplement_only" if is_supp else "main_warning"
            else:  # minor offset etc.
                vis = "supplement_only" if is_supp else "subtle_symbol"
        elif disp == "resolved_with_split_codon":
            vis = "subtle_symbol"
        elif disp == "resolved_phase_not_available":
            vis = "subtle_symbol"
        else:  # robust
            vis = "none"

        short, long = EXPL.get(disp, ("", ""))
        out.append({
            "species": c.get("species_canonical", ""), "isoform": iso,
            "transcript_id": c.get("transcript_id_source", ""), "protein_id": c.get("protein_id", ""),
            "coordinate_resolution_state": cstate, "boundary_precision_state": bstate,
            "protein_evidence_state": pstate, "annotation_review_state": astate,
            "display_uncertainty_class": disp, "plot_visibility_level": vis,
            "uncertainty_explanation_short": short, "uncertainty_explanation_long": long,
        })

    write_tsv(args.outdir / "fgfr2_refined_uncertainty_classes.tsv", out, COLS)
    from collections import Counter
    print(f"[OK] fgfr2_refined_uncertainty_classes.tsv rows={len(out)}")
    print(f"     display_uncertainty_class={dict(Counter(r['display_uncertainty_class'] for r in out))}")
    print(f"     plot_visibility_level={dict(Counter(r['plot_visibility_level'] for r in out))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
