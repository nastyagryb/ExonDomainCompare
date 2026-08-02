#!/usr/bin/env python3
"""
make_pre_interpro_report_package.py  (Task 12)

Assemble the pre-InterProScan report and reproducibility package from the
already-written pipeline outputs. Generates:

  QC_migration_report_tasks_1_to_12_pre_interpro.md
  methods_update_pre_interpro.md
  results_summary_pre_interpro.md
  figure_captions_pre_interpro.md
  run_manifest_pre_interpro.json
  output_file_manifest_pre_interpro.tsv

The package explicitly states that InterProScan / domain annotation has NOT been
executed; no fake domain coordinates are produced anywhere in this pipeline.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

SCRIPT_VERSION = "1.0"


def read_tsv(path: Optional[Path]) -> List[Dict[str, str]]:
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return []
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def read_json(path: Optional[Path]) -> dict:
    if not path or not Path(path).exists():
        return {}
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return {}


def counts_block(d: Dict[str, int]) -> str:
    if not d:
        return "- (none)\n"
    return "".join(f"- `{k}`: {v}\n" for k, v in sorted(d.items(), key=lambda kv: (-int(kv[1]), kv[0])))


def main() -> int:
    ap = argparse.ArgumentParser(description="Build the pre-InterPro report + reproducibility package (Task 12).")
    ap.add_argument("--base", type=Path, required=True, help="Pipeline BASE result directory.")
    ap.add_argument("--master", type=Path, required=True)
    ap.add_argument("--resolver_metadata", type=Path, required=True)
    ap.add_argument("--interpro_summary", type=Path, required=True)
    ap.add_argument("--interpro_dir", type=Path, required=True)
    ap.add_argument("--figures_dir", type=Path, required=True)
    ap.add_argument("--tables_dir", type=Path, required=True)
    ap.add_argument("--paralog_panel_manifest", type=Path, default=None, help="Addendum A panel manifest")
    ap.add_argument("--paralog_screen_metadata", type=Path, default=None, help="Addendum A multi-vertebrate screen metadata json")
    ap.add_argument("--orthology_metadata", type=Path, default=None, help="Addendum B orthology evidence metadata json")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()

    master = read_tsv(args.master)
    resolver_meta = read_json(args.resolver_metadata)
    cds = resolver_meta.get("cds_boundary_precision", {})
    interpro_summary = {r["metric"]: r["value"] for r in read_tsv(args.interpro_summary)}
    panel_manifest = read_tsv(args.paralog_panel_manifest)
    paralog_meta = read_json(args.paralog_screen_metadata)
    orthology_meta = read_json(args.orthology_metadata)

    display_counts = Counter(r.get("final_display_class", "") for r in master)
    eligible_counts = Counter(r.get("main_analysis_eligible", "") for r in master)
    interpro_counts = Counter(r.get("interpro_status", "") for r in master)
    review_species = [r for r in master if r.get("final_display_class") != "main_analysis_high_confidence"]

    interpro_cmd = (
        "interproscan.sh \\\n"
        "  -i fgfr2_interpro_clean_unique.fasta \\\n"
        "  -f TSV,GFF3,JSON \\\n"
        "  -appl Pfam,SMART,PROSITE,PRINTS,CDD \\\n"
        "  -goterms -pa -iprlookup \\\n"
        "  -cpu 4 \\\n"
        "  -o fgfr2_interproscan_results"
    )

    # ---- main migration report ----
    refined_after = {k: int(v) for k, v in cds.get("cds_boundary_precision_refined_counts_after", {}).items()}
    legacy_before = {k: int(v) for k, v in cds.get("legacy_cds_boundary_precision_counts_before", {}).items()}
    resolver_refined = {k: int(v) for k, v in resolver_meta.get("resolver_status_refined_counts", {}).items()}

    md = []
    md.append("# FGFR2 IIIb/IIIc pipeline QC migration report — Tasks 1–12 (pre-InterProScan)\n")
    md.append(f"_Generated: {now}_\n")
    md.append("> **Scope:** This pipeline currently ends at **InterProScan-ready input preparation**. "
              "InterProScan / protein-domain annotation has **NOT** been executed. No InterPro domain "
              "coordinates or domain overlays are produced anywhere in these outputs.\n")

    md.append("## Steps completed up to InterProScan preparation\n")
    md.append(
        "- Steps 1–3: species registry, model collection (dual-source NCBI/Ensembl), initial selection.\n"
        "- Step 4: IIIb/IIIc isoform evidence with **sequence-calibrated direction assignment** "
        "(replaces the provisional exon-order rule; human controls pass; systematic inversion eliminated).\n"
        "- Step 5/5b: role-aware candidate selection with **protein evidence used as QC only** "
        "(never as automatic isoform reassignment); mandatory reference + control checks.\n"
        "- Step 6/6b/6c/6d/6e: paralog screen, III-pair audit, III-region anchoring and final pair audit.\n"
        "- Step 9/10: coordinate resolver with native + normalized III-slot coordinates and "
        "**left/right CDS-boundary precision** (Task 7); pair-level QC.\n"
        "- Step 7: clean, non-redundant **InterProScan-ready FASTA** + mapping/summary/manifest (Task 9).\n"
        "- Task 8/10/11/12: single-source-of-truth species QC master, figure tables, pre-InterPro figures, this package.\n")

    md.append("\n## Task 7 — CDS-boundary precision (before / after)\n")
    md.append("Legacy precision only inspected the transcript-5' (start-phase) boundary. The refined rule "
              "evaluates **both** the left (5') and right (3') CDS boundaries in transcript/translation order "
              "(negative strand respected via GFF3 phase), deriving the right boundary from start phase + CDS length.\n")
    md.append("\n**Before (`legacy_cds_boundary_precision`):**\n")
    md.append(counts_block(legacy_before))
    md.append("\n**After (`cds_boundary_precision_refined`):**\n")
    md.append(counts_block(refined_after))
    md.append("\n**Left boundary:** \n")
    md.append(counts_block({k: int(v) for k, v in cds.get("cds_left_boundary_precision_counts", {}).items()}))
    md.append("\n**Right boundary:** \n")
    md.append(counts_block({k: int(v) for k, v in cds.get("cds_right_boundary_precision_counts", {}).items()}))
    md.append("\n**Phase source:** \n")
    md.append(counts_block({k: int(v) for k, v in cds.get("cds_phase_source_counts", {}).items()}))

    md.append("\n## resolver_status_refined counts\n")
    md.append(counts_block(resolver_refined))

    md.append("\n## main_analysis_eligible counts (species QC master)\n")
    md.append(counts_block(dict(eligible_counts)))
    md.append("\n## final_display_class counts\n")
    md.append(counts_block(dict(display_counts)))

    md.append("\n## Review / supplementary species\n")
    if review_species:
        md.append("| species | display class | review reason |\n|---|---|---|\n")
        for r in review_species:
            md.append(f"| {r.get('species','')} | {r.get('final_display_class','')} | {r.get('review_reason_long','')} |\n")
    else:
        md.append("- None — all species are main-analysis high-confidence.\n")

    md.append("\n## InterPro FASTA preparation summary\n")
    md.append(f"- Total selected proteins: **{interpro_summary.get('total_selected_proteins','?')}**\n")
    md.append(f"- Unique protein sequences: **{interpro_summary.get('unique_sequences','?')}**\n")
    md.append(f"- Duplicates collapsed: **{interpro_summary.get('duplicates_collapsed','?')}**\n")
    md.append(f"- Invalid sequences rejected: **{interpro_summary.get('invalid_sequences_rejected','?')}**\n")
    md.append(f"- Species covered: **{interpro_summary.get('species_covered','?')}** "
              f"(both isoforms: {interpro_summary.get('species_with_both_isoforms','?')})\n")
    md.append(f"- InterPro input status: **{interpro_summary.get('interpro_status','?')}**\n")
    md.append(f"- Per-species `interpro_status`: " + ", ".join(f"{k}={v}" for k, v in interpro_counts.items()) + "\n")

    md.append("\n## Multi-vertebrate FGFR1/2/3/4 paralog panel summary\n")
    if panel_manifest:
        gene_counts = Counter(r.get("gene", "") for r in panel_manifest)
        group_counts = Counter(r.get("taxon_group", "") for r in panel_manifest)
        md.append(f"- Panel reference sequences: **{len(panel_manifest)}** "
                  f"(FGFR1/2/3/4 across {len([g for g in group_counts if g])} vertebrate groups).\n")
        md.append(f"- Per gene: {dict(gene_counts)}\n")
        md.append(f"- Per taxon group: {dict(group_counts)}\n")
    else:
        md.append("- Panel manifest not provided.\n")
    if paralog_meta:
        md.append(f"- Candidate proteins screened: **{paralog_meta.get('n_query_proteins_screened','?')}**\n")
        md.append(f"- Protein-level paralog status: {paralog_meta.get('paralog_status_counts', {})}\n")
        md.append(f"- Species-level paralog status: {paralog_meta.get('species_status_counts', {})}\n")
    md.append("The human-only FGFR1/2/3/4 panel is retained as a legacy control; the multi-vertebrate panel is "
              "the preferred paralog/orthology evidence layer. The screen does not change IIIb/IIIc labels.\n")

    md.append("\n## Orthology evidence counts\n")
    if orthology_meta:
        md.append(f"- Orthology records: **{orthology_meta.get('n_records','?')}**\n")
        md.append("- Per-record orthology status:\n")
        md.append(counts_block({k: int(v) for k, v in orthology_meta.get('orthology_status_counts', {}).items()}))
        md.append("- Per-species orthology status:\n")
        md.append(counts_block({k: int(v) for k, v in orthology_meta.get('species_orthology_status_counts', {}).items()}))
    else:
        md.append("- Orthology evidence metadata not provided.\n")

    md.append("\n## Figure export formats\n")
    md.append("- All pre-InterPro figures are exported in **SVG, PDF and PNG** (not PNG only); see `figure_manifest.tsv`.\n")

    md.append("\n## Exact command to run InterProScan later\n")
    md.append("Run from `" + str(args.interpro_dir) + "`:\n\n```bash\n" + interpro_cmd + "\n```\n")
    md.append("Then map results back with `fgfr2_interpro_id_mapping.tsv` (`unique_id` → species/isoform/role/"
              "transcript/protein/header/sequence_hash) and `fgfr2_interpro_unique_mapping.tsv` for duplicate expansion.\n")

    md.append("\n## Limitations\n")
    md.append("- **InterProScan / domain annotation has not been executed.** No real InterPro domains exist yet; "
              "figure domain tracks are labelled \"InterProScan pending\".\n")
    md.append("- The normalized III-slot axis normalizes cassette shape/length, not absolute biological position; "
              "its confidence is intentionally capped below `high`.\n")
    md.append("- CDS phase is unavailable for some Ensembl-sourced transcripts → `unknown_codon_phase` for those rows.\n")
    (out / "QC_migration_report_tasks_1_to_12_pre_interpro.md").write_text("".join(md), encoding="utf-8")

    # ---- methods update ----
    methods = f"""# Methods update (pre-InterProScan)

_Generated: {now}_

## IIIb/IIIc direction assignment (Step 4)
FGFR2 IIIb/IIIc identity is defined by the mutually exclusive exon/event architecture.
The provisional "first alternative exon = IIIb" order rule was replaced by a
**sequence-calibrated** assignment: each candidate cassette's translated CDS-exon
amino-acid sequence is locally aligned (Smith–Waterman) against curated human IIIb
and IIIc cassette references, and the direction is assigned by identity × query
coverage with an explicit margin. Order-rule labels are preserved as legacy columns.

## Protein evidence as QC only (Step 5b)
Protein-marker comparison is recorded as QC (`validation_status`, bounded
identity/coverage) and never triggers automatic isoform reassignment. Reference
self-control, human RefSeq protein control, human candidate control and a
close-primate control gate the run.

## CDS-boundary precision (Task 7)
CDS boundaries are evaluated in transcript/translation order using GFF3 phase
(negative strand respected). The left (5') boundary is codon-aligned iff phase = 0;
the right (3') boundary is codon-aligned iff (phase + CDS_length) mod 3 = 0.
Each boundary is reported as `exact` / `codon_split`, combined into
`cds_boundary_precision_refined` ∈ {{exact, codon_split_one_side,
codon_split_both_sides, unknown_codon_phase}}. The legacy single-boundary value is
preserved as `legacy_cds_boundary_precision`.

## InterProScan input preparation (Step 7 / Task 9)
Selected FGFR2 proteins are validated (no empty/invalid/stop-containing sequences),
collapsed to unique sequences with stable space-free IDs (e.g. `FGFR2_U0001`), and
exported as `fgfr2_interpro_clean_unique.fasta`. Mapping tables allow exact
reconstruction from any InterProScan result ID back to species/isoform/role/
transcript/protein/original header/sequence hash. **InterProScan is not run here.**
"""
    (out / "methods_update_pre_interpro.md").write_text(methods, encoding="utf-8")

    # ---- results summary ----
    res = f"""# Results summary (pre-InterProScan)

_Generated: {now}_

- Species analysed: **{len(master)}**
- main_analysis_eligible = true: **{eligible_counts.get('true', 0)}**, false: **{eligible_counts.get('false', 0)}**
- Display classes: {dict(display_counts)}
- Direction calibration: systematic order-rule inversions corrected at Step 4; human controls pass.
- CDS-boundary precision (refined): {refined_after}
- InterPro-ready unique sequences: **{interpro_summary.get('unique_sequences','?')}** (from {interpro_summary.get('total_selected_proteins','?')} selected; {interpro_summary.get('duplicates_collapsed','?')} duplicates collapsed)
- Review/supplementary species: **{len(review_species)}** ({', '.join(r.get('species','') for r in review_species) or 'none'})

InterProScan / domain annotation: **pending** (not executed).
"""
    (out / "results_summary_pre_interpro.md").write_text(res, encoding="utf-8")

    # ---- figure captions ----
    caps = f"""# Figure captions (pre-InterProScan)

_Generated: {now}_

**Figure 1 — Framework counts.** Per-category species counts (final display class,
taxon group, CDS-boundary precision, InterPro input status) from the pre-InterPro
species QC master. Green = high-confidence/clean, red = review/supplementary.

**Figure 2 — Exon-to-protein map.** Native protein-axis position of the IIIb (blue)
and IIIc (orange) cassettes per species. Hatching marks split/unknown codon
boundaries. `*` marks review/supplementary species. Protein-domain overlays are
**InterProScan pending** and are not shown.

**Figure 3 — Species evidence matrix.** Per-species support across evidence
dimensions (ortholog, paralog screen, both isoforms, direction, protein QC, native
coordinate, normalized slot, III similarity, CDS boundary, InterPro input).
Green = supported, red = review.

**Figure 4 — Native vs normalized III-slot QC.** Native pair-center distance vs
normalized III-slot pair-center distance, coloured by taxon group; review species
drawn as red-edged squares.

**Supplementary Figure — Review cases.** Table of review/supplementary species with
review reasons and recommended use.

All figures are exported as SVG, PDF and PNG and use a colour-blind-safe palette with
consistent IIIb/IIIc colours. No real InterPro domains are shown.
"""
    (out / "figure_captions_pre_interpro.md").write_text(caps, encoding="utf-8")

    # ---- output file manifest (walk BASE) ----
    manifest_rows: List[Dict[str, object]] = []
    base = args.base
    for p in sorted(base.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(base)
        try:
            size = p.stat().st_size
        except OSError:
            size = 0
        sha = ""
        if size and size < 50_000_000 and p.suffix in (".tsv", ".json", ".md", ".fasta", ".faa", ".txt"):
            h = hashlib.sha256()
            h.update(p.read_bytes())
            sha = h.hexdigest()[:16]
        manifest_rows.append({"path": str(rel), "size_bytes": size, "ext": p.suffix.lstrip("."), "sha256_16": sha})
    with open(out / "output_file_manifest_pre_interpro.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=["path", "size_bytes", "ext", "sha256_16"])
        w.writeheader()
        w.writerows(manifest_rows)

    # ---- run manifest ----
    run_manifest = {
        "script_name": "make_pre_interpro_report_package.py",
        "script_version": SCRIPT_VERSION,
        "generated_utc": now,
        "scope": "pipeline_ends_at_interproscan_preparation",
        "interproscan_executed": False,
        "n_species": len(master),
        "final_display_class_counts": dict(display_counts),
        "main_analysis_eligible_counts": dict(eligible_counts),
        "interpro_status_counts": dict(interpro_counts),
        "cds_boundary_precision": cds,
        "resolver_status_refined_counts": resolver_refined,
        "interpro_prepare_summary": interpro_summary,
        "paralog_panel": {
            "n_panel_sequences": len(panel_manifest),
            "screen_status_counts": paralog_meta.get("paralog_status_counts", {}),
            "species_status_counts": paralog_meta.get("species_status_counts", {}),
        },
        "orthology_evidence": {
            "orthology_status_counts": orthology_meta.get("orthology_status_counts", {}),
            "species_orthology_status_counts": orthology_meta.get("species_orthology_status_counts", {}),
        },
        "figure_export_formats": ["svg", "pdf", "png"],
        "review_species": [r.get("species", "") for r in review_species],
        "interproscan_command": interpro_cmd,
        "n_output_files": len(manifest_rows),
        "figures_dir": str(args.figures_dir),
        "tables_dir": str(args.tables_dir),
    }
    (out / "run_manifest_pre_interpro.json").write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")

    print(f"[OK] pre-InterPro report package written to {out}")
    print(f"     species={len(master)} display={dict(display_counts)} output_files={len(manifest_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
