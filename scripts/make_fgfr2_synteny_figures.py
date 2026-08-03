#!/usr/bin/env python3
"""
Render FGFR2 synteny figures.

Paper-level FGFR2 local gene-neighborhood / synteny figures (SVG/PDF/PNG):
  Figure 9A — FGFR2 local gene-neighborhood map (5 neighbors, representative species)
  Figure 9B — FGFR2 5/10-neighbor conservation matrix
  Figure 9C — Rescue-case locus panels (Gorilla / Canis / Pongo + high-risk)
  Supplement — FGFR2 local gene-neighborhood (10 neighbors, all species)

Synteny validates the FGFR2 locus / orthology context only; it never assigns or relabels IIIb/IIIc.
Neighbor labels use the normalized identity (symbol > curated > RBH > one-way BLAST > raw). BLAST/RBH
inferred names are shown as probable (italic + "?"). MCScanX block-level synteny (and Figure 9D) is
intentionally omitted from this build; the optional-figure slot is reported as unavailable.

Preferred backend: pyGenomeViz (recorded if importable); the comparative multi-row arrow maps use a
colour-blind-safe custom matplotlib renderer for full control of orthology colouring and review/rescue
outlines. DNA Features Viewer availability is also recorded.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402
from exondomaincompare.shared_gene_analysis.strand import strand_sign  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl_syn_"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrow, Patch, Rectangle  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

MAIN_SPECIES = ["homo_sapiens", "pan_troglodytes", "gorilla_gorilla_gorilla", "pongo_abelii",
                "macaca_mulatta", "mus_musculus", "canis_lupus_familiaris", "bos_taurus",
                "gallus_gallus", "xenopus_tropicalis", "danio_rerio"]
FGFR2_COLOR = "#111111"
UNMAPPED_COLOR = "#B0BEC5"
BROAD_COLOR = "#C5B0D5"  # loose human-proteome homology name (LOC... resolved by best hit)
AMBIGUOUS_OUTLINE = "#E69F00"
RESCUE_OUTLINE = "#0072B2"
REVIEW_OUTLINE = "#D55E00"
# colour-blind-safe palette (Okabe-Ito + extensions) for orthology groups
PALETTE = ["#009E73", "#56B4E9", "#E69F00", "#CC79A7", "#0072B2", "#D55E00", "#F0E442",
           "#999933", "#882255", "#44AA99", "#332288", "#AA4499", "#88CCEE", "#117733",
           "#DDCC77", "#661100", "#6699CC", "#888888", "#000000", "#E41A1C"]
MATRIX_STATE_COLOR = {
    "present_same_side_and_order": "#1B7837",
    "present_same_side_reordered": "#7FBF7B",
    "present_opposite_side": "#2166AC",
    "present_elsewhere_in_5neighbor_window": "#5AAE61",
    "present_only_in_10neighbor_supplement": "#A6DBA0",
    "probable_by_blast_rbh": "#9970AB",
    "ambiguous_identity": "#E69F00",
    "missing_or_unmapped": "#E0E0E0",
    "scaffold_unavailable": "#9E9E9E",
    "source_mismatch_review": "#D55E00",
}
SYNTENY_CLASS_COLOR = {
    "synteny_strong": "#1B7837",
    "synteny_supported_with_minor_rearrangement": "#7FBF7B",
    "synteny_partial_blast_supported": "#9970AB",
    "synteny_partial_scaffold_limit": "#80CDC1",
    "synteny_sequence_only_support": "#E69F00",
    "synteny_conflict_review": "#D55E00",
    "synteny_unavailable": "#BDBDBD",
}


def savefig(fig, fig_dir: Path, stem: str) -> None:
    for ext in ("svg", "pdf", "png"):
        fig.savefig(fig_dir / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def backend_report(dirs) -> List[Dict[str, object]]:
    rows = []
    for mod, fig in (("pygenomeviz", "Figure_9A/9C"), ("dna_features_viewer", "Figure_9C"),
                     ("matplotlib", "Figure_9A/9B/9C/Supplement")):
        try:
            m = __import__(mod)
            ver = getattr(m, "__version__", "?")
            status = "available_used" if mod == "matplotlib" else "available_not_required"
            warn = "" if mod == "matplotlib" else "custom matplotlib renderer used for comparative maps"
        except Exception:
            ver, status, warn = "", "missing", "not installed; fallback renderer used"
        rows.append({"backend": mod, "status": status, "version": ver, "figure": fig,
                     "warning": warn})
    rows.append({"backend": "mcscanx", "status": "omitted_by_design", "version": "",
                 "figure": "Figure_9D", "warning": "MCScanX block-level synteny not part of this build"})
    M.write_tsv(dirs["synteny"] / "fgfr2_neighborhood_plotting_backend_report.tsv", rows,
                ["backend", "status", "version", "figure", "warning"])
    return rows


# ---------------------------------------------------------------------------
# data loading
# ---------------------------------------------------------------------------
def load(base, dirs):
    syn = dirs["synteny"]
    master = {(r["species"] or "").lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}
    n5 = M.read_tsv(syn / "fgfr2_local_gene_neighborhood_5neighbors.tsv")
    n10 = M.read_tsv(syn / "fgfr2_local_gene_neighborhood_10neighbors_supplement.tsv")
    ident = {((r["species"] or "").lower(), r["neighbor_gene_id"]): r for r in
             M.read_tsv(syn / "fgfr2_neighbor_identity_resolution.tsv")}
    valid = {(r["species"] or "").lower(): r for r in
             M.read_tsv(syn / "fgfr2_5neighbor_synteny_validation.tsv")}
    matrix = M.read_tsv(syn / "fgfr2_5neighbor_conservation_matrix.tsv")
    href = M.read_tsv(syn / "human_fgfr2_10neighbor_reference.tsv")
    truth = {}
    for r in M.read_tsv(dirs["maps"] / "fgfr2_post_rescue_final_truth_table.tsv"):
        truth.setdefault((r["species"] or "").lower(), []).append(r)
    return master, n5, n10, ident, valid, matrix, href, truth


def neighbor_layout(rows, ident, species, max_n) -> List[Dict[str, object]]:
    """Ordered display slots: upstream rank max..1, FGFR2 (0), downstream 1..max."""
    sp = species.lower()
    spr = [r for r in rows if (r["species"] or "").lower() == sp and r.get("neighbor_rank")]
    by = {(r["neighbor_side"], M.to_int(r["neighbor_rank"], 0)): r for r in spr}
    fr = next((r for r in rows if (r["species"] or "").lower() == sp), {})
    # centered coordinate system: upstream rank r at x=-r, FGFR2 at 0, downstream rank r at x=+r
    slots = []
    for rank in range(max_n, 0, -1):
        slots.append((-rank, "upstream", rank, by.get(("upstream", rank))))
    slots.append((0, "fgfr2", 0, fr))
    for rank in range(1, max_n + 1):
        slots.append((rank, "downstream", rank, by.get(("downstream", rank))))
    out = []
    for x, side, rank, r in slots:
        if side == "fgfr2":
            out.append({"x": x, "side": side, "rank": 0, "symbol": "FGFR2", "method": "anchor",
                        "status": "anchor", "strand": r.get("fgfr2_strand", "+"),
                        "ortho": "FGFR2", "raw": "FGFR2", "is_fgfr2": True})
            continue
        if not r:
            out.append({"x": x, "side": side, "rank": rank, "symbol": "", "method": "missing",
                        "status": "missing", "strand": "", "ortho": "", "raw": "", "is_fgfr2": False,
                        "empty": True})
            continue
        idr = ident.get((sp, r["neighbor_gene_id"]), {})
        out.append({"x": x, "side": side, "rank": rank,
                    "symbol": idr.get("normalized_neighbor_symbol") or r["neighbor_symbol_raw"],
                    "method": idr.get("identity_resolution_method", "raw_annotation_only"),
                    "status": idr.get("identity_resolution_status", "raw_id_only"),
                    "strand": r.get("neighbor_strand", "+"),
                    "ortho": idr.get("normalized_neighbor_orthology_group", ""),
                    "raw": r["neighbor_symbol_raw"], "is_fgfr2": False,
                    "broad_pid": idr.get("broad_homology_percent_identity", ""),
                    "broad_cov": idr.get("broad_homology_query_coverage", ""),
                    "distance": M.to_int(r.get("distance_to_fgfr2"), 0)})
    return out


def color_map(layouts: List[List[Dict[str, object]]]) -> Dict[str, str]:
    groups = []
    for lay in layouts:
        for s in lay:
            g = s.get("ortho")
            if g and g != "FGFR2" and g not in groups:
                groups.append(g)
    cmap = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(sorted(groups))}
    cmap["FGFR2"] = FGFR2_COLOR
    return cmap


def draw_locus_row(ax, y, layout, cmap, fgfr2_strand, label_top=False):
    fstrand = strand_sign(fgfr2_strand)
    for s in layout:
        x = s["x"]
        if s.get("empty"):
            ax.plot([x - 0.3, x + 0.3], [y, y], color="#CCCCCC", lw=0.8, ls=":")
            continue
        ortho = s.get("ortho") or ""
        method = s.get("method", "")
        is_broad = method in ("broad_proteome_best_hit", "broad_proteome_weak_best_hit")
        if s["is_fgfr2"]:
            face = cmap.get("FGFR2")
        elif ortho:
            face = cmap.get(ortho, UNMAPPED_COLOR)
        elif is_broad:
            face = BROAD_COLOR
        else:
            face = UNMAPPED_COLOR
        strand = strand_sign(s.get("strand", "+"))
        rel = strand * fstrand  # orient relative to FGFR2 strand
        dx = 0.62 * rel
        x0 = x - 0.31 * rel
        edge, lw = "black", 0.5
        st = s.get("status", "")
        if st == "ambiguous_neighbor_identity":
            edge, lw = AMBIGUOUS_OUTLINE, 1.8
        named = s["is_fgfr2"] or ortho or st == "symbol_supported_only" or is_broad
        ax.add_patch(FancyArrow(x0, y, dx, 0, width=0.34, head_width=0.34, head_length=0.22,
                                length_includes_head=True, facecolor=face, edgecolor=edge,
                                linewidth=lw, zorder=3, alpha=0.95 if named else 0.55))
        # label
        lab = "FGFR2" if s["is_fgfr2"] else (s["symbol"] or s["raw"])
        italic = method in ("reciprocal_best_hit", "high_confidence_one_way_blast") or is_broad
        if method in ("raw_annotation_only", "unresolved") and not s["is_fgfr2"]:
            lab = s["raw"]
        suffix = "?" if italic else ""
        ty = y + (0.42 if label_top else -0.42)
        lab_color = (FGFR2_COLOR if s["is_fgfr2"] else "#5E3C99" if is_broad
                     else "#222222" if (ortho or st == "symbol_supported_only") else "#888888")
        ax.text(x, ty, lab + suffix, ha="center", va="bottom" if label_top else "top",
                fontsize=6.0, rotation=35 if not s["is_fgfr2"] else 0,
                fontstyle="italic" if italic else "normal", color=lab_color,
                fontweight="bold" if s["is_fgfr2"] else "normal")
        # mark broad-homology coverage (percent identity) directly on the map
        if is_broad and s.get("broad_pid") not in (None, ""):
            cy = y + (0.74 if label_top else -0.78)
            ax.text(x, cy, f"{float(s['broad_pid']):.0f}%", ha="center",
                    va="bottom" if label_top else "top", fontsize=4.6, color="#5E3C99")


# ---------------------------------------------------------------------------
# Figure 9A — local gene-neighborhood map (5 neighbors)
# ---------------------------------------------------------------------------
def fig_9a(base, dirs, master, n5, ident, valid, truth):
    # An ordering preference, not a membership test: the named species lead, then every other
    # species this dataset analysed. Filtering on the thesis panel dropped any species outside
    # it from the figure entirely.
    named = [s for s in MAIN_SPECIES if s in master]
    sps = named + [s for s in master if s not in set(named)]
    layouts = {sp: neighbor_layout(n5, ident, sp, 5) for sp in sps}
    cmap = color_map(list(layouts.values()))
    fig, ax = plt.subplots(figsize=(13, max(5, len(sps) * 0.92)))
    table_rows = []
    for yi, sp in enumerate(sps):
        y = len(sps) - yi
        lay = layouts[sp]
        fstr = next((r.get("fgfr2_strand", "+") for r in n5
                     if (r["species"] or "").lower() == sp), "+")
        draw_locus_row(ax, y, lay, cmap, fstr)
        v = valid.get(sp, {})
        cls = v.get("synteny_validation_class", "")
        claim = ";".join(sorted({r.get("final_claim_status_after_rescue", "")
                                 for r in truth.get(sp, [])} - {""}))
        rescued = any((r.get("rescue_decision") or "").startswith("rescued")
                      for r in truth.get(sp, []))
        # row outline cue for rescue / review
        if rescued:
            ax.add_patch(Rectangle((-5.6, y - 0.5), 11.2, 1.0, fill=False,
                                   edgecolor=RESCUE_OUTLINE, lw=1.4, ls="--", zorder=1))
        elif "supplement" in claim or "excluded" in claim:
            ax.add_patch(Rectangle((-5.6, y - 0.5), 11.2, 1.0, fill=False,
                                   edgecolor=REVIEW_OUTLINE, lw=1.2, ls=":", zorder=1))
        disp = master.get(sp, {}).get("display_species_name", sp)
        ax.text(-6.0, y, disp, ha="right", va="center", fontsize=8, fontstyle="italic")
        ax.add_patch(Rectangle((5.75, y - 0.32), 0.3, 0.64,
                               facecolor=SYNTENY_CLASS_COLOR.get(cls, "#CCCCCC"),
                               edgecolor="none", clip_on=False, zorder=4))
        for s in lay:
            table_rows.append({"species": sp, "slot_x": s["x"], "side": s["side"],
                               "rank": s["rank"], "normalized_symbol": s.get("symbol", ""),
                               "raw_symbol": s.get("raw", ""), "orthology_group": s.get("ortho", ""),
                               "identity_method": s.get("method", ""),
                               "identity_status": s.get("status", ""),
                               "broad_homology_percent_identity": s.get("broad_pid", ""),
                               "broad_homology_query_coverage": s.get("broad_cov", ""),
                               "strand": s.get("strand", ""),
                               "synteny_validation_class": cls,
                               "final_claim_status_after_rescue": claim})
    ax.set_xlim(-7.5, 6.4)
    ax.set_ylim(0.3, len(sps) + 0.8)
    ax.axis("off")
    ax.set_title("Figure 9A — FGFR2 local gene-neighborhood (5 protein-coding neighbors per side)\n"
                 "Synteny validates the FGFR2 locus/orthology context; it does not assign IIIb/IIIc.",
                 fontsize=12, fontweight="bold")
    leg = [Patch(facecolor=FGFR2_COLOR, label="FGFR2 (anchor)"),
           Patch(facecolor="#009E73", label="orthologous neighbor (symbol/curated)"),
           Patch(facecolor="#9970AB", label="RBH/one-way BLAST (probable, italic + ?)"),
           Patch(facecolor=BROAD_COLOR, label="LOC named by human homolog (loose, %id shown)"),
           Patch(facecolor=UNMAPPED_COLOR, label="no human homolog / raw ID"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor="w",
                  markeredgecolor=AMBIGUOUS_OUTLINE, markeredgewidth=2, label="ambiguous identity"),
           Line2D([0], [0], color=RESCUE_OUTLINE, ls="--", label="rescued candidate locus"),
           Line2D([0], [0], color=REVIEW_OUTLINE, ls=":", label="supplement/review species")]
    ax.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, -0.15), ncol=4,
              fontsize=7, frameon=False)
    M.write_tsv(dirs["tables"] / "figure9A_fgfr2_local_gene_neighborhood_5neighbors.tsv",
                table_rows, list(table_rows[0].keys()) if table_rows else ["species"])
    savefig(fig, dirs["figures"], "Figure_9A_FGFR2_local_gene_neighborhood_5neighbors")


# ---------------------------------------------------------------------------
# Figure 9B — conservation matrix
# ---------------------------------------------------------------------------
def fig_9b(base, dirs, master, matrix, href, valid, truth):
    cols = [r["human_gene_symbol"] for r in sorted(
        href, key=lambda r: (r["human_neighbor_side"], M.to_int(r["human_neighbor_rank"], 0)))]
    rows = sorted(matrix, key=lambda r: M.to_int(master.get((r["species"] or "").lower(), {})
                                                 .get("phylo_order"), 999) or 999)
    n, m = len(rows), len(cols)
    fig, ax = plt.subplots(figsize=(max(10, m * 0.7 + 6), max(7, n * 0.34)))
    for yi, r in enumerate(rows):
        for xi, c in enumerate(cols):
            state = r.get(c, "missing_or_unmapped")
            ax.add_patch(Rectangle((xi, yi), 1, 1,
                                   facecolor=MATRIX_STATE_COLOR.get(state, "#E0E0E0"),
                                   edgecolor="white", lw=0.5))
        sp = (r["species"] or "").lower()
        v = valid.get(sp, {})
        # sidebars: taxon, claim, rescue, synteny class
        claim = ";".join(sorted({t.get("final_claim_status_after_rescue", "")
                                 for t in truth.get(sp, [])} - {""}))
        rescued = any((t.get("rescue_decision") or "").startswith("rescued")
                      for t in truth.get(sp, []))
        bars = [("#607D8B" if r.get("taxon_group") else "#ECEFF1"),
                ("#1B7837" if claim.startswith("primary") else "#E69F00"),
                (RESCUE_OUTLINE if rescued else "#ECEFF1"),
                SYNTENY_CLASS_COLOR.get(v.get("synteny_validation_class", ""), "#CCCCCC")]
        for bi, col in enumerate(bars):
            ax.add_patch(Rectangle((m + 0.15 + bi * 0.45, yi), 0.4, 1, facecolor=col,
                                   edgecolor="white", lw=0.4, clip_on=False))
    ax.set_xlim(0, m + 0.15 + 4 * 0.45 + 0.2)
    ax.set_ylim(0, n)
    ax.invert_yaxis()
    ax.set_xticks([i + 0.5 for i in range(m)])
    ax.set_xticklabels(cols, rotation=55, ha="right", fontsize=7)
    ax.set_yticks([i + 0.5 for i in range(n)])
    ax.set_yticklabels([master.get((r["species"] or "").lower(), {})
                        .get("display_species_name", r["species"]) for r in rows], fontsize=6.5)
    sb_labels = ["taxon", "claim", "rescue", "synteny"]
    ax.set_xticks(list(ax.get_xticks()) + [m + 0.35 + bi * 0.45 for bi in range(4)])
    for bi, lab in enumerate(sb_labels):
        ax.text(m + 0.35 + bi * 0.45, -0.4, lab, rotation=55, ha="left", va="bottom", fontsize=6)
    ax.set_title("Figure 9B — FGFR2 5/10-neighbor conservation matrix (human reference neighbors)\n"
                 "Locus/orthology context only; IIIb/IIIc labels are not derived from synteny.",
                 fontsize=12, fontweight="bold", pad=24)
    leg = [Patch(facecolor=c, label=k.replace("_", " ")) for k, c in MATRIX_STATE_COLOR.items()]
    ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, frameon=False)
    M.write_tsv(dirs["tables"] / "figure9B_fgfr2_5neighbor_conservation_matrix.tsv", rows,
                list(rows[0].keys()) if rows else ["species"])
    savefig(fig, dirs["figures"], "Figure_9B_FGFR2_5neighbor_conservation_matrix")


# ---------------------------------------------------------------------------
# Figure 9C — rescue-case locus panels
# ---------------------------------------------------------------------------
def fig_9c(base, dirs, master, n5, ident, valid, truth):
    cases = ["gorilla_gorilla_gorilla", "canis_lupus_familiaris", "pongo_abelii"]
    # add any other high-risk unresolved species (supplement/excluded claim) with a locus
    for sp, ts in truth.items():
        claim = {t.get("final_claim_status_after_rescue", "") for t in ts}
        if sp not in cases and any(c in ("supplement_review", "excluded_from_primary_claim")
                                   for c in claim):
            cases.append(sp)
    cases = [s for s in cases if s in master]
    layouts = {sp: neighbor_layout(n5, ident, sp, 5) for sp in cases}
    cmap = color_map(list(layouts.values()))
    n = len(cases)
    if n == 0:
        # Small custom runs may have no rescue/supplement/review species at all. There is
        # nothing to draw here; write an empty table and skip the panel instead of crashing
        # matplotlib with a 0-row subplot grid (the neighborhood track is JSON/9A-driven).
        M.write_tsv(dirs["tables"] / "figure9C_rescue_case_locus_panels.tsv", [], ["species"])
        return
    fig, axes = plt.subplots(n, 1, figsize=(12, max(4, n * 1.7)))
    axes = axes if hasattr(axes, "__len__") else [axes]
    table_rows = []
    for ax, sp in zip(axes, cases):
        lay = layouts[sp]
        fstr = next((r.get("fgfr2_strand", "+") for r in n5
                     if (r["species"] or "").lower() == sp), "+")
        draw_locus_row(ax, 1.0, lay, cmap, fstr, label_top=True)
        v = valid.get(sp, {})
        ts = truth.get(sp, [])
        claim = ";".join(sorted({t.get("final_claim_status_after_rescue", "") for t in ts} - {""}))
        dec = ";".join(sorted({t.get("rescue_decision", "") for t in ts} - {""}))
        unres = ";".join(sorted({t.get("unresolved_reason_if_any", "") for t in ts} - {""})) or "-"
        disp = master.get(sp, {}).get("display_species_name", sp)
        ax.set_xlim(-6.2, 6.2)
        ax.set_ylim(-0.4, 2.1)
        ax.axis("off")
        ax.text(-6.0, 1.7, disp, fontsize=9, fontweight="bold", fontstyle="italic")
        info = (f"synteny: {v.get('synteny_validation_class','')}  |  "
                f"neighbor support: {v.get('total_neighbor_support_score','')}  |  "
                f"locus: {v.get('rescued_candidate_locus_support','')}\n"
                f"rescue: {dec}  |  claim: {claim}  |  unresolved: {unres}")
        ax.text(-6.0, 0.05, info, fontsize=6.8, va="top", family="monospace", color="#333333")
        table_rows.append({"species": sp, "synteny_validation_class": v.get("synteny_validation_class", ""),
                           "total_neighbor_support_score": v.get("total_neighbor_support_score", ""),
                           "rescued_candidate_locus_support": v.get("rescued_candidate_locus_support", ""),
                           "rescue_decision": dec, "final_claim_status_after_rescue": claim,
                           "unresolved_reason_if_any": unres})
    fig.suptitle("Figure 9C — FGFR2 rescue-case locus panels (Gorilla / Canis / Pongo + high-risk)\n"
                 "Local neighborhood context for rescued/partial cases; synteny does not relabel IIIb/IIIc.",
                 fontsize=12, fontweight="bold", y=1.0)
    M.write_tsv(dirs["tables"] / "figure9C_rescue_case_locus_panels.tsv", table_rows,
                list(table_rows[0].keys()) if table_rows else ["species"])
    savefig(fig, dirs["figures"], "Figure_9C_FGFR2_rescue_case_locus_panels")


# ---------------------------------------------------------------------------
# Supplement — 10-neighbor context, all species
# ---------------------------------------------------------------------------
def fig_supp(base, dirs, master, n10, ident, valid, truth):
    sps = sorted(master.keys(), key=lambda s: M.to_int(master[s].get("phylo_order"), 999) or 999)
    layouts = {sp: neighbor_layout(n10, ident, sp, 10) for sp in sps}
    cmap = color_map(list(layouts.values()))
    fig, ax = plt.subplots(figsize=(20, max(8, len(sps) * 0.62)))
    table_rows = []
    for yi, sp in enumerate(sps):
        y = len(sps) - yi
        fstr = next((r.get("fgfr2_strand", "+") for r in n10
                     if (r["species"] or "").lower() == sp), "+")
        draw_locus_row(ax, y, layouts[sp], cmap, fstr)
        disp = master.get(sp, {}).get("display_species_name", sp)
        ax.text(-11.0, y, disp, ha="right", va="center", fontsize=7, fontstyle="italic")
        v = valid.get(sp, {})
        ax.add_patch(Rectangle((10.7, y - 0.32), 0.4, 0.64,
                               facecolor=SYNTENY_CLASS_COLOR.get(
                                   v.get("synteny_validation_class", ""), "#CCCCCC"),
                               edgecolor="none", clip_on=False, zorder=4))
        for s in layouts[sp]:
            table_rows.append({"species": sp, "slot_x": s["x"], "side": s["side"], "rank": s["rank"],
                               "normalized_symbol": s.get("symbol", ""), "raw_symbol": s.get("raw", ""),
                               "orthology_group": s.get("ortho", ""),
                               "identity_status": s.get("status", "")})
    ax.set_xlim(-13.5, 11.4)
    ax.set_ylim(0.3, len(sps) + 0.8)
    ax.axis("off")
    ax.set_title("Supplement — FGFR2 local gene-neighborhood (10 protein-coding neighbors per side, "
                 "all species)\nShows whether human reference neighbors fall just outside the "
                 "5-neighbor main window; synteny does not relabel IIIb/IIIc.",
                 fontsize=12, fontweight="bold")
    M.write_tsv(dirs["tables"] / "supplement_figure_fgfr2_neighborhood_10neighbors_all_species.tsv",
                table_rows, list(table_rows[0].keys()) if table_rows else ["species"])
    savefig(fig, dirs["figures"],
            "Supplement_Figure_FGFR2_local_gene_neighborhood_10neighbors_all_species")


def main() -> int:
    ap = argparse.ArgumentParser(description="FGFR2 synteny / gene-neighborhood figures (Part H).")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    ok, msgs = M.synteny_gate(base)
    if not ok:
        print("[FAIL] synteny validation gate failed; no synteny figures plotted:", file=sys.stderr)
        for m in msgs:
            print("   - " + m, file=sys.stderr)
        return 2
    backend_report(dirs)
    master, n5, n10, ident, valid, matrix, href, truth = load(base, dirs)
    fig_9a(base, dirs, master, n5, ident, valid, truth)
    fig_9b(base, dirs, master, matrix, href, valid, truth)
    fig_9c(base, dirs, master, n5, ident, valid, truth)
    fig_supp(base, dirs, master, n10, ident, valid, truth)
    print("[OK] synteny figures written (Figure 9A, 9B, 9C + 10-neighbor supplement; svg/pdf/png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
