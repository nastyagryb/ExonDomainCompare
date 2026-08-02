#!/usr/bin/env python3
"""
make_publication_figures_pre_interpro.py  (Sprint Parts 3, 5, 6)

Publication-level pre-InterPro figure engine.

Strict rules:
  * Reads final tables only (species_qc_master.tsv is canonical QC/display).
  * Uses species_phylogenetic_order.tsv / phylo_order for species ordering.
  * Does NOT recompute biological QC, main_analysis_eligible, same_slot_or_adjacent,
    IIIb/IIIc labels, III-region similarity, or review classification.
  * Fails clearly if required fields are missing.
  * No fake InterPro domain calls. Any architecture guide is explicitly labelled
    "canonical FGFR2 architecture guide / InterProScan pending".

This module is normally driven by scripts/make_all_figures.py (--base).
"""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl_pub_preinterpro_"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D


# ---------------------------------------------------------------------------
# Stable, color-blind-safe style (Part 3)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "svg.fonttype": "none",
    "axes.edgecolor": "#404040",
    "axes.linewidth": 0.8,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

C_IIIB = "#0072B2"      # blue/teal  (stable across all figures)
C_IIIC = "#E69F00"      # orange/vermillion (stable across all figures)
C_COMMON = "#D9D9D9"    # common exons / protein body
C_AXIS = "#404040"      # protein backbone / axis
C_REVIEW = "#E69F00"    # amber marker for review (not aggressive red)
C_PENDING = "#BBBBBB"   # pending / missing

# Discrete evidence states for the QC heatmap (Paul Tol color-blind-safe).
STATE_COLORS = {
    "strong": "#117733",          # dark green
    "supported": "#DDCC77",       # sand
    "transcript_only": "#88CCEE", # light blue
    "ambiguous_review": "#CC6677",# muted rose
    "pending_missing": "#BBBBBB", # grey
}
STATE_LABELS = {
    "strong": "strong / pass",
    "supported": "supported / warning",
    "transcript_only": "transcript-level only",
    "ambiguous_review": "ambiguous / review",
    "pending_missing": "pending / missing",
}

# Refined uncertainty visual states (uncertainty-refinement sprint, Part E).
# Minor technical boundary flags (split codon, phase unavailable) get calm,
# non-alarming colors; only true missing data / conflicts / hard fails are loud.
REFINED_STATE_COLORS = {
    "robust_pass": "#2A9D8F",          # blue-green
    "split_codon": "#9ECAE1",          # light blue
    "phase_unavailable": "#C6D8E4",     # pale grey-blue
    "transcript_protein_only": "#B0A8C0",  # muted purple/grey
    "review": "#E6A532",               # amber
    "hard_fail": "#CC3311",            # red (true hard fail only)
    "interpro_pending": "#BBBBBB",      # neutral grey
}
REFINED_STATE_LABELS = {
    "robust_pass": "robust / pass (coordinates resolved, codon-exact)",
    "split_codon": "split codon = known boundary, not an error",
    "phase_unavailable": "phase unavailable = coordinate resolved, phase not inferable",
    "transcript_protein_only": "transcript/protein-level only",
    "review": "review = not used for primary claims, interpreted separately",
    "hard_fail": "hard fail (true coordinate sanity failure)",
    "interpro_pending": "InterPro pending (neutral)",
}


DISPLAY_RANK_FIG4 = {
    "hard_fail_excluded": 6, "protein_overlay_only": 5, "review_protein": 4,
    "review_annotation": 3, "resolved_phase_not_available": 2,
    "resolved_with_split_codon": 1, "robust": 0,
}


def refined_fig4_state(display_class: str, protein_evidence: str) -> str:
    """Map a refined display_uncertainty_class (+ protein evidence) onto a calm,
    color-blind-safe visual state for the evidence matrix."""
    d = display_class or ""
    if d == "hard_fail_excluded":
        return "hard_fail"
    if d in ("review_protein", "review_annotation", "protein_overlay_only"):
        return "review"
    if protein_evidence == "protein_transcript_level_only":
        return "transcript_protein_only"
    if d == "resolved_phase_not_available":
        return "phase_unavailable"
    if d == "resolved_with_split_codon":
        return "split_codon"
    return "robust_pass"

# Subtle, desaturated taxon-group side-bar colors (must not compete with IIIb/IIIc).
TAXON_BAR = {
    "Primates": "#7B68A6",
    "Other mammals": "#5AA469",
    "Birds": "#C97B3C",
    "Reptiles": "#4E89AE",
    "Amphibians": "#5FA8A8",
    "Teleost fish": "#8C8C8C",
}


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


def _to_int(v: object, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(float(str(v).strip()))
    except Exception:
        return default


def is_review_class(fdc: str) -> bool:
    return str(fdc).startswith("supplementary")


def is_minor_review(fdc: str) -> bool:
    return "minor_review" in str(fdc)


# ---------------------------------------------------------------------------
# Data loading (final tables only)
# ---------------------------------------------------------------------------
def _norm_tx(tx: str) -> str:
    t = (tx or "").strip()
    for pref in ("rna-", "transcript:", "transcript-"):
        if t.lower().startswith(pref):
            t = t[len(pref):]
    return t


class FinalData:
    def __init__(self, master: Path, coord_audit: Path, interpro_summary: Path,
                 cds_features: Optional[Path] = None, cassette_map: Optional[Path] = None,
                 refined: Optional[Path] = None):
        for p in (master, coord_audit, interpro_summary):
            if not p or not Path(p).exists() or Path(p).stat().st_size == 0:
                raise RuntimeError(f"Required final table missing or empty: {p}")
        self.master = read_tsv(master)
        # canonical species order by phylo_order
        if not self.master or "phylo_order" not in self.master[0]:
            raise RuntimeError("species_qc_master.tsv lacks phylo_order; run Part 2 first.")
        self.master.sort(key=lambda r: (_to_int(r.get("phylo_order"), 10 ** 6), r.get("species", "")))
        self.species_order = [r["species"] for r in self.master]
        self.master_by_sp = {r["species"]: r for r in self.master}

        self.coords: Dict[Tuple[str, str], Dict[str, str]] = {}
        for r in read_tsv(coord_audit):
            sp = (r.get("species_canonical") or "").strip().lower()
            iso = (r.get("inferred_isoform") or "").strip()
            if sp and iso:
                self.coords[(sp, iso)] = r

        self.interpro = {r["metric"]: r["value"] for r in read_tsv(interpro_summary)}

        # real per-transcript CDS blocks (protein-coordinate projection), if available
        self.cds_by_tx: Dict[str, List[Dict[str, str]]] = {}
        if cds_features and Path(cds_features).exists():
            tmp: Dict[str, List[Dict[str, str]]] = {}
            for c in read_tsv(cds_features):
                tmp.setdefault(_norm_tx(c.get("transcript_id_source")), []).append(c)
            for tx, lst in tmp.items():
                lst.sort(key=lambda c: _to_int(c.get("cds_rank"), 0) or 0)
                self.cds_by_tx[tx] = lst

        # canonical cassette->CDS-block mapping (coordinate-overlap, NOT cds_id);
        # fixes the bug where non-unique NCBI cds_id collapsed cassettes onto block 1.
        self.cassette_map: Dict[Tuple[str, str], Dict[str, str]] = {}
        if cassette_map and Path(cassette_map).exists():
            for m in read_tsv(cassette_map):
                self.cassette_map[((m.get("species") or "").lower(),
                                   (m.get("isoform") or ""))] = m

        # refined uncertainty / display classes (Part A); drives plot visibility
        self.refined_map: Dict[Tuple[str, str], Dict[str, str]] = {}
        if refined and Path(refined).exists():
            for m in read_tsv(refined):
                self.refined_map[((m.get("species") or "").lower(),
                                  (m.get("isoform") or ""))] = m

    def coord(self, species: str, isoform: str) -> Optional[Dict[str, str]]:
        return self.coords.get((species.lower(), isoform))

    def cds_blocks(self, species: str, isoform: str) -> List[Dict[str, str]]:
        """Return real CDS blocks (with protein_start_aa/protein_end_aa) for the
        resolved transcript of (species, isoform); empty if unavailable."""
        c = self.coord(species, isoform)
        if not c:
            return []
        tx = _norm_tx(c.get("transcript_id_source"))
        return self.cds_by_tx.get(tx, [])

    def cassette(self, species: str, isoform: str) -> Dict[str, str]:
        """Canonical coordinate-overlap cassette mapping for (species, isoform)."""
        return self.cassette_map.get((species.lower(), isoform), {})

    def cassette_rank(self, species: str, isoform: str) -> Optional[int]:
        """Unique CDS rank of the resolved cassette block (cds_rank is unique within a
        transcript, unlike the repeated NCBI cds_id)."""
        return _to_int(self.cassette(species, isoform).get("matched_cds_rank"))

    def refined(self, species: str, isoform: str) -> Dict[str, str]:
        """Refined uncertainty/display classes for (species, isoform)."""
        return self.refined_map.get((species.lower(), isoform), {})



def group_spans(master_rows: List[Dict[str, str]]) -> List[Tuple[str, int, int]]:
    """Return [(taxon_group, first_idx, last_idx)] for consecutive species blocks."""
    spans: List[Tuple[str, int, int]] = []
    cur = None
    start = 0
    for i, r in enumerate(master_rows):
        g = r.get("taxon_group", "") or "unknown"
        if cur is None:
            cur, start = g, i
        elif g != cur:
            spans.append((cur, start, i - 1))
            cur, start = g, i
    if cur is not None:
        spans.append((cur, start, len(master_rows) - 1))
    return spans


def draw_taxon_bands(ax, master_rows: List[Dict[str, str]], y_of, x0: float, x1: float,
                     bar_x: Optional[float] = None, label: bool = True) -> None:
    """Subtle horizontal background bands + thin left color bar per taxon group.

    y_of(idx) maps a species index (0=top) to its row center y. Bands span the
    full row height; coloring is desaturated so IIIb/IIIc remain dominant.
    """
    n = len(master_rows)
    row_h = abs(y_of(0) - y_of(1)) if n > 1 else 1.0
    for gi, (g, a, b) in enumerate(group_spans(master_rows)):
        color = TAXON_BAR.get(g, "#999999")
        y_top = y_of(a) + row_h / 2.0
        y_bot = y_of(b) - row_h / 2.0
        # very light background band (alternating subtle)
        ax.add_patch(Rectangle((x0, y_bot), x1 - x0, y_top - y_bot,
                               facecolor=color, edgecolor="none", alpha=0.06, zorder=0))
        # thin left color bar
        bx = bar_x if bar_x is not None else x0
        ax.add_patch(Rectangle((bx, y_bot), (x1 - x0) * 0.006 + 0.0, y_top - y_bot,
                               facecolor=color, edgecolor="none", alpha=0.85, zorder=2,
                               clip_on=False))
        # separator line between groups
        if gi > 0:
            ax.axhline(y_top, color="#CCCCCC", lw=0.6, zorder=1)
        # group label at far left
        if label:
            ax.text(x0 - (x1 - x0) * 0.015, (y_top + y_bot) / 2.0, g,
                    rotation=90, va="center", ha="right", fontsize=7.5,
                    color=color, fontweight="bold", clip_on=False)


def _save(fig, figdir: Path, stem: str) -> Dict[str, str]:
    figdir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 400})):
        p = figdir / f"{stem}.{ext}"
        fig.savefig(p, bbox_inches="tight", **kw)
        paths[ext] = str(p)
    plt.close(fig)
    return paths


def _skip_optional_figure(*, figure_id: str, func_name: str, figdir: Path,
                          reason: str, main_message: str, log_message: str,
                          n_species: int = 0, source_tables: str = "") -> Dict[str, object]:
    """Gracefully skip an OPTIONAL publication/display figure whose drawable
    subset is empty (e.g. a small custom run with zero review-case groups).

    This never calls plt.subplots. It writes a small placeholder note explaining
    the skip and returns a manifest entry marked ``status=skipped_empty`` so the
    pre-InterPro pipeline can continue. Only optional display figures may be
    skipped here — required validation gates are unaffected.
    """
    print(f"[SKIP optional figure] {func_name}: {log_message} (reason={reason})")
    figdir.mkdir(parents=True, exist_ok=True)
    placeholder = figdir / f"{figure_id}_SKIPPED_EMPTY.txt"
    placeholder.write_text(
        f"Figure '{figure_id}' was skipped for this run.\n"
        f"status: skipped_empty\n"
        f"reason: {reason}\n"
        f"rows: 0\n"
        f"function: {func_name}\n"
        f"generated_at: {_now()}\n\n"
        "This is an OPTIONAL publication/display figure. Its drawable subset was "
        "empty for this run (for example, a small custom run with no review cases), "
        "so it was skipped without aborting the pipeline. Required pre-InterPro "
        "validation is unaffected.\n",
        encoding="utf-8")
    return {
        "figure_id": figure_id,
        "paths": {"svg": "", "pdf": "", "png": "", "skip_note": str(placeholder)},
        "source_tables": source_tables,
        "n_species": n_species,
        "n_rows": 0,
        "status": "skipped_empty",
        "skip_reason": reason,
        "main_message": main_message,
    }


# ---------------------------------------------------------------------------
# Discrete evidence-state classification for the QC heatmap (read-only mapping;
# no biological QC is recomputed, only existing master strings are projected
# onto color-blind-safe display states).
# ---------------------------------------------------------------------------
# Ordered (substring, state); first match wins, per master column.
EVIDENCE_LAYERS: List[Tuple[str, str, List[Tuple[str, str]]]] = [
    ("fgfr2_ortholog_status", "FGFR2 orthology", [
        ("high_confidence", "strong"), ("supported", "supported"),
        ("transcript_level_only", "transcript_only"), ("review", "ambiguous_review")]),
    ("paralog_screen_status", "paralog panel", [
        ("high_confidence", "strong"), ("supported", "supported"),
        ("ambiguous", "ambiguous_review"), ("non_fgfr2", "ambiguous_review"),
        ("unavailable", "pending_missing")]),
    ("both_isoforms_detected", "both isoforms", [
        ("true", "strong"), ("false", "ambiguous_review")]),
    ("direction_validation_status", "calibrated direction", [
        ("calibrated", "strong"), ("unresolved_no_sequence", "transcript_only"),
        ("ambiguous", "ambiguous_review"), ("review", "ambiguous_review")]),
    ("protein_validation_summary", "protein QC", [
        ("validated", "strong"), ("conflict", "ambiguous_review"),
        ("ambiguous", "ambiguous_review"), ("review", "ambiguous_review")]),
    ("resolver_status_summary", "resolver status", [
        ("gold_exact", "strong"), ("exact", "strong"),
        ("partial", "supported"), ("review", "ambiguous_review")]),
    ("cds_boundary_precision_summary", "CDS-boundary precision", [
        ("exact", "strong"), ("codon_split_one_side", "supported"),
        ("codon_split_both", "supported"), ("unknown", "pending_missing")]),
    ("native_coordinate_sanity", "native coord sanity", [
        ("same_native", "strong"), ("moderate", "supported"),
        ("major", "ambiguous_review")]),
    ("normalized_slot_sanity", "normalized III-slot", [
        ("same_normalized", "strong"), ("minor", "supported"),
        ("major", "ambiguous_review")]),
    ("iii_region_similarity_class", "III-region similarity", [
        ("full_window_distinct", "strong"), ("window_distinct", "supported"),
        ("nearly_identical", "ambiguous_review"), ("ambiguous", "ambiguous_review")]),
    ("final_display_class", "final display class", [
        ("high_confidence", "strong"), ("minor_review", "supported"),
        ("supplementary", "ambiguous_review")]),
    ("interpro_status", "InterProScan input", [
        ("ready", "strong"), ("pending", "supported"), ("missing", "pending_missing")]),
]


def classify_state(column_rules: List[Tuple[str, str]], value: str) -> str:
    v = (value or "").strip().lower()
    if not v or v in ("na", "none", "unknown", "nan"):
        return "pending_missing"
    for sub, state in column_rules:
        if sub in v:
            return state
    return "supported"


# ---------------------------------------------------------------------------
# FIGURE 1 — framework overview
# ---------------------------------------------------------------------------
def fig1_framework(data: FinalData, figdir: Path, tabledir: Path) -> Dict[str, object]:
    n_species = len(data.species_order)
    n_resolved = len(data.coords)
    n_main = sum(1 for r in data.master if str(r["final_display_class"]).startswith("main_analysis"))
    n_supp = sum(1 for r in data.master if is_review_class(r["final_display_class"]))
    n_selected = _to_int(data.interpro.get("total_selected_proteins"), 0)
    n_unique = _to_int(data.interpro.get("unique_sequences"), 0)
    n_dup = _to_int(data.interpro.get("duplicates_collapsed"), 0)

    # evidence layers (Part F): colored side labels grouping the workflow stages
    layers = [
        ("annotation evidence", "#4C9F70", [
            f"{n_species} vertebrate species \u2192 transcript / CDS extraction",
            "IIIb/IIIc mutually exclusive event detection"]),
        ("sequence evidence", "#0072B2", [
            "FGFR2 orthology / paralog validation",
            "sequence-calibrated IIIb/IIIc direction (human-anchored)"]),
        ("protein QC", "#E69F00", [
            "protein marker QC (never auto-swaps IIIb/IIIc labels)"]),
        ("coordinate resolver", "#7B5EA7", [
            "exon-to-protein coordinate resolver (+ CDS-boundary audit)",
            "main / review classification"]),
        ("InterPro-ready output", "#117733", [
            "non-redundant InterProScan-ready FASTA + mapping"]),
    ]
    stages = [s for _, _, ss in layers for s in ss]
    badges = [
        ("species", n_species),
        ("resolved IIIb/IIIc mappings", n_resolved),
        ("main-analysis species", n_main),
        ("supplement / review species", n_supp),
        ("selected proteins", n_selected),
        ("unique InterProScan-ready sequences", n_unique),
        ("duplicate sequences collapsed", n_dup),
    ]

    # plotting table
    rows: List[Dict[str, object]] = []
    order = 1
    for layer_name, _, ss in layers:
        for s in ss:
            rows.append({"item_type": "stage", "order": order, "label": s,
                         "value": "", "evidence_layer": layer_name})
            order += 1
    for lab, val in badges:
        rows.append({"item_type": "badge", "order": "", "label": lab, "value": val,
                     "evidence_layer": ""})
    table = tabledir / "figure1_framework_counts_pre_interpro.tsv"
    write_tsv(table, rows, ["item_type", "order", "label", "value", "evidence_layer"])

    n_stages = len(stages)
    fig, ax = plt.subplots(figsize=(12.5, 8.5))
    ax.set_xlim(-1.4, 12)
    ax.set_ylim(0, n_stages + 2.2)
    ax.axis("off")

    # message banner
    ax.add_patch(FancyBboxPatch((0.2, n_stages + 1.0), 11.6, 0.95,
                                boxstyle="round,pad=0.02,rounding_size=0.1",
                                facecolor="#F2F7FB", edgecolor=C_AXIS, linewidth=1.0))
    ax.text(6.0, n_stages + 1.47,
            "Annotation-aware exon-to-protein boundary framework up to InterProScan preparation",
            ha="center", va="center", fontsize=12, fontweight="bold", color="#1A1A1A")

    side_x, bar_w = 0.55, 0.42
    box_x, box_w = 1.15, 5.45
    y = n_stages
    for layer_name, col, ss in layers:
        y_top = y + 0.0
        y_bot = y - (len(ss) - 1)
        # colored side label spanning the layer's stages
        ax.add_patch(Rectangle((side_x, y_bot - 0.36), bar_w, (y_top - y_bot) + 0.72,
                               facecolor=col, edgecolor="none", zorder=3))
        ax.text(side_x - 0.18, (y_top + y_bot) / 2, layer_name.replace(" ", "\n", 1),
                rotation=90, ha="center", va="center", fontsize=7.3, fontweight="bold",
                color=col, linespacing=0.9)
        for k, s in enumerate(ss):
            yy = y - k
            ax.add_patch(FancyBboxPatch((box_x, yy - 0.34), box_w, 0.68,
                                        boxstyle="round,pad=0.02,rounding_size=0.12",
                                        facecolor="white", edgecolor=col, linewidth=1.2))
            ax.text(box_x + 0.18, yy, s, ha="left", va="center", fontsize=9, color="#1A1A1A")
        y = y_bot - 1
        if y >= 0.5:  # arrow to next layer
            ax.annotate("", xy=(box_x + box_w / 2, y + 0.34), xytext=(box_x + box_w / 2, y_bot - 0.34),
                        arrowprops=dict(arrowstyle="-|>", color="#999999", lw=1.4))

    # badges column
    bx = 7.4
    ax.text(bx, n_stages + 0.5, "Count badges", fontsize=11, fontweight="bold", color=C_AXIS)
    for j, (lab, val) in enumerate(badges):
        yy = n_stages - 0.1 - j * (n_stages / max(1, len(badges)))
        ax.add_patch(FancyBboxPatch((bx, yy - 0.42), 4.2, 0.84,
                                    boxstyle="round,pad=0.02,rounding_size=0.1",
                                    facecolor="white", edgecolor=C_IIIB, linewidth=1.1))
        ax.text(bx + 0.2, yy, str(val), fontsize=15, fontweight="bold", va="center", color=C_IIIB)
        ax.text(bx + 1.05, yy, lab, fontsize=8.2, va="center", color="#333333")
    ax.text(6.0, 0.1, "InterProScan pending \u2013 no real InterPro domains shown",
            ha="center", va="bottom", fontsize=8.5, color="#7A5A00")
    paths = _save(fig, figdir, "Figure_1_framework_pre_interpro")
    return {"figure_id": "Figure_1_framework", "paths": paths, "source_tables": table.name,
            "n_species": n_species, "n_rows": len(rows),
            "main_message": "Annotation-aware framework, not a single plotting script."}


def _cassette_block_index(blocks: List[Dict[str, str]], cassette_rank: Optional[int],
                          native_start: Optional[int], native_end: Optional[int]) -> Optional[int]:
    """Identify which CDS block is the resolved cassette by its UNIQUE cds_rank
    (from the coordinate-overlap cassette map). cds_rank is unique within a transcript,
    unlike repeated NCBI cds_id values. Fall back to protein-coordinate overlap."""
    if cassette_rank is not None:
        for i, b in enumerate(blocks):
            if (_to_int(b.get("cds_rank")) == cassette_rank):
                return i
    if native_start is None or native_end is None:
        return None
    best_i, best_ov = None, 0
    for i, b in enumerate(blocks):
        ps, pe = _to_int(b.get("protein_start_aa")), _to_int(b.get("protein_end_aa"))
        if ps is None or pe is None:
            continue
        ov = max(0, min(pe, native_end) - max(ps, native_start))
        if ov > best_ov:
            best_i, best_ov = i, ov
    return best_i


def _architecture_table(data: FinalData):
    """Return (track_rows, completeness_rows). One track row per real CDS block per
    (species, isoform); a single fallback row where exon blocks are unavailable."""
    tracks: List[Dict[str, object]] = []
    completeness: List[Dict[str, object]] = []
    for r in data.master:
        sp = r["species"]
        review = is_review_class(r.get("final_display_class", ""))
        rrs = r.get("review_reason_short", "")
        for iso in ("IIIb", "IIIc"):
            c = data.coord(sp, iso) or {}
            blocks = data.cds_blocks(sp, iso)
            cm = data.cassette(sp, iso)
            ns = _to_int(cm.get("matched_protein_start_aa")) or _to_int(c.get("native_protein_start_aa"))
            ne = _to_int(cm.get("matched_protein_end_aa")) or _to_int(c.get("native_protein_end_aa"))
            cass_i = _cassette_block_index(blocks, data.cassette_rank(sp, iso), ns, ne)
            left_p = c.get("cds_left_boundary_precision", "")
            right_p = c.get("cds_right_boundary_precision", "")
            base = {
                "species": sp, "display_species_name": r.get("display_species_name", ""),
                "taxon_group": r.get("taxon_group", ""), "phylo_order": r.get("phylo_order", ""),
                "isoform": iso, "track_type": iso,
                "cassette_overlap_status": cm.get("cassette_overlap_status", ""),
                "is_review": "true" if review else "false", "review_reason_short": rrs,
            }
            if blocks:
                completeness.append({"species": sp, "isoform": iso, "n_cds_blocks": len(blocks),
                                     "has_real_blocks": "true",
                                     "completeness_status": "cds_blocks_reconstructed"})
                for j, b in enumerate(blocks):
                    is_cass = (j == cass_i)
                    iiib_c = is_cass and iso == "IIIb"
                    iiic_c = is_cass and iso == "IIIc"
                    cc = ("IIIb_cassette" if iiib_c else "IIIc_cassette" if iiic_c else "common_cds")
                    tracks.append({**base,
                        "feature_type": "CDS",
                        "exon_or_cds_id": b.get("cds_id_source", ""),
                        "protein_start_aa": b.get("protein_start_aa", ""),
                        "protein_end_aa": b.get("protein_end_aa", ""),
                        "feature_label": (f"{iso} cassette" if is_cass else f"CDS{b.get('cds_rank','')}"),
                        "feature_color_class": cc,
                        "is_IIIb_cassette": "true" if iiib_c else "false",
                        "is_IIIc_cassette": "true" if iiic_c else "false",
                        "boundary_left_precision": left_p if is_cass else "",
                        "boundary_right_precision": right_p if is_cass else "",
                    })
            else:
                completeness.append({"species": sp, "isoform": iso, "n_cds_blocks": 0,
                                     "has_real_blocks": "false",
                                     "completeness_status": "partial_protein_bar_no_cds_blocks"})
                plen = _to_int(c.get("protein_length_aa"), 0) or 0
                tracks.append({**base,
                    "feature_type": "partial_protein_bar",
                    "exon_or_cds_id": "", "protein_start_aa": 1, "protein_end_aa": plen,
                    "feature_label": "CDS architecture partially reconstructed",
                    "feature_color_class": "partial",
                    "is_IIIb_cassette": "false", "is_IIIc_cassette": "false",
                    "boundary_left_precision": "", "boundary_right_precision": "",
                })
                # cassette overlay row so the cassette is still visible
                if ns is not None and ne is not None and ne > ns:
                    tracks.append({**base,
                        "feature_type": "cassette_overlay",
                        "exon_or_cds_id": cm.get("matched_unique_cds_block_id", ""),
                        "protein_start_aa": ns, "protein_end_aa": ne,
                        "feature_label": f"{iso} cassette",
                        "feature_color_class": f"{iso}_cassette",
                        "is_IIIb_cassette": "true" if iso == "IIIb" else "false",
                        "is_IIIc_cassette": "true" if iso == "IIIc" else "false",
                        "boundary_left_precision": left_p, "boundary_right_precision": right_p,
                    })
    return tracks, completeness


BOUNDARY_MARK = {"codon_boundary_exact": ("|", "#117733"),
                 "codon_split_one_side": ("/", "#999999"),
                 "codon_split_both_sides": ("x", "#999999"),
                 "unknown_codon_phase": ("?", "#CC6677")}




# Refined, calm edge symbols (Part E): split/phase flags are subtle and grey, never
# alarming; only true unresolved boundaries get a warning color.
REFINED_EDGE = {
    "codon_boundary_exact": ("", "#117733", 6.0),
    "known_split_codon_boundary": ("/", "#9AA7B0", 6.0),
    "phase_not_available_but_coordinate_resolved": ("\u00b7", "#9DB4C4", 8.0),
    "nucleotide_sequence_unavailable": ("", "#AAAAAA", 6.0),
    "boundary_unresolved": ("?", "#CC6677", 6.5),
}
VIS_MARK = {"main_warning": ("\u26A0", "#B5651D"), "hard_fail": ("\u2715", "#CC3311")}


def _refined_edge_symbol(bstate: str):
    return REFINED_EDGE.get(bstate or "", ("", "#999999", 6.0))


def _render_architecture(data: FinalData, figdir: Path, stem: str, title: str,
                         subtitle: str) -> Dict[str, str]:
    master = data.master
    n = len(master)
    max_len = 900
    for r in master:
        for iso in ("IIIb", "IIIc"):
            c = data.coord(r["species"], iso)
            if c:
                max_len = max(max_len, _to_int(c.get("protein_length_aa"), 0) or 0)

    fig, ax = plt.subplots(figsize=(13.5, max(9.5, 0.56 * n + 2.5)))
    lbl_w = max_len * 0.30          # left label gutter
    qc_w = max_len * 0.16           # right QC strip
    _x0, x1 = 0, max_len
    draw_taxon_bands(ax, master, lambda i: (n - i), -lbl_w * 0.05, x1 + qc_w, label=False)

    th = 0.30
    for i, r in enumerate(master):
        sp = r["species"]
        yc = (n - i)
        review = is_review_class(r.get("final_display_class", ""))
        _minor = is_minor_review(r.get("final_display_class", ""))
        disp = r.get("display_species_name", sp)
        # LEFT species label (Part B/D)
        ax.text(-lbl_w * 0.04, yc, disp, ha="right", va="center", fontsize=7.6,
                style="italic", color=("#7A5A00" if review else "#1A1A1A"))
        for iso, yo in (("IIIb", 0.22), ("IIIc", -0.22)):
            c = data.coord(sp, iso) or {}
            yy = yc + yo
            blocks = data.cds_blocks(sp, iso)
            cm = data.cassette(sp, iso)
            ns = _to_int(cm.get("matched_protein_start_aa")) or _to_int(c.get("native_protein_start_aa"))
            ne = _to_int(cm.get("matched_protein_end_aa")) or _to_int(c.get("native_protein_end_aa"))
            cass_i = _cassette_block_index(blocks, data.cassette_rank(sp, iso), ns, ne)
            ax.text(-lbl_w * 0.005, yy, iso, ha="right", va="center", fontsize=5.6,
                    color=(C_IIIB if iso == "IIIb" else C_IIIC))
            rf = data.refined(sp, iso)
            vis = rf.get("plot_visibility_level", "subtle_symbol")
            bstate = rf.get("boundary_precision_state", "")
            prominent = vis in ("main_warning", "hard_fail")
            if blocks:
                for j, b in enumerate(blocks):
                    ps, pe = _to_int(b.get("protein_start_aa")), _to_int(b.get("protein_end_aa"))
                    if ps is None or pe is None or pe <= ps:
                        continue
                    is_cass = (j == cass_i)
                    col = (C_IIIB if (is_cass and iso == "IIIb") else
                           C_IIIC if (is_cass and iso == "IIIc") else C_COMMON)
                    ax.add_patch(Rectangle((ps, yy - th / 2), pe - ps, th, facecolor=col,
                                           edgecolor="white", linewidth=0.6,
                                           linestyle=("--" if (is_cass and prominent) else "-"),
                                           zorder=4 if is_cass else 3))
                    if is_cass:
                        # subtle edge symbols for boundary precision (split/phase = not errors)
                        sym, scol, sz = _refined_edge_symbol(bstate)
                        if sym:
                            for bx in (ps, pe):
                                ax.text(bx, yy + th / 2 + 0.10, sym, ha="center", va="bottom",
                                        fontsize=sz, color=scol, zorder=6)
                        # strong warning ONLY for true biological/annotation issues
                        if prominent:
                            mk, mcol = VIS_MARK.get(vis, ("\u26A0", "#B5651D"))
                            ax.text((ps + pe) / 2, yy + th / 2 + 0.12, mk, ha="center",
                                    va="bottom", fontsize=8.5, color=mcol, zorder=7)
            else:
                plen = _to_int(c.get("protein_length_aa"), 0) or 0
                ax.add_patch(Rectangle((0, yy - th / 2), plen, th, facecolor="#EDEDED",
                                       edgecolor="#CFCFCF", linewidth=0.4, hatch="///", zorder=3))
                if ns is not None and ne is not None and ne > ns:
                    col = C_IIIB if iso == "IIIb" else C_IIIC
                    ax.add_patch(Rectangle((ns, yy - th / 2), ne - ns, th, facecolor=col,
                                           edgecolor=("#7A5A00" if iso == "IIIc" else "#003E63"),
                                           linewidth=1.0, linestyle="--", zorder=4))
                    # protein-overlay-only = true missing CDS model -> prominent if main
                    if prominent:
                        mk, mcol = VIS_MARK.get(vis, ("\u26A0", "#B5651D"))
                        ax.text((ns + ne) / 2, yy + th / 2 + 0.12, mk, ha="center",
                                va="bottom", fontsize=8.5, color=mcol, zorder=7)
        # right QC strip
        fdc = r.get("final_display_class", "")
        strip = ("HC" if "high_confidence" in fdc else "rev" if review else "minor")
        scol = ("#117733" if strip == "HC" else C_REVIEW if strip == "rev" else "#DDCC77")
        ax.text(x1 + qc_w * 0.55, yc, strip, ha="center", va="center", fontsize=6.5,
                color=scol, fontweight="bold")
        if review:
            ax.text(x1 + qc_w * 0.12, yc, "\u25C6", ha="center", va="center", fontsize=7, color=C_REVIEW)

    # taxon group labels at far left
    for g, a, b in group_spans(master):
        yt, yb = (n - a) + 0.5, (n - b) - 0.5
        ax.text(-lbl_w * 0.85, (yt + yb) / 2, g, rotation=90, va="center", ha="center",
                fontsize=8, fontweight="bold", color=TAXON_BAR.get(g, "#666"))

    ax.set_xlim(-lbl_w, x1 + qc_w * 1.2)
    ax.set_ylim(0.2, n + 1.2)
    ax.set_yticks([])
    ax.set_xlabel("protein coordinate (amino-acid position)", fontsize=10)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=26)
    if subtitle:
        ax.text(0.5, 1.008, subtitle, transform=ax.transAxes, ha="center", va="bottom",
                fontsize=8.5, color="#444")
    legend = [
        Line2D([0], [0], color=C_COMMON, lw=9, label="common CDS/exon blocks"),
        Line2D([0], [0], color=C_IIIB, lw=9, label="IIIb cassette"),
        Line2D([0], [0], color=C_IIIC, lw=9, label="IIIc cassette"),
        Line2D([0], [0], marker="$/$", color="#9AA7B0", markersize=8, lw=0,
               label="split codon boundary (not an error)"),
        Line2D([0], [0], marker="$\u00b7$", color="#9DB4C4", markersize=10, lw=0,
               label="phase unavailable, coordinate resolved"),
        Line2D([0], [0], marker="$\u26A0$", color="#B5651D", markersize=9, lw=0,
               label="true issue (overlay/conflict/offset; dashed cassette)"),
        Line2D([0], [0], color="#EDEDED", lw=9, label="CDS architecture partially reconstructed"),
    ]
    ax.legend(handles=legend, loc="lower right", fontsize=6.6, frameon=True, framealpha=0.92, ncol=1)
    ax.text(0.99, -0.05, "InterProScan pending (no real InterPro domains shown)",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#7A5A00")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    return _save(fig, figdir, stem)


def fig2_architecture(data: FinalData, figdir: Path, tabledir: Path) -> Dict[str, object]:
    tracks, completeness = _architecture_table(data)
    table = tabledir / "figure2_exon_to_protein_architecture_tracks.tsv"
    write_tsv(table, tracks, list(tracks[0].keys()))
    comp = tabledir / "figure2_exon_block_completeness.tsv"
    write_tsv(comp, completeness, list(completeness[0].keys()))
    paths = _render_architecture(
        data, figdir, "Figure_2_exon_to_protein_architecture_pre_interpro",
        "Exon/CDS-to-protein architecture of FGFR2 IIIb/IIIc cassettes",
        "Real CDS-derived protein blocks per transcript; resolved IIIb/IIIc cassette highlighted; "
        "phylogenetic order; InterProScan pending")
    n_real = sum(1 for r in completeness if r["has_real_blocks"] == "true")
    return {"figure_id": "Figure_2_exon_to_protein_architecture", "paths": paths,
            "source_tables": f"{table.name};{comp.name}", "n_species": len(data.species_order),
            "n_rows": len(tracks),
            "main_message": f"Real exon/CDS-derived protein architecture across vertebrate FGFR2 "
                            f"orthologs ({n_real}/60 transcripts with reconstructed CDS blocks); "
                            f"IIIb/IIIc cassettes highlighted with codon-boundary precision."}


def supplement1_all_native(data: FinalData, figdir: Path, tabledir: Path) -> Dict[str, object]:
    tracks, _completeness = _architecture_table(data)
    table = tabledir / "figure2_exon_to_protein_architecture_tracks.tsv"
    if not table.exists():
        write_tsv(table, tracks, list(tracks[0].keys()))
    paths = _render_architecture(
        data, figdir, "Supplement_Figure_1_all_species_native_tracks_pre_interpro",
        "Supplement 1: all-species exon/CDS-to-protein tracks",
        "Same colors / phylogenetic order / review markers as Figure 2; InterProScan pending")
    return {"figure_id": "Supplement_Figure_1_all_species_native_tracks", "paths": paths,
            "source_tables": table.name, "n_species": len(data.species_order),
            "n_rows": len(tracks),
            "main_message": "All 30 species shown with real CDS-block architecture in protein space."}


# ---------------------------------------------------------------------------
# FIGURE 3 — IgIII / D3 cassette zoom (normalized III-slot coordinates)
# ---------------------------------------------------------------------------
# Refined boundary-precision symbols for Figure 3 (calm; split/phase are not errors).
FIG3_SYM = {
    "codon_boundary_exact": ("\u25CF", "#117733"),                       # ● exact
    "known_split_codon_boundary": ("\u25D0", "#6B7B8C"),                 # ◐ split codon
    "phase_not_available_but_coordinate_resolved": ("\u25CB", "#9DB4C4"),# ○ phase unavailable
    "nucleotide_sequence_unavailable": ("\u25A2", "#B5651D"),            # ▢ overlay/missing
    "boundary_unresolved": ("?", "#CC6677"),
}
QC_ICON = {"validated": ("\u2713", "#117733"), "conflict": ("\u2717", "#CC6677"),
           "ambiguous": ("~", "#DDCC77"), "review": ("~", "#DDCC77")}


def _qc_icon(summary: str):
    s = (summary or "").lower()
    for key, val in QC_ICON.items():
        if key in s:
            return val
    return ("", "#444")


def _igIII_segments(c: Dict[str, str]):
    """Derive upstream_conserved / cassette / downstream_conserved within the IgIII
    (III_region) window, normalized so the window starts at 0. Uses real coordinates;
    returns [] if window/cassette coordinates are unavailable."""
    ws, we = _to_int(c.get("III_region_start_aa")), _to_int(c.get("III_region_end_aa"))
    cs, ce = _to_int(c.get("native_protein_start_aa")), _to_int(c.get("native_protein_end_aa"))
    if None in (ws, we, cs, ce) or we <= ws:
        return [], "window_unavailable"
    # clamp cassette into window
    cas_s, cas_e = max(ws, cs), min(we, ce)
    conf = "derived_from_III_region_window"
    if cs < ws or ce > we:
        conf = "cassette_extends_beyond_window_review"
    segs = []
    if cas_s > ws:
        segs.append(("upstream_conserved", 0, cas_s - ws))
    segs.append(("cassette", cas_s - ws, cas_e - ws))
    if we > cas_e:
        segs.append(("downstream_conserved", cas_e - ws, we - ws))
    return segs, conf


def fig3_cassette_zoom(data: FinalData, figdir: Path, tabledir: Path) -> Dict[str, object]:
    master = data.master
    n = len(master)
    rows: List[Dict[str, object]] = []
    max_w = 1
    seg_cache: Dict[Tuple[str, str], list] = {}
    for r in master:
        for iso in ("IIIb", "IIIc"):
            c = data.coord(r["species"], iso) or {}
            segs, conf = _igIII_segments(c)
            seg_cache[(r["species"], iso)] = (segs, conf)
            if segs:
                max_w = max(max_w, max(e for _, _, e in segs))
            for stype, s0, s1 in (segs or [("cassette", 0, 0)]):
                rows.append({
                    "species": r["species"], "isoform": iso, "segment_type": stype,
                    "normalized_start_aa": s0, "normalized_end_aa": s1,
                    "segment_confidence": conf,
                    "cds_boundary_precision_refined": c.get("cds_boundary_precision_refined", ""),
                    "protein_validation_summary": r.get("protein_validation_summary", ""),
                    "iii_region_similarity_class": r.get("iii_region_similarity_class", ""),
                    "is_review": "true" if is_review_class(r.get("final_display_class", "")) else "false",
                    "review_reason_short": r.get("review_reason_short", ""),
                })
    table = tabledir / "figure3_igIII_cassette_zoom_tracks.tsv"
    write_tsv(table, rows, list(rows[0].keys()))

    fig, ax = plt.subplots(figsize=(11, max(9.5, 0.56 * n + 2.5)))
    lbl_w = max_w * 0.42
    th = 0.30
    draw_taxon_bands(ax, master, lambda i: (n - i), -lbl_w * 0.05, max_w, label=False)
    for i, r in enumerate(master):
        sp = r["species"]
        yc = (n - i)
        review = is_review_class(r.get("final_display_class", ""))
        disp = r.get("display_species_name", sp)
        ax.text(-lbl_w * 0.04, yc, disp, ha="right", va="center", fontsize=7.6,
                style="italic", color=("#7A5A00" if review else "#1A1A1A"))
        for iso, yo in (("IIIb", 0.22), ("IIIc", -0.22)):
            c = data.coord(sp, iso) or {}
            segs, conf = seg_cache[(sp, iso)]
            yy = yc + yo
            ax.text(-lbl_w * 0.005, yy, iso, ha="right", va="center", fontsize=5.6,
                    color=(C_IIIB if iso == "IIIb" else C_IIIC))
            cas_col = C_IIIB if iso == "IIIb" else C_IIIC
            for stype, s0, s1 in segs:
                if s1 <= s0:
                    continue
                col = cas_col if stype == "cassette" else C_COMMON
                ax.add_patch(Rectangle((s0, yy - th / 2), s1 - s0, th, facecolor=col,
                                       edgecolor=("#7A5A00" if (stype == "cassette" and iso == "IIIc")
                                                  else "#003E63" if stype == "cassette" else "white"),
                                       linewidth=(1.2 if (stype == "cassette" and review) else 0.6),
                                       linestyle=("--" if (stype == "cassette" and review) else "-"),
                                       zorder=4 if stype == "cassette" else 3))
                if stype == "cassette":
                    # boundary precision symbol from the REFINED boundary state
                    rf = data.refined(sp, iso)
                    fsym, fcol = FIG3_SYM.get(rf.get("boundary_precision_state", ""), ("", "#444"))
                    ax.text(s1 + max_w * 0.012, yy, fsym,
                            ha="left", va="center", fontsize=7.0, color=fcol)
            # protein QC small icon at far right of the track
            sym, scol = _qc_icon(r.get("protein_validation_summary", ""))
            if sym:
                ax.text(max_w * 1.04, yy, sym, ha="center", va="center", fontsize=6.5, color=scol)
        if review:
            ax.text(max_w * 1.10, yc, "\u25C6", ha="center", va="center", fontsize=7, color=C_REVIEW)

    for g, a, b in group_spans(master):
        yt, yb = (n - a) + 0.5, (n - b) - 0.5
        ax.text(-lbl_w * 0.92, (yt + yb) / 2, g, rotation=90, va="center", ha="center",
                fontsize=8, fontweight="bold", color=TAXON_BAR.get(g, "#666"))

    ax.set_xlim(-lbl_w, max_w * 1.16)
    ax.set_ylim(0.2, n + 1.0)
    ax.set_yticks([])
    ax.set_xlabel("normalized IgIII-window coordinate (aa; window start = 0)", fontsize=10)
    ax.set_title("IgIII/D3 IIIb/IIIc event zoom: upstream \u2013 cassette \u2013 downstream",
                 fontsize=12.5, fontweight="bold", pad=26)
    ax.text(0.5, 1.008, "Boundary symbols show codon-precision confidence; split codons are not errors.",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=8.5, color="#444")
    legend = [
        Line2D([0], [0], color=C_COMMON, lw=9, label="upstream/downstream conserved"),
        Line2D([0], [0], color=C_IIIB, lw=9, label="IIIb cassette"),
        Line2D([0], [0], color=C_IIIC, lw=9, label="IIIc cassette"),
        Line2D([0], [0], marker="$\u25CF$", color="#117733", markersize=8, lw=0,
               label="exact boundary"),
        Line2D([0], [0], marker="$\u25D0$", color="#6B7B8C", markersize=8, lw=0,
               label="split codon boundary (not an error)"),
        Line2D([0], [0], marker="$\u25CB$", color="#9DB4C4", markersize=8, lw=0,
               label="phase unavailable but coordinate resolved"),
        Line2D([0], [0], marker="$\u25A2$", color="#B5651D", markersize=8, lw=0,
               label="protein overlay only (no CDS model)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=C_REVIEW, markersize=8,
               label="review case"),
    ]
    ax.legend(handles=legend, loc="lower right", fontsize=6.6, frameon=True, framealpha=0.92)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    paths = _save(fig, figdir, "Figure_3_IgIII_cassette_zoom_pre_interpro")
    return {"figure_id": "Figure_3_IgIII_cassette_zoom", "paths": paths,
            "source_tables": table.name, "n_species": n, "n_rows": len(rows),
            "main_message": "Upstream\u2013cassette\u2013downstream IgIII event structure; normalized "
                            "cassette-internal mapping is stable across vertebrate FGFR2 orthologs."}


# ---------------------------------------------------------------------------
# FIGURE 4 — species evidence matrix / QC heatmap
# ---------------------------------------------------------------------------
def fig4_evidence_matrix(data: FinalData, figdir: Path, tabledir: Path) -> Dict[str, object]:
    master = data.master
    n = len(master)
    layers = EVIDENCE_LAYERS
    rows: List[Dict[str, object]] = []
    grid: List[List[str]] = []
    refined_states: List[str] = []
    for r in master:
        states = []
        out = {"phylo_order": r.get("phylo_order", ""), "species": r["species"],
               "display_species_name": r.get("display_species_name", ""),
               "taxon_group": r.get("taxon_group", "")}
        for col, label, rules in layers:
            st = classify_state(rules, r.get(col, ""))
            states.append(st)
            out[label] = st
        # refined uncertainty (worst-case across IIIb/IIIc); split/phase stay calm
        rfb = data.refined(r["species"], "IIIb")
        rfc = data.refined(r["species"], "IIIc")
        worst = max([rfb, rfc], key=lambda x: DISPLAY_RANK_FIG4.get(
            x.get("display_uncertainty_class", ""), 0))
        rstate = refined_fig4_state(worst.get("display_uncertainty_class", ""),
                                    worst.get("protein_evidence_state", ""))
        refined_states.append(rstate)
        out["display_uncertainty_class"] = worst.get("display_uncertainty_class", "")
        out["refined_uncertainty_state"] = rstate
        out["plot_visibility_level"] = worst.get("plot_visibility_level", "")
        rows.append(out)
        grid.append(states)
    table = tabledir / "figure4_species_evidence_matrix.tsv"
    fields = (["phylo_order", "species", "display_species_name", "taxon_group"]
              + [l for _, l, _ in layers]
              + ["display_uncertainty_class", "refined_uncertainty_state", "plot_visibility_level"])
    write_tsv(table, rows, fields)

    ncol = len(layers)
    # counts for the top annotation (read from canonical master; not recomputed)
    n_main = sum(1 for r in master if str(r["final_display_class"]).startswith("main_analysis"))
    n_supp = sum(1 for r in master if is_review_class(r["final_display_class"]))
    n_resolved = len(data.coords)
    lbl_w = 9.0  # left gutter (matrix x-units) for species + group labels

    fig, ax = plt.subplots(figsize=(max(12, 0.95 * ncol + 6), max(8.5, 0.42 * n + 3)))

    def y_of(i: int) -> float:
        return (n - i)

    refined_col_x = ncol + 0.6  # gap, then refined-uncertainty column
    for i in range(n):
        yc = y_of(i)
        for j in range(ncol):
            st = grid[i][j]
            ax.add_patch(Rectangle((j, yc - 0.45), 0.92, 0.9, facecolor=STATE_COLORS[st],
                                   edgecolor="white", linewidth=0.8, zorder=3))
        # refined-uncertainty column (calm colors; split/phase not alarming)
        rcol = REFINED_STATE_COLORS.get(refined_states[i], "#BBBBBB")
        ax.add_patch(Rectangle((refined_col_x, yc - 0.45), 0.92, 0.9, facecolor=rcol,
                               edgecolor="white", linewidth=0.8, zorder=3))
        # species name on the LEFT (Part D)
        ax.text(-0.4, yc, master[i].get("display_species_name", master[i]["species"]),
                ha="right", va="center", fontsize=7.6, style="italic")
    # taxon separators + group labels further left
    for gi, (g, a, b) in enumerate(group_spans(master)):
        if gi > 0:
            ax.axhline(y_of(a) + 0.5, color="#888888", lw=0.9, zorder=4)
        ax.text(-lbl_w + 0.3, (y_of(a) + y_of(b)) / 2, g, rotation=90, va="center", ha="center",
                fontsize=8, fontweight="bold", color=TAXON_BAR.get(g, "#666"))
    # column labels (mildly angled)
    for j, (_, label, _) in enumerate(layers):
        ax.text(j + 0.46, n + 0.6, label, rotation=35, ha="left", va="bottom", fontsize=8)
    ax.text(refined_col_x + 0.46, n + 0.6, "refined uncertainty", rotation=35, ha="left",
            va="bottom", fontsize=8, fontweight="bold", color="#333")
    # top count annotation
    ax.text(-lbl_w, n + 2.0,
            f"30 species   |   {n_main} main-analysis   |   {n_supp} supplement/review   "
            f"|   {n_resolved} resolved IIIb/IIIc mappings   |   InterProScan pending",
            ha="left", va="bottom", fontsize=9, fontweight="bold", color="#333")
    ax.set_xlim(-lbl_w, refined_col_x + 6.0)
    ax.set_ylim(0.3, n + 2.6)
    ax.axis("off")
    ax.set_title("Species evidence matrix / QC heatmap (pre-InterPro)",
                 fontsize=12.5, fontweight="bold", loc="left")
    legend = [Line2D([0], [0], marker="s", color="w", markerfacecolor=STATE_COLORS[s],
                     markersize=11, label=STATE_LABELS[s]) for s in
              ("strong", "supported", "transcript_only", "ambiguous_review", "pending_missing")]
    leg1 = ax.legend(handles=legend, loc="upper right", fontsize=7.5, frameon=True, framealpha=0.9,
                     ncol=1, bbox_to_anchor=(1.0, 1.0), title="evidence-layer state")
    ax.add_artist(leg1)
    # explicit refined-uncertainty legend (split/phase are NOT errors; review separate)
    rlegend = [Line2D([0], [0], marker="s", color="w", markerfacecolor=REFINED_STATE_COLORS[s],
                      markersize=11, label=REFINED_STATE_LABELS[s]) for s in
               ("robust_pass", "split_codon", "phase_unavailable", "transcript_protein_only",
                "review", "hard_fail", "interpro_pending")]
    ax.legend(handles=rlegend, loc="lower right", fontsize=7.0, frameon=True, framealpha=0.92,
              ncol=1, bbox_to_anchor=(1.0, 0.0), title="refined uncertainty (how to read it)")
    paths = _save(fig, figdir, "Figure_4_species_evidence_matrix_pre_interpro")
    return {"figure_id": "Figure_4_species_evidence_matrix", "paths": paths,
            "source_tables": table.name, "n_species": n, "n_rows": len(rows),
            "main_message": "Integrates multiple independent evidence layers and a refined-uncertainty "
                            "column that separates split-codon, phase-unavailable, transcript/protein-only, "
                            "review and hard-fail states; minor boundary flags stay calm."}


# ---------------------------------------------------------------------------
# FIGURE 5 — native vs normalized coordinate QC
# ---------------------------------------------------------------------------
FDC_STYLE = {
    "main_analysis_high_confidence": ("#117733", "o", "main: high confidence"),
    "main_analysis_with_minor_review": ("#DDCC77", "s", "main: minor review"),
    "supplementary_review_not_primary_claim": (C_REVIEW, "D", "supplement / review"),
}


def fig5_native_vs_normalized(data: FinalData, figdir: Path, tabledir: Path) -> Dict[str, object]:
    rows: List[Dict[str, object]] = []
    for r in data.master:
        cb = data.coord(r["species"], "IIIb") or {}
        cc = data.coord(r["species"], "IIIc") or {}
        nb, nc = _to_int(cb.get("native_protein_center_aa")), _to_int(cc.get("native_protein_center_aa"))
        sb, sc = _to_int(cb.get("iii_slot_center_aa")), _to_int(cc.get("iii_slot_center_aa"))
        native_off = abs(nb - nc) if (nb is not None and nc is not None) else ""
        norm_off = abs(sb - sc) if (sb is not None and sc is not None) else ""
        rows.append({
            "species": r["species"], "display_species_name": r.get("display_species_name", ""),
            "phylo_order": r.get("phylo_order", ""), "taxon_group": r.get("taxon_group", ""),
            "native_offset_aa": native_off, "normalized_offset_aa": norm_off,
            "recommended_use": r.get("recommended_use", ""),
            "final_display_class": r.get("final_display_class", ""),
            "is_review": "true" if is_review_class(r.get("final_display_class", "")) else "false",
        })
    table = tabledir / "figure5_native_vs_normalized_coordinate_qc.tsv"
    write_tsv(table, rows, list(rows[0].keys()))

    fig, ax = plt.subplots(figsize=(9.5, 7.5))
    # threshold bands (Part E): native same/minor/moderate/major; normalized same/minor/review
    xs = [r["native_offset_aa"] for r in rows if r["native_offset_aa"] != ""]
    ys = [r["normalized_offset_aa"] for r in rows if r["normalized_offset_aa"] != ""]
    xmax = max(xs + [10]) * 1.12
    ymax = max(ys + [3]) * 1.25
    NAT_MINOR, NAT_MODERATE = 5, 50
    NORM_MINOR = 2
    band_specs = [(0, NAT_MINOR, "#E8F4EA"), (NAT_MINOR, NAT_MODERATE, "#FBF6E3"),
                  (NAT_MODERATE, xmax, "#FBEAEA")]
    for lo, hi, col in band_specs:
        ax.axvspan(lo, hi, color=col, zorder=0)
    ax.axhspan(NORM_MINOR, ymax, color="#CC6677", alpha=0.05, zorder=0)
    for xv, lab in ((NAT_MINOR, "minor"), (NAT_MODERATE, "moderate/major")):
        ax.axvline(xv, color="#CCCCCC", lw=0.7, ls="--", zorder=1)
        ax.text(xv, ymax * 0.98, lab, rotation=90, fontsize=6.5, color="#999",
                ha="right", va="top")
    ax.axhline(NORM_MINOR, color="#CCCCCC", lw=0.7, ls="--", zorder=1)
    ax.text(xmax * 0.99, NORM_MINOR, "normalized review", fontsize=6.5, color="#999",
            ha="right", va="bottom")

    import random
    rng = random.Random(7)
    seen_labels = set()
    for r in rows:
        x = r["native_offset_aa"]
        y = r["normalized_offset_aa"]
        if x == "" or y == "":
            continue
        fdc = r["final_display_class"]
        color, marker, lab = FDC_STYLE.get(fdc, (C_PENDING, "o", fdc))
        # tiny jitter so overlapping zeros remain visible
        jx = x + rng.uniform(-0.6, 0.6)
        jy = y + rng.uniform(-0.12, 0.12)
        ax.scatter([jx], [jy], c=color, marker=marker, s=70, edgecolors="#333333",
                   linewidths=0.6, zorder=4, label=lab if lab not in seen_labels else None)
        seen_labels.add(lab)
        if r["is_review"] == "true":
            ax.annotate(r["display_species_name"], (jx, jy), fontsize=7, style="italic",
                        xytext=(6, 4), textcoords="offset points", color="#7A5A00")
    ax.set_xlim(-xmax * 0.03, xmax)
    ax.set_ylim(-ymax * 0.06, ymax)
    ax.set_xlabel("native IIIb\u2013IIIc cassette center offset (aa)", fontsize=10)
    ax.set_ylabel("normalized III-slot center offset (aa)", fontsize=10)
    ax.set_title("Native vs normalized coordinate QC (pre-InterPro)",
                 fontsize=12.5, fontweight="bold", pad=22)
    ax.text(0.5, 1.008, "Native offsets identify annotation-dependent cases; normalized III-slot "
            "coordinates remain stable for most species.",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=8.3, color="#444")
    # side count summary
    n_native_same = sum(1 for r in rows if r["native_offset_aa"] != "" and r["native_offset_aa"] <= NAT_MINOR)
    n_norm_same = sum(1 for r in rows if r["normalized_offset_aa"] != "" and r["normalized_offset_aa"] <= NORM_MINOR)
    ax.text(0.99, 0.99, f"native \u2264{NAT_MINOR}aa: {n_native_same}/{len(rows)}\n"
            f"normalized \u2264{NORM_MINOR}aa: {n_norm_same}/{len(rows)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#555",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="#DDD", alpha=0.85))
    ax.grid(True, color="#EEEEEE", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=8, frameon=True, framealpha=0.9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    paths = _save(fig, figdir, "Figure_5_native_vs_normalized_coordinate_qc_pre_interpro")
    return {"figure_id": "Figure_5_native_vs_normalized_coordinate_qc", "paths": paths,
            "source_tables": table.name, "n_species": len(rows), "n_rows": len(rows),
            "main_message": "Native offsets flag annotation/isoform-dependent cases; normalized "
                            "III-slot coordinates remain stable for most species."}


# ---------------------------------------------------------------------------
# SUPPLEMENT 2 — review-case panels (read review species from master)
# ---------------------------------------------------------------------------
def supplement2_review_cases(data: FinalData, figdir: Path, tabledir: Path) -> Dict[str, object]:
    review = [r for r in data.master if is_review_class(r.get("final_display_class", ""))]
    rows: List[Dict[str, object]] = []
    for r in review:
        cb = data.coord(r["species"], "IIIb") or {}
        cc = data.coord(r["species"], "IIIc") or {}
        rows.append({
            "species": r["species"], "display_species_name": r.get("display_species_name", ""),
            "taxon_group": r.get("taxon_group", ""),
            "IIIb_native_start": cb.get("native_protein_start_aa", ""),
            "IIIb_native_end": cb.get("native_protein_end_aa", ""),
            "IIIc_native_start": cc.get("native_protein_start_aa", ""),
            "IIIc_native_end": cc.get("native_protein_end_aa", ""),
            "native_coordinate_sanity": r.get("native_coordinate_sanity", ""),
            "normalized_slot_sanity": r.get("normalized_slot_sanity", ""),
            "protein_validation_summary": r.get("protein_validation_summary", ""),
            "iii_region_similarity_class": r.get("iii_region_similarity_class", ""),
            "cds_boundary_precision_summary": r.get("cds_boundary_precision_summary", ""),
            "review_reason_short": r.get("review_reason_short", ""),
        })
    table = tabledir / "supplement_review_cases.tsv"
    write_tsv(table, rows, list(rows[0].keys()) if rows else ["species"])

    nrev = len(review)
    ncols = 2 if nrev > 1 else 1
    nr = (nrev + ncols - 1) // ncols
    # Guard: small/custom runs can have zero drawable review-case groups. Calling
    # plt.subplots(nr=0, ...) raises "Number of rows must be a positive integer,
    # not 0". This is an optional display figure, so skip it gracefully instead of
    # aborting the pre-InterPro pipeline.
    if nrev <= 0 or nr <= 0:
        return _skip_optional_figure(
            figure_id="Supplement_Figure_2_review_cases",
            func_name="supplement2_review_cases",
            figdir=figdir, reason="no_drawable_review_cases_for_this_run",
            log_message="no drawable review cases for this run",
            n_species=nrev, source_tables=table.name,
            main_message="No review cases in this run; optional review-case panel skipped.")
    fig, axes = plt.subplots(nr, ncols, figsize=(11.5, max(3.5, 3.6 * nr)), squeeze=False)
    for idx, r in enumerate(review):
        ax = axes[idx // ncols][idx % ncols]
        cb = data.coord(r["species"], "IIIb") or {}
        cc = data.coord(r["species"], "IIIc") or {}
        plen = max(_to_int(cb.get("protein_length_aa"), 0) or 0,
                   _to_int(cc.get("protein_length_aa"), 0) or 0, 900)
        # cassette bars occupy the upper part of the panel; QC text the lower part
        for iso, yy, c in (("IIIb", 0.86, cb), ("IIIc", 0.70, cc)):
            ax.add_patch(Rectangle((0, yy - 0.05), plen, 0.10, facecolor=C_COMMON,
                                   edgecolor="#AFAFAF", linewidth=0.4))
            ns, ne = _to_int(c.get("native_protein_start_aa")), _to_int(c.get("native_protein_end_aa"))
            if ns is not None and ne is not None and ne > ns:
                col = C_IIIB if iso == "IIIb" else C_IIIC
                ax.add_patch(Rectangle((ns, yy - 0.05), ne - ns, 0.10, facecolor=col,
                                       edgecolor="#333", linewidth=0.6))
                ax.text(ne + plen * 0.01, yy, f"{ns}\u2013{ne}", ha="left", va="center",
                        fontsize=6.2, color="#555")
            ax.text(-plen * 0.01, yy, iso, ha="right", va="center", fontsize=7,
                    color=(C_IIIB if iso == "IIIb" else C_IIIC))
        # light coordinate ruler at y=0.60
        ax.plot([0, plen], [0.60, 0.60], color="#CCCCCC", lw=0.6)
        for xt in range(0, plen + 1, 200):
            ax.plot([xt, xt], [0.585, 0.615], color="#CCCCCC", lw=0.6)
            ax.text(xt, 0.55, str(xt), ha="center", va="top", fontsize=5.8, color="#999")
        ax.set_xlim(-plen * 0.03, plen * 1.08)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(r.get("display_species_name", r["species"]), fontsize=9.5,
                     style="italic", color="#7A5A00")
        info = (f"native: {r.get('native_coordinate_sanity','')}\n"
                f"normalized: {r.get('normalized_slot_sanity','')}\n"
                f"protein QC: {r.get('protein_validation_summary','')}\n"
                f"III similarity: {r.get('iii_region_similarity_class','')}\n"
                f"CDS precision: {r.get('cds_boundary_precision_summary','')}\n"
                f"review: {r.get('review_reason_short','')}")
        ax.text(0.0, 0.40, info, transform=ax.transAxes, fontsize=6.8, va="top", color="#333")
        for s in ("top", "right", "left", "bottom"):
            ax.spines[s].set_visible(False)
    # hide unused axes
    for k in range(nrev, nr * ncols):
        axes[k // ncols][k % ncols].axis("off")
    fig.suptitle("Supplement 2: review-case panels (retained and interpreted, not hidden)\n"
                 "InterProScan pending", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    paths = _save(fig, figdir, "Supplement_Figure_2_review_cases_pre_interpro")
    return {"figure_id": "Supplement_Figure_2_review_cases", "paths": paths,
            "source_tables": table.name, "n_species": nrev, "n_rows": len(rows),
            "main_message": "Review cases are retained and interpreted, not hidden."}


# ---------------------------------------------------------------------------
# SUPPLEMENT 3 — InterProScan input readiness
# ---------------------------------------------------------------------------
def supplement3_interpro_readiness(data: FinalData, figdir: Path, tabledir: Path) -> Dict[str, object]:
    m = data.interpro
    items = [
        ("selected proteins", _to_int(m.get("total_selected_proteins"), 0)),
        ("unique proteins (after collapse)", _to_int(m.get("unique_sequences"), 0)),
        ("duplicate sequences collapsed", _to_int(m.get("duplicates_collapsed"), 0)),
        ("invalid sequences rejected", _to_int(m.get("invalid_sequences_rejected"), 0)),
        ("species covered", _to_int(m.get("species_covered"), 0)),
        ("species with both isoforms", _to_int(m.get("species_with_both_isoforms"), 0)),
    ]
    rows = [{"metric": k, "value": v} for k, v in items]
    rows.append({"metric": "interpro_status", "value": m.get("interpro_status", "")})
    rows.append({"metric": "interpro_input_fasta", "value": m.get("interpro_input_fasta", "")})
    table = tabledir / "supplement_interproscan_input_readiness.tsv"
    write_tsv(table, rows, ["metric", "value"])

    fig, ax = plt.subplots(figsize=(9, 5.5))
    labels = [k for k, _ in items]
    vals = [v for _, v in items]
    ypos = list(range(len(labels)))[::-1]
    bars = ax.barh(ypos, vals, color=C_IIIB, edgecolor="#003E63", height=0.6)
    bars[3].set_color("#CC6677")  # invalid rejected
    for y, v in zip(ypos, vals):
        ax.text(v + max(vals) * 0.01, y, str(v), va="center", fontsize=9, fontweight="bold")
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("count", fontsize=10)
    ax.set_title("Supplement 3: InterProScan input readiness", fontsize=12.5, fontweight="bold")
    ax.text(0.0, -0.16, "InterProScan not yet executed; this figure summarizes validated input readiness.",
            transform=ax.transAxes, fontsize=9, color="#7A5A00", fontweight="bold")
    ax.grid(True, axis="x", color="#EEEEEE", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    paths = _save(fig, figdir, "Supplement_Figure_3_interproscan_input_readiness_pre_interpro")
    return {"figure_id": "Supplement_Figure_3_interproscan_input_readiness", "paths": paths,
            "source_tables": table.name, "n_species": _to_int(m.get("species_covered"), 0),
            "n_rows": len(rows),
            "main_message": "Validated, non-redundant InterProScan input is prepared; "
                            "InterProScan not yet executed."}


# ---------------------------------------------------------------------------
# Orchestration: render all figures + captions + manifest
# ---------------------------------------------------------------------------
def _locate(base: Path, name: str, hint: str = "") -> Path:
    matches = sorted(base.rglob(name))
    if not matches:
        raise RuntimeError(f"Required input not found under {base}: {name}")
    if hint:
        for m in matches:
            if hint in str(m):
                return m
    return sorted(matches, key=lambda p: len(p.parts))[0]


FIGURE_FUNCS = [
    fig1_framework, fig2_architecture, fig3_cassette_zoom, fig4_evidence_matrix,
    fig5_native_vs_normalized, supplement1_all_native, supplement2_review_cases,
    supplement3_interpro_readiness,
]

CAPTIONS = {
    "Figure_1_framework": (
        "Framework overview of the FGFR2 IIIb/IIIc analysis up to InterProScan input preparation. "
        "The figure shows the annotation-aware workflow from species selection through orthology/paralog "
        "validation, sequence-calibrated IIIb/IIIc direction assignment, protein QC, the exon-to-protein "
        "coordinate resolver, and main/review classification, ending in a non-redundant InterProScan-ready "
        "FASTA. Count badges summarize species, resolved IIIb/IIIc mappings, main vs supplement/review "
        "species, and unique input sequences. "
        "Supported claim: the pipeline is a reproducible, multi-evidence framework. "
        "Not claimed: no InterPro/Pfam domain annotation is shown; InterProScan is pending."),
    "Figure_2_exon_to_protein_architecture": (
        "Exon/CDS-to-protein architecture of resolved IIIb and IIIc cassettes across 30 vertebrate FGFR2 "
        "orthologs. Each species has an IIIb track and an IIIc track. Within each track the real CDS/exon "
        "blocks of the resolved transcript are drawn as separate grey blocks (one block per coding exon, "
        "projected into protein-coordinate space from cds_features.tsv), with thin white separators between "
        "blocks; the resolved IIIb cassette is highlighted in blue and IIIc in orange. The cassette block is "
        "identified by genomic/protein coordinate overlap against a table of unique CDS blocks "
        "(fgfr2_cassette_cds_block_map.tsv), not by a CDS identifier; this corrects an earlier mapping bug in "
        "which non-unique NCBI cds-XP_ ids placed cassettes on the first CDS block (protein_start_aa = 1). "
        "Boundary precision is shown only as small, calm edge symbols (a grey '/' for a known split-codon "
        "boundary, a pale dot for phase-unavailable-but-coordinate-resolved); these minor technical flags "
        "are intentionally not visually alarming because they are not errors. A prominent warning glyph "
        "(\u26A0) and a dashed cassette outline appear ONLY for true issues (protein-overlay-only / no CDS "
        "model, protein conflict, or major native-coordinate offset), and a red \u2715 marks a true hard "
        "coordinate sanity failure; these are driven by plot_visibility_level in "
        "fgfr2_refined_uncertainty_classes.tsv. Where "
        "full per-exon CDS blocks are not available for a transcript, the body is drawn as a lighter hatched "
        "bar labelled 'CDS architecture partially reconstructed' and the cassette is overlaid; missing exon-"
        "block completeness is reported in figure2_exon_block_completeness.tsv. "
        "Supported claim: the IIIb/IIIc cassettes sit within a real, exon-resolved protein architecture that "
        "is consistent across vertebrates. Not claimed: no InterPro/Pfam domain annotation is shown; "
        "InterProScan is pending."),
    "Figure_3_IgIII_cassette_zoom": (
        "IgIII/D3 IIIb/IIIc event zoom. For each species the IIIb and IIIc tracks are shown in a normalized "
        "IgIII-window coordinate (window start = 0) and decomposed into three conceptual segments derived "
        "from real coordinates: an upstream conserved region (grey), the alternative cassette (blue for IIIb, "
        "orange for IIIc) and a downstream conserved region (grey). Cassette edges carry CDS-boundary "
        "precision symbols, a small protein-QC icon is shown per track, and review species use a dashed "
        "cassette outline plus an amber marker. "
        "Supported claim: the upstream\u2013cassette\u2013downstream IgIII event structure and the normalized "
        "cassette-internal localization are stable across vertebrate FGFR2 orthologs even where native "
        "coordinates shift. Not claimed: no domain-level InterPro annotation; InterProScan is pending."),
    "Figure_4_species_evidence_matrix": (
        "Species evidence matrix integrating independent QC layers (orthology, multi-vertebrate paralog "
        "panel, isoform detection, sequence-calibrated direction, protein QC, resolver status, CDS-boundary "
        "precision, native/normalized coordinate sanity, III-region similarity, final display class, and "
        "InterProScan input status). The rightmost column shows the refined uncertainty class with calm, "
        "color-blind-safe states: split-codon boundaries (light blue) and phase-unavailable-but-coordinate-"
        "resolved cases (pale grey-blue) are deliberately NOT alarming because they are known boundary "
        "behaviour, not errors; transcript/protein-level-only cases are muted purple/grey; review cases are "
        "amber and are interpreted separately, not used for primary claims; only a true hard coordinate "
        "sanity failure is red; InterPro status is neutral grey because InterProScan is pending. "
        "Supported claim: multiple evidence layers separate primary-analysis from review/supplement species "
        "and minor boundary-precision flags are distinguished from genuine review/hard-fail cases. "
        "Not claimed: InterPro domain calls are pending."),
    "Figure_5_native_vs_normalized_coordinate_qc": (
        "Native vs normalized coordinate QC. Each point is a species; native IIIb\u2013IIIc cassette offset "
        "is plotted against the normalized III-slot offset. Supported claim: native offsets flag annotation/"
        "isoform-dependent cases while normalized coordinates remain stable. Not claimed: no InterPro domains; "
        "InterProScan pending. Only supplement/review species are labelled."),
    "Supplement_Figure_1_all_species_native_tracks": (
        "All 30 species shown with their real exon/CDS-derived protein architecture and IIIb/IIIc cassette "
        "positions, using the same colors, phylogenetic order, taxon-group bands and review markers as "
        "Figure 2. InterProScan pending."),
    "Supplement_Figure_2_review_cases": (
        "Compact per-species panels for supplement/review species read directly from species_qc_master.tsv. "
        "Each panel shows IIIb/IIIc cassette positions, native/normalized sanity, protein QC, III-region "
        "similarity, CDS-boundary precision and the short review reason. Review cases are retained and "
        "interpreted, not hidden, and are not used for primary claims. InterProScan pending."),
    "Supplement_Figure_3_interproscan_input_readiness": (
        "InterProScan input readiness summary: selected proteins, unique proteins after duplicate collapse, "
        "duplicates collapsed, invalid sequences, and species/isoform coverage. Explicit note: InterProScan "
        "has not yet been executed; the figure summarizes validated input readiness only."),
}

CAPTION_COMMON = (
    "IIIb is shown in blue (#0072B2) and IIIc in orange (#E69F00); these colors are identical across all "
    "figures. Supplement/review species are marked with an amber diamond and dashed cassette outlines and "
    "remain visible. Taxon groups (Primates, Other mammals, Birds, Reptiles, Amphibians, Teleost fish) are "
    "indicated by subtle left-side color bars / light background bands and separator lines; species are "
    "ordered by a reproducible phylogenetic/taxonomic order (species_phylogenetic_order.tsv), not "
    "alphabetically. Domain-aware overlays are a downstream step once InterProScan annotations are available.")


def write_captions(captiondir: Path, manifest: List[Dict[str, object]]) -> Path:
    captiondir.mkdir(parents=True, exist_ok=True)
    out = captiondir / "figure_captions_pre_interpro.md"
    lines = ["# Pre-InterPro figure captions", "",
             "_All figures are pre-InterPro. No real InterPro/Pfam domain annotation is shown; "
             "InterProScan has not yet been executed._", "",
             "## Cassette → CDS-block mapping (corrected)", "",
             "IIIb/IIIc cassettes are mapped onto **unique CDS blocks** by genomic/protein coordinate "
             "overlap (`fgfr2_cassette_cds_block_map.tsv`), not by CDS identifier. This corrects an "
             "earlier bug where non-unique NCBI `cds-XP_` identifiers placed cassettes on the first CDS "
             "block (`protein_start_aa = 1`). All final figures passed the cassette-coordinate sanity "
             "gate: no main-analysis or control-primate IIIb/IIIc cassette is N-terminal in a full-length "
             "FGFR2 protein, and figure-2 cassette positions match the resolved coordinate audit.", "",
             "## CDS-boundary uncertainty (how to read the boundary symbols)", "",
             "The exon-to-protein coordinate resolver records, for every resolved IIIb/IIIc cassette, "
             "whether each cassette boundary falls exactly on a codon boundary, splits a codon, or cannot "
             "be phase-resolved. These states are auditable in `cds_phase_boundary_audit.tsv` and "
             "summarised in `cds_phase_boundary_explainability_summary.tsv`.", "",
             "- **Why unknown codon phase can remain even when CDS coordinates are known:** a CDS feature "
             "has genomic start/end coordinates, but the *reading-frame phase* at the boundary is a separate "
             "annotation. In this dataset all `unknown_codon_phase` cases are explained: the source GFF3 did "
             "not propagate a phase value (`phase_not_propagated_from_source`, all Ensembl-sourced cassettes) "
             "or the cassette transcript was not present in the local CDS model "
             "(`nucleotide_sequence_unavailable`). Coordinates are known; the codon *phase* is simply not "
             "annotated, so the boundary is reported as uncertain rather than guessed.", "",
             "- **Split-codon boundaries** (`codon_split_one_side` / `codon_split_both_sides`) are expected "
             "biology for internal cassette exons whose length is not a multiple of three; they are flagged, "
             "not corrected.", "",
             "- **Uncertain boundaries are evidence-level flags, not errors.** Uncertain cases are never "
             "forced to `exact`. As an independent control, 58/60 resolved cassette transcripts reconstruct "
             "from CDS coordinates to a total length consistent with the selected protein "
             "(`protein_translation_check_status = cds_protein_length_consistent`), confirming the protein "
             "coordinate projection even where codon phase is unannotated.", "",
             "## Uncertainty classes and why they are shown differently", "",
             "The refined uncertainty classes (`fgfr2_refined_uncertainty_classes.tsv`, summarised per "
             "species in `species_qc_master.tsv`) collapse the many low-level flags into a small, "
             "explainable set with an explicit `plot_visibility_level`, so figures stop overstating "
             "uncertainty:", "",
             "- **A known split-codon boundary is not a failure.** Internal cassette exons whose length is "
             "not a multiple of three necessarily split a codon at a boundary; this is expected biology and "
             "is shown only as a small grey edge symbol.", "",
             "- **Phase unavailable does not mean the coordinate is wrong.** Where the source GFF3 did not "
             "propagate a codon phase, the cassette is still coordinate-resolved; for 20/60 such cases the "
             "phase was rescued from a length-consistent CDS reconstruction (the cumulative-CDS reading "
             "frame reproduces the NCBI source-phase split/exact calls exactly), and the remainder are "
             "labelled `phase_not_available_but_coordinate_resolved` rather than wrong.", "",
             "- **True missing data are rare and explicitly counted.** Only 2/60 mappings "
             "(`protein_overlay_no_cds_model` / `nucleotide_sequence_unavailable`) lack a local CDS-block "
             "model; these are the only NCBI-patch candidates.", "",
             "- **NCBI patching is used only for true missing information and is provenance-tracked** "
             "(`fgfr2_ncbi_cds_boundary_patch_report.tsv`); known split/phase flags and locally "
             "reconstructable cases are explicitly NOT patched, and Ensembl/NCBI releases are never mixed "
             "silently.", "",
             "- **Minor boundary-precision flags are shown subtly** (small edge symbols / pale colors) to "
             "avoid implying that coordinate-resolved species are problematic.", "",
             "- **Review cases remain visible and are not hidden.** Protein-conflict, major native-"
             "coordinate-offset and hard coordinate sanity failures keep a prominent marker and are "
             "interpreted separately; they are not used for primary claims.", "",
             "## Figures and their main biological messages", ""]
    for fid, msg in [
        ("Figure 1", "Annotation-aware exon-to-protein boundary framework up to InterProScan preparation."),
        ("Figure 2", "Real exon/CDS-derived protein architecture; IIIb/IIIc cassettes sit in a conserved, exon-resolved context."),
        ("Figure 3", "Upstream\u2013cassette\u2013downstream IgIII event structure; normalized cassette mapping is stable across vertebrates."),
        ("Figure 4", "Multiple independent evidence layers separate main-analysis from review/supplement species."),
        ("Figure 5", "Native offsets flag annotation-dependent cases; normalized III-slot coordinates stay stable for most species."),
        ("Supplements 1\u20133", "All-species architecture, review-case detail, and InterProScan input readiness."),
    ]:
        lines.append(f"- **{fid}:** {msg}")
    lines += ["", "_Pre-InterPro scope: figures show exon-to-protein architecture, resolved IIIb/IIIc "
              "cassette coordinates and QC/evidence layers. They do not show and do not claim real "
              "InterPro/Pfam domain annotations._", ""]
    titles = {
        "Figure_1_framework": "Figure 1. Framework overview up to InterProScan preparation",
        "Figure_2_exon_to_protein_architecture": "Figure 2. Main exon-to-protein architecture map",
        "Figure_3_IgIII_cassette_zoom": "Figure 3. IgIII/D3 IIIb/IIIc cassette zoom",
        "Figure_4_species_evidence_matrix": "Figure 4. Species evidence matrix / QC heatmap",
        "Figure_5_native_vs_normalized_coordinate_qc": "Figure 5. Native vs normalized coordinate QC",
        "Supplement_Figure_1_all_species_native_tracks": "Supplement Figure 1. All-species native tracks",
        "Supplement_Figure_2_review_cases": "Supplement Figure 2. Review-case panels",
        "Supplement_Figure_3_interproscan_input_readiness": "Supplement Figure 3. InterProScan input readiness",
    }
    for entry in manifest:
        fid = entry["figure_id"]
        lines.append(f"## {titles.get(fid, fid)}")
        lines.append("")
        lines.append(CAPTIONS.get(fid, ""))
        lines.append("")
        lines.append(CAPTION_COMMON)
        lines.append("")
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def write_publication_manifest(metadir: Path, manifest: List[Dict[str, object]]) -> Path:
    metadir.mkdir(parents=True, exist_ok=True)
    out = metadir / "publication_figure_manifest.tsv"
    fields = ["figure_id", "status", "skip_reason", "file_svg", "file_pdf", "file_png",
              "source_tables", "n_species", "n_rows", "uses_phylogenetic_order",
              "uses_taxon_group_bands", "canonical_qc_source", "claims_interpro_domains",
              "main_message", "interpro_status", "created_at"]
    bands = {"Figure_2_exon_to_protein_architecture", "Figure_3_IgIII_cassette_zoom",
             "Figure_4_species_evidence_matrix", "Supplement_Figure_1_all_species_native_tracks"}
    rows = []
    for e in manifest:
        rows.append({
            "figure_id": e["figure_id"],
            "status": e.get("status", "rendered"),
            "skip_reason": e.get("skip_reason", ""),
            "file_svg": e["paths"].get("svg", ""),
            "file_pdf": e["paths"].get("pdf", ""),
            "file_png": e["paths"].get("png", ""),
            "source_tables": e.get("source_tables", ""),
            "n_species": e.get("n_species", ""),
            "n_rows": e.get("n_rows", ""),
            "uses_phylogenetic_order": "true",
            "uses_taxon_group_bands": "true" if e["figure_id"] in bands else "false",
            "canonical_qc_source": "species_qc_master.tsv",
            "claims_interpro_domains": "false",
            "main_message": e.get("main_message", ""),
            "interpro_status": "interpro_pending_or_input_prepared",
            "created_at": _now(),
        })
    write_tsv(out, rows, fields)
    return out


CONTROL_SPECIES = ("homo_sapiens", "pan_troglodytes", "macaca_mulatta")


def cassette_sanity_gate(base: Path, data: "FinalData", metadir: Path) -> None:
    """PART G — refuse to generate publication figures if biologically implausible
    cassette positions remain for any MAIN-ANALYSIS species, if a control primate
    cassette is N-terminal, or if figure-2 cassette rows disagree with the coordinate
    audit. On failure, write publication_figure_validation_failed.tsv and raise."""
    failures: List[Dict[str, object]] = []
    for r in data.master:
        sp = r["species"]
        main = str(r.get("final_display_class", "")).startswith("main_analysis")
        for iso in ("IIIb", "IIIc"):
            cm = data.cassette(sp, iso)
            c = data.coord(sp, iso) or {}
            start = _to_int(cm.get("matched_protein_start_aa"))
            rank = _to_int(cm.get("matched_cds_rank"))
            plen = _to_int(c.get("protein_length_aa"))
            full = plen is not None and plen > 500
            # figure2 vs coordinate-audit agreement (matched interval vs native interval)
            ns = _to_int(c.get("native_protein_start_aa"))
            if main and start is not None and ns is not None and abs(start - ns) > 5:
                failures.append({"species": sp, "isoform": iso, "reason": "figure_vs_audit_mismatch",
                                 "detail": f"matched_start={start} vs native_start={ns}"})
            if main and full and start is not None and start < 150:
                failures.append({"species": sp, "isoform": iso, "reason": "cassette_start_lt_150_full_length",
                                 "detail": f"start={start}, protein_length={plen}"})
            if main and full and rank == 1:
                failures.append({"species": sp, "isoform": iso, "reason": "cassette_cds_rank_1_full_length",
                                 "detail": f"protein_length={plen}"})
            if sp in CONTROL_SPECIES and start is not None and start < 150:
                failures.append({"species": sp, "isoform": iso, "reason": "control_primate_cassette_n_terminal",
                                 "detail": f"start={start} (<150) for control species"})
    if failures:
        metadir.mkdir(parents=True, exist_ok=True)
        write_tsv(metadir / "publication_figure_validation_failed.tsv", failures,
                  ["species", "isoform", "reason", "detail"])
        msg = "; ".join(f"{f['species']}/{f['isoform']}:{f['reason']}" for f in failures[:8])
        raise RuntimeError(
            f"Cassette-coordinate sanity gate FAILED for {len(failures)} main-analysis/control "
            f"cassette(s); publication figures NOT generated. See "
            f"publication_figure_validation_failed.tsv. Examples: {msg}")


def render_all(base: Path) -> Tuple[List[Dict[str, object]], Dict[str, Path]]:
    pub = base / "11_publication_figures_pre_interpro"
    figdir, tabledir = pub / "figures", pub / "tables"
    captiondir, metadir = pub / "captions", pub / "metadata"
    for d in (figdir, tabledir, captiondir, metadir):
        d.mkdir(parents=True, exist_ok=True)

    master = _locate(base, "species_qc_master.tsv", "11_pre_interpro_master")
    coord_audit = _locate(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    interpro_summary = _locate(base, "fgfr2_interpro_prepare_summary.tsv")
    cds_features = None
    try:
        cds_features = _locate(base, "cds_features.tsv", "02_models")
    except RuntimeError:
        cds_features = None
    cassette_map = None
    try:
        cassette_map = _locate(base, "fgfr2_cassette_cds_block_map.tsv")
    except RuntimeError:
        cassette_map = None
    refined = None
    try:
        refined = _locate(base, "fgfr2_refined_uncertainty_classes.tsv")
    except RuntimeError:
        refined = None
    data = FinalData(master, coord_audit, interpro_summary, cds_features, cassette_map, refined)

    # PART G — fail-fast biological sanity gate (no figures if implausible cassettes remain)
    cassette_sanity_gate(base, data, metadir)

    manifest: List[Dict[str, object]] = []
    for fn in FIGURE_FUNCS:
        manifest.append(fn(data, figdir, tabledir))

    cap = write_captions(captiondir, manifest)
    man = write_publication_manifest(metadir, manifest)
    return manifest, {"captions": cap, "manifest": man, "figures": figdir,
                      "tables": tabledir, "metadata": metadir}


def main() -> int:
    ap = argparse.ArgumentParser(description="Publication-level pre-InterPro figures (Parts 3,5,6).")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    manifest, paths = render_all(args.base)
    print(f"[OK] rendered {len(manifest)} publication figures (SVG/PDF/PNG)")
    for e in manifest:
        print(f"     {e['figure_id']}: {Path(e['paths']['png']).name}")
    print(f"     captions -> {paths['captions']}")
    print(f"     manifest -> {paths['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
