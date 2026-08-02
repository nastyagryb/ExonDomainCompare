#!/usr/bin/env python3
"""
fgfr2_plot_style.py — shared paper-level plotting style for FGFR2 synteny figures (and others).

Colour-blind-safe (Okabe-Ito based), no rainbow, no red/green-only coding, compact legends, clean
vector export (SVG/PDF/PNG). Deterministic per-orthology-group colours so the same gene keeps the
same colour across every figure. Helpers for rounded gene arrows, evidence badges, panel labels and
taxon ordering.
"""

from __future__ import annotations

import hashlib
import os
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl_fgfr2_"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Patch, Polygon
from matplotlib.lines import Line2D  # noqa: E402

# ---------------------------------------------------------------------------
# global rcParams
# ---------------------------------------------------------------------------
FONT = {"title": 12, "subtitle": 9, "label": 8.5, "gene": 7.2, "small": 6.4, "tick": 7.5,
        "legend": 7.5, "badge": 6.6}
LW = {"thin": 0.5, "gene_edge": 0.6, "ribbon": 1.1, "outline": 1.4, "rule": 0.8}


def apply_rcparams() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": LW["thin"], "figure.dpi": 150, "savefig.dpi": 220,
        "axes.edgecolor": "#444444", "text.color": "#1A1A1A",
        "axes.labelcolor": "#1A1A1A", "xtick.color": "#444444", "ytick.color": "#444444",
    })


# ---------------------------------------------------------------------------
# palette
# ---------------------------------------------------------------------------
# Okabe-Ito + a few extra colour-blind-distinguishable hues (no pure red/green pairing)
PALETTE: List[str] = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442",
    "#117733", "#882255", "#44AA99", "#332288", "#AA4499", "#88CCEE", "#999933",
    "#661100", "#6699CC", "#DDCC77", "#774411", "#117777", "#771155",
]
FGFR2_COLOR = "#1A1A1A"          # anchor
UNRESOLVED_COLOR = "#C9CDD2"     # unresolved LOC / raw
NONREF_COLOR = "#9AA0A6"         # resolved non-reference local gene (neutral, distinct from unresolved)
INK = "#1A1A1A"
MUTED = "#6B7280"

# status / evidence colours (sequential-ish blues→amber, no green/red collision)
STATUS_COLOR = {
    "confirmed_same_order": "#1B6CA8",
    "confirmed_reordered": "#5B9BD5",
    "probable_rbh": "#7B6FB0",
    "probable_blast": "#A9A0D6",
    "non_reference_local_gene": NONREF_COLOR,
    "ambiguous_paralog_family": "#E69F00",
    "unresolved_missing": UNRESOLVED_COLOR,
    "scaffold_unavailable": "#7A7F87",
}
SYNTENY_CLASS_COLOR = {
    "synteny_strong": "#1B6CA8",
    "synteny_supported_with_minor_rearrangement": "#5B9BD5",
    "synteny_partial_blast_supported": "#7B6FB0",
    "synteny_partial_scaffold_limit": "#44AA99",
    "synteny_sequence_only_support": "#E69F00",
    "synteny_conflict_review": "#D55E00",
    "synteny_unavailable": "#B8BCC2",
}
CLAIM_COLOR = {"primary": "#1B6CA8", "supplement": "#E69F00", "excluded": "#D55E00", "": "#C9CDD2"}
BADGE_COLOR = {"ok": "#1B6CA8", "partial": "#E69F00", "fail": "#D55E00", "neutral": "#9AA0A6"}

# stable taxon ordering for grouping
TAXON_ORDER = ["Primates", "Rodentia", "Carnivora", "Cetartiodactyla", "Artiodactyla",
               "Mammalia", "Marsupialia", "Monotremata", "Aves", "Crocodylia", "Testudines",
               "Reptilia", "Squamata", "Amphibia", "Actinopterygii", "Teleostei", "Fish"]


def ortholog_color(group: str, assigned: Dict[str, str]) -> str:
    """Deterministic colour for an orthology group symbol; stable across figures via a hash fallback."""
    if not group:
        return UNRESOLVED_COLOR
    g = group.upper()
    if g == "FGFR2":
        return FGFR2_COLOR
    if g in assigned:
        return assigned[g]
    # assign next free palette colour deterministically by insertion, with hash tiebreak
    idx = len(assigned)
    if idx < len(PALETTE):
        col = PALETTE[idx]
    else:
        h = int(hashlib.md5(g.encode()).hexdigest(), 16)
        col = PALETTE[h % len(PALETTE)]
    assigned[g] = col
    return col


def build_color_map(groups: Sequence[str]) -> Dict[str, str]:
    """Assign palette colours to orthology groups in a stable, frequency-then-alpha order."""
    assigned: Dict[str, str] = {"FGFR2": FGFR2_COLOR}
    for g in sorted({(x or "").upper() for x in groups if x and x.upper() != "FGFR2"}):
        ortholog_color(g, assigned)
    return assigned


def taxon_sort_key(taxon: str, phylo_order) -> Tuple[int, int]:
    try:
        po = int(phylo_order)
    except (TypeError, ValueError):
        po = 999
    t = taxon or ""
    ti = TAXON_ORDER.index(t) if t in TAXON_ORDER else len(TAXON_ORDER)
    return (ti, po)


# ---------------------------------------------------------------------------
# drawing helpers
# ---------------------------------------------------------------------------
def gene_arrow(ax, x_center, y, width, height, strand_rel, facecolor, *, edgecolor=INK,
               lw=LW["gene_edge"], alpha=1.0, zorder=3):
    """Rounded gene arrow (pentagon) centred at x_center, pointing right if strand_rel>=0."""
    w, h = width / 2.0, height / 2.0
    tip = 0.32 * width
    if strand_rel >= 0:
        pts = [(x_center - w, y - h), (x_center + w - tip, y - h), (x_center + w, y),
               (x_center + w - tip, y + h), (x_center - w, y + h)]
    else:
        pts = [(x_center + w, y - h), (x_center - w + tip, y - h), (x_center - w, y),
               (x_center - w + tip, y + h), (x_center + w, y + h)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=facecolor, edgecolor=edgecolor,
                         linewidth=lw, alpha=alpha, zorder=zorder, joinstyle="round"))


def ribbon(ax, x0, y0, x1, y1, color, *, alpha=0.30, lw=LW["ribbon"], zorder=1):
    """Thin connector between conserved neighbors in consecutive rows."""
    ax.plot([x0, x1], [y0, y1], color=color, alpha=alpha, lw=lw, zorder=zorder,
            solid_capstyle="round")


def badge(ax, x, y, label, kind="neutral", *, fontsize=FONT["badge"]):
    """Small rounded evidence badge."""
    col = BADGE_COLOR.get(kind, BADGE_COLOR["neutral"])
    ax.add_patch(FancyBboxPatch((x, y - 0.12), 0.02, 0.24, boxstyle="round,pad=0.16,rounding_size=0.12",
                                facecolor=col, edgecolor="none", alpha=0.16, zorder=4,
                                mutation_aspect=1.0))
    ax.text(x + 0.18, y, label, fontsize=fontsize, va="center", ha="left", color=col,
            fontweight="bold", zorder=5)


def panel_label(ax, letter, *, x=-0.02, y=1.02):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=13, fontweight="bold",
            va="bottom", ha="right", color=INK)


def title(ax, text, subtitle=None):
    ax.set_title(text, fontsize=FONT["title"], fontweight="bold", loc="left", color=INK, pad=10)
    if subtitle:
        ax.text(0.0, 1.005, subtitle, transform=ax.transAxes, fontsize=FONT["subtitle"],
                color=MUTED, va="bottom", ha="left")


def compact_legend(ax, handles, *, loc="lower center", ncol=4, bbox=(0.5, -0.12)):
    leg = ax.legend(handles=handles, loc=loc, bbox_to_anchor=bbox, ncol=ncol,
                    fontsize=FONT["legend"], frameon=False, handlelength=1.1,
                    columnspacing=1.3, handletextpad=0.5)
    return leg


def legend_patch(color, label):
    return Patch(facecolor=color, edgecolor="none", label=label)


def legend_line(color, label, ls="-", lw=LW["outline"]):
    return Line2D([0], [0], color=color, ls=ls, lw=lw, label=label)


def gene_label_style(method: str) -> Tuple[str, str]:
    """Return (suffix, fontstyle) for a neighbor label given its resolution method."""
    if method in ("reciprocal_best_hit", "rbh_supported_neighbor_ortholog"):
        return "", "italic"
    if method in ("high_confidence_one_way_blast", "blast_supported_probable_neighbor"):
        return "?", "italic"
    if method in ("unresolved", "raw_annotation_only", "unresolved_LOC", "raw_id_only_no_sequence"):
        return "", "normal"
    return "", "normal"


def savefig(fig, fig_dir: Path, stem: str) -> None:
    for ext in ("svg", "pdf", "png"):
        fig.savefig(fig_dir / f"{stem}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)
