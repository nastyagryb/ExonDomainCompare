
"""Shared modern publication style for the FGFR2 analysis suite."""
from __future__ import annotations

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

COLORS = {
    "navy": "#102A43",
    "blue": "#2F6B9A",
    "sky": "#64B5D2",
    "teal": "#138A80",
    "green": "#2E8B57",
    "gold": "#D59B22",
    "orange": "#D96B38",
    "red": "#B94A48",
    "purple": "#7561A8",
    "gray": "#6B7785",
    "light": "#EDF2F7",
    "grid": "#DCE4EC",
    "dark": "#243B53",
    "white": "#FFFFFF",
}

def apply_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.titlesize": 18,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10.5,
        "figure.titlesize": 22,
        "axes.linewidth": 0.9,
        "axes.edgecolor": COLORS["dark"],
        "axes.labelcolor": COLORS["dark"],
        "xtick.color": COLORS["dark"],
        "ytick.color": COLORS["dark"],
        "text.color": COLORS["dark"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })

def clean_axis(ax, grid_axis: str | None = "x") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=COLORS["grid"], linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)

def title_block(fig, title: str, subtitle: str) -> None:
    fig.suptitle(title, x=0.06, y=0.985, ha="left", va="top",
                 fontsize=22, fontweight="bold", color=COLORS["navy"])
    fig.text(0.06, 0.925, subtitle, ha="left", va="top",
             fontsize=12.5, color=COLORS["gray"])

def takehome(ax, text: str, loc=(0.98, 0.04), align="right",
             face="#F3F7FA") -> None:
    ax.text(loc[0], loc[1], text, transform=ax.transAxes,
            ha=align, va="bottom", fontsize=11.5, fontweight="bold",
            color=COLORS["navy"],
            bbox=dict(boxstyle="round,pad=0.55", facecolor=face,
                      edgecolor=COLORS["grid"], linewidth=0.8))

def save_figure(fig, outdir: str | Path, stem: str, dpi: int = 400) -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{stem}.png", dpi=dpi)
    fig.savefig(out / f"{stem}.svg")
    fig.savefig(out / f"{stem}.pdf")
    plt.close(fig)
