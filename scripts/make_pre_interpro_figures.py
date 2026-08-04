#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl_preinterpro_"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# Okabe-Ito colour-blind-safe palette.
C_IIIB = "#0072B2"   # blue
C_IIIC = "#E69F00"   # orange
C_OK = "#009E73"     # green
C_REVIEW = "#D55E00" # vermillion
C_NEUTRAL = "#999999"
C_PENDING = "#CCCCCC"
TAXON_COLORS = {
    "mammal": "#0072B2", "bird": "#E69F00", "reptile": "#009E73",
    "amphibian": "#CC79A7", "fish": "#56B4E9", "unknown": "#999999",
}

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9,
    "figure.dpi": 150, "savefig.bbox": "tight", "axes.spines.top": False,
    "axes.spines.right": False,
})


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def to_float(v: str) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def save(fig, outdir: Path, stem: str) -> None:
    for ext in ("svg", "pdf", "png"):
        fig.savefig(outdir / f"{stem}.{ext}")
    plt.close(fig)


def fig1_framework(rows: List[Dict[str, str]], outdir: Path) -> None:
    cats = ["final_display_class", "taxon_group", "cds_boundary_precision_summary", "interpro_status"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, cat in zip(axes.flat, cats):
        sub = [r for r in rows if r["category"] == cat]
        sub.sort(key=lambda r: int(r["count"]), reverse=True)
        labels = [r["level"] for r in sub]
        vals = [int(r["count"]) for r in sub]
        colors = [C_REVIEW if ("review" in l or "supplementary" in l or l in ("false", "unknown_codon_phase",
                  "codon_split_both_sides")) else C_OK for l in labels]
        ax.barh(range(len(labels)), vals, color=colors, edgecolor="black", linewidth=0.4)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels([l[:34] for l in labels], fontsize=7)
        ax.invert_yaxis()
        ax.set_title(cat)
        ax.set_xlabel("species count")
        for i, v in enumerate(vals):
            ax.text(v + 0.1, i, str(v), va="center", fontsize=7)
    fig.suptitle("FGFR2 IIIb/IIIc framework counts (pre-InterProScan)", fontweight="bold")
    legend = [Patch(facecolor=C_OK, edgecolor="black", label="high-confidence / clean"),
              Patch(facecolor=C_REVIEW, edgecolor="black", label="review / supplementary")]
    fig.legend(handles=legend, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    save(fig, outdir, "Figure_1_framework_pre_interpro")


def fig2_exon_to_protein(rows: List[Dict[str, str]], outdir: Path) -> None:
    by_sp: Dict[str, Dict[str, Dict[str, str]]] = {}
    review_sp = set()
    for r in rows:
        by_sp.setdefault(r["species"], {})[r["isoform"]] = r
        if r.get("final_display_class") != "main_analysis_high_confidence":
            review_sp.add(r["species"])
    species = sorted(by_sp, key=lambda s: (s not in review_sp, s))
    fig, ax = plt.subplots(figsize=(11, max(5, 0.34 * len(species) + 1)))
    yticks, ylabels = [], []
    for i, sp in enumerate(species):
        y = len(species) - i
        yticks.append(y)
        disp = by_sp[sp].get("IIIb", by_sp[sp].get("IIIc", {})).get("display_species_name", sp)
        ylabels.append(("* " if sp in review_sp else "") + disp)
        for iso, color in (("IIIb", C_IIIB), ("IIIc", C_IIIC)):
            r = by_sp[sp].get(iso, {})
            s = to_float(r.get("native_start_aa", "")); e = to_float(r.get("native_end_aa", ""))
            if s is None or e is None:
                continue
            off = 0.18 if iso == "IIIb" else -0.18
            split = r.get("cds_boundary_precision_refined", "")
            hatch = "///" if split in ("codon_split_both_sides", "unknown_codon_phase") else None
            ax.barh(y + off, e - s, left=s, height=0.32, color=color, edgecolor="black",
                    linewidth=0.5, hatch=hatch)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("native protein coordinate (aa)")
    ax.set_title("FGFR2 IIIb/IIIc cassette position on the protein axis (pre-InterProScan)\n"
                 "Domain overlays: InterProScan pending (not shown)", fontweight="bold")
    legend = [Patch(facecolor=C_IIIB, edgecolor="black", label="IIIb cassette"),
              Patch(facecolor=C_IIIC, edgecolor="black", label="IIIc cassette"),
              Patch(facecolor="white", edgecolor="black", hatch="///", label="split/unknown codon boundary"),
              Patch(facecolor=C_PENDING, edgecolor="black", label="InterProScan domains pending"),
              Patch(facecolor="none", edgecolor="none", label="* = review/supplementary species")]
    ax.legend(handles=legend, loc="upper right", fontsize=7, frameon=True)
    fig.tight_layout()
    save(fig, outdir, "Figure_2_exon_to_protein_map_pre_interpro")


def fig3_evidence_heatmap(rows: List[Dict[str, str]], outdir: Path) -> None:
    dims = sorted({r["evidence_dimension"] for r in rows})
    species = sorted({r["species"] for r in rows})
    review_sp = {r["species"] for r in rows if r.get("final_display_class") != "main_analysis_high_confidence"}
    species = sorted(species, key=lambda s: (s not in review_sp, s))
    disp = {r["species"]: r.get("display_species_name", r["species"]) for r in rows}
    idx = {(r["species"], r["evidence_dimension"]): int(r["ok_flag"]) for r in rows}
    import numpy as np
    mat = np.array([[idx.get((sp, d), 0) for d in dims] for sp in species], dtype=float)
    fig, ax = plt.subplots(figsize=(max(7, 0.7 * len(dims) + 3), max(5, 0.32 * len(species) + 1)))
    cmap = matplotlib.colors.ListedColormap([C_REVIEW, C_OK])
    ax.imshow(mat, cmap=cmap, aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(dims)))
    ax.set_xticklabels(dims, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(species)))
    ax.set_yticklabels([("* " if s in review_sp else "") + disp[s] for s in species], fontsize=7)
    ax.set_title("FGFR2 per-species evidence matrix (pre-InterProScan)\ngreen = supported, red = review", fontweight="bold")
    legend = [Patch(facecolor=C_OK, edgecolor="black", label="evidence supported"),
              Patch(facecolor=C_REVIEW, edgecolor="black", label="review / not supported")]
    ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, frameon=False)
    fig.tight_layout()
    save(fig, outdir, "Figure_3_species_evidence_heatmap_pre_interpro")


def fig4_native_vs_normalized(rows: List[Dict[str, str]], outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))
    for r in rows:
        x = to_float(r.get("native_pair_center_distance_aa", ""))
        y = to_float(r.get("iii_slot_pair_center_distance_aa", ""))
        if x is None or y is None:
            continue
        review = r.get("final_display_class") != "main_analysis_high_confidence"
        ax.scatter(x, y, s=55, color=TAXON_COLORS.get(r.get("taxon_group", "unknown"), C_NEUTRAL),
                   edgecolor=(C_REVIEW if review else "black"), linewidth=(1.8 if review else 0.6),
                   marker=("s" if review else "o"), zorder=3)
    ax.set_xlabel("native pair-center distance (aa)")
    ax.set_ylabel("normalized III-slot pair-center distance (aa)")
    ax.set_title("Native vs normalized III-slot pair QC (pre-InterProScan)", fontweight="bold")
    ax.axhline(5, color=C_NEUTRAL, ls="--", lw=0.7)
    legend = [Patch(facecolor=c, edgecolor="black", label=t) for t, c in TAXON_COLORS.items()]
    legend.append(Patch(facecolor="white", edgecolor=C_REVIEW, label="review species (square)"))
    ax.legend(handles=legend, loc="best", fontsize=7, frameon=True)
    fig.tight_layout()
    save(fig, outdir, "Figure_4_native_vs_normalized_qc_pre_interpro")


def fig_review(rows: List[Dict[str, str]], outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, max(4, 0.42 * len(rows) + 1.5)))
    ax.axis("off")
    if not rows:
        ax.text(0.5, 0.5, "No review/supplementary species: all main-analysis high-confidence.",
                ha="center", va="center")
        save(fig, outdir, "Supplement_Figure_review_cases_pre_interpro")
        return
    cols = ["display_species_name", "taxon_group", "final_display_class", "review_reason_short", "recommended_use"]
    headers = ["species", "taxon", "display class", "review reason", "recommended use"]
    table_data = [[r.get(c, "")[:34] for c in cols] for r in rows]
    tbl = ax.table(cellText=table_data, colLabels=headers, loc="center", cellLoc="left")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7)
    tbl.scale(1, 1.3)
    for j in range(len(headers)):
        tbl[0, j].set_facecolor(C_NEUTRAL)
        tbl[0, j].set_text_props(color="white", fontweight="bold")
    for i, r in enumerate(rows, start=1):
        c = C_REVIEW if r.get("final_display_class") == "supplementary_review_not_primary_claim" else C_IIIC
        tbl[i, 2].set_text_props(color=c)
    ax.set_title("Review / supplementary FGFR2 species (pre-InterProScan)", fontweight="bold")
    save(fig, outdir, "Supplement_Figure_review_cases_pre_interpro")


def main() -> int:
    ap = argparse.ArgumentParser(description="Render pre-InterPro figures.")
    ap.add_argument("--tables", type=Path, required=True, help="Directory with prepared figure tables.")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    t = args.tables
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    fig1_framework(read_tsv(t / "figure1_framework_counts_pre_interpro.tsv"), out)
    fig2_exon_to_protein(read_tsv(t / "figure2_exon_to_protein_tracks_pre_interpro.tsv"), out)
    fig3_evidence_heatmap(read_tsv(t / "figure3_species_evidence_matrix_pre_interpro.tsv"), out)
    fig4_native_vs_normalized(read_tsv(t / "figure4_native_vs_normalized_qc_pre_interpro.tsv"), out)
    fig_review(read_tsv(t / "figure_review_cases_pre_interpro.tsv"), out)
    print(f"[OK] pre-InterPro figures written to {out} (svg/pdf/png each)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
