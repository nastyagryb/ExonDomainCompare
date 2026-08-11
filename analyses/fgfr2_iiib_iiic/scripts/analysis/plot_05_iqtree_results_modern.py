
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

from publication_style import (
    COLORS, apply_publication_style, clean_axis,
    save_figure, takehome, title_block,
)

HYPOTHESES = {
    1: "Unconstrained ML",
    2: "Isoform monophyly",
    3: "Species-pair topology",
}

def parse_args():
    p = argparse.ArgumentParser(description="Create modern IQ-TREE topology and AU-test plots.")
    p.add_argument("--results-dir", required=True)
    return p.parse_args()

def parse_au_table(path: Path) -> pd.DataFrame:
    text = path.read_text(errors="replace")
    rows = []
    lines = text.splitlines()
    header = next(
        (index for index, line in enumerate(lines) if re.search(r"Tree\s+logL\s+deltaL.*p-AU", line)),
        None,
    )
    if header is None:
        raise ValueError(f"Could not locate the AU-test table in {path}")
    for line in lines[header + 1:]:
        stripped = line.strip()
        if rows and not stripped:
            break
        if not re.match(r"^[123]\s+", stripped):
            continue
        numbers = [
            float(value)
            for value in re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", stripped)
        ]
        if len(numbers) < 3:
            continue
        tree_no = int(numbers[0])
        logl = numbers[1]
        delta = numbers[2]
        p_au = numbers[-1]
        rows.append({
            "tree": tree_no,
            "hypothesis": HYPOTHESES.get(tree_no, f"Tree {tree_no}"),
            "log_likelihood": logl,
            "delta_log_likelihood": delta,
            "p_AU": p_au,
        })
    if not rows:
        raise ValueError(f"Could not parse an AU-test table from {path}")
    table = pd.DataFrame(rows).drop_duplicates("tree").sort_values("tree")
    if table["tree"].tolist() != [1, 2, 3]:
        raise ValueError(f"Expected three tested topologies in {path}")
    return table

def draw_hypothesis_schematic(out: Path):
    fig, ax = plt.subplots(figsize=(12.5, 7.2))
    fig.subplots_adjust(top=0.76, left=0.05, right=0.97, bottom=0.08)
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 1)
    ax.axis("off")
    titles = ["Unconstrained ML", "Isoform monophyly", "Species-pair topology"]
    subtitles = [
        "No imposed grouping",
        "All IIIb sequences together;\nall IIIc sequences together",
        "IIIb and IIIc from each species\nforced to form a pair",
    ]
    for k in range(3):
        x0 = k + 0.08
        ax.add_patch(plt.Rectangle((x0, 0.12), 0.84, 0.70,
                                   facecolor="#F5F8FB",
                                   edgecolor=COLORS["grid"], linewidth=1.2))
        ax.text(x0 + 0.42, 0.75, titles[k], ha="center", va="center",
                fontsize=14, fontweight="bold", color=COLORS["navy"])
        ax.text(x0 + 0.42, 0.62, subtitles[k], ha="center", va="top",
                fontsize=10.5, color=COLORS["gray"])

        if k == 0:
            points = [(x0+0.22,0.26),(x0+0.35,0.42),(x0+0.52,0.30),(x0+0.70,0.45)]
            for px, py in points:
                ax.plot([x0+0.42, px], [0.50, py], color=COLORS["blue"], lw=2)
                ax.scatter(px, py, s=55, color=COLORS["teal"], edgecolor="white")
        elif k == 1:
            ax.plot([x0+0.22,x0+0.22],[0.25,0.48],color=COLORS["blue"],lw=3)
            ax.plot([x0+0.64,x0+0.64],[0.25,0.48],color=COLORS["orange"],lw=3)
            for yy in [0.27,0.34,0.41,0.48]:
                ax.scatter(x0+0.22,yy,s=45,color=COLORS["blue"],edgecolor="white")
                ax.scatter(x0+0.64,yy,s=45,color=COLORS["orange"],edgecolor="white")
            ax.plot([x0+0.22,x0+0.64],[0.52,0.52],color=COLORS["navy"],lw=2)
        else:
            for idx, xx in enumerate([x0+0.20,x0+0.42,x0+0.64]):
                ax.plot([xx-0.05,xx+0.05],[0.34,0.34],color=COLORS["navy"],lw=2)
                ax.scatter(xx-0.05,0.29,s=45,color=COLORS["blue"],edgecolor="white")
                ax.scatter(xx+0.05,0.29,s=45,color=COLORS["orange"],edgecolor="white")
                ax.plot([xx,xx],[0.34,0.49],color=COLORS["navy"],lw=2)
    title_block(
        fig,
        "Three explicit evolutionary hypotheses are compared by maximum likelihood and the AU test",
        "The test asks whether cassette sequences retain isoform identity across species or preferentially pair within species."
    )
    save_figure(fig, out, "modern_IQ_TREE_hypothesis_schematic")

def draw_au_plot(out: Path, table: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(11.2, 7.2))
    fig.subplots_adjust(top=0.78, left=0.24, right=0.94, bottom=0.16)
    plot = table.sort_values("delta_log_likelihood", ascending=False).copy()
    y = np.arange(len(plot))
    accepted = plot["p_AU"] >= 0.05
    colors = np.where(accepted, COLORS["teal"], COLORS["red"])
    ax.barh(y, plot["delta_log_likelihood"], color=colors, alpha=0.82)
    ax.set_yticks(y)
    ax.set_yticklabels(plot["hypothesis"], fontweight="bold")
    ax.set_xlabel("Δ log-likelihood relative to the best topology")
    clean_axis(ax, "x")
    for yi, (_, row) in zip(y, plot.iterrows()):
        ax.text(row["delta_log_likelihood"], yi,
                f'  AU p={row["p_AU"]:.3g}',
                va="center", fontsize=11, fontweight="bold",
                color=COLORS["navy"])
    title_block(
        fig,
        "The approximately unbiased test formally distinguishes competing cassette genealogies",
        "Green hypotheses are not rejected at α=0.05; red hypotheses are incompatible with the observed alignment."
    )
    save_figure(fig, out, "modern_IQ_TREE_AU_topology_test")

def tree_coordinates(tree):
    depths = tree.depths()
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    terminals = tree.get_terminals()
    y = {tip: i for i, tip in enumerate(reversed(terminals))}
    def assign(clade):
        if clade in y:
            return y[clade]
        child_y = [assign(c) for c in clade.clades]
        y[clade] = sum(child_y) / len(child_y)
        return y[clade]
    assign(tree.root)
    return depths, y

def draw_tree(out: Path):
    path = out / "unconstrained.treefile"
    if not path.exists():
        return
    tree = Phylo.read(path, "newick")
    depths, ypos = tree_coordinates(tree)
    tips = tree.get_terminals()
    fig, ax = plt.subplots(figsize=(13.5, max(8.0, 0.24*len(tips))))
    fig.subplots_adjust(top=0.90, left=0.06, right=0.76, bottom=0.06)
    def draw_clade(clade):
        x = depths[clade]
        children = clade.clades
        if children:
            ys = [ypos[c] for c in children]
            ax.plot([x, x], [min(ys), max(ys)], color=COLORS["dark"], lw=0.8)
            for child in children:
                xc = depths[child]
                yc = ypos[child]
                ax.plot([x, xc], [yc, yc], color=COLORS["dark"], lw=0.8)
                draw_clade(child)
    draw_clade(tree.root)
    xmax = max(depths.values())
    for tip in tips:
        label = tip.name.replace("__", " ")
        iso = "IIIb" if "IIIb" in label else "IIIc"
        color = COLORS["blue"] if iso == "IIIb" else COLORS["orange"]
        ax.text(depths[tip] + 0.012*xmax, ypos[tip], label,
                va="center", fontsize=8.8, color=color, fontweight="bold")
    ax.set_ylim(-1, len(tips))
    ax.set_xlim(0, xmax*1.32)
    ax.set_yticks([])
    ax.set_xlabel("Expected substitutions per site")
    clean_axis(ax, "x")
    title_block(
        fig,
        "Maximum-likelihood genealogy of FGFR2 IIIb and IIIc cassette sequences",
        "Tip colors encode isoform identity; topology support should be interpreted together with the formal AU test."
    )
    save_figure(fig, out, "modern_IQ_TREE_unconstrained_phylogram")

def main():
    args = parse_args()
    apply_publication_style()
    out = Path(args.results_dir)
    draw_hypothesis_schematic(out)
    report = out / "topology_AU_test.iqtree"
    if report.exists():
        table = parse_au_table(report)
        table.to_csv(out / "parsed_AU_topology_test.tsv", sep="\t", index=False)
        draw_au_plot(out, table)
    draw_tree(out)

if __name__ == "__main__":
    main()
