#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
from typing import Iterable, List

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

VERSION = "v2.22_final_quality_display_cleanup"

# Import the resolver from the same scripts/ directory.
def _load_resolver():
    here = Path(__file__).resolve().parent
    candidates = [here / "resolve_fgfr2_IIIb_IIIc_exons_exact_v2_22.py", here / "resolve_fgfr2_IIIb_IIIc_exons_exact_v2.py"]
    for p in candidates:
        if p.exists():
            spec = importlib.util.spec_from_file_location("fgfr2_resolver", str(p))
            mod = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(mod)
            return mod
    raise SystemExit("[ERROR] resolve_fgfr2_IIIb_IIIc_exons_exact_v2.py not found next to this script")

R = _load_resolver()

C_BLUE = "#2E86AB"
C_ORANGE = "#F18F01"
C_GREY = "#C9CED6"
C_DARK = "#272A3F"


def read_tsv(path: Path, required: bool = True) -> pd.DataFrame:
    return R.read_tsv(path, required=required)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    R.write_tsv(df, path)


def f(row: pd.Series, names: Iterable[str], default: str = "") -> str:
    return R.first(row, names, default=default)


def ti(x):
    return R.to_int(x)


def col_first(df: pd.DataFrame, names: Iterable[str]):
    return R.col_first(df, names)


def short_species(s: str) -> str:
    parts = str(s).replace("_", " ").split()
    if len(parts) >= 2:
        return f"{parts[0][0]}. {' '.join(parts[1:])}"
    return str(s)


def protein_len_for_species(resolved_rows: pd.DataFrame, cds_rows: pd.DataFrame) -> int:
    vals = []
    for _, r in resolved_rows.iterrows():
        v = ti(f(r, ["protein_length_aa", "length_aa", "protein_length"]))
        if v:
            vals.append(v)
    if not cds_rows.empty:
        p2c = col_first(cds_rows, ["protein_end_aa", "exon_protein_end_aa", "aa_end", "protein_end"])
        if p2c:
            vals.extend([ti(v) for v in cds_rows[p2c].tolist() if ti(v)])
    return max(vals) if vals else 900


def filter_background_cds(cds: pd.DataFrame, selected: pd.DataFrame, species: str, source: str, resolved_rows: pd.DataFrame) -> pd.DataFrame:
    # Choose the transcript/protein from the best resolved IIIb row, otherwise first row.
    if resolved_rows.empty:
        return pd.DataFrame()
    exact = resolved_rows[resolved_rows.get("main_figure_exact_bool", False).astype(str).str.lower().isin(["true", "1", "yes"])] if "main_figure_exact_bool" in resolved_rows.columns else pd.DataFrame()
    seed = exact.iloc[0] if not exact.empty else resolved_rows.iloc[0]
    tx = f(seed, ["transcript_id_source", "transcript_id", "transcript_id_internal"])
    protein = f(seed, ["protein_id", "translation_id_source", "translation_id"])
    sh = R.selected_for_pair(selected, species, source, tx, protein)
    sub, _ = R.cds_subset_by_identity(cds, species, source, tx, protein, sh)
    if sub.empty:
        sub = R.subset_species(cds, species, source)
    return collapse_cds_rows(sub)


def collapse_cds_rows(cds: pd.DataFrame) -> pd.DataFrame:
    if cds.empty:
        return cds
    rows = []
    # Use coding_exon_key when available; otherwise one row per genomic protein interval.
    key_col = col_first(cds, ["coding_exon_key"])
    if key_col:
        groups = cds.groupby(key_col, dropna=False)
    else:
        # Candidate key includes genomic and AA coords; okay for plotting.
        tmp = cds.copy()
        tmp["_plot_key"] = tmp.apply(R.candidate_key, axis=1)
        groups = tmp.groupby("_plot_key", dropna=False)
    for key, g in groups:
        p1s = [ti(f(r, ["protein_start_aa", "exon_protein_start_aa", "aa_start", "protein_start"])) for _, r in g.iterrows()]
        p2s = [ti(f(r, ["protein_end_aa", "exon_protein_end_aa", "aa_end", "protein_end"])) for _, r in g.iterrows()]
        p1s = [x for x in p1s if x is not None]
        p2s = [x for x in p2s if x is not None]
        if not p1s or not p2s:
            continue
        r0 = g.iloc[0]
        rows.append({
            "plot_key": str(key),
            "protein_start_aa": min(p1s),
            "protein_end_aa": max(p2s),
            "matched_exon_rank": f(r0, ["matched_exon_rank", "exon_rank", "display_coding_exon_index", "cds_rank"]),
            "raw_cds_ids": ";".join(sorted(set([f(r, ["raw_cds_ids", "cds_id", "ID", "exon_id_source"]) for _, r in g.iterrows() if f(r, ["raw_cds_ids", "cds_id", "ID", "exon_id_source"])]))),
        })
    out = pd.DataFrame(rows).sort_values(["protein_start_aa", "protein_end_aa"]) if rows else pd.DataFrame()
    if not out.empty:
        out["display_coding_exon_index"] = range(1, len(out)+1)
    return out




def bool_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1", "yes", "y"])


def species_order_from_qc(pair_qc: pd.DataFrame, mode: str = "eligible") -> List[str]:
    if pair_qc.empty:
        return []
    df = pair_qc.copy()
    if mode == "eligible":
        df = df[bool_series(df.get("main_figure_eligible", pd.Series([False]*len(df))))]
    elif mode == "gold":
        df = df[df.get("pair_status", "").astype(str).eq("gold_exact_distinct_CDS_pair")]
    priority = {
        "homo_sapiens": 0, "mus_musculus": 1, "sus_scrofa": 2,
        "bos_taurus": 3, "canis_lupus_familiaris": 4, "felis_catus": 5,
        "gallus_gallus": 6, "xenopus_tropicalis": 7, "danio_rerio": 8,
        "oreochromis_niloticus": 9, "alligator_mississippiensis": 10,
        "monodelphis_domestica": 11, "ornithorhynchus_anatinus": 12,
    }
    df["_prio"] = df["species_canonical"].astype(str).map(lambda x: priority.get(x, 99))
    df = df.sort_values(["_prio", "pair_coordinate_sanity", "species_canonical"])
    return df["species_canonical"].astype(str).tolist()


def representative_species(pair_qc: pd.DataFrame, max_n: int = 16) -> List[str]:
    eligible = species_order_from_qc(pair_qc, mode="eligible")
    if len(eligible) <= max_n:
        return eligible
    fixed = [s for s in [
        "homo_sapiens", "mus_musculus", "sus_scrofa", "bos_taurus", "gallus_gallus",
        "xenopus_tropicalis", "danio_rerio", "oreochromis_niloticus", "alligator_mississippiensis",
        "monodelphis_domestica", "ornithorhynchus_anatinus"
    ] if s in eligible]
    rest = [s for s in eligible if s not in fixed]
    return (fixed + rest)[:max_n]


# Pair-level QC categories are owned by the resolver data layer; plotting only
# visualises its output and must not recompute biological eligibility or labels.
build_pair_qc = R.build_pair_qc


def plot_architecture_v22(resolved: pd.DataFrame, cds: pd.DataFrame, selected: pd.DataFrame, pair_qc: pd.DataFrame, outpath: Path, title: str, species: List[str], show_review_marks: bool = False) -> None:
    species = [s for s in species if s]
    if not species:
        fig, ax = plt.subplots(figsize=(14, 3))
        ax.text(0.5, 0.5, "No species passed the selected QC filter", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        for ext in ["png", "pdf", "svg"]:
            fig.savefig(outpath.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
        plt.close(fig)
        return
    h = max(4.2, 0.58 * len(species) + 1.8)
    fig, ax = plt.subplots(figsize=(15, h))
    sp_col = col_first(resolved, ["species_canonical", "species"])
    qc_idx = pair_qc.set_index("species_canonical") if not pair_qc.empty and "species_canonical" in pair_qc.columns else pd.DataFrame()
    ymap = {sp: len(species)-1-i for i, sp in enumerate(species)}
    max_x = 0
    for sp in species:
        g = resolved[resolved[sp_col].astype(str) == sp]
        if g.empty:
            continue
        source = f(g.iloc[0], ["source_db", "source"])
        bg = filter_background_cds(cds, selected, sp, source, g)
        plen = protein_len_for_species(g, bg)
        max_x = max(max_x, plen)
        y = ymap[sp]
        ax.axhspan(y-0.43, y+0.43, color="#F7F8FA" if y % 2 == 0 else "white", zorder=0)
        ax.plot([0, plen], [y-0.28, y-0.28], lw=3.0, color=C_DARK, solid_capstyle="round", zorder=1)
        # Background coding exons: subtle, no overcrowding in main.
        for _, br in bg.iterrows():
            p1, p2 = ti(br.get("protein_start_aa")), ti(br.get("protein_end_aa"))
            if p1 is None or p2 is None:
                continue
            ax.broken_barh([(p1, max(1, p2-p1+1))], (y-0.37, 0.18), facecolors=C_GREY, edgecolors="#7A828C", linewidth=0.45, zorder=2)
            w = p2-p1+1
            if w >= 22 and len(species) <= 18:
                ax.text((p1+p2)/2, y-0.28, str(br.get("display_coding_exon_index", "")), ha="center", va="center", fontsize=5, color="#31343A", zorder=3)
        # anchors and resolved blocks
        for _, rr in g.iterrows():
            iso = f(rr, ["inferred_isoform", "isoform"])
            color = C_BLUE if iso == "IIIb" else C_ORANGE
            off = 0.19 if iso == "IIIb" else 0.01
            ws, we = ti(f(rr, ["III_region_start_aa", "region_start_aa"])), ti(f(rr, ["III_region_end_aa", "region_end_aa"]))
            if ws and we:
                ax.broken_barh([(ws, max(1, we-ws+1))], (y+off+0.13, 0.18), facecolors="none", edgecolors=color, linewidth=1.05, linestyle="--", zorder=4)
            exact = str(rr.get("main_figure_exact_bool", "False")).lower() in ["true", "1", "yes"]
            ps, pe = ti(rr.get("exon_protein_start_aa")), ti(rr.get("exon_protein_end_aa"))
            if exact and ps and pe:
                ax.broken_barh([(ps, max(1, pe-ps+1))], (y+off, 0.23), facecolors=color, edgecolors="white", linewidth=0.8, zorder=5)
                ax.text((ps+pe)/2, y+off+0.115, iso, ha="center", va="center", fontsize=6.8, color="white", fontweight="bold", zorder=6)
            elif show_review_marks:
                xpos = (we + 2) if we else plen * 0.45
                ax.text(xpos, y+off+0.08, f"{iso} review", ha="left", va="center", fontsize=5.3, color=color, zorder=6)
        # review marker for supplement
        if show_review_marks and not qc_idx.empty and sp in qc_idx.index:
            reason = str(qc_idx.loc[sp].get("review_reason", ""))
            if reason:
                ax.text(max_x + 15, y, "review", ha="left", va="center", fontsize=6, color="#7A3E00")
    ax.set_yticks([ymap[s] for s in species])
    ax.set_yticklabels([short_species(s) for s in species], fontsize=8.5)
    ax.set_xlabel("CDS-derived amino-acid coordinate in reconstructed FGFR2 coding-exon architecture")
    ax.set_title(title, fontsize=14.5, fontweight="bold", pad=14)
    ax.set_xlim(-5, max(900, max_x + (95 if show_review_marks else 50)))
    ax.grid(axis="x", color="#E5E5E5", linewidth=0.7)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    handles = [
        Patch(facecolor=C_GREY, edgecolor="#7A828C", label="coding exons reconstructed from cds_features.tsv"),
        Patch(facecolor=C_BLUE, label="resolved IIIb CDS/exon interval"),
        Patch(facecolor=C_ORANGE, label="resolved IIIc CDS/exon interval"),
        Patch(facecolor="none", edgecolor=C_BLUE, linestyle="--", label="IIIb dynamic protein-region anchor"),
        Patch(facecolor="none", edgecolor=C_ORANGE, linestyle="--", label="IIIc dynamic protein-region anchor"),
        Patch(facecolor=C_DARK, label="selected protein axis"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, frameon=True)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(outpath.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)



def plot_status(resolved: pd.DataFrame, outpath: Path) -> None:
    if resolved is None or resolved.empty:
        fig, ax = plt.subplots(figsize=(10, 2.8))
        ax.text(0.5, 0.5, "No resolved IIIb/IIIc rows available", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        for ext in ["png", "pdf", "svg"]:
            fig.savefig(outpath.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
        plt.close(fig)
        return
    if "exon_coordinate_status" in resolved.columns:
        status = resolved["exon_coordinate_status"].astype(str)
    elif "resolver_status" in resolved.columns:
        status = resolved["resolver_status"].astype(str)
    else:
        status = pd.Series(["status_not_available"] * len(resolved))
    counts = status.value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(11, max(2.7, 0.58 * len(counts) + 1.2)))
    ax.barh(counts.index, counts.values)
    for i, v in enumerate(counts.values):
        ax.text(v + max(0.15, counts.max() * 0.01), i, str(v), va="center", fontsize=9)
    ax.set_xlabel("Number of IIIb/IIIc rows")
    ax.set_title("FGFR2 IIIb/IIIc row-level coordinate status", fontsize=14, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(outpath.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_pair_qc(pair_qc: pd.DataFrame, outpath: Path) -> None:
    if pair_qc.empty:
        return
    counts = pair_qc["pair_status"].astype(str).value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(10, max(2.5, 0.55*len(counts)+1)))
    ax.barh(counts.index, counts.values)
    for i, v in enumerate(counts.values):
        ax.text(v + 0.2, i, str(v), va="center")
    ax.set_xlabel("Number of species")
    ax.set_title("FGFR2 IIIb/IIIc pair-level QC", fontsize=14, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(outpath.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report_v22(resolved: pd.DataFrame, candidates: pd.DataFrame, pair_qc: pd.DataFrame, outdir: Path, prefix: str) -> None:
    exact_rows = int(bool_series(resolved.get("main_figure_exact_bool", pd.Series(dtype=bool))).sum()) if not resolved.empty else 0
    eligible = pair_qc[pair_qc["main_figure_eligible"].astype(bool)] if not pair_qc.empty else pd.DataFrame()
    review = pair_qc[~pair_qc["main_figure_eligible"].astype(bool)] if not pair_qc.empty else pd.DataFrame()
    text = [
        "# FGFR2 IIIb/IIIc final QC and plotting report\n\n",
        f"Version: {VERSION}\n\n",
        "This current-stage analysis performs no InterPro/domain mapping. It uses the v2.21 multi-evidence resolver output and adds final pair-level QC plus publication-oriented display filtering.\n\n",
        f"Resolved IIIb/IIIc rows: {len(resolved)}\n\n",
        f"Exact IIIb/IIIc rows: {exact_rows}\n\n",
        f"Main-figure eligible species: {len(eligible)}\n\n",
        f"Review/supplement species: {len(review)}\n\n",
        "## Main-figure eligible species\n",
    ]
    for sp in eligible["species_canonical"].astype(str).tolist() if not eligible.empty else []:
        text.append(f"- {sp}\n")
    text.append("\n## Review/supplement species and reason\n")
    if review.empty:
        text.append("- none\n")
    else:
        for _, r in review.iterrows():
            text.append(f"- {r['species_canonical']}: {r.get('review_reason','')}\n")
    text.append("\n## Interpretation\n")
    text.append("The compact main figure contains only species with a distinct resolved IIIb/IIIc CDS pair and interpretable same-slot/adjacent CDS-derived amino-acid coordinates. The full supplement figure keeps all selected species for audit and transparency. Full genomic exons including UTRs are not plotted on the amino-acid axis; the grey background tracks are coding-exon intervals reconstructed from CDS features.\n")
    (outdir / f"{prefix}_v2_22_final_qc_report.md").write_text("".join(text), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="FGFR2 IIIb/IIIc v2.22 final QC + paper-oriented plots")
    ap.add_argument("--selected", required=True, type=Path)
    ap.add_argument("--exons", required=False, type=Path)
    ap.add_argument("--cds_features", required=True, type=Path)
    ap.add_argument("--pair_audit", required=True, type=Path)
    ap.add_argument("--alt_exons", required=False, type=Path)
    ap.add_argument("--species_qc", required=False, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--prefix", default="fgfr2")
    ap.add_argument("--max_main_species", type=int, default=16)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    pair = read_tsv(args.pair_audit)
    cds = read_tsv(args.cds_features)
    selected = read_tsv(args.selected)
    alt = read_tsv(args.alt_exons, required=False) if args.alt_exons else pd.DataFrame()

    resolved, candidates = R.build_resolved(pair, cds, selected, alt)
    # QC categories come from the resolver; plotting only visualises them.
    pair_qc = R.build_pair_qc(resolved)

    write_tsv(resolved, args.outdir / f"{args.prefix}_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    write_tsv(resolved, args.outdir / f"{args.prefix}_resolved_IIIb_IIIc_exon_CDS_mapping.tsv")
    write_tsv(candidates, args.outdir / f"{args.prefix}_resolved_IIIb_IIIc_candidate_scores.tsv")
    write_tsv(pair_qc, args.outdir / f"{args.prefix}_pair_level_qc_summary.tsv")

    rep = representative_species(pair_qc, max_n=args.max_main_species)
    all_eligible = species_order_from_qc(pair_qc, mode="eligible")
    all_species = species_order_from_qc(pair_qc, mode="all") if False else pair_qc.sort_values(["main_figure_eligible", "pair_coordinate_sanity", "species_canonical"], ascending=[False, True, True])["species_canonical"].astype(str).tolist()

    plot_architecture_v22(resolved, cds, selected, pair_qc, args.outdir / "Fig_MAIN_REPRESENTATIVE_FGFR2_IIIb_IIIc_exon_architecture", "FGFR2 IIIb/IIIc coding-exon architecture: high-confidence representative species", rep, show_review_marks=False)
    plot_architecture_v22(resolved, cds, selected, pair_qc, args.outdir / "Fig_MAIN_ALL_HIGH_CONFIDENCE_FGFR2_IIIb_IIIc_exon_architecture", "FGFR2 IIIb/IIIc coding-exon architecture: all high-confidence resolved species", all_eligible, show_review_marks=False)
    plot_architecture_v22(resolved, cds, selected, pair_qc, args.outdir / "SuppFig_ALL_SELECTED_FGFR2_IIIb_IIIc_exon_architecture_QC", "FGFR2 IIIb/IIIc coding-exon architecture audit: all selected species", all_species, show_review_marks=True)
    plot_status(resolved, args.outdir / "FigQC_ROW_LEVEL_coordinate_status")
    plot_pair_qc(pair_qc, args.outdir / "FigQC_PAIR_LEVEL_species_status")
    write_report_v22(resolved, candidates, pair_qc, args.outdir, args.prefix)

    print(f"[OK] v2.22 final QC outputs written to {args.outdir}")


if __name__ == "__main__":
    main()
