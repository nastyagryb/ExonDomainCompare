#!/usr/bin/env python3
"""
make_fgfr2_final_closure_figures.py — Parts B/C/D of the closure correction pass.

Regenerates Figures 2, 3, 5 and 6 (paper versions) and their figure-input tables
strictly from the FINAL single source of truth:

    13_final_pre_interpro_closure/final_pre_interpro_truth_table.tsv
    13_final_pre_interpro_closure/MSA/final_cassette_msa_boundary_projection.tsv
    13_final_pre_interpro_closure/MSA/final_human_referenced_residue_agreement.tsv
    13_final_pre_interpro_closure/MSA/final_fgfr2_full_length_protein_msa.aln.faa

Inclusion / review styling is driven ONLY by final_claim_status_after_rescue and
pre_interpro_readiness_class. Stale Step-11 recommended_use / is_review_species are
NOT used as figure logic. Rescued-and-validated primary rows are drawn as accepted
primary rows (no review outline). No domain boxes are drawn (pre-InterProScan).

Outputs (into the closure dir):
    figures/Figure_2_final_exon_to_protein_architecture_pre_interpro.{svg,pdf,png}
    figures/Figure_3_final_IIIb_IIIc_cassette_zoom_pre_interpro.{svg,pdf,png}
    figures/Figure_5_full_length_FGFR2_MSA_integrity_paper.{svg,pdf,png}
    figures/Figure_6_human_referenced_IIIb_IIIc_residue_agreement.{svg,pdf,png}
    tables/figure2_final_exon_to_protein_architecture_pre_interpro.tsv
    tables/figure3_final_IIIb_IIIc_cassette_zoom_pre_interpro.tsv
    tables/figure5_full_length_FGFR2_MSA_integrity_paper.tsv
    tables/figure6_human_referenced_IIIb_IIIc_residue_agreement.tsv
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402
from fgfr2 import human_reference_control as HRC  # noqa: E402
from exondomaincompare.presentation import fgfr2_plot_style as S  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl_closure_"))
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# isoform cassette colours (colour-blind safe, distinct from review amber)
ISO_COLOR = {"IIIb": "#0072B2", "IIIc": "#009E73"}
PROTEIN_BODY = "#D7DBE0"
REVIEW_EDGE = "#E69F00"
EXCLUDED_EDGE = "#D55E00"

AGREE_COLOR = {
    "identical_to_human": "#1B6CA8",
    "conservative_substitution": "#56B4E9",
    "nonconservative_substitution": "#E69F00",
    "gap_or_missing": "#E8EAED",
}

# A cassette position the run carries no observation for. Distinct from
# ``gap_or_missing`` (an observed gap) and never white, so "no data here" and
# "nothing rendered at all" cannot look the same.
NOT_COVERED_COLOR = "#F5F0E1"

# ---- amino-acid property + substitution helpers (descriptive only, no functional claims) ----
AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}
AA_PROPERTY = {
    **{a: "hydrophobic" for a in "AVLIM"},
    **{a: "aromatic" for a in "FWY"},
    **{a: "polar" for a in "STNQ"},
    **{a: "positive" for a in "KRH"},
    **{a: "negative" for a in "DE"},
    **{a: "special_case" for a in "CGP"},
}
PROPERTY_COLOR = {
    "hydrophobic": "#0072B2", "aromatic": "#6A3D9A", "polar": "#117733",
    "positive": "#CC3311", "negative": "#EE7733", "special_case": "#777777", "gap": "#BBBBBB",
}
# Clustal-style conservation groups
_STRONG = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"]
_WEAK = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY"]
POSITION_CLASS_COLOR = {
    "isoform_discriminating_conserved": "#FFF2CC",
    "IIIc_specific_conserved": "#FCE3C8",
    "IIIb_specific_conserved": "#D7E8F5",
    "shared_conserved": "#DCEDE9",
    "variable": "#ECEDEF",
    "gap_rich_review": "#F6F6F6",
}


def _aa_property(aa: str) -> str:
    aa = (aa or "").strip().upper()
    if aa in ("", "-", ".", "X"):
        return "gap"
    return AA_PROPERTY.get(aa, "special_case")


def _substitution_class(a: str, b: str) -> str:
    a = (a or "").strip().upper()
    b = (b or "").strip().upper()
    if a in ("", "-", ".") or b in ("", "-", "."):
        return "gap"
    if a == b:
        return "identical"
    if any(a in g and b in g for g in _STRONG):
        return "conservative"
    if any(a in g and b in g for g in _WEAK):
        return "semi_conservative"
    return "nonconservative"


def _is_primary(row: Dict[str, str]) -> bool:
    return M.claim_is_primary(str(row.get("final_claim_status_after_rescue", "")))


def _plot_status(row: Dict[str, str]) -> str:
    claim = str(row.get("final_claim_status_after_rescue", ""))
    if claim == "primary_claim_supported":
        return "primary_accepted"
    if claim == "primary_claim_supported_with_minor_flags":
        return "primary_minor_flags"
    if claim.startswith("excluded"):
        return "excluded"
    return "supplement_review"


def _review_flag(row: Dict[str, str]) -> str:
    return "false" if _is_primary(row) else "true"


def _short(name: str) -> str:
    """Abbreviate a binomial to 'G. species' for compact y labels."""
    parts = (name or "").replace("_", " ").split()
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    return name or ""


def load_truth(cdir: Path) -> List[Dict[str, str]]:
    rows = M.read_tsv(cdir / "final_pre_interpro_truth_table.tsv")
    if not rows:
        raise RuntimeError("final_pre_interpro_truth_table.tsv missing or empty")
    return rows


def load_boundary(cdir: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    out = {}
    for r in M.read_tsv(cdir / "MSA" / "final_cassette_msa_boundary_projection.tsv"):
        out[((r.get("species") or "").lower(), r.get("isoform") or "")] = r
    return out


def _order_rows(truth: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Primary rows first (truth is already phylo-sorted), then supplement/excluded."""
    primary = [r for r in truth if _is_primary(r)]
    review = [r for r in truth if not _is_primary(r)]
    return primary, review


# ---------------------------------------------------------------------------
# Figure 2 — final exon-to-protein architecture
# ---------------------------------------------------------------------------
def figure2(base: Path, cdir: Path) -> List[Dict[str, str]]:
    truth = load_truth(cdir)
    boundary = load_boundary(cdir)
    primary, review = _order_rows(truth)

    table: List[Dict[str, str]] = []
    for r in truth:
        k = (r["species"].lower(), r["isoform"])
        b = boundary.get(k, {})
        table.append({
            "species": r["species"], "isoform": r["isoform"],
            "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
            "pre_interpro_readiness_class": r.get("pre_interpro_readiness_class", ""),
            "rescue_decision": r.get("rescue_decision", ""),
            "visual_review_flag": _review_flag(r),
            "final_plot_status": _plot_status(r),
            "source_coordinate_table": "MSA/final_cassette_msa_boundary_projection.tsv",
            "final_isoform_label": r.get("final_isoform_label", ""),
            "protein_length": r.get("protein_length", ""),
            "native_cassette_start_aa": b.get("native_cassette_start_aa", ""),
            "native_cassette_end_aa": b.get("native_cassette_end_aa", ""),
        })
    M.write_tsv(cdir / "tables" / "figure2_final_exon_to_protein_architecture_pre_interpro.tsv",
                table, list(table[0].keys()))

    S.apply_rcparams()
    blocks = [("Primary (accepted, incl. rescued)", primary)]
    if review:
        blocks.append(("Supplement / review", review))
    n = len(truth) + len(blocks)
    fig, ax = plt.subplots(figsize=(9.5, max(5.5, n * 0.26 + 1.2)))
    max_len = max((M.to_int(r.get("protein_length"), 0) or 0) for r in truth) or 850

    y = 0
    yticks, ylabels = [], []
    for title_txt, block in blocks:
        ax.text(-0.02 * max_len, y + 0.4, title_txt, fontsize=S.FONT["label"],
                fontweight="bold", color=S.MUTED, va="center", ha="left")
        y += 1
        for r in block:
            k = (r["species"].lower(), r["isoform"])
            b = boundary.get(k, {})
            plen = M.to_int(r.get("protein_length"), 0) or 0
            iso = r.get("final_isoform_label", r["isoform"])
            is_primary = _is_primary(r)
            ax.add_patch(Rectangle((0, y + 0.12), plen, 0.56, facecolor=PROTEIN_BODY,
                                   edgecolor="#B5BAC2", lw=0.4, zorder=2))
            cs = M.to_int(b.get("native_cassette_start_aa"), 0) or 0
            ce = M.to_int(b.get("native_cassette_end_aa"), 0) or 0
            if ce > cs > 0:
                ax.add_patch(Rectangle((cs, y + 0.02), ce - cs, 0.76,
                                       facecolor=ISO_COLOR.get(iso, "#7B6FB0"),
                                       edgecolor="#1A1A1A", lw=0.5, zorder=3))
            if not is_primary:
                ax.add_patch(Rectangle((-2, y), max_len + 4, 0.8, fill=False,
                                       edgecolor=EXCLUDED_EDGE if _plot_status(r) == "excluded"
                                       else REVIEW_EDGE, lw=1.1, ls=":", zorder=4, clip_on=False))
            yticks.append(y + 0.4)
            ylabels.append(f"{_short(r.get('display_species_name', r['species']))} {iso}")
            y += 1

    ax.set_xlim(-0.04 * max_len, max_len * 1.04)
    ax.set_ylim(0, y)
    ax.invert_yaxis()
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=S.FONT["small"])
    ax.set_xlabel("FGFR2 protein position (aa)", fontsize=S.FONT["label"])
    ax.set_title("FGFR2 exon-to-protein architecture (pre-InterPro; no domain calls)",
                 fontsize=S.FONT["title"], fontweight="bold", loc="left", pad=10)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(left=False)
    handles = [Patch(facecolor=PROTEIN_BODY, label="protein body"),
               Patch(facecolor=ISO_COLOR["IIIb"], label="IIIb cassette"),
               Patch(facecolor=ISO_COLOR["IIIc"], label="IIIc cassette"),
               Line2D([0], [0], color=REVIEW_EDGE, ls=":", lw=1.1, label="supplement/review")]
    S.compact_legend(ax, handles, ncol=4, bbox=(0.5, -0.06))
    S.savefig(fig, cdir / "figures", "Figure_2_final_exon_to_protein_architecture_pre_interpro")
    print("[OK] Figure 2 (exon-to-protein architecture)")
    return table


# ---------------------------------------------------------------------------
# Figure 3 — IIIb/IIIc cassette zoom
# ---------------------------------------------------------------------------
def figure3(base: Path, cdir: Path) -> List[Dict[str, str]]:
    truth = load_truth(cdir)
    boundary = load_boundary(cdir)
    primary, review = _order_rows(truth)

    table: List[Dict[str, str]] = []
    for r in truth:
        k = (r["species"].lower(), r["isoform"])
        b = boundary.get(k, {})
        cs = M.to_int(b.get("native_cassette_start_aa"), 0) or 0
        ce = M.to_int(b.get("native_cassette_end_aa"), 0) or 0
        table.append({
            "species": r["species"], "isoform": r["isoform"],
            "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
            "pre_interpro_readiness_class": r.get("pre_interpro_readiness_class", ""),
            "rescue_decision": r.get("rescue_decision", ""),
            "visual_review_flag": _review_flag(r),
            "final_plot_status": _plot_status(r),
            "source_coordinate_table": "MSA/final_cassette_msa_boundary_projection.tsv",
            "final_isoform_label": r.get("final_isoform_label", ""),
            "cassette_length_aa": str(ce - cs + 1) if ce > cs > 0 else "",
            "cassette_msa_start_col": b.get("cassette_msa_start_col", ""),
            "cassette_msa_end_col": b.get("cassette_msa_end_col", ""),
        })
    M.write_tsv(cdir / "tables" / "figure3_final_IIIb_IIIc_cassette_zoom_pre_interpro.tsv",
                table, list(table[0].keys()))

    S.apply_rcparams()
    blocks = [("Primary (accepted, incl. rescued)", primary)]
    if review:
        blocks.append(("Supplement / review", review))
    n = len(truth) + len(blocks)
    fig, ax = plt.subplots(figsize=(8.5, max(5.5, n * 0.26 + 1.2)))
    max_clen = 0
    for r in truth:
        b = boundary.get((r["species"].lower(), r["isoform"]), {})
        cs = M.to_int(b.get("native_cassette_start_aa"), 0) or 0
        ce = M.to_int(b.get("native_cassette_end_aa"), 0) or 0
        max_clen = max(max_clen, (ce - cs + 1) if ce > cs > 0 else 0)
    max_clen = max_clen or 70

    y = 0
    yticks, ylabels = [], []
    for title_txt, block in blocks:
        ax.text(-0.02 * max_clen, y + 0.4, title_txt, fontsize=S.FONT["label"],
                fontweight="bold", color=S.MUTED, va="center", ha="left")
        y += 1
        for r in block:
            b = boundary.get((r["species"].lower(), r["isoform"]), {})
            cs = M.to_int(b.get("native_cassette_start_aa"), 0) or 0
            ce = M.to_int(b.get("native_cassette_end_aa"), 0) or 0
            clen = (ce - cs + 1) if ce > cs > 0 else 0
            iso = r.get("final_isoform_label", r["isoform"])
            is_primary = _is_primary(r)
            if clen > 0:
                ax.add_patch(Rectangle((0, y + 0.12), clen, 0.56,
                                       facecolor=ISO_COLOR.get(iso, "#7B6FB0"),
                                       edgecolor="#1A1A1A", lw=0.5, zorder=3))
                ax.text(clen + 0.5, y + 0.4, f"{clen} aa", fontsize=S.FONT["small"],
                        va="center", ha="left", color=S.MUTED)
            if not is_primary:
                ax.add_patch(Rectangle((-1.2, y), max_clen + 12, 0.8, fill=False,
                                       edgecolor=EXCLUDED_EDGE if _plot_status(r) == "excluded"
                                       else REVIEW_EDGE, lw=1.1, ls=":", zorder=4, clip_on=False))
            yticks.append(y + 0.4)
            ylabels.append(f"{_short(r.get('display_species_name', r['species']))} {iso}")
            y += 1

    ax.set_xlim(-0.04 * max_clen, max_clen * 1.18)
    ax.set_ylim(0, y)
    ax.invert_yaxis()
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=S.FONT["small"])
    ax.set_xlabel("Ig-III alternative cassette length (aa)", fontsize=S.FONT["label"])
    ax.set_title("FGFR2 IIIb/IIIc cassette zoom (pre-InterPro)",
                 fontsize=S.FONT["title"], fontweight="bold", loc="left", pad=10)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(left=False)
    handles = [Patch(facecolor=ISO_COLOR["IIIb"], label="IIIb cassette"),
               Patch(facecolor=ISO_COLOR["IIIc"], label="IIIc cassette"),
               Line2D([0], [0], color=REVIEW_EDGE, ls=":", lw=1.1, label="supplement/review")]
    S.compact_legend(ax, handles, ncol=3, bbox=(0.5, -0.06))
    S.savefig(fig, cdir / "figures", "Figure_3_final_IIIb_IIIc_cassette_zoom_pre_interpro")
    print("[OK] Figure 3 (cassette zoom)")
    return table


# ---------------------------------------------------------------------------
# Figure 5 — full-length MSA integrity (real alignment figure)
# ---------------------------------------------------------------------------
def _per_column_profiles(aln: List[Tuple[str, str]]) -> Tuple[List[float], List[float]]:
    if not aln:
        return [], []
    width = len(aln[0][1])
    occupancy = [0.0] * width
    conservation = [0.0] * width
    n = len(aln)
    for col in range(width):
        counts: Counter = Counter()
        non_gap = 0
        for _, seq in aln:
            aa = seq[col] if col < len(seq) else "-"
            if aa not in ("-", ".", " "):
                counts[aa] += 1
                non_gap += 1
        occupancy[col] = non_gap / n if n else 0.0
        conservation[col] = (counts.most_common(1)[0][1] / non_gap) if non_gap else 0.0
    return occupancy, conservation


def figure5(base: Path, cdir: Path) -> List[Dict[str, str]]:
    truth = load_truth(cdir)
    truth_k = {(r["species"].lower(), r["isoform"]): r for r in truth}
    boundary = load_boundary(cdir)
    integ = {((r.get("species") or "").lower(), r.get("isoform") or ""): r
             for r in M.read_tsv(cdir / "MSA" / "final_full_length_msa_integrity_qc.tsv")}
    aln_all = M.read_fasta(cdir / "MSA" / "final_fgfr2_full_length_protein_msa.aln.faa")

    # map alignment record -> truth key via composite header species|isoform|pid|tag
    def hdr_key(h: str) -> Tuple[str, str]:
        parts = h.split("|")
        return (parts[0].lower() if parts else "", parts[1] if len(parts) > 1 else "")

    primary_aln = [(h, s) for h, s in aln_all
                   if _is_primary(truth_k.get(hdr_key(h), {}))]
    review_aln = [(h, s) for h, s in aln_all
                  if hdr_key(h) in truth_k and not _is_primary(truth_k[hdr_key(h)])]

    occ, cons = _per_column_profiles(aln_all)
    width = len(occ)

    # cassette band (median projected columns)
    starts = [M.to_int(b.get("full_length_msa_start_col"), 0) for b in boundary.values()
              if M.to_int(b.get("full_length_msa_start_col"), 0)]
    ends = [M.to_int(b.get("full_length_msa_end_col"), 0) for b in boundary.values()
            if M.to_int(b.get("full_length_msa_end_col"), 0)]
    cas_start = sorted(starts)[len(starts) // 2] if starts else 0
    cas_end = sorted(ends)[len(ends) // 2] if ends else 0

    # protein-length outliers
    lengths = [M.to_int(r.get("protein_length"), 0) or 0 for r in truth
               if (M.to_int(r.get("protein_length"), 0) or 0) > 0]
    med_len = sorted(lengths)[len(lengths) // 2] if lengths else 0

    table: List[Dict[str, str]] = []
    for r in truth:
        k = (r["species"].lower(), r["isoform"])
        iq = integ.get(k, {})
        plen = M.to_int(r.get("protein_length"), 0) or 0
        dev = (plen - med_len) if med_len else 0
        outlier = "true" if med_len and abs(dev) > max(40, 0.12 * med_len) else "false"
        table.append({
            "species": r["species"], "isoform": r["isoform"],
            "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
            "pre_interpro_readiness_class": r.get("pre_interpro_readiness_class", ""),
            "rescue_decision": r.get("rescue_decision", ""),
            "visual_review_flag": _review_flag(r),
            "final_plot_status": _plot_status(r),
            "source_coordinate_table": "MSA/final_fgfr2_full_length_protein_msa.aln.faa",
            "MSA_full_length_status": r.get("MSA_full_length_status", ""),
            "protein_length": str(plen),
            "length_deviation_from_median": str(dev),
            "length_outlier_flag": outlier,
            "length_outlier_status": iq.get("length_outlier_status", ""),
            "full_length_gap_fraction": r.get("full_length_gap_fraction", ""),
        })
    M.write_tsv(cdir / "tables" / "figure5_full_length_FGFR2_MSA_integrity_paper.tsv",
                table, list(table[0].keys()))

    S.apply_rcparams()
    import numpy as np
    from matplotlib.colors import ListedColormap, BoundaryNorm
    n_pri = len(primary_aln)
    n_rev = len(review_aln)
    # Row labels need vertical room; the matrix must stay readable down to the two
    # sequences per species a comparative run has.
    row_h = 0.42 if (n_pri + n_rev) <= 8 else 0.13
    fig = plt.figure(figsize=(12, max(7.0, (n_pri + n_rev) * row_h + 3.6)))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, max(1.8, n_pri * row_h * 1.6),
                                               max(0.4, n_rev * row_h)],
                          width_ratios=[5.0, 0.9], hspace=0.30, wspace=0.04,
                          left=0.20, right=0.97, top=0.86, bottom=0.10)
    ax_prof = fig.add_subplot(gs[0, 0])
    ax_ras = fig.add_subplot(gs[1, 0], sharex=ax_prof)
    ax_len = fig.add_subplot(gs[1, 1], sharey=ax_ras)
    ax_rev = fig.add_subplot(gs[2, 0], sharex=ax_prof)

    xs = list(range(1, width + 1))
    ax_prof.fill_between(xs, [1 - o for o in occ], color="#C9CDD2", lw=0, label="gap fraction")
    ax_prof.plot(xs, cons, color="#1B6CA8", lw=0.9, label="column conservation")
    if cas_end > cas_start > 0:
        ax_prof.axvspan(cas_start, cas_end, color="#009E73", alpha=0.14, lw=0, zorder=0)
        if n_rev:
            ax_rev.axvspan(cas_start, cas_end, color="#009E73", alpha=0.14, lw=0, zorder=0)
    ax_prof.set_ylim(0, 1.02)
    ax_prof.set_ylabel("frac.", fontsize=S.FONT["small"])
    ax_prof.tick_params(labelbottom=False, labelsize=S.FONT["small"])
    # Keep every legend outside the data area: anchored above its own axes, never
    # floating over the profile or the sequence matrix.
    ax_prof.legend(fontsize=S.FONT["legend"], frameon=False, ncol=2,
                   loc="lower left", bbox_to_anchor=(0.0, 1.02))
    ax_prof.set_title("Full-length FGFR2 protein MSA integrity (primary sequences; green band = IIIb/IIIc cassette)",
                      fontsize=S.FONT["title"], fontweight="bold", loc="left", pad=22)

    # Explicit sequence-state encoding. The previous version wrote residue-present as
    # 0.0 into a "Greys" image, where 0 is *white*: every residue rendered invisible
    # and only the gaps were drawn. States are now named and mapped to opaque colours,
    # and the colormap's "bad" colour is set so no cell can fall through transparent.
    STATE_RESIDUE, STATE_INTERNAL_GAP, STATE_TERMINAL = 0, 1, 2
    state_cmap = ListedColormap(["#1A1A1A", "#FFFFFF", "#E4E7EB"])
    state_cmap.set_bad("#F3B0B0")  # loud, never transparent
    state_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], state_cmap.N)

    def sequence_states(seq: str) -> np.ndarray:
        """Per-column state for one aligned row, separating internal from terminal gaps."""
        row = np.full(width, STATE_TERMINAL, dtype=float)
        residues = [j for j in range(min(len(seq), width)) if seq[j] not in ("-", ".", " ")]
        if not residues:
            return row
        first, last = residues[0], residues[-1]
        for j in range(width):
            if j < first or j > last:
                row[j] = STATE_TERMINAL
            elif j < len(seq) and seq[j] not in ("-", ".", " "):
                row[j] = STATE_RESIDUE
            else:
                row[j] = STATE_INTERNAL_GAP
        return row

    def raster(ax, aln):
        mat = np.vstack([sequence_states(seq) for _, seq in aln]) if aln else np.empty((0, width))
        if mat.size:
            ax.imshow(mat, aspect="auto", cmap=state_cmap, norm=state_norm,
                      extent=[1, width, len(aln), 0], interpolation="nearest", zorder=1)
        return mat

    def row_label(header: str) -> str:
        """species · isoform · accession, as far as the row height allows."""
        parts = header.split("|")
        species = _short((parts[0] if parts else header).replace("_", " ").title())
        isoform = parts[1] if len(parts) > 1 else ""
        accession = parts[2] if len(parts) > 2 else ""
        label = f"{species} · {isoform}" if isoform else species
        return f"{label} · {accession}" if (accession and row_h > 0.2) else label

    primary_matrix = raster(ax_ras, primary_aln)
    if cas_end > cas_start > 0 and n_pri:
        # The cassette is annotated as an outline plus a marker bar, not as a wash over
        # the matrix: a translucent fill over residue cells would dim the very data the
        # cassette region is there to point at.
        ax_ras.add_patch(Rectangle((cas_start, 0), cas_end - cas_start, n_pri, fill=False,
                                   edgecolor="#009E73", lw=1.6, zorder=4))
        ax_ras.add_patch(Rectangle((cas_start, -0.14 * n_pri), cas_end - cas_start, 0.08 * n_pri,
                                   facecolor="#009E73", edgecolor="none", zorder=5,
                                   clip_on=False))
        ax_ras.annotate(f"IIIb/IIIc cassette ({cas_start}–{cas_end})",
                        xy=((cas_start + cas_end) / 2, -0.16 * n_pri), ha="center", va="bottom",
                        fontsize=S.FONT["small"], color="#00614A", zorder=6,
                        annotation_clip=False)
    ax_ras.set_ylabel(f"primary sequences (n={n_pri})", fontsize=S.FONT["small"])
    # Every sequence gets a visible, named row rather than an anonymous band.
    ax_ras.set_yticks([i + 0.5 for i in range(n_pri)])
    ax_ras.set_yticklabels([row_label(h) for h, _ in primary_aln], fontsize=S.FONT["small"])
    ax_ras.tick_params(labelbottom=False, left=False, labelleft=True)

    # length side annotation for primary rows
    for i, (h, _) in enumerate(primary_aln):
        k = hdr_key(h)
        r = truth_k.get(k, {})
        plen = M.to_int(r.get("protein_length"), 0) or 0
        dev = (plen - med_len) if med_len else 0
        is_out = med_len and abs(dev) > max(40, 0.12 * med_len)
        ax_len.barh(i + 0.5, plen, height=0.8,
                    color="#D55E00" if is_out else "#5B9BD5", edgecolor="none")
    ax_len.axvline(med_len, color="#1A1A1A", lw=0.6, ls="--")
    ax_len.set_xlim(0, (max(lengths) if lengths else 900) * 1.05)
    ax_len.set_ylim(n_pri, 0)
    ax_len.set_xlabel("protein\nlen (aa)", fontsize=S.FONT["small"])
    ax_len.tick_params(labelleft=False, labelsize=S.FONT["small"])
    for sp in ("top", "right"):
        ax_len.spines[sp].set_visible(False)

    if n_rev:
        raster(ax_rev, review_aln)
        ax_rev.set_ylabel(f"suppl./review (n={n_rev})", fontsize=S.FONT["small"])
        ax_rev.set_yticks([i + 0.5 for i in range(n_rev)])
        ax_rev.set_yticklabels([row_label(h) for h, _ in review_aln], fontsize=S.FONT["small"])
        ax_rev.tick_params(left=False, labelleft=True, labelsize=S.FONT["small"])
        for s in ax_rev.spines.values():
            s.set_edgecolor(REVIEW_EDGE)
    else:
        ax_rev.axis("off")
    ax_rev.set_xlabel("full-length MSA column", fontsize=S.FONT["label"])

    handles = [Patch(facecolor="#1A1A1A", label="residue present"),
               Patch(facecolor="#FFFFFF", edgecolor="#B5BAC2", label="internal gap / indel"),
               Patch(facecolor="#E4E7EB", edgecolor="#B5BAC2", label="terminal missing"),
               Patch(facecolor="#009E73", alpha=0.3, label="IIIb/IIIc cassette"),
               Patch(facecolor="#D55E00", label="length outlier")]
    # Below the matrix and below the x label: the legend never covers a data row.
    ax_ras.legend(handles=handles, fontsize=S.FONT["legend"], frameon=False, ncol=5,
                  loc="upper left", bbox_to_anchor=(0.0, -0.06), borderaxespad=0.0)

    # The figure states what it contains, so a blank or truncated matrix is a failure
    # here rather than something a reader has to notice.
    if primary_matrix.shape != (n_pri, width):
        raise RuntimeError(f"MSA matrix is {primary_matrix.shape}, expected {(n_pri, width)}")
    if n_pri and not (primary_matrix == STATE_RESIDUE).any():
        raise RuntimeError("MSA matrix contains no residue-present cells")
    if n_pri != len(ax_ras.get_yticklabels()):
        raise RuntimeError(f"{n_pri} sequences but {len(ax_ras.get_yticklabels())} row labels")

    # The rendered counts travel with the figure so the TSV and the image can be
    # checked against each other instead of trusted separately.
    M.write_tsv(cdir / "tables" / "figure5_full_length_MSA_render_qc.tsv", [{
        "sequence_count": str(len(aln_all)),
        "primary_row_count": str(n_pri),
        "review_row_count": str(n_rev),
        "visible_matrix_rows": str(primary_matrix.shape[0]),
        "residue_cells": str(int((primary_matrix == STATE_RESIDUE).sum())),
        "internal_gap_cells": str(int((primary_matrix == STATE_INTERNAL_GAP).sum())),
        "terminal_cells": str(int((primary_matrix == STATE_TERMINAL).sum())),
        "alignment_length": str(width),
        "cassette_msa_start_col": str(cas_start),
        "cassette_msa_end_col": str(cas_end),
        "length_sidebar_bars": str(len(ax_len.patches)),
        "row_labels": " | ".join(row_label(h) for h, _ in primary_aln),
    }], ["sequence_count", "primary_row_count", "review_row_count", "visible_matrix_rows",
         "residue_cells", "internal_gap_cells", "terminal_cells", "alignment_length",
         "cassette_msa_start_col", "cassette_msa_end_col", "length_sidebar_bars", "row_labels"])
    if len(ax_len.patches) != n_pri:
        raise RuntimeError(f"length sidebar has {len(ax_len.patches)} bars for {n_pri} rows")

    S.savefig(fig, cdir / "figures", "Figure_5_full_length_FGFR2_MSA_integrity_paper")
    print(f"[OK] Figure 5 (full-length MSA integrity, real alignment): "
          f"{n_pri} labelled rows × {width} columns, cassette {cas_start}–{cas_end}")
    return table


# ---------------------------------------------------------------------------
# Figure 6 — human-referenced IIIb/IIIc residue agreement
# ---------------------------------------------------------------------------
def figure6(base: Path, cdir: Path) -> List[Dict[str, str]]:
    truth = load_truth(cdir)
    _truth_k = {(r["species"].lower(), r["isoform"]): r for r in truth}
    agree = M.read_tsv(cdir / "MSA" / "final_human_referenced_residue_agreement.tsv")

    # group agreement rows by (species, isoform) -> {residue_index: agreement_class}
    by_si: Dict[Tuple[str, str], Dict[int, str]] = defaultdict(dict)
    ref_iso: Dict[Tuple[str, str], str] = {}
    for a in agree:
        k = ((a.get("species") or "").lower(), a.get("isoform") or "")
        idx = M.to_int(a.get("human_reference_residue_index"), 0) or 0
        by_si[k][idx] = a.get("agreement_class", "")
        ref_iso[k] = a.get("human_reference_isoform", a.get("isoform", ""))

    # figure-input table joined with final truth
    table: List[Dict[str, str]] = []
    for r in truth:
        k = (r["species"].lower(), r["isoform"])
        cls = by_si.get(k, {})
        total = len(cls)
        ident = sum(1 for v in cls.values() if v == "identical_to_human")
        cons = sum(1 for v in cls.values() if v == "conservative_substitution")
        noncons = sum(1 for v in cls.values() if v == "nonconservative_substitution")
        gap = sum(1 for v in cls.values() if v == "gap_or_missing")
        table.append({
            "species": r["species"], "isoform": r["isoform"],
            "final_isoform_label": r.get("final_isoform_label", ""),
            "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
            "pre_interpro_readiness_class": r.get("pre_interpro_readiness_class", ""),
            "rescue_decision": r.get("rescue_decision", ""),
            "final_label_source": r.get("final_label_source", ""),
            "visual_review_flag": _review_flag(r),
            "final_plot_status": _plot_status(r),
            "source_coordinate_table": "MSA/final_human_referenced_residue_agreement.tsv",
            "human_reference_isoform": ref_iso.get(k, ""),
            "n_positions": str(total),
            "n_identical": str(ident), "n_conservative": str(cons),
            "n_nonconservative": str(noncons), "n_gap": str(gap),
            "fraction_identical_or_conservative": (
                f"{(ident + cons) / total:.3f}" if total else ""),
        })
    M.write_tsv(cdir / "tables" / "figure6_human_referenced_IIIb_IIIc_residue_agreement.tsv",
                table, list(table[0].keys()))

    # plot: residue agreement matrix, primary only, faceted IIIb (top) / IIIc (bottom)
    def panel_rows(iso_label: str) -> List[Dict[str, str]]:
        rows = [r for r in truth if r.get("final_isoform_label") == iso_label and _is_primary(r)]
        return rows

    S.apply_rcparams()
    iiib = panel_rows("IIIb")
    iiic = panel_rows("IIIc")
    # IIIb and IIIc are 46 and 48 residues; one shared maximum stretched the shorter
    # panel over positions it does not have. Each panel keeps its own cassette axis.
    panel_width = {
        iso: max([max(v.keys()) for k, v in by_si.items() if k[1] == iso and v]
                 + [HRC.EXPECTED_LENGTHS.get(iso, 48)])
        for iso in ("IIIb", "IIIc")
    }

    fig, axes = plt.subplots(2, 1, figsize=(10, max(6.0, (len(iiib) + len(iiic)) * 0.34 + 2.6)),
                             gridspec_kw={"height_ratios": [max(1, len(iiib)), max(1, len(iiic))],
                                          "hspace": 0.42})

    def draw(ax, rows, iso_label):
        width = panel_width[iso_label]
        for yi, r in enumerate(rows):
            k = (r["species"].lower(), r["isoform"])
            cls = by_si.get(k, {})
            for idx in range(1, width + 1):
                agreement = cls.get(idx, "")
                # A position the run does not cover is drawn as an explicit
                # "not covered" cell. Leaving it white made a populated matrix look
                # like an empty one and hid whether data were missing at all.
                colour = AGREE_COLOR.get(agreement, NOT_COVERED_COLOR) if agreement \
                    else NOT_COVERED_COLOR
                ax.add_patch(Rectangle((idx - 1, yi), 1, 1, facecolor=colour,
                                       edgecolor="#FFFFFF", lw=0.25, zorder=2))
        ax.set_xlim(0, width)
        ax.set_ylim(0, max(1, len(rows)))
        ax.invert_yaxis()
        ax.set_yticks([i + 0.5 for i in range(len(rows))])
        ax.set_yticklabels([_short(r.get("display_species_name", r["species"])) for r in rows],
                           fontsize=S.FONT["small"])
        ticks = sorted({1, 10, 20, 30, 40, width})
        ax.set_xticks([t - 0.5 for t in ticks])
        ax.set_xticklabels([str(t) for t in ticks], fontsize=S.FONT["small"])
        ax.set_title(f"{iso_label} cassette — human-referenced residue agreement "
                     f"(primary, n={len(rows)}; {width} cassette positions)",
                     fontsize=S.FONT["subtitle"], fontweight="bold", loc="left", pad=6)
        for sp in ("top", "right", "left", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(left=False, bottom=False)
        return sum(1 for r in rows if by_si.get((r["species"].lower(), r["isoform"])))

    covered = draw(axes[0], iiib, "IIIb") + draw(axes[1], iiic, "IIIc")
    if (iiib or iiic) and not covered:
        raise RuntimeError("human-referenced agreement figure has no populated rows")
    axes[1].set_xlabel("human reference cassette residue index", fontsize=S.FONT["label"])
    handles = [Patch(facecolor=AGREE_COLOR["identical_to_human"], label="identical"),
               Patch(facecolor=AGREE_COLOR["conservative_substitution"], label="conservative"),
               Patch(facecolor=AGREE_COLOR["nonconservative_substitution"], label="non-conservative"),
               Patch(facecolor=AGREE_COLOR["gap_or_missing"], edgecolor="#B5BAC2", label="gap/missing"),
               Patch(facecolor=NOT_COVERED_COLOR, edgecolor="#B5BAC2", label="not covered")]
    axes[1].legend(handles=handles, fontsize=S.FONT["legend"], frameon=False, ncol=5,
                   loc="upper center", bbox_to_anchor=(0.5, -0.24))
    fig.suptitle("Human-referenced IIIb/IIIc residue agreement (pre-InterPro)",
                 fontsize=S.FONT["title"], fontweight="bold", x=0.02, ha="left")
    S.savefig(fig, cdir / "figures", "Figure_6_human_referenced_IIIb_IIIc_residue_agreement")

    # supplement panel for review/supplement rows (kept honest, separated)
    rev = [r for r in truth if not _is_primary(r)]
    if rev:
        figr, axr = plt.subplots(figsize=(10, max(2.0, len(rev) * 0.3 + 1.2)))
        for yi, r in enumerate(rev):
            k = (r["species"].lower(), r["isoform"])
            cls = by_si.get(k, {})
            for idx in range(1, max_idx + 1):
                c = AGREE_COLOR.get(cls.get(idx, ""), "#FFFFFF")
                axr.add_patch(Rectangle((idx - 1, yi), 1, 1, facecolor=c, edgecolor="#FFFFFF", lw=0.25))
        axr.set_xlim(0, max_idx)
        axr.set_ylim(0, len(rev))
        axr.invert_yaxis()
        axr.set_yticks([i + 0.5 for i in range(len(rev))])
        axr.set_yticklabels([f"{_short(r.get('display_species_name', r['species']))} {r['isoform']}"
                             for r in rev], fontsize=S.FONT["small"])
        axr.set_xlabel("human reference cassette residue index", fontsize=S.FONT["label"])
        axr.set_title("Supplement / review cassettes — residue agreement",
                      fontsize=S.FONT["subtitle"], fontweight="bold", loc="left")
        for s in axr.spines.values():
            s.set_edgecolor(REVIEW_EDGE)
        axr.tick_params(left=False)
        S.savefig(figr, cdir / "figures", "Supplement_Figure_6_residue_agreement_review_cases")
    print("[OK] Figure 6 (human-referenced residue agreement)")
    return table


# ---------------------------------------------------------------------------
# Part A — Pongo IIIb / Canis IIIc review-case explanation
# ---------------------------------------------------------------------------
def review_case_explanation(base: Path, cdir: Path) -> List[Dict[str, str]]:
    truth = load_truth(cdir)
    # focus rows: the two review cases plus their accepted isoform partners for contrast
    focus = {("pongo_abelii", "IIIb"), ("pongo_abelii", "IIIc"),
             ("canis_lupus_familiaris", "IIIb"), ("canis_lupus_familiaris", "IIIc")}
    rows = []
    for r in truth:
        k = (r["species"].lower(), r["isoform"])
        if k not in focus:
            continue
        claim = str(r.get("final_claim_status_after_rescue", ""))
        rescue = str(r.get("rescue_decision", ""))
        unresolved = str(r.get("unresolved_reason_if_any", ""))
        is_primary = M.claim_is_primary(claim)
        ext_rescue = rescue.startswith("rescued")
        confirmed_screen = "confirmed_after_exhaustive_screen" in rescue
        if ext_rescue:
            rescued_found = "yes"
        elif confirmed_screen:
            rescued_found = "not_required_current_confirmed"
        else:
            rescued_found = "no"
        if is_primary and ext_rescue:
            interp = ("Rescued with an external source-compatible validated candidate and accepted as a "
                      "primary pre-InterPro case; shown as an accepted row in all main figures. "
                      "Provenance retained.")
            seqval = "external_validated_candidate_confirmed"
        elif is_primary:
            interp = ("Current candidate confirmed after exhaustive screen / sequence reconciliation; "
                      "accepted as primary (minor flags). Shown as an accepted row in main figures. "
                      "Not an external rescue and not unresolved.")
            seqval = "sequence_reconciliation_confirmed"
        else:
            interp = ("NOT rescued: locus, orthology, synteny, MSA, coordinates and protein integrity "
                      "all pass, but no source-compatible externally validated isoform-specific "
                      "candidate was found (sequence support only). The isoform-specific claim is "
                      "therefore kept as supplement/review with provenance, not asserted as primary. "
                      "This is a genuine unresolved case, not a rescued one.")
            seqval = "sequence_support_only_no_validated_candidate"
        rows.append({
            "species": r["species"], "isoform": r["isoform"],
            "final_claim_status_after_rescue": claim,
            "rescue_decision": rescue,
            "final_label_source": r.get("final_label_source", ""),
            "rescued_candidate_found": rescued_found,
            "sequence_validation_status": seqval,
            "coordinate_validation_status": r.get("coordinate_validation_status", ""),
            "MSA_status": r.get("MSA_full_length_status", ""),
            "synteny_status": r.get("synteny_validation_class", ""),
            "unresolved_reason_if_any": unresolved,
            "final_interpretation": interp,
        })
    rows.sort(key=lambda x: (x["species"], x["isoform"]))
    M.write_tsv(cdir / "tables" / "final_review_case_explanation.tsv", rows,
                ["species", "isoform", "final_claim_status_after_rescue", "rescue_decision",
                 "final_label_source", "rescued_candidate_found", "sequence_validation_status",
                 "coordinate_validation_status", "MSA_status", "synteny_status",
                 "unresolved_reason_if_any", "final_interpretation"])
    print("[OK] Part A review-case explanation table")
    return rows


# ---------------------------------------------------------------------------
# Part B — amino-acid cassette motif map
# ---------------------------------------------------------------------------
def figure3B(base: Path, cdir: Path) -> List[Dict[str, str]]:
    disc = M.read_tsv(cdir / "MSA" / "final_isoform_discriminating_residues.tsv")
    if not disc:
        raise RuntimeError("final_isoform_discriminating_residues.tsv missing")

    table: List[Dict[str, str]] = []
    for d in disc:
        hb, hc = d.get("human_IIIb_aa", ""), d.get("human_IIIc_aa", "")
        bmaj, cmaj = d.get("IIIb_major_aa", ""), d.get("IIIc_major_aa", "")
        pcls = d.get("position_class", "")
        discriminating = pcls in ("isoform_discriminating_conserved", "IIIc_specific_conserved",
                                  "IIIb_specific_conserved")
        table.append({
            "human_reference_residue_index": d.get("human_reference_residue_index", ""),
            "human_IIIb_reference_index": d.get("human_IIIb_reference_index", ""),
            "human_IIIc_reference_index": d.get("human_IIIc_reference_index", ""),
            "human_reference_source": d.get("human_reference_source", ""),
            "MSA_column": d.get("alignment_col", ""),
            "human_IIIb_aa_one_letter": hb, "human_IIIc_aa_one_letter": hc,
            "human_IIIb_aa_three_letter": AA3.get(hb.upper(), "gap" if hb in ("", "-") else hb),
            "human_IIIc_aa_three_letter": AA3.get(hc.upper(), "gap" if hc in ("", "-") else hc),
            "IIIb_major_aa": bmaj, "IIIc_major_aa": cmaj,
            "IIIb_residue_property_class": _aa_property(hb),
            "IIIc_residue_property_class": _aa_property(hc),
            "substitution_class_IIIb_vs_IIIc": _substitution_class(hb, hc),
            "position_class": pcls,
            "is_isoform_discriminating": "true" if discriminating else "false",
            "is_shared_conserved": "true" if pcls == "shared_conserved" else "false",
            "discriminating_score": d.get("discriminating_score", ""),
            "informative_column": d.get("informative_column", ""),
            "gap_rich_excluded_from_main_plot": d.get("gap_rich_excluded_from_main_plot", ""),
            "uses_final_isoform_label": "true",  # built from corrected post-rescue IIIb/IIIc consensus
        })
    M.write_tsv(cdir / "tables" / "figure3B_IIIb_IIIc_cassette_amino_acid_motif_map.tsv",
                table, list(table[0].keys()))

    # The cassette axis is the whole cassette. Dropping gap-rich columns here used to
    # shrink the axis to whatever survived, so a run with sparse columns drew a map
    # narrower than the cassette it claims to show.
    main = table
    n = len(main)

    # "Consensus" over two species is a majority of two, not a vertebrate consensus.
    # The row label says how many species stand behind it.
    n_species = len({r["species"] for r in load_truth(cdir)})
    observed = ("consensus" if n_species >= 10 else f"major state (n={n_species} species)")

    S.apply_rcparams()
    row_defs = [("Human IIIb reference", "human_IIIb_aa_one_letter", ISO_COLOR["IIIb"]),
                (f"IIIb {observed}", "IIIb_major_aa", ISO_COLOR["IIIb"]),
                ("Human IIIc reference", "human_IIIc_aa_one_letter", "#D55E00"),
                (f"IIIc {observed}", "IIIc_major_aa", "#D55E00")]
    nrows = len(row_defs)
    fig, ax = plt.subplots(figsize=(max(11, n * 0.34 + 2.0), 4.6))

    for xi, t in enumerate(main):
        pcls = t["position_class"]
        bg = POSITION_CLASS_COLOR.get(pcls, "#FFFFFF")
        disc = t["is_isoform_discriminating"] == "true"
        # column background spanning the four residue rows
        ax.add_patch(Rectangle((xi, 0), 1, nrows, facecolor=bg, edgecolor="white", lw=0.4, zorder=1))
        if disc:
            ax.add_patch(Rectangle((xi + 0.04, 0.04), 0.92, nrows - 0.08, fill=False,
                                   edgecolor="#B8860B", lw=1.3, zorder=4))
            ax.add_patch(Rectangle((xi, nrows + 0.08), 1, 0.18, facecolor="#B8860B",
                                   edgecolor="none", zorder=4))
        for ri, (_, key, _) in enumerate(row_defs):
            yi = nrows - 1 - ri
            aa = (t.get(key, "") or "").upper()
            letter = aa if aa and aa not in ("-", ".") else "–"
            col = PROPERTY_COLOR.get(_aa_property(aa), "#1A1A1A")
            ax.text(xi + 0.5, yi + 0.5, letter, fontsize=8.6, fontweight="bold",
                    ha="center", va="center", color=col, zorder=5)

    ax.set_xlim(0, n)
    ax.set_ylim(0, nrows + 0.4)
    ax.set_xticks([i + 0.5 for i in range(n)])
    # IIIb and IIIc number their own cassettes, so neither numbering can label a shared
    # axis: where IIIb has a gap its index repeats. The axis is therefore the combined
    # cassette alignment column, and the per-panel indices stay in the source table.
    ax.set_xticklabels([t.get("MSA_column", "") for t in main],
                       fontsize=S.FONT["small"], rotation=0)
    ax.set_yticks([nrows - 1 - ri + 0.5 for ri in range(nrows)])
    ax.set_yticklabels([lab for lab, _, _ in row_defs], fontsize=S.FONT["label"])
    for ri, (_, _, scol) in enumerate(row_defs):
        yi = nrows - 1 - ri
        ax.add_patch(Rectangle((-0.6, yi + 0.1), 0.4, 0.8, facecolor=scol,
                               edgecolor="none", clip_on=False))
    ax.set_xlim(-0.7, n)
    ax.set_xlabel("Combined IIIb/IIIc cassette alignment column "
                  "(gold border / top tick = IIIb/IIIc-discriminating)",
                  fontsize=S.FONT["label"])
    reference_source = (main[0].get("human_reference_source") or "") if main else ""
    reference_note = ("validated human reference control"
                      if reference_source == "canonical_reference_control"
                      else "analysed human reference")
    ax.set_title(f"FGFR2 IIIb/IIIc cassette amino-acid motif map ({reference_note}; pre-InterPro)",
                 fontsize=S.FONT["title"], fontweight="bold", loc="left", pad=10)
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(left=False, bottom=False)

    prop_handles = [Patch(facecolor=PROPERTY_COLOR[p], label=p.replace("_", " "))
                    for p in ("hydrophobic", "aromatic", "polar", "positive", "negative",
                              "special_case", "gap")]
    pos_handles = [Patch(facecolor=POSITION_CLASS_COLOR["shared_conserved"], label="shared conserved"),
                   Patch(facecolor=POSITION_CLASS_COLOR["isoform_discriminating_conserved"],
                         edgecolor="#B8860B", label="discriminating"),
                   Patch(facecolor=POSITION_CLASS_COLOR["variable"], label="variable")]
    leg1 = ax.legend(handles=prop_handles, title="residue property (letter colour)",
                     loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=7,
                     fontsize=S.FONT["legend"], frameon=False, title_fontsize=S.FONT["legend"])
    ax.add_artist(leg1)
    ax.legend(handles=pos_handles, title="column class (background)", loc="upper center",
              bbox_to_anchor=(0.5, -0.42), ncol=3, fontsize=S.FONT["legend"], frameon=False,
              title_fontsize=S.FONT["legend"])
    S.savefig(fig, cdir / "figures", "Figure_3B_IIIb_IIIc_cassette_amino_acid_motif_map")
    print("[OK] Figure 3B (amino-acid cassette motif map)")
    return table


# ---------------------------------------------------------------------------
# Part C — exon-to-protein cassette coordinate map
# ---------------------------------------------------------------------------
def figure3C(base: Path, cdir: Path) -> List[Dict[str, str]]:
    truth = load_truth(cdir)
    _truth_k = {(r["species"].lower(), r["isoform"]): r for r in truth}
    boundary = load_boundary(cdir)
    tracks_path = (base / "11_publication_figures_pre_interpro" / "tables"
                   / "figure2_exon_to_protein_architecture_tracks.tsv")
    tracks = M.read_tsv(tracks_path)
    blocks: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for t in tracks:
        blocks[((t.get("species") or "").lower(), t.get("isoform") or "")].append(t)

    table: List[Dict[str, str]] = []
    for r in truth:
        k = (r["species"].lower(), r["isoform"])
        b = boundary.get(k, {})
        cs = M.to_int(b.get("native_cassette_start_aa"), 0) or 0
        ce = M.to_int(b.get("native_cassette_end_aa"), 0) or 0
        for blk in blocks.get(k, []):
            table.append({
                "species": r["species"], "isoform": r["isoform"],
                "final_isoform_label": r.get("final_isoform_label", ""),
                "validated_exon_type": r.get("validated_exon_type", ""),
                "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
                "visual_review_flag": _review_flag(r),
                "final_plot_status": _plot_status(r),
                "protein_length": r.get("protein_length", ""),
                "feature_type": blk.get("feature_type", ""),
                "exon_or_cds_id": blk.get("exon_or_cds_id", ""),
                "feature_label": blk.get("feature_label", ""),
                "block_start_aa": blk.get("protein_start_aa", ""),
                "block_end_aa": blk.get("protein_end_aa", ""),
                "is_IIIb_cassette": blk.get("is_IIIb_cassette", ""),
                "is_IIIc_cassette": blk.get("is_IIIc_cassette", ""),
                "cassette_start_aa": str(cs) if cs else "",
                "cassette_end_aa": str(ce) if ce else "",
                "boundary_left_precision": blk.get("boundary_left_precision", ""),
                "boundary_right_precision": blk.get("boundary_right_precision", ""),
                "source_coordinate_table": "11_publication_figures_pre_interpro/tables/"
                                           "figure2_exon_to_protein_architecture_tracks.tsv",
            })
    M.write_tsv(cdir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv",
                table, list(table[0].keys()) if table else ["species"])

    primary, review = _order_rows(truth)
    S.apply_rcparams()
    blocks_order = [("Primary (accepted, incl. rescued)", primary)]
    if review:
        blocks_order.append(("Supplement / review", review))
    max_len = max((M.to_int(r.get("protein_length"), 0) or 0) for r in truth) or 850
    nlines = len(truth) + len(blocks_order)
    fig, ax = plt.subplots(figsize=(10, max(6.0, nlines * 0.26 + 1.4)))

    y = 0
    yticks, ylabels = [], []
    for title_txt, block in blocks_order:
        ax.text(-0.02 * max_len, y + 0.4, title_txt, fontsize=S.FONT["label"],
                fontweight="bold", color=S.MUTED, va="center", ha="left")
        y += 1
        for r in block:
            k = (r["species"].lower(), r["isoform"])
            iso = r.get("final_isoform_label", r["isoform"])
            is_primary = _is_primary(r)
            # exon/CDS blocks as thin grey alternating segments
            for bi, blk in enumerate(blocks.get(k, [])):
                bs = M.to_int(blk.get("protein_start_aa"), 0) or 0
                be = M.to_int(blk.get("protein_end_aa"), 0) or 0
                if be <= bs:
                    continue
                shade = "#9AA6B2" if bi % 2 == 0 else "#C2CAD3"
                is_cas = (blk.get("is_IIIb_cassette") == "true" or blk.get("is_IIIc_cassette") == "true")
                if is_cas:
                    shade = ISO_COLOR.get(iso, "#7B6FB0")
                ax.add_patch(Rectangle((bs, y + 0.2), be - bs, 0.4, facecolor=shade,
                                       edgecolor="white", lw=0.3, zorder=3 if not is_cas else 4))
            # cassette boundary ticks
            b = boundary.get(k, {})
            cs = M.to_int(b.get("native_cassette_start_aa"), 0) or 0
            ce = M.to_int(b.get("native_cassette_end_aa"), 0) or 0
            for xc in (cs, ce):
                if xc > 0:
                    ax.plot([xc, xc], [y + 0.08, y + 0.72], color="#1A1A1A", lw=0.8, zorder=5)
            if not is_primary:
                ax.add_patch(Rectangle((-2, y), max_len + 4, 0.8, fill=False,
                                       edgecolor=EXCLUDED_EDGE if _plot_status(r) == "excluded"
                                       else REVIEW_EDGE, lw=1.0, ls=":", zorder=6, clip_on=False))
            yticks.append(y + 0.4)
            ylabels.append(f"{_short(r.get('display_species_name', r['species']))} {iso}")
            y += 1

    ax.set_xlim(-0.04 * max_len, max_len * 1.04)
    ax.set_ylim(0, y)
    ax.invert_yaxis()
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=S.FONT["small"])
    ax.set_xlabel("FGFR2 protein amino-acid coordinate (exon/CDS blocks; black ticks = cassette boundaries)",
                  fontsize=S.FONT["label"])
    ax.set_title("FGFR2 exon/CDS-to-protein cassette coordinate map (pre-InterPro)",
                 fontsize=S.FONT["title"], fontweight="bold", loc="left", pad=10)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(left=False)
    handles = [Patch(facecolor="#9AA6B2", label="exon/CDS block"),
               Patch(facecolor=ISO_COLOR["IIIb"], label="IIIb cassette block"),
               Patch(facecolor=ISO_COLOR["IIIc"], label="IIIc cassette block"),
               Line2D([0], [0], color="#1A1A1A", lw=0.8, label="cassette boundary"),
               Line2D([0], [0], color=REVIEW_EDGE, ls=":", lw=1.0, label="supplement/review")]
    S.compact_legend(ax, handles, ncol=5, bbox=(0.5, -0.05))
    S.savefig(fig, cdir / "figures", "Figure_3C_exon_to_protein_cassette_coordinate_map")
    print("[OK] Figure 3C (exon-to-protein cassette coordinate map)")
    return table


# ---------------------------------------------------------------------------
# Part D — positive label-reconciliation & rescue summary (Figure 4)
# ---------------------------------------------------------------------------
def figure4(base: Path, cdir: Path) -> List[Dict[str, str]]:
    truth = load_truth(cdir)

    def is_corrected(r):
        # genuine label correction = required reconciliation away from the upstream label
        return str(r.get("label_consistency_status", "")) in (
            "swapped_relative_to_upstream", "ambiguous_label_review", "unresolved_no_sequence")

    n_total = len(truth)
    n_consistent = sum(1 for r in truth
                       if str(r.get("label_consistency_status", "")) in ("label_consistent", "consistent"))
    n_corrected_accepted = sum(1 for r in truth if _is_primary(r) and is_corrected(r))
    # rescued = external validated candidate recovered; confirmed = current candidate kept after screen
    n_rescued_validated = sum(1 for r in truth
                              if _is_primary(r)
                              and str(r.get("rescue_decision", "")).startswith("rescued"))
    n_confirmed_screen = sum(1 for r in truth
                             if _is_primary(r)
                             and str(r.get("rescue_decision", ""))
                             == "current_candidate_confirmed_after_exhaustive_screen")
    n_primary = sum(1 for r in truth
                    if str(r.get("final_claim_status_after_rescue", "")) == "primary_claim_supported")
    n_primary_minor = sum(1 for r in truth
                          if str(r.get("final_claim_status_after_rescue", ""))
                          == "primary_claim_supported_with_minor_flags")
    n_supplement = sum(1 for r in truth
                       if str(r.get("final_claim_status_after_rescue", "")) == "supplement_review")
    n_excluded = sum(1 for r in truth
                     if str(r.get("final_claim_status_after_rescue", "")).startswith("excluded")
                     or str(r.get("pre_interpro_readiness_class", "")) == "not_ready_unresolved")
    n_primary_ready = n_primary + n_primary_minor

    table = [
        {"category": "labels_consistent_with_upstream", "count": n_consistent,
         "group": "positive", "description": "upstream IIIb/IIIc label already consistent with sequence evidence"},
        {"category": "labels_corrected_and_accepted", "count": n_corrected_accepted,
         "group": "positive", "description": "label corrected by sequence evidence and accepted as primary"},
        {"category": "rescued_and_validated", "count": n_rescued_validated,
         "group": "positive", "description": "external validated candidate recovered, accepted as primary"},
        {"category": "confirmed_after_exhaustive_screen", "count": n_confirmed_screen,
         "group": "positive", "description": "current candidate confirmed after exhaustive screen, primary"},
        {"category": "primary_ready_total", "count": n_primary_ready,
         "group": "positive", "description": "primary or primary-with-minor-flags (InterPro-ready)"},
        {"category": "primary_claim_supported", "count": n_primary,
         "group": "positive", "description": "primary claim supported"},
        {"category": "primary_claim_supported_with_minor_flags", "count": n_primary_minor,
         "group": "positive", "description": "primary with minor flags (still accepted)"},
        {"category": "supplement_review_only", "count": n_supplement,
         "group": "review", "description": "true review: no source-compatible validated candidate"},
        {"category": "excluded_or_unresolved", "count": n_excluded,
         "group": "review", "description": "excluded / non-recoverable"},
    ]
    M.write_tsv(cdir / "tables" / "figure4_label_reconciliation_and_rescue_summary.tsv",
                table, ["category", "count", "group", "description"])

    S.apply_rcparams()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6),
                                   gridspec_kw={"width_ratios": [1.0, 1.25], "wspace": 0.32})

    # Panel A: final claim composition (mostly primary = positive blue)
    comp = [("Primary", n_primary, "#1B6CA8"),
            ("Primary +\nminor flags", n_primary_minor, "#5B9BD5"),
            ("Supplement /\nreview", n_supplement, "#E69F00"),
            ("Excluded /\nunresolved", n_excluded, "#D55E00")]
    left = 0.0
    for _lab, val, col in comp:
        axA.barh(0, val, left=left, color=col, edgecolor="white", height=0.6)
        if val:
            axA.text(left + val / 2, 0, str(val), ha="center", va="center",
                     color="white", fontweight="bold", fontsize=S.FONT["label"])
        left += val
    axA.set_xlim(0, n_total)
    axA.set_ylim(-0.6, 0.6)
    axA.set_yticks([])
    axA.set_xlabel(f"species/isoform rows (n={n_total})", fontsize=S.FONT["label"])
    axA.set_title(f"Final pre-InterPro claim composition\n{n_primary_ready}/{n_total} primary-ready "
                  f"({100*n_primary_ready/n_total:.0f}%)",
                  fontsize=S.FONT["subtitle"], fontweight="bold", loc="left")
    for sp in ("top", "right", "left"):
        axA.spines[sp].set_visible(False)
    axA.legend(handles=[Patch(facecolor=c, label=l.replace("\n", " ")) for l, _, c in comp],
               loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
               fontsize=S.FONT["legend"], frameon=False)

    # Panel B: reconciliation/rescue outcomes — positive vs true review
    bars = [("Labels consistent", n_consistent, "#1B6CA8"),
            ("Corrected -> accepted", n_corrected_accepted, "#44AA99"),
            ("Confirmed after screen", n_confirmed_screen, "#56B4E9"),
            ("Rescued -> accepted", n_rescued_validated, "#009E73"),
            ("True review /\nunresolved", n_supplement + n_excluded, "#E69F00")]
    ys = list(range(len(bars)))[::-1]
    for yi, (_lab, val, col) in zip(ys, bars):
        axB.barh(yi, val, color=col, edgecolor="white", height=0.62)
        axB.text(val + 0.3, yi, str(val), va="center", ha="left",
                 fontsize=S.FONT["label"], fontweight="bold", color=S.INK)
    axB.set_yticks(ys)
    axB.set_yticklabels([l for l, _, _ in bars], fontsize=S.FONT["label"])
    axB.set_xlim(0, max(n_consistent, n_corrected_accepted, n_rescued_validated,
                        n_supplement + n_excluded) * 1.25 + 1)
    axB.set_xlabel("rows", fontsize=S.FONT["label"])
    axB.set_title("Annotation-aware outcomes\n(corrected & rescued shown as gains, not uncertainty)",
                  fontsize=S.FONT["subtitle"], fontweight="bold", loc="left")
    for sp in ("top", "right"):
        axB.spines[sp].set_visible(False)
    fig.suptitle("FGFR2 IIIb/IIIc label reconciliation and rescue summary (pre-InterPro)",
                 fontsize=S.FONT["title"], fontweight="bold", x=0.02, ha="left")
    S.savefig(fig, cdir / "figures", "Figure_4_label_reconciliation_and_rescue_summary")
    print("[OK] Figure 4 (positive reconciliation/rescue summary)")
    return table


# ---------------------------------------------------------------------------
# Figure 6B — species-by-position IIIb/IIIc cassette residue conservation heatmap
# ---------------------------------------------------------------------------
GREEN_AGREE = {
    "identical_to_human": "#1B7837",
    "conservative_substitution": "#A6DBA0",
    "nonconservative_substitution": "#E69F00",
    "gap_or_missing": "#DADADA",
}


def figure6B(base: Path, cdir: Path) -> List[Dict[str, str]]:
    truth = load_truth(cdir)
    _truth_k = {(r["species"].lower(), r["isoform"]): r for r in truth}
    ag = M.read_tsv(cdir / "MSA" / "final_human_referenced_residue_agreement.tsv")
    disc = M.read_tsv(cdir / "MSA" / "final_isoform_discriminating_residues.tsv")
    # One set per panel. `human_reference_residue_index` is a single collapsed column
    # that repeats a number where IIIc carries its two-residue insertion, so using it
    # for both panels marked IIIb at 16 and 17 — combined columns that exist only in
    # IIIc — and shifted every later IIIc mark. Each column is mapped through the
    # residue index of the panel it is drawn on instead.
    # A closure whose MSA step predates the two per-panel index columns carries only the
    # combined column, so the indices are recovered from that alignment first — the file
    # holds every cassette column, which is what makes counting each panel's residues
    # possible. Without it such a run simply loses its overlay.
    disc_indexed = disc if HRC.has_panel_indices(disc) \
        else HRC.panel_indices_from_combined_alignment(disc)
    discrim_by_panel = HRC.discriminating_positions_by_panel(disc_indexed)

    # agreement keyed by (species, final_isoform_label, position)
    cell: Dict[Tuple[str, str, int], Dict[str, str]] = {}
    href: Dict[Tuple[str, int], str] = {}
    # The axis is the cassette, so it starts at the validated cassette length and only
    # grows if this run observes positions beyond it. Deriving it from the observations
    # alone collapsed the axis to 0–1 whenever the agreement table came out empty.
    maxpos = dict(HRC.EXPECTED_LENGTHS)
    for r in ag:
        iso = r.get("final_isoform_label", "")  # final, not upstream
        if iso not in ("IIIb", "IIIc"):
            continue
        pos = M.to_int(r.get("human_reference_residue_index"), 0) or 0
        if not pos:
            continue
        sp = (r.get("species") or "").lower()
        cell[(sp, iso, pos)] = {"agreement_class": r.get("agreement_class", ""),
                                "species_aa": r.get("species_aa", ""),
                                "human_reference_aa": r.get("human_reference_aa", "")}
        href[(iso, pos)] = r.get("human_reference_aa", "")
        maxpos[iso] = max(maxpos.get(iso, 0), pos)

    # The human reference row is the canonical control wherever the run's own
    # agreement rows do not supply it, so the row is never blank while a validated
    # reference exists.
    try:
        control = HRC.load()
    except HRC.ReferenceControlError:
        control = None
    if control is not None:
        for iso in ("IIIb", "IIIc"):
            for res in HRC.panel_residues(control, iso):
                href.setdefault((iso, int(res["i"])), res.get("aa", ""))

    # figure-input table
    table: List[Dict[str, str]] = []
    for r in truth:
        sp, iso = r["species"].lower(), r.get("final_isoform_label", r["isoform"])
        if iso not in ("IIIb", "IIIc"):
            continue
        panel = "main" if _is_primary(r) else "supplement"
        for pos in range(1, maxpos.get(iso, 0) + 1):
            c = cell.get((sp, iso, pos))
            # A position without an observation is written out as an explicit
            # not-covered row rather than dropped, so the table and the figure agree
            # on the cassette axis instead of the table quietly ending early.
            table.append({
                "species": r["species"], "isoform": r["isoform"],
                "final_isoform_label": iso,
                "validated_exon_type": r.get("validated_exon_type", ""),
                "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
                "visual_review_flag": _review_flag(r),
                "panel": panel,
                "human_reference_residue_index": str(pos),
                "human_reference_aa": c["human_reference_aa"] if c else href.get((iso, pos), ""),
                "species_aa": c["species_aa"] if c else "",
                "agreement_class": c["agreement_class"] if c else "not_covered",
                "is_discriminating_position":
                    "true" if pos in discrim_by_panel.get(iso, set()) else "false",
                "source_coordinate_table": "MSA/final_human_referenced_residue_agreement.tsv",
            })
    M.write_tsv(cdir / "tables" / "figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv",
                table, list(table[0].keys()) if table else ["species"])

    def species_rows(iso_label: str, primary: bool) -> List[Dict[str, str]]:
        out = [r for r in truth if r.get("final_isoform_label") == iso_label
               and (_is_primary(r) == primary)]
        return out  # truth already phylo-sorted

    def draw_panel(ax, iso_label, rows, title):
        npos = maxpos.get(iso_label, 0)
        panel_disc = discrim_by_panel.get(iso_label, set())
        # human reference letter row ABOVE the heatmap (negative y; y-axis is inverted)
        for pos in range(1, npos + 1):
            hr = href.get((iso_label, pos), "")
            disc = pos in panel_disc
            ax.add_patch(Rectangle((pos - 1, -1.0), 1, 0.85, facecolor="#F2F2F2",
                                   edgecolor="#D0D0D0", lw=0.3))
            if hr:
                ax.text(pos - 0.5, -0.57, hr, fontsize=6.0, ha="center", va="center",
                        color="#1A1A1A", fontweight="bold")
            if disc:
                ax.add_patch(Rectangle((pos - 1, -1.18), 1, 0.16, facecolor="#B8860B",
                                       edgecolor="none"))
        ax.text(-0.4, -0.57, "human ref", fontsize=4.8, ha="right", va="center", color="#555555")
        # taxon group separators
        prev_tax = None
        for yi, r in enumerate(rows):
            tax = r.get("taxon_group", "")
            if tax != prev_tax and yi > 0:
                ax.axhline(yi, color="#FFFFFF", lw=1.4)
            prev_tax = tax
        for yi, r in enumerate(rows):
            sp = r["species"].lower()
            for pos in range(1, npos + 1):
                c = cell.get((sp, iso_label, pos))
                cls = c["agreement_class"] if c else "not_covered"
                col = GREEN_AGREE.get(cls, NOT_COVERED_COLOR)
                ax.add_patch(Rectangle((pos - 1, yi), 1, 1, facecolor=col,
                                       edgecolor="#FFFFFF", lw=0.25))
                if c and cls in ("conservative_substitution", "nonconservative_substitution"):
                    aa = (c["species_aa"] or "").strip()
                    if aa and aa not in ("-", "."):
                        ax.text(pos - 0.5, yi + 0.5, aa, fontsize=4.6, ha="center", va="center",
                                color="#1A1A1A")
                if pos in panel_disc:
                    ax.add_patch(Rectangle((pos - 1, yi), 1, 1, fill=False,
                                           edgecolor="#B8860B", lw=0.4))
        ax.set_xlim(0, npos)
        ax.set_ylim(-1.3, len(rows))
        ax.invert_yaxis()
        ax.set_yticks([i + 0.5 for i in range(len(rows))])
        ax.set_yticklabels([_short(r.get("display_species_name", r["species"])) for r in rows],
                           fontsize=5.8)
        ticks = sorted({1, 10, 20, 30, 40, npos})
        ax.set_xticks([t - 0.5 for t in ticks])
        ax.set_xticklabels([str(t) for t in ticks], fontsize=S.FONT["small"])
        ax.set_title(title, fontsize=S.FONT["subtitle"], fontweight="bold", loc="left", pad=4)
        for s in ("top", "right", "left", "bottom"):
            ax.spines[s].set_visible(False)
        ax.tick_params(left=False, bottom=False)

    S.apply_rcparams()
    iiib = species_rows("IIIb", True)
    iiic = species_rows("IIIc", True)
    fig, axes = plt.subplots(2, 1, figsize=(11, max(7, (len(iiib) + len(iiic)) * 0.17 + 2.6)),
                             gridspec_kw={"height_ratios": [len(iiib) + 1.5, len(iiic) + 1.5],
                                          "hspace": 0.18})
    draw_panel(axes[0], "IIIb", iiib,
               f"IIIb cassette — species-resolved residue conservation (primary, n={len(iiib)})")
    draw_panel(axes[1], "IIIc", iiic,
               f"IIIc cassette — species-resolved residue conservation (primary, n={len(iiic)})")
    axes[1].set_xlabel("Human reference cassette residue position "
                       "(top row = human reference aa; gold = IIIb/IIIc-discriminating)",
                       fontsize=S.FONT["label"])
    handles = [Patch(facecolor=GREEN_AGREE["identical_to_human"], label="identical to human"),
               Patch(facecolor=GREEN_AGREE["conservative_substitution"], label="conservative"),
               Patch(facecolor=GREEN_AGREE["nonconservative_substitution"], label="non-conservative"),
               Patch(facecolor=GREEN_AGREE["gap_or_missing"], label="gap / missing")]
    axes[1].legend(handles=handles, fontsize=S.FONT["legend"], frameon=False, ncol=4,
                   loc="upper center", bbox_to_anchor=(0.5, -0.12))
    fig.suptitle("Species-resolved IIIb/IIIc cassette residue conservation map (pre-InterPro)",
                 fontsize=S.FONT["title"], fontweight="bold", x=0.02, ha="left")
    S.savefig(fig, cdir / "figures", "Figure_6B_species_resolved_IIIb_IIIc_cassette_residue_map")

    # supplement panel: review/supplement rows only
    rev = [r for r in truth if not _is_primary(r) and r.get("final_isoform_label") in ("IIIb", "IIIc")]
    if rev:
        figr, axr = plt.subplots(figsize=(11, max(2.2, len(rev) * 0.4 + 1.4)))
        # use the larger position axis
        npos = max(maxpos.values())
        for yi, r in enumerate(rev):
            sp = r["species"].lower()
            iso = r.get("final_isoform_label")
            for pos in range(1, (maxpos.get(iso, 0)) + 1):
                c = cell.get((sp, iso, pos))
                cls = c["agreement_class"] if c else "gap_or_missing"
                axr.add_patch(Rectangle((pos - 1, yi), 1, 1,
                                        facecolor=GREEN_AGREE.get(cls, "#FFFFFF"),
                                        edgecolor="#FFFFFF", lw=0.25))
                if c and cls in ("conservative_substitution", "nonconservative_substitution"):
                    aa = (c["species_aa"] or "").strip()
                    if aa and aa not in ("-", "."):
                        axr.text(pos - 0.5, yi + 0.5, aa, fontsize=5.0, ha="center", va="center")
        axr.set_xlim(0, npos)
        axr.set_ylim(0, len(rev))
        axr.invert_yaxis()
        axr.set_yticks([i + 0.5 for i in range(len(rev))])
        axr.set_yticklabels([f"{_short(r.get('display_species_name', r['species']))} "
                             f"{r.get('final_isoform_label')}" for r in rev], fontsize=S.FONT["small"])
        axr.set_xlabel("Human reference cassette residue position", fontsize=S.FONT["label"])
        axr.set_title("Supplement / review cassettes — species-resolved residue conservation",
                      fontsize=S.FONT["subtitle"], fontweight="bold", loc="left")
        for s in axr.spines.values():
            s.set_edgecolor(REVIEW_EDGE)
        axr.tick_params(left=False)
        S.savefig(figr, cdir / "figures", "Supplement_Figure_6B_review_rows_cassette_residue_map")
    print("[OK] Figure 6B (species-resolved cassette residue conservation map)")
    return table


def main() -> int:
    ap = argparse.ArgumentParser(description="Final closure figures 2/3/3B/3C/4/5/6/6B from truth table.")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    cdir = M.closure_dir(base)
    (cdir / "figures").mkdir(parents=True, exist_ok=True)
    (cdir / "tables").mkdir(parents=True, exist_ok=True)
    review_case_explanation(base, cdir)
    figure2(base, cdir)
    figure3(base, cdir)
    figure3B(base, cdir)
    figure3C(base, cdir)
    figure4(base, cdir)
    figure5(base, cdir)
    figure6(base, cdir)
    figure6B(base, cdir)
    print("[DONE] closure figures 2/3/3B/3C/4/5/6/6B + review explanation regenerated from final truth table")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
