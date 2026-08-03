#!/usr/bin/env python3
"""
Render FGFR2 MSA reference figures.

Paper-level reference-guided figures (SVG/PDF/PNG):
  6C  Human-referenced IIIb/IIIc residue agreement map
  6D  MSA cassette boundary map, local zoom (alignment-column space)
  7C  Isoform-discriminating residues, informative positions
  8C  Alignment evidence stack
  Supplement  Per-species cassette difference panels

All figures use the FINAL (sequence-calibrated) isoform labels. A validation gate enforces
final_isoform_label == validated_exon_type for all plotted rows and that human/mouse controls
pass; otherwise NO figure is generated. MSA does not relabel IIIb/IIIc; no InterPro domains.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl_msaref_"))

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402


# ---- Part H colour mapping (color-blind safe; stable across figures) ----
AGREE_COLORS = {
    "identical_to_human": "#1B9E77",          # strong blue-green
    "conservative_substitution": "#A6CEE3",   # light blue
    "nonconservative_substitution": "#D55E00",  # orange/vermillion
    "gap_or_missing": "#E0E0E0",              # light grey
    "insertion_relative_to_human": "#C7A9D9",  # pale purple
    "unmapped_review": "#F0E0A0",
}
C_REVIEW_OUTLINE = "#E69F00"   # amber
C_FAIL = "#CC0000"
TAXON_BANDS = {
    "primate": "#4C72B0", "rodent": "#55A868", "other_mammal": "#8172B2",
    "bird": "#CCB974", "reptile": "#64B5CD", "amphibian": "#937860", "fish": "#DA8BC3",
    "other": "#BBBBBB",
}


def taxon_group(sp: str) -> str:
    s = sp.lower()
    prim = ("homo", "pan_trog", "gorilla", "pongo", "macaca", "callithrix")
    rod = ("mus_mus", "rattus")
    bird = ("gallus", "meleagris", "taeniopygia")
    rept = ("anolis", "chrysemys", "alligator")
    amph = ("xenopus", "ambystoma")
    fish = ("danio", "takifugu", "gasterosteus", "oreochromis", "oryzias")
    for keys, name in ((prim, "primate"), (rod, "rodent"), (bird, "bird"),
                       (rept, "reptile"), (amph, "amphibian"), (fish, "fish")):
        if any(s.startswith(k) or k in s for k in keys):
            return name
    return "other_mammal"


def save(fig, outdir: Path, stem: str) -> None:
    for ext in ("svg", "pdf", "png"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def phylo_order(base: Path, master: Dict[str, Dict[str, str]]) -> Dict[str, int]:
    order = {}
    p = M.locate(base, "species_phylogenetic_order.tsv")
    if p:
        for r in M.read_tsv(p):
            sp = (r.get("species") or r.get("species_canonical") or "").lower()
            o = M.to_int(r.get("phylo_order") or r.get("order") or r.get("rank"))
            if sp and o is not None:
                order[sp] = o
    if not order:
        for sp, r in master.items():
            o = M.to_int(r.get("phylo_order"))
            if o is not None:
                order[sp] = o
    return order


def is_review(rec_row: Dict[str, str], recommended_use: str, claim: str = "") -> bool:
    status = (rec_row or {}).get("label_consistency_status", "")
    return (status in ("ambiguous_label_review", "unresolved_no_sequence")
            or (claim and not M.claim_is_primary(claim))
            or not M.is_main_use(recommended_use))


def is_excluded(claim: str) -> bool:
    return (claim or "") == "excluded_from_primary_claim"


def load_agreement(cons: Path, iso: str) -> Dict[Tuple[str, int], Dict[str, str]]:
    out = {}
    for r in M.read_tsv(cons / f"fgfr2_{iso}_human_reference_residue_agreement.tsv"):
        idx = M.to_int(r.get("human_reference_residue_index"))
        if idx is not None:
            out[(r["species"].lower(), idx)] = r
    return out


def species_in_order(species: List[str], order: Dict[str, int]) -> List[str]:
    return sorted(set(species), key=lambda s: (order.get(s, 999), s))


def iso_primary(claims: Dict[Tuple[str, str], Dict[str, str]], sp: str, iso: str) -> bool:
    """Whether the (species, isoform) cassette holds a primary claim after maximal rescue.
    Part G: primary figures include only primary_claim_supported(_with_minor_flags) rows."""
    return M.claim_is_primary(M.claim_value(claims.get(((sp or "").lower(), iso), {})))


def fig_6c(base, dirs, master, order, recon, recommended, claims):
    cons, figd, tabd = dirs["conservation"], dirs["figures"], dirs["tables"]
    pos_sum = defaultdict(dict)  # (iso) -> {idx: pic}
    for r in M.read_tsv(cons / "fgfr2_reference_agreement_summary_by_position.tsv"):
        idx = M.to_int(r.get("human_reference_residue_index"))
        if idx is not None:
            pos_sum[r["human_reference_isoform"]][idx] = M.to_float(
                r.get("percent_identical_or_conservative"), 0.0)
    table_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(15, 8.6),
                             gridspec_kw={"width_ratios": [51, 42], "wspace": 0.16})
    for ax, iso in zip(axes, ("IIIb", "IIIc")):
        agree = load_agreement(cons, iso)
        # Part G: primary figure shows only species whose (species,iso) cassette is a primary claim
        sps = [s for s in species_in_order([k[0] for k in agree], order)
               if iso_primary(claims, s, iso)]
        maxidx = max([k[1] for k in agree], default=1)
        # top conservation track axis
        for yi, sp in enumerate(sps):
            claim_sp = M.claim_value(claims.get((sp, iso), {})) or M.species_claim(claims, sp)
            for x in range(1, maxidx + 1):
                r = agree.get((sp, x))
                if r is None:
                    color = AGREE_COLORS["gap_or_missing"]
                    acl = "gap_or_missing"
                else:
                    acl = r["agreement_class"]
                    color = AGREE_COLORS.get(acl, "#FFFFFF")
                ax.add_patch(Rectangle((x - 0.5, yi - 0.5), 1, 1, facecolor=color,
                                       edgecolor="white", linewidth=0.25))
                if r is not None:
                    table_rows.append({"isoform": iso, "species": sp,
                                       "human_reference_residue_index": x,
                                       "human_reference_aa": r.get("human_reference_aa", ""),
                                       "species_aa": r.get("species_aa", ""),
                                       "agreement_class": acl,
                                       "label_consistency_status": r.get("label_consistency_status", ""),
                                       "final_claim_status_after_rescue": claim_sp,
                                       "is_review_species": r.get("is_review_species", "")})
        # boundary markers (cassette start/end)
        ax.axvline(0.5, color="#222222", lw=1.4, ls="-")
        ax.axvline(maxidx + 0.5, color="#222222", lw=1.4, ls="-")
        # taxon sidebar + review markers
        for yi, sp in enumerate(sps):
            ax.add_patch(Rectangle((-1.8, yi - 0.5), 1.4, 1,
                                    facecolor=TAXON_BANDS[taxon_group(sp)], edgecolor="white",
                                    linewidth=0.3, clip_on=False))
            rrow = recon.get((sp, None))
            claim = M.species_claim(claims, sp)
            if is_review(rrow, recommended.get(sp, ""), claim):
                ax.add_patch(Rectangle((-0.5, yi - 0.5), maxidx + 1, 1, fill=False,
                                        edgecolor=C_REVIEW_OUTLINE, linewidth=1.3))
            if is_excluded(claim):
                ax.add_patch(Rectangle((-0.5, yi - 0.5), maxidx + 1, 1, fill=False,
                                        edgecolor=C_FAIL, linewidth=1.0, hatch="///", alpha=0.9))
        ax.set_xlim(-2.0, maxidx + 1.0)
        ax.set_ylim(len(sps) - 0.5, -2.4)
        ax.set_yticks(range(len(sps)))
        ax.set_yticklabels([master.get(s, {}).get("display_species_name", s) for s in sps],
                           fontsize=7)
        ax.set_xlabel(f"human {iso} reference residue position", fontsize=9)
        ax.set_title(f"{iso} cassette", fontsize=11, fontweight="bold",
                     color=M.C_IIIB if iso == "IIIb" else M.C_IIIC)
        # conservation track above
        for x in range(1, maxidx + 1):
            pic = pos_sum.get(iso, {}).get(x, 0.0)
            ax.add_patch(Rectangle((x - 0.5, -2.3), 1, 0.9,
                                    facecolor=plt.cm.Greys(0.25 + 0.6 * pic),
                                    edgecolor="white", linewidth=0.2, clip_on=False))
        ax.text(maxidx / 2, -2.5, "per-position identical-or-conservative (dark = high)",
                ha="center", va="bottom", fontsize=7, color="#444444")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
    handles = [Patch(facecolor=c, edgecolor="white", label=k.replace("_", " "))
               for k, c in AGREE_COLORS.items() if k != "unmapped_review"]
    handles.append(Patch(facecolor="white", edgecolor=C_REVIEW_OUTLINE,
                         label="review / supplement (not primary)"))
    handles.append(Patch(facecolor="white", edgecolor=C_FAIL, hatch="///",
                         label="excluded from primary claim"))
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Figure 6C — Human-referenced IIIb/IIIc residue agreement map",
                 fontsize=13, fontweight="bold", y=0.99)
    fig.text(0.5, 0.945, "Per-position agreement vs the curated human FGFR2 (UniProt P21802) "
             "IIIb/IIIc cassette; final sequence-calibrated labels. Primary-claim cassettes only "
             "(rescued candidates included); review/excluded cases in the supplement. "
             "InterProScan pending.", ha="center", fontsize=8.5, color="#444444")
    save(fig, figd, "Figure_6C_human_referenced_IIIb_IIIc_residue_agreement_map")
    M.write_tsv(tabd / "figure6C_human_referenced_residue_agreement_map.tsv", table_rows,
                ["isoform", "species", "human_reference_residue_index", "human_reference_aa",
                 "species_aa", "agreement_class", "label_consistency_status",
                 "final_claim_status_after_rescue", "is_review_species"])


def fig_6d(base, dirs, master, order, recon, recommended, claims):
    _cons, figd = dirs["conservation"], dirs["figures"]
    fig, axes = plt.subplots(1, 2, figsize=(15, 8.0), gridspec_kw={"wspace": 0.16})
    for ax, iso in zip(axes, ("IIIb", "IIIc")):
        rows = M.read_tsv(dirs["maps"] / f"fgfr2_{iso}_human_reference_msa_coordinate_map.tsv")
        bycol = defaultdict(dict)
        cols = set()
        for r in rows:
            c = M.to_int(r.get("alignment_col"))
            cols.add(c)
            bycol[r["species"].lower()][c] = r
        sps = [s for s in species_in_order(list(bycol), order) if iso_primary(claims, s, iso)]
        cols = sorted(cols)
        # gap fraction per column for shading
        gapfrac = {}
        for c in cols:
            vals = [bycol[sp].get(c, {}) for sp in sps]
            ng = sum(1 for v in vals if v.get("is_gap") == "true" or not v.get("residue_aa"))
            gapfrac[c] = ng / max(1, len(sps))
        for yi, sp in enumerate(sps):
            for ci, c in enumerate(cols):
                r = bycol[sp].get(c, {})
                mst = r.get("mapping_status", "")
                if r.get("is_gap") == "true" or not r.get("residue_aa"):
                    color = AGREE_COLORS["gap_or_missing"]
                elif mst == "insertion_relative_to_human":
                    color = AGREE_COLORS["insertion_relative_to_human"]
                else:
                    color = M.C_IIIB if iso == "IIIb" else M.C_IIIC
                ax.add_patch(Rectangle((ci - 0.5, yi - 0.5), 1, 1, facecolor=color,
                                       edgecolor="white", linewidth=0.2))
            rrow = recon.get((sp, None))
            if is_review(rrow, recommended.get(sp, ""), M.species_claim(claims, sp)):
                ax.add_patch(Rectangle((-0.5, yi - 0.5), len(cols), 1, fill=False,
                                        edgecolor=C_REVIEW_OUTLINE, linewidth=1.2))
        # gap-rich columns lightly shaded on top strip
        for ci, c in enumerate(cols):
            if gapfrac[c] >= 0.5:
                ax.add_patch(Rectangle((ci - 0.5, -1.6), 1, 0.8, facecolor="#BBBBBB",
                                        edgecolor="white", linewidth=0.2, clip_on=False))
        # human-reference boundary ticks (first/last mapped human position columns)
        human_cols = [M.to_int(r.get("alignment_col")) for r in rows
                      if r["species"].lower() == "homo_sapiens"
                      and r.get("human_reference_residue_index")]
        if human_cols:
            for hc in (min(human_cols), max(human_cols)):
                if hc in cols:
                    ax.axvline(cols.index(hc), color="#222222", lw=1.3)
        ax.set_xlim(-0.6, len(cols) - 0.4)
        ax.set_ylim(len(sps) - 0.5, -1.8)
        ax.set_yticks(range(len(sps)))
        ax.set_yticklabels([master.get(s, {}).get("display_species_name", s) for s in sps],
                           fontsize=7)
        ax.set_xlabel(f"{iso} cassette MSA column (L-INS-i)", fontsize=9)
        ax.set_title(f"{iso} cassette", fontsize=11, fontweight="bold",
                     color=M.C_IIIB if iso == "IIIb" else M.C_IIIC)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
    handles = [Patch(facecolor=M.C_IIIB, label="IIIb cassette residue"),
               Patch(facecolor=M.C_IIIC, label="IIIc cassette residue"),
               Patch(facecolor=AGREE_COLORS["insertion_relative_to_human"],
                     label="insertion vs human"),
               Patch(facecolor=AGREE_COLORS["gap_or_missing"], label="gap/missing"),
               Patch(facecolor="#BBBBBB", label="gap-rich column"),
               Patch(facecolor="white", edgecolor=C_REVIEW_OUTLINE, label="review species")]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Figure 6D — MSA cassette boundary map, local zoom", fontsize=13,
                 fontweight="bold", y=0.99)
    fig.text(0.5, 0.95, "Validated exon/cassette residues project to comparable local alignment "
             "columns in robust species; black ticks = human-reference cassette boundaries.",
             ha="center", fontsize=8.5, color="#444444")
    save(fig, figd, "Figure_6D_MSA_projected_boundary_local_zoom")


def fig_7c(base, dirs):
    cons, figd = dirs["conservation"], dirs["figures"]
    rows = M.read_tsv(cons / "fgfr2_IIIb_IIIc_discriminating_positions_informative.tsv")
    rows = [r for r in rows if M.to_int(r.get("human_reference_residue_index")) is not None]
    rows.sort(key=lambda r: M.to_int(r.get("human_reference_residue_index"), 0))
    if not rows:
        return
    x = list(range(len(rows)))
    labels = [r.get("human_reference_residue_index", "") for r in rows]
    fig, (axt, axs) = plt.subplots(2, 1, figsize=(max(9, len(rows) * 0.42), 6.2),
                                   gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.32})
    tracks = [("human IIIb", "human_IIIb_aa", M.C_IIIB),
              ("human IIIc", "human_IIIc_aa", M.C_IIIC),
              ("IIIb major", "IIIb_major_aa", M.C_IIIB),
              ("IIIc major", "IIIc_major_aa", M.C_IIIC)]
    for ti, (name, key, col) in enumerate(tracks):
        y = len(tracks) - ti
        for xi, r in zip(x, rows):
            disc = r.get("position_class") == "isoform_discriminating_conserved"
            axt.text(xi, y, r.get(key, "") or "-", ha="center", va="center", fontsize=8,
                     fontweight="bold" if disc else "normal", color=col)
        axt.text(-1.2, y, name, ha="right", va="center", fontsize=8.5)
    for xi, r in zip(x, rows):
        if r.get("position_class") == "isoform_discriminating_conserved":
            axt.add_patch(Rectangle((xi - 0.5, 0.4), 1, len(tracks) + 1.1, fill=False,
                                     edgecolor="#222222", linewidth=0.8, linestyle=":"))
    axt.set_xlim(-2.5, len(rows) - 0.5)
    axt.set_ylim(0.2, len(tracks) + 1.4)
    axt.axis("off")
    axt.set_title("Figure 7C — Isoform-discriminating residues (informative positions)",
                  fontsize=12, fontweight="bold")
    axs.bar(x, [M.to_float(r.get("discriminating_score"), 0.0) for r in rows],
            color=["#222222" if r.get("position_class") == "isoform_discriminating_conserved"
                   else "#AAAAAA" for r in rows], width=0.8)
    axs.set_xticks(x)
    axs.set_xticklabels(labels, fontsize=6, rotation=90)
    axs.set_ylabel("discriminating\nscore", fontsize=8)
    axs.set_xlabel("human reference cassette residue position", fontsize=9)
    axs.set_ylim(0, 1)
    for s in ("top", "right"):
        axs.spines[s].set_visible(False)
    fig.text(0.5, -0.02, "Dotted boxes / black bars: positions conserved within each isoform but "
             "different between IIIb and IIIc. Gap-rich columns excluded (see supplement). "
             "MSA supports but does not relabel IIIb/IIIc.", ha="center", fontsize=8, color="#444444")
    save(fig, figd, "Figure_7C_IIIb_IIIc_isoform_discriminating_residues_informative")


GROUP_COLORS = {
    "human_curated_positive_control": "#1B4F72", "close_primate_control": "#2E86C1",
    "known_label_risk_mammal": "#B9770E", "global_review_outlier": "#7D6608",
    "standard_species": "#BBBBBB",
}
CLAIM_COLORS = {
    "primary_claim_supported": "#1B9E77",
    "primary_claim_supported_with_minor_flags": "#A6CEE3",
    "supplement_review": "#F0E0A0", "excluded_from_primary_claim": "#D55E00",
}
RESCUE_COLORS = {
    "keep_current_candidate": "#1B9E77", "use_rescued_candidate": "#2E86C1",
    "manual_review_required": "#F0E0A0", "exclude_from_primary_claim": "#D55E00",
    "not_required": "#E0E0E0", "": "#FFFFFF",
}


def fig_8c(base, dirs, master, order):
    figd = dirs["figures"]
    scores = M.read_tsv(dirs["robustness"] / "fgfr2_boundary_robustness_scores.tsv")
    if not scores:
        return
    claims = M.load_claim_status(base)
    # Part G: primary evidence stack shows only primary-claim (species, isoform) rows
    rows = sorted([r for r in scores
                   if iso_primary(claims, r["species"].lower(), r.get("final_isoform_label", ""))],
                  key=lambda r: (order.get(r["species"].lower(), 999),
                                 r.get("final_isoform_label", "")))
    cols = [("IIIb ref agr", "reference_agreement_percent_identical_or_conservative", "frac"),
            ("IIIc ref agr", "reference_agreement_percent_identical_or_conservative", "frac"),
            ("boundary proj", "component_msa_projection_score", "frac"),
            ("left bound cons", "left_boundary_reference_agreement", "frac"),
            ("right bound cons", "right_boundary_reference_agreement", "frac"),
            ("core cons", "cassette_core_reference_agreement", "frac"),
            ("gap burden", "gap_rich_penalty", "gap"),
            ("disc support", "discriminating_residue_support", "frac")]
    # build label rows
    ylabels = [f"{master.get(r['species'].lower(), {}).get('display_species_name', r['species'])} "
               f"[{r.get('final_isoform_label','')}]" for r in rows]
    fig, ax = plt.subplots(figsize=(13, max(6, len(rows) * 0.32)))
    # columns: [group sidebar] + evidence cols + [final evidence] + [rescue status] + [final claim]
    ncol = 1 + len(cols) + 3
    for yi, r in enumerate(rows):
        sp = r["species"].lower()
        cr = claims.get((sp, r.get("final_isoform_label", "")), {})
        # validation-group sidebar
        ax.add_patch(Rectangle((0, yi), 1, 1,
                               facecolor=GROUP_COLORS.get(cr.get("validation_group", ""), "#FFFFFF"),
                               edgecolor="white", lw=0.5))
        for ci, (_name, key, kind) in enumerate(cols):
            v = M.to_float(r.get(key))
            if v is None:
                color = "#FFFFFF"
            elif kind == "gap":
                color = plt.cm.Oranges(min(1.0, max(0.0, v)))
            else:
                color = plt.cm.GnBu(0.15 + 0.8 * min(1.0, max(0.0, v)))
            ax.add_patch(Rectangle((1 + ci, yi), 1, 1, facecolor=color, edgecolor="white", lw=0.5))
        # final evidence class as discrete cell
        evid = r.get("overall_alignment_evidence_class", "")
        ec = {"alignment_supports_boundary": "#1B9E77",
              "alignment_supports_boundary_with_minor_variation": "#A6CEE3",
              "alignment_gap_rich_review": "#E0E0E0",
              "alignment_shift_review": "#D55E00",
              "alignment_unresolved": "#F0E0A0"}.get(evid, "#FFFFFF")
        ax.add_patch(Rectangle((1 + len(cols), yi), 1, 1, facecolor=ec, edgecolor="white", lw=0.5))
        # rescue status + final claim status
        ax.add_patch(Rectangle((2 + len(cols), yi), 1, 1,
                               facecolor=RESCUE_COLORS.get(cr.get("rescue_status", ""), "#FFFFFF"),
                               edgecolor="white", lw=0.5))
        cc = CLAIM_COLORS.get(cr.get("final_claim_status", ""), "#FFFFFF")
        ax.add_patch(Rectangle((3 + len(cols), yi), 1, 1, facecolor=cc, edgecolor="white", lw=0.5,
                               hatch="///" if cr.get("final_claim_status") == "excluded_from_primary_claim" else None))
    ax.set_xlim(0, ncol)
    ax.set_ylim(len(rows), 0)
    ax.set_xticks([i + 0.5 for i in range(ncol)])
    ax.set_xticklabels(["valid. group"] + [c[0] for c in cols] +
                       ["final evidence", "rescue status", "final claim"],
                       rotation=40, ha="right", fontsize=8)
    ax.set_yticks([i + 0.5 for i in range(len(rows))])
    ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Figure 8C — Alignment evidence stack (primary-claim cassettes)",
                 fontsize=12, fontweight="bold")
    handles = [Patch(facecolor="#1B9E77", label="supports / primary"),
               Patch(facecolor="#A6CEE3", label="minor variation / minor flags"),
               Patch(facecolor="#F0E0A0", label="review / supplement"),
               Patch(facecolor="#D55E00", label="shift / excluded"),
               Patch(facecolor=GROUP_COLORS["close_primate_control"], label="primate control"),
               Patch(facecolor=GROUP_COLORS["known_label_risk_mammal"], label="known-risk mammal"),
               Patch(facecolor=GROUP_COLORS["global_review_outlier"], label="review outlier"),
               Patch(facecolor=GROUP_COLORS["standard_species"], label="standard species")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.08))
    fig.text(0.5, -0.11, "Columns: validation group, reference agreement, boundary projection and "
             "conservation (blue-green = stronger), gap burden (orange), discriminating support, "
             "final evidence, rescue status and final claim status. Hatched = excluded from primary "
             "claim. InterProScan pending.", ha="center", fontsize=8, color="#444444")
    save(fig, figd, "Figure_8C_alignment_evidence_stack")


def fig_supplement(base, dirs, master, order, recon, recommended, claims):
    cons, figd = dirs["conservation"], dirs["figures"]
    diag = {r["species"].lower() for r in
            M.read_tsv(dirs["review_diagnostics"] / "fgfr2_msa_review_case_diagnostics.tsv")}
    agree = {iso: load_agreement(cons, iso) for iso in ("IIIb", "IIIc")}
    all_sps = species_in_order([k[0] for iso in agree for k in agree[iso]], order)
    review_sps = [s for s in all_sps
                  if s in diag or is_review(recon.get((s, None)), recommended.get(s, ""),
                                            M.species_claim(claims, s))]
    # always include known-risk / close-primate review cases explicitly
    must = [s for s in ("canis_lupus_familiaris", "gorilla_gorilla_gorilla", "pongo_abelii")
            if s in all_sps]
    reps = [s for s in ("homo_sapiens", "mus_musculus", "gallus_gallus", "danio_rerio") if s in all_sps]
    panel_sps = list(dict.fromkeys(must + review_sps + reps))[:14] or all_sps[:8]
    n = len(panel_sps)
    fig, axes = plt.subplots(n, 1, figsize=(11, max(4, n * 1.0)), squeeze=False)
    for i, sp in enumerate(panel_sps):
        ax = axes[i][0]
        offset = 0
        for iso in ("IIIb", "IIIc"):
            ag = agree[iso]
            idxs = sorted([k[1] for k in ag if k[0] == sp])
            for x in idxs:
                r = ag[(sp, x)]
                color = AGREE_COLORS.get(r["agreement_class"], "#FFFFFF")
                ax.add_patch(Rectangle((offset + x - 0.5, 0), 1, 1, facecolor=color,
                                       edgecolor="white", linewidth=0.2))
            mx = max(idxs, default=0)
            ax.text(offset + mx / 2 if mx else offset, 1.15, iso, ha="center", va="bottom",
                    fontsize=8, color=M.C_IIIB if iso == "IIIb" else M.C_IIIC, fontweight="bold")
            offset += (mx + 4)
        rrow = recon.get((sp, None)) or {}
        claim = M.species_claim(claims, sp)
        cr = claims.get((sp, "IIIb"), {}) or claims.get((sp, "IIIc"), {})
        dec = (cr.get("maximal_rescue_decision", "")
               or claims.get((sp, "IIIc"), {}).get("maximal_rescue_decision", ""))
        reason = (f"{rrow.get('label_consistency_status','') or 'main'} | "
                  f"local+external rescue: {dec or 'n/a'} | {claim or 'primary'}")
        ax.text(offset + 1, 0.5, reason, va="center", fontsize=6.5,
                color=C_FAIL if is_excluded(claim) else "#666666")
        ax.set_xlim(-1, offset + 24)
        ax.set_ylim(-0.2, 1.5)
        ax.set_yticks([0.5])
        ax.set_yticklabels([master.get(sp, {}).get("display_species_name", sp)], fontsize=7)
        ax.set_xticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    handles = [Patch(facecolor=c, label=k.replace("_", " ")) for k, c in AGREE_COLORS.items()
               if k != "unmapped_review"]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Supplement — Per-species cassette difference panels", fontsize=12,
                 fontweight="bold", y=0.995)
    fig.text(0.5, -0.01, "Right-hand text: label-consistency | rescue status | final claim status. "
             "Includes known-risk and review species (Canis, Gorilla, Pongo) with explanation.",
             ha="center", fontsize=7.5, color="#444444")
    save(fig, figd, "Supplement_Figure_per_species_cassette_difference_panels")


def main() -> int:
    ap = argparse.ArgumentParser(description="Reference-guided MSA figures (Parts G,H).")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)

    # ---- validation gates: label gate + general rescue gate + maximal rescue gate (Part H) ----
    ok, msgs = M.label_gate(base)
    ok2, msgs2 = M.general_rescue_gate(base)
    ok3, msgs3 = M.maximal_rescue_gate(base)
    if not ok or not ok2 or not ok3:
        print("[FAIL] validation gate failed; no reference figures generated:", file=sys.stderr)
        for m in (msgs + msgs2 + msgs3):
            print("   - " + m, file=sys.stderr)
        return 2

    master = {r["species"].lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}
    recommended = {sp: r.get("recommended_use", "") for sp, r in master.items()}
    order = phylo_order(base, master)
    rec_rows = M.read_tsv(dirs["maps"] / "fgfr2_exon_type_label_reconciliation.tsv")
    recon = {}
    for r in rec_rows:
        recon.setdefault((r["species"].lower(), None), r)  # any one row per species for status
    claims = M.load_claim_status(base)

    fig_6c(base, dirs, master, order, recon, recommended, claims)
    fig_6d(base, dirs, master, order, recon, recommended, claims)
    fig_7c(base, dirs)
    fig_8c(base, dirs, master, order)
    fig_supplement(base, dirs, master, order, recon, recommended, claims)

    # Part E: write the authoritative cross-table consistency gate now that ALL final figure
    # tables (figure6, figure8, figure6C) exist; fail loudly if the corrected state is inconsistent.
    okc, msgsc = M.post_rescue_consistency_gate(base)
    if not okc:
        print("[FAIL] post-rescue cross-table consistency gate FAILED after figures:", file=sys.stderr)
        for m in msgsc:
            print("   - " + m, file=sys.stderr)
        return 2
    print("[OK] reference-guided figures 6C, 6D, 7C, 8C + supplement written (SVG/PDF/PNG); "
          "cross-table consistency gate PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
