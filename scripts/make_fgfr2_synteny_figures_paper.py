#!/usr/bin/env python3
"""
make_fgfr2_synteny_figures_paper.py  (synteny paper-level figures, Parts E-H)

Redesigned, paper-level FGFR2 local-synteny figures using the shared style module and the v2
neighbor-identity resolution (LOC resolver v2):

  Figure 9A  — FGFR2 local synteny map (equal-spacing, representative species, orthology ribbons,
               compact rescue/claim badges).            ...local_synteny_map_paper.{svg,pdf,png}
  Supplement — FGFR2 local synteny, TRUE genomic scale, all species (coordinate-level evidence).
  Figure 9B  — FGFR2 neighbor conservation matrix (clean, all species, 10 reference anchors).
  Figure 9C  — FGFR2 synteny rescue-case cards (narrative badges + one-line interpretation).

Synteny validates the FGFR2 locus / orthology context only; it never assigns or relabels IIIb/IIIc.
BLAST/RBH-inferred names are shown as probable (italic + "?"); unresolved LOCs are shown as "LOC?"
and never forced into a false ortholog name.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fgfr2_msa_common as M  # noqa: E402
from shared_gene_analysis.strand import same_strand, strand_sign  # noqa: E402
import fgfr2_plot_style as S  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

MAIN_SPECIES = ["homo_sapiens", "mus_musculus", "pan_troglodytes", "gorilla_gorilla_gorilla",
                "canis_lupus_familiaris", "bos_taurus", "gallus_gallus", "xenopus_tropicalis",
                "danio_rerio"]
RESOLVED_REF = {"ortholog_supported", "curated_or_symbol_supported", "rbh_supported_neighbor_ortholog"}
PROBABLE = {"blast_supported_probable_neighbor", "high_confidence_one_way_blast",
            "probable_ortholog_supported"}
UNRESOLVED = {"unresolved_LOC", "raw_id_only", "unmapped_neighbor", "raw_id_only_no_sequence"}
NONREF = {"blast_supported_non_reference_gene", "symbol_supported_only"}
AMBIG = {"ambiguous_paralog_family", "ambiguous_neighbor_identity"}

MATRIX2_COLOR = {
    "resolved_reference_neighbor_same_order": S.STATUS_COLOR["confirmed_same_order"],
    "resolved_reference_neighbor_reordered": S.STATUS_COLOR["confirmed_reordered"],
    "probable_reference_neighbor_by_rbh": S.STATUS_COLOR["probable_rbh"],
    "probable_reference_neighbor_by_blast": S.STATUS_COLOR["probable_blast"],
    "ambiguous_paralog_family": S.STATUS_COLOR["ambiguous_paralog_family"],
    "missing_or_unmapped": S.STATUS_COLOR["unresolved_missing"],
    "scaffold_unavailable": S.STATUS_COLOR["scaffold_unavailable"],
}
MATRIX2_LABEL = {
    "resolved_reference_neighbor_same_order": "confirmed (same order)",
    "resolved_reference_neighbor_reordered": "confirmed (reordered)",
    "probable_reference_neighbor_by_rbh": "probable (RBH)",
    "probable_reference_neighbor_by_blast": "probable (BLAST)",
    "ambiguous_paralog_family": "ambiguous paralog",
    "missing_or_unmapped": "missing / unmapped",
    "scaffold_unavailable": "scaffold/source N/A",
}


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load(base, dirs):
    syn = dirs["synteny"]
    master = {(r["species"] or "").lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}
    idv2 = {((r["species"] or "").lower(), r["neighbor_gene_id"]): r for r in
            M.read_tsv(syn / "fgfr2_neighbor_identity_resolution_v2.tsv")}
    n5, n10 = {}, {}
    for r in M.read_tsv(syn / "fgfr2_local_gene_neighborhood_5neighbors.tsv"):
        n5.setdefault((r["species"] or "").lower(), []).append(r)
    for r in M.read_tsv(syn / "fgfr2_local_gene_neighborhood_10neighbors_supplement.tsv"):
        n10.setdefault((r["species"] or "").lower(), []).append(r)
    valid = {(r["species"] or "").lower(): r for r in
             M.read_tsv(syn / "fgfr2_5neighbor_synteny_validation.tsv")}
    matrix = M.read_tsv(syn / "fgfr2_5neighbor_conservation_matrix.tsv")
    href = M.read_tsv(syn / "human_fgfr2_10neighbor_reference.tsv")
    truth = {}
    for r in M.read_tsv(dirs["maps"] / "fgfr2_post_rescue_final_truth_table.tsv"):
        truth.setdefault((r["species"] or "").lower(), []).append(r)
    return master, idv2, n5, n10, valid, matrix, href, truth


def ref_symbols(href) -> List[str]:
    return [r["human_gene_symbol"] for r in
            sorted(href, key=lambda r: (r["human_neighbor_side"],
                                        M.to_int(r["human_neighbor_rank"], 0)))]


def disp_for(idr, raw) -> Tuple[str, str, str]:
    """Return (label, style_kind, group) for a neighbor. style_kind in {ref,probable,nonref,ambig,unresolved,fgfr2}."""
    status = (idr or {}).get("identity_resolution_status", "")
    method = (idr or {}).get("identity_resolution_method", "")
    sym = (idr or {}).get("normalized_neighbor_symbol", "") or raw
    group = ((idr or {}).get("normalized_neighbor_orthology_group", "") or "").upper()
    if status in UNRESOLVED or not sym:
        return "LOC?", "unresolved", ""
    if status in AMBIG:
        return (sym + "?" if sym and not sym.upper().startswith("LOC") else "paralog?"), "ambig", ""
    suf, _ = S.gene_label_style(method)
    if status in RESOLVED_REF:
        return sym, "ref", group
    if status in PROBABLE:
        return sym + "?", "probable", group
    if status in NONREF:
        return sym + (suf or ""), "nonref", group
    return sym, "ref", group


def neighbor_slots(neigh, idv2, species, max_n, equal=True):
    sp = species.lower()
    by = {(r["neighbor_side"], M.to_int(r["neighbor_rank"], 0)): r for r in neigh
          if r.get("neighbor_rank")}
    fr = neigh[0] if neigh else {}
    fstart = M.to_int(fr.get("fgfr2_start"), 0)
    fend = M.to_int(fr.get("fgfr2_end"), 0)
    fmid = (fstart + fend) / 2 if fstart else 0
    fstrand = fr.get("fgfr2_strand", "+")
    slots = [{"x": 0, "side": "fgfr2", "rank": 0, "label": "FGFR2", "kind": "fgfr2", "group": "FGFR2",
              "strand_rel": 1, "is_fgfr2": True, "offset": 0.0}]
    for side, sign in (("upstream", -1), ("downstream", 1)):
        for rank in range(1, max_n + 1):
            r = by.get((side, rank))
            x = sign * rank
            if not r:
                slots.append({"x": x, "side": side, "rank": rank, "label": "", "kind": "empty",
                              "group": "", "strand_rel": 1, "is_fgfr2": False, "offset": None})
                continue
            idr = idv2.get((sp, r["neighbor_gene_id"]), {})
            label, kind, group = disp_for(idr, r.get("neighbor_symbol_raw", ""))
            nstrand = r.get("neighbor_strand", "+")
            rel = 1 if same_strand(nstrand, fstrand) else -1
            mid = (M.to_int(r.get("neighbor_start"), 0) + M.to_int(r.get("neighbor_end"), 0)) / 2
            off = (mid - fmid) * strand_sign(fstrand)
            slots.append({"x": x, "side": side, "rank": rank, "label": label, "kind": kind,
                          "group": group, "strand_rel": rel, "is_fgfr2": False, "offset": off,
                          "raw": r.get("neighbor_symbol_raw", "")})
    return slots, fstrand


def dataset_species(master):
    """Every species this dataset actually analysed, in the thesis reading order first.

    ``MAIN_SPECIES`` is the nine-species reading order of the 30-species thesis panel, and
    using it as a *filter* meant a dataset containing any other species drew that species
    nowhere. On the two-species human/cat run that left a single human row presented as a
    cross-species comparison. The list is therefore an ordering preference, not a membership
    test: species it names come first, and every other species the run analysed follows.
    """
    named = [s for s in MAIN_SPECIES if s in master]
    return named + [s for s in master if s not in set(named)]


def kind_color(kind, group, cmap):
    if kind == "fgfr2":
        return S.FGFR2_COLOR
    if kind in ("ref", "probable") and group in cmap:
        return cmap[group]
    if kind == "nonref":
        return S.NONREF_COLOR
    if kind == "ambig":
        return S.STATUS_COLOR["ambiguous_paralog_family"]
    return S.UNRESOLVED_COLOR


# ---------------------------------------------------------------------------
# Figure 9A — equal-spacing paper map
# ---------------------------------------------------------------------------
def fig_9a(base, dirs, data):
    master, idv2, n5, _n10, valid, _matrix, href, truth = data
    refs = [s.upper() for s in ref_symbols(href)]
    cmap = S.build_color_map(refs)
    sps = dataset_species(master)
    sps.sort(key=lambda s: S.taxon_sort_key(master[s].get("taxon_group", ""),
                                             master[s].get("phylo_order")))
    layouts = {sp: neighbor_slots(n5.get(sp, []), idv2, sp, 5)[0] for sp in sps}

    fig, ax = S.plt.subplots(figsize=(12.5, max(5, len(sps) * 0.86)))
    S.apply_rcparams()
    table_rows = []
    ys = {sp: len(sps) - i for i, sp in enumerate(sps)}
    gh = 0.34  # half gene height
    # ribbons between consecutive rows for shared reference orthology groups
    for a, b in zip(sps, sps[1:]):
        la = {s["group"]: s["x"] for s in layouts[a] if s["kind"] in ("ref", "probable") and s["group"]}
        lb = {s["group"]: s["x"] for s in layouts[b] if s["kind"] in ("ref", "probable") and s["group"]}
        for g in set(la) & set(lb):
            S.ribbon(ax, la[g], ys[a] - gh, lb[g], ys[b] + gh, cmap.get(g, S.MUTED), alpha=0.22)
    for sp in sps:
        y = ys[sp]
        for s in layouts[sp]:
            if s["kind"] == "empty":
                ax.plot([s["x"] - 0.18, s["x"] + 0.18], [y, y], color="#D5D8DC", lw=1.0, ls=(0, (1, 1)),
                        zorder=2)
                continue
            col = kind_color(s["kind"], s["group"], cmap)
            S.gene_arrow(ax, s["x"], y, 0.74, 0.5, s["strand_rel"], col,
                         alpha=1.0 if s["kind"] in ("fgfr2", "ref") else
                         (0.92 if s["kind"] == "probable" else 0.7))
            italic = s["kind"] in ("probable", "ambig")
            color = S.INK if s["kind"] in ("fgfr2", "ref", "probable", "nonref") else S.MUTED
            ax.text(s["x"], y - gh - 0.12, s["label"], ha="center", va="top",
                    fontsize=S.FONT["gene"], fontstyle="italic" if italic else "normal",
                    color=color, fontweight="bold" if s["is_fgfr2"] else "normal")
            table_rows.append({"species": sp, "slot": s["x"], "side": s["side"], "rank": s["rank"],
                               "label": s["label"], "kind": s["kind"], "orthology_group": s["group"]})
        # badges
        v = valid.get(sp, {})
        ts = truth.get(sp, [])
        claim = ";".join(sorted({t.get("final_claim_status_after_rescue", "") for t in ts} - {""}))
        rescued = any((t.get("rescue_decision") or "").startswith("rescued") for t in ts)
        bx = 5.55
        ax.add_patch(Rectangle((bx, y - 0.16), 0.26, 0.32, facecolor=S.SYNTENY_CLASS_COLOR.get(
            v.get("synteny_validation_class", ""), "#CCCCCC"), edgecolor="none", clip_on=False))
        ckind = "primary" if claim.startswith("primary") else ("excluded" if "excluded" in claim
                                                               else "supplement")
        ax.add_patch(Rectangle((bx + 0.32, y - 0.16), 0.26, 0.32,
                               facecolor=S.CLAIM_COLOR.get(ckind, "#CCCCCC"), edgecolor="none",
                               clip_on=False))
        if rescued:
            ax.text(bx + 0.78, y, "rescued", fontsize=S.FONT["badge"], va="center",
                    color=S.RESCUE_OUTLINE if hasattr(S, "RESCUE_OUTLINE") else S.BADGE_COLOR["ok"],
                    fontweight="bold")
        disp = master.get(sp, {}).get("display_species_name", sp)
        ax.text(-5.5, y, disp, ha="right", va="center", fontsize=S.FONT["label"], fontstyle="italic")
    # subtle taxon group labels
    _taxon_bands(ax, sps, ys, master, x=-7.0)
    ax.set_xlim(-7.4, 7.2)
    ax.set_ylim(0.3, len(sps) + 0.9)
    ax.axis("off")
    S.title(ax, "FGFR2 local gene-neighborhood synteny",
            "Equal-spacing order map · 5 protein-coding neighbors each side · locus/orthology context only")
    handles = [S.legend_patch(S.FGFR2_COLOR, "FGFR2 anchor"),
               S.legend_patch(cmap.get(refs[0], S.PALETTE[0]) if refs else S.PALETTE[0],
                              "conserved reference neighbor"),
               S.legend_patch(S.STATUS_COLOR["probable_blast"], "probable (italic + ?)"),
               S.legend_patch(S.NONREF_COLOR, "non-reference local gene"),
               S.legend_patch(S.UNRESOLVED_COLOR, "unresolved (LOC?)"),
               S.legend_line(S.MUTED, "orthology ribbon", lw=1.1)]
    S.compact_legend(ax, handles, ncol=6, bbox=(0.5, -0.06))
    M.write_tsv(dirs["tables"] / "figure9A_fgfr2_local_synteny_map_paper.tsv", table_rows,
                list(table_rows[0].keys()) if table_rows else ["species"])
    S.savefig(fig, dirs["figures"], "Figure_9A_FGFR2_local_synteny_5neighbor_paper")


def _taxon_bands(ax, sps, ys, master, x):
    groups = {}
    for sp in sps:
        groups.setdefault(master[sp].get("taxon_group", ""), []).append(ys[sp])
    for g, yy in groups.items():
        if not g:
            continue
        ax.text(x, sum(yy) / len(yy), g, rotation=90, ha="center", va="center",
                fontsize=S.FONT["small"], color=S.MUTED, fontweight="bold")


# ---------------------------------------------------------------------------
# Supplement — true genomic scale, all species
# ---------------------------------------------------------------------------
def fig_true_scale(base, dirs, data):
    master, idv2, _n5, n10, _valid, _matrix, href, _truth = data
    refs = [s.upper() for s in ref_symbols(href)]
    cmap = S.build_color_map(refs)
    sps = sorted(master.keys(), key=lambda s: S.taxon_sort_key(master[s].get("taxon_group", ""),
                                                               master[s].get("phylo_order")))
    WIN = 2_000_000.0  # +/- 2 Mb window
    fig, ax = S.plt.subplots(figsize=(16, max(8, len(sps) * 0.6)))
    S.apply_rcparams()
    table_rows = []
    for i, sp in enumerate(sps):
        y = len(sps) - i
        slots, _fstrand = neighbor_slots(n10.get(sp, []), idv2, sp, 10)
        ax.plot([-1, 1], [y, y], color="#EAECEE", lw=0.8, zorder=1)
        for s in slots:
            if s["kind"] == "empty" or s["offset"] is None:
                continue
            off = s["offset"]
            if abs(off) > WIN:
                xb = 1.02 if off > 0 else -1.02
                ax.text(xb, y, "//", fontsize=7, va="center", ha="center", color=S.MUTED)
                continue
            x = off / WIN
            col = kind_color(s["kind"], s["group"], cmap)
            S.gene_arrow(ax, x, y, 0.045, 0.42, s["strand_rel"], col,
                         alpha=1.0 if s["kind"] in ("fgfr2", "ref") else 0.8, lw=0.4)
            if s["is_fgfr2"] or s["kind"] == "ref":
                ax.text(x, y + 0.28, s["label"], ha="center", va="bottom", fontsize=5.4,
                        color=S.INK, rotation=0, fontweight="bold" if s["is_fgfr2"] else "normal")
            table_rows.append({"species": sp, "label": s["label"], "kind": s["kind"],
                               "offset_bp": int(s["offset"]) if s["offset"] is not None else "",
                               "orthology_group": s["group"]})
        disp = master[sp].get("display_species_name", sp)
        ax.text(-1.06, y, disp, ha="right", va="center", fontsize=6.6, fontstyle="italic")
    ax.axhline
    ax.set_xlim(-1.35, 1.12)
    ax.set_ylim(0.3, len(sps) + 1.0)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xticklabels(["-2 Mb", "-1 Mb", "FGFR2", "+1 Mb", "+2 Mb"], fontsize=S.FONT["tick"])
    ax.set_yticks([])
    for sp_ in ("top", "right", "left"):
        ax.spines[sp_].set_visible(False)
    S.title(ax, "FGFR2 local synteny — true genomic scale (all species)",
            "Neighbors placed by real genomic offset from FGFR2 (±2 Mb); // = neighbor beyond window")
    handles = [S.legend_patch(S.FGFR2_COLOR, "FGFR2"),
               S.legend_patch(cmap.get(refs[0], S.PALETTE[0]), "reference neighbor"),
               S.legend_patch(S.NONREF_COLOR, "non-reference gene"),
               S.legend_patch(S.UNRESOLVED_COLOR, "unresolved LOC")]
    S.compact_legend(ax, handles, ncol=4, bbox=(0.5, -0.05))
    M.write_tsv(dirs["tables"] / "supplement_fgfr2_local_synteny_true_scale_all_species.tsv",
                table_rows, list(table_rows[0].keys()) if table_rows else ["species"])
    S.savefig(fig, dirs["figures"], "Supplement_Figure_FGFR2_local_synteny_10neighbor_all_species")


# ---------------------------------------------------------------------------
# Figure 9B — clean conservation matrix
# ---------------------------------------------------------------------------
def fig_9b(base, dirs, data):
    master, _idv2, _n5, _n10, valid, matrix, href, truth = data
    cols = ref_symbols(href)
    rows = sorted(matrix, key=lambda r: S.taxon_sort_key(
        master.get((r["species"] or "").lower(), {}).get("taxon_group", ""),
        master.get((r["species"] or "").lower(), {}).get("phylo_order")))
    n, m = len(rows), len(cols)
    fig, ax = S.plt.subplots(figsize=(max(11, m * 0.62 + 6), max(7.5, n * 0.32)))
    S.apply_rcparams()
    for yi, r in enumerate(rows):
        sp = (r["species"] or "").lower()
        for xi, c in enumerate(cols):
            st = r.get(c, "missing_or_unmapped")
            ax.add_patch(Rectangle((xi, yi), 0.92, 0.92, facecolor=MATRIX2_COLOR.get(
                st, S.STATUS_COLOR["unresolved_missing"]), edgecolor="none"))
        v = valid.get(sp, {})
        claim = ";".join(sorted({t.get("final_claim_status_after_rescue", "")
                                 for t in truth.get(sp, [])} - {""}))
        ckind = "primary" if claim.startswith("primary") else ("excluded" if "excluded" in claim
                                                               else "supplement")
        q = M.to_float(v.get("neighbor_label_quality_score"), 0.0)
        side = [S.SYNTENY_CLASS_COLOR.get(v.get("synteny_validation_class", ""), "#CCCCCC"),
                S.plt.cm.Blues(0.25 + 0.6 * q), S.CLAIM_COLOR.get(ckind, "#CCCCCC")]
        for bi, col in enumerate(side):
            ax.add_patch(Rectangle((m + 0.25 + bi * 0.55, yi), 0.5, 0.92, facecolor=col,
                                   edgecolor="white", lw=0.4, clip_on=False))
    ax.set_xlim(0, m + 0.25 + 3 * 0.55 + 0.2)
    ax.set_ylim(0, n)
    ax.invert_yaxis()
    ax.set_xticks([i + 0.46 for i in range(m)])
    ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=S.FONT["small"])
    ax.set_yticks([i + 0.46 for i in range(n)])
    ax.set_yticklabels([master.get((r["species"] or "").lower(), {})
                        .get("display_species_name", r["species"]) for r in rows],
                       fontsize=S.FONT["small"], fontstyle="italic")
    for bi, lab in enumerate(["synteny", "label quality", "claim"]):
        ax.text(m + 0.5 + bi * 0.55, -0.5, lab, rotation=45, ha="left", va="bottom",
                fontsize=S.FONT["small"], color=S.MUTED)
    for sp_ in ("top", "right", "left", "bottom"):
        ax.spines[sp_].set_visible(False)
    ax.tick_params(length=0)
    S.title(ax, "FGFR2 neighbor conservation across species",
            "Rows = species (taxon order) · columns = human FGFR2 reference neighbors")
    handles = [S.legend_patch(c, MATRIX2_LABEL[k]) for k, c in MATRIX2_COLOR.items()]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.005, 1.0),
              fontsize=S.FONT["legend"], frameon=False, handlelength=1.1)
    M.write_tsv(dirs["tables"] / "figure9B_fgfr2_neighbor_conservation_matrix_paper.tsv", rows,
                list(rows[0].keys()) if rows else ["species"])
    S.savefig(fig, dirs["figures"], "Figure_9B_FGFR2_5neighbor_conservation_matrix_paper")


# ---------------------------------------------------------------------------
# Figure 9C — narrative rescue cards
# ---------------------------------------------------------------------------
def _interpret(sp_disp, v, claim, rescued, unres):
    cls = v.get("synteny_validation_class", "")
    strong = cls in ("synteny_strong", "synteny_supported_with_minor_rearrangement")
    if rescued and claim.startswith("primary") and strong:
        return f"{sp_disp} FGFR2 is retained as a primary locus after sequence rescue and local synteny support."
    if rescued and "supplement" in claim:
        return (f"{sp_disp} is sequence-rescued in part; the unresolved isoform stays supplement-only "
                f"because no source-compatible validated candidate was found.")
    if "supplement" in claim or "excluded" in claim:
        return f"{sp_disp} remains supplement/review; local synteny ({cls.replace('synteny_','')}) is reported but does not relabel IIIb/IIIc."
    return f"{sp_disp} sits in the conserved FGFR2 neighborhood ({cls.replace('synteny_','')})."


def fig_9c(base, dirs, data):
    master, idv2, n5, _n10, valid, _matrix, href, truth = data
    refs = [s.upper() for s in ref_symbols(href)]
    cmap = S.build_color_map(refs)
    cases = ["gorilla_gorilla_gorilla", "canis_lupus_familiaris", "pongo_abelii"]
    for sp, ts in truth.items():
        claim = {t.get("final_claim_status_after_rescue", "") for t in ts}
        if sp not in cases and any(c in ("supplement_review", "excluded_from_primary_claim")
                                   for c in claim):
            cases.append(sp)
    cases = [s for s in cases if s in master]
    n = len(cases)
    if n == 0:
        # A dataset where no species needs review has no review panel. Saying so is the
        # honest outcome; asking matplotlib for a nought-row figure is a crash.
        print("[skip] Figure 9C: no review-case species in this dataset")
        return
    fig, axes = S.plt.subplots(n, 1, figsize=(11.5, max(4, n * 1.85)))
    S.apply_rcparams()
    axes = axes if hasattr(axes, "__len__") else [axes]
    table_rows = []
    for ax, sp in zip(axes, cases):
        slots, _fstrand = neighbor_slots(n5.get(sp, []), idv2, sp, 5)
        for s in slots:
            if s["kind"] == "empty":
                continue
            col = kind_color(s["kind"], s["group"], cmap)
            S.gene_arrow(ax, s["x"], 1.25, 0.74, 0.5, s["strand_rel"], col,
                         alpha=1.0 if s["kind"] in ("fgfr2", "ref") else 0.85)
            ax.text(s["x"], 1.25 - 0.42, s["label"], ha="center", va="top", fontsize=S.FONT["gene"],
                    fontstyle="italic" if s["kind"] in ("probable", "ambig") else "normal",
                    color=S.INK if s["kind"] != "unresolved" else S.MUTED,
                    fontweight="bold" if s["is_fgfr2"] else "normal")
        v = valid.get(sp, {})
        ts = truth.get(sp, [])
        claim = ";".join(sorted({t.get("final_claim_status_after_rescue", "") for t in ts} - {""}))
        rescued = any((t.get("rescue_decision") or "").startswith("rescued") for t in ts)
        unres = ";".join(sorted({t.get("unresolved_reason_if_any", "") for t in ts} - {"", "-"}))
        disp = master[sp].get("display_species_name", sp)
        ax.set_xlim(-6.0, 6.4)
        ax.set_ylim(-0.2, 2.2)
        ax.axis("off")
        ax.text(-5.9, 1.95, disp, fontsize=10, fontweight="bold", fontstyle="italic")
        # badges
        bx, by = -5.9, 0.55
        _badge(ax, bx, by, "sequence rescue", "ok" if rescued else "neutral")
        _badge(ax, bx + 2.05, by, f"synteny: {v.get('synteny_validation_class','').replace('synteny_','')}",
               "ok" if v.get("synteny_validation_class") in ("synteny_strong",
               "synteny_supported_with_minor_rearrangement") else "partial")
        _badge(ax, bx + 4.6, by, "primary" if claim.startswith("primary") else "supplement",
               "ok" if claim.startswith("primary") else "partial")
        if unres:
            _badge(ax, bx, by - 0.5, "unresolved: " + unres[:40], "fail")
        ax.text(-5.9, -0.05, _interpret(disp, v, claim, rescued, unres), fontsize=S.FONT["label"],
                va="top", color=S.INK, style="italic")
        table_rows.append({"species": sp, "synteny_validation_class": v.get("synteny_validation_class", ""),
                           "neighbor_label_quality_score": v.get("neighbor_label_quality_score", ""),
                           "rescued": "true" if rescued else "false",
                           "final_claim_status_after_rescue": claim,
                           "unresolved_reason": unres or "-",
                           "interpretation": _interpret(disp, v, claim, rescued, unres)})
    fig.suptitle("FGFR2 synteny rescue-case cards", fontsize=S.FONT["title"], fontweight="bold",
                 x=0.07, ha="left", y=0.995)
    M.write_tsv(dirs["tables"] / "figure9C_fgfr2_synteny_rescue_case_cards.tsv", table_rows,
                list(table_rows[0].keys()) if table_rows else ["species"])
    S.savefig(fig, dirs["figures"], "Figure_9C_FGFR2_synteny_review_cases_paper")


def _badge(ax, x, y, label, kind):
    col = S.BADGE_COLOR.get(kind, S.BADGE_COLOR["neutral"])
    ax.add_patch(FancyBboxPatch((x, y - 0.14), 1.9, 0.28, boxstyle="round,pad=0.02,rounding_size=0.08",
                                facecolor=col, edgecolor="none", alpha=0.16, zorder=3))
    ax.text(x + 0.1, y, label, fontsize=S.FONT["badge"], va="center", ha="left", color=col,
            fontweight="bold", zorder=4)


def main() -> int:
    ap = argparse.ArgumentParser(description="FGFR2 paper-level synteny figures (Parts E-H).")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    ok, msgs = M.synteny_gate(base)
    if not ok:
        print("[FAIL] synteny gate failed; no paper figures plotted:", file=sys.stderr)
        for mm in msgs:
            print("   - " + mm, file=sys.stderr)
        return 2
    S.apply_rcparams()
    data = load(base, dirs)
    fig_9a(base, dirs, data)
    fig_true_scale(base, dirs, data)
    fig_9b(base, dirs, data)
    fig_9c(base, dirs, data)
    print("[OK] paper synteny figures written (Figure 9A/9B/9C + true-scale supplement; svg/pdf/png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
