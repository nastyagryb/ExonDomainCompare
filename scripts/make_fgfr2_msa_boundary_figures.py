#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl_msa_"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


C_IIIB, C_IIIC, C_REVIEW = M.C_IIIB, M.C_IIIC, M.C_REVIEW
ROBUST_COLORS = {
    "robust_boundary": "#009E73",
    "supported_boundary_with_minor_flags": "#56B4E9",
    "review_boundary": "#E69F00",
    "unresolved_or_annotation_dependent_boundary": "#D55E00",
    "no_data": "#CCCCCC",
}
TAXON_BAND = ["#F2F2F2", "#E6EEF5"]


def savefig(fig, fig_dir: Path, stem: str) -> None:
    for ext in ("svg", "pdf", "png"):
        fig.savefig(fig_dir / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def load_common(base: Path, dirs: Dict[str, Path]):
    master = {r["species"].lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}
    proj = M.read_tsv(dirs["maps"] / "fgfr2_exon_boundary_msa_projection.tsv")
    scores = {(r["species"].lower(), r["isoform"]): r for r in
              M.read_tsv(dirs["robustness"] / "fgfr2_boundary_robustness_scores.tsv")}
    integ = M.read_tsv(dirs["protein_integrity"] / "fgfr2_pre_interpro_protein_integrity_qc.tsv")
    diag = M.read_tsv(dirs["review_diagnostics"] / "fgfr2_msa_review_case_diagnostics.tsv")
    disc = M.read_tsv(dirs["conservation"] / "fgfr2_IIIb_IIIc_discriminating_positions_main_only.tsv")
    full_cons = M.read_tsv(dirs["conservation"] / "fgfr2_full_length_msa_column_conservation.tsv")
    return master, proj, scores, integ, diag, disc, full_cons




SYNC_COLS = ["rescue_decision", "final_label_source", "final_claim_status_after_rescue",
             "recommended_use_pre_rescue", "recommended_use_post_rescue", "rescue_evidence_summary",
             "unresolved_reason_if_any"]


def build_fig6_table(base, dirs, master, proj, scores, truth) -> List[Dict[str, object]]:
    rows = []
    for p in proj:
        sp, iso = p["species"], p["isoform"]
        mr = master.get(sp.lower(), {})
        sc = scores.get((sp.lower(), iso), {})
        tr = truth.get((sp.lower(), iso), {})
        row = {
            "species": sp, "isoform": iso, "phylo_order": mr.get("phylo_order", ""),
            "taxon_group": mr.get("taxon_group_display", mr.get("taxon_group", "")),
            "alignment_col_start": p.get("full_length_msa_start_col", ""),
            "alignment_col_end": p.get("full_length_msa_end_col", ""),
            "boundary_left_col": p.get("full_length_msa_start_col", ""),
            "boundary_right_col": p.get("full_length_msa_end_col", ""),
            "boundary_projection_status": p.get("boundary_projection_status", ""),
            "boundary_robustness_class": sc.get("boundary_robustness_class", "no_data"),
            "conservation_score_region": sc.get("cassette_conservation_score", ""),
            "gap_fraction_region": p.get("internal_region_gap_fraction", ""),
            "review_reason_short": mr.get("review_reason_short", ""),
        }
        # post-rescue truth columns (single source of truth)
        for c in SYNC_COLS:
            row[c] = tr.get(c, p.get(c, ""))
        rows.append(row)
    cols = ["species", "isoform", "phylo_order", "taxon_group", "alignment_col_start",
            "alignment_col_end", "boundary_left_col", "boundary_right_col",
            "boundary_projection_status", "boundary_robustness_class",
            "conservation_score_region", "gap_fraction_region",
            "review_reason_short"] + SYNC_COLS
    M.write_tsv(dirs["tables"] / "figure6_msa_projected_boundary_map.tsv", rows, cols)
    return rows


def _is_primary(claim: str) -> bool:
    return (claim or "").startswith("primary_claim")


def fig6(dirs, master, fig6_rows, full_cons):
    aln_len = max((M.to_int(c["alignment_col"], 0) for c in full_cons), default=1)
    cons_by_col = {M.to_int(c["alignment_col"]): M.to_float(c["conservation_score"], 0.0)
                   for c in full_cons}
    fig, axes = plt.subplots(1, 2, figsize=(15, 10), sharey=True,
                             gridspec_kw={"wspace": 0.06})
    for ax, iso, col in zip(axes, ("IIIb", "IIIc"), (C_IIIB, C_IIIC)):
        # Part D/G: primary figure shows only post-rescue primary-claim rows
        sub = sorted([r for r in fig6_rows if r["isoform"] == iso
                      and _is_primary(r.get("final_claim_status_after_rescue", ""))],
                     key=lambda r: M.to_int(r["phylo_order"], 999) or 999)
        y = 0
        yticks, ylabels = [], []
        last_tax = None
        band = 0
        for r in sub:
            tax = r["taxon_group"]
            if tax != last_tax:
                band ^= 1
                last_tax = tax
            ax.axhspan(y - 0.5, y + 0.5, color=TAXON_BAND[band], zorder=0)
            s = M.to_int(r["alignment_col_start"])
            e = M.to_int(r["alignment_col_end"])
            is_review = r.get("final_claim_status_after_rescue") == \
                "primary_claim_supported_with_minor_flags"
            if s and e and e >= s:
                ax.add_patch(Rectangle((s, y - 0.32), e - s, 0.64, facecolor=col,
                                       edgecolor="black", linewidth=0.4,
                                       alpha=0.55 if is_review else 0.9, zorder=2))
                ax.plot([s, s], [y - 0.4, y + 0.4], color="black", lw=1.1, zorder=3)
                ax.plot([e, e], [y - 0.4, y + 0.4], color="black", lw=1.1, zorder=3)
            stat = r["boundary_projection_status"]
            if "review" in stat:
                ax.plot(aln_len * 1.01, y, marker=">", color=C_REVIEW, ms=7,
                        clip_on=False, zorder=4)
            # robustness side strip
            rc = ROBUST_COLORS.get(r["boundary_robustness_class"], "#CCCCCC")
            ax.add_patch(Rectangle((-aln_len * 0.035, y - 0.45), aln_len * 0.03, 0.9,
                                   facecolor=rc, edgecolor="none", clip_on=False, zorder=4))
            disp = master.get(r["species"].lower(), {}).get("display_species_name", r["species"])
            yticks.append(y)
            ylabels.append(("• " if is_review else "") + disp)
            y += 1
        # conservation track (top)
        top = y + 0.5
        xs = sorted(cons_by_col)
        cons_vals = [cons_by_col[c] for c in xs]
        ax.plot(xs, [top + 0.8 * v for v in cons_vals], color="#444444", lw=0.7, zorder=2)
        ax.text(aln_len * 0.5, top + 1.05, "per-column conservation", ha="center",
                fontsize=8, color="#444444")
        ax.set_xlim(-aln_len * 0.04, aln_len * 1.03)
        ax.set_ylim(-1, top + 1.4)
        ax.invert_yaxis()
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, fontsize=7)
        ax.set_xlabel("full-length MSA alignment column")
        ax.set_title(f"{iso} cassette projected onto MSA", color=col, fontweight="bold")
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
    leg = [Patch(facecolor=ROBUST_COLORS[k], label=k.replace("_", " ")) for k in
           ("robust_boundary", "supported_boundary_with_minor_flags", "review_boundary",
            "unresolved_or_annotation_dependent_boundary")]
    leg.append(Line2D([0], [0], marker=">", color="w", markerfacecolor=C_REVIEW,
                      label="MSA boundary projection review", markersize=8))
    leg.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#777",
                      label="• primary with minor flags", markersize=7))
    fig.legend(handles=leg, loc="lower center", ncol=3, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Figure 6 — MSA-projected IIIb/IIIc cassette boundaries (primary-claim cassettes)\n"
                 "Exon-defined boundaries (validated/rescued coordinates) projected onto the MSA; "
                 "post-rescue primary claims only; MSA does not relabel IIIb/IIIc.",
                 fontsize=12, y=0.99)
    savefig(fig, dirs["figures"], "Figure_6_MSA_projected_IIIb_IIIc_boundary_map")


def fig7(dirs, disc):
    cols = [M.to_int(d["alignment_col"]) for d in disc]
    iiib = [M.to_float(d["IIIb_major_aa_fraction"], 0.0) for d in disc]
    iiic = [M.to_float(d["IIIc_major_aa_fraction"], 0.0) for d in disc]
    dsc = [M.to_float(d["discriminating_score"], 0.0) for d in disc]
    pclass = [d["position_class"] for d in disc]
    M.write_tsv(dirs["tables"] / "figure7_isoform_discriminating_residues.tsv", disc,
                list(disc[0].keys()) if disc else ["alignment_col"])
    fig, axes = plt.subplots(3, 1, figsize=(13, 8), sharex=True,
                             gridspec_kw={"hspace": 0.18})
    axes[0].bar(cols, iiib, color=C_IIIB, width=1.0)
    axes[0].set_ylabel("IIIb\nconservation", fontsize=9)
    axes[0].set_title("Figure 7 — Isoform-discriminating cassette residues "
                      "(main-analysis species; combined IIIb+IIIc MSA)", fontweight="bold")
    axes[1].bar(cols, iiic, color=C_IIIC, width=1.0)
    axes[1].set_ylabel("IIIc\nconservation", fontsize=9)
    bar_colors = ["#CC79A7" if pc == "isoform_discriminating_conserved" else "#999999"
                  for pc in pclass]
    axes[2].bar(cols, dsc, color=bar_colors, width=1.0)
    axes[2].set_ylabel("discriminating\nscore", fontsize=9)
    axes[2].set_xlabel("combined cassette MSA alignment column")
    for d, c, s in zip(disc, cols, dsc):
        if d["position_class"] == "isoform_discriminating_conserved":
            hb = d.get("human_IIIb_aa_if_available", "")
            hc = d.get("human_IIIc_aa_if_available", "")
            if hb or hc:
                axes[2].annotate(f"{hb}/{hc}", (c, s), textcoords="offset points",
                                 xytext=(0, 3), ha="center", fontsize=6, color="#333333")
    for ax in axes:
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
    fig.legend(handles=[Patch(facecolor="#CC79A7", label="isoform-discriminating conserved"),
                        Patch(facecolor="#999999", label="other positions")],
               loc="lower center", ncol=2, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.03))
    savefig(fig, dirs["figures"], "Figure_7_IIIb_IIIc_isoform_discriminating_residues")


def _fig8_sorted(master, scores):
    return sorted(scores.values(),
                  key=lambda r: (M.to_int(master.get(r["species"].lower(), {}).get("phylo_order"), 999),
                                 r["isoform"]))


def fig8_table_only(dirs, master, scores):
    all_rows = _fig8_sorted(master, scores)
    M.write_tsv(dirs["tables"] / "figure8_boundary_robustness_evidence_stack.tsv", all_rows,
                list(all_rows[0].keys()) if all_rows else ["species"])


def fig8(dirs, master, scores):
    comp_cols = [("component_annotation_score", "coord\nresolution"),
                 ("component_codon_phase_score", "codon/\nboundary"),
                 ("component_protein_qc_score", "protein\nQC"),
                 ("component_msa_projection_score", "MSA\nprojection"),
                 ("component_conservation_score", "conserv./\ngap"),
                 ("component_integrity_score", "protein\nintegrity")]
    # Part D/G: primary evidence stack plots only post-rescue primary-claim rows
    rows = [r for r in _fig8_sorted(master, scores)
            if _is_primary(r.get("final_claim_status_after_rescue", ""))]
    n = len(rows)
    fig, ax = plt.subplots(figsize=(12, max(8, n * 0.32)))
    cmap = plt.get_cmap("RdYlGn")
    for yi, r in enumerate(rows):
        for xi, (key, _) in enumerate(comp_cols):
            v = M.to_float(r.get(key), 0.0)
            ax.add_patch(Rectangle((xi, yi), 1, 1, facecolor=cmap(v), edgecolor="white", lw=0.5))
            ax.text(xi + 0.5, yi + 0.5, f"{v:.2f}", ha="center", va="center", fontsize=6)
        # final class swatch + score
        rc = ROBUST_COLORS.get(r["boundary_robustness_class"], "#CCCCCC")
        ax.add_patch(Rectangle((len(comp_cols), yi), 1, 1, facecolor=rc, edgecolor="white", lw=0.5))
        ax.text(len(comp_cols) + 0.5, yi + 0.5, r.get("boundary_robustness_score", ""),
                ha="center", va="center", fontsize=6, fontweight="bold")
    ax.set_xlim(0, len(comp_cols) + 1)
    ax.set_ylim(0, n)
    ax.invert_yaxis()
    ax.set_xticks([i + 0.5 for i in range(len(comp_cols) + 1)])
    ax.set_xticklabels([c[1] for c in comp_cols] + ["robustness\nclass/score"], fontsize=8)
    ax.xaxis.tick_top()
    disp = [f"{master.get(r['species'].lower(), {}).get('display_species_name', r['species'])} "
            f"[{r['isoform']}]" for r in rows]
    ax.set_yticks([i + 0.5 for i in range(n)])
    ax.set_yticklabels(disp, fontsize=6)
    ax.set_title("Figure 8 — Boundary robustness evidence stack (primary-claim cassettes)\n"
                 "Independent evidence layers per species/isoform (not a single coordinate source)",
                 fontweight="bold", pad=30)
    leg = [Patch(facecolor=ROBUST_COLORS[k], label=k.replace("_", " ")) for k in
           ("robust_boundary", "supported_boundary_with_minor_flags", "review_boundary",
            "unresolved_or_annotation_dependent_boundary")]
    ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, frameon=False)
    savefig(fig, dirs["figures"], "Figure_8_boundary_robustness_evidence_stack")


def supp_review(dirs, master, diag):
    if not diag:
        diag = [{"species": "none", "isoform": "", "final_interpretation": "no review cases"}]
    n = len(diag)
    ncol = 3
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(15, max(4, nrow * 2.4)))
    axes = (axes.flatten() if hasattr(axes, "flatten") else [axes])
    for ax, d in zip(axes, diag):
        ax.axis("off")
        disp = master.get(d["species"].lower(), {}).get("display_species_name", d["species"])
        rc = ROBUST_COLORS.get(d.get("boundary_robustness_class", ""), "#888888")
        ax.add_patch(Rectangle((0.02, 0.85), 0.96, 0.13, transform=ax.transAxes,
                               facecolor="#EFEFEF", edgecolor="none"))
        ax.text(0.05, 0.91, f"{disp} [{d.get('isoform','')}]", transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="center")
        lines = [
            f"recommended_use: {d.get('recommended_use','')}",
            f"native coord: {d.get('native_coordinate_sanity','')}",
            f"normalized slot: {d.get('normalized_slot_sanity','')}",
            f"protein evidence: {d.get('protein_evidence_state','')}",
            f"MSA projection: {d.get('msa_boundary_projection_status','')}",
            f"bw gap frac: {d.get('boundary_window_gap_fraction','')}",
            f"cassette cons: {d.get('cassette_conservation_score','')}",
            f"robustness: {d.get('boundary_robustness_score','')}",
            f"interpretation: {d.get('final_interpretation','')}",
            f"display: {d.get('suggested_display_location','')}",
        ]
        ax.text(0.05, 0.80, "\n".join(lines), transform=ax.transAxes, fontsize=7,
                va="top", family="monospace")
        ax.add_patch(Rectangle((0.9, 0.02), 0.08, 0.1, transform=ax.transAxes, facecolor=rc))
    for ax in axes[len(diag):]:
        ax.axis("off")
    fig.suptitle("Supplement — MSA review-case diagnostics (cases retained and explained, not hidden)",
                 fontsize=12, y=1.0)
    savefig(fig, dirs["figures"], "Supplement_Figure_MSA_review_case_diagnostics")


def supp_integrity(dirs, master, integ, dirs_full_manifest):
    lengths = [(r["species"], r["isoform"], M.to_int(r["sequence_length"], 0),
                r["protein_integrity_status"], r["interpro_ready"]) for r in integ]
    run_man = M.read_tsv(dirs["metadata"] / "msa_run_manifest.tsv")
    full = next((r for r in run_man if r["msa_name"] == "full_length_protein"), {})
    # per-sequence gap fraction from full-length alignment
    aln = [(i, M.clean_alignment_seq(s)) for i, s in
           M.read_fasta(dirs["alignments"] / "fgfr2_full_length_protein_msa.aln.faa")]
    gap_frac = []
    for _sid, s in aln:
        g = sum(1 for c in s if c in M.GAP_CHARS) / len(s) if s else 1.0
        gap_frac.append(g)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].hist([l[2] for l in lengths], bins=20, color="#0072B2", edgecolor="white")
    axes[0].set_title("protein length distribution")
    axes[0].set_xlabel("length (aa)")
    axes[1].hist(gap_frac, bins=20, color="#E69F00", edgecolor="white")
    axes[1].set_title("per-sequence gap fraction (full-length MSA)")
    axes[1].set_xlabel("gap fraction")
    from collections import Counter
    cc = Counter(l[3] for l in lengths)
    axes[2].barh(list(cc.keys()), list(cc.values()), color="#009E73")
    axes[2].set_title("protein integrity status")
    n_ready = sum(1 for l in lengths if l[4] == "true")
    fig.suptitle(f"Supplement — Full-length MSA / protein integrity overview "
                 f"(InterProScan-ready: {n_ready}/{len(lengths)}; aligned length "
                 f"{full.get('aligned_length','?')})", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    savefig(fig, dirs["figures"], "Supplement_Figure_full_length_MSA_protein_integrity")


def main() -> int:
    ap = argparse.ArgumentParser(description="Render MSA boundary figures.")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    # validation gates: final_isoform_label==validated_exon_type + general rescue gate (Part J)
    ok, msgs = M.label_gate(base)
    ok2, msgs2 = M.general_rescue_gate(base)
    ok3, msgs3 = M.maximal_rescue_gate(base)
    if not ok or not ok2 or not ok3:
        print("[FAIL] validation gate failed; no figures generated:", file=sys.stderr)
        for m in (msgs + msgs2 + msgs3):
            print("   - " + m, file=sys.stderr)
        return 2
    master, proj, scores, integ, diag, disc, full_cons = load_common(base, dirs)
    truth = M.load_post_rescue_truth(base)

    # build synchronized figure tables first (single source of truth = post-rescue truth table)
    fig6_rows = build_fig6_table(base, dirs, master, proj, scores, truth)
    fig8_table_only(dirs, master, scores)

    # Part E: cross-table consistency gate must pass before any final figure is plotted
    okc, msgsc = M.post_rescue_consistency_gate(base)
    if not okc:
        print("[FAIL] post-rescue cross-table consistency gate failed; no figures plotted:",
              file=sys.stderr)
        for m in msgsc:
            print("   - " + m, file=sys.stderr)
        return 2

    fig6(dirs, master, fig6_rows, full_cons)
    fig7(dirs, disc)
    fig8(dirs, master, scores)
    supp_review(dirs, master, diag)
    supp_integrity(dirs, master, integ, None)
    print("[OK] MSA figures written (Figure 6, 7, 8 + 2 supplements; svg/pdf/png each)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
