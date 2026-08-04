#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402
from exondomaincompare.presentation import fgfr2_plot_style as S  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl_fw_"))

# Only evidence layers that were actually performed are shown as main columns.
# Not-performed analyses (e.g. MCScanX block synteny) are reported in text, NOT as a
# confusing fully-grey main column.
COLS = [
    ("label_reconciliation", "Label recon."),
    ("rescue_status", "Rescue"),
    ("orthology_status", "Orthology/paralog"),
    ("coordinate_validation_status", "Coordinates"),
    ("MSA_full_length_status", "Full-length MSA"),
    ("boundary_robustness_class", "Cassette MSA"),
    ("reference_agreement", "Ref. agreement"),
    ("synteny_validation_class", "Local synteny"),
    ("pre_interpro_readiness_class", "InterPro ready"),
]

# discrete colour-blind-safe state colours
STATE = {
    "pass": "#1B6CA8", "primary": "#1B6CA8", "robust": "#1B6CA8", "strong": "#1B6CA8",
    "minor": "#5B9BD5", "supported": "#5B9BD5", "probable": "#7B6FB0",
    "review": "#E69F00", "supplement": "#E69F00", "partial": "#7B6FB0",
    "fail": "#D55E00", "excluded": "#D55E00", "unavailable": "#B8BCC2", "na": "#E8EAED",
    "confirmed": "#009E73", "rescued_ok": "#44AA99", "not_run": "#B8BCC2",
}
REVIEW_EDGE = "#E69F00"


def _bucket(val: str, col: str) -> str:
    v = (val or "").lower()
    if col == "label_reconciliation":
        if "swapped" in v:
            return "review"
        if "consistent" in v or "matches" in v or v in ("", "consistent_with_reference"):
            return "pass"
        return "minor" if v else "na"
    if col == "rescue_status":
        if v.startswith("rescued") or "confirmed" in v:
            return "rescued_ok"
        if "not_suspicious" in v or v in ("", "none"):
            return "pass"
        if "manual" in v or "exclude" in v:
            return "review"
        return "pass"
    if col == "orthology_status":
        if "pass" in v or "supported" in v or "ortholog" in v:
            return "pass"
        if "review" in v or "ambiguous" in v:
            return "review"
        return "minor"
    if col == "coordinate_validation_status":
        if "major" in v or "fail" in v:
            return "fail"
        if "moderate" in v or "review" in v:
            return "review"
        return "pass"
    if col == "MSA_full_length_status":
        if "unavailable" in v:
            return "unavailable"
        if "outlier" in v or "review" in v:
            return "review"
        if "minor" in v:
            return "minor"
        return "pass"
    if col == "boundary_robustness_class":
        if "robust" in v:
            return "robust"
        if "review" in v or "unresolved" in v:
            return "review"
        return "minor"
    if col == "reference_agreement":
        try:
            pct = float(v.replace("%", ""))
            if pct >= 0.9:
                return "pass"
            if pct >= 0.7:
                return "minor"
            return "review"
        except Exception:
            return "na"
    if col == "synteny_validation_class":
        if "strong" in v:
            return "strong"
        if "supported" in v or "partial" in v:
            return "minor"
        if "conflict" in v or "unavailable" in v:
            return "review"
        return "minor"
    if col == "pre_interpro_readiness_class":
        if "primary" in v and "minor" not in v:
            return "primary"
        if "minor" in v:
            return "minor"
        if "supplement" in v:
            return "supplement"
        if "excluded" in v or "not_ready" in v:
            return "fail"
        return "na"
    return "na"


def _color(bucket: str) -> str:
    return STATE.get(bucket, STATE["na"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    cdir = M.closure_dir(base)
    cdir.mkdir(parents=True, exist_ok=True)
    (cdir / "figures").mkdir(exist_ok=True)
    (cdir / "tables").mkdir(exist_ok=True)

    truth_path = cdir / "final_pre_interpro_truth_table.tsv"
    if not truth_path.exists():
        # allow running before closure — build from post-rescue
        import run_fgfr2_final_pre_interpro_closure as C  # noqa: E402
        truth = C.build_truth_table(base)
    else:
        truth = M.read_tsv(truth_path)

    md = M.module_dir(base)
    ref_pct = {}
    for r in M.read_tsv(md / "robustness" / "fgfr2_boundary_robustness_scores.tsv"):
        ref_pct[(r["species"], r.get("isoform", ""))] = r.get(
            "reference_agreement_percent_identical_or_conservative", "")

    table_rows = []
    for r in truth:
        sp, iso = r["species"], r["isoform"]
        claim = r.get("final_claim_status_after_rescue", "")
        is_review = not M.claim_is_primary(str(claim))
        row_data = {
            "species": sp, "isoform": iso,
            "display_species_name": r.get("display_species_name", sp),
            "taxon_group": r.get("taxon_group", ""),
            "final_claim_status_after_rescue": claim,
            "visual_review_flag": "true" if is_review else "false",
        }
        cells = {}
        vals = {
            "label_reconciliation": r.get("label_consistency_status", ""),
            "rescue_status": r.get("rescue_decision", ""),
            "orthology_status": r.get("orthology_status", ""),
            "coordinate_validation_status": r.get("coordinate_validation_status", ""),
            "MSA_full_length_status": r.get("MSA_full_length_status", ""),
            "boundary_robustness_class": r.get("boundary_robustness_class", ""),
            "reference_agreement": ref_pct.get((sp, iso), ""),
            "synteny_validation_class": r.get("synteny_validation_class", ""),
            "pre_interpro_readiness_class": r.get("pre_interpro_readiness_class", ""),
        }
        for col_id, _ in COLS:
            b = _bucket(str(vals.get(col_id, "")), col_id)
            cells[col_id] = b
            row_data[col_id] = b
            row_data[f"{col_id}_raw"] = vals.get(col_id, "")
        table_rows.append(row_data)

    M.write_tsv(cdir / "tables" / "figure_final_framework_evidence_stack.tsv", table_rows,
                list(table_rows[0].keys()) if table_rows else ["species"])

    # order: primary block first, then a labelled supplement/review block (Part E readability)
    primary = [r for r in table_rows if r.get("visual_review_flag") != "true"]
    review = [r for r in table_rows if r.get("visual_review_flag") == "true"]

    S.apply_rcparams()
    import matplotlib.pyplot as plt
    m = len(COLS)
    n_rows = len(table_rows)
    fig_h = max(6.5, n_rows * 0.26 + 2.4)
    fig, ax = plt.subplots(figsize=(max(11.5, m * 1.0 + 4.5), fig_h))

    blocks = [("Primary (accepted, incl. rescued)", primary)]
    if review:
        blocks.append(("Supplement / review only", review))

    y = 0.0
    yticks, ylabels = [], []
    group_spans = []  # (taxon, y_start, y_end) for sidebars
    cur_taxon, span_start = None, None
    for title_txt, block in blocks:
        # block header band
        ax.add_patch(Rectangle((-0.02, y), m + 0.5, 0.7, facecolor="#F0F1F3",
                               edgecolor="none", clip_on=False))
        ax.text(0.05, y + 0.35, title_txt, fontsize=S.FONT["label"], fontweight="bold",
                color=S.INK, va="center", ha="left")
        y += 0.85
        if cur_taxon is not None:
            group_spans.append((cur_taxon, span_start, y - 0.85))
            cur_taxon, span_start = None, None
        for row in block:
            is_review = row.get("visual_review_flag") == "true"
            tax = row.get("taxon_group", "")
            if tax != cur_taxon:
                if cur_taxon is not None:
                    group_spans.append((cur_taxon, span_start, y))
                cur_taxon, span_start = tax, y
            for xi, (col_id, _) in enumerate(COLS):
                b = row[col_id]
                ax.add_patch(Rectangle((xi, y), 1, 1, facecolor=_color(b),
                                       edgecolor="white", lw=0.5))
            # claim sidebar
            ax.add_patch(Rectangle((m + 0.12, y), 0.32, 1,
                                   facecolor=S.CLAIM_COLOR.get(
                                       "supplement" if is_review else "primary", "#C9CDD2"),
                                   edgecolor="none", clip_on=False))
            if is_review:
                ax.add_patch(Rectangle((-0.06, y), m + 0.55, 1, fill=False,
                                       edgecolor=REVIEW_EDGE, lw=1.0, ls=":", clip_on=False))
            yticks.append(y + 0.5)
            ylabels.append(f"{row['display_species_name']} {row['isoform']}")
            y += 1
        if cur_taxon is not None:
            group_spans.append((cur_taxon, span_start, y))
            cur_taxon, span_start = None, None

    # taxon group sidebars + labels on the far left
    for tax, y0, y1 in group_spans:
        if not tax or y1 - y0 < 0.5:
            continue
        ax.plot([-0.95, -0.95], [y0 + 0.06, y1 - 0.06], color=S.MUTED, lw=2.0,
                solid_capstyle="round", clip_on=False)
        ax.text(-1.05, (y0 + y1) / 2, tax, rotation=90, fontsize=S.FONT["small"],
                color=S.MUTED, va="center", ha="right")

    ax.set_xlim(-1.1, m + 0.6)
    ax.set_ylim(0, y)
    ax.invert_yaxis()
    # column labels both top and bottom for readability
    ax.set_xticks([i + 0.5 for i in range(m)])
    ax.set_xticklabels([lab for _, lab in COLS], rotation=40, ha="left", fontsize=S.FONT["tick"])
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(axis="x", length=0)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=S.FONT["small"])
    ax.tick_params(axis="y", length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Integrated annotation-aware evidence stack (pre-InterPro)",
                 fontsize=S.FONT["title"], fontweight="bold", loc="left", pad=26)
    leg = [Patch(facecolor=c, label=k) for k, c in
           [("pass / primary", STATE["pass"]), ("minor / supported", STATE["minor"]),
            ("confirmed / rescued", STATE["confirmed"]),
            ("review / supplement", STATE["review"]), ("fail / excluded", STATE["fail"]),
            ("not applicable", STATE["na"])]]
    ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -0.04),
              ncol=6, fontsize=S.FONT["legend"], frameon=False, handlelength=1.1,
              columnspacing=1.2, handletextpad=0.5)
    stem = "Figure_Final_Framework_Evidence_Stack"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(cdir / "figures" / f"{stem}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] framework evidence stack -> {cdir / 'figures' / stem}.{{svg,pdf,png}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
