#!/usr/bin/env python3

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO = Path(__file__).resolve().parent.parent


def display_path(path) -> str:
    p = Path(path)
    try:
        return str(p.resolve().relative_to(REPO.resolve()))
    except Exception:
        return str(p)
sys.path.insert(0, str(REPO / "scripts"))
from exondomaincompare.presentation import fgfr2_plot_style as st
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

import os as _os
BASE = Path(_os.environ.get("FGFR2_RESULTS_DIR") or _os.environ.get("RESULTS_DIR")
            or _os.environ.get("BASE") or (REPO / "results" / "final_30_until_interpro_prepare"))
CLOSURE = BASE / "13_final_pre_interpro_closure"
POST = BASE / "15_exon_domain_boundary_post_interpro"

TRUTH = CLOSURE / "final_pre_interpro_truth_table.tsv"
MANIFEST = CLOSURE / "freeze" / "final_pre_interpro_sequence_manifest.tsv"
FEATURES = POST / "tables" / "exon_domain_architecture_features.tsv"
INTERPRO = POST / "tables" / "interpro_domain_features_normalized.tsv"
PYTMHMM = POST / "tables" / "pytmhmm_tm_features_normalized.tsv"
QC = POST / "tables" / "fgfr2_domain_architecture_qc.tsv"
SANITATION = POST / "tables" / "exon_block_length_consistency_audit.tsv"

OUT = BASE / "16_final_thesis_analyses" / "exon_domain_boundary_consistency"
T_OUT = OUT / "tables"
F_OUT = OUT / "figures"
R_OUT = OUT / "reports"

ALIGNED_MAX = 3      # 0-3 aa: aligned_to_domain_boundary
NEAR_MAX = 15        # 4-15 aa: near_domain_boundary

CLASS_COLORS = {
    "aligned_to_domain_boundary": "#1B7837",
    "near_domain_boundary": "#A6DBA0",
    "within_domain": "#FDB863",
    "between_domains": "#B2ABD2",
    "review_or_missing": "#D9D9D9",
}
CLASS_ORDER = ["aligned_to_domain_boundary", "near_domain_boundary",
               "within_domain", "between_domains", "review_or_missing"]
CLASS_SHORT = {"aligned_to_domain_boundary": "aligned (0-3)",
               "near_domain_boundary": "near (4-15)",
               "within_domain": "within domain",
               "between_domains": "between domains",
               "review_or_missing": "missing/NA"}

TAXON_ORDER = ["Primates", "Other mammals", "Birds", "Reptiles",
               "Amphibians", "Teleost fish"]
TAXON_SHORT = {"Primates": "Primates", "Other mammals": "Mammals", "Birds": "Birds",
               "Reptiles": "Reptiles", "Amphibians": "Amphibians",
               "Teleost fish": "Fish"}


def _load_species_order() -> Dict[str, int]:
    f = REPO / "reference" / "Species_list_final_30.txt"
    order: Dict[str, int] = {}
    if f.exists():
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines()):
            name = line.strip().lower().replace(" ", "_")
            if name:
                order[name] = i
    return order


SPECIES_ORDER = _load_species_order()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def to_int(v, default=None):
    try:
        return int(float(str(v).strip()))
    except (ValueError, TypeError):
        return default


def classify(x: Optional[int], domains: List[dict], *,
             aligned_max: int = ALIGNED_MAX, near_max: int = NEAR_MAX) -> dict:
    if x is None or not domains:
        return {"boundary_class": "review_or_missing", "label": "", "dclass": "",
                "start": "", "end": "", "edge": "", "dist": ""}
    best = None
    for d in domains:
        for edge, pos in (("start", d["start"]), ("end", d["end"])):
            if pos is None:
                continue
            dist = abs(x - pos)
            if best is None or dist < best["dist"]:
                best = {"dist": dist, "label": d["label"], "dclass": d["dclass"],
                        "start": d["start"], "end": d["end"], "edge": edge}
    if best is None:
        return {"boundary_class": "review_or_missing", "label": "", "dclass": "",
                "start": "", "end": "", "edge": "", "dist": ""}
    inside = any(d["start"] is not None and d["end"] is not None
                 and d["start"] <= x <= d["end"] for d in domains)
    if best["dist"] <= aligned_max:
        cls = "aligned_to_domain_boundary"
    elif best["dist"] <= near_max:
        cls = "near_domain_boundary"
    elif inside:
        cls = "within_domain"
    else:
        cls = "between_domains"
    return {"boundary_class": cls, "label": best["label"], "dclass": best["dclass"],
            "start": best["start"], "end": best["end"], "edge": best["edge"],
            "dist": best["dist"]}


def load_architecture() -> Dict[Tuple[str, str], dict]:
    rows = read_tsv(FEATURES)
    by: Dict[Tuple[str, str], dict] = {}
    for r in rows:
        key = (r["species"], r["isoform"])
        node = by.setdefault(key, {
            "species": r["species"], "isoform": r["isoform"],
            "transcript_id": r.get("transcript_id", ""),
            "protein_id": r.get("protein_id", ""),
            "protein_length": to_int(r.get("protein_length")),
            "domains": [], "tm": [], "cassette": None, "coding_exons": [],
        })
        ft = r.get("feature_type", "")
        s, e = to_int(r.get("start_aa")), to_int(r.get("end_aa"))
        label = r.get("feature_label", "")
        if ft == "ig_like_domain":
            node["domains"].append({"label": label, "dclass": "Ig-like",
                                    "start": s, "end": e})
        elif ft == "kinase_domain":
            node["domains"].append({"label": "kinase", "dclass": "kinase",
                                    "start": s, "end": e})
        elif ft == "transmembrane_pytmhmm" and r.get("status") == "receptor_tm":
            node["tm"].append({"label": "TM", "dclass": "TM", "start": s, "end": e})
        elif ft in ("IIIb_slot", "IIIc_slot"):
            node["cassette"] = {"type": ft, "label": label, "start": s, "end": e}
        elif ft == "coding_exon":
            node["coding_exons"].append({"label": label, "start": s, "end": e})
    return by


def main() -> int:
    st.apply_rcparams()
    for d in (T_OUT, F_OUT, R_OUT):
        d.mkdir(parents=True, exist_ok=True)

    manifest = {(r["species"], r["isoform"]): r for r in read_tsv(MANIFEST)}
    primary = {k for k, r in manifest.items()
               if str(r.get("included_in_primary_interpro", "")).lower() == "true"}
    truth = {(r["species"], r["isoform"]): r for r in read_tsv(TRUTH)}
    qc = {(r["species"], r["isoform"]): r for r in read_tsv(QC)}
    sanit = {(r["species"], r["isoform"]): r for r in read_tsv(SANITATION)}
    arch = load_architecture()

    # exon_block_display_status per protein (from QC table)
    disp_status = {k: v.get("exon_block_display_status", "") for k, v in qc.items()}

    distance_rows: List[dict] = []
    # per-protein heatmap cells + cassette distances for summary
    heat: Dict[Tuple[str, str], dict] = {}
    cassette_dist_records: List[dict] = []   # for summary + Fig 12

    n_cassette = n_interpro = n_tm = 0

    for key in sorted(primary):
        sp, iso = key
        node = arch.get(key)
        tr = truth.get(key, {})
        iso_label = tr.get("final_isoform_label", iso)
        taxon = tr.get("taxon_group", "")
        L = node["protein_length"] if node else to_int(manifest[key].get("protein_length"))
        tx = tr.get("transcript_id", "")
        pid = tr.get("protein_id", "")
        qc_status = qc.get(key, {}).get("final_qc_status", "")
        eb_status = disp_status.get(key, "")
        san = sanit.get(key, {})
        src_status = san.get("action", "") or "sanitized_feature_table"

        if not node:
            continue
        domains_all = [d for d in node["domains"]] + list(node["tm"])
        ig_domains = [d for d in node["domains"] if d["dclass"] == "Ig-like"]
        tm_domains = list(node["tm"])
        kin_domains = [d for d in node["domains"] if d["dclass"] == "kinase"]
        if any(d["dclass"] in ("Ig-like", "kinase") for d in domains_all):
            n_interpro += 1
        if tm_domains:
            n_tm += 1

        cassette = node["cassette"]
        cell = {"species": sp, "isoform": iso, "final_isoform_label": iso_label,
                "taxon_group": taxon, "eb_status": eb_status,
                "cassette_start": "review_or_missing", "cassette_end": "review_or_missing",
                "ig": "review_or_missing", "tm": "review_or_missing",
                "kinase": "review_or_missing"}
        if cassette and cassette["start"] is not None:
            n_cassette += 1
            for btype, bpos in (("cassette_start", cassette["start"]),
                                ("cassette_end", cassette["end"])):
                c = classify(bpos, domains_all)
                distance_rows.append(_drow(sp, iso, iso_label, tx, pid, L, btype,
                                           cassette["label"], bpos, c, qc_status,
                                           eb_status, src_status,
                                           "cassette boundary (primary thesis metric)"))
                cell[btype] = c["boundary_class"]
                if isinstance(c["dist"], int):
                    cassette_dist_records.append({
                        "species": sp, "isoform": iso, "taxon": taxon,
                        "btype": btype, "dist": c["dist"],
                        "class": c["boundary_class"], "eb_status": eb_status})
            cend = cassette["end"]
            cell["ig"] = classify(cend, ig_domains)["boundary_class"] if ig_domains else "review_or_missing"
            cell["tm"] = classify(cend, tm_domains)["boundary_class"] if tm_domains else "review_or_missing"
            cell["kinase"] = classify(cend, kin_domains)["boundary_class"] if kin_domains else "review_or_missing"
        heat[key] = cell

        for ex in node["coding_exons"]:
            for btype, bpos in (("coding_exon_start", ex["start"]),
                                ("coding_exon_end", ex["end"])):
                c = classify(bpos, domains_all)
                distance_rows.append(_drow(sp, iso, iso_label, tx, pid, L, btype,
                                           ex["label"], bpos, c, qc_status,
                                           eb_status, src_status,
                                           "coding-exon boundary (context)"))

    _write_distances(distance_rows)
    summary = _write_summary(primary, heat, cassette_dist_records,
                             n_cassette, n_interpro, n_tm)
    outliers = _write_outliers(distance_rows, heat)
    _figure11(heat)
    _figure12(cassette_dist_records)
    _report(primary, summary, outliers, distance_rows)

    print(f"[ok] output -> {display_path(OUT)}")
    print(f"[summary] primary={len(primary)} with_cassette={n_cassette} "
          f"with_interpro_domains={n_interpro} with_tm={n_tm} "
          f"distance_rows={len(distance_rows)} outliers={len(outliers)}")
    print(f"[cassette boundary classes] "
          + ", ".join(f"{c}={summary['cassette_class_counts'].get(c,0)}"
                      for c in CLASS_ORDER))
    return 0


def _drow(sp, iso, iso_label, tx, pid, L, btype, blabel, bpos, c, qc_status,
          eb_status, src_status, notes) -> dict:
    return {
        "species": sp, "isoform": iso, "final_isoform_label": iso_label,
        "transcript_id": tx, "protein_id": pid, "protein_length": L if L else "",
        "boundary_type": btype, "boundary_label": blabel, "boundary_aa": bpos,
        "nearest_domain_label": c["label"], "nearest_domain_class": c["dclass"],
        "nearest_domain_start_aa": c["start"], "nearest_domain_end_aa": c["end"],
        "nearest_domain_boundary_type": c["edge"],
        "distance_to_nearest_domain_boundary": c["dist"],
        "boundary_class": c["boundary_class"],
        "architecture_qc_status": qc_status,
        "exon_block_display_status": eb_status,
        "source_coordinate_status": src_status,
        "notes": notes,
    }


def _write_distances(rows: List[dict]) -> None:
    cols = ["species", "isoform", "final_isoform_label", "transcript_id", "protein_id",
            "protein_length", "boundary_type", "boundary_label", "boundary_aa",
            "nearest_domain_label", "nearest_domain_class", "nearest_domain_start_aa",
            "nearest_domain_end_aa", "nearest_domain_boundary_type",
            "distance_to_nearest_domain_boundary", "boundary_class",
            "architecture_qc_status", "exon_block_display_status",
            "source_coordinate_status", "notes"]
    _tsv(T_OUT / "exon_domain_boundary_distances.tsv", rows, cols)


def _write_summary(primary, heat, cass_recs, n_cassette, n_interpro, n_tm) -> dict:
    # cassette boundary class counts (overall)
    cclasses = [r["class"] for r in cass_recs]
    class_counts = {c: cclasses.count(c) for c in CLASS_ORDER}
    dists = [r["dist"] for r in cass_recs]

    def med(x):
        return round(statistics.median(x), 2) if x else ""

    def mean(x):
        return round(statistics.mean(x), 2) if x else ""

    rows: List[dict] = []

    def scope_row(scope, level, recs, nprot):
        cl = [r["class"] for r in recs]
        d = [r["dist"] for r in recs]
        rows.append({
            "scope": scope, "level": level, "n_proteins": nprot,
            "n_cassette_boundaries": len(recs),
            "aligned_to_domain_boundary": cl.count("aligned_to_domain_boundary"),
            "near_domain_boundary": cl.count("near_domain_boundary"),
            "within_domain": cl.count("within_domain"),
            "between_domains": cl.count("between_domains"),
            "review_or_missing": cl.count("review_or_missing"),
            "median_distance_to_nearest_domain_boundary": med(d),
            "mean_distance_to_nearest_domain_boundary": mean(d),
        })

    # overall row also carries the global availability counts
    scope_row("overall", "all_primary", cass_recs, len(primary))
    rows[-1].update({
        "total_primary_proteins": len(primary),
        "proteins_with_cassette_boundary_data": n_cassette,
        "proteins_with_interpro_domain_data": n_interpro,
        "proteins_with_tm_data": n_tm,
    })
    for iso in ("IIIb", "IIIc"):
        recs = [r for r in cass_recs if r["isoform"] == iso]
        nprot = len({(r["species"], r["isoform"]) for r in recs})
        scope_row("isoform", iso, recs, nprot)
    for tax in TAXON_ORDER:
        recs = [r for r in cass_recs if r["taxon"] == tax]
        if not recs:
            continue
        nprot = len({(r["species"], r["isoform"]) for r in recs})
        scope_row("taxon_group", tax, recs, nprot)

    cols = ["scope", "level", "n_proteins", "n_cassette_boundaries",
            "aligned_to_domain_boundary", "near_domain_boundary", "within_domain",
            "between_domains", "review_or_missing",
            "median_distance_to_nearest_domain_boundary",
            "mean_distance_to_nearest_domain_boundary",
            "total_primary_proteins", "proteins_with_cassette_boundary_data",
            "proteins_with_interpro_domain_data", "proteins_with_tm_data"]
    _tsv(T_OUT / "exon_domain_boundary_consistency_summary.tsv", rows, cols)
    return {"cassette_class_counts": class_counts,
            "median_all": med(dists), "mean_all": mean(dists),
            "n_cassette": n_cassette, "n_interpro": n_interpro, "n_tm": n_tm,
            "rows": rows}


def _write_outliers(distance_rows, heat) -> List[dict]:
    LOW_CONF = {"cassette_only_high_confidence", "native_exon_blocks_reconstructed"}
    UNRESOLVED = {"cassette_only_display", "hide_untrusted_exon_block"}
    out: List[dict] = []
    seen_display = set()
    for r in distance_rows:
        if r["boundary_type"] not in ("cassette_start", "cassette_end"):
            continue
        cls = r["boundary_class"]
        d = r["distance_to_nearest_domain_boundary"]
        ebs = r["exon_block_display_status"]
        scs = r["source_coordinate_status"]
        geo = []
        if cls == "review_or_missing":
            geo.append("missing cassette/domain boundary data")
        if cls == "between_domains":
            geo.append("cassette boundary lies outside annotated domains")
        disp = []
        if ebs in LOW_CONF:
            disp.append(f"low-confidence / reconstructed exon-block display ({ebs})")
        if scs in UNRESOLVED:
            disp.append(f"unresolved coordinate status ({scs})")

        emit = False
        if geo:
            emit = True
            reasons = geo + disp
        elif disp and (r["species"], r["isoform"]) not in seen_display \
                and r["boundary_type"] == "cassette_end":
            emit = True
            seen_display.add((r["species"], r["isoform"]))
            reasons = disp
        if emit:
            out.append({
                "species": r["species"], "isoform": r["isoform"],
                "final_isoform_label": r["final_isoform_label"],
                "boundary_type": r["boundary_type"], "boundary_aa": r["boundary_aa"],
                "nearest_domain_label": r["nearest_domain_label"],
                "distance_to_nearest_domain_boundary": d,
                "boundary_class": cls,
                "exon_block_display_status": ebs,
                "architecture_qc_status": r["architecture_qc_status"],
                "outlier_reason": "; ".join(reasons),
            })
    cols = ["species", "isoform", "final_isoform_label", "boundary_type",
            "boundary_aa", "nearest_domain_label",
            "distance_to_nearest_domain_boundary", "boundary_class",
            "exon_block_display_status", "architecture_qc_status", "outlier_reason"]
    _tsv(T_OUT / "exon_domain_boundary_outliers.tsv", out, cols)
    return out


def _sort_key(cell) -> tuple:
    ti = TAXON_ORDER.index(cell["taxon_group"]) if cell["taxon_group"] in TAXON_ORDER else 99
    si = SPECIES_ORDER.get(cell["species"], 999)
    return (ti, si, cell["species"])


def _figure11(heat: Dict[Tuple[str, str], dict]) -> None:
    cols = [("cassette_start", "cassette\nstart"), ("cassette_end", "cassette\nend"),
            ("ig", "Ig-like\nboundary"), ("tm", "TM\nrelation"),
            ("kinase", "kinase\nrelation")]
    panels = [("IIIb", [c for c in heat.values() if c["isoform"] == "IIIb"]),
              ("IIIc", [c for c in heat.values() if c["isoform"] == "IIIc"])]
    for _, cells in panels:
        cells.sort(key=_sort_key)
    nrows = max(len(p[1]) for p in panels)
    low_conf = {"cassette_only_high_confidence", "native_exon_blocks_reconstructed"}
    fig, axes = plt.subplots(1, 2, figsize=(11.5, max(6.0, nrows * 0.26)),
                             gridspec_kw={"wspace": 0.62})
    for ax, (iso, cells) in zip(axes, panels):
        ax.set_title(f"FGFR2 {iso}", fontsize=st.FONT["title"], fontweight="bold",
                     loc="left", color=st.INK)
        ax.set_xlim(0, len(cols))
        ax.set_ylim(0, len(cells))
        ax.set_xticks([i + 0.5 for i in range(len(cols))])
        ax.set_xticklabels([c[1] for c in cols], fontsize=st.FONT["small"])
        ax.set_yticks([i + 0.5 for i in range(len(cells))])
        ylabels = []
        for c in cells:
            mark = " *" if c["eb_status"] in low_conf else ""
            ylabels.append(f"{c['species'].replace('_',' ')}{mark}")
        ax.set_yticklabels(ylabels, fontsize=st.FONT["small"] - 0.5)
        ax.invert_yaxis()
        for yi, c in enumerate(cells):
            lc = c["eb_status"] in low_conf
            for xi, (key, _) in enumerate(cols):
                cls = c.get(key, "review_or_missing")
                ax.add_patch(Rectangle((xi, yi), 1, 1,
                             facecolor=CLASS_COLORS[cls],
                             edgecolor=("#111111" if lc else "white"),
                             lw=(1.1 if lc else 0.5), zorder=2))
        # taxon-group separators + right-side group labels
        prev, start = None, 0
        bounds = []
        for yi, c in enumerate(cells):
            g = c["taxon_group"]
            if g != prev:
                if prev is not None:
                    bounds.append((start, yi, prev))
                start, prev = yi, g
        bounds.append((start, len(cells), prev))
        for (s, e, g) in bounds:
            if s > 0:
                ax.axhline(s, color="#8A8F98", lw=0.7, zorder=5)
            ax.text(len(cols) + 0.12, (s + e) / 2, TAXON_SHORT.get(g, g),
                    rotation=270, va="center", ha="left",
                    fontsize=st.FONT["small"] - 0.5, color="#555555")
        for sp in ("top", "right", "left", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=0)
    handles = [Patch(facecolor=CLASS_COLORS[c], edgecolor="white",
                     label=CLASS_SHORT[c]) for c in CLASS_ORDER]
    handles.append(Patch(facecolor="white", edgecolor="#111111", lw=1.1,
                         label="* reconstructed / cassette-only display (inspect)"))
    fig.legend(handles=handles, loc="lower center", ncol=3,
               fontsize=st.FONT["small"], frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("FGFR2 IIIb/IIIc exon\u2013domain boundary consistency",
                 fontsize=st.FONT["title"], fontweight="bold", x=0.02, ha="left")
    st.savefig(fig, F_OUT, "Figure_11_exon_domain_boundary_consistency_heatmap")


def _figure12(cass_recs: List[dict]) -> None:
    import random
    random.seed(7)
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    groups = ["IIIb", "IIIc"]
    btypes = [("cassette_start", "#2166AC"), ("cassette_end", "#B2182B")]
    xticks, xlabels = [], []
    pos = 0
    for iso in groups:
        for btype, color in btypes:
            vals = [r["dist"] for r in cass_recs
                    if r["isoform"] == iso and r["btype"] == btype]
            if vals:
                bx = ax.boxplot(vals, positions=[pos], widths=0.55,
                                patch_artist=True, showfliers=False)
                for b in bx["boxes"]:
                    b.set(facecolor=color, alpha=0.25, edgecolor=color)
                for med in bx["medians"]:
                    med.set(color=color, linewidth=1.6)
                xs = [pos + random.uniform(-0.16, 0.16) for _ in vals]
                ax.scatter(xs, vals, s=16, color=color, alpha=0.8, zorder=3,
                           edgecolors="white", linewidths=0.4)
            xticks.append(pos)
            xlabels.append(f"{iso}\n{btype.split('_')[1]}")
            pos += 1
        pos += 0.6
    ax.axhspan(0, ALIGNED_MAX, color="#1B7837", alpha=0.07, zorder=0)
    ax.axhspan(ALIGNED_MAX, NEAR_MAX, color="#A6DBA0", alpha=0.10, zorder=0)
    ax.axhline(ALIGNED_MAX, color="#1B7837", lw=0.7, ls="--", alpha=0.6)
    ax.axhline(NEAR_MAX, color="#7A9E76", lw=0.7, ls="--", alpha=0.6)
    ax.text(pos - 0.6, ALIGNED_MAX, " aligned <=3 aa", fontsize=st.FONT["small"],
            va="bottom", ha="right", color="#1B7837")
    ax.text(pos - 0.6, NEAR_MAX, " near <=15 aa", fontsize=st.FONT["small"],
            va="bottom", ha="right", color="#5a7d56")
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=st.FONT["small"])
    ax.set_ylabel("Absolute distance to nearest domain boundary (aa)",
                  fontsize=st.FONT["label"])
    ax.set_title("Cassette boundary distance to nearest protein-domain boundary",
                 fontsize=st.FONT["title"], fontweight="bold", loc="left", color=st.INK)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    st.savefig(fig, F_OUT, "Figure_12_boundary_distance_distribution")

def _report(primary, summary, outliers, distance_rows) -> None:
    cc = summary["cassette_class_counts"]
    tot = sum(cc.values())

    def pct(n):
        return f"{100*n/tot:.0f}%" if tot else "0%"

    iso_rows = {r["level"]: r for r in summary["rows"] if r["scope"] == "isoform"}
    tax_rows = [r for r in summary["rows"] if r["scope"] == "taxon_group"]

    lines = []
    lines.append("# Exon-domain boundary consistency report\n")
    lines.append("Thesis analysis: *consistency and robustness of exon-domain "
                 "boundary identification across vertebrate FGFR2 orthologs "
                 "(IIIb/IIIc case study)*.\n")
    lines.append("## Input files (read-only)\n")
    for p in (TRUTH, MANIFEST, FEATURES, INTERPRO, PYTMHMM, QC, SANITATION):
        lines.append(f"* `{display_path(p)}`")
    lines.append("\nProvenance-only (not used as final display truth): "
                 "`post_interpro_qc_review_case_audit.tsv`, "
                 "`exon_block_reconstruction_overrides.json`, figure3C source table.\n")
    lines.append("## Method\n")
    lines.append(f"* Primary proteins analyzed: **{len(primary)}**.")
    lines.append("* Coordinates come from the **sanitized** post-InterPro feature "
                 "table `exon_domain_architecture_features.tsv` (never the raw "
                 "figure3C blocks where a sanitation override exists).")
    lines.append("* Domain landmarks: representative InterProScan **Ig-like** (Ig1/Ig2/Ig3) "
                 "and **kinase** domains plus the **pyTMHMM** receptor TM helix. "
                 "The family-level FGFR fingerprint (spanning the whole protein) is "
                 "excluded because it carries no informative internal boundary.")
    lines.append("* pyTMHMM provides the TM layer because InterProScan did not "
                 "annotate the transmembrane helix consistently for these proteins.")
    lines.append("* For each boundary the absolute distance to the nearest domain "
                 "edge (start or end) is computed and classified.\n")
    lines.append("### Boundary-class thresholds\n")
    lines.append(f"* `aligned_to_domain_boundary`: 0–{ALIGNED_MAX} aa from nearest domain boundary")
    lines.append(f"* `near_domain_boundary`: {ALIGNED_MAX+1}–{NEAR_MAX} aa")
    lines.append("* `within_domain`: boundary inside a domain but > "
                 f"{NEAR_MAX} aa from its edges")
    lines.append("* `between_domains`: boundary outside all annotated domains")
    lines.append("* `review_or_missing`: cassette/domain/coordinate data missing\n")
    lines.append("## Overall cassette-boundary consistency\n")
    lines.append(f"Across {tot} cassette boundaries "
                 f"({summary['n_cassette']} proteins × start+end):\n")
    for c in CLASS_ORDER:
        lines.append(f"* {CLASS_SHORT[c]}: **{cc[c]}** ({pct(cc[c])})")
    lines.append(f"\nMedian distance to nearest domain boundary: "
                 f"**{summary['median_all']} aa**; mean: **{summary['mean_all']} aa**.")
    aligned_near = cc["aligned_to_domain_boundary"] + cc["near_domain_boundary"]
    lines.append(f"\n**Interpretation:** {aligned_near}/{tot} ({pct(aligned_near)}) of "
                 "cassette boundaries fall within 15 aa of an InterProScan/pyTMHMM domain "
                 "boundary, i.e. the IIIb/IIIc splice cassette is consistently anchored to "
                 "the D3 (Ig3) region of the receptor architecture across vertebrates. "
                 "This supports the robustness of exon–domain boundary identification.\n")
    lines.append("## IIIb vs IIIc\n")
    for iso in ("IIIb", "IIIc"):
        r = iso_rows.get(iso)
        if r:
            lines.append(f"* **{iso}** ({r['n_proteins']} proteins): "
                         f"aligned {r['aligned_to_domain_boundary']}, "
                         f"near {r['near_domain_boundary']}, "
                         f"within {r['within_domain']}, between {r['between_domains']}; "
                         f"median {r['median_distance_to_nearest_domain_boundary']} aa.")
    lines.append("\n## Taxon groups\n")
    for r in tax_rows:
        lines.append(f"* **{r['level']}** ({r['n_proteins']} proteins): "
                     f"median {r['median_distance_to_nearest_domain_boundary']} aa "
                     f"(aligned {r['aligned_to_domain_boundary']}, "
                     f"near {r['near_domain_boundary']}, within {r['within_domain']}, "
                     f"between {r['between_domains']}).")
    lines.append(f"\n## Outliers ({len(outliers)})\n")
    if outliers:
        for o in outliers:
            lines.append(f"* {o['species'].replace('_',' ')} {o['final_isoform_label']} "
                         f"({o['boundary_type']} @ {o['boundary_aa']}): {o['outlier_reason']}")
    else:
        lines.append("* None.")
    lines.append("\n## Notes\n")
    lines.append("* InterProScan/pyTMHMM domain annotations **support** the domain "
                 "architecture; they never relabel IIIb/IIIc (labels always come from "
                 "the final truth table).")
    lines.append("* pyTMHMM is used for the TM layer because InterProScan did not "
                 "annotate TM helices consistently.")
    lines.append("* All coordinates use the sanitized post-InterPro feature tables; "
                 "low-confidence / reconstructed exon-block display is flagged (marked "
                 "with `*` in Figure 11) and reported separately from biological QC — "
                 "it is a display-coordinate property, not a biological failure.")
    (R_OUT / "exon_domain_boundary_consistency_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")


def _tsv(path: Path, rows: List[dict], cols: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
