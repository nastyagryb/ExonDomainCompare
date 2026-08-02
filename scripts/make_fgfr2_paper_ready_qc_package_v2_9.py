#!/usr/bin/env python3
"""
Make a paper-ready QC package for FGFR2 IIIb/IIIc analysis.

This script does NOT change biological calls. It classifies species/records into:
  - main_analysis_eligible: robust enough for primary claims/figures
  - supplementary_review: shown as QC/supplement, not used for strong conclusions
  - excluded_no_pair: cannot support IIIb/IIIc pair claims

Inputs are the v5.9 marker-validated outputs plus paralog screen outputs.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import textwrap
import pandas as pd
import numpy as np


def read_tsv(path: str | Path, required: bool = True) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        if required:
            raise FileNotFoundError(str(p))
        return pd.DataFrame()
    if p.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(p, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()




def safe_str(v) -> str:
    if pd.isna(v):
        return ""
    return str(v)


def classify_species(row: pd.Series, anchor_map: pd.DataFrame, paralog_map: dict) -> dict:
    species = safe_str(row.get("species"))
    reasons = []
    warnings = []

    pair_status = safe_str(row.get("pair_audit_status"))
    human_status = safe_str(row.get("human_control_status"))
    regional_identity = row.get("regional_local_identity", np.nan)

    if species == "homo_sapiens" and human_status not in {"pass", "not_applicable"}:
        reasons.append("human_positive_control_not_passed")

    if pair_status == "III_region_sequence_distinct":
        pass
    elif pair_status == "III_region_nearly_identical":
        reasons.append("sequence_uninformative_nearly_identical_IIIb_IIIc_windows")
    elif pair_status == "missing_pair_member":
        reasons.append("missing_IIIb_or_IIIc_pair_member")
    else:
        reasons.append(f"pair_status_{pair_status or 'missing'}")

    sub = anchor_map[anchor_map["species"].astype(str).eq(species)].copy() if not anchor_map.empty and "species" in anchor_map.columns else pd.DataFrame()
    if sub.empty:
        reasons.append("no_anchor_map_records")
    else:
        # Conservative record-level warnings.
        for col, badvals in [
            ("full_region_status", {"full_region_low_confidence_or_missing", "missing", "not_available"}),
            ("anchor_position_status", {"position_outlier"}),
            ("candidate_window_mode", {"fallback_fixed_window_low_confidence_anchor"}),
            ("human_anchor_warning", {"anchor_weak", "full_region_low_confidence_or_missing"}),
        ]:
            if col in sub.columns:
                vals = set(sub[col].dropna().astype(str))
                hits = sorted(v for v in vals if v in badvals)
                if hits:
                    warnings.append(f"{col}:{','.join(hits)}")
        if "final_anchor_status_simplified" in sub.columns:
            vals = set(sub["final_anchor_status_simplified"].dropna().astype(str))
            if any("region not reliably detected" in v for v in vals):
                reasons.append("region_not_reliably_detected")
            elif any("sequence uninformative" == v for v in vals):
                reasons.append("sequence_uninformative_record")
            elif any("exon supported; sequence uninformative" == v for v in vals):
                # keep as review reason only if pair is not sequence-distinct
                if pair_status != "III_region_sequence_distinct":
                    reasons.append("exon_supported_sequence_uninformative")

    paralog_status = paralog_map.get(species, "missing_paralog_screen")
    if paralog_status != "all_high_confidence_FGFR2":
        warnings.append(f"paralog_screen:{paralog_status}")
        # Not necessarily fatal if probable, but not primary eligible.
        if "high_confidence" not in paralog_status:
            reasons.append("paralog_screen_not_high_confidence")

    # Strict main figure criteria.
    fatal = [r for r in reasons if r]
    has_hard_warning = any(
        w.startswith("full_region_status:")
        or w.startswith("anchor_position_status:")
        or w.startswith("candidate_window_mode:")
        or w.startswith("human_anchor_warning:")
        for w in warnings
    )

    main = (not fatal) and (not has_hard_warning) and pair_status == "III_region_sequence_distinct" and paralog_status == "all_high_confidence_FGFR2"

    if main:
        qc_class = "main_analysis_eligible"
    elif pair_status in {"III_region_sequence_distinct", "III_region_nearly_identical"}:
        qc_class = "supplementary_review_not_primary_claim"
    else:
        qc_class = "excluded_no_reliable_pair_claim"

    return {
        "species": species,
        "main_analysis_eligible": int(main),
        "qc_class": qc_class,
        "primary_exclusion_reasons": ";".join(fatal) if fatal else "none",
        "qc_warnings": ";".join(warnings) if warnings else "none",
        "pair_audit_status": pair_status,
        "human_control_status": human_status,
        "regional_local_identity": regional_identity,
        "paralog_screen_status": paralog_status,
        "IIIb_transcript": row.get("IIIb_transcript", ""),
        "IIIc_transcript": row.get("IIIc_transcript", ""),
        "IIIb_protein": row.get("IIIb_protein", ""),
        "IIIc_protein": row.get("IIIc_protein", ""),
    }


def write_if_not_empty(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def make_markdown(outdir: Path, species_qc: pd.DataFrame, anchor_map: pd.DataFrame, pair_audit: pd.DataFrame):
    counts = species_qc["qc_class"].value_counts().to_dict()
    pair_counts = pair_audit["pair_audit_status"].value_counts().to_dict() if "pair_audit_status" in pair_audit else {}
    main_species = species_qc.loc[species_qc["main_analysis_eligible"].eq(1), "species"].tolist()
    review = species_qc.loc[species_qc["main_analysis_eligible"].ne(1), ["species", "qc_class", "primary_exclusion_reasons", "qc_warnings"]]

    human = species_qc[species_qc["species"].eq("homo_sapiens")]
    human_status = human["human_control_status"].iloc[0] if not human.empty else "missing"

    md = []
    md.append("# FGFR2 paper-ready QC package\n")
    md.append("## Summary\n")
    md.append(f"- Species in pair audit: `{len(species_qc)}`")
    md.append(f"- Main-analysis eligible species: `{len(main_species)}`")
    md.append(f"- Human positive-control status: `{human_status}`")
    md.append("- Pair-audit status counts:")
    for k, v in pair_counts.items():
        md.append(f"  - `{k}`: `{v}`")
    md.append("- QC class counts:")
    for k, v in counts.items():
        md.append(f"  - `{k}`: `{v}`")

    md.append("\n## Main-analysis eligible species\n")
    if main_species:
        for s in main_species:
            md.append(f"- `{s}`")
    else:
        md.append("None under the strict filter.")

    md.append("\n## Review / supplementary / excluded cases\n")
    if review.empty:
        md.append("None.")
    else:
        for _, r in review.iterrows():
            md.append(f"- `{r['species']}` — `{r['qc_class']}`; reasons: `{r['primary_exclusion_reasons']}`; warnings: `{r['qc_warnings']}`")

    md.append("\n## Recommended wording\n")
    md.append(textwrap.dedent("""
    Use cautious, QC-aware wording:

    - The human FGFR2-IIIb/IIIc positive-control pair passed protein-level validation using human-calibrated isoform-specific peptide references.
    - FGFR2-IIIb/IIIc windows were interpreted as sequence-supported only for species passing exon support, dynamic III-region anchoring, sequence distinction, and FGFR2 paralog-screen checks.
    - Species marked as sequence-uninformative, anchor-weak, fallback-window, or position-outlier were retained as review/supplementary cases and not used for strong primary biological claims.
    - The FGFR1-4 paralog identity screen supports FGFR2-like identity, but should be complemented by external orthology evidence or gene-tree/orthogroup inference for paper-level orthology wording.
    """).strip())

    (outdir / "fgfr2_paper_ready_qc_report.md").write_text("\n".join(md) + "\n")


def make_claims_file(outdir: Path, species_qc: pd.DataFrame):
    n_total = len(species_qc)
    n_main = int(species_qc["main_analysis_eligible"].sum())
    n_distinct = int(species_qc["pair_audit_status"].eq("III_region_sequence_distinct").sum())
    n_near = int(species_qc["pair_audit_status"].eq("III_region_nearly_identical").sum())

    txt = f"""# Figure bullet descriptions — QC-vetted draft

## Safe main-text bullets

- The human FGFR2-IIIb/IIIc positive-control pair passed protein-level validation using human-calibrated isoform-specific peptide references.
- Across the {n_total}-species panel, {n_distinct} species had sequence-distinct FGFR2-IIIb/IIIc candidate windows; {n_near} species were retained as sequence-uninformative review cases.
- Under the strict primary-analysis filter, {n_main} species passed the combined QC criteria: FGFR2 paralog-screen support, exon support, dynamic III-region anchoring, sequence distinction, and no major anchor/position warning.
- Records failing the strict filter were not interpreted as negative biological evidence; they were retained as annotation/QC review cases.

## Do NOT claim

- Do not claim that FGFR2 orthology is proven solely by the FGFR1-4 paralog screen.
- Do not claim that all species have validated sequence-distinct IIIb/IIIc isoforms.
- Do not describe the current IIIb/IIIc references as independent external curated references if they were calibrated from the human positive-control pair.

## Orthology wording after adding external evidence

- If Ensembl Compara / OMA / OrthoDB / OrthoFinder supports the same species as FGFR2 orthologs, you may write: “FGFR2-like identity was supported by an FGFR1-4 paralog screen and was consistent with external orthology evidence.”
- If OrthoFinder/gene-tree support is added, stronger wording is possible: “Candidate FGFR2 sequences clustered within the FGFR2 orthogroup/gene-tree clade, supporting orthology in addition to paralog discrimination.”
"""
    (outdir / "figure_bullet_descriptions_QC_vetted.md").write_text(txt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair_audit", required=True)
    ap.add_argument("--anchor_map", required=True)
    ap.add_argument("--paralog_species_summary", required=False)
    ap.add_argument("--pair_difference_positions", required=False)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", default="fgfr2")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pair = read_tsv(args.pair_audit)
    anchor = read_tsv(args.anchor_map)
    para = read_tsv(args.paralog_species_summary, required=False) if args.paralog_species_summary else pd.DataFrame()

    paralog_map = {}
    if not para.empty and {"species", "species_fgfr2_screen_status"}.issubset(para.columns):
        paralog_map = dict(zip(para["species"].astype(str), para["species_fgfr2_screen_status"].astype(str)))

    rows = [classify_species(r, anchor, paralog_map) for _, r in pair.iterrows()]
    qc = pd.DataFrame(rows).sort_values(["main_analysis_eligible", "species"], ascending=[False, True])
    write_if_not_empty(qc, outdir / f"{args.prefix}_paper_ready_species_qc.tsv")

    main_species = set(qc.loc[qc["main_analysis_eligible"].eq(1), "species"].astype(str))
    review_species = set(qc.loc[qc["main_analysis_eligible"].ne(1), "species"].astype(str))

    write_if_not_empty(pair[pair["species"].astype(str).isin(main_species)], outdir / f"{args.prefix}_pair_audit_MAIN.tsv")
    write_if_not_empty(pair[pair["species"].astype(str).isin(review_species)], outdir / f"{args.prefix}_pair_audit_REVIEW.tsv")

    if not anchor.empty and "species" in anchor.columns:
        write_if_not_empty(anchor[anchor["species"].astype(str).isin(main_species)], outdir / f"{args.prefix}_anchor_map_MAIN.tsv")
        write_if_not_empty(anchor[anchor["species"].astype(str).isin(review_species)], outdir / f"{args.prefix}_anchor_map_REVIEW.tsv")

    if args.pair_difference_positions:
        diff = read_tsv(args.pair_difference_positions, required=False)
        if not diff.empty and "species" in diff.columns:
            write_if_not_empty(diff[diff["species"].astype(str).isin(main_species)], outdir / f"{args.prefix}_pair_difference_positions_MAIN.tsv")
            write_if_not_empty(diff[diff["species"].astype(str).isin(review_species)], outdir / f"{args.prefix}_pair_difference_positions_REVIEW.tsv")
        else:
            # Keep an empty but parseable file.
            pd.DataFrame({"note": ["No parseable pair-difference positions available."]}).to_csv(outdir / f"{args.prefix}_pair_difference_positions_MAIN.tsv", sep="\t", index=False)

    make_markdown(outdir, qc, anchor, pair)
    make_claims_file(outdir, qc)

    print("Wrote paper-ready QC package to", outdir)
    print(qc["qc_class"].value_counts().to_string())


if __name__ == "__main__":
    main()
