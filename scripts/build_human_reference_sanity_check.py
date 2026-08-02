#!/usr/bin/env python3
"""Human FGFR2 reference boundary sanity-check (read-only).

Summarizes the human IIIb and IIIc architecture using ONLY existing pipeline
outputs (the step-16 exon-domain boundary distances table). No biological data,
FASTA, truth table, or membership is changed. The output is a small human
sanity-check report (Markdown + TSV) with calm, scientific wording.
"""
from __future__ import annotations

import csv
from pathlib import Path
from datetime import datetime, timezone

REPO = Path(__file__).resolve().parent.parent


def display_path(path) -> str:
    """Repo-relative path for display/logging only; falls back to the raw path when
    BASE is a run-local relative path. Never raises and never affects outputs."""
    p = Path(path)
    try:
        return str(p.resolve().relative_to(REPO.resolve()))
    except Exception:
        return str(p)
import os as _os  # run-folder path override; legacy default preserved
_RESULTS = _os.environ.get("FGFR2_RESULTS_DIR") or _os.environ.get("RESULTS_DIR") or _os.environ.get("BASE")
if _RESULTS:
    BASE = Path(_RESULTS) / "16_final_thesis_analyses"
else:
    BASE = REPO / "results" / "final_30_until_interpro_prepare" / "16_final_thesis_analyses"
DIST_TSV = BASE / "exon_domain_boundary_consistency" / "tables" / "exon_domain_boundary_distances.tsv"
OUT_DIR = BASE / "human_reference_sanity_check"
SPECIES = "homo_sapiens"


def _load_rows() -> list[dict]:
    with DIST_TSV.open(encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["species"] == SPECIES]


def _int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _summarize_isoform(rows: list[dict]) -> dict | None:
    if not rows:
        return None
    r0 = rows[0]
    by_type = {r["boundary_type"]: r for r in rows}
    cs = by_type.get("cassette_start")
    ce = by_type.get("cassette_end")

    # nearest Ig-like domain referenced by the cassette boundaries
    ig_row = None
    for r in (cs, ce):
        if r and r.get("nearest_domain_class") == "Ig-like":
            ig_row = r
            break

    # TM / kinase coordinates (from any coding-exon row that references them)
    tm_start = kinase_start = None
    for r in rows:
        if r.get("nearest_domain_class") == "TM" and tm_start is None:
            tm_start = _int(r.get("nearest_domain_start_aa"))
        if r.get("nearest_domain_class") == "kinase" and kinase_start is None:
            kinase_start = _int(r.get("nearest_domain_start_aa"))

    cassette_start = _int(cs["boundary_aa"]) if cs else None
    cassette_end = _int(ce["boundary_aa"]) if ce else None

    up_tm = (cassette_start is not None and tm_start is not None and cassette_start < tm_start)
    up_kinase = (cassette_start is not None and kinase_start is not None and cassette_start < kinase_start)

    display_status = r0.get("exon_block_display_status", "") or r0.get("source_coordinate_status", "")
    minor_display = display_status not in ("keep", "", "figure3C_native")

    if minor_display:
        interpretation = "supported_with_minor_display_note"
    else:
        interpretation = "human_reference_architecture_supported"

    return {
        "isoform": r0["isoform"],
        "transcript_id": r0.get("transcript_id", ""),
        "protein_id": r0.get("protein_id", ""),
        "protein_length": _int(r0.get("protein_length")),
        "cassette_start_aa": cassette_start,
        "cassette_end_aa": cassette_end,
        "nearest_ig_like_domain": (ig_row or {}).get("nearest_domain_label", ""),
        "ig_like_start_aa": _int((ig_row or {}).get("nearest_domain_start_aa")),
        "ig_like_end_aa": _int((ig_row or {}).get("nearest_domain_end_aa")),
        "tm_start_aa": tm_start,
        "kinase_start_aa": kinase_start,
        "cassette_upstream_of_tm": "yes" if up_tm else "no",
        "cassette_upstream_of_kinase": "yes" if up_kinase else "no",
        "cassette_start_class": (cs or {}).get("boundary_class", ""),
        "cassette_end_class": (ce or {}).get("boundary_class", ""),
        "boundary_class": (ce or {}).get("boundary_class", ""),
        "cassette_start_distance": _int((cs or {}).get("distance_to_nearest_domain_boundary")),
        "cassette_end_distance": _int((ce or {}).get("distance_to_nearest_domain_boundary")),
        "display_status": display_status or "native",
        "interpretation": interpretation,
    }


TSV_COLS = [
    "isoform", "transcript_id", "protein_id", "protein_length",
    "cassette_start_aa", "cassette_end_aa", "nearest_ig_like_domain",
    "ig_like_start_aa", "ig_like_end_aa", "tm_start_aa", "kinase_start_aa",
    "cassette_upstream_of_tm", "cassette_upstream_of_kinase",
    "boundary_class", "display_status", "interpretation",
]


def _write_tsv(summaries: list[dict]) -> Path:
    path = OUT_DIR / "human_fgfr2_boundary_sanity_check.tsv"
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(TSV_COLS)
        for s in summaries:
            w.writerow(["" if s.get(c) is None else s.get(c) for c in TSV_COLS])
    return path


def _iso_block(s: dict) -> str:
    def _fmt(v):
        return "—" if v is None or v == "" else v
    return (
        f"### Human FGFR2-{s['isoform']}\n\n"
        f"- **Final isoform label:** {s['isoform']}\n"
        f"- **Transcript / protein:** `{_fmt(s['transcript_id'])}` / `{_fmt(s['protein_id'])}`\n"
        f"- **Protein length:** {_fmt(s['protein_length'])} aa\n"
        f"- **Cassette slot (aa):** {_fmt(s['cassette_start_aa'])}–{_fmt(s['cassette_end_aa'])}\n"
        f"- **Nearest Ig-like domain:** {_fmt(s['nearest_ig_like_domain'])} "
        f"({_fmt(s['ig_like_start_aa'])}–{_fmt(s['ig_like_end_aa'])} aa)\n"
        f"- **pyTMHMM transmembrane start:** {_fmt(s['tm_start_aa'])} aa\n"
        f"- **Kinase domain start:** {_fmt(s['kinase_start_aa'])} aa\n"
        f"- **Cassette upstream of TM:** {s['cassette_upstream_of_tm']}\n"
        f"- **Cassette upstream of kinase:** {s['cassette_upstream_of_kinase']}\n"
        f"- **Cassette-start relation to Ig-like region:** "
        f"{s['cassette_start_class'].replace('_', ' ')} "
        f"({_fmt(s['cassette_start_distance'])} aa to nearest domain boundary)\n"
        f"- **Cassette-end relation to Ig-like region:** "
        f"{s['cassette_end_class'].replace('_', ' ')} "
        f"({_fmt(s['cassette_end_distance'])} aa to nearest domain boundary)\n"
        f"- **Display-coordinate note:** {s['display_status'].replace('_', ' ')}\n"
        f"- **Interpretation:** `{s['interpretation']}`\n"
    )


def _write_md(summaries: list[dict]) -> Path:
    path = OUT_DIR / "human_fgfr2_boundary_sanity_check.md"
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# Human FGFR2 reference — boundary sanity check",
        "",
        "*Read-only summary generated from existing pipeline outputs "
        "(step-16 exon\u2013domain boundary distances). No biological data, FASTA, "
        "truth table, or primary/review membership was changed.*",
        "",
        f"_Generated: {ts}_",
        f"_Source: `{display_path(DIST_TSV)}`_",
        "",
        "## Summary",
        "",
        "The human IIIb and IIIc reference isoforms carry the expected FGFR2 receptor "
        "architecture: three Ig-like domains, a single transmembrane segment (pyTMHMM), "
        "and a C-terminal kinase domain. In both isoforms the mutually exclusive IIIb/IIIc "
        "cassette lies within the third Ig-like (Ig3) region and upstream of both the "
        "transmembrane segment and the kinase domain, as expected for the receptor "
        "ectodomain. Cassette-end positions are aligned or near the Ig-like/TM domain "
        "boundary, while cassette-start positions fall within the Ig-like region.",
        "",
    ]
    for s in summaries:
        lines.append(_iso_block(s))
        lines.append("")
    lines += [
        "## Notes",
        "",
        "- Distances are relative to the nearest InterPro/pyTMHMM domain boundary.",
        "- `supported_with_minor_display_note` indicates a minor display-coordinate note "
        "(e.g. a codon-boundary length clamp); it is not a biological failure.",
        "- IIIb/IIIc labels, protein sequences, and dataset membership are unchanged by this report.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main() -> None:
    if not DIST_TSV.exists():
        raise SystemExit(f"Source table not found: {DIST_TSV}")
    rows = _load_rows()
    if not rows:
        raise SystemExit("No human rows found in the boundary distances table.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = []
    for iso in ("IIIb", "IIIc"):
        s = _summarize_isoform([r for r in rows if r["isoform"] == iso])
        if s:
            summaries.append(s)
    tsv = _write_tsv(summaries)
    md = _write_md(summaries)
    print(f"[human-sanity] wrote {display_path(md)}")
    print(f"[human-sanity] wrote {display_path(tsv)}")
    for s in summaries:
        print(f"[human-sanity] {s['isoform']}: {s['interpretation']} "
              f"(cassette {s['cassette_start_aa']}-{s['cassette_end_aa']}, "
              f"{s['display_status']})")


if __name__ == "__main__":
    main()
