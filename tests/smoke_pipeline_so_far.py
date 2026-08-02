#!/usr/bin/env python3
"""Synthetic end-to-end smoke test for the FGFR2 pipeline steps completed so far.

This does not contact Ensembl or NCBI. It uses a tiny artificial FGFR2-like data set
with two mutually exclusive internal exons so the local selection and IIIb/IIIc
classification logic can be checked quickly and reproducibly.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


def find_file(root: Path, name: str) -> Path:
    for candidate in (root / "scripts" / name, root / name):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Cannot find {name} in scripts/ or repository root")


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(cmd: list[str]) -> None:
    print("+", " ".join(str(x) for x in cmd))
    subprocess.run(cmd, check=True)


def make_exons(tx_id: str, alt_start: int, alt_end: int) -> list[dict[str, str]]:
    # Eight exons are used so the classifier does not treat the transcript as too short.
    coords = [
        (1, 100, 150),
        (2, 200, 250),
        (3, 300, 350),
        (4, alt_start, alt_end),       # mutually exclusive alternative exon
        (5, 800, 850),                 # identical right flank
        (6, 900, 950),
        (7, 1000, 1050),
        (8, 1100, 1150),
    ]
    return [
        {
            "transcript_id_internal": tx_id,
            "exon_rank": str(rank),
            "chrom": "chr10",
            "start": str(start),
            "end": str(end),
            "strand": "+",
        }
        for rank, start, end in coords
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    args = parser.parse_args()

    root = args.repo_root.resolve()
    work = root / "tmp" / "smoke_pipeline_so_far"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    py = sys.executable
    script_registry = find_file(root, "build_species_registry_improved.py")
    script_select = find_file(root, "select_fgfr2_transcripts_annotation_aware_v2.py")
    script_classify = find_file(root, "classify_fgfr2_IIIb_IIIc_by_exon_structure_v2.py")

    # Step 1: species registry
    species_list = work / "species_list.txt"
    species_list.write_text("Homo sapiens\nMus musculus\n", encoding="utf-8")
    registry_out = work / "01_registry"
    run([py, str(script_registry), "--species_list", str(species_list), "--outdir", str(registry_out), "--strict"])
    assert (registry_out / "species_registry.tsv").exists(), "species_registry.tsv was not created"

    # Synthetic outputs equivalent to what the collector would provide for downstream scripts.
    transcripts = [
        {
            "species_input": "Homo sapiens", "species_canonical": "Homo sapiens", "source_db": "synthetic",
            "gene_id_internal": "gene1", "gene_symbol": "FGFR2", "internal_transcript_id": "tx_ref",
            "transcript_id_source": "TX_REF", "transcript_name": "FGFR2 reference", "transcript_biotype": "protein_coding",
            "translation_id_source": "PROT_REF", "protein_length_aa": "821", "is_canonical_source": "1",
            "support_level": "tsl1", "completeness_flags": "complete",
        },
        {
            "species_input": "Homo sapiens", "species_canonical": "Homo sapiens", "source_db": "synthetic",
            "gene_id_internal": "gene1", "gene_symbol": "FGFR2", "internal_transcript_id": "tx_iiib",
            "transcript_id_source": "TX_IIIB", "transcript_name": "FGFR2 IIIb", "transcript_biotype": "protein_coding",
            "translation_id_source": "PROT_B", "protein_length_aa": "820", "is_canonical_source": "0",
            "support_level": "tsl1", "completeness_flags": "complete",
        },
        {
            "species_input": "Homo sapiens", "species_canonical": "Homo sapiens", "source_db": "synthetic",
            "gene_id_internal": "gene1", "gene_symbol": "FGFR2", "internal_transcript_id": "tx_iiic",
            "transcript_id_source": "TX_IIIC", "transcript_name": "FGFR2 IIIc", "transcript_biotype": "protein_coding",
            "translation_id_source": "PROT_C", "protein_length_aa": "820", "is_canonical_source": "0",
            "support_level": "tsl1", "completeness_flags": "complete",
        },
    ]
    tx_fields = list(transcripts[0].keys())
    transcript_path = work / "collector_like" / "transcripts.tsv"
    write_tsv(transcript_path, transcripts, tx_fields)

    exons = []
    exons.extend(make_exons("tx_ref", 400, 550))
    exons.extend(make_exons("tx_iiib", 400, 550))
    exons.extend(make_exons("tx_iiic", 600, 750))
    exon_fields = ["transcript_id_internal", "exon_rank", "chrom", "start", "end", "strand"]
    exon_path = work / "collector_like" / "exons.tsv"
    write_tsv(exon_path, exons, exon_fields)

    # Step 3 initial selection without structural isoform evidence.
    select_initial_out = work / "03_selection_initial"
    run([py, str(script_select), "--transcripts", str(transcript_path), "--exons", str(exon_path), "--outdir", str(select_initial_out), "--strict"])
    assert (select_initial_out / "selected_transcripts.tsv").exists(), "initial selected_transcripts.tsv missing"
    assert (select_initial_out / "run_metadata.json").exists(), "initial run_metadata.json missing"

    # Step 4 structural IIIb/IIIc classification.
    iso_out = work / "04_isoform_evidence"
    run([
        py, str(script_classify), "--transcripts", str(transcript_path), "--exons", str(exon_path),
        "--selected_transcripts", str(select_initial_out / "selected_transcripts.tsv"), "--outdir", str(iso_out), "--strict",
    ])
    evidence_path = iso_out / "fgfr2_isoform_evidence.tsv"
    assert evidence_path.exists(), "fgfr2_isoform_evidence.tsv missing"
    evidence = read_tsv(evidence_path)
    classes = {r["internal_transcript_id"]: r["isoform_class"] for r in evidence}
    assert classes.get("tx_iiib") == "IIIb", f"Expected tx_iiib as IIIb, got {classes.get('tx_iiib')}"
    assert classes.get("tx_iiic") == "IIIc", f"Expected tx_iiic as IIIc, got {classes.get('tx_iiic')}"

    # Step 5 selection repeated with structural isoform evidence.
    select_final_out = work / "05_selection_with_isoforms"
    run([
        py, str(script_select), "--transcripts", str(transcript_path), "--exons", str(exon_path),
        "--isoform_evidence", str(evidence_path), "--outdir", str(select_final_out), "--strict",
    ])
    selected = read_tsv(select_final_out / "selected_transcripts.tsv")
    roles = {r["selection_role"] for r in selected}
    assert "reference" in roles, f"reference role missing; roles={roles}"
    assert "FGFR2_IIIb_candidate" in roles, f"IIIb candidate missing; roles={roles}"
    assert "FGFR2_IIIc_candidate" in roles, f"IIIc candidate missing; roles={roles}"
    assert (select_final_out / "transcript_selection_report.md").exists(), "final Markdown report missing"
    assert (select_final_out / "transcript_selection_report.html").exists(), "final HTML report missing"

    print(f"\nSmoke pipeline passed. Outputs written to: {work}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
