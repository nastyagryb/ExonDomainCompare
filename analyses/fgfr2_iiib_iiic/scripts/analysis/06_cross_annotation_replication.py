
from __future__ import annotations

import argparse
import json
import logging
import re
import shutil
import subprocess
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from common import ensure_dir, read_tsv, setup_logging, write_json

ENSEMBL_BASE = "https://rest.ensembl.org"


def make_ensembl_session() -> requests.Session:
    retry = Retry(
        total=6,
        connect=6,
        read=6,
        status=6,
        backoff_factor=1.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "POST"}),
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    session = requests.Session()
    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=4,
        pool_maxsize=4,
    )
    session.mount("https://", adapter)
    session.headers.update({"Accept": "application/json"})
    return session


ENSEMBL_SESSION = make_ensembl_session()
ENSEMBL_TIMEOUT = (15, 180)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Acquire NCBI and Ensembl FGFR2 protein annotations, classify IIIb/IIIc candidates against curated cassette references, and prepare cross-annotation replication tables."
    )
    p.add_argument("--manifest", required=True, help="TSV: species, ncbi_taxon, ensembl_species, include")
    p.add_argument("--references", required=True, help="Curated human IIIb/IIIc cassette FASTA")
    p.add_argument("--outdir", required=True)
    p.add_argument("--sources", default="ncbi,ensembl")
    p.add_argument("--datasets-executable", default="datasets")
    p.add_argument("--sleep", type=float, default=0.15)
    p.add_argument("--skip-fetch", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


def ensembl_get(path: str, params: dict | None = None):
    response = ENSEMBL_SESSION.get(
        ENSEMBL_BASE + path,
        params=params,
        headers={"Content-Type": "application/json"},
        timeout=ENSEMBL_TIMEOUT,
    )
    response.raise_for_status()
    return response.json()


def ensembl_post(
    path: str,
    payload: dict,
    params: dict | None = None,
):
    response = ENSEMBL_SESSION.post(
        ENSEMBL_BASE + path,
        params=params,
        json=payload,
        headers={"Content-Type": "application/json"},
        timeout=ENSEMBL_TIMEOUT,
    )
    response.raise_for_status()
    return response.json()


def chunks(items: list[str], size: int = 50):
    for start in range(0, len(items), size):
        yield items[start:start + size]


def fetch_ensembl(
    species: str,
    outdir: Path,
    sleep: float,
) -> list[SeqRecord]:
    logging.info("Ensembl %s: FGFR2 lookup", species)

    data = ensembl_get(
        f"/lookup/symbol/{species}/FGFR2",
        {"expand": 1},
    )

    (outdir / "lookup.json").write_text(
        json.dumps(data, indent=2) + "\n"
    )

    transcript_by_translation = {}

    for transcript in data.get("Transcript", []):
        translation = transcript.get("Translation")

        if translation and translation.get("id"):
            transcript_by_translation[translation["id"]] = transcript

    translation_ids = sorted(transcript_by_translation)

    logging.info(
        "Ensembl %s: fetching %d translations in batches",
        species,
        len(translation_ids),
    )

    sequences = {}
    failures = []

    for batch_number, batch in enumerate(
        chunks(translation_ids, 50),
        start=1,
    ):
        try:
            result = ensembl_post(
                "/sequence/id",
                {"ids": batch},
                {"type": "protein"},
            )

            if isinstance(result, dict):
                result = [result]

            for item in result:
                sequence_id = item.get("id") or item.get("query")
                sequence = item.get("seq", "")

                if sequence_id and sequence:
                    sequences[sequence_id] = sequence

            logging.info(
                "Ensembl %s: batch %d complete (%d/%d)",
                species,
                batch_number,
                len(sequences),
                len(translation_ids),
            )

        except requests.RequestException as exc:
            logging.warning(
                "Batch request failed for %s: %s",
                species,
                exc,
            )
            logging.warning(
                "Falling back to individual translation requests"
            )

            for translation_id in batch:
                try:
                    item = ensembl_get(
                        f"/sequence/id/{translation_id}",
                        {"type": "protein"},
                    )

                    sequence = item.get("seq", "")

                    if sequence:
                        sequences[translation_id] = sequence
                    else:
                        failures.append({
                            "translation_id": translation_id,
                            "error": "empty_sequence",
                        })

                except requests.RequestException as single_exc:
                    failures.append({
                        "translation_id": translation_id,
                        "error": str(single_exc),
                    })

        if sleep:
            time.sleep(sleep)

    records = []
    metadata = []

    for translation_id in translation_ids:
        sequence = sequences.get(translation_id, "")

        if not sequence:
            continue

        transcript = transcript_by_translation[translation_id]

        records.append(
            SeqRecord(
                Seq(sequence),
                id=translation_id,
                description=(
                    f"transcript={transcript.get('id', '')}"
                ),
            )
        )

        metadata.append({
            "source": "ensembl",
            "ensembl_species": species,
            "transcript_id": transcript.get("id"),
            "protein_id": translation_id,
            "protein_length": len(sequence),
            "is_canonical": transcript.get("is_canonical"),
            "biotype": transcript.get("biotype"),
        })

    SeqIO.write(
        records,
        outdir / "protein.faa",
        "fasta",
    )

    pd.DataFrame(metadata).to_csv(
        outdir / "metadata.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        failures,
        columns=["translation_id", "error"],
    ).to_csv(
        outdir / "fetch_failures.tsv",
        sep="\t",
        index=False,
    )

    logging.info(
        "Ensembl %s: wrote %d protein sequences",
        species,
        len(records),
    )

    return records


def fetch_ncbi(taxon: str, outdir: Path, executable: str) -> list[SeqRecord]:
    exe = shutil.which(executable)
    if not exe:
        raise FileNotFoundError(
            f"NCBI Datasets executable '{executable}' not found. Install it or use --skip-fetch with a pre-downloaded protein.faa."
        )
    package = outdir / "ncbi_gene_package.zip"
    subprocess.run(
        [
            exe,
            "download",
            "gene",
            "symbol",
            "FGFR2",
            "--taxon",
            taxon,
            "--include",
            "gene,rna,cds,protein",
            "--filename",
            str(package),
        ],
        check=True,
    )
    with zipfile.ZipFile(package) as zf:
        zf.extractall(outdir / "package")
    candidates = list((outdir / "package").rglob("protein.faa"))
    if not candidates:
        raise FileNotFoundError(f"No protein.faa in NCBI package for {taxon}")
    records = list(SeqIO.parse(candidates[0], "fasta"))
    SeqIO.write(records, outdir / "protein.faa", "fasta")
    return records


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner(mode="local")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def score_reference(aligner: PairwiseAligner, protein: str, reference: str) -> tuple[float, int, int, float]:
    alignment = aligner.align(protein, reference)[0]
    target_blocks, query_blocks = alignment.aligned
    aligned_query = int(sum(end - start for start, end in query_blocks))
    if len(target_blocks):
        start = int(target_blocks[0][0]) + 1
        end = int(target_blocks[-1][1])
    else:
        start = end = -1
    coverage = aligned_query / len(reference) if reference else 0.0
    return float(alignment.score), start, end, coverage


def classify_records(records: list[SeqRecord], refs: dict[str, str], source: str, species: str) -> pd.DataFrame:
    aligner = make_aligner()
    rows = []
    for record in records:
        protein = str(record.seq).replace("*", "").upper()
        scores = {}
        for isoform, ref in refs.items():
            scores[isoform] = score_reference(aligner, protein, ref)
        ranked = sorted(scores.items(), key=lambda kv: kv[1][0], reverse=True)
        best_label, best = ranked[0]
        second = ranked[1][1]
        rows.append(
            {
                "source": source,
                "species": species,
                "protein_id": record.id,
                "protein_length": len(protein),
                "predicted_isoform": best_label,
                "best_local_alignment_score": best[0],
                "second_best_score": second[0],
                "score_margin": best[0] - second[0],
                "cassette_start_aa": best[1],
                "cassette_end_aa": best[2],
                "reference_coverage": best[3],
                "IIIb_score": scores["IIIb"][0],
                "IIIc_score": scores["IIIc"][0],
            }
        )
    return pd.DataFrame(rows)


def select_candidates(classified: pd.DataFrame) -> pd.DataFrame:
    if classified.empty:
        return classified
    ranked = classified.sort_values(
        ["species", "source", "predicted_isoform", "reference_coverage", "score_margin", "best_local_alignment_score", "protein_length"],
        ascending=[True, True, True, False, False, False, False],
    )
    return ranked.groupby(["species", "source", "predicted_isoform"], as_index=False).head(1)


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)
    out = ensure_dir(args.outdir)
    manifest = read_tsv(args.manifest, required=["species", "ncbi_taxon", "ensembl_species"])
    if "include" in manifest.columns:
        manifest = manifest[manifest["include"].astype(str).str.lower().isin(["1", "true", "yes", "y"])]
    ref_records = list(SeqIO.parse(args.references, "fasta"))
    refs = {}
    for r in ref_records:
        label = "IIIb" if "iiib" in r.id.lower() or "iiib" in r.description.lower() else "IIIc" if "iiic" in r.id.lower() or "iiic" in r.description.lower() else None
        if label:
            refs[label] = str(r.seq).replace("-", "").upper()
    if set(refs) != {"IIIb", "IIIc"}:
        raise ValueError("Reference FASTA must contain identifiable IIIb and IIIc records")

    sources = {x.strip().lower() for x in args.sources.split(",")}
    classified_all = []
    fasta_index = {}
    for _, row in manifest.iterrows():
        species = row["species"]
        if "ensembl" in sources:
            d = ensure_dir(out / "raw" / safe(species) / "ensembl")
            protein_file = d / "protein.faa"
            records = (list(SeqIO.parse(protein_file, "fasta")) if protein_file.exists() and protein_file.stat().st_size > 0 else fetch_ensembl(row["ensembl_species"], d, args.sleep))
            fasta_index[(species, "ensembl")] = {r.id: r for r in records}
            classified_all.append(classify_records(records, refs, "ensembl", species))
        if "ncbi" in sources:
            d = ensure_dir(out / "raw" / safe(species) / "ncbi")
            protein_file = d / "protein.faa"
            records = (list(SeqIO.parse(protein_file, "fasta")) if protein_file.exists() and protein_file.stat().st_size > 0 else fetch_ncbi(row["ncbi_taxon"], d, args.datasets_executable))
            fasta_index[(species, "ncbi")] = {r.id: r for r in records}
            classified_all.append(classify_records(records, refs, "ncbi", species))

    classified = pd.concat(classified_all, ignore_index=True) if classified_all else pd.DataFrame()
    classified.to_csv(out / "all_annotation_candidates_classified.tsv", sep="\t", index=False)
    selected = select_candidates(classified)
    selected.to_csv(out / "selected_cross_annotation_candidates.tsv", sep="\t", index=False)

    selected_records = []
    for _, row in selected.iterrows():
        rec = fasta_index[(row["species"], row["source"])][row["protein_id"]]
        new_id = f"{safe(row['species'])}|{row['source']}|{row['predicted_isoform']}|{row['protein_id']}"
        selected_records.append(SeqRecord(rec.seq, id=new_id, description=""))
    SeqIO.write(selected_records, out / "selected_cross_annotation_proteins_for_interpro.faa", "fasta")

    comparison_rows = []
    for (species, isoform), group in selected.groupby(["species", "predicted_isoform"]):
        if set(group["source"]) >= {"ncbi", "ensembl"}:
            n = group[group["source"] == "ncbi"].iloc[0]
            e = group[group["source"] == "ensembl"].iloc[0]
            comparison_rows.append(
                {
                    "species": species,
                    "isoform": isoform,
                    "ncbi_protein_id": n["protein_id"],
                    "ensembl_protein_id": e["protein_id"],
                    "protein_length_delta_ncbi_minus_ensembl": n["protein_length"] - e["protein_length"],
                    "cassette_start_delta_ncbi_minus_ensembl": n["cassette_start_aa"] - e["cassette_start_aa"],
                    "cassette_end_delta_ncbi_minus_ensembl": n["cassette_end_aa"] - e["cassette_end_aa"],
                    "same_isoform_assignment": n["predicted_isoform"] == e["predicted_isoform"],
                    "minimum_reference_coverage": min(n["reference_coverage"], e["reference_coverage"]),
                    "minimum_score_margin": min(n["score_margin"], e["score_margin"]),
                }
            )
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(out / "ncbi_vs_ensembl_coordinate_comparison.tsv", sep="\t", index=False)
    write_json(
        {
            "n_species_requested": int(manifest["species"].nunique()),
            "sources": sorted(sources),
            "n_candidates_classified": int(len(classified)),
            "n_selected_candidates": int(len(selected)),
            "n_ncbi_ensembl_pairs": int(len(comparison)),
            "candidate_selection": {
                "grouping": ["species", "source", "predicted_isoform"],
                "ranking": ["reference_coverage", "score_margin", "best_local_alignment_score", "protein_length"],
                "direction": "descending for every ranking field",
                "identity_times_coverage_used": False,
            },
            "next_step": "Run InterProScan on selected_cross_annotation_proteins_for_interpro.faa, then compare D3 topology using the same annotation-aware boundary rules.",
        },
        out / "cross_annotation_manifest.json",
    )


if __name__ == "__main__":
    main()
