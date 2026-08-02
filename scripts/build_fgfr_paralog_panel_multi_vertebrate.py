#!/usr/bin/env python3
"""
build_fgfr_paralog_panel_multi_vertebrate.py  (Addendum Task A)

Build a multi-vertebrate FGFR1/2/3/4 paralog reference panel for the FGFR2
orthology / paralog screen. Sequences are fetched with NCBI ``datasets`` for a
small, curated set of vertebrate groups (mammal, bird, amphibian, teleost fish,
reptile). The human-only panel remains the legacy control; this panel is the
preferred orthology evidence layer.

Outputs:
  references/fgfr_paralog_panel_multi_vertebrate.fasta
      headers: ">FGFR2|species=gallus_gallus|accession=NP_...|source=NCBI_RefSeq"
  <outdir>/fgfr2_paralog_reference_panel_manifest.tsv

A per-gene/per-species download cache avoids repeated network calls. If a
download fails, that entry is skipped and recorded as a manifest warning so the
panel is still usable (and reproducible) from whatever was retrievable.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SCRIPT_VERSION = "1.0"
GENES = ("FGFR1", "FGFR2", "FGFR3", "FGFR4")

# Curated representative vertebrate groups. (scientific_name, taxon_group)
DEFAULT_SPECIES: List[Tuple[str, str]] = [
    ("Homo sapiens", "mammal"),
    ("Mus musculus", "mammal"),
    ("Gallus gallus", "bird"),
    ("Xenopus tropicalis", "amphibian"),
    ("Danio rerio", "teleost_fish"),
    ("Anolis carolinensis", "reptile"),
]


def canon(name: str) -> str:
    return name.strip().lower().replace(" ", "_")


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    hdr = None
    seq: List[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(seq)))
            hdr = line[1:].strip()
            seq = []
        elif line.strip():
            seq.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(seq)))
    return out


def fetch_gene_protein(gene: str, species: str, datasets_bin: str, cache_dir: Path,
                       timeout: int) -> Tuple[Optional[str], Optional[str], str]:
    """Return (accession, sequence, note). Picks the longest protein as canonical."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"{gene}_{canon(species)}.faa"
    if cache.exists() and cache.stat().st_size > 0:
        recs = parse_fasta(cache.read_text(encoding="utf-8"))
        if recs:
            hdr, seq = max(recs, key=lambda r: len(r[1]))
            return hdr.split()[0], seq, "from_cache"
    with tempfile.TemporaryDirectory() as td:
        zip_path = Path(td) / "pkg.zip"
        cmd = [datasets_bin, "download", "gene", "symbol", gene, "--taxon", species,
               "--include", "protein", "--filename", str(zip_path)]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=timeout)
        except subprocess.CalledProcessError as exc:
            return None, None, f"datasets_error:{exc.returncode}"
        except subprocess.TimeoutExpired:
            return None, None, "datasets_timeout"
        if not zip_path.exists():
            return None, None, "no_zip_downloaded"
        try:
            with zipfile.ZipFile(zip_path) as zf:
                faa_names = [n for n in zf.namelist() if n.endswith("protein.faa")]
                if not faa_names:
                    return None, None, "no_protein_faa_in_package"
                text = zf.read(faa_names[0]).decode("utf-8", errors="replace")
        except zipfile.BadZipFile:
            return None, None, "bad_zip"
    recs = parse_fasta(text)
    if not recs:
        return None, None, "empty_protein_faa"
    cache.write_text("".join(f">{h}\n{s}\n" for h, s in recs), encoding="utf-8")
    hdr, seq = max(recs, key=lambda r: len(r[1]))
    return hdr.split()[0], seq, "downloaded"


def main() -> int:
    ap = argparse.ArgumentParser(description="Build multi-vertebrate FGFR1-4 paralog panel (Addendum A).")
    ap.add_argument("--panel_fasta", type=Path, default=Path("references/fgfr_paralog_panel_multi_vertebrate.fasta"))
    ap.add_argument("--outdir", type=Path, required=True, help="Where to write the panel manifest.")
    ap.add_argument("--cache_dir", type=Path, default=Path("references/_paralog_panel_cache"))
    ap.add_argument("--datasets_bin", default="datasets")
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--legacy_human_panel", type=Path, default=Path("references/human_FGFR1_2_3_4.fasta"),
                    help="Legacy human-only panel; merged in as the human_only fallback layer.")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.panel_fasta.parent.mkdir(parents=True, exist_ok=True)

    fasta_records: List[Tuple[str, str]] = []
    manifest_rows: List[Dict[str, object]] = []
    warnings: List[Dict[str, str]] = []

    for gene in GENES:
        for species, group in DEFAULT_SPECIES:
            acc, seq, note = fetch_gene_protein(gene, species, args.datasets_bin, args.cache_dir, args.timeout)
            if acc and seq:
                header = f"{gene}|species={canon(species)}|accession={acc}|source=NCBI_RefSeq"
                fasta_records.append((header, seq))
                manifest_rows.append({
                    "gene": gene, "species": canon(species), "taxon_group": group,
                    "accession": acc, "source": "NCBI_RefSeq", "sequence_length": len(seq),
                    "retrieval_note": note,
                })
            else:
                warnings.append({"gene": gene, "species": canon(species), "taxon_group": group,
                                 "warning": f"panel_entry_unavailable:{note}"})

    # Merge legacy human-only panel as an explicit human_only fallback layer.
    if args.legacy_human_panel.exists():
        for hdr, seq in parse_fasta(args.legacy_human_panel.read_text(encoding="utf-8")):
            gene = next((g for g in GENES if g in hdr.upper()), "FGFR")
            acc = hdr.split("_")[-1] if "_" in hdr else hdr
            header = f"{gene}|species=homo_sapiens|accession={acc}|source=legacy_human_panel"
            # Avoid duplicating if we already have a human RefSeq entry for this gene.
            if not any(r[0].startswith(f"{gene}|species=homo_sapiens|") and "NCBI_RefSeq" in r[0]
                       for r in fasta_records):
                fasta_records.append((header, seq))
                manifest_rows.append({
                    "gene": gene, "species": "homo_sapiens", "taxon_group": "mammal",
                    "accession": acc, "source": "legacy_human_panel", "sequence_length": len(seq),
                    "retrieval_note": "legacy_control",
                })

    # De-duplicate identical headers.
    seen = set()
    uniq: List[Tuple[str, str]] = []
    for hdr, seq in fasta_records:
        if hdr in seen:
            continue
        seen.add(hdr)
        uniq.append((hdr, seq))

    args.panel_fasta.write_text("".join(f">{h}\n{s}\n" for h, s in uniq), encoding="utf-8")

    with open(args.outdir / "fgfr2_paralog_reference_panel_manifest.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["gene", "species", "taxon_group", "accession", "source",
                                       "sequence_length", "retrieval_note"])
        w.writeheader()
        w.writerows(manifest_rows)

    with open(args.outdir / "fgfr2_paralog_reference_panel_warnings.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=["gene", "species", "taxon_group", "warning"])
        w.writeheader()
        w.writerows(warnings)

    from collections import Counter
    gene_counts = Counter(r["gene"] for r in manifest_rows)
    group_counts = Counter(r["taxon_group"] for r in manifest_rows)
    meta = {
        "script_version": SCRIPT_VERSION,
        "panel_fasta": str(args.panel_fasta),
        "n_panel_sequences": len(uniq),
        "gene_counts": dict(gene_counts),
        "taxon_group_counts": dict(group_counts),
        "n_warnings": len(warnings),
        "multi_vertebrate": len(group_counts) >= 3,
    }
    (args.outdir / "fgfr2_paralog_reference_panel_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"[OK] panel sequences: {len(uniq)}  genes={dict(gene_counts)}  groups={dict(group_counts)}  warnings={len(warnings)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
