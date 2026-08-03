#!/usr/bin/env python3
"""
collect_fgfr2_models_dual_source_v3.py

Improved dual-source FGFR2 model collector.

Main improvements over earlier versions
---------------------------------------
1. Separates model status from source-comparison conflict status.
2. Uses deterministic NCBI Datasets workflow:
   summary genome taxon -> choose best assembly -> download selected accession.
3. Treats chromosome/seqid naming differences as minor label mismatches, not as
   automatic moderate conflicts.
4. Only computes gene_overlap_fraction when coordinate systems are plausibly comparable.
5. Produces richer source_comparison.tsv and selection_decisions.tsv outputs.
6. Writes compact ncbi_assembly_selected.tsv in addition to the full selection table.
7. Adds internal_consistency_checks.tsv for end-of-run structural validation.
8. Performs stricter validation of NCBI GFF3-derived models.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import time
import zipfile
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import sys as _sys

if str(Path(__file__).resolve().parent) not in _sys.path:
    _sys.path.insert(0, str(Path(__file__).resolve().parent))

from exondomaincompare.shared_gene_analysis.strand import MINUS, is_reverse, normalize_strand  # noqa: E402
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode, unquote
from urllib.request import Request, urlopen

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.shared_gene_analysis import assembly_selection as asel  # noqa: E402
from exondomaincompare.shared_gene_analysis import gene_identification as gid  # noqa: E402
from exondomaincompare.shared_gene_analysis import model_recovery as recovery  # noqa: E402

ENSEMBL_REST = "https://rest.ensembl.org"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


# ---------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------
@dataclass
class SpeciesRecord:
    input_name: str
    ensembl_species: str
    ncbi_species: str
    taxid: str
    common_name: str = ""
    preferred_source: str = "ensembl_first"
    assembly_preference: str = "RefSeq"


@dataclass
class GeneModel:
    source_db: str
    species_input: str
    species_canonical: str
    taxid: str
    assembly_accession: str
    assembly_name: str
    gene_symbol_requested: str
    gene_symbol_found: str
    gene_id_source: str
    gene_biotype: str
    chrom: str
    start: str
    end: str
    strand: str
    model_status: str
    model_confidence: str
    lookup_method: str = ""
    internal_gene_id: str = ""


@dataclass
class TranscriptModel:
    source_db: str
    species_canonical: str
    gene_id_internal: str
    transcript_id_source: str
    transcript_name: str
    transcript_biotype: str
    translation_id_source: str
    protein_length_aa: str
    protein_length_source: str
    is_canonical_source: str
    support_level: str
    completeness_flags: str
    transcript_model_confidence: str
    chrom: str
    start: str
    end: str
    strand: str
    internal_transcript_id: str = ""


@dataclass
class ExonModel:
    source_db: str
    species_canonical: str
    transcript_id_internal: str
    exon_id_source: str
    exon_rank: str
    chrom: str
    start: str
    end: str
    strand: str
    phase: str
    end_phase: str
    internal_exon_id: str = ""


@dataclass
class CDSFeature:
    source_db: str
    species_input: str
    species_canonical: str
    transcript_id_internal: str
    transcript_id_source: str
    translation_id_source: str
    cds_id_source: str
    cds_rank: str
    chrom: str
    start: str
    end: str
    strand: str
    phase: str
    cds_length_bp: str
    cds_offset_start_0based: str
    cds_offset_end_0based: str
    protein_start_aa: str
    protein_end_aa: str
    parent_feature_id: str
    internal_cds_id: str = ""
    coordinate_source: str = ""
    confidence: str = ""
    warning: str = ""
    # How this part came to sit where it does, recorded rather than inferred. A reader (or a
    # regression test) can see whether the transcript order was taken from the annotation's
    # own ranks or derived from the normalised strand, and the source spelling is kept beside
    # the normalised one so a spelling difference is visible instead of silent.
    normalized_strand: str = ""
    genomic_order: str = ""
    transcript_order: str = ""
    coding_exon_order: str = ""
    source_ordering_method: str = ""


# ---------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------
def to_str(x) -> str:
    return "" if x is None else str(x)


def safe_sleep(seconds: float) -> None:
    if seconds > 0:
        time.sleep(seconds)


def strip_ensembl_version(value: str) -> str:
    value = to_str(value).strip()
    return value.split(".", 1)[0] if value else ""


def to_int(value: str) -> Optional[int]:
    try:
        if value is None or str(value).strip() == "":
            return None
        return int(value)
    except Exception:
        return None


def safe_float_str(x: Optional[float]) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    return f"{x:.6f}"


def read_species_registry(path: Path) -> List[SpeciesRecord]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Species registry is empty.")

    required = {"input_name", "ensembl_species", "ncbi_species", "taxid"}
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Species registry missing required columns: {sorted(missing)}")

    out = []
    for row in rows:
        out.append(
            SpeciesRecord(
                input_name=row["input_name"].strip(),
                ensembl_species=row["ensembl_species"].strip(),
                ncbi_species=row["ncbi_species"].strip(),
                taxid=row["taxid"].strip(),
                common_name=row.get("common_name", "").strip(),
                preferred_source=row.get("preferred_source", "ensembl_first").strip() or "ensembl_first",
                assembly_preference=row.get("assembly_preference", "RefSeq").strip() or "RefSeq",
            )
        )
    return out


def write_tsv(rows: List[Dict[str, str]], path: Path, fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def write_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


# ---------------------------------------------------------------------
# HTTP / JSON
# ---------------------------------------------------------------------
def fetch_json(url: str, headers: Optional[Dict[str, str]] = None, sleep_s: float = 0.25, retries: int = 2) -> Dict:
    """Fetch JSON with a small retry loop.

    Failures are re-raised after retries so caller functions can record
    per-species warnings and continue.
    """
    req_headers = {"User-Agent": "Mozilla/5.0"}
    if headers:
        req_headers.update(headers)
    last_exc = None
    for attempt in range(retries + 1):
        try:
            req = Request(url, headers=req_headers)
            with urlopen(req, timeout=60) as resp:
                data = resp.read().decode("utf-8")
            safe_sleep(sleep_s)
            return json.loads(data)
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_exc = exc
            if attempt < retries:
                safe_sleep(min(2.0, sleep_s + 0.5 * (attempt + 1)))
                continue
            raise last_exc
    raise RuntimeError("unreachable fetch_json retry state")


def ensembl_get_json(path: str, params: Optional[Dict[str, str]] = None, sleep_s: float = 0.25) -> Dict:
    url = ENSEMBL_REST + path
    if params:
        url = f"{url}?{urlencode(params)}"
    return fetch_json(
        url,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        sleep_s=sleep_s,
    )


def esearch(db: str, term: str, retmax: int = 20, sleep_s: float = 0.34) -> List[str]:
    url = EUTILS_BASE + "esearch.fcgi?" + urlencode(
        {"db": db, "term": term, "retmode": "json", "retmax": str(retmax)}
    )
    data = fetch_json(url, sleep_s=sleep_s)
    return data.get("esearchresult", {}).get("idlist", [])


def esearch_safe(db: str, term: str, retmax: int = 20, sleep_s: float = 0.34) -> Tuple[List[str], str]:
    try:
        return esearch(db, term, retmax=retmax, sleep_s=sleep_s), "ok"
    except Exception as exc:
        return [], f"failed:{type(exc).__name__}"


def esummary(db: str, ids: List[str], sleep_s: float = 0.34) -> Dict:
    if not ids:
        return {}
    url = EUTILS_BASE + "esummary.fcgi?" + urlencode(
        {"db": db, "id": ",".join(ids), "retmode": "json"}
    )
    return fetch_json(url, sleep_s=sleep_s)


def esummary_safe(db: str, ids: List[str], sleep_s: float = 0.34) -> Tuple[Dict, str]:
    try:
        return esummary(db, ids, sleep_s=sleep_s), "ok"
    except Exception as exc:
        return {}, f"failed:{type(exc).__name__}"




def build_cds_features_from_parts(
    *,
    source_db: str,
    species_input: str,
    species_canonical: str,
    transcript_id_internal: str,
    transcript_id_source: str,
    translation_id_source: str,
    parts: List[Dict[str, str]],
    coordinate_source: str,
    protein_length_aa: Optional[int] = None,
) -> List[CDSFeature]:
    """Build an auditable CDS feature table and AA intervals from ordered CDS parts.

    CDS parts must already be true CDS intervals, not whole exons. The function puts them in
    transcript 5'→3' order, accumulates CDS offsets, and converts those offsets to 1-based
    protein AA intervals. Codons split across CDS-part boundaries are represented by inclusive
    AA intervals; downstream scripts can flag phase/split-codon cases rather than pretending
    the part is a clean codon block.

    Two properties are worth stating because getting either wrong is silent.

    *Order.* Transcript order comes from the annotation's own ranks when the parts carry them,
    and only otherwise from the normalised strand. Ranks are what the annotation asserts; a
    coordinate sort is a reconstruction of it. Deriving the order from a raw strand comparison
    is what reversed every Ensembl minus-strand transcript, because Ensembl spells that strand
    ``-1`` and the comparison tested for ``-``. Parts already in transcript order are also
    left alone rather than sorted twice.

    *The stop codon.* An annotation's CDS features usually include the terminator, so the
    nucleotide length is three more than the protein needs. Projected without allowance, the
    last three bases become a residue the protein does not have — a phantom final position
    one past the end, and, where a source splits the stop across two CDS features, a
    degenerate trailing part covering that phantom alone. Pass ``protein_length_aa`` (from the
    translated sequence) and the projection is truncated to the real protein; parts that fall
    entirely in the terminator are marked and carry no protein interval.
    """
    clean = []
    for r in parts or []:
        st, en = to_int(r.get("start", "")), to_int(r.get("end", ""))
        if st is None or en is None or st > en:
            continue
        rr = dict(r)
        rr["start"] = str(st)
        rr["end"] = str(en)
        clean.append(rr)
    if not clean:
        return []

    strand_value = next((r.get("strand") for r in clean if r.get("strand") not in (None, "")),
                        "+")
    normalized = normalize_strand(strand_value)
    clean, ordering_method = _in_transcript_order(clean, normalized)
    genomic_rank = {id(r): i + 1 for i, r in
                    enumerate(sorted(clean, key=lambda x: int(x["start"])))}

    out: List[CDSFeature] = []
    cds_offset = 0
    coding_rank = 0
    for rank, r in enumerate(clean, start=1):
        st, en = int(r["start"]), int(r["end"])
        length = en - st + 1
        offset_start = cds_offset
        offset_end = cds_offset + length - 1
        aa_start = offset_start // 3 + 1
        aa_end = offset_end // 3 + 1
        cds_id = r.get("cds_id_source", "") or r.get("exon_id", "") or r.get("ID", "") or f"{transcript_id_source}_CDS_{rank}"
        warning = ""
        if length % 3 != 0 or offset_start % 3 != 0 or (offset_end + 1) % 3 != 0:
            warning = "CDS feature boundary may split codons; AA interval is inclusive"

        # Truncate at the real protein. Everything past it is the terminator.
        terminator_only = False
        if protein_length_aa and protein_length_aa > 0:
            if aa_start > protein_length_aa:
                terminator_only = True
                warning = ("CDS part encodes only the translation terminator; it has no "
                           "protein interval")
            elif aa_end > protein_length_aa:
                aa_end = protein_length_aa
        if not terminator_only:
            coding_rank += 1
        out.append(CDSFeature(
            source_db=source_db, species_input=species_input, species_canonical=species_canonical,
            transcript_id_internal=transcript_id_internal, transcript_id_source=transcript_id_source,
            translation_id_source=translation_id_source, cds_id_source=cds_id, cds_rank=str(rank),
            chrom=r.get("seqid", r.get("chrom", "")), start=str(st), end=str(en),
            strand=r.get("strand", strand_value),
            phase=r.get("phase", ""), cds_length_bp=str(length),
            cds_offset_start_0based=str(offset_start), cds_offset_end_0based=str(offset_end),
            protein_start_aa="" if terminator_only else str(aa_start),
            protein_end_aa="" if terminator_only else str(aa_end),
            parent_feature_id=r.get("parent_feature_id", transcript_id_source),
            internal_cds_id=f"{transcript_id_internal}|cds|{rank}|{cds_id}",
            coordinate_source=coordinate_source, confidence="exact_from_CDS_features",
            warning=warning,
            normalized_strand="" if normalized is None else str(normalized),
            genomic_order=str(genomic_rank[id(r)]),
            transcript_order=str(rank),
            coding_exon_order="" if terminator_only else str(coding_rank),
            source_ordering_method=ordering_method,
        ))
        cds_offset += length
    return out


def _in_transcript_order(parts: List[Dict[str, str]],
                         normalized_strand: Optional[int]) -> Tuple[List[Dict[str, str]], str]:
    """CDS parts in transcript 5'→3' order, plus the name of the method that got them there.

    The annotation's own ranks win when present, because they are the assertion and a
    coordinate sort is only a reconstruction of it. Parts that already ascend or descend
    consistently with the strand are returned untouched, so a table produced in transcript
    order is not reversed a second time.
    """
    ranks = [to_int(r.get("exon_rank") or r.get("cds_rank") or r.get("rank") or "")
             for r in parts]
    if all(v is not None for v in ranks) and len(set(ranks)) == len(ranks):
        return ([p for _, p in sorted(zip(ranks, parts), key=lambda pair: pair[0])],
                "annotation_rank")

    if normalized_strand is None:
        # Unknown strand: ascending coordinates, said out loud rather than assumed correct.
        return (sorted(parts, key=lambda x: int(x["start"])), "genomic_ascending_strand_unknown")

    starts = [int(p["start"]) for p in parts]
    descending = starts == sorted(starts, reverse=True)
    ascending = starts == sorted(starts)
    if normalized_strand == MINUS:
        if descending:
            return parts, "already_transcript_order_minus_strand"
        return sorted(parts, key=lambda x: int(x["start"]), reverse=True), "normalized_strand"
    if ascending:
        return parts, "already_transcript_order_plus_strand"
    return sorted(parts, key=lambda x: int(x["start"])), "normalized_strand"


def build_ensembl_cds_features_from_expanded_transcript(
    species: SpeciesRecord, gene_internal: str, tx: Dict, tx_id: str, translation_id: str
) -> List[CDSFeature]:
    """Derive CDS pieces for an Ensembl expanded transcript when Translation start/end are present.

    Ensembl expanded lookup returns exon intervals and, for translated transcripts,
    a Translation object. When genomic translation start/end are available, the
    CDS feature pieces are the intersections of transcript exons with that interval.
    If these coordinates are absent, no CDS rows are emitted rather than guessing.
    """
    translation = tx.get("Translation") or {}
    tstart, tend = to_int(translation.get("start", "")), to_int(translation.get("end", ""))
    if tstart is None or tend is None:
        return []
    lo, hi = min(tstart, tend), max(tstart, tend)
    _strand = to_str(tx.get("strand", "")) or "+"
    parts=[]
    for rank, ex in enumerate(tx.get("Exon", []) or [], start=1):
        est, een = to_int(ex.get("start", "")), to_int(ex.get("end", ""))
        if est is None or een is None:
            continue
        ov_s, ov_e = max(est, lo), min(een, hi)
        if ov_s <= ov_e:
            parts.append({
                "cds_id_source": strip_ensembl_version(ex.get("id")) or f"{tx_id}_CDS_{rank}",
                "seqid": to_str(ex.get("seq_region_name", tx.get("seq_region_name", ""))),
                "start": str(ov_s), "end": str(ov_e),
                "strand": to_str(ex.get("strand", tx.get("strand", ""))),
                "phase": to_str(ex.get("phase", "")),
                "parent_feature_id": tx_id,
            })
    tx_internal = f"{gene_internal}|tx|{tx_id}"
    return build_cds_features_from_parts(
        source_db="Ensembl", species_input=species.input_name, species_canonical=species.ensembl_species,
        transcript_id_internal=tx_internal, transcript_id_source=tx_id, translation_id_source=translation_id,
        parts=parts, coordinate_source="Ensembl_expanded_lookup_translation_exon_intersection",
    )

# ---------------------------------------------------------------------
# Ensembl adapter
# ---------------------------------------------------------------------
def try_ensembl_lookup(species: SpeciesRecord, gene_symbol: str, sleep_s: float) -> Tuple[Optional[Dict], str]:
    try:
        data = ensembl_get_json(
            f"/lookup/symbol/{quote(species.ensembl_species)}/{quote(gene_symbol)}",
            params={"expand": 1},
            sleep_s=sleep_s,
        )
        if isinstance(data, dict) and data.get("id"):
            return data, "symbol_lookup"
        return None, "symbol_lookup_empty"
    except (HTTPError, URLError, TimeoutError) as exc:
        return None, f"symbol_lookup_failed:{type(exc).__name__}"


def try_ensembl_xref_then_lookup(species: SpeciesRecord, gene_symbol: str, sleep_s: float) -> Tuple[Optional[Dict], str]:
    try:
        xrefs = ensembl_get_json(
            f"/xrefs/symbol/{quote(species.ensembl_species)}/{quote(gene_symbol)}",
            sleep_s=sleep_s,
        )
        if not isinstance(xrefs, list) or not xrefs:
            return None, "xref_empty"
        gene_id = next((rec.get("id") for rec in xrefs if rec.get("type") == "gene"), None) or xrefs[0].get("id")
        if not gene_id:
            return None, "xref_no_gene_id"
        data = ensembl_get_json(
            f"/lookup/id/{quote(gene_id)}",
            params={"expand": 1},
            sleep_s=sleep_s,
        )
        if isinstance(data, dict) and data.get("id"):
            return data, "xref_lookup"
        return None, "xref_lookup_empty"
    except (HTTPError, URLError, TimeoutError) as exc:
        return None, f"xref_lookup_failed:{type(exc).__name__}"


def fetch_ensembl_model(species: SpeciesRecord, gene_symbol: str, sleep_s: float) -> Tuple[Optional[Dict], str]:
    data, mode = try_ensembl_lookup(species, gene_symbol, sleep_s)
    if data is not None:
        return data, mode
    data2, mode2 = try_ensembl_xref_then_lookup(species, gene_symbol, sleep_s)
    if data2 is not None:
        return data2, mode2
    return None, f"{mode}|{mode2}"


def validate_ensembl_model(gene: Dict) -> Tuple[str, List[str]]:
    warnings: List[str] = []
    if gene.get("biotype") != "protein_coding":
        warnings.append(f"gene biotype is {gene.get('biotype', '')}, not protein_coding")

    transcripts = gene.get("Transcript", []) or []
    if not transcripts:
        return "gene_found_but_no_transcripts", warnings + ["no transcripts returned"]

    protein_coding_txs = [tx for tx in transcripts if tx.get("biotype") == "protein_coding"]
    if not protein_coding_txs:
        warnings.append("no protein_coding transcripts")

    with_translation = [tx for tx in protein_coding_txs if (tx.get("Translation") or {}).get("id")]
    if not with_translation:
        warnings.append("no translated protein_coding transcripts")

    exons_present = any((tx.get("Exon") or []) for tx in transcripts)
    if not exons_present:
        warnings.append("no exons returned")

    if with_translation and exons_present and gene.get("biotype") == "protein_coding":
        return "model_found_high_confidence", warnings
    if protein_coding_txs:
        return "model_found_medium_confidence", warnings
    return "model_found_partial", warnings


def normalize_ensembl_model(
    species: SpeciesRecord,
    gene_symbol: str,
    gene: Dict,
    lookup_method: str,
) -> Tuple[GeneModel, List[TranscriptModel], List[ExonModel], List[CDSFeature], List[str]]:
    model_status, warnings = validate_ensembl_model(gene)
    model_confidence = (
        "high" if model_status == "model_found_high_confidence"
        else "medium" if model_status == "model_found_medium_confidence"
        else "low"
    )

    canonical_transcript = strip_ensembl_version(gene.get("canonical_transcript"))
    gene_internal = f"{species.ensembl_species}|Ensembl|{strip_ensembl_version(gene.get('id'))}"

    gene_model = GeneModel(
        source_db="Ensembl",
        species_input=species.input_name,
        species_canonical=species.ensembl_species,
        taxid=species.taxid,
        assembly_accession="",
        assembly_name=to_str(gene.get("assembly_name")),
        gene_symbol_requested=gene_symbol,
        gene_symbol_found=to_str(gene.get("display_name")) or gene_symbol,
        gene_id_source=strip_ensembl_version(gene.get("id")),
        gene_biotype=to_str(gene.get("biotype")),
        chrom=to_str(gene.get("seq_region_name")),
        start=to_str(gene.get("start")),
        end=to_str(gene.get("end")),
        strand=to_str(gene.get("strand")),
        model_status=model_status,
        model_confidence=model_confidence,
        lookup_method=lookup_method,
        internal_gene_id=gene_internal,
    )

    txs: List[TranscriptModel] = []
    exs: List[ExonModel] = []
    cds_features: List[CDSFeature] = []
    for tx in gene.get("Transcript", []) or []:
        tx_id = strip_ensembl_version(tx.get("id"))
        translation = tx.get("Translation") or {}
        translation_id = strip_ensembl_version(translation.get("id"))
        tx_internal = f"{gene_internal}|tx|{tx_id}"

        completeness_flags: List[str] = []
        if not translation_id:
            completeness_flags.append("no_translation")
        if not (tx.get("Exon") or []):
            completeness_flags.append("no_exons")

        tx_conf = "high" if translation_id and (tx.get("Exon") or []) else "medium"

        txs.append(
            TranscriptModel(
                source_db="Ensembl",
                species_canonical=species.ensembl_species,
                gene_id_internal=gene_internal,
                transcript_id_source=tx_id,
                transcript_name=to_str(tx.get("display_name")),
                transcript_biotype=to_str(tx.get("biotype")),
                translation_id_source=translation_id,
                protein_length_aa=to_str(translation.get("length")),
                protein_length_source="Ensembl_translation_length" if to_str(translation.get("length")) else "missing",
                is_canonical_source="1" if tx_id and tx_id == canonical_transcript else "0",
                support_level="",
                completeness_flags=";".join(completeness_flags),
                transcript_model_confidence=tx_conf,
                chrom=to_str(tx.get("seq_region_name", gene.get("seq_region_name"))),
                start=to_str(tx.get("start")),
                end=to_str(tx.get("end")),
                strand=to_str(tx.get("strand", gene.get("strand"))),
                internal_transcript_id=tx_internal,
            )
        )

        for rank, ex in enumerate(tx.get("Exon", []) or [], start=1):
            exon_id = strip_ensembl_version(ex.get("id"))
            exs.append(
                ExonModel(
                    source_db="Ensembl",
                    species_canonical=species.ensembl_species,
                    transcript_id_internal=tx_internal,
                    exon_id_source=exon_id,
                    exon_rank=str(rank),
                    chrom=to_str(ex.get("seq_region_name", tx.get("seq_region_name", gene.get("seq_region_name")))),
                    start=to_str(ex.get("start")),
                    end=to_str(ex.get("end")),
                    strand=to_str(ex.get("strand", tx.get("strand", gene.get("strand")))),
                    phase=to_str(ex.get("phase")),
                    end_phase=to_str(ex.get("end_phase")),
                    internal_exon_id=f"{tx_internal}|exon|{rank}|{exon_id}",
                )
            )
        cds_features.extend(build_ensembl_cds_features_from_expanded_transcript(species, gene_internal, tx, tx_id, translation_id))

    return gene_model, txs, exs, cds_features, warnings


# ---------------------------------------------------------------------
# NCBI evidence adapter
# ---------------------------------------------------------------------
def fetch_ncbi_gene_evidence(species: SpeciesRecord, gene_symbol: str, sleep_s: float) -> Dict[str, str]:
    term = f'{gene_symbol}[Gene Name] AND "{species.ncbi_species}"[Organism]'
    ids, search_status = esearch_safe("gene", term, retmax=5, sleep_s=sleep_s)
    if not ids:
        return {
            "ncbi_gene_found": "0",
            "ncbi_gene_id": "",
            "ncbi_gene_name": "",
            "ncbi_gene_description": "",
            "ncbi_gene_taxname": "",
            "ncbi_eutils_status": search_status,
        }
    summary, summary_status = esummary_safe("gene", ids, sleep_s=sleep_s)
    result = summary.get("result", {})
    uids = result.get("uids", [])
    if not uids:
        return {
            "ncbi_gene_found": "0",
            "ncbi_gene_id": "",
            "ncbi_gene_name": "",
            "ncbi_gene_description": "",
            "ncbi_gene_taxname": "",
            "ncbi_eutils_status": summary_status,
        }
    uid = uids[0]
    rec = result.get(uid, {})
    return {
        "ncbi_gene_found": "1",
        "ncbi_gene_id": to_str(uid),
        "ncbi_gene_name": to_str(rec.get("name")),
        "ncbi_gene_description": to_str(rec.get("description")),
        "ncbi_gene_taxname": to_str((rec.get("organism") or {}).get("scientificname")),
        "ncbi_eutils_status": summary_status,
    }


def fetch_ncbi_refseq_protein_evidence(species: SpeciesRecord, gene_symbol: str, sleep_s: float) -> Dict[str, str]:
    term = f'{gene_symbol}[Gene] AND "{species.ncbi_species}"[Organism] AND srcdb_refseq[PROP]'
    ids, status = esearch_safe("protein", term, retmax=10, sleep_s=sleep_s)
    return {
        "refseq_protein_hits": str(len(ids)),
        "refseq_protein_example_ids": ",".join(ids[:5]),
        "ncbi_protein_eutils_status": status,
    }


# ---------------------------------------------------------------------
# NCBI Datasets helpers
# ---------------------------------------------------------------------
def gene_product_name(gene_config: Optional[Path], gene_symbol: str) -> str:
    """The spelled-out product name from the gene config, as a fallback.

    Preferred source is the annotation service's own description of the requested
    gene, which needs no configuration. This is the offline fallback, so a run whose
    E-utils lookup was unavailable can still tell a description of this gene's product
    from a description of a paralog's.
    """
    if not gene_config:
        return ""
    try:
        text = Path(gene_config).read_text(encoding="utf-8")
    except OSError:
        return ""
    match = re.search(r"^\s*display_name:\s*(.+)$", text, re.MULTILINE)
    while match:
        value = match.group(1).strip().strip('"\'')
        # The gene block's display_name is the product name; other blocks have their
        # own. Accept only a value that ends in this symbol's number.
        if value.lower().endswith(gid.paralog_family(gene_symbol)[1]):
            return value
        text = text[match.end():]
        match = re.search(r"^\s*display_name:\s*(.+)$", text, re.MULTILINE)
    return ""


def datasets_available(user_bin: Optional[str]) -> Optional[str]:
    if user_bin:
        # accept either an explicit path or a bare command name resolvable on PATH
        if Path(user_bin).exists():
            return user_bin
        return shutil.which(user_bin)
    return shutil.which("datasets")


def run_command(cmd: List[str], cwd: Optional[Path] = None, timeout: Optional[int] = None) -> Tuple[int, str]:
    proc = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout,
    )
    return proc.returncode, proc.stdout




def pick_best_ncbi_assembly_from_summary(candidates: List[Dict[str, str]], assembly_preference: str) -> Tuple[Optional[Dict[str, str]], List[Dict[str, str]]]:
    if not candidates:
        return None, []

    level_rank = {
        "complete genome": 40,
        "chromosome": 35,
        "scaffold": 20,
        "contig": 10,
    }

    scored_rows = []
    for row in candidates:
        accession = row.get("accession", "")
        is_refseq = 1 if accession.startswith("GCF_") else 0
        annotated = 1 if row.get("annotated", "") == "1" else 0

        refseq_category = row.get("refseq_category", "").lower()
        reference_or_rep = 1 if ("reference genome" in refseq_category or "representative genome" in refseq_category) else 0
        level = level_rank.get(row.get("assembly_level", "").lower(), 0)
        pref = (assembly_preference or "RefSeq").lower()
        is_genbank = 1 if accession.startswith("GCA_") else 0
        pref_match = 1 if ((pref == "refseq" and is_refseq) or (pref == "genbank" and is_genbank) or pref in {"any", "none", ""}) else 0

        score = annotated * 1000 + is_refseq * 500 + pref_match * 200 + reference_or_rep * 100 + level
        notes = []
        if annotated:
            notes.append("annotated")
        if is_refseq:
            notes.append("refseq")
        if pref_match:
            notes.append("matches_species_registry_assembly_preference")
        if reference_or_rep:
            notes.append("reference_or_representative")
        notes.append(f"level_score={level}")

        scored = dict(row)
        scored["assembly_score"] = str(score)
        scored["assembly_decision_notes"] = ";".join(notes)
        scored_rows.append(scored)

    best = sorted(scored_rows, key=lambda r: int(r["assembly_score"]), reverse=True)[0]
    return best, scored_rows


def parse_gff3_attributes(attr_text: str) -> Dict[str, str]:
    """Parse and URL-decode a GFF3 attribute column."""
    out: Dict[str, str] = {}
    for item in attr_text.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[unquote(k.strip())] = unquote(v.strip())
    return out


def split_gff3_multi_value(value: str) -> List[str]:
    return [v.strip() for v in to_str(value).replace("|", ",").split(",") if v.strip()]


def symbol_matches_gene(attr: Dict[str, str], gene_symbol: str) -> bool:
    target = gene_symbol.strip().lower()
    fields = ["gene", "Name", "gene_synonym", "gene_name", "product", "description"]
    for field in fields:
        raw = attr.get(field, "")
        for token in split_gff3_multi_value(raw):
            if token.strip().lower() == target:
                return True
    return False


def normalize_biotype(value: str, has_cds: bool = False) -> str:
    v = to_str(value).strip().lower().replace("-", "_").replace(" ", "_")
    if v in {"protein_coding", "mrna", "gene"} or has_cds:
        return "coding_gene"
    if "pseudo" in v:
        return "pseudogene"
    if v in {"lncrna", "ncrna", "mirna", "snorna", "snrna", "rrna", "trna"}:
        return "noncoding_gene"
    return v or "unknown"


def classify_ncbi_model_status(tx_count: int, translated_tx_count: int, exon_count: int, warnings: List[str]) -> Tuple[str, str]:
    if tx_count >= 1 and translated_tx_count >= 1 and exon_count >= 1 and not warnings:
        return "model_found_high_confidence", "high"
    if tx_count >= 1 and exon_count >= 1:
        return "model_found_medium_confidence", "medium"
    if tx_count >= 1:
        return "model_found_partial", "low"
    return "no_model_found", "low"


#: GFF3 feature types read as a transcript. Shared by the candidate tally and the
#: parser so that a candidate's reported transcript count is the number the parser
#: would actually read from that locus.
TRANSCRIPT_FEATURE_TYPES = frozenset({"mRNA", "transcript"})


def collect_gene_candidates(gff_path: Path, gene_symbol: str,
                            expected_gene_ids: Sequence[str] = (),
                            product_name: str = "",
                            ) -> List[gid.GeneCandidate]:
    """Every locus in the annotation that could be the requested gene.

    A first pass over the file, before any transcript is read, because which locus to
    read depends on a decision that needs all the candidates in hand. The old parser
    decided on the spot from a token match and could not see a LOC-labelled locus at
    all, which is how an annotated gene stays invisible in an annotated genome.
    """
    seen: Dict[str, gid.GeneCandidate] = {}
    # Transcripts and distinct proteins per locus, tallied in the same pass. These are
    # part of what makes a candidate assessable — a locus annotated with no transcript
    # cannot yield a model, and the difference has to be visible in the candidate
    # inventory rather than inferred later from an empty transcript table.
    transcripts_of: Dict[str, int] = {}
    proteins_of: Dict[str, set] = {}
    tx_parent: Dict[str, str] = {}

    with gff_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attrs = parts
            if ftype in ("gene", "pseudogene"):
                attr = parse_gff3_attributes(attrs)
                cand = gid.candidate_from_attributes(
                    attr, seqid=seqid, strand=strand, start=start, end=end)
                route, _evidence = gid.classify_route(
                    cand, gene_symbol, expected_gene_ids, product_name)
                if not route:
                    continue
                key = cand.source_gene_id or f"{seqid}:{start}-{end}"
                seen.setdefault(key, cand)
            elif ftype in TRANSCRIPT_FEATURE_TYPES:
                attr = parse_gff3_attributes(attrs)
                parent = (attr.get("Parent") or "").split(",")[0]
                if parent:
                    transcripts_of[parent] = transcripts_of.get(parent, 0) + 1
                    if attr.get("ID"):
                        tx_parent[attr["ID"]] = parent
            elif ftype == "CDS":
                attr = parse_gff3_attributes(attrs)
                protein = attr.get("protein_id") or ""
                parent = tx_parent.get((attr.get("Parent") or "").split(",")[0], "")
                if protein and parent:
                    proteins_of.setdefault(parent, set()).add(protein)

    for key, cand in seen.items():
        cand.transcript_count = transcripts_of.get(key, 0)
        cand.protein_count = len(proteins_of.get(key, ()))
    return list(seen.values())


def read_archive_proteins(gff_path: Path, candidates: Sequence[gid.GeneCandidate],
                          ) -> Dict[str, str]:
    """One representative protein per candidate locus, from the package's protein FASTA.

    Only needed for candidates whose annotation evidence is too weak to name the
    paralog. NCBI protein FASTA headers carry the gene symbol in brackets, which is
    what links a sequence back to its locus without re-deriving the translation.
    """
    package = gff_path.parent
    faa = next((p for p in sorted(package.rglob("protein.faa"))), None)
    if faa is None:
        return {}
    by_symbol: Dict[str, str] = {}
    header = ""
    chunks: List[str] = []

    def flush() -> None:
        if not header or not chunks:
            return
        match = re.search(r"\[gene=([^\]]+)\]", header)
        symbol = (match.group(1) if match else "").strip()
        if symbol and symbol not in by_symbol:
            by_symbol[symbol] = "".join(chunks)

    try:
        with faa.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    flush()
                    header, chunks = line[1:].strip(), []
                else:
                    chunks.append(line.strip())
        flush()
    except OSError:
        return {}

    out: Dict[str, str] = {}
    for cand in candidates:
        sequence = by_symbol.get(cand.source_symbol, "")
        if sequence:
            out[cand.source_gene_id] = sequence
    return out


def parse_ncbi_gff3_for_gene(
    gff_path: Path,
    species: SpeciesRecord,
    gene_symbol: str,
    expected_gene_ids: Sequence[str] = (),
    paralog_panel: Optional[Dict[str, str]] = None,
    identification_out: Optional[List[gid.Identification]] = None,
    product_name: str = "",
) -> Tuple[Optional[Tuple[GeneModel, List[TranscriptModel], List[ExonModel], List[CDSFeature]]], List[str]]:
    """
    Conservative parser for NCBI GFF3.

    Strategy:
    - find exactly one gene locus matching FGFR2
    - collect mRNA/transcript features when available
    - collect exon features if present; otherwise reconstruct exon-like intervals from CDS
    - estimate protein length from CDS span sum / 3 if translation ids are unavailable
    """
    warnings: List[str] = []

    # Decide which locus is the requested gene before reading a single transcript.
    # Deciding while scanning is what let a token match in a paralog's synonym list
    # win, and what made a LOC-labelled locus unreachable.
    candidates = collect_gene_candidates(gff_path, gene_symbol, expected_gene_ids,
                                         product_name)
    proteins = (read_archive_proteins(gff_path, candidates)
                if any(gid.classify_route(c, gene_symbol, expected_gene_ids,
                                          product_name)[0]
                       in gid.ROUTES_NEEDING_DISCRIMINATION for c in candidates)
                else {})
    identification = gid.identify(candidates, gene_symbol,
                                  expected_gene_ids=expected_gene_ids,
                                  proteins=proteins, panel=paralog_panel or {},
                                  product_name=product_name)
    if identification_out is not None:
        identification_out.append(identification)
    if identification.accepted is None:
        warnings.append(f"{identification.status}: {identification.detail}")
        return None, warnings
    accepted_gene_id = identification.accepted.source_gene_id

    genes = []
    transcripts = []
    transcript_to_exons: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    transcript_to_cds: Dict[str, List[Dict[str, str]]] = defaultdict(list)

    current_gene_ids = set()

    with gff_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _source, feature_type, start, end, _score, strand, phase, attrs = parts
            attr = parse_gff3_attributes(attrs)

            if feature_type in ("gene", "pseudogene") and attr.get("ID", "") == accepted_gene_id:
                gene_id = attr.get("ID", "")
                genes.append(
                    {
                        "seqid": seqid,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "gene_id": gene_id,
                        # The symbol the annotation actually uses, which for a
                        # LOC-labelled locus is not the requested symbol. Recording the
                        # requested one here would erase how the gene was identified.
                        "symbol": identification.accepted.source_symbol or gene_symbol,
                        "biotype": attr.get("gene_biotype", attr.get("gbkey", "")),
                    }
                )
                if gene_id:
                    current_gene_ids.add(gene_id)

            elif feature_type in TRANSCRIPT_FEATURE_TYPES:
                parents = split_gff3_multi_value(attr.get("Parent", ""))
                if any(parent in current_gene_ids for parent in parents):
                    tx_id = attr.get("ID", "")
                    transcripts.append(
                        {
                            "tx_id": tx_id,
                            "name": attr.get("Name", ""),
                            "biotype": attr.get("transcript_biotype", feature_type),
                            "seqid": seqid,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "protein_id": attr.get("protein_id", ""),
                        }
                    )

            elif feature_type == "exon":
                parents = split_gff3_multi_value(attr.get("Parent", ""))
                for parent in parents:
                    transcript_to_exons[parent].append(
                        {
                            "exon_id": attr.get("ID", ""),
                            "cds_id_source": attr.get("ID", ""),
                            "parent_feature_id": parent,
                            "seqid": seqid,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "phase": "",
                            "end_phase": "",
                        }
                    )

            elif feature_type == "CDS":
                parents = split_gff3_multi_value(attr.get("Parent", ""))
                for parent in parents:
                    transcript_to_cds[parent].append(
                        {
                            "exon_id": attr.get("ID", ""),
                            "cds_id_source": attr.get("ID", ""),
                            "parent_feature_id": parent,
                            "seqid": seqid,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "phase": phase if phase != "." else "",
                            "end_phase": "",
                            "protein_id": attr.get("protein_id", ""),
                        }
                    )

    if len(genes) != 1:
        warnings.append(
            f"the accepted {gene_symbol} locus {accepted_gene_id!r} was identified but "
            f"{len(genes)} gene features carry that ID on re-scan; the annotation file "
            "may be inconsistent")
        return None, warnings

    gene = genes[0]
    gene_internal = f"{species.ensembl_species}|NCBI|{gene['gene_id'] or gene['symbol']}"

    # If no explicit transcripts but CDS parents exist, synthesize transcript entries.
    if not transcripts and transcript_to_cds:
        warnings.append("no explicit mRNA/transcript features; synthesizing transcripts from CDS parents")
        for tx_id, cds_parts in transcript_to_cds.items():
            starts = [to_int(c["start"]) for c in cds_parts if to_int(c["start"]) is not None]
            ends = [to_int(c["end"]) for c in cds_parts if to_int(c["end"]) is not None]
            if not starts or not ends:
                continue
            transcripts.append(
                {
                    "tx_id": tx_id,
                    "name": tx_id,
                    "biotype": "mRNA",
                    "seqid": cds_parts[0]["seqid"],
                    "start": str(min(starts)),
                    "end": str(max(ends)),
                    "strand": cds_parts[0]["strand"],
                    "protein_id": next((c.get("protein_id", "") for c in cds_parts if c.get("protein_id", "")), ""),
                }
            )

    tx_models: List[TranscriptModel] = []
    ex_models: List[ExonModel] = []
    cds_feature_models: List[CDSFeature] = []
    translated_tx_count = 0

    for tx in transcripts:
        tx_id = tx["tx_id"]
        tx_internal = f"{gene_internal}|tx|{tx_id}"

        exs = transcript_to_exons.get(tx_id, [])
        cds_parts = transcript_to_cds.get(tx_id, [])

        use_parts = exs if exs else cds_parts
        if not use_parts:
            # maybe Parent references differ; skip but warn
            warnings.append(f"transcript {tx_id} has no exon or CDS child features")
            continue

        strand = tx["strand"]
        sorted_parts = sorted(use_parts, key=lambda r: int(r["start"]),
                              reverse=is_reverse(strand))

        # Estimate protein length from CDS spans.
        protein_len = ""
        translation_id = tx.get("protein_id", "")
        if cds_parts:
            cds_bp = sum((to_int(c["end"]) or 0) - (to_int(c["start"]) or 0) + 1 for c in cds_parts)
            if cds_bp > 0:
                # RefSeq/Gnomon CDS features include the terminator, so a complete CDS holds
                # one codon more than the protein has residues. Counting it made every
                # estimated length one too long, which in turn let the projection place a
                # residue where the protein only has a stop.
                codons = cds_bp // 3
                complete = cds_bp % 3 == 0
                protein_len = str(max(0, codons - 1 if complete else codons))
        protein_length_source = ("NCBI_GFF3_CDS_span_estimate_stop_excluded" if protein_len
                                 else "missing")
        if translation_id or protein_len:
            translated_tx_count += 1

        completeness_flags = []
        if not translation_id:
            completeness_flags.append("no_translation_id_from_gff3")
        if not cds_parts:
            completeness_flags.append("no_cds_features")
        tx_conf = "high" if cds_parts and use_parts else "medium"

        tx_models.append(
            TranscriptModel(
                source_db="NCBI",
                species_canonical=species.ensembl_species,
                gene_id_internal=gene_internal,
                transcript_id_source=tx_id,
                transcript_name=tx["name"],
                transcript_biotype=tx["biotype"],
                translation_id_source=translation_id,
                protein_length_aa=protein_len,
                protein_length_source=protein_length_source,
                is_canonical_source="",
                support_level="",
                completeness_flags=";".join(completeness_flags),
                transcript_model_confidence=tx_conf,
                chrom=tx["seqid"],
                start=tx["start"],
                end=tx["end"],
                strand=tx["strand"],
                internal_transcript_id=tx_internal,
            )
        )

        for rank, ex in enumerate(sorted_parts, start=1):
            exon_id = ex.get("exon_id", "") or f"{tx_id}_part{rank}"
            ex_models.append(
                ExonModel(
                    source_db="NCBI",
                    species_canonical=species.ensembl_species,
                    transcript_id_internal=tx_internal,
                    exon_id_source=exon_id,
                    exon_rank=str(rank),
                    chrom=ex["seqid"],
                    start=ex["start"],
                    end=ex["end"],
                    strand=ex["strand"],
                    phase=ex.get("phase", ""),
                    end_phase=ex.get("end_phase", ""),
                    internal_exon_id=f"{tx_internal}|exon|{rank}|{exon_id}",
                )
            )
        cds_feature_models.extend(build_cds_features_from_parts(
            source_db="NCBI", species_input=species.input_name, species_canonical=species.ensembl_species,
            transcript_id_internal=tx_internal, transcript_id_source=tx_id, translation_id_source=translation_id,
            parts=cds_parts, coordinate_source="NCBI_GFF3_CDS_features",
        ))

    model_status, model_confidence = classify_ncbi_model_status(
        tx_count=len(tx_models),
        translated_tx_count=translated_tx_count,
        exon_count=len(ex_models),
        warnings=warnings,
    )

    if not tx_models:
        warnings.append("no parseable FGFR2 transcripts found in NCBI GFF3")
        return None, warnings

    gene_model = GeneModel(
        source_db="NCBI",
        species_input=species.input_name,
        species_canonical=species.ensembl_species,
        taxid=species.taxid,
        assembly_accession="",
        assembly_name="",
        gene_symbol_requested=gene_symbol,
        gene_symbol_found=gene["symbol"],
        gene_id_source=gene["gene_id"],
        gene_biotype=normalize_biotype(gene["biotype"], has_cds=bool(ex_models)),
        chrom=gene["seqid"],
        start=gene["start"],
        end=gene["end"],
        strand=gene["strand"],
        model_status=model_status,
        model_confidence=model_confidence,
        lookup_method="ncbi_datasets_gff3",
        internal_gene_id=gene_internal,
    )

    return (gene_model, tx_models, ex_models, cds_feature_models), warnings


def choose_best_gff3_file(gff_candidates: List[Path]) -> Optional[Path]:
    if not gff_candidates:
        return None

    def score(path: Path) -> Tuple[int, int, str]:
        name = path.name.lower()
        s = 0
        if name in {"genomic.gff", "genomic.gff3"}:
            s += 100
        if name.endswith(".gff3"):
            s += 20
        if "ncbi_dataset" in str(path).lower():
            s += 10
        try:
            size = path.stat().st_size
        except OSError:
            size = 0
        return (s, size, str(path))

    return sorted(gff_candidates, key=score, reverse=True)[0]


def fetch_ncbi_model_from_datasets(
    species: SpeciesRecord,
    gene_symbol: str,
    workdir: Path,
    datasets_bin: Optional[str],
    datasets_include: str,
    datasets_summary_timeout: int,
    datasets_download_timeout: int,
    sleep_s: float,
    expected_gene_ids: Sequence[str] = (),
    paralog_panel: Optional[Dict[str, str]] = None,
    identification_out: Optional[List[gid.Identification]] = None,
    product_name: str = "",
) -> Tuple[Optional[Tuple[GeneModel, List[TranscriptModel], List[ExonModel], List[CDSFeature]]], List[str], Dict[str, str], List[Dict[str, str]]]:
    warnings: List[str] = []
    provenance = {
        "datasets_used": "0",
        "datasets_summary_status": "",
        "datasets_download_status": "",
        "datasets_package": "",
        "taxid": "",
        "assembly_accession": "",
        "assembly_name": "",
        "assembly_level": "",
        "assembly_status": "",
        "assembly_score": "",
        "assembly_decision_notes": "",
    }
    assembly_selection_rows: List[Dict[str, str]] = []

    if not datasets_bin:
        warnings.append("NCBI Datasets CLI not available")
        provenance["datasets_summary_status"] = "datasets_cli_not_found"
        return None, warnings, provenance, assembly_selection_rows

    # The numeric taxid, or the *accepted* scientific name — never the submitted slug.
    # Passing the slug is what produced "The taxonomy name 'equus_quagga' is not
    # recognized" and ended the Equus quagga run before any assembly was listed. The
    # registry now resolves the name, so an empty query term here means resolution
    # failed and there is nothing worth asking the service.
    taxon = (species.taxid or "").strip() or (species.ncbi_species or "").strip()
    if not taxon:
        warnings.append(
            "species identity was not resolved against NCBI Taxonomy, so no assembly "
            "query was issued; a query under the unverified name would be rejected")
        provenance["datasets_summary_status"] = asel.NO_QUERY_TERM
        assembly_selection_rows.append(asel.failure_row(
            asel.NO_QUERY_TERM,
            "no verified taxid or accepted scientific name for this species",
            species.input_name, species.ensembl_species, species.taxid, ""))
        return None, warnings, provenance, assembly_selection_rows

    summary_cmd = [
        datasets_bin,
        "summary",
        "genome",
        "taxon",
        taxon,
        "--as-json-lines",
    ]
    try:
        code, output = run_command(summary_cmd, timeout=datasets_summary_timeout)
    except subprocess.TimeoutExpired:
        warnings.append("datasets summary timed out")
        provenance["datasets_summary_status"] = "timeout"
        return None, warnings, provenance, assembly_selection_rows

    provenance["datasets_used"] = "1"
    provenance["datasets_summary_status"] = f"exit_{code}"
    if code != 0:
        # Two very different problems arrive on this branch: a query term the service
        # will never accept, and a service that is temporarily unhappy. The first needs
        # the name corrected, the second needs a retry, and the original run reported
        # both as an empty table.
        rejected = "not recognized" in output.lower() or "not found" in output.lower()
        status = asel.TAXON_REJECTED if rejected else asel.SERVICE_FAILED
        provenance["datasets_summary_status"] = status
        warnings.append(f"datasets summary failed ({status}): {output[:500]}")
        assembly_selection_rows.append(asel.failure_row(
            status, output.strip()[:300], species.input_name,
            species.ensembl_species, species.taxid, taxon))
        return None, warnings, provenance, assembly_selection_rows

    reports = []
    for line in output.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            reports.append(json.loads(line))
        except Exception:
            pass

    if not reports:
        provenance["datasets_summary_status"] = asel.SERVICE_FAILED
        warnings.append("datasets summary returned no parseable JSON")
        assembly_selection_rows.append(asel.failure_row(
            asel.SERVICE_FAILED, "the summary response could not be parsed as JSON",
            species.input_name, species.ensembl_species, species.taxid, taxon))
        return None, warnings, provenance, assembly_selection_rows

    merged: List[Dict[str, Any]] = []
    for obj in reports:
        if isinstance(obj, dict):
            if "reports" in obj and isinstance(obj["reports"], list):
                merged.extend(obj["reports"])
            else:
                merged.append(obj)

    selection = asel.select(
        asel.parse_summary(merged),
        preference=species.assembly_preference,
        requested_taxid=species.taxid,
        requested_name=species.ncbi_species,
        query_term=taxon,
    )
    provenance["assembly_selection_status"] = selection.status
    provenance["assembly_selection_detail"] = selection.detail
    provenance["n_assembly_candidates"] = str(len(selection.candidates))
    provenance["n_annotated_candidates"] = str(selection.n_annotated)

    resolved_taxid = ((species.taxid or "").strip()
                      or ((selection.selected or {}).get("tax_id") or "").strip())
    provenance["taxid"] = resolved_taxid

    assembly_selection_rows.extend(asel.selection_rows(
        selection, species.input_name, species.ensembl_species, resolved_taxid))

    if selection.selected is None:
        # An honest absence. The candidate rows above already say which assemblies exist
        # and why none of them can yield a gene model, which is what the empty table
        # could not say.
        warnings.append(f"{selection.status}: {selection.detail}")
        return None, warnings, provenance, assembly_selection_rows

    best = selection.selected

    dir_key = resolved_taxid or re.sub(r"[^a-z0-9]+", "_", taxon.lower()).strip("_")
    species_dir = workdir / f"ncbi_{dir_key}"
    species_dir.mkdir(parents=True, exist_ok=True)

    accession = best["accession"]
    provenance["assembly_accession"] = accession
    provenance["assembly_name"] = best.get("assembly_name", "")
    provenance["assembly_level"] = best.get("assembly_level", "")
    provenance["assembly_status"] = best.get("refseq_category", "") or best.get("assembly_status", "")
    provenance["assembly_score"] = best.get("assembly_score", "")
    provenance["assembly_decision_notes"] = best.get("assembly_decision_notes", "")

    zip_dir = species_dir / accession
    zip_dir.mkdir(parents=True, exist_ok=True)
    zip_path = zip_dir / "ncbi_dataset.zip"

    download_cmd = [
        datasets_bin,
        "download",
        "genome",
        "accession",
        accession,
        "--include",
        datasets_include,
        "--filename",
        str(zip_path),
    ]
    try:
        code, output = run_command(download_cmd, timeout=datasets_download_timeout)
    except subprocess.TimeoutExpired:
        warnings.append("datasets accession download timed out")
        provenance["datasets_download_status"] = "timeout"
        provenance["assembly_selection_status"] = asel.DOWNLOAD_FAILED
        return None, warnings, provenance, assembly_selection_rows

    provenance["datasets_download_status"] = f"exit_{code}"
    provenance["datasets_package"] = str(zip_path)

    if code != 0 or not zip_path.exists():
        warnings.append(f"datasets accession download failed: {output[:500]}")
        provenance["assembly_selection_status"] = asel.DOWNLOAD_FAILED
        provenance["assembly_selection_detail"] = output.strip()[:300]
        return None, warnings, provenance, assembly_selection_rows

    extract_dir = zip_dir / "unzipped"
    extract_dir.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(extract_dir)
    except Exception as exc:
        # An archive that arrived and cannot be opened is a processing fault. Reporting
        # it as a parser failure rather than as missing data is what keeps a user from
        # being told their species has no annotation when in fact it has.
        warnings.append(f"datasets zip extraction failed: {exc}")
        provenance["assembly_selection_status"] = asel.PARSE_FAILED
        provenance["assembly_selection_detail"] = f"zip extraction failed: {exc}"
        return None, warnings, provenance, assembly_selection_rows

    gff_candidates = sorted(set(list(extract_dir.rglob("*.gff")) + list(extract_dir.rglob("*.gff3"))))
    gff_path = choose_best_gff3_file(gff_candidates)
    if gff_path is None:
        warnings.append("no GFF3 found in downloaded NCBI assembly package")
        provenance["assembly_selection_status"] = asel.PARSE_FAILED
        provenance["assembly_selection_detail"] = (
            "the downloaded package contained no GFF3 annotation file; the Datasets "
            "package layout may have changed")
        return None, warnings, provenance, assembly_selection_rows
    provenance["gff3_path"] = str(gff_path)

    try:
        parsed, parse_warnings = parse_ncbi_gff3_for_gene(
            gff_path, species, gene_symbol,
            expected_gene_ids=expected_gene_ids,
            paralog_panel=paralog_panel,
            identification_out=identification_out,
            product_name=product_name)
    except Exception as exc:
        # A parser exception must never reach the caller as "no model". That conversion
        # is how a genuine defect ends up presented as an empty table, which is what
        # made the original failure unreadable.
        warnings.append(f"GFF3 parsing raised {type(exc).__name__}: {exc}")
        provenance["assembly_selection_status"] = asel.PARSE_FAILED
        provenance["assembly_selection_detail"] = f"{type(exc).__name__}: {exc}"
        provenance["parser_traceback_in_log"] = "1"
        return None, warnings, provenance, assembly_selection_rows
    warnings.extend(parse_warnings)
    if parsed is None:
        warnings.append(
            f"the annotation for {accession} was read but contained no locus "
            f"identifiable as {gene_symbol}")
        provenance["assembly_selection_status"] = "gene_not_in_annotation"
        return None, warnings, provenance, assembly_selection_rows

    # Stamp assembly metadata onto gene model.
    gene_model, tx_models, ex_models, cds_features = parsed
    gene_model.assembly_accession = provenance["assembly_accession"]
    gene_model.assembly_name = provenance["assembly_name"]

    safe_sleep(sleep_s)
    return (gene_model, tx_models, ex_models, cds_features), warnings, provenance, assembly_selection_rows


# ---------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------
def seqid_normalization_status(chrom_a: str, chrom_b: str) -> Tuple[str, str]:
    a = to_str(chrom_a).strip()
    b = to_str(chrom_b).strip()
    if not a or not b:
        return "unknown", "missing_seqid"
    if a == b:
        return "same", "exact_match"

    # Common case: Ensembl uses bare chromosome names; NCBI uses accessions.
    if a.startswith("NC_") or b.startswith("NC_"):
        return "different_label", "seqid_naming_difference_only"
    return "different_label", "different_seqid_labels"


def count_translated_transcripts(txs: List[TranscriptModel]) -> int:
    return sum(1 for tx in txs if tx.translation_id_source or to_int(tx.protein_length_aa) is not None)


def max_protein_length(txs: List[TranscriptModel]) -> Optional[int]:
    values = [to_int(tx.protein_length_aa) for tx in txs]
    values = [v for v in values if v is not None]
    return max(values) if values else None


def plausible_long_count(txs: List[TranscriptModel], cutoff: int = 750) -> int:
    return sum(1 for tx in txs if (to_int(tx.protein_length_aa) or -1) >= cutoff)


def compute_gene_overlap_fraction_if_comparable(a: Optional[GeneModel], b: Optional[GeneModel]) -> Tuple[str, str]:
    if not a or not b:
        return "NA", "missing_model"
    if a.chrom != b.chrom:
        return "NA", "different_coordinate_system_or_seqid_labels"
    a_start, a_end = to_int(a.start), to_int(a.end)
    b_start, b_end = to_int(b.start), to_int(b.end)
    if None in {a_start, a_end, b_start, b_end}:
        return "NA", "missing_coordinates"
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    if overlap_start > overlap_end:
        return safe_float_str(0.0), "same_seqid_no_overlap"
    ov = overlap_end - overlap_start + 1
    denom = min(a_end - a_start + 1, b_end - b_start + 1)
    if denom <= 0:
        return "NA", "invalid_coordinate_span"
    return safe_float_str(ov / denom), "same_coordinate_system"


def compare_models(
    ensembl_gene: Optional[GeneModel],
    ensembl_txs: List[TranscriptModel],
    ensembl_exons: List[ExonModel],
    ncbi_gene: Optional[GeneModel],
    ncbi_txs: List[TranscriptModel],
    ncbi_exons: List[ExonModel],
    preferred_source: str = "ensembl_first",
) -> Tuple[str, str, str, Dict[str, str], List[str]]:
    """
    Returns:
      selected_source, selected_confidence, detailed_reason, comparison_row, warnings
    """
    warnings: List[str] = []

    row = {
        "ensembl_available": "1" if ensembl_gene else "0",
        "ncbi_available": "1" if ncbi_gene else "0",
        "same_gene_symbol": "0",
        "same_biotype": "0",
        "same_chrom_label": "NA",
        "seqid_relation": "NA",
        "gene_overlap_fraction": "NA",
        "gene_overlap_reason": "NA",
        "translation_support_ensembl": str(count_translated_transcripts(ensembl_txs)),
        "translation_support_ncbi": str(count_translated_transcripts(ncbi_txs)),
        "max_protein_len_ensembl": to_str(max_protein_length(ensembl_txs)),
        "max_protein_len_ncbi": to_str(max_protein_length(ncbi_txs)),
        "long_transcript_count_ensembl": str(plausible_long_count(ensembl_txs)),
        "long_transcript_count_ncbi": str(plausible_long_count(ncbi_txs)),
        "transcript_count_ensembl": str(len(ensembl_txs)),
        "transcript_count_ncbi": str(len(ncbi_txs)),
        "exon_count_ensembl": str(len(ensembl_exons)),
        "exon_count_ncbi": str(len(ncbi_exons)),
        "transcript_count_diff": "NA",
        "exon_count_diff": "NA",
        "conflict_level": "NA",
        "ensembl_source_score": "NA",
        "ncbi_source_score": "NA",
        "source_score_difference": "NA",
        "preferred_source_applied": preferred_source,
        "comparison_decision_reason": "",
    }

    if ensembl_gene and ncbi_gene:
        row["same_gene_symbol"] = "1" if ensembl_gene.gene_symbol_found == ncbi_gene.gene_symbol_found else "0"
        row["same_biotype"] = "1" if normalize_biotype(ensembl_gene.gene_biotype) == normalize_biotype(ncbi_gene.gene_biotype, has_cds=bool(ncbi_exons)) else "0"
        seqid_state, seqid_reason = seqid_normalization_status(ensembl_gene.chrom, ncbi_gene.chrom)
        row["same_chrom_label"] = "1" if seqid_state == "same" else "0"
        row["seqid_relation"] = seqid_reason
        row["transcript_count_diff"] = str(abs(len(ensembl_txs) - len(ncbi_txs)))
        row["exon_count_diff"] = str(abs(len(ensembl_exons) - len(ncbi_exons)))
        overlap_frac, overlap_reason = compute_gene_overlap_fraction_if_comparable(ensembl_gene, ncbi_gene)
        row["gene_overlap_fraction"] = overlap_frac
        row["gene_overlap_reason"] = overlap_reason

    # Determine conflict level.
    if ensembl_gene and not ncbi_gene:
        row["conflict_level"] = "no_conflict"
        row["comparison_decision_reason"] = "only_ensembl_structured_model_available"
        return "Ensembl", ensembl_gene.model_confidence, row["comparison_decision_reason"], row, warnings
    if ncbi_gene and not ensembl_gene:
        row["conflict_level"] = "no_conflict"
        row["comparison_decision_reason"] = "only_ncbi_structured_model_available"
        return "NCBI", ncbi_gene.model_confidence, row["comparison_decision_reason"], row, warnings
    if not ensembl_gene and not ncbi_gene:
        row["conflict_level"] = "major_conflict"
        row["comparison_decision_reason"] = "no_structured_model_found"
        warnings.append("no structured model found in either source")
        return "NONE", "none", row["comparison_decision_reason"], row, warnings

    # Both available.
    same_symbol = row["same_gene_symbol"] == "1"
    same_biotype = row["same_biotype"] == "1"
    tx_diff = abs(len(ensembl_txs) - len(ncbi_txs))
    ex_diff = abs(len(ensembl_exons) - len(ncbi_exons))
    seqid_relation = row["seqid_relation"]
    max_len_e = max_protein_length(ensembl_txs) or -1
    max_len_n = max_protein_length(ncbi_txs) or -1
    len_diff = abs(max_len_e - max_len_n) if (max_len_e > 0 and max_len_n > 0) else 999

    if same_symbol and same_biotype and tx_diff <= 5 and ex_diff <= 80 and len_diff <= 100:
        if seqid_relation == "seqid_naming_difference_only":
            row["conflict_level"] = "minor_conflict"
            warnings.append("seqid label difference only between Ensembl and NCBI")
        else:
            row["conflict_level"] = "no_conflict"
    elif same_symbol and same_biotype:
        row["conflict_level"] = "minor_conflict"
    else:
        row["conflict_level"] = "moderate_conflict"

    # Source decision.
    ensembl_score = 0
    ncbi_score = 0

    conf_rank = {"high": 3, "medium": 2, "low": 1, "none": 0}
    ensembl_score += conf_rank.get(ensembl_gene.model_confidence, 0) * 100
    ncbi_score += conf_rank.get(ncbi_gene.model_confidence, 0) * 100

    ensembl_score += count_translated_transcripts(ensembl_txs) * 10
    ncbi_score += count_translated_transcripts(ncbi_txs) * 10

    ensembl_score += plausible_long_count(ensembl_txs) * 15
    ncbi_score += plausible_long_count(ncbi_txs) * 15

    ensembl_score += (max_protein_length(ensembl_txs) or 0) // 100
    ncbi_score += (max_protein_length(ncbi_txs) or 0) // 100

    preferred = (preferred_source or "ensembl_first").lower()
    if preferred in {"ensembl_first", "strict_ensembl"}:
        ensembl_score += 10 if preferred == "ensembl_first" else 50
    elif preferred in {"ncbi_first", "strict_ncbi"}:
        ncbi_score += 10 if preferred == "ncbi_first" else 50

    row["ensembl_source_score"] = str(ensembl_score)
    row["ncbi_source_score"] = str(ncbi_score)
    row["source_score_difference"] = str(abs(ensembl_score - ncbi_score))

    if preferred == "strict_ensembl" and ensembl_gene and ensembl_gene.model_confidence in {"high", "medium"}:
        reason = "strict_ensembl_preference_applied"
        row["comparison_decision_reason"] = reason
        return "Ensembl", ensembl_gene.model_confidence, reason, row, warnings
    if preferred == "strict_ncbi" and ncbi_gene and ncbi_gene.model_confidence in {"high", "medium"}:
        reason = "strict_ncbi_preference_applied"
        row["comparison_decision_reason"] = reason
        return "NCBI", ncbi_gene.model_confidence, reason, row, warnings

    if ensembl_score >= ncbi_score:
        reason = (
            "ensembl_preferred_same_locus_minor_label_differences"
            if seqid_relation == "seqid_naming_difference_only" and row["conflict_level"] in {"no_conflict", "minor_conflict"}
            else "ensembl_preferred_higher_confidence_or_better_translation_support"
        )
        row["comparison_decision_reason"] = reason
        return "Ensembl", ensembl_gene.model_confidence, reason, row, warnings
    else:
        reason = "ncbi_preferred_higher_confidence_or_better_translation_support"
        row["comparison_decision_reason"] = reason
        return "NCBI", ncbi_gene.model_confidence, reason, row, warnings


# ---------------------------------------------------------------------
# Internal consistency checks
# ---------------------------------------------------------------------
def run_internal_consistency_checks(
    genes_rows: List[Dict[str, str]],
    transcripts_rows: List[Dict[str, str]],
    exons_rows: List[Dict[str, str]],
) -> List[Dict[str, str]]:
    checks: List[Dict[str, str]] = []

    gene_ids = {r["internal_gene_id"] for r in genes_rows}
    tx_ids = {r["internal_transcript_id"] for r in transcripts_rows}

    orphan_txs = [r for r in transcripts_rows if r.get("gene_id_internal", "") not in gene_ids]
    checks.append(
        {
            "check_name": "transcripts_reference_existing_genes",
            "status": "PASS" if not orphan_txs else "FAIL",
            "affected_species": ",".join(sorted({r.get("species_canonical", "") for r in orphan_txs})),
            "details": f"orphan_transcripts={len(orphan_txs)}",
        }
    )

    orphan_exons = [r for r in exons_rows if r.get("transcript_id_internal", "") not in tx_ids]
    checks.append(
        {
            "check_name": "exons_reference_existing_transcripts",
            "status": "PASS" if not orphan_exons else "FAIL",
            "affected_species": ",".join(sorted({r.get("species_canonical", "") for r in orphan_exons})),
            "details": f"orphan_exons={len(orphan_exons)}",
        }
    )

    invalid_coords = [r for r in exons_rows if (to_int(r.get("start", "")) is None or to_int(r.get("end", "")) is None or (to_int(r.get("start", "")) or 0) > (to_int(r.get("end", "")) or 0))]
    checks.append(
        {
            "check_name": "exon_coordinates_valid",
            "status": "PASS" if not invalid_coords else "FAIL",
            "affected_species": ",".join(sorted({r.get("species_canonical", "") for r in invalid_coords})),
            "details": f"invalid_exons={len(invalid_coords)}",
        }
    )

    by_tx = defaultdict(list)
    for row in exons_rows:
        by_tx[row["transcript_id_internal"]].append(row)

    gap_species = set()
    order_species = set()
    for _tx_id, rows in by_tx.items():
        species = rows[0].get("species_canonical", "")
        ranks = sorted([to_int(r["exon_rank"]) for r in rows if to_int(r["exon_rank"]) is not None])
        if ranks and ranks != list(range(1, len(ranks) + 1)):
            gap_species.add(species)

        strand = rows[0].get("strand", "")
        sorted_by_rank = sorted(rows, key=lambda r: int(r["exon_rank"]))
        coords = [(to_int(r["start"]) or 0, to_int(r["end"]) or 0) for r in sorted_by_rank]
        # Normalised, because this check only ever recognised the Ensembl spellings "1" and
        # "-1". A RefSeq transcript spelling the same strand "-" matched neither branch, so
        # the one check that would have caught a mis-ordered transcript silently skipped it.
        normalized = normalize_strand(strand)
        if normalized is not None:
            starts = [c[0] for c in coords]
            if starts != sorted(starts, reverse=normalized == MINUS):
                order_species.add(species)

    checks.append(
        {
            "check_name": "exon_rank_gaps_absent",
            "status": "PASS" if not gap_species else "FAIL",
            "affected_species": ",".join(sorted(gap_species)),
            "details": f"species_with_rank_gaps={len(gap_species)}",
        }
    )
    checks.append(
        {
            "check_name": "exon_coordinate_order_matches_rank",
            "status": "PASS" if not order_species else "FAIL",
            "affected_species": ",".join(sorted(order_species)),
            "details": f"species_with_order_issues={len(order_species)}",
        }
    )

    empty_models = [r for r in genes_rows if not any(tx.get("gene_id_internal") == r["internal_gene_id"] for tx in transcripts_rows)]
    checks.append(
        {
            "check_name": "selected_genes_have_transcripts",
            "status": "PASS" if not empty_models else "FAIL",
            "affected_species": ",".join(sorted({r.get("species_canonical", "") for r in empty_models})),
            "details": f"genes_without_transcripts={len(empty_models)}",
        }
    )

    return checks


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(description="Improved dual-source FGFR2 collector.")
    parser.add_argument("--species_registry", required=True, type=Path)
    parser.add_argument("--gene_symbol", default="FGFR2")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sleep_between_requests", type=float, default=0.25)
    parser.add_argument("--use_ncbi_datasets", action="store_true")
    parser.add_argument("--datasets_bin", default="")
    parser.add_argument("--datasets_include", default="gff3,protein")
    parser.add_argument("--gene_config", type=Path, default=None,
                        help="Gene configuration used only as an offline fallback for "
                             "the gene's spelled-out product name.")
    parser.add_argument("--paralog_reference_fasta", type=Path,
                        default=Path("references/human_FGFR1_2_3_4.fasta"),
                        help="Reference panel used to tell the requested gene apart "
                             "from its close paralogs when the annotation alone "
                             "cannot.")
    parser.add_argument("--datasets_summary_timeout", type=int, default=180)
    parser.add_argument("--datasets_download_timeout", type=int, default=900)
    args = parser.parse_args()

    species_registry = read_species_registry(args.species_registry)
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    datasets_bin = datasets_available(args.datasets_bin or None) if args.use_ncbi_datasets else None
    workdir = outdir / "_ncbi_datasets_cache"
    workdir.mkdir(parents=True, exist_ok=True)

    genes_rows: List[Dict[str, str]] = []
    transcripts_rows: List[Dict[str, str]] = []
    exons_rows: List[Dict[str, str]] = []
    cds_rows: List[Dict[str, str]] = []
    species_status_rows: List[Dict[str, str]] = []
    selection_rows: List[Dict[str, str]] = []
    source_conflict_rows: List[Dict[str, str]] = []
    source_comparison_rows: List[Dict[str, str]] = []
    annotation_warning_rows: List[Dict[str, str]] = []
    raw_source_map_rows: List[Dict[str, str]] = []
    ncbi_rescue_rows: List[Dict[str, str]] = []
    ncbi_assembly_selection_rows: List[Dict[str, str]] = []
    ncbi_assembly_selected_rows: List[Dict[str, str]] = []
    gene_candidate_rows: List[Dict[str, str]] = []
    outcomes: List[recovery.SpeciesOutcome] = []

    # Loaded once. Used only for candidates whose annotation evidence cannot say which
    # paralog was found; a locus with an official symbol never reaches it.
    paralog_panel = gid.read_panel(args.paralog_reference_fasta)
    if not paralog_panel:
        print(f"[WARN] paralog reference panel not readable at "
              f"{args.paralog_reference_fasta}; weak-route candidates cannot be "
              f"discriminated by sequence and will be reported as ambiguous")

    per_species_ncbi_provenance: Dict[str, Dict[str, str]] = {}

    run_meta = {
        "run_timestamp_epoch": int(time.time()),
        "script_name": "collect_fgfr2_models_dual_source_v3.py",
        "gene_symbol": args.gene_symbol,
        "species_count": len(species_registry),
        "use_ncbi_datasets": bool(args.use_ncbi_datasets),
        "datasets_bin": datasets_bin or "",
        "ncbi_datasets_strategy": "summary genome taxon -> deterministic best assembly -> download genome accession",
        "parameters": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
        },
        "per_species_ncbi_provenance": {},
    }

    for sp in species_registry:
        print(f"[RUN] {sp.input_name}")

        # Ensembl
        ensembl_gene_model: Optional[GeneModel] = None
        ensembl_txs: List[TranscriptModel] = []
        ensembl_exons: List[ExonModel] = []
        ensembl_cds: List[CDSFeature] = []
        ensembl_lookup_method = ""

        ensembl_model_raw, ensembl_lookup_method = fetch_ensembl_model(sp, args.gene_symbol, args.sleep_between_requests)
        if ensembl_model_raw is not None:
            ensembl_gene_model, ensembl_txs, ensembl_exons, ensembl_cds, ensembl_warnings = normalize_ensembl_model(
                sp, args.gene_symbol, ensembl_model_raw, ensembl_lookup_method
            )
            for w in ensembl_warnings:
                annotation_warning_rows.append(
                    {
                        "species_input": sp.input_name,
                        "species_canonical": sp.ensembl_species,
                        "source_db": "Ensembl",
                        "warning": w,
                    }
                )
        else:
            annotation_warning_rows.append(
                {
                    "species_input": sp.input_name,
                    "species_canonical": sp.ensembl_species,
                    "source_db": "Ensembl",
                    "warning": f"lookup failed: {ensembl_lookup_method}",
                }
            )

        # NCBI evidence
        ncbi_gene_ev = fetch_ncbi_gene_evidence(sp, args.gene_symbol, args.sleep_between_requests)
        ncbi_prot_ev = fetch_ncbi_refseq_protein_evidence(sp, args.gene_symbol, args.sleep_between_requests)
        ncbi_rescue_rows.append(
            {
                "species_input": sp.input_name,
                "species_canonical": sp.ensembl_species,
                "taxid": sp.taxid,
                "gene_symbol": args.gene_symbol,
                "ncbi_gene_found": ncbi_gene_ev["ncbi_gene_found"],
                "ncbi_gene_id": ncbi_gene_ev["ncbi_gene_id"],
                "ncbi_gene_name": ncbi_gene_ev["ncbi_gene_name"],
                "ncbi_gene_description": ncbi_gene_ev["ncbi_gene_description"],
                "ncbi_gene_taxname": ncbi_gene_ev["ncbi_gene_taxname"],
                "ncbi_eutils_status": ncbi_gene_ev.get("ncbi_eutils_status", ""),
                "refseq_protein_hits": ncbi_prot_ev["refseq_protein_hits"],
                "refseq_protein_example_ids": ncbi_prot_ev["refseq_protein_example_ids"],
                "ncbi_protein_eutils_status": ncbi_prot_ev.get("ncbi_protein_eutils_status", ""),
            }
        )

        # NCBI structured fallback
        ncbi_gene_model: Optional[GeneModel] = None
        ncbi_txs: List[TranscriptModel] = []
        ncbi_exons: List[ExonModel] = []
        ncbi_cds: List[CDSFeature] = []
        datasets_provenance = {
            "datasets_used": "0",
            "datasets_summary_status": "not_requested",
            "datasets_download_status": "not_requested",
            "datasets_package": "",
            "assembly_accession": "",
            "assembly_name": "",
            "assembly_level": "",
            "assembly_status": "",
            "assembly_score": "",
            "assembly_decision_notes": "",
        }

        # The rescue evidence gathered a moment ago is an *input* to the structured
        # route, not a note beside it. The Equus quagga run held the right GeneID
        # (124236178, FGFR2, Equus quagga) in ncbi_rescue_candidates.tsv and never used
        # it; passing it here is what lets a locus be recognised by identifier when its
        # symbol is still LOC-labelled, before any step that assumes a non-empty table.
        expected_gene_ids = [ncbi_gene_ev["ncbi_gene_id"]] if ncbi_gene_ev.get("ncbi_gene_id") else []
        # The source's own spelled-out product name, e.g. "fibroblast growth factor
        # receptor 2". It makes the description route specific to this gene instead of
        # to every locus whose description ends in the same number.
        product_name = ncbi_gene_ev.get("ncbi_gene_description", "") or gene_product_name(
            args.gene_config, args.gene_symbol)
        identifications: List[gid.Identification] = []

        if args.use_ncbi_datasets:
            ncbi_model_result, ncbi_model_warnings, datasets_provenance, assembly_selection_rows = fetch_ncbi_model_from_datasets(
                sp,
                args.gene_symbol,
                workdir,
                datasets_bin,
                args.datasets_include,
                args.datasets_summary_timeout,
                args.datasets_download_timeout,
                args.sleep_between_requests,
                expected_gene_ids=expected_gene_ids,
                paralog_panel=paralog_panel,
                identification_out=identifications,
                product_name=product_name,
            )
            ncbi_assembly_selection_rows.extend(assembly_selection_rows)
            if datasets_provenance.get("assembly_accession", ""):
                ncbi_assembly_selected_rows.append(
                    {
                        "species_input": sp.input_name,
                        "species_canonical": sp.ensembl_species,
                        "taxid": datasets_provenance.get("taxid") or sp.taxid,
                        "assembly_accession": datasets_provenance.get("assembly_accession", ""),
                        "assembly_name": datasets_provenance.get("assembly_name", ""),
                        "assembly_level": datasets_provenance.get("assembly_level", ""),
                        "assembly_status": datasets_provenance.get("assembly_status", ""),
                        "assembly_score": datasets_provenance.get("assembly_score", ""),
                        "assembly_decision_notes": datasets_provenance.get("assembly_decision_notes", ""),
                    }
                )
            for w in ncbi_model_warnings:
                annotation_warning_rows.append(
                    {
                        "species_input": sp.input_name,
                        "species_canonical": sp.ensembl_species,
                        "source_db": "NCBI",
                        "warning": w,
                    }
                )
            if ncbi_model_result is not None:
                ncbi_gene_model, ncbi_txs, ncbi_exons, ncbi_cds = ncbi_model_result

        per_species_ncbi_provenance[sp.input_name] = datasets_provenance

        # The recovery story for this species, in the order it happened. This is what
        # replaces the traceback as the run's explanation of itself.
        outcome = recovery.SpeciesOutcome(
            species_id=sp.ensembl_species, species_input=sp.input_name,
            gene_symbol=args.gene_symbol, taxid=sp.taxid,
            accepted_scientific_name=sp.ncbi_species)
        outcome.assembly_accession = datasets_provenance.get("assembly_accession", "")
        outcome.assembly_status = datasets_provenance.get(
            "assembly_selection_status", datasets_provenance.get(
                "datasets_summary_status", ""))

        outcome.record("ensembl_preferred_source",
                       "model_found" if ensembl_gene_model else "no_model",
                       ensembl_lookup_method, len(ensembl_txs))
        outcome.record("ncbi_annotated_assembly",
                       "model_found" if ncbi_gene_model else "no_model",
                       datasets_provenance.get("assembly_selection_detail",
                                               outcome.assembly_status),
                       len(ncbi_txs))
        if expected_gene_ids:
            outcome.record("gene_identifier_lookup", "gene_id_found",
                           f"NCBI GeneID {expected_gene_ids[0]} "
                           f"({ncbi_gene_ev.get('ncbi_gene_name')}, "
                           f"{ncbi_gene_ev.get('ncbi_gene_taxname')})")
        for ident in identifications:
            outcome.record("annotation_gene_identification", ident.status, ident.detail)
            for cand in ident.candidates:
                gene_candidate_rows.append({
                    "species_input": sp.input_name,
                    "species_canonical": sp.ensembl_species,
                    "gene_symbol": args.gene_symbol,
                    **cand.as_row(),
                })
            if ident.accepted is not None:
                outcome.accepted_route = ident.accepted.route

        # Compare and resolve
        selected_source, selected_confidence, detailed_reason, comparison_row, compare_warnings = compare_models(
            ensembl_gene_model, ensembl_txs, ensembl_exons,
            ncbi_gene_model, ncbi_txs, ncbi_exons,
            preferred_source=sp.preferred_source,
        )

        source_comparison_rows.append(
            {
                "species_input": sp.input_name,
                "species_canonical": sp.ensembl_species,
                "gene_symbol": args.gene_symbol,
                **comparison_row,
                "selected_source": selected_source,
                "selected_confidence": selected_confidence,
                "selection_reason_detailed": detailed_reason,
            }
        )

        for w in compare_warnings:
            source_conflict_rows.append(
                {
                    "species_input": sp.input_name,
                    "species_canonical": sp.ensembl_species,
                    "gene_symbol": args.gene_symbol,
                    "selected_source": selected_source,
                    "conflict_level": comparison_row.get("conflict_level", ""),
                    "warning": w,
                }
            )

        selected_gene: Optional[GeneModel] = None
        selected_txs: List[TranscriptModel] = []
        selected_exons: List[ExonModel] = []
        selected_cds: List[CDSFeature] = []
        alt_source = "yes" if (ensembl_gene_model and ncbi_gene_model) else "no"

        if selected_source == "Ensembl":
            selected_gene, selected_txs, selected_exons, selected_cds = ensembl_gene_model, ensembl_txs, ensembl_exons, ensembl_cds
        elif selected_source == "NCBI":
            selected_gene, selected_txs, selected_exons, selected_cds = ncbi_gene_model, ncbi_txs, ncbi_exons, ncbi_cds

        final_model_status = "no_model_found"
        selected_gene_conf = "none"
        if selected_gene:
            genes_rows.append(asdict(selected_gene))
            transcripts_rows.extend({**asdict(tx), "species_input": sp.input_name} for tx in selected_txs)
            exons_rows.extend({**asdict(ex), "species_input": sp.input_name} for ex in selected_exons)
            cds_rows.extend(asdict(cds) for cds in selected_cds)
            final_model_status = selected_gene.model_status
            selected_gene_conf = selected_gene.model_confidence

        # Conclude this species. Deciding here — before any step that assumes a
        # non-empty table — is the point: the run states what it recovered and why,
        # instead of leaving a later step to infer it from an empty file.
        outcome.n_genes = 1 if selected_gene else 0
        outcome.n_transcripts = len(selected_txs)
        outcome.n_exons = len(selected_exons)
        outcome.n_cds_features = len(selected_cds)
        outcome.n_translated_proteins = count_translated_transcripts(selected_txs)

        if selected_gene and outcome.n_transcripts:
            if not outcome.n_translated_proteins:
                outcome.conclude(
                    recovery.NO_VALID_TRANSLATED_CDS,
                    f"the locus has {outcome.n_transcripts} transcript(s) but none "
                    "carries a translatable CDS")
            elif outcome.accepted_route in gid.ROUTES_NEEDING_DISCRIMINATION:
                # Recovered, but on evidence weaker than an official symbol. Analysed
                # as a reviewable model rather than silently as an equal.
                outcome.conclude(
                    recovery.REVIEW_REQUIRED,
                    f"the locus was identified via {outcome.accepted_route} rather than "
                    "an official gene symbol")
            else:
                outcome.conclude(recovery.MODELS_AVAILABLE)
        else:
            identification_statuses = [i.status for i in identifications]
            if identification_statuses:
                outcome.conclude(
                    recovery.status_from_identification(identification_statuses[-1]),
                    identifications[-1].detail)
            else:
                status = recovery.status_from_assembly(outcome.assembly_status)
                if status == recovery.RESCUE_REQUIRED:
                    status = (recovery.SOURCE_UNAVAILABLE
                              if "HTTPError" in ensembl_lookup_method
                              else recovery.ANNOTATION_NOT_FOUND)
                outcome.conclude(
                    status,
                    datasets_provenance.get("assembly_selection_detail", ""))
        outcomes.append(outcome)

        species_status_rows.append(
            {
                "species_input": sp.input_name,
                "species_canonical": sp.ensembl_species,
                "taxid": sp.taxid,
                "gene_symbol": args.gene_symbol,
                "selected_source": selected_source,
                "final_model_status": final_model_status,
                "selected_gene_confidence": selected_gene_conf,
                "comparison_conflict_level": comparison_row.get("conflict_level", ""),
                "alternative_source_available": alt_source,
                "ensembl_lookup_method": ensembl_lookup_method,
                "ncbi_gene_found": ncbi_gene_ev["ncbi_gene_found"],
                "ncbi_refseq_protein_hits": ncbi_prot_ev["refseq_protein_hits"],
                "datasets_summary_status": datasets_provenance.get("datasets_summary_status", ""),
                "datasets_download_status": datasets_provenance.get("datasets_download_status", ""),
                "selected_ncbi_accession": datasets_provenance.get("assembly_accession", ""),
                "selected_ncbi_assembly_name": datasets_provenance.get("assembly_name", ""),
            }
        )

        selection_rows.append(
            {
                "species_input": sp.input_name,
                "species_canonical": sp.ensembl_species,
                "gene_symbol": args.gene_symbol,
                "preferred_source": sp.preferred_source,
                "selected_source": selected_source,
                "selected_confidence": selected_confidence,
                "ensembl_model_confidence": ensembl_gene_model.model_confidence if ensembl_gene_model else "none",
                "ncbi_model_confidence": ncbi_gene_model.model_confidence if ncbi_gene_model else "none",
                "ensembl_translation_count": str(count_translated_transcripts(ensembl_txs)),
                "ncbi_translation_count": str(count_translated_transcripts(ncbi_txs)),
                "ensembl_max_protein_length": to_str(max_protein_length(ensembl_txs)),
                "ncbi_max_protein_length": to_str(max_protein_length(ncbi_txs)),
                "ensembl_available": "1" if ensembl_gene_model else "0",
                "ncbi_structured_available": "1" if ncbi_gene_model else "0",
                "selection_reason_detailed": detailed_reason,
            }
        )

        raw_source_map_rows.append(
            {
                "species_input": sp.input_name,
                "species_canonical": sp.ensembl_species,
                "taxid": sp.taxid,
                "ensembl_species": sp.ensembl_species,
                "ncbi_species": sp.ncbi_species,
                "preferred_source": sp.preferred_source,
                "assembly_preference": sp.assembly_preference,
                "selected_source": selected_source,
            }
        )

    # Final metadata
    run_meta["per_species_ncbi_provenance"] = per_species_ncbi_provenance

    # Internal consistency checks on selected outputs
    # Structural checks over the recovered tables, plus checks that cannot pass
    # vacuously. The failed run reported six of six PASS over four empty tables, so the
    # one artefact meant to catch an inconsistent result endorsed an empty one.
    internal_checks_rows = recovery.consistency_checks(
        len(genes_rows), len(transcripts_rows), len(exons_rows), len(cds_rows))
    if genes_rows or transcripts_rows:
        internal_checks_rows += run_internal_consistency_checks(
            genes_rows, transcripts_rows, exons_rows)
    else:
        internal_checks_rows += [{
            "check_name": name,
            "status": "NOT_APPLICABLE",
            "affected_species": "",
            "details": "no rows to check",
        } for name in ("transcripts_reference_existing_genes",
                       "exons_reference_existing_transcripts",
                       "exon_coordinates_valid",
                       "exon_rank_gaps_absent",
                       "exon_coordinate_order_matches_rank",
                       "selected_genes_have_transcripts")]

    contract = recovery.CollectionContract(gene_symbol=args.gene_symbol,
                                          outcomes=outcomes)
    contract_path = contract.write(outdir)

    # Write outputs
    write_tsv(
        genes_rows,
        outdir / "genes.tsv",
        [
            "source_db", "species_input", "species_canonical", "taxid", "assembly_accession", "assembly_name",
            "gene_symbol_requested", "gene_symbol_found", "gene_id_source", "gene_biotype",
            "chrom", "start", "end", "strand", "model_status", "model_confidence", "lookup_method", "internal_gene_id",
        ],
    )

    write_tsv(
        transcripts_rows,
        outdir / "transcripts.tsv",
        [
            "species_input", "source_db", "species_canonical", "gene_id_internal", "transcript_id_source", "transcript_name",
            "transcript_biotype", "translation_id_source", "protein_length_aa", "protein_length_source", "is_canonical_source",
            "support_level", "completeness_flags", "transcript_model_confidence",
            "chrom", "start", "end", "strand", "internal_transcript_id",
        ],
    )

    write_tsv(
        exons_rows,
        outdir / "exons.tsv",
        [
            "species_input", "source_db", "species_canonical", "transcript_id_internal", "exon_id_source", "exon_rank",
            "chrom", "start", "end", "strand", "phase", "end_phase", "internal_exon_id",
        ],
    )


    write_tsv(
        cds_rows,
        outdir / "cds_features.tsv",
        [
            "species_input", "source_db", "species_canonical", "transcript_id_internal", "transcript_id_source",
            "translation_id_source", "cds_id_source", "cds_rank", "chrom", "start", "end", "strand", "phase",
            "cds_length_bp", "cds_offset_start_0based", "cds_offset_end_0based", "protein_start_aa", "protein_end_aa",
            "parent_feature_id", "internal_cds_id", "coordinate_source", "confidence", "warning",
        ],
    )

    write_tsv(
        species_status_rows,
        outdir / "species_status.tsv",
        [
            "species_input", "species_canonical", "taxid", "gene_symbol", "selected_source",
            "final_model_status", "selected_gene_confidence", "comparison_conflict_level",
            "alternative_source_available", "ensembl_lookup_method",
            "ncbi_gene_found", "ncbi_refseq_protein_hits",
            "datasets_summary_status", "datasets_download_status",
            "selected_ncbi_accession", "selected_ncbi_assembly_name",
        ],
    )

    write_tsv(
        selection_rows,
        outdir / "selection_decisions.tsv",
        [
            "species_input", "species_canonical", "gene_symbol", "preferred_source", "selected_source",
            "selected_confidence", "ensembl_model_confidence", "ncbi_model_confidence",
            "ensembl_translation_count", "ncbi_translation_count",
            "ensembl_max_protein_length", "ncbi_max_protein_length",
            "ensembl_available", "ncbi_structured_available",
            "selection_reason_detailed",
        ],
    )

    write_tsv(
        source_conflict_rows,
        outdir / "source_conflicts.tsv",
        [
            "species_input", "species_canonical", "gene_symbol", "selected_source", "conflict_level", "warning",
        ],
    )

    write_tsv(
        source_comparison_rows,
        outdir / "source_comparison.tsv",
        [
            "species_input", "species_canonical", "gene_symbol",
            "ensembl_available", "ncbi_available",
            "same_gene_symbol", "same_biotype",
            "same_chrom_label", "seqid_relation",
            "gene_overlap_fraction", "gene_overlap_reason",
            "translation_support_ensembl", "translation_support_ncbi",
            "max_protein_len_ensembl", "max_protein_len_ncbi",
            "long_transcript_count_ensembl", "long_transcript_count_ncbi",
            "transcript_count_ensembl", "transcript_count_ncbi",
            "exon_count_ensembl", "exon_count_ncbi",
            "transcript_count_diff", "exon_count_diff",
            "conflict_level", "ensembl_source_score", "ncbi_source_score", "source_score_difference",
            "preferred_source_applied", "comparison_decision_reason",
            "selected_source", "selected_confidence", "selection_reason_detailed",
        ],
    )

    write_tsv(
        annotation_warning_rows,
        outdir / "annotation_warnings.tsv",
        ["species_input", "species_canonical", "source_db", "warning"],
    )

    write_tsv(
        raw_source_map_rows,
        outdir / "raw_source_map.tsv",
        [
            "species_input", "species_canonical", "taxid", "ensembl_species", "ncbi_species",
            "preferred_source", "assembly_preference", "selected_source",
        ],
    )

    write_tsv(
        ncbi_rescue_rows,
        outdir / "ncbi_rescue_candidates.tsv",
        [
            "species_input", "species_canonical", "taxid", "gene_symbol",
            "ncbi_gene_found", "ncbi_gene_id", "ncbi_gene_name",
            "ncbi_gene_description", "ncbi_gene_taxname", "ncbi_eutils_status",
            "refseq_protein_hits", "refseq_protein_example_ids", "ncbi_protein_eutils_status",
        ],
    )

    write_tsv(
        ncbi_assembly_selection_rows,
        outdir / "ncbi_assembly_selection.tsv",
        [
            "species_input", "species_canonical", "taxid",
            # The outcome of the selection as a whole, repeated on every candidate row,
            # so that a table with rows still says why it has no selected assembly.
            "selection_status", "query_term",
            "accession", "assembly_name", "assembly_level", "refseq_category",
            "assembly_source", "assembly_status", "annotated", "organism_name",
            "annotation_name", "annotation_release_date", "paired_accession",
            "tax_id", "taxon_match", "decision", "rejection_reason",
            "assembly_score", "assembly_decision_notes",
        ],
    )

    write_tsv(
        ncbi_assembly_selected_rows,
        outdir / "ncbi_assembly_selected.tsv",
        [
            "species_input", "species_canonical", "taxid",
            "assembly_accession", "assembly_name", "assembly_level",
            "assembly_status", "assembly_score", "assembly_decision_notes",
        ],
    )

    write_tsv(
        internal_checks_rows,
        outdir / "internal_consistency_checks.tsv",
        ["check_name", "status", "affected_species", "details"],
    )

    write_tsv(
        gene_candidate_rows,
        outdir / "gene_candidates.tsv",
        ["species_input", "species_canonical", "gene_symbol", "source_gene_id",
         "source_symbol", "description", "seqid", "strand", "start", "end", "biotype",
         "dbxrefs", "synonyms", "transcript_count", "protein_count",
         "identification_route", "orthology_evidence", "similarity_evidence",
         "paralog_discrimination_evidence", "decision", "reason"],
    )

    write_json(run_meta, outdir / "run_metadata.json")

    print(f"[OK] genes.tsv rows: {len(genes_rows)}")
    print(f"[OK] transcripts.tsv rows: {len(transcripts_rows)}")
    print(f"[OK] exons.tsv rows: {len(exons_rows)}")
    print(f"[OK] cds_features.tsv rows: {len(cds_rows)}")
    print(f"[OK] species_status.tsv rows: {len(species_status_rows)}")
    print(f"[OK] selection_decisions.tsv rows: {len(selection_rows)}")
    print(f"[OK] source_conflicts.tsv rows: {len(source_conflict_rows)}")
    print(f"[OK] source_comparison.tsv rows: {len(source_comparison_rows)}")
    print(f"[OK] annotation_warnings.tsv rows: {len(annotation_warning_rows)}")
    print(f"[OK] raw_source_map.tsv rows: {len(raw_source_map_rows)}")
    print(f"[OK] ncbi_rescue_candidates.tsv rows: {len(ncbi_rescue_rows)}")
    print(f"[OK] ncbi_assembly_selection.tsv rows: {len(ncbi_assembly_selection_rows)}")
    print(f"[OK] ncbi_assembly_selected.tsv rows: {len(ncbi_assembly_selected_rows)}")
    print(f"[OK] internal_consistency_checks.tsv rows: {len(internal_checks_rows)}")
    print(f"[OK] gene_candidates.tsv rows: {len(gene_candidate_rows)}")
    print(f"[OK] run_metadata.json written to: {outdir / 'run_metadata.json'}")
    print(f"[OK] collection_status.json: {contract.status} -> {contract_path}")
    for outcome in outcomes:
        print(f"[STATUS] {outcome.species_id}: {outcome.status} — {outcome.message()}")
    if not contract.usable:
        # Non-zero, because a run that recovered no model has not succeeded. The reason
        # is on the contract, so the next step can report it instead of raising.
        print(f"[FAIL] {contract.message()}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
