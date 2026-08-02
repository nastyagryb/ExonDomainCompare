#!/usr/bin/env python3
"""
protein_lookup.py

Shared protein-sequence retrieval logic for the FGFR2 IIIb/IIIc pipeline.

Why this module exists
----------------------
Historically, Step 5b (protein validation) and Step 6 (protein export) used two
different, incompatible protein-retrieval strategies. Step 6 robustly resolved
NCBI protein sequences from the NCBI ``datasets`` cache (GFF3 transcript->protein
mapping + protein.faa indexing + product/species/length rescue), while Step 5b
only looked up sequences by translation ID in a flat FASTA index and otherwise
fell back to Ensembl REST. As a result Step 5b reported hundreds of
``protein_sequence_unavailable`` rows for NCBI candidates that Step 6 could
export without problems.

This module centralises the retrieval so that both steps use *exactly* the same
logic and report consistent, informative lookup statuses.

Design
------
``ProteinCacheIndex`` is built once from a cache directory and reused for every
lookup. ``lookup_protein`` resolves a single transcript/protein row to a
``ProteinLookupResult`` describing the sequence, the database it came from, the
method used, and a confidence level.

The lookup order is deterministic and conservative:
  1. Exact protein-accession match (NCBI XP_/NP_ or Ensembl translation ID).
  2. Transcript-linked match (GFF3 CDS Parent -> protein accession -> FASTA).
  3. Cache product/species/length rescue (clearly reported, lower confidence).
  4. Ensembl REST fetch for ENSP* translation IDs (optional, network).

It never silently "finds" an ambiguous sequence: rescue and REST matches are
flagged with explicit lower-confidence method labels.
"""
from __future__ import annotations

import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import unquote
from urllib.request import Request, urlopen


# Protein-sequence status values.
STATUS_FOUND_EXACT = "protein_found_exact_accession"
STATUS_FOUND_TX_LINKED = "protein_found_transcript_linked"
STATUS_FOUND_CACHE_FALLBACK = "protein_found_cache_fallback"
STATUS_FOUND_ENSEMBL_REST = "protein_found_ensembl_rest"
STATUS_UNAVAIL_NO_ACCESSION = "protein_unavailable_no_accession"
STATUS_UNAVAIL_NOT_IN_CACHE = "protein_unavailable_not_in_cache"
STATUS_UNAVAIL_AMBIGUOUS = "protein_unavailable_ambiguous_mapping"

FOUND_STATUSES = {
    STATUS_FOUND_EXACT,
    STATUS_FOUND_TX_LINKED,
    STATUS_FOUND_CACHE_FALLBACK,
    STATUS_FOUND_ENSEMBL_REST,
}

CONFIDENCE_BY_STATUS = {
    STATUS_FOUND_EXACT: "high",
    STATUS_FOUND_TX_LINKED: "high",
    STATUS_FOUND_ENSEMBL_REST: "high",
    STATUS_FOUND_CACHE_FALLBACK: "medium",
    STATUS_UNAVAIL_NO_ACCESSION: "none",
    STATUS_UNAVAIL_NOT_IN_CACHE: "none",
    STATUS_UNAVAIL_AMBIGUOUS: "none",
}


@dataclass
class ProteinLookupResult:
    sequence: str = ""
    source_db: str = ""
    protein_sequence_status: str = STATUS_UNAVAIL_NO_ACCESSION
    lookup_method: str = "none"
    lookup_confidence: str = "none"
    matched_accession: str = ""
    detail: str = ""

    @property
    def found(self) -> bool:
        return bool(self.sequence) and self.protein_sequence_status in FOUND_STATUSES

    @property
    def length_aa(self) -> int:
        return len(self.sequence.replace("*", "")) if self.sequence else 0


@dataclass
class _FastaRecord:
    accession: str
    accession_no_version: str
    header: str
    sequence: str
    path: str
    product_lower: str
    bracket_species_lower: str


@dataclass
class _GffProteinMap:
    protein_accession: str = ""
    product: str = ""
    source_gff: str = ""


# ----------------------------- generic helpers -----------------------------

def clean_accession(x: object) -> str:
    s = str(x or "").strip()
    if not s or s.lower() in {"nan", "none", "null", "na"}:
        return ""
    s = s.split()[0]
    for prefix in ("rna-", "cds-", "protein-", "transcript:", "protein:", "RefSeq:", "Genbank:", "NCBI:"):
        if s.startswith(prefix):
            s = s[len(prefix):]
    return s.strip()


def no_version(acc: object) -> str:
    a = clean_accession(acc)
    return re.sub(r"\.\d+$", "", a)


def parse_int_maybe(x: object) -> Optional[int]:
    s = str(x or "").strip()
    if not s or s.lower() in {"nan", "none", "null", "na"}:
        return None
    try:
        return int(float(s))
    except Exception:
        return None


def _parse_attrs(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in attr.strip().split(";"):
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k] = unquote(v)
    return out


def _parse_fasta(path: Path):
    header = None
    seq_parts: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).replace(" ", "")
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts).replace(" ", "")


def _bracket_species(header: str) -> str:
    m = re.search(r"\[([^\]]+)\]\s*$", header)
    return m.group(1).lower() if m else ""


def _product_text(header: str) -> str:
    h = re.sub(r"\s*\[[^\]]+\]\s*$", "", header).strip()
    parts = h.split(maxsplit=1)
    return parts[1].lower() if len(parts) > 1 else ""


def _choose_gff_paths(cache: Path) -> List[Path]:
    paths = list(cache.rglob("*.gff")) + list(cache.rglob("*.gff3"))

    def score(p: Path) -> Tuple[int, int, str]:
        name = p.name.lower()
        s = 0
        if name == "genomic.gff":
            s += 100
        if "ncbi_dataset" in str(p):
            s += 20
        try:
            size = p.stat().st_size
        except OSError:
            size = 0
        return (-s, -size, str(p))

    return sorted(paths, key=score)


# ----------------------------- cache index -----------------------------

_FGFR2_PRODUCT_HINTS = ("fibroblast growth factor receptor 2", "fgfr2")
_FGFR2_LINE_HINTS = ("FGFR2", "fibroblast growth factor receptor 2", "rowth factor receptor 2")


class ProteinCacheIndex:
    """Index protein sequences from an NCBI ``datasets`` cache directory.

    Provides exact-accession lookup, transcript->protein mapping (from GFF3),
    and species/product/length rescue identical to the Step 6 export logic.

    Performance
    -----------
    NCBI ``datasets`` caches are large (tens of GB of genomic GFF3). To keep the
    pipeline fast and reproducible we:
      * index only protein FASTA records that are either explicitly requested
        (``wanted_accessions``) or are FGFR2 products (needed for rescue);
      * build the transcript->protein GFF3 mapping lazily and only scan GFF lines
        that mention FGFR2 (a cheap substring pre-filter), since exact-accession
        lookups resolve the overwhelming majority of candidates already.
    """

    def __init__(self, cache: Optional[Path], wanted_accessions: Optional[set] = None,
                 wanted_transcripts: Optional[set] = None):
        self.cache = Path(cache) if cache else None
        self.tx_to_protein: Dict[str, _GffProteinMap] = {}
        self.by_acc: Dict[str, _FastaRecord] = {}
        self.by_species: Dict[str, List[_FastaRecord]] = {}
        self._gff_built = False
        self._fasta_records = 0
        self.wanted_accessions = {clean_accession(a) for a in (wanted_accessions or set()) if a}
        self.wanted_accessions |= {no_version(a) for a in list(self.wanted_accessions)}
        self.wanted_transcripts = {clean_accession(t) for t in (wanted_transcripts or set()) if t}
        self.wanted_transcripts |= {no_version(t) for t in list(self.wanted_transcripts)}
        if self.cache and self.cache.exists():
            self._index_fastas()

    def _keep_fasta(self, accession: str, accession_no_version: str, product_lower: str) -> bool:
        if not self.wanted_accessions:
            return True
        if accession in self.wanted_accessions or accession_no_version in self.wanted_accessions:
            return True
        return any(h in product_lower for h in _FGFR2_PRODUCT_HINTS)

    def _ensure_gff(self) -> None:
        if self._gff_built or not self.cache:
            return
        self._gff_built = True
        self._index_gff()

    def _index_gff(self) -> None:
        assert self.cache is not None
        for gff in _choose_gff_paths(self.cache):
            try:
                with gff.open("r", encoding="utf-8", errors="replace") as f:
                    for line in f:
                        if not line or line.startswith("#"):
                            continue
                        # Cheap substring pre-filter: only FGFR2-related lines matter.
                        if not any(h in line for h in _FGFR2_LINE_HINTS):
                            continue
                        cols = line.rstrip("\n").split("\t")
                        if len(cols) < 9 or cols[2].lower() not in {"cds", "protein"}:
                            continue
                        attrs = _parse_attrs(cols[8])
                        raw_protein = attrs.get("protein_id") or attrs.get("Name") or attrs.get("ID") or attrs.get("Dbxref", "")
                        protein = clean_accession(raw_protein)
                        m = re.search(r"(?:Genbank|RefSeq|NCBI):([^,;]+)", raw_protein)
                        if m:
                            protein = clean_accession(m.group(1))
                        if not re.match(r"^[XN]P_", protein):
                            m2 = re.search(r"([XN]P_\d+(?:\.\d+)?)", cols[8])
                            if m2:
                                protein = m2.group(1)
                        if not protein:
                            continue
                        product = attrs.get("product", "")
                        for parent in attrs.get("Parent", "").split(","):
                            parent = clean_accession(parent)
                            for key in {parent, no_version(parent)}:
                                if key and key not in self.tx_to_protein:
                                    self.tx_to_protein[key] = _GffProteinMap(protein, product, str(gff))
            except OSError:
                continue

    def _index_fastas(self) -> None:
        assert self.cache is not None
        fasta_paths = [p for p in self.cache.rglob("protein.faa")]
        if not fasta_paths:
            fasta_paths = [p for p in self.cache.rglob("*") if p.is_file() and p.suffix.lower() in {".faa", ".fa", ".fasta"}]
        for fasta in sorted(fasta_paths):
            for header, seq in _parse_fasta(fasta):
                first = clean_accession(header.split()[0] if header else "")
                if not first:
                    continue
                acc_nv = no_version(first)
                product_lower = _product_text(header)
                if not self._keep_fasta(first, acc_nv, product_lower):
                    continue
                rec = _FastaRecord(
                    accession=first,
                    accession_no_version=acc_nv,
                    header=header,
                    sequence=seq,
                    path=str(fasta),
                    product_lower=product_lower,
                    bracket_species_lower=_bracket_species(header),
                )
                for k in {rec.accession, rec.accession_no_version}:
                    if k and k not in self.by_acc:
                        self.by_acc[k] = rec
                # Only FGFR2 product records are useful for species/product rescue.
                if any(h in product_lower for h in _FGFR2_PRODUCT_HINTS):
                    self.by_species.setdefault(rec.bracket_species_lower, []).append(rec)
                self._fasta_records += 1

    # ---- composition reporting ----
    def composition(self) -> Dict[str, object]:
        return {
            "cache_path": str(self.cache) if self.cache else "",
            "fasta_records_indexed": self._fasta_records,
            "protein_accession_keys": len(self.by_acc),
            "fgfr2_species_groups": len(self.by_species),
            "gff_built": self._gff_built,
            "transcript_protein_mappings": len(self.tx_to_protein),
        }

    # ---- matching primitives ----
    def by_accession(self, acc: str) -> Optional[_FastaRecord]:
        for key in (clean_accession(acc), no_version(acc)):
            if key and key in self.by_acc:
                return self.by_acc[key]
        return None

    def transcript_protein(self, tx: str) -> _GffProteinMap:
        self._ensure_gff()
        for key in (clean_accession(tx), no_version(tx)):
            if key and key in self.tx_to_protein:
                return self.tx_to_protein[key]
        return _GffProteinMap()

    def product_rescue(self, species_input: str, species_canonical: str, product: str, expected_len: Optional[int]) -> Optional[_FastaRecord]:
        prod_l = (product or "").lower().strip()
        candidates: List[_FastaRecord] = []
        for sp in _species_bracket_candidates(species_input, species_canonical):
            for bracket, recs in self.by_species.items():
                if bracket and (bracket == sp or bracket.startswith(sp) or sp.startswith(bracket)):
                    candidates.extend(recs)
        if not candidates:
            return None
        scored: List[Tuple[int, _FastaRecord]] = []
        for rec in candidates:
            p = rec.product_lower
            if "fibroblast growth factor receptor 2" not in p and "fgfr2" not in p:
                continue
            score = 40
            if prod_l:
                if p == prod_l:
                    score += 100
                elif prod_l in p or p in prod_l:
                    score += 70
            iso = re.search(r"isoform\s+([A-Za-z0-9]+)", prod_l)
            if iso:
                if re.search(rf"isoform\s+{re.escape(iso.group(1))}(?:\s|$)", p):
                    score += 50
                else:
                    score -= 40
            score += _length_score(len(rec.sequence), expected_len)
            scored.append((score, rec))
        if not scored:
            return None
        scored.sort(key=lambda x: (x[0], -abs(len(x[1].sequence) - expected_len) if expected_len else 0), reverse=True)
        return scored[0][1] if scored[0][0] >= 50 else None


def _species_bracket_candidates(species_input: str, species_canonical: str) -> List[str]:
    vals: List[str] = []
    if species_input:
        vals.append(species_input.lower())
    if species_canonical:
        vals.append(species_canonical.replace("_", " ").lower())
    for v in list(vals):
        toks = v.split()
        if len(toks) >= 2:
            vals.append(" ".join(toks[:2]))
    out, seen = [], set()
    for v in vals:
        if v and v not in seen:
            out.append(v)
            seen.add(v)
    return out


def _length_score(seq_len: int, expected_len: Optional[int]) -> int:
    if expected_len is None:
        return 0
    diff = abs(seq_len - expected_len)
    if diff == 0:
        return 60
    if diff <= 2:
        return 50
    if diff <= 10:
        return 35
    if diff <= 30:
        return 15
    return -50


# ----------------------------- Ensembl REST -----------------------------

def fetch_ensembl_protein(translation_id: str, retries: int = 3, sleep: float = 0.4, timeout: int = 30,
                          user_agent: str = "FGFR2-boundary-mapping-bachelor-thesis/1.0") -> Tuple[str, str]:
    tid = clean_accession(translation_id)
    if not tid:
        return "", "missing_translation_id"
    url = f"https://rest.ensembl.org/sequence/id/{tid}?type=protein"
    headers = {"Content-Type": "text/plain", "User-Agent": user_agent}
    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers=headers)
            with urlopen(req, timeout=timeout) as resp:
                text = resp.read().decode("utf-8", errors="replace").strip()
            seq = re.sub(r"\s+", "", text)
            if seq and not seq.startswith("{") and re.fullmatch(r"[A-Za-z*]+", seq) and len(seq) > 20:
                return seq.replace("*", ""), "ensembl_rest_matched"
            return "", f"ensembl_rest_unexpected_response:{text[:80]}"
        except HTTPError as e:
            if e.code in {429, 500, 502, 503, 504} and attempt < retries:
                time.sleep(sleep * attempt)
                continue
            return "", f"ensembl_rest_http_{e.code}"
        except URLError as e:
            if attempt < retries:
                time.sleep(sleep * attempt)
                continue
            return "", f"ensembl_rest_urlerror:{e.reason}"
        except Exception as e:  # pragma: no cover - network edge cases
            if attempt < retries:
                time.sleep(sleep * attempt)
                continue
            return "", f"ensembl_rest_error:{type(e).__name__}:{e}"
    return "", "ensembl_rest_failed"


# ----------------------------- unified lookup -----------------------------

def lookup_protein(
    index: ProteinCacheIndex,
    *,
    source_db: str,
    transcript_id: str,
    translation_id: str,
    species_input: str = "",
    species_canonical: str = "",
    expected_length_aa: Optional[int] = None,
    allow_ensembl_rest: bool = True,
    allow_cache_rescue: bool = True,
    ensembl_sleep: float = 0.4,
    ensembl_timeout: int = 30,
) -> ProteinLookupResult:
    """Resolve a single transcript/protein to a sequence using shared logic.

    The same retrieval path is used by Step 5b (validation) and Step 6 (export)
    so that protein availability is consistent across the pipeline.
    """
    src = (source_db or "").strip()
    src_upper = src.upper()
    tx = clean_accession(transcript_id)
    translation = clean_accession(translation_id)

    # NCBI path: cache-based resolution (exact accession, transcript-linked, rescue).
    if src_upper == "NCBI" or (not src and (tx.startswith("XM_") or tx.startswith("NM_") or translation.startswith("XP_") or translation.startswith("NP_"))):
        # 1) exact accession (explicit translation ID) - does NOT need the GFF.
        if translation:
            rec = index.by_accession(translation)
            if rec:
                return ProteinLookupResult(
                    sequence=rec.sequence, source_db="NCBI",
                    protein_sequence_status=STATUS_FOUND_EXACT,
                    lookup_method="ncbi_exact_accession",
                    lookup_confidence=CONFIDENCE_BY_STATUS[STATUS_FOUND_EXACT],
                    matched_accession=rec.accession, detail=rec.header,
                )
        # 2) transcript-linked protein accession from GFF3 (lazily triggers GFF parse).
        mapped = index.transcript_protein(tx)
        requested = translation or mapped.protein_accession
        requested_clean = clean_accession(requested)
        product = mapped.product
        if mapped.protein_accession:
            rec = index.by_accession(mapped.protein_accession)
            if rec:
                return ProteinLookupResult(
                    sequence=rec.sequence, source_db="NCBI",
                    protein_sequence_status=STATUS_FOUND_TX_LINKED,
                    lookup_method="ncbi_transcript_linked_gff3",
                    lookup_confidence=CONFIDENCE_BY_STATUS[STATUS_FOUND_TX_LINKED],
                    matched_accession=rec.accession, detail=rec.header,
                )
        # 3) product/species/length rescue
        if allow_cache_rescue:
            rec = index.product_rescue(species_input, species_canonical, product, expected_length_aa)
            if rec:
                return ProteinLookupResult(
                    sequence=rec.sequence, source_db="NCBI",
                    protein_sequence_status=STATUS_FOUND_CACHE_FALLBACK,
                    lookup_method="ncbi_product_species_length_rescue",
                    lookup_confidence=CONFIDENCE_BY_STATUS[STATUS_FOUND_CACHE_FALLBACK],
                    matched_accession=rec.accession,
                    detail="exact_accession_not_found_product_rescue_used",
                )
        if not requested_clean:
            return ProteinLookupResult(
                source_db="NCBI",
                protein_sequence_status=STATUS_UNAVAIL_NO_ACCESSION,
                lookup_method="ncbi_no_accession",
                lookup_confidence="none",
                detail="no_translation_id_and_no_gff3_protein_accession",
            )
        return ProteinLookupResult(
            source_db="NCBI",
            protein_sequence_status=STATUS_UNAVAIL_NOT_IN_CACHE,
            lookup_method="ncbi_not_in_cache",
            lookup_confidence="none",
            detail=f"requested={requested_clean}_not_found_in_cache",
        )

    # Ensembl / other path.
    # Exact match in cache first (rare, but supports custom local FASTA caches).
    if translation:
        rec = index.by_accession(translation)
        if rec:
            return ProteinLookupResult(
                sequence=rec.sequence, source_db=src or "Ensembl",
                protein_sequence_status=STATUS_FOUND_EXACT,
                lookup_method="cache_exact_accession",
                lookup_confidence=CONFIDENCE_BY_STATUS[STATUS_FOUND_EXACT],
                matched_accession=rec.accession, detail=rec.header,
            )
    if translation and allow_ensembl_rest:
        seq, detail = fetch_ensembl_protein(translation, sleep=ensembl_sleep, timeout=ensembl_timeout)
        if seq:
            return ProteinLookupResult(
                sequence=seq, source_db=src or "Ensembl",
                protein_sequence_status=STATUS_FOUND_ENSEMBL_REST,
                lookup_method="ensembl_rest_translation",
                lookup_confidence=CONFIDENCE_BY_STATUS[STATUS_FOUND_ENSEMBL_REST],
                matched_accession=translation, detail=detail,
            )
        return ProteinLookupResult(
            source_db=src or "Ensembl",
            protein_sequence_status=STATUS_UNAVAIL_NOT_IN_CACHE,
            lookup_method="ensembl_rest_failed",
            lookup_confidence="none", detail=detail,
        )
    if not translation:
        return ProteinLookupResult(
            source_db=src or "Ensembl",
            protein_sequence_status=STATUS_UNAVAIL_NO_ACCESSION,
            lookup_method="ensembl_missing_translation",
            lookup_confidence="none",
            detail="no_translation_id",
        )
    return ProteinLookupResult(
        source_db=src or "Ensembl",
        protein_sequence_status=STATUS_UNAVAIL_NOT_IN_CACHE,
        lookup_method="ensembl_rest_disabled",
        lookup_confidence="none",
        detail="ensembl_rest_disabled_and_not_in_cache",
    )
