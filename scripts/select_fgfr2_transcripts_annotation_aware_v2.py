#!/usr/bin/env python3
"""
select_fgfr2_transcripts_annotation_aware_v2.py

Annotation-aware transcript prioritisation for FGFR2 exon--domain boundary studies.

Designed for use after collect_fgfr2_models_dual_source_v3.py. The script selects a
cross-species reference transcript while retaining FGFR2 IIIb/IIIc isoform candidates
separately. It writes full audit, warnings, run metadata, a Markdown/HTML report and
optional per-species top-score plots.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))


class EmptyModelInput(RuntimeError):
    """No transcripts to select from, with the recorded reason attached.

    A distinct type so the pipeline can tell "the upstream step recovered nothing, and
    here is why" from a genuine defect in this step. The message is written for the
    person who started the run; the traceback stays in the log.
    """


#: Set by ``run()`` so the schema check can quote the collection status without being
#: handed the argument namespace.
_TRANSCRIPTS_PATH: Optional[Path] = None


def _explain_no_transcripts() -> str:
    from exondomaincompare.shared_gene_analysis import model_recovery as recovery

    contract = None
    if _TRANSCRIPTS_PATH is not None:
        contract = recovery.read_contract(
            _TRANSCRIPTS_PATH.parent / "collection_status.json")
    return recovery.explain_empty_input(contract, "FGFR2", "Transcript selection")

SCRIPT_NAME = "select_fgfr2_transcripts_annotation_aware_v2.py"
SCRIPT_VERSION = "2.0.0"

# ----------------------------- generic I/O -----------------------------
def read_tsv(path: Optional[Path]) -> List[Dict[str, str]]:
    if path is None:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(rows: List[Dict[str, str]], path: Path, fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def clean_id(value: str) -> str:
    value = (value or "").strip()
    return value.split(".", 1)[0] if value else ""


def as_int(value: str) -> Optional[int]:
    try:
        if value is None or str(value).strip() == "":
            return None
        return int(float(str(value).strip()))
    except Exception:
        return None


def get_any(row: Dict[str, str], names: Iterable[str], default: str = "") -> str:
    lower = {k.lower(): k for k in row.keys()}
    for name in names:
        if name in row:
            return row.get(name, default)
        lk = name.lower()
        if lk in lower:
            return row.get(lower[lk], default)
    return default


def bool_flag(value: str) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y", "t"}


def row_key(row: Dict[str, str]) -> Tuple[str, str]:
    internal = get_any(row, ["internal_transcript_id", "transcript_id_internal", "parent_transcript_id", "transcript_internal_id"])
    source = clean_id(get_any(row, ["transcript_id_source", "transcript_id", "ensembl_transcript_id", "refseq_mrna", "mrna_id"]))
    return internal, source


# ----------------------------- validation and warnings -----------------------------
def columns(rows: List[Dict[str, str]]) -> set[str]:
    out = set()
    for r in rows:
        out.update(r.keys())
    return out


def missing_columns(rows: List[Dict[str, str]], required: Sequence[str]) -> List[str]:
    existing = columns(rows)
    return [c for c in required if c not in existing]


def add_warning(warnings: List[Dict[str, str]], severity: str, code: str, message: str, species: str = "", transcript: str = "") -> None:
    warnings.append({
        "severity": severity,
        "warning_code": code,
        "species_canonical": species,
        "transcript_id": transcript,
        "message": message,
    })


def validate_input_schemas(transcripts: List[Dict[str, str]], exons: List[Dict[str, str]], strict: bool, warnings: List[Dict[str, str]]) -> None:
    required_tx = ["species_input", "species_canonical", "source_db", "gene_id_internal", "transcript_id_source", "internal_transcript_id"]
    recommended_tx = ["transcript_biotype", "translation_id_source", "protein_length_aa", "is_canonical_source", "support_level", "completeness_flags"]
    required_ex = ["exon_rank", "chrom", "start", "end", "strand"]
    if not transcripts:
        # Not this step's failure, and not this step's story to invent. An empty
        # transcript table means model collection recovered nothing, and it recorded
        # why. Raising here produced the traceback that became the Equus quagga run's
        # user-facing explanation — pointing at the wrong file and the wrong stage,
        # and never mentioning the species, the assembly, or the misspelled taxonomy
        # query four stages upstream.
        raise EmptyModelInput(_explain_no_transcripts())
    missing_tx = missing_columns(transcripts, required_tx)
    if missing_tx:
        msg = f"Missing required transcript columns: {', '.join(missing_tx)}"
        add_warning(warnings, "critical", "missing_required_transcript_columns", msg)
        if strict:
            raise ValueError(msg)
    missing_tx_recommended = missing_columns(transcripts, recommended_tx)
    if missing_tx_recommended:
        add_warning(warnings, "moderate", "missing_recommended_transcript_columns", f"Missing recommended transcript columns: {', '.join(missing_tx_recommended)}")
    if exons:
        missing_ex = missing_columns(exons, required_ex)
        id_cols_present = columns(exons).intersection({"transcript_id_internal", "internal_transcript_id", "parent_transcript_id", "transcript_id", "transcript_id_source"})
        if missing_ex:
            msg = f"Missing required exon feature columns: {', '.join(missing_ex)}"
            add_warning(warnings, "critical", "missing_required_exon_columns", msg)
            if strict:
                raise ValueError(msg)
        if not id_cols_present:
            msg = "Exon table lacks a supported transcript ID column. Supported names: transcript_id_internal, internal_transcript_id, parent_transcript_id, transcript_id, transcript_id_source."
            add_warning(warnings, "critical", "missing_exon_transcript_id_column", msg)
            if strict:
                raise ValueError(msg)
    else:
        add_warning(warnings, "critical" if strict else "moderate", "empty_exon_table", "No exon rows supplied.")
        if strict:
            raise ValueError("No exon rows supplied.")


def detect_duplicate_transcripts(transcripts: List[Dict[str, str]], warnings: List[Dict[str, str]]) -> Dict[str, str]:
    internal_counter = Counter(r.get("internal_transcript_id", "") for r in transcripts if r.get("internal_transcript_id", ""))
    source_counter = Counter(clean_id(r.get("transcript_id_source", "")) for r in transcripts if clean_id(r.get("transcript_id_source", "")))
    duplicate_map: Dict[str, str] = {}
    for r in transcripts:
        internal = r.get("internal_transcript_id", "")
        source = clean_id(r.get("transcript_id_source", ""))
        flags = []
        if internal and internal_counter[internal] > 1:
            flags.append("duplicate_internal_transcript_id")
            add_warning(warnings, "moderate", "duplicate_internal_transcript_id", f"Duplicate internal_transcript_id: {internal}", r.get("species_canonical", ""), internal)
        if source and source_counter[source] > 1:
            flags.append("duplicate_source_transcript_id")
            add_warning(warnings, "moderate", "duplicate_source_transcript_id", f"Duplicate transcript_id_source: {source}", r.get("species_canonical", ""), source)
        duplicate_map[internal or source] = ";".join(sorted(set(flags)))
    return duplicate_map


# ----------------------------- evidence normalisation -----------------------------
def normalize_biotype(value: str) -> str:
    v = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    if v in {"protein_coding", "mrna", "m_rna", "messenger_rna", "coding", "cds", "gene_with_cds", "protein_coding_cds_not_defined"}:
        return "coding"
    if "pseudogene" in v:
        return "pseudogene"
    if v in {"nmd", "nonsense_mediated_decay", "non_stop_decay", "retained_intron"} or "nonsense" in v or "retained_intron" in v:
        return "decay_or_retained_intron"
    if v in {"lncrna", "lincrna", "mirna", "snorna", "snrna", "rrna", "noncoding", "non_coding"}:
        return "noncoding"
    return "unknown" if not v else "other"


def is_refseq_select_annotation(row: Dict[str, str]) -> bool:
    val = get_any(row, ["refseq_select", "is_refseq_select", "refseq_select_status", "refseq_select_transcript"])
    return str(val or "").strip().lower() in {"1", "true", "yes", "y", "select", "refseq select", "refseq_select"}


def is_mane_select_annotation(row: Dict[str, str]) -> bool:
    val = get_any(row, ["mane_select", "is_mane_select", "mane_select_status", "mane"])
    return str(val or "").strip().lower() in {"1", "true", "yes", "y", "mane_select", "mane select", "select"}


def is_ensembl_canonical_annotation(tx: Dict[str, str], ann: Dict[str, str]) -> bool:
    val = get_any(ann, ["ensembl_canonical", "is_ensembl_canonical", "canonical_transcript"])
    return tx.get("is_canonical_source", "") == "1" or str(val or "").strip().lower() in {"1", "true", "yes", "y", "canonical", "ensembl_canonical"}


def is_appris_principal_annotation(row: Dict[str, str]) -> bool:
    text = get_any(row, ["appris", "appris_annotation", "appris_principal", "appris_principal_status"])
    t = str(text or "").strip().lower()
    return bool(re.search(r"\bprincipal(?:[1-5])?\b", t))


def classify_protein_length(protein_len: Optional[int]) -> str:
    if protein_len is None:
        return "unavailable"
    if 650 <= protein_len <= 950:
        return "plausible_fgfr2_length"
    if 450 <= protein_len < 650 or 950 < protein_len <= 1200:
        return "borderline_plausible_length"
    return "implausible_or_partial_length"


# ----------------------------- annotation/domain/exon indexes -----------------------------
def build_annotation_index(rows: List[Dict[str, str]]) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
    by_internal: Dict[str, Dict[str, str]] = {}
    by_source: Dict[str, Dict[str, str]] = {}
    for r in rows:
        internal, source = row_key(r)
        if internal:
            by_internal[internal] = r
        if source:
            by_source[source] = r
    return by_internal, by_source


def annotation_for(tx: Dict[str, str], by_internal: Dict[str, Dict[str, str]], by_source: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    internal = tx.get("internal_transcript_id", "")
    source = clean_id(tx.get("transcript_id_source", ""))
    ann: Dict[str, str] = {}
    if source and source in by_source:
        ann.update(by_source[source])
    if internal and internal in by_internal:
        ann.update(by_internal[internal])
    return ann


def build_domain_index(domain_rows: List[Dict[str, str]]) -> Tuple[Dict[str, List[Dict[str, str]]], Dict[str, List[Dict[str, str]]]]:
    by_internal: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    by_source: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for r in domain_rows:
        internal, source = row_key(r)
        if internal:
            by_internal[internal].append(r)
        if source:
            by_source[source].append(r)
    return by_internal, by_source


def domains_for(tx: Dict[str, str], by_internal: Dict[str, List[Dict[str, str]]], by_source: Dict[str, List[Dict[str, str]]]) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    internal = tx.get("internal_transcript_id", "")
    source = clean_id(tx.get("transcript_id_source", ""))
    out.extend(by_source.get(source, []))
    out.extend(by_internal.get(internal, []))
    seen = set()
    uniq = []
    for r in out:
        marker = tuple(sorted(r.items()))
        if marker not in seen:
            uniq.append(r)
            seen.add(marker)
    return uniq


def infer_domain_evidence(domain_rows: List[Dict[str, str]]) -> Tuple[bool, bool, str, str]:
    any_ig_hits: List[str] = []
    iii_hits: List[str] = []
    for d in domain_rows:
        text = " ".join(str(v) for v in d.values()).lower()
        label = get_any(d, ["domain_name", "interpro_id", "pfam_id", "domain_id"], "domain")
        if any(p in text for p in ["ig-like", "immunoglobulin", "ig domain", "i-set", "ig-like fold", "ig-like domain"]):
            any_ig_hits.append(label)
        if any(p in text for p in ["igiii", "ig iii", "domain iii", "d3", "third immunoglobulin", "immunoglobulin-like domain iii"]):
            iii_hits.append(label)
    return bool(any_ig_hits), bool(iii_hits), ";".join(sorted(set(any_ig_hits))), ";".join(sorted(set(iii_hits)))


def build_exon_index(exon_rows: List[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
    by_tx: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    id_cols = ["transcript_id_internal", "internal_transcript_id", "parent_transcript_id", "transcript_id", "transcript_id_source"]
    for ex in exon_rows:
        ids: List[str] = []
        for col in id_cols:
            raw = ex.get(col, "")
            if raw:
                ids.extend([x.strip() for x in str(raw).split(",") if x.strip()])
        for tx_id in set(ids):
            by_tx[tx_id].append(ex)
            cleaned = clean_id(tx_id)
            if cleaned and cleaned != tx_id:
                by_tx[cleaned].append(ex)
    for tx_id in by_tx:
        by_tx[tx_id].sort(key=lambda r: as_int(r.get("exon_rank", "")) or 10**9)
    return by_tx


def exons_for(tx: Dict[str, str], exon_index: Dict[str, List[Dict[str, str]]]) -> List[Dict[str, str]]:
    candidates = [tx.get("internal_transcript_id", ""), tx.get("transcript_id_source", ""), clean_id(tx.get("transcript_id_source", ""))]
    seen = set()
    out: List[Dict[str, str]] = []
    for key in candidates:
        for ex in exon_index.get(key, []):
            marker = tuple(sorted(ex.items()))
            if marker not in seen:
                out.append(ex)
                seen.add(marker)
    out.sort(key=lambda r: as_int(r.get("exon_rank", "")) or 10**9)
    return out


def exon_signature(exons: List[Dict[str, str]]) -> str:
    return "|".join(f"{e.get('chrom','')}:{e.get('start','')}-{e.get('end','')}:{e.get('strand','')}" for e in exons)


# ----------------------------- isoform classification -----------------------------
def normalize_isoform(value: str) -> str:
    v = (value or "").strip().lower().replace("-", "").replace("_", "").replace(" ", "")
    if v in {"iiib", "fgfr2iiib", "b", "kgfr"}:
        return "IIIb"
    if v in {"iiic", "fgfr2iiic", "c"}:
        return "IIIc"
    if v in {"ambiguous", "unknown", "mixed", "both"}:
        return "ambiguous"
    return "unclassified"


def build_isoform_index(rows: List[Dict[str, str]]) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    by_internal: Dict[str, str] = {}
    by_source: Dict[str, str] = {}
    evidence: Dict[str, str] = {}
    for r in rows:
        internal, source = row_key(r)
        cls = normalize_isoform(get_any(r, ["isoform_class", "iii_isoform_assignment", "fgfr2_isoform", "isoform"]))
        ev = get_any(r, ["evidence", "selection_reason", "reason"], "provided_isoform_evidence")
        if internal:
            by_internal[internal] = cls
            evidence[internal] = ev
        if source:
            by_source[source] = cls
            evidence[source] = ev
    return by_internal, by_source, evidence


def classify_isoform(tx: Dict[str, str], ann: Dict[str, str], iso_by_internal: Dict[str, str], iso_by_source: Dict[str, str], iso_ev: Dict[str, str]) -> Tuple[str, str, str]:
    internal = tx.get("internal_transcript_id", "")
    source = clean_id(tx.get("transcript_id_source", ""))
    if internal in iso_by_internal:
        return iso_by_internal[internal], "high", iso_ev.get(internal, "provided_isoform_evidence")
    if source in iso_by_source:
        return iso_by_source[source], "high", iso_ev.get(source, "provided_isoform_evidence")
    text = " ".join([tx.get("transcript_name", ""), ann.get("transcript_name", ""), ann.get("product", ""), ann.get("note", "")]).lower()
    if "iiib" in text or "isoform b" in text or "kgfr" in text:
        return "IIIb", "medium", "transcript/product name hint only"
    if "iiic" in text or "isoform c" in text:
        return "IIIc", "medium", "transcript/product name hint only"
    return "unclassified", "none", "no explicit IIIb/IIIc evidence supplied"


# ----------------------------- scoring -----------------------------
def score_transcript(tx: Dict[str, str], ann: Dict[str, str], tx_exons: List[Dict[str, str]], tx_domains: List[Dict[str, str]]) -> Tuple[int, List[str], Dict[str, str]]:
    annotation_score = 0
    structure_score = 0
    protein_score = 0
    domain_score = 0
    support_score = 0
    penalty_score = 0
    reasons: List[str] = []

    biotype_raw = tx.get("transcript_biotype", "") or ann.get("transcript_biotype", "") or ann.get("biotype", "")
    biotype_norm = normalize_biotype(biotype_raw)
    protein_len = as_int(tx.get("protein_length_aa", ""))
    has_translation = bool(tx.get("translation_id_source", "")) or (protein_len is not None and protein_len > 0)
    flags = (tx.get("completeness_flags", "") or "").lower()
    support = (tx.get("support_level", "") or ann.get("tsl", "") or ann.get("transcript_support_level", "")).lower()

    is_mane = is_mane_select_annotation(ann)
    is_refseq_select = is_refseq_select_annotation(ann)
    is_ensembl_canonical = is_ensembl_canonical_annotation(tx, ann)
    is_appris_principal = is_appris_principal_annotation(ann)
    is_ccds = bool(get_any(ann, ["ccds", "ccds_id"]))
    is_uniprot_canonical = bool_flag(get_any(ann, ["uniprot_canonical", "is_uniprot_canonical"]))
    has_any_ig, has_igiii, ig_hits, igiii_hits = infer_domain_evidence(tx_domains)

    if is_mane:
        annotation_score += 1000; reasons.append("MANE Select")
    if is_refseq_select:
        annotation_score += 900; reasons.append("RefSeq Select")
    if is_ensembl_canonical:
        annotation_score += 800; reasons.append("Ensembl Canonical")
    if is_appris_principal:
        annotation_score += 700; reasons.append("APPRIS Principal")
    if is_uniprot_canonical:
        annotation_score += 120; reasons.append("UniProt canonical concordance")
    if is_ccds:
        annotation_score += 80; reasons.append("CCDS annotated")

    if biotype_norm == "coding":
        protein_score += 200; reasons.append("coding transcript biotype")
    elif biotype_norm in {"pseudogene", "decay_or_retained_intron", "noncoding"}:
        penalty_score -= 250; reasons.append(f"penalty:biotype={biotype_norm}")
    else:
        penalty_score -= 100; reasons.append(f"non-standard biotype={biotype_raw or 'NA'}")

    if has_translation:
        protein_score += 200; reasons.append("translation/protein present")
    else:
        penalty_score -= 500; reasons.append("no translation/protein")

    n_exons = len(tx_exons)
    if n_exons > 0:
        structure_score += 100; reasons.append(f"exons present ({n_exons})")
    else:
        penalty_score -= 300; reasons.append("no exons")

    protein_length_class = classify_protein_length(protein_len)
    if protein_length_class == "plausible_fgfr2_length":
        protein_score += 120; reasons.append(f"FGFR2-plausible protein length ({protein_len} aa)")
    elif protein_length_class == "borderline_plausible_length":
        protein_score += 30; reasons.append(f"borderline plausible protein length ({protein_len} aa)")
    elif protein_length_class == "implausible_or_partial_length":
        penalty_score -= 100; reasons.append(f"implausible/partial protein length ({protein_len} aa)")
    else:
        reasons.append("protein length unavailable")

    if has_igiii:
        domain_score += 200; reasons.append("IgIII/D3-specific domain label present")
    elif has_any_ig:
        domain_score += 100; reasons.append("Ig-like domain evidence present, not D3-specific")
    elif tx_domains:
        domain_score += 25; reasons.append("domain evidence present but no Ig-like label detected")
    else:
        reasons.append("domain evidence unavailable")

    m = re.search(r"([1-5])", support)
    if m:
        tsl = int(m.group(1))
        support_score += {1: 80, 2: 50, 3: 20, 4: -20, 5: -50}.get(tsl, 0)
        reasons.append(f"TSL/support level {tsl}")

    bad_terms = ["no_translation", "no_exons", "partial", "retained_intron", "nonsense_mediated_decay", "non_stop_decay"]
    for term in bad_terms:
        if term in flags:
            penalty = -250 if term in {"no_translation", "no_exons"} else -120
            penalty_score += penalty
            reasons.append(f"penalty:{term}")

    total = annotation_score + structure_score + protein_score + domain_score + support_score + penalty_score
    features = {
        "score": str(total),
        "annotation_score": str(annotation_score),
        "structure_score": str(structure_score),
        "protein_score": str(protein_score),
        "domain_score": str(domain_score),
        "support_score": str(support_score),
        "penalty_score": str(penalty_score),
        "normalized_biotype": biotype_norm,
        "protein_length_class": protein_length_class,
        "is_mane_select": "1" if is_mane else "0",
        "is_refseq_select": "1" if is_refseq_select else "0",
        "is_ensembl_canonical": "1" if is_ensembl_canonical else "0",
        "is_appris_principal": "1" if is_appris_principal else "0",
        "is_ccds": "1" if is_ccds else "0",
        "is_uniprot_canonical": "1" if is_uniprot_canonical else "0",
        "has_translation": "1" if has_translation else "0",
        "exon_count": str(n_exons),
        "has_any_ig_like_domain_evidence": "1" if has_any_ig else "0",
        "has_igIII_or_D3_specific_domain_evidence": "1" if has_igiii else "0",
        "ig_like_domain_hits": ig_hits,
        "igIII_D3_domain_hits": igiii_hits,
    }
    return total, reasons, features


def selection_mode(row: Optional[Dict[str, str]]) -> str:
    if not row:
        return "none"
    if row.get("is_mane_select") == "1":
        return "MANE_based"
    if row.get("is_refseq_select") == "1":
        return "RefSeq_Select_based"
    if row.get("is_ensembl_canonical") == "1":
        return "Ensembl_Canonical_based"
    if row.get("is_appris_principal") == "1":
        return "APPRIS_based"
    if row.get("has_translation") == "1" and as_int(row.get("exon_count", "")) and (as_int(row.get("score", "")) or 0) > 0:
        return "rule_based_fallback"
    return "low_confidence_fallback"


def choose_best(rows: List[Dict[str, str]]) -> Optional[Dict[str, str]]:
    if not rows:
        return None
    def key(r: Dict[str, str]) -> Tuple[int, int, int, str]:
        return (
            as_int(r.get("score", "")) or -10**9,
            as_int(r.get("protein_length_aa", "")) or -1,
            as_int(r.get("exon_count", "")) or -1,
            r.get("transcript_id_source", ""),
        )
    return sorted(rows, key=key, reverse=True)[0]


# ----------------------------- reports and plots -----------------------------
def write_reports(summary_rows: List[Dict[str, str]], selected_rows: List[Dict[str, str]], warnings: List[Dict[str, str]], outdir: Path, gene_symbol: str) -> None:
    md_lines = [
        f"# Transcript selection report for {gene_symbol}",
        "",
        f"Generated by `{SCRIPT_NAME}` version `{SCRIPT_VERSION}`.",
        "",
        "## Overview",
        "",
        f"- Species processed: {len(summary_rows)}",
        f"- Selected transcript rows: {len(selected_rows)}",
        f"- Warnings: {len(warnings)}",
        "",
        "## Per-species summary",
        "",
        "| Species | Transcripts | Reference | Score | Selection mode | Isoform status | Review note |",
        "|---|---:|---|---:|---|---|---|",
    ]
    for r in summary_rows:
        md_lines.append("| {species} | {n} | {ref} | {score} | {mode} | {iso} | {note} |".format(
            species=r.get("species_canonical", ""), n=r.get("transcript_count", ""), ref=r.get("reference_transcript_id", ""),
            score=r.get("reference_score", ""), mode=r.get("selection_mode", ""), iso=r.get("isoform_status", ""), note=r.get("review_note", "")
        ))
    md_lines.extend(["", "## Warning summary", ""])
    if warnings:
        counts = Counter(w.get("warning_code", "") for w in warnings)
        for code, n in counts.most_common():
            md_lines.append(f"- `{code}`: {n}")
    else:
        md_lines.append("No warnings were generated.")
    md_lines.extend(["", "## Method note", "", "Reference transcript selection and IIIb/IIIc candidate retention are separated. Medium-confidence name hints are retained as provisional isoform candidates, whereas high-confidence candidates require explicit isoform evidence."])
    report_md = "\n".join(md_lines) + "\n"
    (outdir / "transcript_selection_report.md").write_text(report_md, encoding="utf-8")

    html_rows = []
    for r in summary_rows:
        html_rows.append("<tr>" + "".join(f"<td>{escape(r.get(k, ''))}</td>" for k in ["species_canonical", "transcript_count", "reference_transcript_id", "reference_score", "selection_mode", "isoform_status", "review_note"]) + "</tr>")
    html = f"""<!doctype html>
<html><head><meta charset='utf-8'><title>Transcript selection report</title>
<style>body{{font-family:Arial,sans-serif;max-width:1100px;margin:2rem auto;line-height:1.45}}table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #ccc;padding:0.35rem}}th{{background:#eee}}</style></head>
<body><h1>Transcript selection report for {escape(gene_symbol)}</h1><p>Generated by <code>{SCRIPT_NAME}</code> version <code>{SCRIPT_VERSION}</code>.</p>
<p><b>Species:</b> {len(summary_rows)} &nbsp; <b>Selected rows:</b> {len(selected_rows)} &nbsp; <b>Warnings:</b> {len(warnings)}</p>
<h2>Per-species summary</h2><table><thead><tr><th>Species</th><th>Transcripts</th><th>Reference</th><th>Score</th><th>Selection mode</th><th>Isoform status</th><th>Review note</th></tr></thead><tbody>{''.join(html_rows)}</tbody></table>
</body></html>"""
    (outdir / "transcript_selection_report.html").write_text(html, encoding="utf-8")


def write_score_plots(audit_rows: List[Dict[str, str]], outdir: Path, top_n: int = 10) -> List[str]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return []
    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    written: List[str] = []
    by_species: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for r in audit_rows:
        by_species[r.get("species_canonical", "unknown")].append(r)
    for species, rows in by_species.items():
        rows_sorted = sorted(rows, key=lambda r: as_int(r.get("score", "")) or -10**9, reverse=True)[:top_n]
        if not rows_sorted:
            continue
        labels = [r.get("transcript_id_source", "") or r.get("internal_transcript_id", "") for r in rows_sorted]
        scores = [as_int(r.get("score", "")) or 0 for r in rows_sorted]
        fig_width = max(7, min(14, len(labels) * 1.0))
        plt.figure(figsize=(fig_width, 4.8))
        plt.bar(range(len(scores)), scores)
        plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
        plt.ylabel("Transcript selection score")
        plt.title(f"Top FGFR2 transcript scores: {species}")
        plt.tight_layout()
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", species).strip("_") or "unknown"
        path = plot_dir / f"top_scores_{safe}.png"
        plt.savefig(path, dpi=160)
        plt.close()
        written.append(str(path))
    return written


# ----------------------------- main workflow -----------------------------
def run(args: argparse.Namespace) -> Dict[str, object]:
    warnings: List[Dict[str, str]] = []
    global _TRANSCRIPTS_PATH
    _TRANSCRIPTS_PATH = Path(args.transcripts)
    transcripts = read_tsv(args.transcripts)
    exons = read_tsv(args.exons)
    annotations = read_tsv(args.annotations)
    domains = read_tsv(args.domains)
    isoform_evidence = read_tsv(args.isoform_evidence)

    validate_input_schemas(transcripts, exons, args.strict, warnings)
    duplicate_flags = detect_duplicate_transcripts(transcripts, warnings)

    ann_by_internal, ann_by_source = build_annotation_index(annotations)
    dom_by_internal, dom_by_source = build_domain_index(domains)
    exon_index = build_exon_index(exons)
    iso_by_internal, iso_by_source, iso_ev = build_isoform_index(isoform_evidence)

    audit_rows: List[Dict[str, str]] = []
    iso_rows: List[Dict[str, str]] = []
    selected_rows: List[Dict[str, str]] = []
    summary_rows: List[Dict[str, str]] = []

    for tx in transcripts:
        ann = annotation_for(tx, ann_by_internal, ann_by_source)
        tx_exons = exons_for(tx, exon_index)
        tx_domains = domains_for(tx, dom_by_internal, dom_by_source)
        _score, reasons, feat = score_transcript(tx, ann, tx_exons, tx_domains)
        iso_cls, iso_conf, iso_reason = classify_isoform(tx, ann, iso_by_internal, iso_by_source, iso_ev)
        ident = tx.get("internal_transcript_id", "") or clean_id(tx.get("transcript_id_source", ""))
        if not tx_exons:
            add_warning(warnings, "moderate", "transcript_without_exon_links", "No exon rows linked to transcript.", tx.get("species_canonical", ""), ident)
        if feat.get("has_translation") != "1":
            add_warning(warnings, "moderate", "transcript_without_translation", "Transcript has no detectable translation/protein length.", tx.get("species_canonical", ""), ident)
        base = {
            "species_input": tx.get("species_input", ""),
            "species_canonical": tx.get("species_canonical", ""),
            "source_db": tx.get("source_db", ""),
            "gene_id_internal": tx.get("gene_id_internal", ""),
            "transcript_id_source": tx.get("transcript_id_source", ""),
            "internal_transcript_id": tx.get("internal_transcript_id", ""),
            "transcript_name": tx.get("transcript_name", "") or get_any(ann, ["transcript_name"]),
            "transcript_biotype": tx.get("transcript_biotype", ""),
            "translation_id_source": tx.get("translation_id_source", ""),
            "protein_length_aa": tx.get("protein_length_aa", ""),
            "support_level": tx.get("support_level", "") or get_any(ann, ["tsl", "transcript_support_level"]),
            "completeness_flags": tx.get("completeness_flags", ""),
            **feat,
            "duplicate_transcript_flags": duplicate_flags.get(ident, ""),
            "iii_isoform_assignment": iso_cls,
            "iii_isoform_confidence": iso_conf,
            "iii_isoform_evidence": iso_reason,
            "exon_signature": exon_signature(tx_exons),
            "selection_reason": "; ".join(reasons),
        }
        audit_rows.append(base)
        iso_rows.append({k: base.get(k, "") for k in [
            "species_input", "species_canonical", "source_db", "transcript_id_source", "internal_transcript_id",
            "iii_isoform_assignment", "iii_isoform_confidence", "iii_isoform_evidence", "exon_count",
            "protein_length_aa", "has_any_ig_like_domain_evidence", "has_igIII_or_D3_specific_domain_evidence", "exon_signature"
        ]})

    by_species: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for r in audit_rows:
        by_species[(r.get("species_input", ""), r.get("species_canonical", ""))].append(r)

    for (species_input, species_canonical), rows in sorted(by_species.items()):
        eligible_ref = [r for r in rows if r.get("has_translation") == "1" and (as_int(r.get("score", "")) or 0) >= args.min_reference_score]
        ref = choose_best(eligible_ref) or choose_best(rows)
        high_iiib = choose_best([r for r in rows if r.get("iii_isoform_assignment") == "IIIb" and r.get("iii_isoform_confidence") == "high" and r.get("has_translation") == "1"])
        high_iiic = choose_best([r for r in rows if r.get("iii_isoform_assignment") == "IIIc" and r.get("iii_isoform_confidence") == "high" and r.get("has_translation") == "1"])
        med_iiib = choose_best([r for r in rows if r.get("iii_isoform_assignment") == "IIIb" and r.get("iii_isoform_confidence") == "medium" and r.get("has_translation") == "1"])
        med_iiic = choose_best([r for r in rows if r.get("iii_isoform_assignment") == "IIIc" and r.get("iii_isoform_confidence") == "medium" and r.get("has_translation") == "1"])

        selected_by_id: Dict[str, List[str]] = defaultdict(list)
        staged: List[Tuple[str, Optional[Dict[str, str]], str, str]] = [
            ("reference", ref, "1", "highest-ranked representative transcript for cross-species comparison"),
            ("FGFR2_IIIb_candidate", high_iiib, "1", "best transcript with high-confidence IIIb evidence"),
            ("FGFR2_IIIc_candidate", high_iiic, "1", "best transcript with high-confidence IIIc evidence"),
            ("FGFR2_IIIb_provisional", med_iiib if not high_iiib else None, "provisional", "best transcript with medium-confidence IIIb name evidence only"),
            ("FGFR2_IIIc_provisional", med_iiic if not high_iiic else None, "provisional", "best transcript with medium-confidence IIIc name evidence only"),
        ]
        for role, r, rank, _reason in staged:
            if not r:
                continue
            selected_by_id[r.get("internal_transcript_id", "") or r.get("transcript_id_source", "")].append(role)
        used_ids = set(selected_by_id.keys())
        alt_count = 0
        for r in sorted(rows, key=lambda x: as_int(x.get("score", "")) or -10**9, reverse=True):
            ident = r.get("internal_transcript_id", "") or r.get("transcript_id_source", "")
            if ident in used_ids:
                continue
            if r.get("has_translation") == "1" and (as_int(r.get("score", "")) or 0) >= args.min_reference_score:
                if alt_count < args.max_alternatives_per_species:
                    staged.append(("alternative_complete", r, "secondary", "complete high-scoring alternative retained for sensitivity analysis"))
                    selected_by_id[ident].append("alternative_complete")
                    alt_count += 1

        for role, r, rank, extra_reason in staged:
            if not r:
                continue
            ident = r.get("internal_transcript_id", "") or r.get("transcript_id_source", "")
            roles = sorted(set(selected_by_id.get(ident, [role])))
            overlapping = "1" if len(roles) > 1 else "0"
            if overlapping == "1":
                add_warning(warnings, "info", "overlapping_selection_roles", f"Transcript selected for multiple roles: {', '.join(roles)}", species_canonical, ident)
            selected_rows.append({
                **r,
                "gene_symbol": args.gene_symbol,
                "selection_role": role,
                "selection_rank": rank,
                "selection_mode": selection_mode(ref) if role == "reference" else selection_mode(r),
                "overlapping_roles": overlapping,
                "all_roles_for_transcript": ";".join(roles),
                "final_selection_reason": extra_reason + "; " + r.get("selection_reason", ""),
            })

        if not ref:
            add_warning(warnings, "critical", "no_reference_transcript", "No reference transcript could be selected.", species_canonical)
        elif (as_int(ref.get("score", "")) or 0) < args.min_reference_score:
            add_warning(warnings, "moderate", "reference_below_threshold", "Reference transcript is below min_reference_score and requires review.", species_canonical, ref.get("internal_transcript_id", ""))
        if not high_iiib and med_iiib:
            add_warning(warnings, "info", "only_medium_IIIb_evidence", "Only medium-confidence IIIb evidence available.", species_canonical, med_iiib.get("internal_transcript_id", ""))
        if not high_iiic and med_iiic:
            add_warning(warnings, "info", "only_medium_IIIc_evidence", "Only medium-confidence IIIc evidence available.", species_canonical, med_iiic.get("internal_transcript_id", ""))
        if not high_iiib and not high_iiic:
            add_warning(warnings, "moderate", "no_high_confidence_IIIb_IIIc_pair", "No high-confidence IIIb/IIIc pair detected. Provide --isoform_evidence for final analysis.", species_canonical)

        summary_rows.append({
            "species_input": species_input,
            "species_canonical": species_canonical,
            "transcript_count": str(len(rows)),
            "translated_transcript_count": str(sum(1 for r in rows if r.get("has_translation") == "1")),
            "reference_transcript_id": ref.get("transcript_id_source", "") if ref else "",
            "reference_score": ref.get("score", "") if ref else "",
            "reference_quality": "usable" if ref and (as_int(ref.get("score", "")) or 0) >= args.min_reference_score else "low_confidence_review_required",
            "selection_mode": selection_mode(ref),
            "IIIb_candidate_transcript_id": high_iiib.get("transcript_id_source", "") if high_iiib else "",
            "IIIc_candidate_transcript_id": high_iiic.get("transcript_id_source", "") if high_iiic else "",
            "IIIb_provisional_transcript_id": med_iiib.get("transcript_id_source", "") if med_iiib and not high_iiib else "",
            "IIIc_provisional_transcript_id": med_iiic.get("transcript_id_source", "") if med_iiic and not high_iiic else "",
            "isoform_status": (
                "both_high_confidence_IIIb_and_IIIc_detected" if high_iiib and high_iiic else
                "only_high_confidence_IIIb_detected" if high_iiib else
                "only_high_confidence_IIIc_detected" if high_iiic else
                "provisional_name_hints_only" if med_iiib or med_iiic else
                "no_explicit_IIIb_IIIc_evidence"
            ),
            "review_note": "Provide --isoform_evidence from exon/sequence/domain mapping for high-confidence IIIb/IIIc calls." if not (high_iiib and high_iiic) else "",
        })

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    common_fields = [
        "species_input", "species_canonical", "source_db", "gene_symbol", "gene_id_internal",
        "selection_role", "selection_rank", "selection_mode", "overlapping_roles", "all_roles_for_transcript",
        "transcript_id_source", "internal_transcript_id", "transcript_name", "transcript_biotype", "normalized_biotype",
        "translation_id_source", "protein_length_aa", "protein_length_class", "score", "annotation_score", "structure_score",
        "protein_score", "domain_score", "support_score", "penalty_score", "is_mane_select", "is_refseq_select",
        "is_ensembl_canonical", "is_appris_principal", "is_ccds", "is_uniprot_canonical", "has_translation",
        "exon_count", "support_level", "completeness_flags", "has_any_ig_like_domain_evidence",
        "has_igIII_or_D3_specific_domain_evidence", "ig_like_domain_hits", "igIII_D3_domain_hits", "duplicate_transcript_flags",
        "iii_isoform_assignment", "iii_isoform_confidence", "iii_isoform_evidence", "final_selection_reason", "exon_signature",
    ]
    write_tsv(selected_rows, outdir / "selected_transcripts.tsv", common_fields)
    audit_fields = [f for f in common_fields if f not in {"gene_symbol", "selection_role", "selection_rank", "selection_mode", "overlapping_roles", "all_roles_for_transcript", "final_selection_reason"}]
    audit_fields.append("selection_reason")
    write_tsv(audit_rows, outdir / "transcript_selection_audit.tsv", audit_fields)
    write_tsv(iso_rows, outdir / "isoform_classification.tsv", [
        "species_input", "species_canonical", "source_db", "transcript_id_source", "internal_transcript_id",
        "iii_isoform_assignment", "iii_isoform_confidence", "iii_isoform_evidence", "exon_count",
        "protein_length_aa", "has_any_ig_like_domain_evidence", "has_igIII_or_D3_specific_domain_evidence", "exon_signature"
    ])
    write_tsv(summary_rows, outdir / "transcript_selection_summary.tsv", [
        "species_input", "species_canonical", "transcript_count", "translated_transcript_count", "reference_transcript_id",
        "reference_score", "reference_quality", "selection_mode", "IIIb_candidate_transcript_id", "IIIc_candidate_transcript_id",
        "IIIb_provisional_transcript_id", "IIIc_provisional_transcript_id", "isoform_status", "review_note"
    ])
    write_tsv(warnings, outdir / "transcript_selection_warnings.tsv", ["severity", "warning_code", "species_canonical", "transcript_id", "message"])

    methods = f"""Transcript selection method for {args.gene_symbol}\n\nTranscript selection was performed with an annotation-aware hierarchical strategy. For each species, one primary reference transcript was selected for cross-species comparability. The ranking separated annotation, structure, protein, domain, support and penalty sub-scores. It prioritized MANE Select, RefSeq Select, Ensembl Canonical and APPRIS Principal annotations, followed by rule-based evidence including coding biotype, translation availability, exon structure, plausible FGFR2 protein length, transcript support level and Ig-like domain evidence. Because FGFR2 IIIb/IIIc isoforms are produced by mutually exclusive exon usage in the IgIII/D3 ligand-binding region, isoform-specific IIIb and IIIc candidates were retained separately from the reference transcript. High-confidence IIIb/IIIc calls require explicit isoform evidence; name-based hints are retained only as provisional candidates. All scores, penalties, warnings and final decisions were written to audit tables for reproducibility and manual review.\n"""
    (outdir / "transcript_selection_methods.txt").write_text(methods, encoding="utf-8")
    write_reports(summary_rows, selected_rows, warnings, outdir, args.gene_symbol)
    plot_paths: List[str] = []
    if not args.no_plots:
        plot_paths = write_score_plots(audit_rows, outdir, top_n=args.plot_top_n)

    metadata = {
        "script_name": SCRIPT_NAME,
        "script_version": SCRIPT_VERSION,
        "run_datetime_utc": datetime.now(timezone.utc).isoformat(),
        "parameters": {
            "gene_symbol": args.gene_symbol,
            "min_reference_score": args.min_reference_score,
            "max_alternatives_per_species": args.max_alternatives_per_species,
            "strict": args.strict,
            "plot_top_n": args.plot_top_n,
            "no_plots": args.no_plots,
        },
        "input_files": {k: str(getattr(args, k)) if getattr(args, k) else "" for k in ["transcripts", "exons", "annotations", "domains", "isoform_evidence"]},
        "input_row_counts": {"transcripts": len(transcripts), "exons": len(exons), "annotations": len(annotations), "domains": len(domains), "isoform_evidence": len(isoform_evidence)},
        "output_row_counts": {"audit": len(audit_rows), "selected": len(selected_rows), "summary": len(summary_rows), "warnings": len(warnings)},
        "plot_files": plot_paths,
    }
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    if args.strict and any(w.get("severity") == "critical" for w in warnings):
        raise ValueError("Critical warnings were generated in strict mode. See transcript_selection_warnings.tsv.")
    return metadata


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Annotation-aware FGFR2 transcript selection with audit trail, warnings, reports and plots.")
    parser.add_argument("--transcripts", required=True, type=Path)
    parser.add_argument("--exons", required=True, type=Path)
    parser.add_argument("--annotations", type=Path, default=None)
    parser.add_argument("--domains", type=Path, default=None)
    parser.add_argument("--isoform_evidence", type=Path, default=None)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min_reference_score", type=int, default=350)
    parser.add_argument("--max_alternatives_per_species", type=int, default=3)
    parser.add_argument("--gene_symbol", default="FGFR2")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--no_plots", action="store_true")
    parser.add_argument("--plot_top_n", type=int, default=10)
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    try:
        metadata = run(args)
    except EmptyModelInput as exc:
        # A precise, user-facing reason on stderr and a non-zero exit; no traceback,
        # because the traceback describes this step and the cause is upstream.
        print(f"[BLOCKED] {exc}", file=sys.stderr)
        return 2
    print(f"[OK] transcripts processed: {metadata['input_row_counts']['transcripts']}")
    print(f"[OK] species processed: {metadata['output_row_counts']['summary']}")
    print(f"[OK] selected rows: {metadata['output_row_counts']['selected']}")
    print(f"[OK] warnings: {metadata['output_row_counts']['warnings']}")
    print(f"[OK] output directory: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
