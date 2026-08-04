#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

SCRIPT_NAME = "validate_fgfr2_paralog_screen_and_paralogy.py"
SCRIPT_VERSION = "2.0-final"
FGFR_GENES = ("FGFR1", "FGFR2", "FGFR3", "FGFR4")


@dataclass
class FastaRecord:
    id: str
    header: str
    seq: str


@dataclass
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int
    ref_gene: str = "UNKNOWN"
    ref_species: str = "unknown"

    @property
    def qcov(self) -> float:
        return self.length / self.qlen if self.qlen else 0.0

    @property
    def scov(self) -> float:
        return self.length / self.slen if self.slen else 0.0


def read_fasta(path: Path) -> List[FastaRecord]:
    records: List[FastaRecord] = []
    cur_header: Optional[str] = None
    cur_seq: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_header is not None:
                    rid = cur_header.split()[0]
                    records.append(FastaRecord(rid, cur_header, "".join(cur_seq).upper().replace(" ", "")))
                cur_header = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(re.sub(r"\s+", "", line))
    if cur_header is not None:
        rid = cur_header.split()[0]
        records.append(FastaRecord(rid, cur_header, "".join(cur_seq).upper().replace(" ", "")))
    return records


def write_tsv(path: Path, rows: List[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(fieldnames), extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def read_tsv_optional(path: Optional[Path]) -> Dict[str, dict]:
    if not path or not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    out: Dict[str, dict] = {}
    for r in rows:
        for key in ("output_id", "fasta_id", "query_id", "protein_fasta_id"):
            if r.get(key):
                out[r[key]] = r
                break
    return out


def extract_header_field(text: str, key: str) -> str:
    for part in str(text).split("|"):
        if part.startswith(key + "="):
            return part.split("=", 1)[1]
    m = re.search(rf"(?:^|[\s|;]){re.escape(key)}=([^\s|;]+)", str(text))
    return m.group(1) if m else ""


def extract_species_from_query_id(query_id: str) -> str:
    return extract_header_field(query_id, "species") or "unknown"


def extract_role_from_query_id(query_id: str) -> str:
    return extract_header_field(query_id, "role") or ""


def extract_isoform_from_query_id(query_id: str) -> str:
    return extract_header_field(query_id, "isoform") or ""


def infer_ref_gene(record: FastaRecord) -> str:
    text = f"{record.id} {record.header}".upper()
    # prefer exact FGFR1-4 tokens
    for gene in FGFR_GENES:
        if re.search(rf"(^|[^A-Z0-9]){gene}([^A-Z0-9]|$)", text):
            return gene
    # fallback for receptor names with spaces
    for n, gene in [("1", "FGFR1"), ("2", "FGFR2"), ("3", "FGFR3"), ("4", "FGFR4")]:
        if re.search(rf"FIBROBLAST GROWTH FACTOR RECEPTOR\s*{n}", text):
            return gene
    return "UNKNOWN"


def infer_ref_species(record: FastaRecord) -> str:
    header = record.header
    for key in ("species", "organism", "taxon"):
        val = extract_header_field(header, key)
        if val:
            return val
    # UniProt OS=Homo sapiens OX=9606 style
    m = re.search(r"\bOS=([^=]+?)(?:\sOX=|\sGN=|\sPE=|\sSV=|$)", header)
    if m:
        return re.sub(r"\s+", "_", m.group(1).strip().lower())
    # common simplified names
    text = header.lower()
    known = [
        ("homo_sapiens", "homo_sapiens"), ("human", "homo_sapiens"),
        ("mus_musculus", "mus_musculus"), ("mouse", "mus_musculus"),
        ("gallus_gallus", "gallus_gallus"), ("chicken", "gallus_gallus"),
        ("danio_rerio", "danio_rerio"), ("zebrafish", "danio_rerio"),
        ("xenopus", "xenopus"),
    ]
    for token, species in known:
        if token in text:
            return species
    return "unknown"


def check_reference_fasta(refs: List[FastaRecord]) -> Tuple[Dict[str, str], Dict[str, str], List[dict]]:
    warnings: List[dict] = []
    sseq_to_gene: Dict[str, str] = {}
    sseq_to_species: Dict[str, str] = {}
    genes_seen: Counter[str] = Counter()
    for r in refs:
        gene = infer_ref_gene(r)
        species = infer_ref_species(r)
        sseq_to_gene[r.id] = gene
        sseq_to_species[r.id] = species
        genes_seen[gene] += 1
        if gene == "UNKNOWN":
            warnings.append({
                "warning_code": "reference_gene_unknown",
                "severity": "error",
                "message": f"Could not infer FGFR gene name from reference FASTA header: {r.header}",
                "query_id": "",
                "subject_id": r.id,
            })
    for gene in FGFR_GENES:
        if genes_seen[gene] == 0:
            warnings.append({
                "warning_code": "reference_gene_missing",
                "severity": "error",
                "message": f"Reference FASTA does not contain a detectable {gene} sequence.",
                "query_id": "",
                "subject_id": "",
            })
    return sseq_to_gene, sseq_to_species, warnings


def run_blastp(query_fasta: Path, db_fasta: Path, out_tsv: Path, outdir: Path, db_name: str, threads: int, evalue: str, max_targets: int = 1000) -> None:
    blastp = shutil.which("blastp")
    makeblastdb = shutil.which("makeblastdb")
    if not blastp or not makeblastdb:
        raise RuntimeError("BLAST+ not found: blastp and makeblastdb must be available in PATH, or use --fallback_pairwise.")
    db_prefix = outdir / db_name
    subprocess.run([makeblastdb, "-in", str(db_fasta), "-dbtype", "prot", "-out", str(db_prefix)], check=True)
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    cmd = [
        blastp, "-query", str(query_fasta), "-db", str(db_prefix), "-out", str(out_tsv),
        "-outfmt", outfmt, "-evalue", str(evalue), "-num_threads", str(threads),
        "-max_target_seqs", str(max_targets), "-seg", "yes",
    ]
    subprocess.run(cmd, check=True)


def simple_similarity_score(a: str, b: str) -> Tuple[float, float, float, int]:
    a = a.replace("*", "")
    b = b.replace("*", "")
    if not a or not b:
        return 0.0, 0.0, 0.0, 0
    if len(a) <= len(b):
        short, long = a, b
    else:
        short, long = b, a
    step = max(1, len(short) // 250)
    best = 0
    best_len = len(short)
    for i in range(0, max(1, len(long) - len(short) + 1), step):
        window = long[i:i + len(short)]
        matches = sum(1 for x, y in zip(short, window) if x == y)
        best = max(best, matches)
    pident = 100.0 * best / best_len if best_len else 0.0
    qcov = min(len(a), len(b)) / len(a) if a else 0.0
    scov = min(len(a), len(b)) / len(b) if b else 0.0
    return pident, qcov, scov, int(best)


def run_fallback_pairwise(queries: List[FastaRecord], refs: List[FastaRecord], sseq_to_gene: Dict[str, str], sseq_to_species: Dict[str, str], out_tsv: Path) -> None:
    rows: List[dict] = []
    for q in queries:
        for s in refs:
            pident, qcov, _scov, score = simple_similarity_score(q.seq, s.seq)
            align_len = int(min(len(q.seq), len(s.seq)) * qcov)
            rows.append({
                "qseqid": q.id, "sseqid": s.id, "pident": f"{pident:.3f}",
                "length": str(align_len), "mismatch": "NA", "gapopen": "NA",
                "qstart": "1", "qend": str(align_len), "sstart": "1", "send": str(align_len),
                "evalue": "NA", "bitscore": f"{score:.3f}", "qlen": str(len(q.seq)),
                "slen": str(len(s.seq)), "ref_gene": sseq_to_gene.get(s.id, "UNKNOWN"),
                "ref_species": sseq_to_species.get(s.id, "unknown"),
            })
    fields = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "ref_gene", "ref_species"]
    write_tsv(out_tsv, rows, fields)


def parse_blast_tsv(path: Path, sseq_to_gene: Dict[str, str], sseq_to_species: Dict[str, str]) -> List[BlastHit]:
    hits: List[BlastHit] = []
    if not path.exists():
        return hits
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            c = line.split("\t")
            if len(c) < 14:
                continue
            try:
                sseqid = c[1]
                hits.append(BlastHit(
                    qseqid=c[0], sseqid=sseqid, pident=float(c[2]), length=int(float(c[3])),
                    mismatch=int(float(c[4])) if c[4] != "NA" else 0,
                    gapopen=int(float(c[5])) if c[5] != "NA" else 0,
                    qstart=int(float(c[6])), qend=int(float(c[7])),
                    sstart=int(float(c[8])), send=int(float(c[9])),
                    evalue=float(c[10]) if c[10] != "NA" else math.nan,
                    bitscore=float(c[11]), qlen=int(float(c[12])), slen=int(float(c[13])),
                    ref_gene=sseq_to_gene.get(sseqid, c[14] if len(c) > 14 else "UNKNOWN"),
                    ref_species=sseq_to_species.get(sseqid, c[15] if len(c) > 15 else "unknown"),
                ))
            except Exception:
                continue
    return hits


def best_hit_per_gene(hits: List[BlastHit]) -> Dict[str, BlastHit]:
    best: Dict[str, BlastHit] = {}
    for h in hits:
        if h.ref_gene not in FGFR_GENES:
            continue
        old = best.get(h.ref_gene)
        if old is None or (h.bitscore, h.pident, h.qcov) > (old.bitscore, old.pident, old.qcov):
            best[h.ref_gene] = h
    return best


def classify(best_gene: str, pident: float, qcov: float, margin_fraction: float, args: argparse.Namespace) -> str:
    if not best_gene or best_gene == "NO_HIT":
        return "no_detectable_fgfr_similarity"
    if qcov < args.min_qcov_for_any_call or pident < args.min_pident_for_any_call:
        return "insufficient_sequence_similarity"
    if best_gene != "FGFR2":
        if margin_fraction >= args.paralog_margin_fraction:
            return "possible_paralog_misassignment"
        return "ambiguous_fgfr_family_member"
    if margin_fraction >= args.high_conf_margin_fraction and qcov >= args.high_conf_qcov and pident >= args.high_conf_pident:
        return "high_confidence_FGFR2"
    if margin_fraction >= args.probable_margin_fraction and qcov >= args.probable_qcov and pident >= args.probable_pident:
        return "probable_FGFR2"
    return "FGFR2_best_hit_but_low_margin"


def recommendation_from_status(status: str) -> str:
    if status == "high_confidence_FGFR2":
        return "include_primary"
    if status == "probable_FGFR2":
        return "include_primary_with_note"
    if status == "FGFR2_best_hit_but_low_margin":
        return "keep_as_uncertain_case"
    if status in {"ambiguous_fgfr_family_member", "possible_paralog_misassignment"}:
        return "exclude_or_manual_review"
    return "exclude_from_primary_analysis"


def make_validation_rows(queries: List[FastaRecord], hits: List[BlastHit], protein_report: Dict[str, dict], args: argparse.Namespace) -> Tuple[List[dict], List[dict]]:
    by_query: Dict[str, List[BlastHit]] = defaultdict(list)
    for h in hits:
        by_query[h.qseqid].append(h)
    rows: List[dict] = []
    risks: List[dict] = []
    for q in queries:
        qhits = by_query.get(q.id, [])
        per_gene = best_hit_per_gene(qhits)
        ranked = sorted(per_gene.values(), key=lambda h: (h.bitscore, h.pident, h.qcov), reverse=True)
        if ranked:
            best = ranked[0]
            second = ranked[1] if len(ranked) > 1 else None
            best_gene = best.ref_gene
            best_score = best.bitscore
            second_gene = second.ref_gene if second else "NA"
            second_score = second.bitscore if second else 0.0
            margin_abs = best_score - second_score
            margin_fraction = margin_abs / best_score if best_score > 0 else 0.0
            status = classify(best_gene, best.pident, best.qcov, margin_fraction, args)
            best_pident, best_qcov, best_scov = best.pident, best.qcov, best.scov
            best_evalue, best_subject = best.evalue, best.sseqid
            best_ref_species = best.ref_species
        else:
            best_gene, best_subject, best_ref_species, second_gene = "NO_HIT", "", "", "NA"
            best_score = second_score = margin_abs = margin_fraction = 0.0
            best_pident = best_qcov = best_scov = 0.0
            best_evalue = math.nan
            status = "no_detectable_fgfr_similarity"

        gene_scores = {g: per_gene[g].bitscore for g in per_gene}
        gene_subjects = {g: per_gene[g].sseqid for g in per_gene}
        gene_species = {g: per_gene[g].ref_species for g in per_gene}
        gene_pident = {g: per_gene[g].pident for g in per_gene}
        gene_qcov = {g: per_gene[g].qcov for g in per_gene}
        n_refs_hit = Counter(h.ref_gene for h in qhits if h.ref_gene in FGFR_GENES)
        n_ref_species_hit = {g: len({h.ref_species for h in qhits if h.ref_gene == g and h.bitscore > 0}) for g in FGFR_GENES}

        meta = protein_report.get(q.id, {})
        species = meta.get("species_canonical") or meta.get("species") or extract_species_from_query_id(q.id)
        role = meta.get("selection_role") or meta.get("role") or extract_role_from_query_id(q.id)
        isoform = meta.get("isoform") or extract_isoform_from_query_id(q.id)
        source_db = meta.get("source_db") or extract_header_field(q.id, "source")
        transcript_id = meta.get("transcript_id_source") or meta.get("transcript_id") or extract_header_field(q.id, "transcript")
        protein_id = meta.get("translation_id_source") or meta.get("protein_id") or extract_header_field(q.id, "protein")

        row = {
            "query_id": q.id,
            "species_canonical": species,
            "species_input": meta.get("species_input", ""),
            "source_db": source_db,
            "selection_role": role,
            "isoform": isoform,
            "transcript_id_source": transcript_id,
            "translation_id_source": protein_id,
            "protein_length_aa": str(len(q.seq.replace("*", ""))),
            "best_ref_gene": best_gene,
            "best_ref_subject": best_subject,
            "best_ref_species": best_ref_species,
            "best_bitscore": f"{best_score:.3f}",
            "second_ref_gene": second_gene,
            "second_bitscore": f"{second_score:.3f}",
            "fgfr1_bitscore": f"{gene_scores.get('FGFR1', 0.0):.3f}",
            "fgfr2_bitscore": f"{gene_scores.get('FGFR2', 0.0):.3f}",
            "fgfr3_bitscore": f"{gene_scores.get('FGFR3', 0.0):.3f}",
            "fgfr4_bitscore": f"{gene_scores.get('FGFR4', 0.0):.3f}",
            "fgfr1_best_subject": gene_subjects.get("FGFR1", ""),
            "fgfr2_best_subject": gene_subjects.get("FGFR2", ""),
            "fgfr3_best_subject": gene_subjects.get("FGFR3", ""),
            "fgfr4_best_subject": gene_subjects.get("FGFR4", ""),
            "fgfr1_best_ref_species": gene_species.get("FGFR1", ""),
            "fgfr2_best_ref_species": gene_species.get("FGFR2", ""),
            "fgfr3_best_ref_species": gene_species.get("FGFR3", ""),
            "fgfr4_best_ref_species": gene_species.get("FGFR4", ""),
            "fgfr1_best_pident": f"{gene_pident.get('FGFR1', 0.0):.3f}",
            "fgfr2_best_pident": f"{gene_pident.get('FGFR2', 0.0):.3f}",
            "fgfr3_best_pident": f"{gene_pident.get('FGFR3', 0.0):.3f}",
            "fgfr4_best_pident": f"{gene_pident.get('FGFR4', 0.0):.3f}",
            "fgfr1_best_qcov": f"{gene_qcov.get('FGFR1', 0.0):.5f}",
            "fgfr2_best_qcov": f"{gene_qcov.get('FGFR2', 0.0):.5f}",
            "fgfr3_best_qcov": f"{gene_qcov.get('FGFR3', 0.0):.5f}",
            "fgfr4_best_qcov": f"{gene_qcov.get('FGFR4', 0.0):.5f}",
            "fgfr1_n_reference_hits": str(n_refs_hit.get("FGFR1", 0)),
            "fgfr2_n_reference_hits": str(n_refs_hit.get("FGFR2", 0)),
            "fgfr3_n_reference_hits": str(n_refs_hit.get("FGFR3", 0)),
            "fgfr4_n_reference_hits": str(n_refs_hit.get("FGFR4", 0)),
            "fgfr1_n_reference_species_hit": str(n_ref_species_hit.get("FGFR1", 0)),
            "fgfr2_n_reference_species_hit": str(n_ref_species_hit.get("FGFR2", 0)),
            "fgfr3_n_reference_species_hit": str(n_ref_species_hit.get("FGFR3", 0)),
            "fgfr4_n_reference_species_hit": str(n_ref_species_hit.get("FGFR4", 0)),
            "score_margin_abs": f"{margin_abs:.3f}",
            "score_margin_fraction": f"{margin_fraction:.5f}",
            "best_pident": f"{best_pident:.3f}",
            "best_query_coverage": f"{best_qcov:.5f}",
            "best_subject_coverage": f"{best_scov:.5f}",
            "best_evalue": "NA" if math.isnan(best_evalue) else f"{best_evalue:.3g}",
            "fgfr2_screen_status": status,
            "primary_analysis_recommendation": recommendation_from_status(status),
        }
        rows.append(row)
        if status not in {"high_confidence_FGFR2", "probable_FGFR2"}:
            risks.append(row)
    return rows, risks


def safe_float(x) -> float:
    try:
        if x in (None, "", "NA"):
            return math.nan
        return float(x)
    except Exception:
        return math.nan


def load_rows(path: Path) -> List[dict]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def species_summary(validation_rows: List[dict]) -> List[dict]:
    by_species: Dict[str, List[dict]] = defaultdict(list)
    for r in validation_rows:
        sp = r.get("species_canonical") or extract_species_from_query_id(r.get("query_id", "")) or "unknown"
        by_species[sp].append(r)
    out: List[dict] = []
    for sp, rows in sorted(by_species.items()):
        statuses = sorted({r.get("fgfr2_screen_status", "") for r in rows})
        min_margin = min([safe_float(r.get("score_margin_fraction")) for r in rows if not math.isnan(safe_float(r.get("score_margin_fraction")))] or [math.nan])
        min_qcov = min([safe_float(r.get("best_query_coverage")) for r in rows if not math.isnan(safe_float(r.get("best_query_coverage")))] or [math.nan])
        best_genes = Counter(r.get("best_ref_gene", "") for r in rows)
        next_genes = Counter(r.get("second_ref_gene", "") for r in rows)
        n_high = sum(1 for r in rows if r.get("fgfr2_screen_status") == "high_confidence_FGFR2")
        n_prob = sum(1 for r in rows if r.get("fgfr2_screen_status") == "probable_FGFR2")
        if n_high == len(rows):
            sp_status = "all_high_confidence_FGFR2"
        elif n_high + n_prob == len(rows):
            sp_status = "all_supported_FGFR2"
        elif any(r.get("best_ref_gene") != "FGFR2" for r in rows):
            sp_status = "paralog_or_ambiguous_case_present"
        else:
            sp_status = "mixed_or_uncertain"
        out.append({
            "species": sp,
            "n_proteins": str(len(rows)),
            "n_high_confidence_FGFR2": str(n_high),
            "n_probable_FGFR2": str(n_prob),
            "min_score_margin_fraction": "NA" if math.isnan(min_margin) else f"{min_margin:.5f}",
            "min_query_coverage": "NA" if math.isnan(min_qcov) else f"{min_qcov:.5f}",
            "best_ref_gene_counts": ";".join(f"{k}:{v}" for k, v in best_genes.most_common()),
            "second_ref_gene_counts": ";".join(f"{k}:{v}" for k, v in next_genes.most_common()),
            "statuses": ";".join(statuses),
            "species_fgfr2_screen_status": sp_status,
        })
    return out


def write_annotated_blast(hits: List[BlastHit], path: Path) -> None:
    rows = []
    for h in hits:
        rows.append({
            "qseqid": h.qseqid, "sseqid": h.sseqid, "ref_gene": h.ref_gene, "ref_species": h.ref_species,
            "pident": f"{h.pident:.3f}", "length": h.length, "mismatch": h.mismatch,
            "gapopen": h.gapopen, "qstart": h.qstart, "qend": h.qend, "sstart": h.sstart, "send": h.send,
            "evalue": f"{h.evalue:.3g}" if not math.isnan(h.evalue) else "NA", "bitscore": f"{h.bitscore:.3f}",
            "qlen": h.qlen, "slen": h.slen, "qcov": f"{h.qcov:.5f}", "scov": f"{h.scov:.5f}",
        })
    fields = ["qseqid", "sseqid", "ref_gene", "ref_species", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcov", "scov"]
    write_tsv(path, rows, fields)


def reciprocal_panel_test(refs: List[FastaRecord], queries: List[FastaRecord], reference_fasta: Path, query_fasta: Path, outdir: Path, prefix: str, threads: int, evalue: str, fallback: bool, sseq_to_gene: Dict[str, str], sseq_to_species: Dict[str, str]) -> Tuple[Path, Path, List[dict]]:
    raw = outdir / f"{prefix}_reciprocal_panel_blastp.tsv"
    annotated = outdir / f"{prefix}_reciprocal_panel_blastp_annotated.tsv"
    warnings: List[dict] = []

    # For reciprocal parsing, qseqid is reference, sseqid is candidate. Need maps for reference query metadata.
    ref_id_to_gene = {r.id: infer_ref_gene(r) for r in refs}
    ref_id_to_species = {r.id: infer_ref_species(r) for r in refs}
    cand_id_to_species = {q.id: extract_species_from_query_id(q.id) for q in queries}
    cand_id_to_role = {q.id: extract_role_from_query_id(q.id) for q in queries}
    cand_id_to_isoform = {q.id: extract_isoform_from_query_id(q.id) for q in queries}

    if shutil.which("blastp") and shutil.which("makeblastdb"):
        run_blastp(reference_fasta, query_fasta, raw, outdir, f"{prefix}_candidate_blastdb", threads, evalue, max_targets=max(1000, len(queries)))
        # parse manually because subject IDs are candidates, not refs
        rows = []
        with raw.open("r", encoding="utf-8") as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 14:
                    continue
                qid, sid = c[0], c[1]
                rows.append({
                    "reference_id": qid,
                    "reference_gene": ref_id_to_gene.get(qid, "UNKNOWN"),
                    "reference_species": ref_id_to_species.get(qid, "unknown"),
                    "candidate_id": sid,
                    "candidate_species": cand_id_to_species.get(sid, "unknown"),
                    "candidate_role": cand_id_to_role.get(sid, ""),
                    "candidate_isoform": cand_id_to_isoform.get(sid, ""),
                    "pident": c[2], "length": c[3], "evalue": c[10], "bitscore": c[11],
                    "reference_coverage": f"{(float(c[3]) / float(c[12])) if float(c[12]) else 0:.5f}",
                    "candidate_coverage": f"{(float(c[3]) / float(c[13])) if float(c[13]) else 0:.5f}",
                })
    elif fallback:
        rows = []
        for ref in refs:
            for cand in queries:
                pident, qcov, scov, score = simple_similarity_score(ref.seq, cand.seq)
                rows.append({
                    "reference_id": ref.id, "reference_gene": ref_id_to_gene.get(ref.id, "UNKNOWN"),
                    "reference_species": ref_id_to_species.get(ref.id, "unknown"), "candidate_id": cand.id,
                    "candidate_species": cand_id_to_species.get(cand.id, "unknown"),
                    "candidate_role": cand_id_to_role.get(cand.id, ""), "candidate_isoform": cand_id_to_isoform.get(cand.id, ""),
                    "pident": f"{pident:.3f}", "length": int(min(len(ref.seq), len(cand.seq)) * qcov),
                    "evalue": "NA", "bitscore": f"{score:.3f}", "reference_coverage": f"{qcov:.5f}", "candidate_coverage": f"{scov:.5f}",
                })
        warnings.append({"warning_code": "reciprocal_fallback_pairwise_used", "severity": "warning", "message": "Reciprocal panel test used fallback pairwise mode.", "query_id": "", "subject_id": ""})
    else:
        raise RuntimeError("BLAST+ unavailable for reciprocal panel test.")

    rec_fields = ["reference_id", "reference_gene", "reference_species", "candidate_id", "candidate_species", "candidate_role", "candidate_isoform", "pident", "length", "evalue", "bitscore", "reference_coverage", "candidate_coverage"]
    write_tsv(annotated, rows, rec_fields)

    by_species_gene: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for r in rows:
        if r.get("reference_gene") in FGFR_GENES and r.get("candidate_species") != "unknown":
            by_species_gene[(r["candidate_species"], r["reference_gene"])].append(r)

    species_set = sorted({r["candidate_species"] for r in rows if r.get("candidate_species") != "unknown"})
    summary_rows: List[dict] = []
    for sp in species_set:
        gene_best = {}
        for gene in FGFR_GENES:
            candidates = by_species_gene.get((sp, gene), [])
            if candidates:
                best = max(candidates, key=lambda r: (safe_float(r["bitscore"]), safe_float(r["pident"])))
                gene_best[gene] = best
        scores = {g: safe_float(gene_best[g]["bitscore"]) for g in gene_best}
        ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
        top_gene, top_score = ranked[0] if ranked else ("NO_HIT", 0.0)
        second_gene, second_score = ranked[1] if len(ranked) > 1 else ("NA", 0.0)
        margin = (top_score - second_score) / top_score if top_score > 0 else 0.0
        fgfr2_top_candidate = gene_best.get("FGFR2", {}).get("candidate_id", "")
        if top_gene == "FGFR2" and margin >= 0.03:
            status = "reciprocal_FGFR2_supported"
        elif top_gene == "FGFR2":
            status = "reciprocal_FGFR2_low_margin"
        elif top_gene == "NO_HIT":
            status = "reciprocal_no_hit"
        else:
            status = "reciprocal_non_FGFR2_top"
        summary_rows.append({
            "species": sp,
            "reciprocal_top_reference_gene": top_gene,
            "reciprocal_second_reference_gene": second_gene,
            "reciprocal_margin_fraction": f"{margin:.5f}",
            "fgfr1_best_candidate_bitscore": f"{scores.get('FGFR1', 0.0):.3f}",
            "fgfr2_best_candidate_bitscore": f"{scores.get('FGFR2', 0.0):.3f}",
            "fgfr3_best_candidate_bitscore": f"{scores.get('FGFR3', 0.0):.3f}",
            "fgfr4_best_candidate_bitscore": f"{scores.get('FGFR4', 0.0):.3f}",
            "fgfr2_best_candidate_id": fgfr2_top_candidate,
            "reciprocal_status": status,
        })
    summary_path = outdir / f"{prefix}_reciprocal_species_summary.tsv"
    summary_fields = ["species", "reciprocal_top_reference_gene", "reciprocal_second_reference_gene", "reciprocal_margin_fraction", "fgfr1_best_candidate_bitscore", "fgfr2_best_candidate_bitscore", "fgfr3_best_candidate_bitscore", "fgfr4_best_candidate_bitscore", "fgfr2_best_candidate_id", "reciprocal_status"]
    write_tsv(summary_path, summary_rows, summary_fields)
    return annotated, summary_path, warnings


def pca_2d(matrix: List[List[float]]) -> Tuple[List[float], List[float]]:
    try:
        import numpy as np
        X = np.array(matrix, dtype=float)
        if X.ndim != 2 or X.shape[0] == 0:
            return [], []
        X = X - X.mean(axis=0)
        U, S, _Vt = np.linalg.svd(X, full_matrices=False)
        coords = U[:, :2] * S[:2]
        if coords.shape[1] == 1:
            return coords[:, 0].tolist(), [0.0] * coords.shape[0]
        return coords[:, 0].tolist(), coords[:, 1].tolist()
    except Exception:
        return [], []


def make_plots(validation_tsv: Path, species_summary_tsv: Path, outdir: Path, prefix: str, reciprocal_summary_tsv: Optional[Path] = None) -> List[str]:
    created: List[str] = []
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return created
    rows = load_rows(validation_tsv)
    if not rows:
        return created

    status_order = [
        "high_confidence_FGFR2", "probable_FGFR2", "FGFR2_best_hit_but_low_margin",
        "ambiguous_fgfr_family_member", "possible_paralog_misassignment",
        "insufficient_sequence_similarity", "no_detectable_fgfr_similarity",
    ]

    # 1) Protein-level scatter plot: coverage vs paralog margin.
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    markers = {"high_confidence_FGFR2": "o", "probable_FGFR2": "s", "FGFR2_best_hit_but_low_margin": "^", "ambiguous_fgfr_family_member": "D", "possible_paralog_misassignment": "X", "insufficient_sequence_similarity": "v", "no_detectable_fgfr_similarity": "P"}
    for status in status_order:
        xs = [safe_float(r.get("best_query_coverage")) for r in rows if r.get("fgfr2_screen_status") == status]
        ys = [safe_float(r.get("score_margin_fraction")) for r in rows if r.get("fgfr2_screen_status") == status]
        xs2, ys2 = zip(*[(x, y) for x, y in zip(xs, ys) if not math.isnan(x) and not math.isnan(y)]) if any(not math.isnan(x) and not math.isnan(y) for x, y in zip(xs, ys)) else ([], [])
        if xs2:
            ax.scatter(xs2, ys2, label=status.replace("_", " "), marker=markers.get(status, "o"), s=46, alpha=0.85)
    ax.axhline(0.08, linestyle="--", linewidth=1)
    ax.axhline(0.03, linestyle=":", linewidth=1)
    ax.set_xlabel("Best-hit query coverage")
    ax.set_ylabel("FGFR2 vs next-best paralog bitscore margin")
    ax.set_title("Protein-level FGFR2 paralog-screen support")
    ax.set_xlim(-0.02, 1.02)
    ax.legend(fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = outdir / f"{prefix}_paralog_screen_margin_scatter.{ext}"
        fig.savefig(p, dpi=300)
        created.append(str(p))
    plt.close(fig)

    # 2) Status count bar plot.
    counts = Counter(r.get("fgfr2_screen_status", "") for r in rows)
    labels = [s for s in status_order if counts.get(s, 0)]
    values = [counts[s] for s in labels]
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.bar(range(len(labels)), values)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([x.replace("_", "\n") for x in labels], fontsize=8)
    ax.set_ylabel("Number of selected proteins")
    ax.set_title("Paralog-screen/paralogy validation outcome")
    for i, v in enumerate(values):
        ax.text(i, v + max(values + [1]) * 0.01, str(v), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = outdir / f"{prefix}_fgfr2_screen_status_counts.{ext}"
        fig.savefig(p, dpi=300)
        created.append(str(p))
    plt.close(fig)

    # 3) Paralog score heatmap for all proteins.
    genes = list(FGFR_GENES)
    sorted_rows = sorted(rows, key=lambda r: (r.get("species_canonical", ""), r.get("isoform", ""), r.get("selection_role", ""), r.get("query_id", "")))
    matrix, ylabels = [], []
    for r in sorted_rows:
        scores = [safe_float(r.get(f"{g.lower()}_bitscore", "0")) for g in genes]
        scores = [0.0 if math.isnan(s) else s for s in scores]
        max_score = max(scores) if scores else 0.0
        matrix.append([(s / max_score if max_score > 0 else 0.0) for s in scores])
        sp = r.get("species_canonical") or extract_species_from_query_id(r.get("query_id", ""))
        role = (r.get("selection_role") or "").replace("FGFR2_", "").replace("_candidate", "")
        iso = r.get("isoform") or extract_isoform_from_query_id(r.get("query_id", ""))
        ylabels.append(f"{sp} | {role} | {iso}"[:95])
    if matrix:
        fig_h = max(8, 0.20 * len(matrix) + 2.0)
        fig, ax = plt.subplots(figsize=(7.2, fig_h))
        im = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, interpolation="nearest")
        ax.set_xticks(range(len(genes)))
        ax.set_xticklabels(genes)
        ax.set_yticks(range(len(ylabels)))
        ax.set_yticklabels(ylabels, fontsize=5.5 if len(ylabels) > 70 else 6.5)
        ax.set_xlabel("Reference paralog")
        ax.set_ylabel("Candidate protein")
        ax.set_title("Relative BLASTP support across FGFR paralogs")
        for i, row in enumerate(matrix):
            if row:
                j = max(range(len(row)), key=lambda k: row[k])
                ax.text(j, i, "●", ha="center", va="center", fontsize=3.8)
        cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
        cbar.set_label("Bitscore / best paralog bitscore")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            p = outdir / f"{prefix}_paralog_score_heatmap.{ext}"
            fig.savefig(p, dpi=300)
            created.append(str(p))
        plt.close(fig)

    # 4) Species-level minimum margin plot.
    sp_rows = load_rows(species_summary_tsv) if species_summary_tsv.exists() else []
    if sp_rows:
        sp_rows = sorted(sp_rows, key=lambda r: safe_float(r.get("min_score_margin_fraction")))
        species = [r.get("species", "unknown") for r in sp_rows]
        margins = [safe_float(r.get("min_score_margin_fraction")) for r in sp_rows]
        fig, ax = plt.subplots(figsize=(8.2, max(5.5, 0.28 * len(species))))
        ax.barh(species, margins)
        ax.axvline(0.08, linestyle="--", linewidth=1)
        ax.axvline(0.03, linestyle=":", linewidth=1)
        ax.set_xlabel("Minimum FGFR2 margin vs next-best FGFR paralog")
        ax.set_ylabel("Species")
        ax.set_title("Species-level FGFR2 paralog-screen robustness")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            p = outdir / f"{prefix}_species_min_margin_plot.{ext}"
            fig.savefig(p, dpi=300)
            created.append(str(p))
        plt.close(fig)

    # 5) Paralog support PCA-like plot from FGFR1-4 score vectors.
    score_matrix, labels_for_pca, status_for_pca = [], [], []
    for r in rows:
        scores = [safe_float(r.get(f"{g.lower()}_bitscore", "0")) for g in genes]
        scores = [0.0 if math.isnan(s) else s for s in scores]
        total = sum(scores)
        score_matrix.append([s / total if total > 0 else 0.0 for s in scores])
        labels_for_pca.append(r.get("species_canonical") or extract_species_from_query_id(r.get("query_id", "")))
        status_for_pca.append(r.get("fgfr2_screen_status", ""))
    pc1, pc2 = pca_2d(score_matrix)
    if pc1 and pc2:
        fig, ax = plt.subplots(figsize=(7.2, 5.8))
        for status in status_order:
            xs = [x for x, s in zip(pc1, status_for_pca) if s == status]
            ys = [y for y, s in zip(pc2, status_for_pca) if s == status]
            if xs:
                ax.scatter(xs, ys, label=status.replace("_", " "), s=46, alpha=0.85, marker=markers.get(status, "o"))
        ax.axhline(0, linewidth=0.7)
        ax.axvline(0, linewidth=0.7)
        ax.set_xlabel("PC1 of normalized FGFR paralog support")
        ax.set_ylabel("PC2 of normalized FGFR paralog support")
        ax.set_title("Candidate separation by FGFR paralog-support profile")
        ax.legend(fontsize=7, frameon=False, loc="best")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            p = outdir / f"{prefix}_paralog_support_pca.{ext}"
            fig.savefig(p, dpi=300)
            created.append(str(p))
        plt.close(fig)

    # 6) Reciprocal panel plot, if available.
    if reciprocal_summary_tsv and reciprocal_summary_tsv.exists():
        rec_rows = load_rows(reciprocal_summary_tsv)
        if rec_rows:
            rec_rows = sorted(rec_rows, key=lambda r: safe_float(r.get("reciprocal_margin_fraction")))
            species = [r.get("species", "unknown") for r in rec_rows]
            margins = [safe_float(r.get("reciprocal_margin_fraction")) for r in rec_rows]
            fig, ax = plt.subplots(figsize=(8.2, max(5.5, 0.28 * len(species))))
            ax.barh(species, margins)
            ax.axvline(0.03, linestyle=":", linewidth=1)
            ax.set_xlabel("Reciprocal FGFR2 margin vs next-best reference paralog")
            ax.set_ylabel("Species")
            ax.set_title("Reciprocal panel consistency by species")
            fig.tight_layout()
            for ext in ("png", "pdf"):
                p = outdir / f"{prefix}_reciprocal_species_margin_plot.{ext}"
                fig.savefig(p, dpi=300)
                created.append(str(p))
            plt.close(fig)

            # reciprocal heatmap genes x species
            rec_genes = ["fgfr1_best_candidate_bitscore", "fgfr2_best_candidate_bitscore", "fgfr3_best_candidate_bitscore", "fgfr4_best_candidate_bitscore"]
            mat = []
            for r in rec_rows:
                vals = [safe_float(r.get(c)) for c in rec_genes]
                vals = [0.0 if math.isnan(v) else v for v in vals]
                mx = max(vals) if vals else 0.0
                mat.append([v / mx if mx > 0 else 0.0 for v in vals])
            fig, ax = plt.subplots(figsize=(7.2, max(5.5, 0.28 * len(species))))
            im = ax.imshow(mat, aspect="auto", vmin=0, vmax=1, interpolation="nearest")
            ax.set_xticks(range(4))
            ax.set_xticklabels(["FGFR1", "FGFR2", "FGFR3", "FGFR4"])
            ax.set_yticks(range(len(species)))
            ax.set_yticklabels(species, fontsize=6)
            ax.set_xlabel("Reference paralog queried against candidate set")
            ax.set_ylabel("Species")
            ax.set_title("Reciprocal relative support by species")
            cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
            cbar.set_label("Bitscore / best reciprocal reference")
            fig.tight_layout()
            for ext in ("png", "pdf"):
                p = outdir / f"{prefix}_reciprocal_support_heatmap.{ext}"
                fig.savefig(p, dpi=300)
                created.append(str(p))
            plt.close(fig)

    return created


def write_markdown_report(path: Path, validation_rows: List[dict], risk_rows: List[dict], species_rows: List[dict], reciprocal_rows: List[dict], metadata: dict, plots: List[str]) -> None:
    counts = Counter(r.get("fgfr2_screen_status", "") for r in validation_rows)
    species_counts = Counter(r.get("species_fgfr2_screen_status", "") for r in species_rows)
    rec_counts = Counter(r.get("reciprocal_status", "") for r in reciprocal_rows) if reciprocal_rows else Counter()
    margins = [safe_float(r.get("score_margin_fraction")) for r in validation_rows]
    margins = [m for m in margins if not math.isnan(m)]
    min_margin = min(margins) if margins else math.nan
    next_best = Counter(r.get("second_ref_gene", "") for r in validation_rows)
    lines = [
        "# FGFR2 paralog identity screening report", "",
        f"**Script:** `{SCRIPT_NAME}` v{SCRIPT_VERSION}",
        f"**Run time UTC:** {metadata.get('run_datetime_utc')}",
        f"**Query proteins:** {metadata.get('query_count')}",
        f"**Candidate species:** {metadata.get('species_count')}",
        f"**Reference proteins:** {metadata.get('reference_count')}",
        f"**Reference genes:** {metadata.get('reference_gene_counts')}",
        f"**Reference species:** {metadata.get('reference_species_count')}",
        f"**Risk/uncertain proteins:** {len(risk_rows)}", "",
        "## Paralog-screen status counts", "", "| Status | Count |", "|---|---:|",
    ]
    for status, count in counts.most_common():
        lines.append(f"| {status} | {count} |")
    lines += ["", "## Species-level status counts", "", "| Species status | Count |", "|---|---:|"]
    for status, count in species_counts.most_common():
        lines.append(f"| {status} | {count} |")
    if rec_counts:
        lines += ["", "## Reciprocal panel status counts", "", "| Reciprocal status | Count |", "|---|---:|"]
        for status, count in rec_counts.most_common():
            lines.append(f"| {status} | {count} |")
    lines += ["", "## Key quantitative checks", ""]
    lines.append(f"- Minimum protein-level FGFR2 margin fraction: `{min_margin:.5f}`" if not math.isnan(min_margin) else "- Minimum protein-level FGFR2 margin fraction: `NA`")
    if next_best:
        lines.append(f"- Most common next-best paralog: `{next_best.most_common(1)[0][0]}`")
    lines += ["", "## Interpretation", ""]
    lines.append("Candidate proteins classified as high-confidence or probable FGFR2 can be used for the primary exon-domain analysis. Low-margin, ambiguous or non-FGFR2-best-hit cases should be displayed with uncertainty markers or moved to manual review.")
    if metadata.get("reference_species_count", 0) <= 1:
        lines.append("The current reference panel appears to contain one reference species. For the final thesis run, a multi-vertebrate FGFR1-4 panel is recommended to reduce human-reference bias.")
    if plots:
        lines += ["", "## Generated plots", ""]
        for p in plots:
            lines.append(f"- `{Path(p).name}`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Validate selected FGFR2 proteins against a paralog-aware FGFR reference panel using BLASTP.")
    ap.add_argument("--query_fasta", required=True, type=Path, help="Selected FGFR2 candidate proteins FASTA")
    ap.add_argument("--reference_fasta", required=True, type=Path, help="Reference FASTA containing FGFR1, FGFR2, FGFR3 and FGFR4 proteins; multi-vertebrate panel recommended")
    ap.add_argument("--protein_report", type=Path, default=None, help="Optional protein export report TSV")
    ap.add_argument("--outdir", required=True, type=Path, help="Output directory")
    ap.add_argument("--prefix", default="fgfr2", help="Output file prefix")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--evalue", default="1e-20")
    ap.add_argument("--fallback_pairwise", action="store_true", help="Debug fallback if BLAST+ is unavailable; not recommended for final analysis")
    ap.add_argument("--reciprocal_panel_test", action="store_true", help="Run reference-panel-vs-candidate-set reciprocal consistency test")

    ap.add_argument("--high_conf_margin_fraction", type=float, default=0.08)
    ap.add_argument("--probable_margin_fraction", type=float, default=0.03)
    ap.add_argument("--paralog_margin_fraction", type=float, default=0.03)
    ap.add_argument("--high_conf_qcov", type=float, default=0.70)
    ap.add_argument("--probable_qcov", type=float, default=0.55)
    ap.add_argument("--high_conf_pident", type=float, default=45.0)
    ap.add_argument("--probable_pident", type=float, default=35.0)
    ap.add_argument("--min_qcov_for_any_call", type=float, default=0.25)
    ap.add_argument("--min_pident_for_any_call", type=float, default=20.0)
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    queries = read_fasta(args.query_fasta)
    refs = read_fasta(args.reference_fasta)
    if not queries:
        raise ValueError(f"No query sequences found in {args.query_fasta}")
    if not refs:
        raise ValueError(f"No reference sequences found in {args.reference_fasta}")

    sseq_to_gene, sseq_to_species, warnings = check_reference_fasta(refs)
    warnings_path = args.outdir / f"{args.prefix}_paralog_screen_warnings.tsv"
    if any(w["severity"] == "error" for w in warnings):
        write_tsv(warnings_path, warnings, ["warning_code", "severity", "message", "query_id", "subject_id"])
        raise ValueError(f"Reference FASTA failed validation. See {warnings_path}")

    ref_gene_counts = Counter(sseq_to_gene.values())
    ref_species = sorted({sseq_to_species.get(r.id, "unknown") for r in refs})
    if len(ref_species) <= 1:
        warnings.append({
            "warning_code": "single_species_reference_panel",
            "severity": "note",
            "message": "Reference panel contains only one inferred species. A multi-vertebrate FGFR1-4 panel is recommended for the final thesis run.",
            "query_id": "", "subject_id": "",
        })

    raw_blast = args.outdir / f"{args.prefix}_fgfr_paralog_blastp.tsv"
    blast_with_gene = args.outdir / f"{args.prefix}_fgfr_paralog_blastp_with_gene.tsv"
    if shutil.which("blastp") and shutil.which("makeblastdb"):
        run_blastp(args.query_fasta, args.reference_fasta, raw_blast, args.outdir, f"{args.prefix}_reference_fgfr_blastdb", args.threads, args.evalue, max_targets=max(1000, len(refs)))
        hits = parse_blast_tsv(raw_blast, sseq_to_gene, sseq_to_species)
        write_annotated_blast(hits, blast_with_gene)
        method = "blastp"
    elif args.fallback_pairwise:
        run_fallback_pairwise(queries, refs, sseq_to_gene, sseq_to_species, blast_with_gene)
        hits = parse_blast_tsv(blast_with_gene, sseq_to_gene, sseq_to_species)
        warnings.append({"warning_code": "fallback_pairwise_used", "severity": "warning", "message": "BLAST+ was unavailable; used fallback pairwise mode. Do not use for final evidence.", "query_id": "", "subject_id": ""})
        method = "fallback_pairwise"
    else:
        raise RuntimeError("BLAST+ is unavailable. Install blastp/makeblastdb or rerun with --fallback_pairwise for debugging only.")

    protein_report = read_tsv_optional(args.protein_report)
    validation_rows, risk_rows = make_validation_rows(queries, hits, protein_report, args)

    validation_fields = [
        "query_id", "species_canonical", "species_input", "source_db", "selection_role", "isoform",
        "transcript_id_source", "translation_id_source", "protein_length_aa",
        "best_ref_gene", "best_ref_subject", "best_ref_species", "best_bitscore", "second_ref_gene", "second_bitscore",
        "fgfr1_bitscore", "fgfr2_bitscore", "fgfr3_bitscore", "fgfr4_bitscore",
        "fgfr1_best_subject", "fgfr2_best_subject", "fgfr3_best_subject", "fgfr4_best_subject",
        "fgfr1_best_ref_species", "fgfr2_best_ref_species", "fgfr3_best_ref_species", "fgfr4_best_ref_species",
        "fgfr1_best_pident", "fgfr2_best_pident", "fgfr3_best_pident", "fgfr4_best_pident",
        "fgfr1_best_qcov", "fgfr2_best_qcov", "fgfr3_best_qcov", "fgfr4_best_qcov",
        "fgfr1_n_reference_hits", "fgfr2_n_reference_hits", "fgfr3_n_reference_hits", "fgfr4_n_reference_hits",
        "fgfr1_n_reference_species_hit", "fgfr2_n_reference_species_hit", "fgfr3_n_reference_species_hit", "fgfr4_n_reference_species_hit",
        "score_margin_abs", "score_margin_fraction", "best_pident", "best_query_coverage", "best_subject_coverage", "best_evalue",
        "fgfr2_screen_status", "primary_analysis_recommendation",
    ]
    validation_path = args.outdir / f"{args.prefix}_paralog_screen_validation.tsv"
    risk_path = args.outdir / f"{args.prefix}_paralog_risk.tsv"
    species_path = args.outdir / f"{args.prefix}_paralog_screen_species_summary.tsv"
    metadata_path = args.outdir / f"{args.prefix}_paralog_screen_metadata.json"
    md_path = args.outdir / f"{args.prefix}_paralog_screen_report.md"

    write_tsv(validation_path, validation_rows, validation_fields)
    write_tsv(risk_path, risk_rows, validation_fields)
    species_rows = species_summary(validation_rows)
    species_fields = ["species", "n_proteins", "n_high_confidence_FGFR2", "n_probable_FGFR2", "min_score_margin_fraction", "min_query_coverage", "best_ref_gene_counts", "second_ref_gene_counts", "statuses", "species_fgfr2_screen_status"]
    write_tsv(species_path, species_rows, species_fields)

    reciprocal_annotated: Optional[Path] = None
    reciprocal_summary: Optional[Path] = None
    reciprocal_rows: List[dict] = []
    if args.reciprocal_panel_test:
        reciprocal_annotated, reciprocal_summary, rec_warnings = reciprocal_panel_test(refs, queries, args.reference_fasta, args.query_fasta, args.outdir, args.prefix, args.threads, args.evalue, args.fallback_pairwise, sseq_to_gene, sseq_to_species)
        warnings.extend(rec_warnings)
        reciprocal_rows = load_rows(reciprocal_summary) if reciprocal_summary and reciprocal_summary.exists() else []

    write_tsv(warnings_path, warnings, ["warning_code", "severity", "message", "query_id", "subject_id"])
    plots = make_plots(validation_path, species_path, args.outdir, args.prefix, reciprocal_summary)

    metadata = {
        "script_name": SCRIPT_NAME,
        "script_version": SCRIPT_VERSION,
        "run_datetime_utc": datetime.now(timezone.utc).isoformat(),
        "query_fasta": str(args.query_fasta),
        "reference_fasta": str(args.reference_fasta),
        "protein_report": str(args.protein_report or ""),
        "outdir": str(args.outdir),
        "query_count": len(queries),
        "species_count": len({r.get("species_canonical") or extract_species_from_query_id(r.get("query_id", "")) for r in validation_rows}),
        "reference_count": len(refs),
        "reference_gene_counts": dict(ref_gene_counts),
        "reference_species_count": len(ref_species),
        "reference_species": ref_species,
        "method": method,
        "reciprocal_panel_test": bool(args.reciprocal_panel_test),
        "thresholds": {
            "high_conf_margin_fraction": args.high_conf_margin_fraction,
            "probable_margin_fraction": args.probable_margin_fraction,
            "paralog_margin_fraction": args.paralog_margin_fraction,
            "high_conf_qcov": args.high_conf_qcov,
            "probable_qcov": args.probable_qcov,
            "high_conf_pident": args.high_conf_pident,
            "probable_pident": args.probable_pident,
            "min_qcov_for_any_call": args.min_qcov_for_any_call,
            "min_pident_for_any_call": args.min_pident_for_any_call,
        },
        "outputs": {
            "validation_tsv": str(validation_path),
            "risk_tsv": str(risk_path),
            "species_summary_tsv": str(species_path),
            "warnings_tsv": str(warnings_path),
            "blast_tsv": str(raw_blast if raw_blast.exists() else blast_with_gene),
            "blast_with_gene_tsv": str(blast_with_gene),
            "reciprocal_blast_annotated_tsv": str(reciprocal_annotated or ""),
            "reciprocal_species_summary_tsv": str(reciprocal_summary or ""),
            "plots": plots,
        },
        "status_counts": dict(Counter(r.get("fgfr2_screen_status", "") for r in validation_rows)),
        "species_status_counts": dict(Counter(r.get("species_fgfr2_screen_status", "") for r in species_rows)),
        "reciprocal_status_counts": dict(Counter(r.get("reciprocal_status", "") for r in reciprocal_rows)),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    write_markdown_report(md_path, validation_rows, risk_rows, species_rows, reciprocal_rows, metadata, plots)

    print(f"Wrote validation table: {validation_path}")
    print(f"Wrote species summary: {species_path}")
    print(f"Wrote paralog-risk table: {risk_path}")
    if reciprocal_summary:
        print(f"Wrote reciprocal summary: {reciprocal_summary}")
    print(f"Wrote report: {md_path}")
    if plots:
        print("Wrote plots:")
        for p in plots:
            print(f"  {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
