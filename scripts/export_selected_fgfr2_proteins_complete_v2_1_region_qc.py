#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
import re
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import unquote
from urllib.request import Request, urlopen

SCRIPT_NAME = "export_selected_fgfr2_proteins_complete_v2.py"
SCRIPT_VERSION = "2.1"
DEFAULT_ROLES = "reference,FGFR2_IIIb_candidate,FGFR2_IIIc_candidate"
AA_ALLOWED = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")

REQUIRED_SELECTED_COLUMNS = {
    "species_canonical",
    "source_db",
    "selection_role",
    "transcript_id_source",
}
OPTIONAL_BUT_IMPORTANT_COLUMNS = {
    "species_input",
    "translation_id_source",
    "protein_length_aa",
    "iii_isoform_assignment",
    "internal_transcript_id",
}


@dataclass
class FastaRecord:
    accession: str
    accession_no_version: str
    header: str
    sequence: str
    path: str
    product_lower: str
    bracket_species_lower: str


@dataclass
class GffProteinMap:
    protein_accession: str = ""
    product: str = ""
    source_gff: str = ""


@dataclass
class SequenceCheck:
    status: str
    warnings: List[str]


def read_tsv(path: Path) -> List[dict]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return list(reader)


def tsv_fieldnames(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader.fieldnames or [])


def write_tsv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def add_warning(rows: List[dict], code: str, severity: str, message: str, row: Optional[dict] = None, **extra) -> None:
    base = {
        "warning_code": code,
        "severity": severity,
        "message": message,
        "species_canonical": "",
        "source_db": "",
        "selection_role": "",
        "transcript_id_source": "",
        "output_id": "",
    }
    if row:
        base.update({
            "species_canonical": row.get("species_canonical", ""),
            "source_db": row.get("source_db", row.get("source", "")),
            "selection_role": row.get("selection_role", row.get("role", "")),
            "transcript_id_source": row.get("transcript_id_source", row.get("transcript_id", "")),
        })
    base.update({k: str(v) for k, v in extra.items()})
    rows.append(base)


def validate_columns(path: Path, required: set[str], strict: bool, warnings: List[dict]) -> List[str]:
    fields = tsv_fieldnames(path)
    missing = sorted(required - set(fields))
    if missing:
        msg = f"Missing required columns in {path.name}: {', '.join(missing)}"
        add_warning(warnings, "missing_required_column", "error", msg, missing_columns=",".join(missing), file=str(path))
        if strict:
            raise ValueError(msg)
    optional_missing = sorted(OPTIONAL_BUT_IMPORTANT_COLUMNS - set(fields))
    for col in optional_missing:
        add_warning(warnings, "missing_optional_column", "info", f"Optional column {col} is missing from {path.name}.", file=str(path), column=col)
    return fields


def clean_accession(x: str) -> str:
    x = str(x or "").strip()
    if not x or x.lower() in {"nan", "none", "null", "na"}:
        return ""
    x = x.split()[0]
    for prefix in ("rna-", "cds-", "protein-", "transcript:", "protein:", "RefSeq:", "Genbank:", "NCBI:"):
        if x.startswith(prefix):
            x = x[len(prefix):]
    return x.strip()


def no_version(acc: str) -> str:
    acc = clean_accession(acc)
    return re.sub(r"\.\d+$", "", acc)


def parse_int_maybe(x: str) -> Optional[int]:
    s = str(x or "").strip()
    if not s or s.lower() in {"nan", "none", "null", "na"}:
        return None
    try:
        return int(float(s))
    except Exception:
        return None


def sanitize_header_value(x: str) -> str:
    s = str(x or "NA").strip()
    if not s:
        s = "NA"
    return re.sub(r"[^A-Za-z0-9_.:=,\-]+", "_", s)


def parse_attrs(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in attr.strip().split(";"):
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k] = unquote(v)
    return out


def choose_gff_paths(cache: Path) -> List[Path]:
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


def parse_gff3_transcript_to_protein(cache: Path) -> Dict[str, GffProteinMap]:
    mapping: Dict[str, GffProteinMap] = {}
    for gff in choose_gff_paths(cache):
        try:
            with gff.open("r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    if not line or line.startswith("#"):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 9 or cols[2].lower() not in {"cds", "protein"}:
                        continue
                    attrs = parse_attrs(cols[8])
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
                            if key and key not in mapping:
                                mapping[key] = GffProteinMap(protein, product, str(gff))
        except OSError:
            continue
    return mapping


def parse_fasta(path: Path) -> Iterable[Tuple[str, str]]:
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


def bracket_species(header: str) -> str:
    m = re.search(r"\[([^\]]+)\]\s*$", header)
    return m.group(1).lower() if m else ""


def product_text(header: str) -> str:
    h = re.sub(r"\s*\[[^\]]+\]\s*$", "", header).strip()
    parts = h.split(maxsplit=1)
    return parts[1].lower() if len(parts) > 1 else ""


def index_protein_fastas(cache: Path) -> Tuple[Dict[str, FastaRecord], Dict[str, List[FastaRecord]]]:
    by_acc: Dict[str, FastaRecord] = {}
    by_species: Dict[str, List[FastaRecord]] = {}
    fasta_paths = [p for p in cache.rglob("protein.faa")]
    if not fasta_paths:
        fasta_paths = [p for p in cache.rglob("*") if p.is_file() and p.suffix.lower() in {".faa", ".fa", ".fasta"}]
    for fasta in sorted(fasta_paths):
        for header, seq in parse_fasta(fasta):
            first = clean_accession(header.split()[0] if header else "")
            if not first:
                continue
            rec = FastaRecord(
                accession=first,
                accession_no_version=no_version(first),
                header=header,
                sequence=seq,
                path=str(fasta),
                product_lower=product_text(header),
                bracket_species_lower=bracket_species(header),
            )
            for k in {rec.accession, rec.accession_no_version}:
                if k and k not in by_acc:
                    by_acc[k] = rec
            by_species.setdefault(rec.bracket_species_lower, []).append(rec)
    return by_acc, by_species


def species_to_bracket_candidates(species_input: str, species_canonical: str) -> List[str]:
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
            out.append(v); seen.add(v)
    return out


def length_score(seq_len: int, expected_len: Optional[int]) -> int:
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


def find_by_product_fallback(by_species: Dict[str, List[FastaRecord]], species_input: str, species_canonical: str, product: str, expected_len: Optional[int]) -> Optional[FastaRecord]:
    prod_l = (product or "").lower().strip()
    candidates: List[FastaRecord] = []
    for sp in species_to_bracket_candidates(species_input, species_canonical):
        for bracket, recs in by_species.items():
            if bracket == sp or bracket.startswith(sp) or sp.startswith(bracket):
                candidates.extend(recs)
    if not candidates:
        return None
    scored: List[Tuple[int, FastaRecord]] = []
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
        score += length_score(len(rec.sequence), expected_len)
        scored.append((score, rec))
    if not scored:
        return None
    scored.sort(key=lambda x: (x[0], -abs(len(x[1].sequence) - expected_len) if expected_len else 0), reverse=True)
    return scored[0][1] if scored[0][0] >= 50 else None


def fetch_ensembl_protein(translation_id: str, retries: int = 3, sleep: float = 0.4, timeout: int = 30) -> Tuple[str, str]:
    tid = clean_accession(translation_id)
    if not tid:
        return "", "missing_translation_id"
    url = f"https://rest.ensembl.org/sequence/id/{tid}?type=protein"
    headers = {"Content-Type": "text/plain", "User-Agent": "FGFR2-boundary-mapping-bachelor-thesis/1.0"}
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
            detail = f"ensembl_rest_http_{e.code}"
            if e.code in {429, 500, 502, 503, 504} and attempt < retries:
                time.sleep(sleep * attempt)
                continue
            return "", detail
        except URLError as e:
            if attempt < retries:
                time.sleep(sleep * attempt)
                continue
            return "", f"ensembl_rest_urlerror:{e.reason}"
        except Exception as e:
            if attempt < retries:
                time.sleep(sleep * attempt)
                continue
            return "", f"ensembl_rest_error:{type(e).__name__}:{e}"
    return "", "ensembl_rest_failed"


def check_sequence(seq: str, expected_len: Optional[int], min_len: int, max_len: int) -> SequenceCheck:
    warnings: List[str] = []
    if not seq:
        return SequenceCheck("missing", ["missing_sequence"])
    seq_u = seq.upper().replace("*", "")
    invalid = sorted(set(seq_u) - AA_ALLOWED)
    if invalid:
        warnings.append("invalid_amino_acid_characters:" + "".join(invalid))
    dna_chars = set("ACGTN")
    if len(seq_u) > 0 and sum(1 for c in seq_u if c in dna_chars) / len(seq_u) > 0.9:
        warnings.append("sequence_looks_like_dna")
    if "*" in seq[:-1]:
        warnings.append("internal_stop_codon")
    if len(seq_u) < min_len:
        warnings.append("protein_too_short")
    if len(seq_u) > max_len:
        warnings.append("protein_too_long")
    if expected_len is not None:
        diff = abs(len(seq_u) - expected_len)
        frac = diff / expected_len if expected_len else 0
        if diff == 0:
            pass
        elif diff <= 2:
            warnings.append("minor_length_difference")
        elif frac <= 0.05:
            warnings.append("moderate_length_difference")
        else:
            warnings.append("major_length_difference")
    return SequenceCheck("pass" if not warnings else "warning", warnings)


def length_check_status(expected: Optional[int], observed: Optional[int]) -> Tuple[str, str, str]:
    if observed is None:
        return "missing_observed_length", "", ""
    if expected is None:
        return "missing_expected_length", "", ""
    diff = observed - expected
    frac = abs(diff) / expected if expected else 0
    if diff == 0:
        status = "exact"
    elif abs(diff) <= 2:
        status = "minor_difference"
    elif frac <= 0.05:
        status = "moderate_difference"
    else:
        status = "major_difference"
    return status, str(diff), f"{frac:.6f}"


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_fasta(path: Path, records: List[Tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for h, s in records:
            f.write(f">{h}\n{wrap(s)}\n")


def write_markdown_report(path: Path, report_rows: List[dict], warnings: List[dict], metadata: dict) -> None:
    counts: Dict[str, int] = {}
    methods: Dict[str, int] = {}
    for r in report_rows:
        counts[r.get("match_status", "")] = counts.get(r.get("match_status", ""), 0) + 1
        methods[r.get("match_method", "")] = methods.get(r.get("match_method", ""), 0) + 1
    lines = [
        "# FGFR2 protein export report",
        "",
        f"Script: `{metadata.get('script_name')}` version `{metadata.get('script_version')}`",
        f"Run time: {metadata.get('run_datetime_utc')}",
        "",
        "## Summary",
        "",
        f"Selected rows considered: **{metadata.get('selected_rows_considered')}**",
        f"FASTA records written: **{metadata.get('fasta_records_written')}**",
        f"Warnings: **{len(warnings)}**",
        "",
        "## Match status counts",
        "",
    ]
    for k, v in sorted(counts.items()):
        lines.append(f"- `{k}`: {v}")
    lines += ["", "## Match method counts", ""]
    for k, v in sorted(methods.items()):
        lines.append(f"- `{k}`: {v}")
    lines += ["", "## Scientific interpretation", "", "Exact accession matches are treated as high-confidence sequence exports. Product/species/length rescue matches are retained for completeness but flagged as medium-confidence and should be inspected before final biological interpretation.", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_html_report(path: Path, md_path: Path) -> None:
    text = md_path.read_text(encoding="utf-8")
    html_lines = ["<html><head><meta charset='utf-8'><title>FGFR2 protein export report</title></head><body>"]
    for line in text.splitlines():
        if line.startswith("# "):
            html_lines.append(f"<h1>{line[2:]}</h1>")
        elif line.startswith("## "):
            html_lines.append(f"<h2>{line[3:]}</h2>")
        elif line.startswith("- "):
            html_lines.append(f"<p>{line}</p>")
        elif not line.strip():
            html_lines.append("<br/>")
        else:
            html_lines.append(f"<p>{line}</p>")
    html_lines.append("</body></html>")
    path.write_text("\n".join(html_lines), encoding="utf-8")


def _norm_isoform(value: object) -> str:
    v = str(value or "").strip().lower()
    if "iiib" in v or v == "3b":
        return "IIIb"
    if "iiic" in v or v == "3c":
        return "IIIc"
    return "unclassified"


def _extract_header_field(header: str, key: str) -> str:
    for part in str(header).split("|"):
        if part.startswith(key + "="):
            return part.split("=", 1)[1]
    return ""


def _hamming_identity(a: str, b: str) -> Tuple[float, int, int]:
    n = min(len(a), len(b))
    if n <= 0:
        return 0.0, 0, 0
    matches = sum(1 for i in range(n) if a[i] == b[i])
    mismatches = n - matches
    return matches / n, matches, mismatches


def make_iii_region_export_qc(report_rows: List[dict], seq_by_output_id: Dict[str, str], start_1based: int, end_1based: int) -> List[dict]:
    w0 = max(0, int(start_1based) - 1)
    w1 = max(w0, int(end_1based))
    rows: List[dict] = []

    # First pass: one row per exported protein.
    for r in report_rows:
        output_id = str(r.get("output_id", ""))
        seq = seq_by_output_id.get(output_id, "")
        role = str(r.get("selection_role", ""))
        iso = _norm_isoform(role or r.get("output_header", ""))
        if iso == "unclassified":
            iso = _norm_isoform(_extract_header_field(str(r.get("output_header", "")), "isoform"))
        window = seq[w0:w1] if seq else ""
        expected_len = w1 - w0
        coverage = len(window) / expected_len if expected_len else 0.0
        if not seq:
            window_status = "no_sequence"
        elif coverage >= 0.95:
            window_status = "window_complete"
        elif coverage >= 0.50:
            window_status = "window_partial"
        else:
            window_status = "window_missing_or_too_short"
        rows.append({
            "output_id": output_id,
            "species_canonical": r.get("species_canonical", ""),
            "source_db": r.get("source_db", ""),
            "selection_role": role,
            "inferred_isoform_from_role": iso,
            "transcript_id_source": r.get("transcript_id_source", ""),
            "protein_accession_matched": r.get("protein_accession_matched", ""),
            "protein_length": len(seq) if seq else "",
            "iii_region_window_start_1based": start_1based,
            "iii_region_window_end_1based": end_1based,
            "iii_region_window_length": len(window),
            "iii_region_window_coverage": f"{coverage:.4f}",
            "iii_region_window_status": window_status,
            "iii_region_window_sequence": window,
            "same_as_reference": "",
            "same_as_other_isoform": "",
            "pair_window_identity": "",
            "pair_window_mismatches": "",
            "pair_distinguishability_status": "not_evaluated",
        })

    # Second pass: add species-local pair information and reference identity.
    by_species: Dict[str, List[dict]] = defaultdict(list)
    for row in rows:
        by_species[str(row.get("species_canonical", ""))].append(row)

    for _sp, srows in by_species.items():
        refs = [x for x in srows if str(x.get("selection_role", "")).lower() == "reference"]
        ref_seq = refs[0].get("iii_region_window_sequence", "") if refs else ""
        iiib = next((x for x in srows if x.get("inferred_isoform_from_role") == "IIIb" and "candidate" in str(x.get("selection_role", "")).lower()), None)
        iiic = next((x for x in srows if x.get("inferred_isoform_from_role") == "IIIc" and "candidate" in str(x.get("selection_role", "")).lower()), None)
        pair_identity = ""
        pair_mismatches = ""
        pair_status = "missing_pair_member"
        if iiib and iiic:
            idv, _matches, mismatches = _hamming_identity(str(iiib.get("iii_region_window_sequence", "")), str(iiic.get("iii_region_window_sequence", "")))
            pair_identity = f"{idv:.4f}"
            pair_mismatches = str(mismatches + abs(len(str(iiib.get("iii_region_window_sequence", ""))) - len(str(iiic.get("iii_region_window_sequence", "")))))
            if idv >= 0.97 and int(pair_mismatches) <= 5:
                pair_status = "III_region_nearly_identical"
            else:
                pair_status = "III_region_sequence_distinct"
        for row in srows:
            wseq = str(row.get("iii_region_window_sequence", ""))
            row["same_as_reference"] = "1" if ref_seq and wseq == ref_seq else "0"
            if row.get("inferred_isoform_from_role") == "IIIb" and iiic:
                row["same_as_other_isoform"] = "1" if wseq == str(iiic.get("iii_region_window_sequence", "")) else "0"
            elif row.get("inferred_isoform_from_role") == "IIIc" and iiib:
                row["same_as_other_isoform"] = "1" if wseq == str(iiib.get("iii_region_window_sequence", "")) else "0"
            else:
                row["same_as_other_isoform"] = ""
            row["pair_window_identity"] = pair_identity
            row["pair_window_mismatches"] = pair_mismatches
            row["pair_distinguishability_status"] = pair_status
    return rows

def export_selected(args: argparse.Namespace) -> Tuple[List[dict], List[dict], dict]:
    warnings: List[dict] = []
    fields = validate_columns(args.selected, REQUIRED_SELECTED_COLUMNS, args.strict, warnings)
    selected = read_tsv(args.selected)
    roles = {x.strip() for x in args.roles.split(",") if x.strip()}
    selected_rows = [r for r in selected if (r.get("selection_role") or r.get("role") or "") in roles]
    if not selected_rows:
        add_warning(warnings, "no_selected_rows_for_roles", "error" if args.strict else "warning", "No rows matched requested selection roles.", roles=",".join(sorted(roles)))
        if args.strict:
            raise ValueError("No selected rows matched requested roles")

    tx_to_protein = parse_gff3_transcript_to_protein(args.cache) if args.cache else {}
    by_acc, by_species = index_protein_fastas(args.cache) if args.cache else ({}, {})

    fasta_records: List[Tuple[str, str]] = []
    report_rows: List[dict] = []
    seq_by_output_id: Dict[str, str] = {}
    seen_headers: set[str] = set()

    for idx, row in enumerate(selected_rows, start=1):
        species_input = row.get("species_input") or row.get("species") or row.get("species_name") or row.get("species_canonical", "")
        species_canonical = row.get("species_canonical") or species_input.replace(" ", "_").lower()
        source_db = row.get("source_db") or row.get("source") or ""
        role = row.get("selection_role") or row.get("role") or ""
        tx = clean_accession(row.get("transcript_id_source") or row.get("transcript_id") or "")
        translation = clean_accession(row.get("translation_id_source") or row.get("protein_id") or row.get("translation_id") or "")
        expected_len = parse_int_maybe(row.get("protein_length_aa") or row.get("expected_length_aa"))

        output_id = f"fgfr2prot_{idx:06d}"
        mapped = GffProteinMap()
        for key in [tx, no_version(tx)]:
            if key in tx_to_protein:
                mapped = tx_to_protein[key]
                break
        requested_protein = translation or mapped.protein_accession
        requested_clean = clean_accession(requested_protein)
        product = mapped.product

        status = "not_found"
        status_detail = ""
        seq = ""
        source_header = ""
        matched_accession = ""
        source_method = ""
        match_confidence = "failed"

        if source_db.upper() == "NCBI":
            if not requested_clean:
                add_warning(warnings, "ncbi_missing_protein_accession", "warning", "NCBI row has no translation ID and no GFF3-derived protein accession.", row, output_id=output_id)
            for key in [requested_clean, no_version(requested_clean)]:
                rec = by_acc.get(key)
                if rec:
                    seq = rec.sequence
                    source_header = rec.header
                    matched_accession = rec.accession
                    status = "matched"
                    source_method = "ncbi_exact_accession"
                    match_confidence = "high"
                    break
            if not seq and not args.disable_ncbi_rescue:
                rec = find_by_product_fallback(by_species, species_input, species_canonical, product, expected_len)
                if rec:
                    seq = rec.sequence
                    source_header = rec.header
                    matched_accession = rec.accession
                    status = "matched"
                    source_method = "ncbi_product_species_length_rescue"
                    match_confidence = "medium"
                    status_detail = "exact_accession_not_found_but_product_rescue_used"
                    add_warning(warnings, "ncbi_rescue_match_used", "warning", "NCBI product/species/length rescue was used instead of exact accession matching.", row, output_id=output_id, matched_accession=matched_accession)
            if not seq:
                status = "not_found"
                source_method = "ncbi_failed"
                status_detail = "protein_not_found_in_fasta"
                add_warning(warnings, "protein_not_found", "error" if args.strict else "warning", "Protein sequence could not be found for NCBI row.", row, output_id=output_id, requested_protein=requested_clean)
        else:
            if translation and not args.no_ensembl_rest:
                seq, detail = fetch_ensembl_protein(translation, sleep=args.ensembl_sleep, timeout=args.ensembl_timeout)
                if seq:
                    status = "matched"
                    source_method = "ensembl_rest_translation"
                    status_detail = detail
                    matched_accession = translation
                    match_confidence = "high"
                    source_header = f"{translation} Ensembl REST protein sequence"
                else:
                    status = "not_found"
                    source_method = "ensembl_rest_failed"
                    status_detail = detail
                    match_confidence = "failed"
                    add_warning(warnings, "ensembl_rest_failed", "error" if args.strict else "warning", f"Ensembl REST sequence retrieval failed: {detail}", row, output_id=output_id)
            else:
                status = "not_found"
                source_method = "ensembl_no_rest_or_missing_translation"
                status_detail = "no_ensembl_rest_enabled_or_missing_translation_id"
                add_warning(warnings, "ensembl_no_rest_or_missing_translation", "error" if args.strict else "warning", "Ensembl sequence was not retrieved because REST is disabled or translation ID is missing.", row, output_id=output_id)

        observed_len = len(seq.replace("*", "")) if seq else None
        len_status, len_diff, len_frac = length_check_status(expected_len, observed_len)
        seq_check = check_sequence(seq, expected_len, args.min_protein_len, args.max_protein_len)
        if seq and seq_check.warnings:
            for w in seq_check.warnings:
                severity = "warning"
                if w in {"invalid_amino_acid_characters", "sequence_looks_like_dna", "internal_stop_codon", "major_length_difference"}:
                    severity = "error" if args.strict else "warning"
                add_warning(warnings, w, severity, f"Sequence validation warning: {w}", row, output_id=output_id, observed_length_aa=observed_len, expected_length_aa=expected_len or "")

        output_header = ""
        if seq:
            metadata_header = (
                f"{output_id}|species={sanitize_header_value(species_canonical)}|source={sanitize_header_value(source_db)}"
                f"|role={sanitize_header_value(role)}|transcript={sanitize_header_value(tx)}"
                f"|protein={sanitize_header_value(matched_accession or requested_clean)}"
                f"|isoform={sanitize_header_value(row.get('iii_isoform_assignment',''))}"
            )
            if metadata_header in seen_headers:
                add_warning(warnings, "duplicate_fasta_header", "error" if args.strict else "warning", "Duplicate FASTA header detected.", row, output_id=output_id)
            seen_headers.add(metadata_header)
            output_header = metadata_header
            clean_seq_for_output = seq.replace("*", "")
            fasta_records.append((output_header, clean_seq_for_output))
            seq_by_output_id[output_id] = clean_seq_for_output

        report_rows.append({
            "output_id": output_id,
            "row_number": str(idx),
            "species_input": species_input,
            "species_canonical": species_canonical,
            "source_db": source_db,
            "selection_role": role,
            "internal_transcript_id": row.get("internal_transcript_id", ""),
            "transcript_id_source": tx,
            "translation_id_source": translation,
            "protein_accession_mapped": mapped.protein_accession,
            "protein_accession_requested": requested_clean,
            "protein_accession_matched": matched_accession,
            "protein_product_from_gff": product,
            "expected_length_aa": str(expected_len or ""),
            "observed_length_aa": str(observed_len or ""),
            "length_difference_aa": len_diff,
            "length_difference_fraction": len_frac,
            "length_check_status": len_status,
            "sequence_check_status": seq_check.status,
            "sequence_check_warnings": ";".join(seq_check.warnings),
            "match_status": status,
            "match_method": source_method,
            "match_confidence": match_confidence,
            "status_detail": status_detail,
            "source_fasta_or_rest_header": source_header,
            "output_header": output_header,
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(args.out, fasta_records)
    report_fields = [
        "output_id", "row_number", "species_input", "species_canonical", "source_db", "selection_role",
        "internal_transcript_id", "transcript_id_source", "translation_id_source", "protein_accession_mapped",
        "protein_accession_requested", "protein_accession_matched", "protein_product_from_gff",
        "expected_length_aa", "observed_length_aa", "length_difference_aa", "length_difference_fraction",
        "length_check_status", "sequence_check_status", "sequence_check_warnings", "match_status", "match_method",
        "match_confidence", "status_detail", "source_fasta_or_rest_header", "output_header",
    ]
    write_tsv(args.report, report_rows, report_fields)
    warning_fields = ["warning_code", "severity", "message", "species_canonical", "source_db", "selection_role", "transcript_id_source", "output_id", "file", "column", "missing_columns", "roles", "requested_protein", "matched_accession", "observed_length_aa", "expected_length_aa"]
    write_tsv(args.warnings, warnings, warning_fields)
    region_qc_rows = make_iii_region_export_qc(report_rows, seq_by_output_id, args.iii_region_start, args.iii_region_end)
    region_qc_fields = [
        "output_id", "species_canonical", "source_db", "selection_role", "inferred_isoform_from_role",
        "transcript_id_source", "protein_accession_matched", "protein_length",
        "iii_region_window_start_1based", "iii_region_window_end_1based", "iii_region_window_length",
        "iii_region_window_coverage", "iii_region_window_status", "iii_region_window_sequence",
        "same_as_reference", "same_as_other_isoform", "pair_window_identity", "pair_window_mismatches", "pair_distinguishability_status",
    ]
    write_tsv(args.region_qc, region_qc_rows, region_qc_fields)

    metadata = {
        "script_name": SCRIPT_NAME,
        "script_version": SCRIPT_VERSION,
        "run_datetime_utc": datetime.now(timezone.utc).isoformat(),
        "selected": str(args.selected),
        "cache": str(args.cache),
        "output_fasta": str(args.out),
        "report": str(args.report),
        "warnings": str(args.warnings),
        "region_qc": str(args.region_qc),
        "roles": sorted(roles),
        "no_ensembl_rest": bool(args.no_ensembl_rest),
        "disable_ncbi_rescue": bool(args.disable_ncbi_rescue),
        "min_protein_len": args.min_protein_len,
        "max_protein_len": args.max_protein_len,
        "selected_rows_total": len(selected),
        "selected_rows_considered": len(selected_rows),
        "fasta_records_written": len(fasta_records),
        "report_rows": len(report_rows),
        "region_qc_rows": len(region_qc_rows),
        "warning_rows": len(warnings),
        "ncbi_transcript_protein_mappings": len(tx_to_protein),
        "ncbi_fasta_accession_keys": len(by_acc),
        "ncbi_fasta_species_groups": len(by_species),
        "input_columns": fields,
    }
    args.metadata.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    write_markdown_report(args.md_report, report_rows, warnings, metadata)
    write_html_report(args.html_report, args.md_report)

    if args.strict and any(w.get("severity") == "error" for w in warnings):
        raise RuntimeError("Strict mode failed because error-level warnings were produced. See protein_export_warnings.tsv.")
    return report_rows, warnings, metadata


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Export selected FGFR2 proteins for InterProScan with explicit provenance and quality checks.")
    ap.add_argument("--selected", required=True, type=Path, help="selected_transcripts.tsv from transcript selection step")
    ap.add_argument("--cache", required=True, type=Path, help="NCBI datasets cache directory from model collection step")
    ap.add_argument("--out", required=True, type=Path, help="Output FASTA")
    ap.add_argument("--report", required=True, type=Path, help="Protein export mapping/report TSV")
    ap.add_argument("--warnings", type=Path, default=None, help="Warnings TSV; default next to report")
    ap.add_argument("--region_qc", type=Path, default=None, help="III-region fixed-window protein QC TSV; default next to report")
    ap.add_argument("--iii_region_start", type=int, default=250, help="1-based start of fixed protein III-region QC window")
    ap.add_argument("--iii_region_end", type=int, default=430, help="1-based end of fixed protein III-region QC window")
    ap.add_argument("--metadata", type=Path, default=None, help="Run metadata JSON; default next to report")
    ap.add_argument("--md_report", type=Path, default=None, help="Markdown report; default next to report")
    ap.add_argument("--html_report", type=Path, default=None, help="HTML report; default next to report")
    ap.add_argument("--roles", default=DEFAULT_ROLES)
    ap.add_argument("--no_ensembl_rest", action="store_true")
    ap.add_argument("--disable_ncbi_rescue", action="store_true")
    ap.add_argument("--ensembl_sleep", type=float, default=0.4)
    ap.add_argument("--ensembl_timeout", type=int, default=30)
    ap.add_argument("--min_protein_len", type=int, default=500)
    ap.add_argument("--max_protein_len", type=int, default=1200)
    ap.add_argument("--strict", action="store_true")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = build_arg_parser()
    args = ap.parse_args(argv)
    base = args.report.parent
    if args.warnings is None:
        args.warnings = base / "protein_export_warnings.tsv"
    if args.region_qc is None:
        args.region_qc = base / "fgfr2_exported_protein_region_qc.tsv"
    if args.metadata is None:
        args.metadata = base / "run_metadata.json"
    if args.md_report is None:
        args.md_report = base / "protein_export_report.md"
    if args.html_report is None:
        args.html_report = base / "protein_export_report.html"
    try:
        report_rows, _warnings, metadata = export_selected(args)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2
    counts: Dict[str, int] = {}
    methods: Dict[str, int] = {}
    for r in report_rows:
        counts[r["match_status"]] = counts.get(r["match_status"], 0) + 1
        methods[r["match_method"]] = methods.get(r["match_method"], 0) + 1
    print(f"Written FASTA records: {metadata['fasta_records_written']}")
    print(f"Output FASTA: {args.out}")
    print(f"Report TSV: {args.report}")
    print(f"Warnings TSV: {args.warnings}")
    print("Match status counts:")
    for k, v in sorted(counts.items()):
        print(f"  {k}: {v}")
    print("Match method counts:")
    for k, v in sorted(methods.items()):
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
