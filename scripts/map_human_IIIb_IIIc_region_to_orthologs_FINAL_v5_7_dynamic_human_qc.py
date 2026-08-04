#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Iterable
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------- data models -----------------------------

@dataclass
class FastaRecord:
    header: str
    seq: str
    species: str
    role: str
    isoform: str
    transcript: str
    protein: str
    source: str

@dataclass
class LocalAlignment:
    score: int
    q_start: int
    q_end: int
    t_start: int
    t_end: int
    q_aln: str
    t_aln: str
    identity: float
    coverage_q: float
    coverage_t: float
    mismatches: int
    gaps: int


# ----------------------------- FASTA / parsing -----------------------------

def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = None
    seq: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def field(header: str, key: str) -> str:
    for p in str(header).split("|"):
        if p.startswith(key + "="):
            return p.split("=", 1)[1]
    return "unknown"


def parse_fasta_records(path: Path) -> List[FastaRecord]:
    out: List[FastaRecord] = []
    for h, s in read_fasta(path):
        out.append(FastaRecord(
            header=h,
            seq=s,
            species=field(h, "species"),
            role=field(h, "role"),
            isoform=field(h, "isoform"),
            transcript=field(h, "transcript"),
            protein=field(h, "protein"),
            source=field(h, "source"),
        ))
    return out


def norm_iso(x: object) -> str:
    s = str(x or "").strip().lower()
    if s in {"iiib", "fgfr2_iiib", "iii_b", "b"}:
        return "IIIb"
    if s in {"iiic", "fgfr2_iiic", "iii_c", "c"}:
        return "IIIc"
    if "iiib" in s:
        return "IIIb"
    if "iiic" in s:
        return "IIIc"
    return "unclassified"


def safe_float(x, default=float("nan")) -> float:
    try:
        if x is None or str(x).strip() == "" or str(x).lower() == "nan":
            return default
        return float(x)
    except Exception:
        return default


def safe_int(x, default=None):
    try:
        if x is None or str(x).strip() == "" or str(x).lower() == "nan":
            return default
        return int(float(str(x)))
    except Exception:
        return default


# ----------------------------- alignment helpers -----------------------------

def smith_waterman(q: str, t: str, match: int = 3, mismatch: int = -2, gap: int = -3) -> LocalAlignment:
    n, m = len(q), len(t)
    if n == 0 or m == 0:
        return LocalAlignment(0, 0, 0, 0, 0, "", "", 0.0, 0.0, 0.0, 0, 0)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    P = [[0] * (m + 1) for _ in range(n + 1)]
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = H[i-1][j-1] + (match if q[i-1] == t[j-1] else mismatch)
            up = H[i-1][j] + gap
            left = H[i][j-1] + gap
            val = max(0, diag, up, left)
            H[i][j] = val
            if val == 0:
                P[i][j] = 0
            elif val == diag:
                P[i][j] = 1
            elif val == up:
                P[i][j] = 2
            else:
                P[i][j] = 3
            if val > best:
                best, bi, bj = val, i, j
    i, j = bi, bj
    qa, ta = [], []
    while i > 0 and j > 0 and H[i][j] > 0:
        p = P[i][j]
        if p == 1:
            qa.append(q[i-1]); ta.append(t[j-1]); i -= 1; j -= 1
        elif p == 2:
            qa.append(q[i-1]); ta.append("-"); i -= 1
        elif p == 3:
            qa.append("-"); ta.append(t[j-1]); j -= 1
        else:
            break
    q_start, t_start = i, j
    q_end, t_end = bi, bj
    qa = "".join(reversed(qa)); ta = "".join(reversed(ta))
    pairs = [(a, b) for a, b in zip(qa, ta) if a != "-" and b != "-"]
    matches = sum(1 for a, b in pairs if a == b)
    mism = sum(1 for a, b in pairs if a != b)
    gaps = sum(1 for a, b in zip(qa, ta) if a == "-" or b == "-")
    identity = matches / len(pairs) if pairs else 0.0
    cov_q = (q_end - q_start) / n if n else 0.0
    cov_t = (t_end - t_start) / m if m else 0.0
    return LocalAlignment(best, q_start, q_end, t_start, t_end, qa, ta, identity, cov_q, cov_t, mism, gaps)


def fixed_identity(a: str, b: str) -> Tuple[float, int, int]:
    n = min(len(a), len(b))
    if n == 0:
        return 0.0, 0, 0
    matches = sum(1 for i in range(n) if a[i] == b[i])
    mism = n - matches + abs(len(a) - len(b))
    return matches / n, matches, mism


def classify_iii_region_similarity(fixed_id: float, local: "LocalAlignment", wb: str, wc: str,
                                   identity_threshold: float,
                                   full_window_id: float = 0.90,
                                   high_local_coverage: float = 0.80,
                                   ambiguous_local_coverage: float = 0.50):
    if not wb or not wc:
        return "unresolved", float("nan"), float("nan"), 0, "missing_window_for_one_isoform"
    cov_min = min(local.coverage_q, local.coverage_t)
    cov_max = max(local.coverage_q, local.coverage_t)
    local_len = len(local.q_aln)
    local_id = local.identity
    if fixed_id >= full_window_id:
        return ("full_window_nearly_identical", cov_min, cov_max, local_len,
                "IIIb_IIIc_windows_positionally_near_identical_isoforms_not_distinguishable")
    if local_id >= identity_threshold:
        if cov_min >= high_local_coverage:
            return ("window_distinct_but_local_subregion_identical", cov_min, cov_max, local_len,
                    "long_identical_subregion_spans_most_of_both_windows_isoforms_poorly_distinguishable")
        if cov_min >= ambiguous_local_coverage:
            return ("ambiguous_similarity_review", cov_min, cov_max, local_len,
                    "substantial_local_identity_below_full_coverage_manual_review_recommended")
        return ("window_distinct_but_short_local_identity", cov_min, cov_max, local_len,
                "only_a_short_local_subregion_is_identical_overall_windows_distinct")
    return ("full_window_distinct", cov_min, cov_max, local_len, "")


def similarity_to_window(candidate_window: str, reference_window: str) -> float:
    n = min(len(candidate_window), len(reference_window))
    if n == 0:
        return 0.0
    matches = sum(1 for i in range(n) if candidate_window[i] == reference_window[i])
    denom = max(len(candidate_window), len(reference_window))
    return matches / denom if denom else 0.0


def unique_kmer_support(query: str, ref: str, other_ref: str, k: int = 8) -> float:
    q = (query or "").upper()
    r = (ref or "").upper()
    o = (other_ref or "").upper()
    if not q or not r or len(r) < k:
        return 0.0
    ref_kmers = {r[i:i+k] for i in range(0, len(r)-k+1)}
    other_kmers = {o[i:i+k] for i in range(0, max(0, len(o)-k+1))}
    unique = sorted(x for x in ref_kmers if x not in other_kmers and "X" not in x)
    if not unique:
        return 0.0
    return sum(1 for x in unique if x in q) / len(unique)


def reference_preference_by_local_alignment(candidate_window: str, seg_b: str, seg_c: str) -> Tuple[str, float, float, float, float]:
    aln_b = smith_waterman(seg_b, candidate_window)
    aln_c = smith_waterman(seg_c, candidate_window)
    marker_b = unique_kmer_support(candidate_window, seg_b, seg_c, k=8)
    marker_c = unique_kmer_support(candidate_window, seg_c, seg_b, k=8)
    support_b = marker_b + 0.05 * (aln_b.identity * aln_b.coverage_q)
    support_c = marker_c + 0.05 * (aln_c.identity * aln_c.coverage_q)
    cov_b = marker_b if marker_b > 0 else aln_b.coverage_q
    cov_c = marker_c if marker_c > 0 else aln_c.coverage_q
    label = "IIIb" if support_b >= support_c else "IIIc"
    return label, support_b, support_c, cov_b, cov_c


# ----------------------------- evidence loading -----------------------------

def read_tsv_optional(path: Optional[Path]) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def get_col(df: pd.DataFrame, aliases: Iterable[str]) -> Optional[str]:
    for a in aliases:
        if a in df.columns:
            return a
    return None


def norm_key(x: object) -> str:
    x = str(x or "").strip()
    if not x or x.lower() in {"nan", "none", "unknown"}:
        return ""
    # Remove FASTA/TSV whitespace payload and common version suffixes.
    x = x.split()[0]
    x = re.sub(r"\.\d+$", "", x)  # XM_... .1, XP_... .2, ENS... .5
    return x


def source_norm(x: object) -> str:
    x = str(x or "").strip().lower()
    if "ensembl" in x:
        return "ensembl"
    if "ncbi" in x or "refseq" in x:
        return "ncbi"
    return x


def species_norm(x: object) -> str:
    return str(x or "").strip().lower().replace(" ", "_")


class EvidenceIndex:
    def __init__(self, evidence: pd.DataFrame):
        self.evidence = evidence.copy() if evidence is not None else pd.DataFrame()
        self.by_species_tx: Dict[Tuple[str, str], dict] = {}
        self.by_tx: Dict[str, List[dict]] = defaultdict(list)
        self.by_species_source_iso: Dict[Tuple[str, str, str], List[dict]] = defaultdict(list)
        self.loaded = not self.evidence.empty
        if self.evidence.empty:
            return

        sp_col = get_col(self.evidence, ["species_canonical", "species", "species_input"])
        src_col = get_col(self.evidence, ["source_db", "source", "database"])
        tx_cols = [c for c in [
            get_col(self.evidence, ["transcript_id_source", "transcript_id", "accession"]),
            get_col(self.evidence, ["internal_transcript_id", "transcript_id_internal"]),
            get_col(self.evidence, ["transcript_name", "name"]),
        ] if c]

        for _, r in self.evidence.iterrows():
            row = dict(r)
            sp = species_norm(r.get(sp_col, "")) if sp_col else ""
            src = source_norm(r.get(src_col, "")) if src_col else ""
            iso = norm_iso(row.get("iii_isoform_assignment") or row.get("isoform_class"))
            row["__species_norm"] = sp
            row["__source_norm"] = src
            row["__iso_norm"] = iso
            # Add normalized transcript IDs for auditability.
            tx_norms = []
            for c in tx_cols:
                val = norm_key(r.get(c, ""))
                if val:
                    tx_norms.append(val)
                    self.by_species_tx[(sp, val)] = row
                    self.by_tx[val].append(row)
            row["__tx_norms"] = ";".join(sorted(set(tx_norms)))
            if sp and src and iso in {"IIIb", "IIIc", "unclassified"}:
                self.by_species_source_iso[(sp, src, iso)].append(row)

    def lookup(self, record: "FastaRecord") -> Tuple[dict, str, str]:
        if not self.loaded:
            return {}, "not_available", "isoform_evidence_not_supplied_or_empty"

        sp = species_norm(record.species)
        src = source_norm(record.source)
        tx = norm_key(record.transcript)
        expected = norm_iso(record.isoform)

        # 1) Best: same species + normalized transcript ID.
        if sp and tx and (sp, tx) in self.by_species_tx:
            return self.by_species_tx[(sp, tx)], "species_transcript_normalized", "matched_by_species_and_transcript_id"

        # 2) Transcript ID globally unique in evidence.
        if tx and tx in self.by_tx:
            rows = self.by_tx[tx]
            if len(rows) == 1:
                return rows[0], "transcript_normalized_unique", "matched_by_unique_transcript_id_without_species"
            same_sp = [r for r in rows if r.get("__species_norm") == sp]
            if len(same_sp) == 1:
                return same_sp[0], "transcript_normalized_species_resolved", "matched_by_transcript_id_and_resolved_species"
            return {}, "not_available", "transcript_id_ambiguous_in_isoform_evidence"

        # 3) Conservative fallback: species + source + expected isoform is unique.
        # This is useful for reference/candidate duplicates where transcript IDs were
        # rewritten upstream but only one transcript of that isoform exists in evidence.
        if expected in {"IIIb", "IIIc"}:
            rows = self.by_species_source_iso.get((sp, src, expected), [])
            informative = [r for r in rows if norm_iso(r.get("iii_isoform_assignment") or r.get("isoform_class")) == expected]
            if len(informative) == 1:
                return informative[0], "species_source_isoform_unique", "matched_by_unique_species_source_isoform_fallback"
            if len(informative) > 1:
                return {}, "not_available", "multiple_isoform_evidence_rows_for_species_source_isoform"

        # 4) Species-only isoform fallback, only if unique and source was unavailable/mismatched.
        if expected in {"IIIb", "IIIc"}:
            candidates = []
            for (k_sp, _k_src, k_iso), rows in self.by_species_source_iso.items():
                if k_sp == sp and k_iso == expected:
                    candidates.extend(rows)
            if len(candidates) == 1:
                return candidates[0], "species_isoform_unique", "matched_by_unique_species_isoform_fallback"
            if len(candidates) > 1:
                return {}, "not_available", "multiple_isoform_evidence_rows_for_species_isoform"

        return {}, "not_available", "transcript_not_found_in_isoform_evidence"


def make_evidence_lookup(evidence: pd.DataFrame) -> EvidenceIndex:
    return EvidenceIndex(evidence)

def exon_signature(row: pd.Series) -> str:
    chrom = str(row.get("chrom", ""))
    start = str(row.get("start", ""))
    end = str(row.get("end", ""))
    strand = str(row.get("strand", ""))
    return f"{chrom}:{start}-{end}:{strand}"




def normalize_exons(exons: pd.DataFrame) -> pd.DataFrame:
    if exons.empty:
        return exons
    df = exons.copy()
    for c in ["start", "end", "exon_rank"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "exon_sig" not in df.columns:
        df["exon_sig"] = df.apply(exon_signature, axis=1)
    return df


# ----------------------------- core analysis -----------------------------

def choose_candidate(records: List[FastaRecord], species: str, isoform: str, ev_index: Optional[EvidenceIndex] = None) -> Optional[FastaRecord]:
    cands = [r for r in records if r.species == species and norm_iso(r.isoform) == isoform]
    if not cands:
        return None

    def candidate_score(r: FastaRecord):
        exon_score = 0
        match_level = "not_available"
        exon_call = "not_available"
        if ev_index is not None and getattr(ev_index, "loaded", False):
            ev, level, _reasonn = ev_index.lookup(r)
            match_level = level
            exon_call = norm_iso(ev.get("iii_isoform_assignment") or ev.get("isoform_class")) if ev else "not_available"
            if exon_call == isoform and level in {"species_transcript_normalized", "transcript_global_normalized"}:
                exon_score = 1000
            elif exon_call == isoform:
                exon_score = 700
            elif exon_call in {"unclassified", "not_available"}:
                exon_score = 0
            else:
                exon_score = -1000
        role_score = 100 if "candidate" in r.role.lower() else 0
        # prefer non-unknown IDs for auditability
        id_score = (20 if r.transcript != "unknown" else 0) + (20 if r.protein != "unknown" else 0)
        return (exon_score, role_score, id_score, len(r.seq), r.transcript, r.protein, match_level, exon_call)

    return sorted(cands, key=candidate_score, reverse=True)[0]


def human_candidate_protein_qc(record: FastaRecord, expected: str, seg_b: str, seg_c: str, fallback_w0: int, fallback_w1: int) -> dict:
    window, w0, w1, mode, anch = extract_dynamic_window(record.seq, seg_b, seg_c, fallback_w0, fallback_w1)
    pref, to_b, to_c, cov_b, cov_c = reference_preference_by_local_alignment(window, seg_b, seg_c)
    expected = norm_iso(expected)
    expected_support = to_b if expected == "IIIb" else to_c
    opposite_support = to_c if expected == "IIIb" else to_b
    expected_cov = cov_b if expected == "IIIb" else cov_c
    full_status = str(anch.get("full_region_status", ""))
    best_segment = str(anch.get("human_full_best_segment", ""))
    delta = expected_support - opposite_support
    passed = (
        expected in {"IIIb", "IIIc"}
        and pref == expected
        and full_status == "full_region_detected"
        and best_segment == expected
        and expected_support >= 0.25
        and expected_cov >= 0.25
        and delta > 0.05
    )
    if passed:
        reason = "pass"
    elif full_status != "full_region_detected":
        reason = f"fail_{expected}_candidate_full_region_anchor_weak_or_missing"
    elif best_segment != expected:
        reason = f"fail_{expected}_candidate_anchor_best_segment_is_{best_segment}"
    elif pref != expected:
        reason = f"fail_{expected}_candidate_window_prefers_{pref}"
    elif expected_support < 0.25 or expected_cov < 0.25:
        reason = f"fail_{expected}_candidate_low_curated_reference_support"
    else:
        reason = f"fail_{expected}_candidate_insufficient_margin"
    return {
        "window": window, "w0": w0, "w1": w1, "mode": mode, "anchor": anch,
        "prefers": pref, "support_to_IIIb_ref": to_b, "support_to_IIIc_ref": to_c,
        "coverage_IIIb_ref": cov_b, "coverage_IIIc_ref": cov_c,
        "expected_support": expected_support, "opposite_support": opposite_support,
        "ref_delta_expected_minus_opposite": delta,
        "full_region_status": full_status, "anchor_best_segment": best_segment,
        "status": "pass" if passed else "fail", "failure_reason": reason,
    }


def anchor_human_region(seq: str, seg_b: str, seg_c: str, plausible_min: int, plausible_max: int) -> dict:
    aln_b = smith_waterman(seg_b, seq)
    aln_c = smith_waterman(seg_c, seq)
    best_label, best = ("IIIb", aln_b) if aln_b.score >= aln_c.score else ("IIIc", aln_c)
    midpoint = (best.t_start + best.t_end) / 2.0 + 1
    if best.coverage_q >= 0.70 and best.identity >= 0.35:
        full_status = "full_region_detected"
    else:
        full_status = "full_region_low_confidence_or_missing"
    if best.t_end == 0:
        pos_status = "position_missing"
    elif plausible_min <= midpoint <= plausible_max:
        pos_status = "position_plausible"
    elif plausible_min - 50 <= midpoint <= plausible_max + 50:
        pos_status = "position_shifted_warning"
    else:
        pos_status = "position_outlier"
    return {
        "human_full_best_segment": best_label,
        "human_full_score": best.score,
        "human_full_identity": best.identity,
        "human_full_coverage": best.coverage_q,
        "human_full_anchor_start": best.t_start + 1 if best.t_end else "",
        "human_full_anchor_end": best.t_end if best.t_end else "",
        "human_full_anchor_midpoint": midpoint if best.t_end else float("nan"),
        "full_region_status": full_status,
        "anchor_position_status": pos_status,
    }




def dynamic_region_bounds(seq: str, seg_b: str, seg_c: str, fallback_w0: int, fallback_w1: int, flank: int = 12,
                          min_identity: float = 0.35, min_coverage: float = 0.55) -> Tuple[int, int, str, dict]:
    anch = anchor_human_region(seq, seg_b, seg_c, plausible_min=1, plausible_max=max(len(seq), 1))
    start = safe_int(anch.get("human_full_anchor_start"), None)
    end = safe_int(anch.get("human_full_anchor_end"), None)
    ident = safe_float(anch.get("human_full_identity"), 0.0)
    cov = safe_float(anch.get("human_full_coverage"), 0.0)
    if start is not None and end is not None and end > start and ident >= min_identity and cov >= min_coverage:
        w0 = max(0, start - 1 - flank)
        w1 = min(len(seq), end + flank)
        mode = "dynamic_alignment_anchor"
    else:
        w0 = max(0, min(len(seq), fallback_w0))
        w1 = max(w0, min(len(seq), fallback_w1))
        mode = "fallback_fixed_window_low_confidence_anchor"
    return w0, w1, mode, anch


def extract_dynamic_window(seq: str, seg_b: str, seg_c: str, fallback_w0: int, fallback_w1: int) -> Tuple[str, int, int, str, dict]:
    w0, w1, mode, anch = dynamic_region_bounds(seq, seg_b, seg_c, fallback_w0, fallback_w1)
    return seq[w0:w1], w0, w1, mode, anch

def build_pair_audit(records: List[FastaRecord], seg_b: str, seg_c: str, w0: int, w1: int, identity_threshold: float, diff_threshold: int, ev_index: Optional[EvidenceIndex] = None) -> Tuple[pd.DataFrame, Dict[str, dict], List[str], List[dict]]:
    species_list = sorted({r.species for r in records})
    rows = []
    pair_lookup: Dict[str, dict] = {}
    fasta_lines: List[str] = []
    diff_rows: List[dict] = []
    for sp in species_list:
        rb = choose_candidate(records, sp, "IIIb", ev_index)
        rc = choose_candidate(records, sp, "IIIc", ev_index)
        row = {"species": sp, "has_IIIb": rb is not None, "has_IIIc": rc is not None}
        if rb is None or rc is None:
            row["pair_audit_status"] = "missing_pair_member"
            rows.append(row); pair_lookup[sp] = row
            continue
        wb, b0, b1, bmode, _banch = extract_dynamic_window(rb.seq, seg_b, seg_c, w0, w1)
        wc, c0, c1, cmode, _canch = extract_dynamic_window(rc.seq, seg_b, seg_c, w0, w1)
        fid, _fm, fdiff = fixed_identity(wb, wc)
        local = smith_waterman(wb, wc)
        n_local_diff = local.mismatches + local.gaps
        # Refined similarity is independent of the legacy status.
        similarity_class, local_cov_min, local_cov_max, local_len, similarity_warning = classify_iii_region_similarity(
            fid, local, wb, wc, identity_threshold)
        if rb.protein == rc.protein and rb.protein != "unknown":
            status = "same_protein_accession_warning"
        elif rb.transcript == rc.transcript and rb.transcript != "unknown":
            status = "same_transcript_accession_warning"
        elif rb.seq == rc.seq:
            status = "full_proteins_identical"
        elif local.identity >= identity_threshold and n_local_diff <= diff_threshold:
            status = "III_region_nearly_identical"
        else:
            status = "III_region_sequence_distinct"
        # Preserve the exact legacy label before any human-control override below.
        legacy_pair_audit_status = status

        # Built-in positive control: in Homo sapiens, the selected IIIb candidate
        # should be more similar to the curated human IIIb segment, and IIIc should
        # be more similar to the curated human IIIc segment. If this fails, the
        # upstream exon/selection layer probably chose the wrong alternative slot.
        human_control_status = "not_applicable"
        human_IIIb_ref_delta = float("nan")
        human_IIIc_ref_delta = float("nan")
        human_IIIb_validation_reason = ""
        human_IIIc_validation_reason = ""
        human_IIIb_anchor_status = ""
        human_IIIc_anchor_status = ""
        human_IIIb_anchor_best_segment = ""
        human_IIIc_anchor_best_segment = ""
        if sp.lower() in {"homo_sapiens", "human"}:
            bq = human_candidate_protein_qc(rb, "IIIb", seg_b, seg_c, w0, w1)
            cq = human_candidate_protein_qc(rc, "IIIc", seg_b, seg_c, w0, w1)
            b_pref = bq["prefers"]; c_pref = cq["prefers"]
            b_to_b = bq["support_to_IIIb_ref"]; b_to_c = bq["support_to_IIIc_ref"]
            c_to_b = cq["support_to_IIIb_ref"]; c_to_c = cq["support_to_IIIc_ref"]
            human_IIIb_ref_delta = bq["ref_delta_expected_minus_opposite"]
            human_IIIc_ref_delta = cq["ref_delta_expected_minus_opposite"]
            human_IIIb_validation_reason = bq["failure_reason"]
            human_IIIc_validation_reason = cq["failure_reason"]
            human_IIIb_anchor_status = bq["full_region_status"]
            human_IIIc_anchor_status = cq["full_region_status"]
            human_IIIb_anchor_best_segment = bq["anchor_best_segment"]
            human_IIIc_anchor_best_segment = cq["anchor_best_segment"]
            if wb == wc:
                human_control_status = "fail_selected_human_IIIb_IIIc_windows_are_identical"
                status = "human_positive_control_failed_identical_windows"
            elif bq["status"] == "pass" and cq["status"] == "pass":
                human_control_status = "pass"
                if status == "III_region_nearly_identical":
                    status = "III_region_sequence_distinct"
            elif bq["status"] != "pass" and cq["status"] != "pass":
                human_control_status = "fail_selected_human_IIIb_and_IIIc_candidates_not_protein_validated"
                status = "human_positive_control_failed_both_candidates_not_protein_validated"
            elif cq["status"] != "pass":
                human_control_status = "fail_selected_human_IIIc_candidate_not_protein_IIIc_like"
                status = "human_positive_control_failed_selected_IIIc_not_protein_IIIc_like"
            elif bq["status"] != "pass":
                human_control_status = "fail_selected_human_IIIb_candidate_not_protein_IIIb_like"
                status = "human_positive_control_failed_selected_IIIb_not_protein_IIIb_like"
        row.update({
            "IIIb_transcript": rb.transcript, "IIIc_transcript": rc.transcript,
            "IIIb_protein": rb.protein, "IIIc_protein": rc.protein,
            "IIIb_length": len(rb.seq), "IIIc_length": len(rc.seq),
            "length_difference_IIIc_minus_IIIb": len(rc.seq) - len(rb.seq),
            "IIIb_window": wb, "IIIc_window": wc,
            "IIIb_window_start_1based": b0 + 1 if wb else "",
            "IIIb_window_end_1based": b1 if wb else "",
            "IIIc_window_start_1based": c0 + 1 if wc else "",
            "IIIc_window_end_1based": c1 if wc else "",
            "IIIb_window_mode": bmode, "IIIc_window_mode": cmode,
            "human_control_status": human_control_status,
            "human_IIIb_ref_delta": human_IIIb_ref_delta,
            "human_IIIc_ref_delta": human_IIIc_ref_delta,
            "human_IIIb_validation_reason": human_IIIb_validation_reason,
            "human_IIIc_validation_reason": human_IIIc_validation_reason,
            "human_IIIb_anchor_status": human_IIIb_anchor_status,
            "human_IIIc_anchor_status": human_IIIc_anchor_status,
            "human_IIIb_anchor_best_segment": human_IIIb_anchor_best_segment,
            "human_IIIc_anchor_best_segment": human_IIIc_anchor_best_segment,
            "human_IIIb_prefers": b_pref if sp.lower() in {"homo_sapiens", "human"} else "",
            "human_IIIc_prefers": c_pref if sp.lower() in {"homo_sapiens", "human"} else "",
            "human_IIIb_support_to_IIIb_ref": b_to_b if sp.lower() in {"homo_sapiens", "human"} else "",
            "human_IIIb_support_to_IIIc_ref": b_to_c if sp.lower() in {"homo_sapiens", "human"} else "",
            "human_IIIc_support_to_IIIb_ref": c_to_b if sp.lower() in {"homo_sapiens", "human"} else "",
            "human_IIIc_support_to_IIIc_ref": c_to_c if sp.lower() in {"homo_sapiens", "human"} else "",
            "fixed_window_identity": fid,
            "fixed_window_mismatches_plus_length_delta": fdiff,
            "regional_local_identity": local.identity,
            "regional_local_mismatches_plus_gaps": n_local_diff,
            "local_alignment_coverage_min": local_cov_min,
            "local_alignment_coverage_max": local_cov_max,
            "local_alignment_length_aa": local_len,
            "iii_region_similarity_class": similarity_class,
            "iii_region_similarity_warning": similarity_warning,
            "legacy_pair_audit_status": legacy_pair_audit_status,
            "pair_audit_status": status,
        })
        rows.append(row); pair_lookup[sp] = row
        fasta_lines += [
            f">{sp}|isoform=IIIb|window={b0+1}-{b1}|mode={bmode}|transcript={rb.transcript}|protein={rb.protein}", wb,
            f">{sp}|isoform=IIIc|window={c0+1}-{c1}|mode={cmode}|transcript={rc.transcript}|protein={rc.protein}", wc,
        ]
        n = min(len(wb), len(wc))
        for i in range(n):
            if wb[i] != wc[i]:
                diff_rows.append({"species": sp, "relative_window_position": i+1, "IIIb_protein_position_1based": b0+i+1, "IIIc_protein_position_1based": c0+i+1, "IIIb_aa": wb[i], "IIIc_aa": wc[i], "difference_type": "substitution"})
        if len(wb) != len(wc):
            for i in range(n, max(len(wb), len(wc))):
                diff_rows.append({"species": sp, "relative_window_position": i+1, "IIIb_protein_position_1based": b0+i+1 if i < len(wb) else "", "IIIc_protein_position_1based": c0+i+1 if i < len(wc) else "", "IIIb_aa": wb[i] if i < len(wb) else "-", "IIIc_aa": wc[i] if i < len(wc) else "-", "difference_type": "gap_or_length_difference"})
    return pd.DataFrame(rows), pair_lookup, fasta_lines, diff_rows


def species_window_score(record: FastaRecord, pair: dict, seg_b: str, seg_c: str, w0: int, w1: int, delta_threshold: float) -> dict:
    status = pair.get("pair_audit_status", "missing_pair_member")
    if status != "III_region_sequence_distinct":
        return {"species_motif_call": "uninformative", "species_IIIb_support": float("nan"), "species_IIIc_support": float("nan"), "species_motif_delta": float("nan")}
    cand_w, cw0, cw1, cmode, _canch = extract_dynamic_window(record.seq, seg_b, seg_c, w0, w1)
    wb, wc = pair.get("IIIb_window", ""), pair.get("IIIc_window", "")
    sb = similarity_to_window(cand_w, wb)
    sc = similarity_to_window(cand_w, wc)
    delta = sb - sc
    if abs(delta) < delta_threshold:
        call = "ambiguous"
    elif delta > 0:
        call = "IIIb"
    else:
        call = "IIIc"
    return {"species_motif_call": call, "species_IIIb_support": sb, "species_IIIc_support": sc, "species_motif_delta": delta, "candidate_window_start_1based": cw0 + 1 if cand_w else "", "candidate_window_end_1based": cw1 if cand_w else "", "candidate_window_mode": cmode}


def join_exon_evidence(record: FastaRecord, ev_index: EvidenceIndex) -> dict:
    ev, level, reason = ev_index.lookup(record)
    if not ev:
        return {
            "exon_isoform_call": "not_available",
            "exon_evidence_confidence": "",
            "exon_iiib_sig": "",
            "exon_iiic_sig": "",
            "exon_matched_alt_exons": "",
            "exon_slot_id": "",
            "exon_evidence_match_level": level,
            "exon_evidence_missing_reason": reason,
            "exon_evidence_transcript_ids": "",
            "exon_evidence_source": "",
        }
    call = norm_iso(ev.get("iii_isoform_assignment") or ev.get("isoform_class"))
    txs = ev.get("__tx_norms", "")
    return {
        "exon_isoform_call": call,
        "exon_evidence_confidence": ev.get("confidence", ""),
        "exon_iiib_sig": ev.get("iiib_exon_sig", ""),
        "exon_iiic_sig": ev.get("iiic_exon_sig", ""),
        "exon_matched_alt_exons": ev.get("matched_alt_exons", ""),
        "exon_slot_id": ev.get("slot_id", ""),
        "exon_evidence_match_level": level,
        "exon_evidence_missing_reason": reason,
        "exon_evidence_transcript_ids": txs,
        "exon_evidence_source": ev.get("source_db", ""),
    }

def final_status(row: dict) -> str:
    call = row.get("species_motif_call", "")
    exon_call = row.get("exon_isoform_call", "not_available")
    exp = row.get("expected_isoform", "unclassified")
    agreement = row.get("exon_sequence_agreement", "not_available")
    full = row.get("full_region_status", "")
    pos = row.get("anchor_position_status", "")
    pair = row.get("pair_audit_status", "")

    # Highest-confidence integrated evidence: exon and sequence agree.
    if agreement == "agree" and call in {"IIIb", "IIIc"} and exon_call == call:
        if pos == "position_outlier":
            return f"exon_sequence_supported_position_outlier_{call}"
        if full != "full_region_detected":
            return f"exon_sequence_supported_anchor_weak_{call}"
        return f"exon_sequence_supported_{call}"

    # Pair is locally non-informative: do not force IIIb/IIIc from sequence.
    if call == "uninformative" or pair != "III_region_sequence_distinct":
        if exon_call in {"IIIb", "IIIc"}:
            return "sequence_uninformative_exon_supported"
        if full == "full_region_detected":
            return "region_detected_species_pair_uninformative"
        return "region_not_reliably_detected"

    # Position outlier without integrated support remains a region-level warning.
    if pos == "position_outlier":
        return "region_detected_position_outlier"

    # Region-level sequence support without exon evidence.
    if call in {"IIIb", "IIIc"}:
        if exp == call:
            return "region_detected_species_pair_supported"
        return "region_detected_sequence_label_discordance"

    if call == "ambiguous":
        if exon_call in {"IIIb", "IIIc"}:
            return "exon_supported_but_sequence_ambiguous"
        return "region_detected_isoform_ambiguous"

    return "region_not_reliably_detected"


def simplify_final_status(status: str) -> str:
    status = str(status)
    if status.startswith("exon_sequence_supported_anchor_weak"):
        return "exon+sequence supported; human anchor weak"
    if status.startswith("exon_sequence_supported_position_outlier"):
        return "exon+sequence supported; position outlier"
    if status.startswith("exon_sequence_supported"):
        return "exon+sequence supported"
    if status == "sequence_uninformative_exon_supported":
        return "exon supported; sequence uninformative"
    if "position_outlier" in status:
        return "position outlier"
    if "ambiguous" in status:
        return "sequence ambiguous"
    if "uninformative" in status:
        return "sequence uninformative"
    if "not_reliably" in status:
        return "region not reliably detected"
    if "discordance" in status:
        return "discordance"
    return status


# ----------------------------- plots -----------------------------

def savefig(outdir: Path, name: str):
    plt.tight_layout()
    plt.savefig(outdir / f"{name}.png", dpi=300)
    plt.savefig(outdir / f"{name}.pdf")
    plt.close()


def plot_counts(series: pd.Series, outdir: Path, name: str, title: str, xlabel: str):
    counts = series.value_counts().sort_values(ascending=True)
    plt.figure(figsize=(8, max(3.5, 0.45 * len(counts))))
    plt.barh(counts.index, counts.values)
    for i, v in enumerate(counts.values):
        plt.text(v + 0.15, i, str(v), va="center")
    plt.xlabel(xlabel); plt.ylabel(""); plt.title(title)
    savefig(outdir, name)


def ensure_pair_audit_plot_columns(audit: pd.DataFrame) -> pd.DataFrame:
    if audit is None or audit.empty:
        return audit
    audit = audit.copy()
    if "regional_local_identity" not in audit.columns:
        fallback_cols = [
            "local_identity",
            "III_region_local_identity",
            "dynamic_window_identity",
            "fixed_window_identity",
        ]
        for col in fallback_cols:
            if col in audit.columns:
                audit["regional_local_identity"] = pd.to_numeric(audit[col], errors="coerce")
                break
        else:
            audit["regional_local_identity"] = float("nan")
    else:
        audit["regional_local_identity"] = pd.to_numeric(audit["regional_local_identity"], errors="coerce")

    if "regional_local_mismatches_plus_gaps" not in audit.columns:
        fallback_cols = [
            "local_mismatches_plus_gaps",
            "III_region_mismatches_plus_gaps",
            "fixed_window_mismatches_plus_length_delta",
        ]
        for col in fallback_cols:
            if col in audit.columns:
                audit["regional_local_mismatches_plus_gaps"] = pd.to_numeric(audit[col], errors="coerce")
                break
        else:
            audit["regional_local_mismatches_plus_gaps"] = float("nan")
    return audit


def plot_pair_identity(audit: pd.DataFrame, outdir: Path, prefix: str, threshold: float):
    audit = ensure_pair_audit_plot_columns(audit)
    if audit is None or audit.empty:
        return
    df = audit.copy().sort_values("regional_local_identity", ascending=True, na_position="first")
    plt.figure(figsize=(8, max(5, 0.28 * len(df))))
    plt.barh(df["species"], pd.to_numeric(df["regional_local_identity"], errors="coerce"))
    plt.axvline(threshold, linestyle="--", linewidth=1)
    plt.xlabel("Local identity between exported IIIb and IIIc candidate windows")
    plt.ylabel("Species")
    plt.title("Pair-level IIIb/IIIc distinguishability")
    savefig(outdir, f"{prefix}_v5_pair_identity_barplot")


def plot_matrix(df: pd.DataFrame, row_col: str, col_col: str, outdir: Path, name: str, title: str):
    if df.empty or row_col not in df.columns or col_col not in df.columns:
        return
    tab = pd.crosstab(df[row_col], df[col_col])
    if tab.empty:
        return
    plt.figure(figsize=(1.4 * len(tab.columns) + 3.5, 1.0 * len(tab.index) + 3))
    im = plt.imshow(tab.values, aspect="auto", interpolation="nearest")
    plt.xticks(range(len(tab.columns)), tab.columns, rotation=35, ha="right")
    plt.yticks(range(len(tab.index)), tab.index)
    for i in range(tab.shape[0]):
        for j in range(tab.shape[1]):
            plt.text(j, i, str(tab.iloc[i, j]), ha="center", va="center")
    plt.colorbar(im, label="Number of proteins")
    plt.xlabel(col_col.replace("_", " "))
    plt.ylabel(row_col.replace("_", " "))
    plt.title(title)
    savefig(outdir, name)


def plot_support_scatter(df: pd.DataFrame, outdir: Path, prefix: str):
    if df.empty or "species_IIIb_support" not in df or "species_IIIc_support" not in df:
        return
    plt.figure(figsize=(7, 6))
    for status in sorted(df["final_anchor_status"].dropna().unique()):
        sub = df[df["final_anchor_status"] == status]
        plt.scatter(pd.to_numeric(sub["species_IIIb_support"], errors="coerce"), pd.to_numeric(sub["species_IIIc_support"], errors="coerce"), label=status, alpha=0.8)
    plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1)
    plt.xlabel("Species-calibrated IIIb window support")
    plt.ylabel("Species-calibrated IIIc window support")
    plt.title("Species-specific IIIb vs IIIc sequence support")
    plt.legend(fontsize=7, loc="best")
    savefig(outdir, f"{prefix}_v5_species_support_scatter")


def plot_difference_map(diff_df: pd.DataFrame, audit: pd.DataFrame, outdir: Path, prefix: str, w0: int, w1: int):
    audit = ensure_pair_audit_plot_columns(audit)
    if audit is None or audit.empty or "species" not in audit.columns:
        return
    species = audit.sort_values("regional_local_identity", ascending=True, na_position="first")["species"].tolist()
    width = w1 - w0
    mat = []
    diff_lookup = defaultdict(dict)
    for _, r in diff_df.iterrows():
        diff_lookup[r["species"]][int(r["relative_window_position"]) - 1] = 2 if r.get("difference_type") == "gap_or_length_difference" else 1
    for sp in species:
        row = [0] * width
        for pos, val in diff_lookup.get(sp, {}).items():
            if 0 <= pos < width:
                row[pos] = val
        mat.append(row)
    plt.figure(figsize=(12, max(5, 0.26 * len(species))))
    im = plt.imshow(mat, aspect="auto", interpolation="nearest", vmin=0, vmax=2)
    plt.yticks(range(len(species)), species, fontsize=6)
    xticks = list(range(0, width, max(1, width // 10)))
    plt.xticks(xticks, [str(w0 + x + 1) for x in xticks], rotation=0)
    plt.xlabel("Aligned relative position in dynamically anchored III-region window")
    plt.ylabel("Species")
    plt.title("IIIb/IIIc candidate-pair difference map")
    cbar = plt.colorbar(im, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(["same", "substitution", "gap/length"])
    savefig(outdir, f"{prefix}_v5_pair_difference_map")


def plot_species_stacked(df: pd.DataFrame, outdir: Path, prefix: str):
    if df.empty:
        return
    tab = pd.crosstab(df["species"], df["final_anchor_status"])
    tab = tab.loc[sorted(tab.index)]
    plt.figure(figsize=(10, max(6, 0.28 * len(tab))))
    left = [0] * len(tab)
    y = range(len(tab))
    for col in tab.columns:
        vals = tab[col].values
        plt.barh(y, vals, left=left, label=col)
        left = [a + b for a, b in zip(left, vals)]
    plt.yticks(y, tab.index)
    plt.xlabel("Number of candidate proteins")
    plt.ylabel("Species")
    plt.title("Species-level final III-region evidence composition")
    plt.legend(fontsize=7, loc="best")
    savefig(outdir, f"{prefix}_v5_species_status_stacked")


def plot_exon_architecture(records: List[FastaRecord], exons: pd.DataFrame, ev_lookup: EvidenceIndex, outdir: Path, prefix: str, max_species: int):
    if exons.empty or not getattr(ev_lookup, "loaded", False):
        return
    tx_col = get_col(exons, ["transcript_id_internal", "internal_transcript_id", "tx_internal_id"])
    rank_col = get_col(exons, ["exon_rank", "rank"])
    start_col = get_col(exons, ["start"])
    end_col = get_col(exons, ["end"])
    _sp_col = get_col(exons, ["species_canonical", "species", "species_input"])
    if not all([tx_col, rank_col, start_col, end_col]):
        return
    ex = normalize_exons(exons)
    chosen_species = sorted({r.species for r in records})[:max_species]
    rows = []
    for sp in chosen_species:
        for iso in ["IIIb", "IIIc"]:
            rec = choose_candidate(records, sp, iso)
            if rec:
                rows.append((sp, iso, rec))
    if not rows:
        return
    fig_h = max(6, 0.45 * len(rows))
    plt.figure(figsize=(13, fig_h))
    ylabels = []
    for y, (sp, iso, rec) in enumerate(rows):
        ev, _match_level, _missing_reason = ev_lookup.lookup(rec)
        # Find exons by internal transcript if possible. If evidence has internal id, use it.
        internal_id = ev.get("internal_transcript_id", "") if ev else ""
        sub = ex[ex[tx_col].astype(str) == str(internal_id)] if internal_id else pd.DataFrame()
        if sub.empty:
            # fallback: try source transcript id in transcript_id_internal rarely impossible; skip schematic line with protein only
            ylabels.append(f"{sp} | {iso}")
            plt.text(0.02, y, "exon rows not matched", va="center", fontsize=6)
            continue
        sub = sub.sort_values(rank_col)
        min_s = pd.to_numeric(sub[start_col], errors="coerce").min()
        max_e = pd.to_numeric(sub[end_col], errors="coerce").max()
        span = max(1, max_e - min_s)
        iiib_sig = ev.get("iiib_exon_sig", "")
        iiic_sig = ev.get("iiic_exon_sig", "")
        matched = set(str(ev.get("matched_alt_exons", "")).split(";")) if ev else set()
        for _, er in sub.iterrows():
            s = float(er[start_col]); e = float(er[end_col]); rank = int(er[rank_col]) if not pd.isna(er[rank_col]) else 0
            x = (s - min_s) / span
            w = max(0.004, (e - s + 1) / span)
            sig = er.get("exon_sig", exon_signature(er))
            height = 0.28
            # default rectangles are uncolored by explicit matplotlib defaults; use hatching/edge semantics
            if sig == iiib_sig:
                hatch = "///"; lw = 1.8
            elif sig == iiic_sig:
                hatch = "\\\\\\"; lw = 1.8
            elif sig in matched:
                hatch = "xx"; lw = 1.4
            else:
                hatch = ""; lw = 0.8
            plt.gca().add_patch(plt.Rectangle((x, y - height/2), w, height, fill=False, linewidth=lw, hatch=hatch))
            if rank % 2 == 1 or sig in {iiib_sig, iiic_sig}:
                plt.text(x + w/2, y + 0.18, str(rank), ha="center", va="bottom", fontsize=5)
        ylabels.append(f"{sp} | {iso} | exon={norm_iso(ev.get('iii_isoform_assignment',''))}")
    plt.ylim(-1, len(rows))
    plt.xlim(-0.02, 1.02)
    plt.yticks(range(len(rows)), ylabels, fontsize=6)
    plt.xlabel("Normalized genomic span of transcript model")
    plt.ylabel("Candidate transcript")
    plt.title("Exon architecture overview with IIIb/IIIc alternative-exon signatures\nHatched rectangles mark inferred IIIb/IIIc slot exons when available")
    savefig(outdir, f"{prefix}_v5_exon_architecture_overview")



def plot_simplified_status_counts(df: pd.DataFrame, outdir: Path, prefix: str):
    if df.empty or "final_anchor_status_simplified" not in df.columns:
        return
    counts = df["final_anchor_status_simplified"].value_counts().sort_values(ascending=True)
    plt.figure(figsize=(8, max(3.5, 0.45 * len(counts))))
    plt.barh(counts.index, counts.values)
    for i, v in enumerate(counts.values):
        plt.text(v + 0.15, i, str(v), va="center")
    plt.xlabel("Number of candidate proteins")
    plt.ylabel("")
    plt.title("Simplified final III-region evidence classes")
    savefig(outdir, f"{prefix}_v5_final_status_simplified_counts")


def plot_evidence_agreement_storyboard(df: pd.DataFrame, audit: pd.DataFrame, outdir: Path, prefix: str):
    _fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # A simplified classes
    counts = df["final_anchor_status_simplified"].value_counts().sort_values(ascending=True)
    axes[0,0].barh(counts.index, counts.values)
    axes[0,0].set_title("A. Integrated evidence classes")
    axes[0,0].set_xlabel("Proteins")

    # B exon vs sequence agreement
    agree_counts = df["exon_sequence_agreement"].value_counts().sort_values(ascending=True)
    axes[0,1].barh(agree_counts.index, agree_counts.values)
    axes[0,1].set_title("B. Exon–sequence agreement")
    axes[0,1].set_xlabel("Proteins")

    # C exon call vs sequence call matrix
    tab = pd.crosstab(df["exon_isoform_call"], df["species_motif_call"])
    _im = axes[1,0].imshow(tab.values, aspect="auto", interpolation="nearest")
    axes[1,0].set_xticks(range(len(tab.columns)))
    axes[1,0].set_xticklabels(tab.columns, rotation=35, ha="right")
    axes[1,0].set_yticks(range(len(tab.index)))
    axes[1,0].set_yticklabels(tab.index)
    for i in range(tab.shape[0]):
        for j in range(tab.shape[1]):
            axes[1,0].text(j, i, str(tab.iloc[i, j]), ha="center", va="center")
    axes[1,0].set_title("C. Exon call vs sequence call")
    axes[1,0].set_xlabel("Species-calibrated sequence call")
    axes[1,0].set_ylabel("Exon isoform call")

    # D pair audit
    pc = audit["pair_audit_status"].value_counts().sort_values(ascending=True)
    axes[1,1].barh(pc.index, pc.values)
    axes[1,1].set_title("D. Species-level pair audit")
    axes[1,1].set_xlabel("Species")

    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_v5_evidence_agreement_storyboard.png", dpi=300)
    plt.savefig(outdir / f"{prefix}_v5_evidence_agreement_storyboard.pdf")
    plt.close()


def parse_species_list(species_arg: str) -> List[str]:
    if not species_arg:
        return []
    return [x.strip() for x in species_arg.split(",") if x.strip()]


def plot_representative_exon_motif_tracks(records: List[FastaRecord], exons: pd.DataFrame, ev_lookup: EvidenceIndex, df: pd.DataFrame, outdir: Path, prefix: str, representative_species: List[str]):
    if exons.empty or not getattr(ev_lookup, "loaded", False) or not representative_species:
        return
    tx_col = get_col(exons, ["transcript_id_internal", "internal_transcript_id", "tx_internal_id"])
    rank_col = get_col(exons, ["exon_rank", "rank"])
    start_col = get_col(exons, ["start"])
    end_col = get_col(exons, ["end"])
    if not all([tx_col, rank_col, start_col, end_col]):
        return
    ex = normalize_exons(exons)
    status_lookup = {(r["species"], r["transcript"]): r for _, r in df.iterrows()}

    rows = []
    for sp in representative_species:
        for iso in ["IIIb", "IIIc"]:
            rec = choose_candidate(records, sp, iso)
            if rec:
                rows.append((sp, iso, rec))
    if not rows:
        return

    fig_h = max(5, 0.7 * len(rows))
    _fig, ax = plt.subplots(figsize=(14, fig_h))
    ylabels = []
    for y, (sp, iso, rec) in enumerate(rows):
        ev, _match_level, _missing_reason = ev_lookup.lookup(rec)
        internal_id = ev.get("internal_transcript_id", "") if ev else ""
        sub = ex[ex[tx_col].astype(str) == str(internal_id)] if internal_id else pd.DataFrame()
        row = status_lookup.get((sp, rec.transcript), {})
        seq_call = row.get("species_motif_call", "NA") if isinstance(row, dict) else row.get("species_motif_call", "NA")
        final = row.get("final_anchor_status_simplified", "NA") if isinstance(row, dict) else row.get("final_anchor_status_simplified", "NA")
        exon_call = row.get("exon_isoform_call", "NA") if isinstance(row, dict) else row.get("exon_isoform_call", "NA")

        if sub.empty:
            ax.text(0.02, y, "exons not matched", va="center", fontsize=8)
            ylabels.append(f"{sp} | {iso}")
            continue
        sub = sub.sort_values(rank_col)
        min_s = pd.to_numeric(sub[start_col], errors="coerce").min()
        max_e = pd.to_numeric(sub[end_col], errors="coerce").max()
        span = max(1, max_e - min_s)
        iiib_sig = ev.get("iiib_exon_sig", "") if ev else ""
        iiic_sig = ev.get("iiic_exon_sig", "") if ev else ""
        alt_ranks = []
        for _, er in sub.iterrows():
            s = float(er[start_col]); e = float(er[end_col])
            rank = int(er[rank_col]) if not pd.isna(er[rank_col]) else 0
            x = (s - min_s) / span
            w = max(0.006, (e - s + 1) / span)
            sig = er.get("exon_sig", exon_signature(er))
            height = 0.34
            hatch = ""
            lw = 0.9
            if sig == iiib_sig:
                hatch = "///"; lw = 2.0; alt_ranks.append(f"IIIb exon {rank}")
            elif sig == iiic_sig:
                hatch = "\\\\\\"; lw = 2.0; alt_ranks.append(f"IIIc exon {rank}")
            ax.add_patch(plt.Rectangle((x, y - height/2), w, height, fill=False, linewidth=lw, hatch=hatch))
            if sig in {iiib_sig, iiic_sig}:
                ax.text(x + w/2, y + 0.24, str(rank), ha="center", va="bottom", fontsize=7)
        # right-side evidence badge as text
        ax.text(1.02, y, f"exon={exon_call}; seq={seq_call}; {final}", va="center", fontsize=7)
        ylabels.append(f"{sp} | {iso}")

    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xlim(-0.02, 1.55)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Normalized transcript span")
    ax.set_ylabel("Representative candidate transcript")
    ax.set_title("Representative exon–motif tracks for FGFR2 IIIb/IIIc candidates\nHatched exons mark inferred IIIb/IIIc alternative-exon signatures; right labels show exon/sequence evidence")
    savefig(outdir, f"{prefix}_v5_representative_exon_motif_tracks")


def plot_storyboard(df: pd.DataFrame, audit: pd.DataFrame, diff_df: pd.DataFrame, outdir: Path, prefix: str):
    # A compact multi-panel summary for thesis figures.
    _fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    # A counts
    counts = df["final_anchor_status"].value_counts().sort_values(ascending=True)
    axes[0,0].barh(counts.index, counts.values)
    axes[0,0].set_title("A. Final III-region evidence status")
    axes[0,0].set_xlabel("Proteins")
    # B pair audit
    pc = audit["pair_audit_status"].value_counts().sort_values(ascending=True)
    axes[0,1].barh(pc.index, pc.values)
    axes[0,1].set_title("B. Species-level pair audit")
    axes[0,1].set_xlabel("Species")
    # C matrix
    tab = pd.crosstab(df["expected_isoform"], df["species_motif_call"])
    _im = axes[1,0].imshow(tab.values, aspect="auto", interpolation="nearest")
    axes[1,0].set_xticks(range(len(tab.columns))); axes[1,0].set_xticklabels(tab.columns, rotation=35, ha="right")
    axes[1,0].set_yticks(range(len(tab.index))); axes[1,0].set_yticklabels(tab.index)
    for i in range(tab.shape[0]):
        for j in range(tab.shape[1]):
            axes[1,0].text(j, i, str(tab.iloc[i,j]), ha="center", va="center")
    axes[1,0].set_title("C. Expected label vs sequence call")
    # D support scatter
    axes[1,1].scatter(pd.to_numeric(df["species_IIIb_support"], errors="coerce"), pd.to_numeric(df["species_IIIc_support"], errors="coerce"), alpha=0.75)
    axes[1,1].plot([0,1],[0,1], linestyle="--", linewidth=1)
    axes[1,1].set_xlabel("IIIb support")
    axes[1,1].set_ylabel("IIIc support")
    axes[1,1].set_title("D. Species-calibrated sequence support")
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_v5_exon_motif_storyboard.png", dpi=300)
    plt.savefig(outdir / f"{prefix}_v5_exon_motif_storyboard.pdf")
    plt.close()


# ----------------------------- main -----------------------------




def recompute_species_call_for_threshold(row: pd.Series, threshold: float) -> str:
    if str(row.get("pair_audit_status", "")) != "III_region_sequence_distinct":
        return "uninformative"
    try:
        b = float(row.get("species_IIIb_support"))
        c = float(row.get("species_IIIc_support"))
    except Exception:
        return "uninformative"
    if not (math.isfinite(b) and math.isfinite(c)):
        return "uninformative"
    delta = b - c
    if delta >= threshold:
        return "IIIb"
    if delta <= -threshold:
        return "IIIc"
    return "ambiguous"



def main():
    ap = argparse.ArgumentParser(description="Exon/CDS-aware and sequence-aware FGFR2 IIIb/IIIc-region anchoring v5.8 structure-guided.")
    ap.add_argument("--query_fasta", required=True, type=Path)
    ap.add_argument("--human_iiib_segment_fasta", required=True, type=Path)
    ap.add_argument("--human_iiic_segment_fasta", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--prefix", default="fgfr2")
    ap.add_argument("--exons_tsv", type=Path, default=None, help="Optional exons.tsv from collect_fgfr2_models_dual_source_v3.py")
    ap.add_argument("--isoform_evidence_tsv", type=Path, default=None, help="Optional fgfr2_isoform_evidence.tsv from classify_fgfr2_IIIb_IIIc_by_exon_structure.py")
    ap.add_argument("--region_start", type=int, default=250, help="1-based fallback protein window start; used only when dynamic alignment anchoring fails")
    ap.add_argument("--region_end", type=int, default=430, help="1-based fallback protein window end; used only when dynamic alignment anchoring fails")
    ap.add_argument("--pair_identity_threshold", type=float, default=0.97)
    ap.add_argument("--pair_diff_threshold", type=int, default=5)
    ap.add_argument("--species_delta_threshold", type=float, default=0.12)
    ap.add_argument("--plausible_anchor_min", type=int, default=320)
    ap.add_argument("--plausible_anchor_max", type=int, default=430)
    ap.add_argument("--max_exon_plot_species", type=int, default=30)
    ap.add_argument("--representative_species", default="homo_sapiens,mus_musculus,canis_lupus_familiaris,gallus_gallus,anolis_carolinensis,danio_rerio,xenopus_tropicalis,takifugu_rubripes", help="Comma-separated species for readable representative exon-motif tracks")
    ap.add_argument("--sensitivity_thresholds", default="0.02,0.05,0.10,0.15,0.20,0.30", help="Comma-separated species-motif delta thresholds for robustness/sensitivity summary")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    records = parse_fasta_records(args.query_fasta)
    seg_b = read_fasta(args.human_iiib_segment_fasta)[0][1]
    seg_c = read_fasta(args.human_iiic_segment_fasta)[0][1]
    w0 = max(0, args.region_start - 1)
    w1 = max(w0, args.region_end)

    iso_ev = read_tsv_optional(args.isoform_evidence_tsv)
    exons = read_tsv_optional(args.exons_tsv)
    ev_lookup = make_evidence_lookup(iso_ev)

    audit, pair_lookup, window_fasta, diff_rows = build_pair_audit(records, seg_b, seg_c, w0, w1, args.pair_identity_threshold, args.pair_diff_threshold, ev_lookup)
    diff_df = pd.DataFrame(diff_rows)

    rows = []
    for rec in records:
        expected = norm_iso(rec.isoform)
        pair = pair_lookup.get(rec.species, {})
        row = {
            "query_id": rec.header,
            "species": rec.species,
            "role": rec.role,
            "source": rec.source,
            "transcript": rec.transcript,
            "protein": rec.protein,
            "expected_isoform": expected,
            "protein_length": len(rec.seq),
            "pair_audit_status": pair.get("pair_audit_status", "missing_pair_member"),
            "human_control_status": pair.get("human_control_status", "not_applicable"),
        }
        row.update(anchor_human_region(rec.seq, seg_b, seg_c, args.plausible_anchor_min, args.plausible_anchor_max))
        row.update(species_window_score(rec, pair, seg_b, seg_c, w0, w1, args.species_delta_threshold))
        row.update(join_exon_evidence(rec, ev_lookup))
        row["exon_sequence_agreement"] = "not_available"
        if row.get("exon_isoform_call") in {"IIIb", "IIIc"} and row.get("species_motif_call") in {"IIIb", "IIIc"}:
            row["exon_sequence_agreement"] = "agree" if row["exon_isoform_call"] == row["species_motif_call"] else "discordant"
        elif row.get("exon_isoform_call") in {"IIIb", "IIIc"} and row.get("species_motif_call") == "ambiguous":
            row["exon_sequence_agreement"] = "exon_supported_sequence_ambiguous"
        elif row.get("species_motif_call") == "uninformative":
            row["exon_sequence_agreement"] = "sequence_uninformative"
        row["final_anchor_status"] = final_status(row)
        row["final_anchor_status_simplified"] = simplify_final_status(row["final_anchor_status"])
        # Human-anchor weakness is tracked separately, so status labels can remain biologically interpretable.
        row["human_anchor_warning"] = "anchor_weak" if row.get("full_region_status") != "full_region_detected" else "none"
        rows.append(row)

    df = pd.DataFrame(rows)
    species_summary = df.groupby("species").agg(
        n_proteins=("query_id", "count"),
        final_statuses=("final_anchor_status", lambda x: ";".join(sorted(set(map(str, x))))),
        n_sequence_supported=("final_anchor_status", lambda x: sum(str(v).startswith("exon_sequence_supported") or str(v) == "region_detected_species_pair_supported" for v in x)),
        n_ambiguous=("final_anchor_status", lambda x: sum("ambiguous" in str(v) for v in x)),
        n_uninformative=("final_anchor_status", lambda x: sum("uninformative" in str(v) for v in x)),
        pair_audit_status=("pair_audit_status", lambda x: ";".join(sorted(set(map(str, x))))),
    ).reset_index()

    # Output tables
    df.to_csv(args.outdir / f"{args.prefix}_III_region_anchor_map.tsv", sep="\t", index=False)
    species_summary.to_csv(args.outdir / f"{args.prefix}_III_region_species_summary.tsv", sep="\t", index=False)
    audit.to_csv(args.outdir / f"{args.prefix}_III_pair_audit.tsv", sep="\t", index=False)
    diff_df.to_csv(args.outdir / f"{args.prefix}_III_pair_difference_positions.tsv", sep="\t", index=False)
    (args.outdir / f"{args.prefix}_III_species_pair_windows.fasta").write_text("\n".join(window_fasta) + "\n", encoding="utf-8")
    review = df[df["final_anchor_status"].astype(str).str.contains("ambiguous|uninformative|outlier|not_reliably|discordance", regex=True, na=False) | df["pair_audit_status"].astype(str).str.contains("human_positive_control_failed", regex=False, na=False)]
    review.to_csv(args.outdir / f"{args.prefix}_III_region_review_cases.tsv", sep="\t", index=False)
    agreement_cols = [c for c in ["species", "role", "source", "transcript", "protein", "expected_isoform", "exon_isoform_call", "exon_evidence_confidence", "exon_evidence_match_level", "exon_evidence_missing_reason", "exon_evidence_transcript_ids", "exon_evidence_source", "species_motif_call", "exon_sequence_agreement", "final_anchor_status"] if c in df.columns]
    df[agreement_cols].to_csv(args.outdir / f"{args.prefix}_III_exon_sequence_agreement.tsv", sep="\t", index=False)

    # Plots
    plot_counts(df["final_anchor_status"], args.outdir, f"{args.prefix}_v5_final_status_counts", "Detailed final III-region evidence status", "Number of candidate proteins")
    plot_simplified_status_counts(df, args.outdir, args.prefix)
    plot_counts(audit["pair_audit_status"], args.outdir, f"{args.prefix}_v5_pair_audit_status_counts", "IIIb/IIIc candidate pair audit", "Number of species")
    plot_pair_identity(audit, args.outdir, args.prefix, args.pair_identity_threshold)
    plot_matrix(df, "expected_isoform", "species_motif_call", args.outdir, f"{args.prefix}_v5_expected_vs_sequence_call_matrix", "Expected isoform vs species-calibrated sequence call")
    plot_matrix(df, "exon_isoform_call", "species_motif_call", args.outdir, f"{args.prefix}_v5_exon_vs_sequence_call_matrix", "Exon-structure call vs species-calibrated sequence call")
    if "exon_evidence_match_level" in df.columns:
        plot_counts(df["exon_evidence_match_level"], args.outdir, f"{args.prefix}_v5_exon_evidence_match_level_counts", "Exon evidence match level", "Number of candidate proteins")
    if "exon_evidence_missing_reason" in df.columns:
        missing_reasons = df.loc[df["exon_evidence_match_level"] == "not_available", "exon_evidence_missing_reason"]
        if len(missing_reasons):
            plot_counts(missing_reasons, args.outdir, f"{args.prefix}_v5_exon_missing_reason_counts", "Why exon evidence could not be joined", "Number of candidate proteins")
    plot_difference_map(diff_df, audit, args.outdir, args.prefix, w0, w1)
    plot_support_scatter(df, args.outdir, args.prefix)
    plot_species_stacked(df, args.outdir, args.prefix)
    plot_exon_architecture(records, exons, ev_lookup, args.outdir, args.prefix, args.max_exon_plot_species)
    plot_representative_exon_motif_tracks(records, exons, ev_lookup, df, args.outdir, args.prefix, parse_species_list(args.representative_species))
    plot_storyboard(df, audit, diff_df, args.outdir, args.prefix)
    plot_evidence_agreement_storyboard(df, audit, args.outdir, args.prefix)

    # Report
    with (args.outdir / f"{args.prefix}_III_region_anchor_report.md").open("w", encoding="utf-8") as f:
        f.write("# FGFR2 IIIb/IIIc-region anchoring v5.7 report\n\n")
        f.write("## Summary\n\n")
        f.write(f"- Candidate proteins analysed: {len(df)}\n")
        f.write(f"- Species represented: {df['species'].nunique()}\n")
        f.write(f"- Pair-audit status counts: {audit['pair_audit_status'].value_counts().to_dict()}\n")
        if 'human_control_status' in audit.columns:
            f.write(f"- Human positive-control status: {audit[audit['species'].astype(str).str.lower().isin(['homo_sapiens','human'])]['human_control_status'].tolist()}\n")
        f.write(f"- Species-calibrated sequence calls: {df['species_motif_call'].value_counts().to_dict()}\n")
        f.write(f"- Final integrated status counts: {df['final_anchor_status'].value_counts().to_dict()}\n")
        if not iso_ev.empty:
            f.write(f"- Exon evidence rows loaded: {len(iso_ev)}\n")
        else:
            f.write("- Exon evidence was not supplied; exon-aware agreement plots are limited.\n")
        if not exons.empty:
            f.write(f"- Exon rows loaded: {len(exons)}\n")
        else:
            f.write("- Exon table was not supplied; exon architecture plot was skipped.\n")
        f.write("\n## Interpretation\n\n")
        f.write("The v5.7 workflow treats the human-derived full-region anchor as a conservative localization QC rather than as a hard isoform classifier. The main IIIb/IIIc evidence comes from exon-structure calls, species-local IIIb/IIIc candidate-pair distinguishability and species-calibrated sequence support. Cases with weak human anchors but concordant exon/sequence evidence are retained as supported with an anchor-warning flag. Threshold-sensitivity tables and plots document the robustness of sequence calls to motif-delta cutoffs.\n")
    print(f"[OK] wrote outputs to {args.outdir}")
    print(df["final_anchor_status"].value_counts())
    print("\nSimplified final evidence classes:")
    print(df["final_anchor_status_simplified"].value_counts())


if __name__ == "__main__":
    main()
