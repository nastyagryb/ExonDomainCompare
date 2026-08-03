#!/usr/bin/env python3
"""
resolve_fgfr2_IIIb_IIIc_exons_exact_v1.py

Multi-evidence resolver for current-stage FGFR2 IIIb/IIIc exon-to-protein mapping.
No InterPro/domain mapping is performed.

Goal:
  Resolve the biologically most consistent IIIb/IIIc CDS/exon interval per species
  by combining:
    1) human-calibrated IIIb/IIIc protein-region anchors from pair audit,
    2) selected transcript/protein identity,
    3) CDS features exported by collect_fgfr2_models_dual_source_v3.py,
    4) optional mutually-exclusive exon structure metadata,
    5) pair-level distinctness of IIIb and IIIc genomic/CDS intervals.

This resolver intentionally DOES NOT use the structure exon as the only source of
truth, because annotation tables and protein-anchor tables can refer to different
transcript/protein versions. Instead, structure metadata contributes evidence;
the final candidate must still be transcript/protein compatible and overlap the
expected dynamic IIIb/IIIc protein window.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd

VERSION = "v2.22_refined_evidence_native_normalized_resolver"



# Refined resolver evidence labels.
# The single broad legacy label `gold_exact_multi_evidence_CDS_pair` hid whether
# the CDS feature was matched using transcript+protein evidence or only one of
# them. These helpers split that into explicit evidence levels / match types /
# confidence / refined status, while keeping the legacy label in a column.

def classify_resolver_evidence(ov: int, txm: bool, prm: bool, same_pair_key: bool,
                               identity_incomplete: bool, has_candidate: bool) -> Dict[str, str]:
    """Return refined resolver evidence labels for a single IIIb/IIIc CDS choice.

    Rules:
      * transcript AND protein match -> gold_exact_transcript_and_protein_CDS_pair
      * transcript match only         -> gold_exact_transcript_specific_CDS_pair
      * protein match only            -> silver_exact_protein_specific_CDS_pair
      * candidate CDS in window only  -> bronze_candidate_CDS_in_window
      * identity incomplete           -> review_identity_incomplete
      * otherwise                     -> unresolved_no_CDS_match
    """
    if not has_candidate:
        return {
            "resolver_evidence_level": "unresolved",
            "resolver_match_type": "unresolved",
            "resolver_confidence": "unresolved",
            "resolver_status_refined": "unresolved_no_CDS_match",
        }
    if same_pair_key:
        return {
            "resolver_evidence_level": "incomplete_identity",
            "resolver_match_type": "species_level_candidate_CDS",
            "resolver_confidence": "low",
            "resolver_status_refined": "review_identity_incomplete",
        }
    if ov > 0 and txm and prm:
        return {
            "resolver_evidence_level": "transcript_and_protein",
            "resolver_match_type": "exact_transcript_exact_protein_CDS",
            "resolver_confidence": "high",
            "resolver_status_refined": "gold_exact_transcript_and_protein_CDS_pair",
        }
    if ov > 0 and txm:
        return {
            "resolver_evidence_level": "transcript_only",
            "resolver_match_type": "exact_transcript_CDS",
            "resolver_confidence": "moderate",
            "resolver_status_refined": "gold_exact_transcript_specific_CDS_pair",
        }
    if ov > 0 and prm:
        return {
            "resolver_evidence_level": "protein_only",
            "resolver_match_type": "exact_protein_CDS",
            "resolver_confidence": "moderate",
            "resolver_status_refined": "silver_exact_protein_specific_CDS_pair",
        }
    if ov > 0 and not identity_incomplete:
        return {
            "resolver_evidence_level": "species_level_fallback",
            "resolver_match_type": "candidate_CDS_in_expected_window",
            "resolver_confidence": "low",
            "resolver_status_refined": "bronze_candidate_CDS_in_window",
        }
    return {
        "resolver_evidence_level": "incomplete_identity",
        "resolver_match_type": "candidate_CDS_in_expected_window" if ov > 0 else "unresolved",
        "resolver_confidence": "low",
        "resolver_status_refined": "review_identity_incomplete" if ov > 0 else "unresolved_no_CDS_match",
    }


def codon_boundary_precision(phase: str, coordinate_source: str, warning: str) -> Tuple[str, str]:
    """Return explicit CDS codon-phase and boundary precision.

    NCBI/Ensembl CDS phase indicates how many bases of the first codon lie in the
    previous exon. Phase 0 means the exon boundary coincides with a codon boundary;
    phase 1/2 means the IIIb/IIIc protein-axis start falls inside a split codon, so
    the amino-acid boundary is approximate.
    """
    p = str(phase or "").strip()
    warn = str(warning or "").strip()
    if p in ("", "nan", "None", "-1", "."):
        precision = "unknown_codon_phase"
    elif p in ("0", "0.0"):
        precision = "codon_boundary_exact"
    elif p in ("1", "2", "1.0", "2.0"):
        precision = "codon_split_phase_offset"
    else:
        precision = "unknown_codon_phase"
    note_bits = []
    if precision == "codon_split_phase_offset":
        note_bits.append(f"cds_phase={p}_amino_acid_boundary_approximate")
    if str(coordinate_source or "").strip():
        note_bits.append(f"coordinate_source={coordinate_source}")
    if warn:
        note_bits.append(f"cds_warning={warn}")
    return precision, ";".join(note_bits)


def _phase_int(phase: str) -> Optional[int]:
    p = str(phase or "").strip()
    if p in ("", "nan", "None", "-1", "."):
        return None
    try:
        v = int(float(p))
    except Exception:
        return None
    return v if v in (0, 1, 2) else None


def cds_boundary_precision_lr(phase: str, cds_len_bp: Optional[int], coordinate_source: str,
                              warning: str, strand: str) -> Dict[str, str]:
    """Return left/right CDS-boundary precision from GFF3 phase.

    GFF3 phase is defined in transcript/translation order (not raw genomic order),
    so it is correct on the negative strand without flipping. ``phase`` is the
    number of bases of the first codon that lie in the previous exon, i.e. the
    transcript-5' (left) boundary is codon-aligned iff phase == 0. The transcript-3'
    (right) boundary is codon-aligned iff (phase + coding_length) % 3 == 0. A single
    exon's genomic length equals its coding length, so both boundaries are derived
    from the start phase plus the CDS-exon length.
    """
    p = _phase_int(phase)
    out = {
        "cds_left_phase": "", "cds_right_phase": "",
        "cds_left_boundary_precision": "unknown_codon_phase",
        "cds_right_boundary_precision": "unknown_codon_phase",
        "cds_boundary_precision_refined": "unknown_codon_phase",
        "cds_boundary_confidence": "low",
        "cds_phase_source": "phase_missing_or_uninterpretable",
        "cds_phase_warning": "",
    }
    warn_bits: List[str] = []
    if str(warning or "").strip():
        warn_bits.append(f"cds_warning={str(warning).strip()}")
    if p is None or not isinstance(cds_len_bp, int) or cds_len_bp <= 0:
        if p is None:
            warn_bits.append("cds_phase_missing_or_uninterpretable")
        if not isinstance(cds_len_bp, int) or cds_len_bp <= 0:
            warn_bits.append("cds_length_unavailable")
        out["cds_phase_warning"] = ";".join(warn_bits)
        return out
    right_offset = (p + cds_len_bp) % 3  # 0 == codon-aligned 3' boundary
    left_exact = (p == 0)
    right_exact = (right_offset == 0)
    left_prec = "exact" if left_exact else "codon_split"
    right_prec = "exact" if right_exact else "codon_split"
    if left_exact and right_exact:
        refined, conf = "exact", "high"
    elif (not left_exact) and (not right_exact):
        refined, conf = "codon_split_both_sides", "low"
    else:
        refined, conf = "codon_split_one_side", "medium"
    if not str(strand or "").strip():
        warn_bits.append("strand_unknown_phase_assumed_transcript_order")
    out.update({
        "cds_left_phase": str(p),
        "cds_right_phase": str(right_offset),
        "cds_left_boundary_precision": left_prec,
        "cds_right_boundary_precision": right_prec,
        "cds_boundary_precision_refined": refined,
        "cds_boundary_confidence": conf,
        "cds_phase_source": "gff3_cds_phase_transcript_order",
        "cds_phase_warning": ";".join(warn_bits),
    })
    return out


def read_tsv(path: Path, required: bool = True) -> pd.DataFrame:
    if path is None or str(path) == "" or not Path(path).exists() or Path(path).stat().st_size == 0:
        if required:
            raise SystemExit(f"[ERROR] Missing or empty input: {path}")
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def norm_col(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d.columns = [str(c).strip() for c in d.columns]
    return d


def first(row: pd.Series, names: Iterable[str], default: str = "") -> str:
    for n in names:
        if n in row.index:
            v = row.get(n, "")
            if pd.notna(v) and str(v).strip() not in ("", "nan", "None"):
                return str(v).strip()
    return default


def col_first(df: pd.DataFrame, names: Iterable[str]) -> Optional[str]:
    for n in names:
        if n in df.columns:
            return n
    return None


def to_int(x: Any) -> Optional[int]:
    try:
        if x is None or str(x).strip() == "" or str(x).lower() == "nan":
            return None
        return int(float(str(x)))
    except Exception:
        return None


def clean_id(x: Any) -> str:
    s = str(x or "").strip()
    if not s or s.lower() == "nan":
        return ""
    # Remove common GFF3 prefixes but keep version suffix.
    s = re.sub(r"^(cds-|rna-|gene-|exon-|protein-)", "", s)
    return s


def strip_version(x: Any) -> str:
    s = clean_id(x)
    return re.sub(r"\.\d+$", "", s)


def ids_match(a: Any, b: Any) -> bool:
    aa, bb = clean_id(a), clean_id(b)
    if not aa or not bb:
        return False
    return aa == bb or strip_version(aa) == strip_version(bb) or aa in bb or bb in aa


def any_match(value: str, candidates: Iterable[str]) -> bool:
    return any(ids_match(value, c) for c in candidates if str(c).strip())


def overlap_len(a1: Optional[int], a2: Optional[int], b1: Optional[int], b2: Optional[int]) -> int:
    if None in (a1, a2, b1, b2):
        return 0
    lo = max(int(a1), int(b1))
    hi = min(int(a2), int(b2))
    return max(0, hi - lo + 1)


def interval_len(a1: Optional[int], a2: Optional[int]) -> int:
    if a1 is None or a2 is None:
        return 0
    return max(0, int(a2) - int(a1) + 1)


def center_distance(a1: Optional[int], a2: Optional[int], b1: Optional[int], b2: Optional[int]) -> Optional[float]:
    if None in (a1, a2, b1, b2):
        return None
    return abs(((int(a1) + int(a2)) / 2.0) - ((int(b1) + int(b2)) / 2.0))


def get_numeric(row: pd.Series, names: Iterable[str]) -> Optional[int]:
    return to_int(first(row, names))


def candidate_key(row: pd.Series) -> str:
    chrom = first(row, ["chrom", "seqid", "seq_id", "contig", "accession"])
    start = first(row, ["coding_exon_cds_start", "cds_start", "start", "genomic_start", "exon_start"])
    end = first(row, ["coding_exon_cds_end", "cds_end", "end", "genomic_end", "exon_end"])
    strand = first(row, ["strand"])
    p1 = first(row, ["protein_start_aa", "exon_protein_start_aa", "aa_start", "protein_start"])
    p2 = first(row, ["protein_end_aa", "exon_protein_end_aa", "aa_end", "protein_end"])
    raw = first(row, ["coding_exon_key", "exon_id_source", "exon_id", "raw_cds_ids", "cds_id", "ID"])
    return "|".join([clean_id(raw), chrom, start, end, strand, p1, p2])


def genomic_overlap(cds_row: pd.Series, exon_row: pd.Series) -> int:
    cchrom = first(cds_row, ["chrom", "seqid", "seq_id", "contig", "accession"])
    echrom = first(exon_row, ["chrom", "seqid", "seq_id", "contig", "accession"])
    if cchrom and echrom and cchrom != echrom:
        return 0
    c1 = get_numeric(cds_row, ["coding_exon_cds_start", "cds_start", "start", "genomic_start", "exon_start"])
    c2 = get_numeric(cds_row, ["coding_exon_cds_end", "cds_end", "end", "genomic_end", "exon_end"])
    e1 = get_numeric(exon_row, ["exon_start", "start", "genomic_start", "coding_exon_cds_start"])
    e2 = get_numeric(exon_row, ["exon_end", "end", "genomic_end", "coding_exon_cds_end"])
    return overlap_len(c1, c2, e1, e2)


def subset_species(df: pd.DataFrame, species: str, source: str = "") -> pd.DataFrame:
    if df.empty:
        return df
    d = df.copy()
    sp_col = col_first(d, ["species_canonical", "species", "species_name"])
    if sp_col:
        d = d[d[sp_col].astype(str) == species]
    src_col = col_first(d, ["source_db", "source"])
    if source and src_col:
        dd = d[d[src_col].astype(str) == source]
        if not dd.empty:
            d = dd
    return d



def normalize_pair_audit(pair: pd.DataFrame) -> pd.DataFrame:
    """Normalize either long or wide 6e pair-audit/anchor tables to one row per species+isoform.

    The original 6e fgfr2_III_pair_audit.tsv is WIDE: one row per species with
    IIIb_transcript, IIIc_transcript, IIIb_protein, IIIc_protein,
    IIIb_window_start_1based, IIIc_window_start_1based, etc.
    Earlier plotting code expected LONG rows. v2.20 returned an empty resolver
    table when it received the wide table. This function fixes that at the
    mapping layer, not only at plotting.
    """
    pair = norm_col(pair)
    if pair.empty:
        return pair
    iso_col = col_first(pair, ["inferred_isoform", "isoform", "expected_isoform"])
    # Already long: standardize a few column names and return.
    if iso_col is not None:
        d = pair.copy()
        if iso_col != "inferred_isoform":
            d["inferred_isoform"] = d[iso_col]
        sp_col = col_first(d, ["species_canonical", "species", "species_name"])
        if sp_col and sp_col != "species_canonical":
            d["species_canonical"] = d[sp_col]
        src_col = col_first(d, ["source_db", "source"])
        if src_col and src_col != "source_db":
            d["source_db"] = d[src_col]
        tx_col = col_first(d, ["transcript_id_source", "transcript", "transcript_id", "transcript_id_internal"])
        if tx_col and "transcript_id_source" not in d.columns:
            d["transcript_id_source"] = d[tx_col]
        pr_col = col_first(d, ["protein_id", "protein", "translation_id_source", "translation_id"])
        if pr_col and "protein_id" not in d.columns:
            d["protein_id"] = d[pr_col]
        pl_col = col_first(d, ["protein_length_aa", "protein_length", "length_aa", "IIIb_length", "IIIc_length"])
        if pl_col and "protein_length_aa" not in d.columns:
            d["protein_length_aa"] = d[pl_col]
        # Normalize anchor coordinates if possible.
        s_col = col_first(d, ["III_region_start_aa", "region_start_aa", "protein_region_start_aa", "start_aa", "window_start_1based"])
        e_col = col_first(d, ["III_region_end_aa", "region_end_aa", "protein_region_end_aa", "end_aa", "window_end_1based"])
        if s_col and "III_region_start_aa" not in d.columns:
            d["III_region_start_aa"] = d[s_col]
        if e_col and "III_region_end_aa" not in d.columns:
            d["III_region_end_aa"] = d[e_col]
        return d

    # Wide 6e pair audit: one species row -> two rows.
    has_wide = any(c.startswith("IIIb_") for c in pair.columns) and any(c.startswith("IIIc_") for c in pair.columns)
    if not has_wide:
        return pair
    out = []
    for _, r in pair.iterrows():
        species = first(r, ["species_canonical", "species", "species_name"])
        source = first(r, ["source_db", "source"])
        pair_status = first(r, ["pair_audit_status"])
        human_status = first(r, ["human_control_status"])
        # Carry refined similarity fields when present.
        similarity_class = first(r, ["iii_region_similarity_class"])
        similarity_warning = first(r, ["iii_region_similarity_warning"])
        legacy_pair_status = first(r, ["legacy_pair_audit_status"], pair_status)
        for iso in ["IIIb", "IIIc"]:
            rec = {
                "species_canonical": species,
                "source_db": source,
                "inferred_isoform": iso,
                "transcript_id_source": first(r, [f"{iso}_transcript", f"{iso}_transcript_id", f"{iso}_transcript_id_source"]),
                "protein_id": first(r, [f"{iso}_protein", f"{iso}_protein_id", f"{iso}_translation_id_source"]),
                "translation_id_source": first(r, [f"{iso}_protein", f"{iso}_protein_id", f"{iso}_translation_id_source"]),
                "protein_length_aa": first(r, [f"{iso}_length", f"{iso}_protein_length", f"{iso}_protein_length_aa"]),
                "III_region_start_aa": first(r, [f"{iso}_window_start_1based", f"{iso}_region_start_aa", f"{iso}_protein_region_start_aa"]),
                "III_region_end_aa": first(r, [f"{iso}_window_end_1based", f"{iso}_region_end_aa", f"{iso}_protein_region_end_aa"]),
                "III_region_source": first(r, [f"{iso}_window_mode", f"{iso}_anchor_status", "III_region_source"]),
                "pair_audit_status": pair_status,
                "legacy_pair_audit_status": legacy_pair_status,
                "iii_region_similarity_class": similarity_class,
                "iii_region_similarity_warning": similarity_warning,
                "human_control_status": human_status,
                "pair_table_format": "wide_6e_pair_audit_normalized_to_long",
            }
            out.append(rec)
    return pd.DataFrame(out)


def selected_for_pair(selected: pd.DataFrame, species: str, source: str, tx: str, protein: str) -> pd.DataFrame:
    d = subset_species(selected, species, source)
    if d.empty:
        return d
    tx_cols = [c for c in ["transcript_id_source", "transcript_id", "transcript", "transcript_id_internal"] if c in d.columns]
    pr_cols = [c for c in ["protein_id", "translation_id_source", "translation_id"] if c in d.columns]
    masks = []
    if tx and tx_cols:
        masks.append(pd.Series(False, index=d.index))
        for c in tx_cols:
            masks[-1] = masks[-1] | d[c].map(lambda v: ids_match(v, tx))
    if protein and pr_cols:
        masks.append(pd.Series(False, index=d.index))
        for c in pr_cols:
            masks[-1] = masks[-1] | d[c].map(lambda v: ids_match(v, protein))
    if masks:
        m = masks[0]
        for mm in masks[1:]:
            m = m | mm
        dd = d[m]
        if not dd.empty:
            return dd
    return d.head(1)


def alt_candidates(alt: pd.DataFrame, species: str, source: str, iso: str) -> pd.DataFrame:
    if alt.empty:
        return alt
    d = subset_species(alt, species, source)
    if d.empty:
        return d
    # Keep rows that explicitly refer to the isoform if possible.
    iso_cols = [c for c in ["inferred_isoform", "isoform", "iii_isoform", "slot", "slot_label", "fgfr2_iii_slot", "candidate_isoform"] if c in d.columns]
    if iso_cols:
        mask = pd.Series(False, index=d.index)
        for c in iso_cols:
            mask = mask | d[c].astype(str).str.contains(iso, case=False, na=False)
        dd = d[mask]
        if not dd.empty:
            d = dd
    if "selected_as_fgfr2_iii_slot" in d.columns:
        dd = d[d["selected_as_fgfr2_iii_slot"].astype(str).str.lower().isin(["1", "true", "yes", "y"])]
        if not dd.empty:
            d = dd
    return d


def cds_subset_by_identity(cds: pd.DataFrame, species: str, source: str, tx: str, protein: str, selected_hits: pd.DataFrame) -> Tuple[pd.DataFrame, str]:
    d = subset_species(cds, species, source)
    if d.empty:
        return d, "no_species_cds"
    tx_candidates = {tx}
    pr_candidates = {protein}
    if selected_hits is not None and not selected_hits.empty:
        for _, r in selected_hits.iterrows():
            for c in ["transcript_id_source", "transcript_id", "transcript", "transcript_id_internal"]:
                if c in r.index:
                    tx_candidates.add(str(r.get(c, "")))
            for c in ["protein_id", "translation_id_source", "translation_id"]:
                if c in r.index:
                    pr_candidates.add(str(r.get(c, "")))
    tx_cols = [c for c in ["transcript_id_source", "transcript_id", "transcript", "transcript_id_internal", "parent_mrna", "Parent"] if c in d.columns]
    pr_cols = [c for c in ["protein_id", "translation_id_source", "translation_id", "product", "protein_product", "raw_cds_ids", "cds_id", "ID"] if c in d.columns]
    def tx_match_row(r):
        return any(any_match(str(r.get(c,"")), tx_candidates) for c in tx_cols)
    def pr_match_row(r):
        return any(any_match(str(r.get(c,"")), pr_candidates) for c in pr_cols)
    if tx_cols and pr_cols:
        m = d.apply(lambda r: tx_match_row(r) and pr_match_row(r), axis=1)
        dd = d[m]
        if not dd.empty:
            return dd, "transcript_and_protein_specific"
    if tx_cols:
        m = d.apply(tx_match_row, axis=1)
        dd = d[m]
        if not dd.empty:
            return dd, "transcript_specific"
    if pr_cols:
        m = d.apply(pr_match_row, axis=1)
        dd = d[m]
        if not dd.empty:
            return dd, "protein_specific"
    return d, "species_only_fallback"


def score_candidates(cds_df: pd.DataFrame, pr: pd.Series, selected_hits: pd.DataFrame, alt_rows: pd.DataFrame) -> pd.DataFrame:
    if cds_df.empty:
        return pd.DataFrame()
    species = first(pr, ["species_canonical", "species"])
    iso = first(pr, ["inferred_isoform", "isoform"])
    tx = first(pr, ["transcript_id_source", "transcript_id", "transcript_id_internal"])
    protein = first(pr, ["protein_id", "translation_id_source", "translation_id"])
    win_s = get_numeric(pr, ["III_region_start_aa", "region_start_aa", "protein_region_start_aa", "start_aa"])
    win_e = get_numeric(pr, ["III_region_end_aa", "region_end_aa", "protein_region_end_aa", "end_aa"])

    selected_tx = set([tx])
    selected_pr = set([protein])
    if selected_hits is not None and not selected_hits.empty:
        for _, r in selected_hits.iterrows():
            for c in ["transcript_id_source", "transcript_id", "transcript", "transcript_id_internal"]:
                if c in r.index:
                    selected_tx.add(str(r.get(c, "")))
            for c in ["protein_id", "translation_id_source", "translation_id"]:
                if c in r.index:
                    selected_pr.add(str(r.get(c, "")))

    rows = []
    for _, c in cds_df.iterrows():
        p1 = get_numeric(c, ["protein_start_aa", "exon_protein_start_aa", "aa_start", "protein_start"])
        p2 = get_numeric(c, ["protein_end_aa", "exon_protein_end_aa", "aa_end", "protein_end"])
        if p1 is None or p2 is None:
            continue
        ov = overlap_len(p1, p2, win_s, win_e)
        clen = interval_len(p1, p2)
        wlen = interval_len(win_s, win_e)
        ov_frac_cds = ov / clen if clen else 0
        ov_frac_win = ov / wlen if wlen else 0
        dist = center_distance(p1, p2, win_s, win_e)
        # transcript/protein match details
        tx_cols = [cc for cc in ["transcript_id_source", "transcript_id", "transcript", "transcript_id_internal", "parent_mrna", "Parent"] if cc in c.index]
        pr_cols = [cc for cc in ["protein_id", "translation_id_source", "translation_id", "product", "protein_product", "raw_cds_ids", "cds_id", "ID"] if cc in c.index]
        tx_match = any(any_match(str(c.get(cc,"")), selected_tx) for cc in tx_cols)
        pr_match = any(any_match(str(c.get(cc,"")), selected_pr) for cc in pr_cols)
        alt_ov = 0
        if alt_rows is not None and not alt_rows.empty:
            for _, ar in alt_rows.iterrows():
                alt_ov = max(alt_ov, genomic_overlap(c, ar))
        # Scoring: protein window and transcript/protein identity are primary;
        # alternative exon structure is a bonus, not a hard constraint.
        score = 0.0
        if tx_match:
            score += 30
        if pr_match:
            score += 35
        if ov > 0:
            score += 80 + 20 * ov_frac_cds + 20 * ov_frac_win
        else:
            # Allow candidate but penalize strongly; useful for diagnostics.
            if dist is not None:
                score -= min(80, dist / 2)
            else:
                score -= 50
        if alt_ov > 0:
            score += 25
        # Prefer a CDS interval length consistent with a mutually exclusive III exon/core (roughly 35-80 aa).
        if 30 <= clen <= 90:
            score += 10
        elif clen > 120:
            score -= 20
        strict_key = candidate_key(c)
        status = "candidate"
        if ov <= 0:
            status = "outside_expected_III_window"
        elif not (tx_match or pr_match):
            status = "window_match_but_not_tx_or_protein_specific"
        elif not tx_match:
            status = "window_match_protein_only"
        elif not pr_match:
            status = "window_match_transcript_only"
        elif alt_rows is not None and not alt_rows.empty and alt_ov <= 0:
            status = "window_tx_protein_match_no_structure_overlap"
        else:
            status = "window_tx_protein_structure_supported" if alt_ov > 0 else "window_tx_protein_match"
        rec = {k: first(c, [k]) for k in c.index}
        rec.update({
            "species_canonical": species,
            "inferred_isoform": iso,
            "candidate_score": round(score, 3),
            "candidate_status": status,
            "candidate_window_overlap_aa": ov,
            "candidate_window_overlap_fraction_cds": round(ov_frac_cds, 3),
            "candidate_window_overlap_fraction_window": round(ov_frac_win, 3),
            "candidate_center_distance_to_window": "" if dist is None else round(dist, 3),
            "candidate_transcript_match": int(bool(tx_match)),
            "candidate_protein_match": int(bool(pr_match)),
            "candidate_structure_exon_overlap_bp": alt_ov,
            "strict_cds_exon_key": strict_key,
            "candidate_protein_start_aa": p1,
            "candidate_protein_end_aa": p2,
        })
        rows.append(rec)
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows).sort_values(["candidate_score", "candidate_window_overlap_aa"], ascending=[False, False])
    return out


def resolve_pair_for_species(pair_rows: pd.DataFrame, candidates_by_iso: Dict[str, pd.DataFrame]) -> Dict[str, Optional[pd.Series]]:
    # Choose the best distinct pair if both isoforms exist.
    choices: Dict[str, Optional[pd.Series]] = {"IIIb": None, "IIIc": None}
    b = candidates_by_iso.get("IIIb", pd.DataFrame())
    c = candidates_by_iso.get("IIIc", pd.DataFrame())
    if not b.empty and not c.empty:
        best = None
        for _, rb in b.head(20).iterrows():
            for _, rc in c.head(20).iterrows():
                kb = first(rb, ["strict_cds_exon_key"])
                kc = first(rc, ["strict_cds_exon_key"])
                # Distinct genomic/CDS interval required for gold. If same key, penalize heavily but keep as fallback.
                penalty = -10000 if kb and kc and kb == kc else 0
                # Both must overlap expected windows for exact.
                exact_bonus = 100 if to_int(rb.get("candidate_window_overlap_aa")) and to_int(rc.get("candidate_window_overlap_aa")) else 0
                score = float(rb.get("candidate_score", 0)) + float(rc.get("candidate_score", 0)) + penalty + exact_bonus
                if best is None or score > best[0]:
                    best = (score, rb, rc, penalty)
        if best is not None:
            choices["IIIb"] = best[1]
            choices["IIIc"] = best[2]
            return choices
    for iso in ["IIIb", "IIIc"]:
        d = candidates_by_iso.get(iso, pd.DataFrame())
        if not d.empty:
            choices[iso] = d.iloc[0]
    return choices


def build_resolved(pair: pd.DataFrame, cds: pd.DataFrame, selected: pd.DataFrame, alt: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    pair = normalize_pair_audit(pair)
    cds = norm_col(cds)
    selected = norm_col(selected)
    alt = norm_col(alt) if not alt.empty else pd.DataFrame()

    all_candidates = []
    unresolved = []
    resolved_rows = []

    # prepare group by species source
    for species, grp in pair.groupby(pair[col_first(pair, ["species_canonical", "species"])]):
        # Group by source if available, otherwise within species.
        src_col = col_first(grp, ["source_db", "source"])
        subgroups = grp.groupby(grp[src_col]) if src_col else [("", grp)]
        for source, sg in subgroups:
            candidates_by_iso: Dict[str, pd.DataFrame] = {}
            pr_by_iso: Dict[str, pd.Series] = {}
            for _, pr in sg.iterrows():
                iso = first(pr, ["inferred_isoform", "isoform"])
                if iso not in ("IIIb", "IIIc"):
                    continue
                pr_by_iso[iso] = pr
                tx = first(pr, ["transcript_id_source", "transcript_id", "transcript_id_internal"])
                protein = first(pr, ["protein_id", "translation_id_source", "translation_id"])
                sh = selected_for_pair(selected, str(species), str(source), tx, protein)
                cds_sub, filter_status = cds_subset_by_identity(cds, str(species), str(source), tx, protein, sh)
                ar = alt_candidates(alt, str(species), str(source), iso)
                cand = score_candidates(cds_sub, pr, sh, ar)
                if not cand.empty:
                    cand["identity_filter_status"] = filter_status
                    candidates_by_iso[iso] = cand
                    all_candidates.append(cand)
                else:
                    candidates_by_iso[iso] = pd.DataFrame()
                    unresolved.append({
                        "species_canonical": species, "source_db": source, "inferred_isoform": iso,
                        "transcript_id_source": tx, "protein_id": protein,
                        "reason": f"no_CDS_candidate_after_{filter_status}",
                    })
            choices = resolve_pair_for_species(sg, candidates_by_iso)
            keys = {iso: first(ch, ["strict_cds_exon_key"]) if ch is not None else "" for iso, ch in choices.items()}
            same_pair_key = bool(keys.get("IIIb") and keys.get("IIIc") and keys["IIIb"] == keys["IIIc"])
            pair_has_both = choices.get("IIIb") is not None and choices.get("IIIc") is not None
            for iso in ["IIIb", "IIIc"]:
                pr = pr_by_iso.get(iso)
                ch = choices.get(iso)
                if pr is None:
                    continue
                base = {c: first(pr, [c]) for c in pr.index}
                if ch is None:
                    base.update({
                        "resolved_status": "anchor_only_no_CDS_candidate",
                        "exon_coordinate_status": "anchor_only_no_CDS_candidate",
                        "strict_pair_class": "not_resolved_pair",
                        "main_figure_exact_bool": False,
                        "legacy_resolver_status": "anchor_only_no_CDS_candidate",
                        "resolver_evidence_level": "unresolved",
                        "resolver_match_type": "unresolved",
                        "resolver_confidence": "unresolved",
                        "resolver_status_refined": "unresolved_no_CDS_match",
                        "native_protein_start_aa": "",
                        "native_protein_end_aa": "",
                        "native_protein_center_aa": "",
                        "native_protein_length_aa": "",
                        "iii_slot_start_aa": "",
                        "iii_slot_end_aa": "",
                        "iii_slot_center_aa": "",
                        "iii_slot_length_aa": "",
                        "iii_slot_coordinate_method": "not_normalized_insufficient_evidence",
                        "iii_slot_coordinate_confidence": "unresolved",
                        "iii_slot_coordinate_note": "insufficient_evidence_for_normalization",
                        "cds_phase_value": "",
                        "cds_boundary_precision": "unknown_codon_phase",
                        "legacy_cds_boundary_precision": "unknown_codon_phase",
                        "cds_boundary_precision_note": "",
                        "cds_left_phase": "",
                        "cds_right_phase": "",
                        "cds_left_boundary_precision": "unknown_codon_phase",
                        "cds_right_boundary_precision": "unknown_codon_phase",
                        "cds_boundary_precision_refined": "unknown_codon_phase",
                        "cds_boundary_confidence": "low",
                        "cds_phase_source": "unresolved_no_CDS_match",
                        "cds_phase_warning": "no_resolved_CDS_exon",
                    })
                    resolved_rows.append(base)
                    continue
                # Determine status.
                ov = to_int(ch.get("candidate_window_overlap_aa")) or 0
                txm = str(ch.get("candidate_transcript_match", "0")) in ("1", "True", "true")
                prm = str(ch.get("candidate_protein_match", "0")) in ("1", "True", "true")
                structure_ov = to_int(ch.get("candidate_structure_exon_overlap_bp")) or 0
                identity_incomplete = False
                if same_pair_key:
                    status = "ambiguous_same_CDS_feature"
                    pair_class = "ambiguous_same_CDS_feature"
                    exact = False
                    identity_incomplete = True
                elif ov > 0 and (txm or prm) and pair_has_both:
                    # Legacy gold label retained verbatim; refined labels below.
                    status = "gold_exact_multi_evidence_CDS_pair"
                    pair_class = "gold_exact_distinct_IIIb_IIIc_CDS_pair"
                    exact = True
                elif ov > 0:
                    status = "candidate_CDS_in_window_but_identity_incomplete"
                    pair_class = "not_gold_identity_incomplete"
                    exact = False
                    identity_incomplete = True
                else:
                    status = "candidate_CDS_outside_expected_III_window"
                    pair_class = "not_gold_outside_expected_window"
                    exact = False

                # Keep the legacy resolver status beside refined evidence labels.
                refined = classify_resolver_evidence(
                    ov=ov, txm=txm, prm=prm, same_pair_key=same_pair_key,
                    identity_incomplete=identity_incomplete, has_candidate=ch is not None,
                )

                # Native protein coordinates use the transcript-specific axis.
                p1 = to_int(ch.get("candidate_protein_start_aa"))
                p2 = to_int(ch.get("candidate_protein_end_aa"))
                if p1 is not None and p2 is not None:
                    native_center = round((p1 + p2) / 2.0, 1)
                    native_length = p2 - p1 + 1
                else:
                    native_center = ""
                    native_length = ""

                # Normalized III-slot coordinates use an exon-internal axis.
                # IMPORTANT: this is NOT an anchor/human-alignment biological
                # normalization. Each mutually-exclusive cassette's own start is set
                # to origin 1, so the axis normalizes cassette SHAPE/LENGTH, not the
                # absolute biological position of the slot in the protein. Confidence
                # is therefore capped below "high" until a true anchor-based or
                # human-alignment-based normalization is implemented.
                if isinstance(native_length, int) and native_length > 0:
                    iii_slot_start = 1
                    iii_slot_end = native_length
                    iii_slot_center = round((1 + native_length) / 2.0, 1)
                    iii_slot_len = native_length
                    iii_slot_method = "exon_internal_relative"
                    # Cap confidence: best evidence -> "moderate", weaker -> "low".
                    conf_map = {"high": "moderate", "moderate": "low", "low": "low", "unresolved": "unresolved"}
                    iii_slot_conf = conf_map.get(refined["resolver_confidence"], "low")
                    iii_slot_note = "normalizes_cassette_shape_length_not_absolute_biological_position"
                else:
                    iii_slot_start = iii_slot_end = iii_slot_center = iii_slot_len = ""
                    iii_slot_method = "not_normalized_insufficient_evidence"
                    iii_slot_conf = "unresolved"
                    iii_slot_note = "insufficient_evidence_for_normalization"

                # Retain the legacy single-value CDS boundary precision.
                cds_precision, cds_precision_note = codon_boundary_precision(
                    first(ch, ["phase", "phase_values"]),
                    first(ch, ["coordinate_source"]),
                    first(ch, ["warning"]),
                )
                # Store explicit left and right CDS-boundary precision.
                _cds_start = to_int(first(ch, ["coding_exon_cds_start", "cds_start", "start", "genomic_start", "exon_start"]))
                _cds_end = to_int(first(ch, ["coding_exon_cds_end", "cds_end", "end", "genomic_end", "exon_end"]))
                _cds_len_bp = to_int(first(ch, ["cds_length_bp", "coding_exon_cds_length_bp"]))
                if _cds_len_bp is None and _cds_start is not None and _cds_end is not None:
                    _cds_len_bp = abs(_cds_end - _cds_start) + 1
                cds_lr = cds_boundary_precision_lr(
                    first(ch, ["phase", "phase_values"]), _cds_len_bp,
                    first(ch, ["coordinate_source"]), first(ch, ["warning"]),
                    first(ch, ["strand"]),
                )

                # Flatten chosen candidate columns with resolver_ prefix where needed.
                chosen = {f"resolver_{k}": ch.get(k, "") for k in ch.index}
                # Standardized coordinate columns for plotting/downstream.
                base.update(chosen)
                base.update({
                    "legacy_resolver_status": status,
                    "resolver_evidence_level": refined["resolver_evidence_level"],
                    "resolver_match_type": refined["resolver_match_type"],
                    "resolver_confidence": refined["resolver_confidence"],
                    "resolver_status_refined": refined["resolver_status_refined"],
                    "native_protein_start_aa": p1 if p1 is not None else "",
                    "native_protein_end_aa": p2 if p2 is not None else "",
                    "native_protein_center_aa": native_center,
                    "native_protein_length_aa": native_length,
                    "iii_slot_start_aa": iii_slot_start,
                    "iii_slot_end_aa": iii_slot_end,
                    "iii_slot_center_aa": iii_slot_center,
                    "iii_slot_length_aa": iii_slot_len,
                    "iii_slot_coordinate_method": iii_slot_method,
                    "iii_slot_coordinate_confidence": iii_slot_conf,
                    "iii_slot_coordinate_note": iii_slot_note,
                    "cds_phase_value": first(ch, ["phase", "phase_values"]),
                    "cds_boundary_precision": cds_precision,
                    "legacy_cds_boundary_precision": cds_precision,
                    "cds_boundary_precision_note": cds_precision_note,
                    "cds_left_phase": cds_lr["cds_left_phase"],
                    "cds_right_phase": cds_lr["cds_right_phase"],
                    "cds_left_boundary_precision": cds_lr["cds_left_boundary_precision"],
                    "cds_right_boundary_precision": cds_lr["cds_right_boundary_precision"],
                    "cds_boundary_precision_refined": cds_lr["cds_boundary_precision_refined"],
                    "cds_boundary_confidence": cds_lr["cds_boundary_confidence"],
                    "cds_phase_source": cds_lr["cds_phase_source"],
                    "cds_phase_warning": cds_lr["cds_phase_warning"],
                    "exon_id_source": first(ch, ["coding_exon_key", "exon_id_source", "exon_id", "raw_cds_ids", "cds_id", "ID"]),
                    "exon_rank": first(ch, ["matched_exon_rank", "exon_rank", "display_coding_exon_index", "cds_rank"]),
                    "chrom": first(ch, ["chrom", "seqid", "seq_id", "contig", "accession"]),
                    "exon_start": first(ch, ["coding_exon_cds_start", "cds_start", "start", "genomic_start", "exon_start"]),
                    "exon_end": first(ch, ["coding_exon_cds_end", "cds_end", "end", "genomic_end", "exon_end"]),
                    "strand": first(ch, ["strand"]),
                    "phase": first(ch, ["phase", "phase_values"]),
                    "exon_protein_start_aa": ch.get("candidate_protein_start_aa", ""),
                    "exon_protein_end_aa": ch.get("candidate_protein_end_aa", ""),
                    "exon_coordinate_status": status,
                    "coordinate_source": first(ch, ["coordinate_source"]),
                    "coordinate_warning": first(ch, ["warning"]) + ("; " if first(ch, ["warning"]) else "") + f"resolver_score={ch.get('candidate_score','')}; candidate_status={ch.get('candidate_status','')}",
                    "strict_cds_exon_key": first(ch, ["strict_cds_exon_key"]),
                    "strict_pair_class": pair_class,
                    "main_figure_exact_bool": bool(exact),
                    "resolver_version": VERSION,
                    "resolver_pair_same_key": int(same_pair_key),
                    "resolver_pair_has_both_isoforms": int(pair_has_both),
                    "resolver_structure_overlap_bp": structure_ov,
                    "resolver_window_overlap_aa": ov,
                    "resolver_transcript_match": int(txm),
                    "resolver_protein_match": int(prm),
                })
                resolved_rows.append(base)
    resolved = pd.DataFrame(resolved_rows)
    candidates = pd.concat(all_candidates, ignore_index=True) if all_candidates else pd.DataFrame()
    if unresolved:
        _unresolved_df = pd.DataFrame(unresolved)
    else:
        _unresolved_df = pd.DataFrame()
    return resolved, candidates


# Conservative pair-level coordinate sanity and figure eligibility.
# These QC categories used to live in the plotting script. They are defined here
# (data layer) so plotting only visualises existing QC fields.

_EVIDENCE_RANK = {
    "transcript_and_protein": 4,
    "transcript_only": 3,
    "protein_only": 3,
    "species_level_fallback": 2,
    "incomplete_identity": 1,
    "unresolved": 0,
    "": 0,
}

_CONF_RANK = {"high": 3, "moderate": 2, "low": 1, "unresolved": 0, "": 0}


def _native_coordinate_sanity(dist: Optional[float]) -> str:
    if dist is None:
        return "native_coordinate_unresolved"
    if dist <= 5:
        return "same_native_coordinate"
    if dist <= 20:
        return "minor_native_offset"
    if dist <= 60:
        return "moderate_native_offset_review"
    return "major_native_offset_review"


def _iii_slot_coordinate_sanity(dist: Optional[float], confidence_ok: bool) -> str:
    if dist is None or not confidence_ok:
        return "normalized_III_slot_unresolved"
    if dist <= 5:
        return "same_normalized_III_slot"
    if dist <= 15:
        return "minor_normalized_III_slot_offset"
    return "normalized_III_slot_review"




def build_pair_qc(resolved: pd.DataFrame) -> pd.DataFrame:
    """Return one row per species with native and normalized III-slot pair QC."""
    if resolved.empty:
        return pd.DataFrame()
    sp_col = col_first(resolved, ["species_canonical", "species"])
    rows: List[dict] = []
    for sp, g in resolved.groupby(resolved[sp_col].astype(str)):
        row: Dict[str, Any] = {"species_canonical": str(sp)}
        row["n_rows"] = len(g)
        row["statuses"] = ";".join(sorted(set(g.get("legacy_resolver_status", pd.Series(dtype=str)).astype(str))))
        row["resolver_status_refined_set"] = ";".join(sorted(set(g.get("resolver_status_refined", pd.Series(dtype=str)).astype(str))))
        row["has_both_isoforms"] = set(g.get("inferred_isoform", pd.Series(dtype=str)).astype(str)) >= {"IIIb", "IIIc"}
        row["pair_audit_status"] = first(g.iloc[0], ["pair_audit_status"]) if len(g) else ""
        row["iii_region_similarity_class"] = first(g.iloc[0], ["iii_region_similarity_class"]) if len(g) else ""
        per_iso: Dict[str, pd.Series] = {}
        for iso in ["IIIb", "IIIc"]:
            gi = g[g.get("inferred_isoform", pd.Series(dtype=str)).astype(str).eq(iso)]
            if gi.empty:
                for suffix in ["status", "refined", "evidence_level", "confidence",
                               "native_start_aa", "native_end_aa", "native_center_aa", "native_len_aa",
                               "iii_slot_center_aa", "iii_slot_len_aa", "iii_slot_confidence",
                               "strict_cds_exon_key", "transcript_id_source", "protein_id",
                               "window_overlap_aa", "cds_boundary_precision",
                               "cds_boundary_precision_refined", "cds_left_boundary_precision",
                               "cds_right_boundary_precision", "cds_boundary_confidence"]:
                    row[f"{iso}_{suffix}"] = "missing" if suffix in ("status", "refined") else ""
                continue
            # Prefer the strongest-evidence row.
            gi = gi.assign(_rank=gi.get("resolver_evidence_level", pd.Series(dtype=str)).astype(str).map(lambda v: _EVIDENCE_RANK.get(v, 0)))
            r = gi.sort_values("_rank", ascending=False).iloc[0]
            per_iso[iso] = r
            row[f"{iso}_status"] = str(r.get("legacy_resolver_status", ""))
            row[f"{iso}_refined"] = str(r.get("resolver_status_refined", ""))
            row[f"{iso}_evidence_level"] = str(r.get("resolver_evidence_level", ""))
            row[f"{iso}_confidence"] = str(r.get("resolver_confidence", ""))
            row[f"{iso}_native_start_aa"] = first(r, ["native_protein_start_aa"])
            row[f"{iso}_native_end_aa"] = first(r, ["native_protein_end_aa"])
            row[f"{iso}_native_center_aa"] = first(r, ["native_protein_center_aa"])
            row[f"{iso}_native_len_aa"] = first(r, ["native_protein_length_aa"])
            row[f"{iso}_iii_slot_center_aa"] = first(r, ["iii_slot_center_aa"])
            row[f"{iso}_iii_slot_len_aa"] = first(r, ["iii_slot_length_aa"])
            row[f"{iso}_iii_slot_confidence"] = first(r, ["iii_slot_coordinate_confidence"])
            row[f"{iso}_strict_cds_exon_key"] = first(r, ["strict_cds_exon_key"])
            row[f"{iso}_transcript_id_source"] = first(r, ["transcript_id_source"])
            row[f"{iso}_protein_id"] = first(r, ["protein_id"])
            row[f"{iso}_window_overlap_aa"] = first(r, ["resolver_window_overlap_aa"])
            row[f"{iso}_cds_boundary_precision"] = first(r, ["cds_boundary_precision"])
            row[f"{iso}_cds_boundary_precision_refined"] = first(r, ["cds_boundary_precision_refined"])
            row[f"{iso}_cds_left_boundary_precision"] = first(r, ["cds_left_boundary_precision"])
            row[f"{iso}_cds_right_boundary_precision"] = first(r, ["cds_right_boundary_precision"])
            row[f"{iso}_cds_boundary_confidence"] = first(r, ["cds_boundary_confidence"])
            # Back-compat coordinate columns used by existing figures.
            row[f"{iso}_start_aa"] = first(r, ["native_protein_start_aa"])
            row[f"{iso}_end_aa"] = first(r, ["native_protein_end_aa"])

        both_resolved = ("IIIb" in per_iso) and ("IIIc" in per_iso)
        # Native pair center distance.
        try:
            nb = float(row.get("IIIb_native_center_aa")) if row.get("IIIb_native_center_aa") not in ("", None) else None
            nc = float(row.get("IIIc_native_center_aa")) if row.get("IIIc_native_center_aa") not in ("", None) else None
        except Exception:
            nb = nc = None
        native_dist = abs(nb - nc) if (nb is not None and nc is not None) else None
        # Normalized III-slot pair center distance.
        try:
            sb = float(row.get("IIIb_iii_slot_center_aa")) if row.get("IIIb_iii_slot_center_aa") not in ("", None) else None
            sc = float(row.get("IIIc_iii_slot_center_aa")) if row.get("IIIc_iii_slot_center_aa") not in ("", None) else None
        except Exception:
            sb = sc = None
        slot_dist = abs(sb - sc) if (sb is not None and sc is not None) else None

        row["native_pair_center_distance_aa"] = round(native_dist, 1) if native_dist is not None else ""
        row["iii_slot_pair_center_distance_aa"] = round(slot_dist, 1) if slot_dist is not None else ""

        # The normalized III-slot is an exon_internal_relative (cassette shape/length)
        # axis with intentionally capped confidence (<= moderate). conf_ok therefore
        # only requires that BOTH cassettes were normalized at all (>= low, i.e. not
        # unresolved); the absolute-position guard is native_coordinate_sanity.
        conf_ok = (_CONF_RANK.get(str(row.get("IIIb_iii_slot_confidence", "")), 0) >= 1 and
                   _CONF_RANK.get(str(row.get("IIIc_iii_slot_confidence", "")), 0) >= 1)
        row["native_coordinate_sanity"] = _native_coordinate_sanity(native_dist)
        row["iii_slot_coordinate_sanity"] = _iii_slot_coordinate_sanity(slot_dist, conf_ok)
        # Legacy column retained for existing plotting sort.
        row["pair_coordinate_sanity"] = row["native_coordinate_sanity"]

        distinct = bool(row.get("IIIb_strict_cds_exon_key")) and bool(row.get("IIIc_strict_cds_exon_key")) and row.get("IIIb_strict_cds_exon_key") != row.get("IIIc_strict_cds_exon_key")
        row["distinct_CDS_pair"] = bool(distinct)

        # Pair precision uses the weaker boundary result.
        _prec_rank = {"exact": 3, "codon_split_one_side": 2, "codon_split_both_sides": 1, "unknown_codon_phase": 0, "": 0, "missing": 0}
        pb = str(row.get("IIIb_cds_boundary_precision_refined", "") or "")
        pc = str(row.get("IIIc_cds_boundary_precision_refined", "") or "")
        if both_resolved and pb and pc and pb != "missing" and pc != "missing":
            worst = min((pb, pc), key=lambda v: _prec_rank.get(v, 0))
            row["cds_boundary_precision_summary"] = worst
            if pb == "exact" and pc == "exact":
                row["cds_boundary_confidence_summary"] = "high"
            elif _prec_rank.get(worst, 0) == 0:
                row["cds_boundary_confidence_summary"] = "low"
            else:
                row["cds_boundary_confidence_summary"] = "medium"
        else:
            row["cds_boundary_precision_summary"] = "unknown_codon_phase"
            row["cds_boundary_confidence_summary"] = "low"

        # Pair-level evidence level (weakest of the two).
        eb = _EVIDENCE_RANK.get(str(row.get("IIIb_evidence_level", "")), 0)
        ec = _EVIDENCE_RANK.get(str(row.get("IIIc_evidence_level", "")), 0)
        row["pair_min_evidence_rank"] = min(eb, ec) if both_resolved else 0
        both_transcript_specific = both_resolved and eb >= 3 and ec >= 3

        # Legacy-style pair_status for backward compatibility with plots.
        both_gold = both_resolved and str(row.get("IIIb_refined", "")).startswith(("gold", "silver")) and str(row.get("IIIc_refined", "")).startswith(("gold", "silver"))
        if both_gold and distinct:
            row["pair_status"] = "gold_exact_distinct_CDS_pair"
        elif both_resolved and (str(row.get("IIIb_refined", "")).startswith(("gold", "silver", "bronze")) or str(row.get("IIIc_refined", "")).startswith(("gold", "silver", "bronze"))):
            row["pair_status"] = "partial_or_review"
        else:
            row["pair_status"] = "review_no_exact_pair"

        # Anchor contradiction uses the refined
        # similarity class when present: only FULL-window near-identity (or a
        # failed human control / identical proteins) is treated as contradictory.
        # Window-distinct-but-local-subregion-identical pairs are NOT contradictory.
        pa_status = str(row.get("pair_audit_status", ""))
        sim_class = str(row.get("iii_region_similarity_class", "")).strip()
        if sim_class and sim_class not in ("", "nan"):
            anchor_contradictory = (
                sim_class in ("full_window_nearly_identical",
                              "window_distinct_but_local_subregion_identical",
                              "unresolved")
                or "human_positive_control_failed" in pa_status
                or "full_proteins_identical" in pa_status
            )
        else:
            anchor_contradictory = ("human_positive_control_failed" in pa_status) or ("identical" in pa_status) or pa_status.startswith("same_")

        # Review flags.
        flags: List[str] = []
        if not both_resolved:
            flags.append("missing_resolved_isoform")
        if not both_transcript_specific:
            flags.append("evidence_below_transcript_specific")
        if not distinct:
            flags.append("non_distinct_CDS_pair")
        if anchor_contradictory:
            flags.append(f"anchor_contradictory:{sim_class or pa_status}")
        elif sim_class == "ambiguous_similarity_review":
            flags.append("iii_region_ambiguous_similarity_review")
        elif sim_class == "window_distinct_but_short_local_identity":
            flags.append("iii_region_short_local_identity_review")
        if row["native_coordinate_sanity"] == "moderate_native_offset_review":
            flags.append("moderate_native_offset")
        if row["native_coordinate_sanity"] == "major_native_offset_review":
            flags.append("major_native_offset")
        if row["iii_slot_coordinate_sanity"] == "normalized_III_slot_review":
            flags.append("normalized_slot_offset")
        if row["iii_slot_coordinate_sanity"] == "normalized_III_slot_unresolved":
            flags.append("normalized_slot_unresolved")
        if "review_identity_incomplete" in row.get("resolver_status_refined_set", ""):
            flags.append("identity_incomplete")

        # Main-figure eligibility.
        slot_ok = row["iii_slot_coordinate_sanity"] in ("same_normalized_III_slot", "minor_normalized_III_slot_offset")
        native_not_major = row["native_coordinate_sanity"] != "major_native_offset_review"
        eligible = (
            both_transcript_specific
            and distinct
            and not anchor_contradictory
            and conf_ok
            and slot_ok
            and native_not_major
        )
        row["main_figure_eligible"] = bool(eligible)
        if eligible:
            reason_bits = ["both_transcript_specific", "distinct_CDS_pair",
                           row["iii_slot_coordinate_sanity"], row["native_coordinate_sanity"]]
            if row["native_coordinate_sanity"] == "moderate_native_offset_review":
                reason_bits.append("moderate_native_offset_explained_by_normalized_slot")
            row["main_figure_eligibility_reason"] = ";".join(reason_bits)
            row["supplement_only_reason"] = ""
        else:
            row["main_figure_eligibility_reason"] = ""
            supp = []
            if not both_resolved:
                supp.append("missing_resolved_isoform")
            elif not both_transcript_specific:
                supp.append("evidence_below_transcript_specific")
            if not distinct:
                supp.append("non_distinct_CDS_pair")
            if anchor_contradictory:
                supp.append("anchor_contradictory")
            if not conf_ok or not slot_ok:
                supp.append(row["iii_slot_coordinate_sanity"])
            if not native_not_major:
                supp.append("major_native_offset_review")
            row["supplement_only_reason"] = ";".join(supp) if supp else "review"
        row["review_flags"] = ";".join(flags)
        # Legacy review_reason column (plots reference it).
        row["review_reason"] = row["supplement_only_reason"]
        rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["main_figure_eligible", "native_coordinate_sanity", "species_canonical"],
                           ascending=[False, True, True]).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser(description="Resolve FGFR2 IIIb/IIIc exact CDS/exon mapping using multi-evidence scoring.")
    ap.add_argument("--pair_audit", required=True, type=Path)
    ap.add_argument("--cds_features", required=True, type=Path)
    ap.add_argument("--selected", required=True, type=Path)
    ap.add_argument("--alt_exons", required=False, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--prefix", default="fgfr2")
    args = ap.parse_args()

    pair = read_tsv(args.pair_audit)
    cds = read_tsv(args.cds_features)
    selected = read_tsv(args.selected)
    alt = read_tsv(args.alt_exons, required=False) if args.alt_exons else pd.DataFrame()
    resolved, candidates = build_resolved(pair, cds, selected, alt)

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    # Deterministic ordering by species, isoform, transcript.
    if not resolved.empty:
        sp_col = col_first(resolved, ["species_canonical", "species"]) or "species_canonical"
        sort_cols = [c for c in [sp_col, "inferred_isoform", "transcript_id_source"] if c in resolved.columns]
        if sort_cols:
            resolved = resolved.sort_values(sort_cols).reset_index(drop=True)
    write_tsv(resolved, out / f"{args.prefix}_resolved_IIIb_IIIc_exon_CDS_mapping.tsv")
    write_tsv(candidates, out / f"{args.prefix}_resolved_IIIb_IIIc_candidate_scores.tsv")
    # Compatibility name for existing plotting/QC expectations.
    write_tsv(resolved, out / f"{args.prefix}_current_stage_IIIb_IIIc_coordinate_audit.tsv")

    pair_qc = build_pair_qc(resolved)
    write_tsv(pair_qc, out / f"{args.prefix}_pair_level_qc_summary.tsv")

    status_counts = resolved.get("exon_coordinate_status", pd.Series(dtype=str)).value_counts().rename_axis("status").reset_index(name="n")
    write_tsv(status_counts, out / f"{args.prefix}_resolved_status_counts.tsv")
    refined_counts = resolved.get("resolver_status_refined", pd.Series(dtype=str)).value_counts().rename_axis("resolver_status_refined").reset_index(name="n")
    write_tsv(refined_counts, out / f"{args.prefix}_resolver_status_refined_counts.tsv")
    def _vc(col: str) -> Dict[str, int]:
        if resolved.empty or col not in resolved.columns:
            return {}
        return {str(k): int(v) for k, v in resolved[col].value_counts().to_dict().items()}

    def _vc_pair(col: str) -> Dict[str, int]:
        if pair_qc.empty or col not in pair_qc.columns:
            return {}
        return {str(k): int(v) for k, v in pair_qc[col].value_counts().to_dict().items()}

    meta = {
        "version": VERSION,
        "rows_resolved": int(len(resolved)),
        "rows_candidates": int(len(candidates)),
        "species_in_pair_qc": int(len(pair_qc)),
        "main_figure_eligible": int(pair_qc.get("main_figure_eligible", pd.Series(dtype=bool)).sum()) if not pair_qc.empty else 0,
        "resolver_status_refined_counts": refined_counts.set_index("resolver_status_refined")["n"].to_dict() if not refined_counts.empty else {},
        "cds_boundary_precision": {
            "legacy_cds_boundary_precision_counts_before": _vc("legacy_cds_boundary_precision"),
            "cds_boundary_precision_refined_counts_after": _vc("cds_boundary_precision_refined"),
            "cds_left_boundary_precision_counts": _vc("cds_left_boundary_precision"),
            "cds_right_boundary_precision_counts": _vc("cds_right_boundary_precision"),
            "cds_boundary_confidence_counts": _vc("cds_boundary_confidence"),
            "cds_phase_source_counts": _vc("cds_phase_source"),
            "pair_cds_boundary_precision_summary_counts": _vc_pair("cds_boundary_precision_summary"),
        },
        "inputs": {"pair_audit": str(args.pair_audit), "cds_features": str(args.cds_features), "selected": str(args.selected), "alt_exons": str(args.alt_exons or "")},
    }
    (out / f"{args.prefix}_resolver_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"[OK] wrote resolved mapping: {out / (args.prefix + '_resolved_IIIb_IIIc_exon_CDS_mapping.tsv')}")
    print(f"[OK] wrote pair-level QC: {out / (args.prefix + '_pair_level_qc_summary.tsv')}")
    if not refined_counts.empty:
        print(refined_counts.to_string(index=False))


if __name__ == "__main__":
    main()
