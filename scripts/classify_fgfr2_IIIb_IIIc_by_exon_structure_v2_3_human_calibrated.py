#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.shared_gene_analysis.strand import is_reverse, same_strand  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import protein_lookup as PL  # noqa: E402

SCRIPT_NAME = "classify_fgfr2_IIIb_IIIc_by_exon_structure_v2.py"
SCRIPT_VERSION = "2.4.0_sequence_calibrated_direction"

@dataclass
class AlignMetrics:
    score: int = 0
    identity: float = 0.0        # matched columns / aligned columns, [0,1]
    coverage_query: float = 0.0  # aligned query span / query length, [0,1]
    aln_len: int = 0
    n_match: int = 0


def smith_waterman_local(query: str, ref: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> AlignMetrics:
    q = (query or "").upper()
    t = (ref or "").upper()
    n, m = len(q), len(t)
    if n == 0 or m == 0:
        return AlignMetrics()
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best = bi = bj = 0
    for i in range(1, n + 1):
        qi = q[i - 1]
        Hi = H[i]
        Hp = H[i - 1]
        for j in range(1, m + 1):
            s = match if qi == t[j - 1] else mismatch
            v = Hp[j - 1] + s
            d = Hp[j] + gap
            if d > v:
                v = d
            e = Hi[j - 1] + gap
            if e > v:
                v = e
            if v < 0:
                v = 0
            Hi[j] = v
            if v > best:
                best = v
                bi, bj = i, j
    if best == 0:
        return AlignMetrics()
    i, j = bi, bj
    q_end = i
    n_match = aln_len = 0
    while i > 0 and j > 0 and H[i][j] > 0:
        cur = H[i][j]
        diag = H[i - 1][j - 1] + (match if q[i - 1] == t[j - 1] else mismatch)
        if cur == diag:
            aln_len += 1
            if q[i - 1] == t[j - 1]:
                n_match += 1
            i -= 1
            j -= 1
        elif cur == H[i - 1][j] + gap:
            i -= 1
            aln_len += 1
        else:
            j -= 1
            aln_len += 1
    q_start = i + 1
    identity = (n_match / aln_len) if aln_len else 0.0
    coverage_query = ((q_end - q_start + 1) / n) if n else 0.0
    return AlignMetrics(
        score=best,
        identity=round(min(1.0, max(0.0, identity)), 6),
        coverage_query=round(min(1.0, max(0.0, coverage_query)), 6),
        aln_len=aln_len,
        n_match=n_match,
    )


def read_single_fasta(path: Path) -> str:
    seq_parts: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if seq_parts:
                    break
                continue
            seq_parts.append(re.sub(r"[^A-Za-z]", "", line))
    return "".join(seq_parts).upper()


def clean_id(value: str) -> str:
    v = (value or "").strip()
    for prefix in ("rna-", "cds-", "gene-", "id-"):
        if v.startswith(prefix):
            v = v[len(prefix):]
    return v


def read_tsv(path: Optional[Path]) -> List[Dict[str, str]]:
    if path is None:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [dict(r) for r in reader]


def write_tsv(rows: List[Dict[str, object]], path: Path, fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def get_any(row: Dict[str, str], names: Iterable[str], default: str = "") -> str:
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return str(value)
    return default


def to_int(value: object) -> Optional[int]:
    try:
        if value is None or value == "" or str(value).lower() == "nan":
            return None
        return int(float(str(value)))
    except Exception:
        return None


def to_float(value: object) -> Optional[float]:
    try:
        if value is None or value == "" or str(value).lower() == "nan":
            return None
        return float(str(value))
    except Exception:
        return None


def truthy(value: object) -> bool:
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "selected", "canonical"}


TRANSCRIPT_REQUIRED_ALIASES: Dict[str, Sequence[str]] = {
    "internal_transcript_id": ("internal_transcript_id", "transcript_id_internal", "tx_internal_id"),
    "transcript_id_source": ("transcript_id_source", "transcript_id", "id", "accession"),
    "species": ("species_canonical", "species", "species_input"),
    "source_db": ("source_db", "source"),
    "gene_id": ("gene_id_internal", "gene_id", "gene"),
}

EXON_REQUIRED_ALIASES: Dict[str, Sequence[str]] = {
    "transcript_id": ("transcript_id_internal", "internal_transcript_id", "parent_transcript_id", "transcript_id", "transcript_id_source", "tx_internal_id"),
    "rank": ("exon_rank", "rank", "exon_number"),
    "chrom": ("chrom", "seqid", "chromosome"),
    "start": ("start",),
    "end": ("end",),
    "strand": ("strand",),
}


def validate_alias_schema(rows: List[Dict[str, str]], aliases: Dict[str, Sequence[str]], label: str) -> List[Dict[str, object]]:
    warnings: List[Dict[str, object]] = []
    if not rows:
        warnings.append({"level": "error", "warning_type": "empty_input", "file": label, "message": f"{label} contains no rows"})
        return warnings
    fieldnames = set(rows[0].keys())
    for logical_name, possible_names in aliases.items():
        if not any(name in fieldnames for name in possible_names):
            warnings.append({
                "level": "error",
                "warning_type": "missing_required_column",
                "file": label,
                "logical_column": logical_name,
                "accepted_aliases": ",".join(possible_names),
                "message": f"{label} lacks required logical column '{logical_name}'",
            })
    return warnings



@dataclass(frozen=True)
class Exon:
    tx_internal: str
    exon_id_source: str
    rank: int
    chrom: str
    start: int
    end: int
    strand: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def sig(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}:{self.strand}"


def transcript_internal_id(row: Dict[str, str]) -> str:
    return get_any(row, ["internal_transcript_id", "transcript_id_internal", "tx_internal_id"])


def transcript_source_id(row: Dict[str, str]) -> str:
    return get_any(row, ["transcript_id_source", "transcript_id", "id", "accession"])


def exon_transcript_id(row: Dict[str, str]) -> str:
    return get_any(row, ["transcript_id_internal", "internal_transcript_id", "parent_transcript_id", "transcript_id", "transcript_id_source", "tx_internal_id"])


def group_key_from_tx(row: Dict[str, str]) -> Tuple[str, str, str, str]:
    species = get_any(row, ["species_canonical", "species", "species_input"], "unknown_species")
    source = get_any(row, ["source_db", "source"], "unknown_source")
    gene = get_any(row, ["gene_id_internal", "gene_id", "gene"], "unknown_gene")
    symbol = get_any(row, ["gene_symbol", "symbol"], "FGFR2")
    return species, source, gene, symbol


def load_exons(exon_rows: List[Dict[str, str]], warnings: Optional[List[Dict[str, object]]] = None) -> Dict[str, List[Exon]]:
    by_tx: Dict[str, List[Exon]] = defaultdict(list)
    for idx, r in enumerate(exon_rows, start=1):
        tx = exon_transcript_id(r)
        if not tx:
            if warnings is not None:
                warnings.append({"level": "warning", "warning_type": "exon_without_transcript_id", "row_number": idx, "message": "Exon row lacks a transcript identifier"})
            continue
        rank = to_int(get_any(r, ["exon_rank", "rank", "exon_number"]))
        start = to_int(r.get("start"))
        end = to_int(r.get("end"))
        chrom = get_any(r, ["chrom", "seqid", "chromosome"])
        strand = get_any(r, ["strand"], "+") or "+"
        if rank is None or start is None or end is None or not chrom:
            if warnings is not None:
                warnings.append({"level": "warning", "warning_type": "invalid_exon_row", "row_number": idx, "transcript_id": tx, "message": "Exon row lacks rank/start/end/chrom and was ignored"})
            continue
        if start > end:
            start, end = end, start
        if strand not in {"+", "-"}:
            strand = "+"
        by_tx[tx].append(Exon(tx, get_any(r, ["exon_id_source", "exon_id", "id"]), rank, chrom, start, end, strand))
    for tx in list(by_tx):
        by_tx[tx] = sorted(by_tx[tx], key=lambda e: e.rank)
    return by_tx


def normalize_biotype(value: str) -> str:
    text = (value or "").lower().replace("-", "_").replace(" ", "_")
    if not text:
        return "unknown"
    if any(x in text for x in ["pseudogene"]):
        return "pseudogene"
    if any(x in text for x in ["nonsense_mediated_decay", "nmd", "retained_intron"]):
        return "noncoding_or_nmd"
    if any(x in text for x in ["protein_coding", "mrna", "messenger_rna", "coding", "protein"]):
        return "coding"
    if any(x in text for x in ["lncrna", "ncrna", "mirna", "rrna", "trna"]):
        return "noncoding_or_nmd"
    return "unknown"


def is_complete_candidate_tx(
    tx: Dict[str, str],
    exons: List[Exon],
    min_protein_aa: int,
    max_protein_aa: int,
) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    if len(exons) < 8:
        reasons.append("too_few_exons")
    biotype = normalize_biotype(get_any(tx, ["transcript_biotype", "biotype"]))
    if biotype not in {"coding", "unknown"}:
        reasons.append(f"noncoding_biotype:{biotype}")
    length = to_float(get_any(tx, ["protein_length_aa", "protein_length", "translation_length"]))
    if length is not None and not (min_protein_aa <= length <= max_protein_aa):
        reasons.append("implausible_protein_length")
    return len(reasons) == 0, reasons


def transcript_order_position(exon: Exon) -> int:
    return -exon.start if is_reverse(exon.strand) else exon.start


def build_slot_candidates(
    tx_rows: List[Dict[str, str]],
    exons_by_tx: Dict[str, List[Exon]],
    min_alt_exon_len: int,
    max_alt_exon_len: int,
    min_tx_fraction_with_flanks: float,
) -> List[Dict[str, object]]:
    tx_ids = [transcript_internal_id(t) for t in tx_rows]
    tx_ids = [t for t in tx_ids if t in exons_by_tx]
    n_tx = max(1, len(tx_ids))
    slot_map: Dict[Tuple[str, str], Dict[str, object]] = {}
    for tx_id in tx_ids:
        exons = exons_by_tx.get(tx_id, [])
        for i in range(1, len(exons) - 1):
            alt = exons[i]
            if not (min_alt_exon_len <= alt.length <= max_alt_exon_len):
                continue
            left = exons[i - 1]
            right = exons[i + 1]
            if not (left.chrom == alt.chrom == right.chrom
                and same_strand(left.strand, alt.strand, right.strand)):
                continue
            key = (left.sig, right.sig)
            if key not in slot_map:
                slot_map[key] = {"left_flank": left, "right_flank": right, "alt_exons": {}, "tx_to_alt": defaultdict(set), "tx_seen": set()}
            slot_map[key]["alt_exons"].setdefault(alt.sig, alt)
            slot_map[key]["tx_to_alt"][tx_id].add(alt.sig)
            slot_map[key]["tx_seen"].add(tx_id)

    candidates: List[Dict[str, object]] = []
    for idx, ((left_sig, right_sig), obj) in enumerate(slot_map.items(), start=1):
        alt_exons: Dict[str, Exon] = obj["alt_exons"]
        tx_to_alt: Dict[str, set] = obj["tx_to_alt"]
        tx_seen: set = obj["tx_seen"]
        if len(alt_exons) < 2:
            continue
        tx_fraction = len(tx_seen) / n_tx
        if tx_fraction < min_tx_fraction_with_flanks:
            continue
        multi_alt_tx = sum(1 for s in tx_to_alt.values() if len(s) > 1)
        mutual_exclusivity = 1.0 - (multi_alt_tx / max(1, len(tx_to_alt)))
        lengths = [e.length for e in alt_exons.values()]
        length_spread = max(lengths) - min(lengths)
        score = len(tx_seen) * 10 + len(alt_exons) * 25 + int(mutual_exclusivity * 100) - length_spread - max(0, len(alt_exons) - 2) * 20
        candidates.append({
            "slot_id": f"slot_{idx}",
            "left_flank_sig": left_sig,
            "right_flank_sig": right_sig,
            "alt_count": len(alt_exons),
            "tx_with_slot": len(tx_seen),
            "tx_total": n_tx,
            "tx_fraction": round(tx_fraction, 4),
            "mutual_exclusivity": round(mutual_exclusivity, 4),
            "length_spread": length_spread,
            "score": score,
            "alt_exons": alt_exons,
            "tx_to_alt": tx_to_alt,
        })
    candidates.sort(key=lambda c: (c["score"], c["tx_with_slot"], -c["length_spread"]), reverse=True)
    return candidates


def slot_geometry_features(slot: Dict[str, object]) -> Dict[str, object]:
    alt_exons: Dict[str, Exon] = slot.get("alt_exons", {})
    alts = sorted(alt_exons.values(), key=transcript_order_position)
    out: Dict[str, object] = {
        "alt_pair_nonoverlapping": 0,
        "alt_pair_min_gap_bp": "",
        "alt_pair_order_distance_bp": "",
        "alt_pair_geometry_class": "not_exact_two_alt_exons",
        "fgfr2_iii_geometry_score": 0,
    }
    if len(alts) != 2:
        return out
    a, b = alts
    if a.chrom != b.chrom or not same_strand(a.strand, b.strand):
        out["alt_pair_geometry_class"] = "different_chrom_or_strand"
        return out
    overlap = max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
    if overlap > 0:
        gap = -overlap
    else:
        gap = max(a.start, b.start) - min(a.end, b.end) - 1
    order_distance = abs(transcript_order_position(b) - transcript_order_position(a))
    length_similarity_penalty = abs(a.length - b.length)
    nonoverlap = 1 if gap >= 1 else 0
    if gap >= 50:
        geom_class = "separated_cassette_exons_consistent_with_FGFR2_IIIb_IIIc"
        geom_score = 1000 + min(gap, 5000) - length_similarity_penalty
    elif gap >= 1:
        geom_class = "nonoverlapping_but_close_cassette_exons"
        geom_score = 500 + gap - length_similarity_penalty
    else:
        geom_class = "overlapping_alternative_splice_site_like_exons_not_preferred_for_IIIb_IIIc"
        geom_score = -1000 + gap - length_similarity_penalty
    out.update({
        "alt_pair_nonoverlapping": nonoverlap,
        "alt_pair_min_gap_bp": gap,
        "alt_pair_order_distance_bp": order_distance,
        "alt_pair_geometry_class": geom_class,
        "fgfr2_iii_geometry_score": geom_score,
    })
    return out


def choose_fgfr2_slot(candidates: List[Dict[str, object]]) -> Optional[Dict[str, object]]:
    if not candidates:
        return None
    for c in candidates:
        c.update(slot_geometry_features(c))
    exact_two = [c for c in candidates if c.get("alt_count") == 2 and float(c.get("mutual_exclusivity", 0)) >= 0.95]
    if exact_two:
        exact_two.sort(key=lambda c: (
            int(c.get("fgfr2_iii_geometry_score", 0)),
            int(c.get("alt_pair_nonoverlapping", 0)),
            int(c.get("tx_with_slot", 0)),
            float(c.get("score", 0)),
        ), reverse=True)
        return exact_two[0]
    candidates.sort(key=lambda c: (
        int(c.get("fgfr2_iii_geometry_score", 0)) if "fgfr2_iii_geometry_score" in c else 0,
        float(c.get("score", 0)),
        int(c.get("tx_with_slot", 0)),
    ), reverse=True)
    return candidates[0]


def slot_confidence_and_reason(slot: Dict[str, object]) -> Tuple[str, str]:
    alt_count = int(slot.get("alt_count", 0))
    mutual = float(slot.get("mutual_exclusivity", 0))
    if alt_count == 2 and mutual >= 0.95:
        return "high", "exactly_two_alternative_exons_and_high_mutual_exclusivity"
    if alt_count >= 2 and mutual >= 0.75:
        return "medium", "candidate_slot_detected_but_alt_count_or_mutual_exclusivity_not_ideal"
    return "low", "weak_candidate_slot_detected;manual_review_recommended"


LEGACY_DIRECTION_RULE = ("order_rule:first_alternative_exon=IIIb;second_alternative_exon=IIIc "
                         "(transcript_order; provisional, requires sequence/domain validation)")


@dataclass
class DirectionContext:
    cache_index: object
    protein_cache: Dict[str, object]
    cds_by_tx: Dict[str, List[Dict[str, int]]]
    tx_by_internal: Dict[str, Dict[str, str]]
    iiib_ref: str
    iiic_ref: str
    args: argparse.Namespace
    enabled: bool


def load_cds_by_tx(cds_rows: List[Dict[str, str]]) -> Dict[str, List[Dict[str, int]]]:
    out: Dict[str, List[Dict[str, int]]] = defaultdict(list)
    for r in cds_rows:
        tx = get_any(r, ["transcript_id_internal", "internal_transcript_id", "transcript_id_source"])
        if not tx:
            continue
        start = to_int(r.get("start"))
        end = to_int(r.get("end"))
        ps = to_int(r.get("protein_start_aa"))
        pe = to_int(r.get("protein_end_aa"))
        if start is None or end is None or ps is None or pe is None:
            continue
        if start > end:
            start, end = end, start
        if ps > pe:
            ps, pe = pe, ps
        out[tx].append({"start": start, "end": end, "ps": ps, "pe": pe})
    return out


def get_protein_for_tx(ctx: DirectionContext, tx_internal: str) -> str:
    tx = ctx.tx_by_internal.get(tx_internal)
    if not tx:
        return ""
    translation = clean_id(get_any(tx, ["translation_id_source", "protein_id", "translation_id"]))
    transcript = clean_id(get_any(tx, ["transcript_id_source", "transcript_id"]))
    source_db = get_any(tx, ["source_db", "source"])
    key = "|".join([source_db, translation, transcript])
    if key in ctx.protein_cache:
        res = ctx.protein_cache[key]
    else:
        expected_len = to_int(get_any(tx, ["protein_length_aa"]))
        res = PL.lookup_protein(
            ctx.cache_index,
            source_db=source_db,
            transcript_id=transcript,
            translation_id=translation,
            species_input=get_any(tx, ["species_input"]),
            species_canonical=get_any(tx, ["species_canonical", "species"]),
            expected_length_aa=expected_len if expected_len else None,
            allow_ensembl_rest=not ctx.args.no_ensembl_rest,
            ensembl_sleep=ctx.args.ensembl_sleep,
            ensembl_timeout=ctx.args.ensembl_timeout,
        )
        ctx.protein_cache[key] = res
    return (getattr(res, "sequence", "") or "").upper().replace("*", "")


def extract_exon_aa(ctx: DirectionContext, tx_internal: str, exon: Exon) -> str:
    cds = ctx.cds_by_tx.get(tx_internal, [])
    if not cds:
        return ""
    best = None
    best_ov = 0
    for c in cds:
        ov = min(exon.end, c["end"]) - max(exon.start, c["start"]) + 1
        if ov > best_ov:
            best_ov = ov
            best = c
    if best is None or best_ov <= 0:
        return ""
    protein = get_protein_for_tx(ctx, tx_internal)
    if not protein:
        return ""
    ps = max(1, best["ps"])
    pe = min(len(protein), best["pe"])
    if pe < ps:
        return ""
    return protein[ps - 1:pe]


def _direction_score(m: AlignMetrics) -> float:
    return round(m.identity * m.coverage_query, 6)


def calibrate_exon(ctx: DirectionContext, exon: Exon, tx_ids: List[str]) -> Dict[str, object]:
    out: Dict[str, object] = {
        "prefers": "unresolved",
        "id_b": "", "id_c": "", "cov_b": "", "cov_c": "",
        "score_b": "", "score_c": "",
        "method": "unresolved", "n_evaluated": 0,
    }
    if not ctx.enabled:
        out["method"] = "order_rule_provisional"
        return out
    ordered_tx = sorted(set(tx_ids), key=lambda t: 0 if str(get_any(ctx.tx_by_internal.get(t, {}),
                        ["source_db", "source"])).upper().startswith("N") else 1)
    candidates = []
    n_eval = 0
    for tx_id in ordered_tx:
        if n_eval >= int(ctx.args.max_proteins_per_exon):
            break
        aa = extract_exon_aa(ctx, tx_id, exon)
        if not aa or len(aa) < 10:
            continue
        n_eval += 1
        mb = smith_waterman_local(aa, ctx.iiib_ref)
        mc = smith_waterman_local(aa, ctx.iiic_ref)
        sb, sc = _direction_score(mb), _direction_score(mc)
        candidates.append((abs(sb - sc), aa, mb, mc, sb, sc))
        if max(sb, sc) >= 0.60 and abs(sb - sc) >= float(ctx.args.direction_margin):
            break
    out["n_evaluated"] = n_eval
    if not candidates:
        out["method"] = "unresolved"
        return out
    _, aa, mb, mc, sb, sc = max(candidates, key=lambda x: x[0])
    margin = float(ctx.args.direction_margin)
    floor = float(ctx.args.direction_min_identity)
    if max(sb, sc) < floor:
        prefers = "unresolved"
    elif sb - sc >= margin:
        prefers = "IIIb"
    elif sc - sb >= margin:
        prefers = "IIIc"
    else:
        prefers = "ambiguous"
    out.update({
        "prefers": prefers,
        "id_b": round(mb.identity, 4), "id_c": round(mc.identity, 4),
        "cov_b": round(mb.coverage_query, 4), "cov_c": round(mc.coverage_query, 4),
        "score_b": round(sb, 4), "score_c": round(sc, 4),
        "method": "translated_CDS_sequence_calibrated",
    })
    return out


def complement_iso(iso: str) -> str:
    return "IIIc" if iso == "IIIb" else ("IIIb" if iso == "IIIc" else "")


def resolve_pair_direction(metA: Dict[str, object], metB: Dict[str, object]) -> Dict[str, str]:
    pa = str(metA.get("prefers", "unresolved"))
    pb = str(metB.get("prefers", "unresolved"))
    legacy_a, legacy_b = "IIIb", "IIIc"  # order rule
    res = {"A_iso": legacy_a, "B_iso": legacy_b,
           "method": "order_rule_provisional", "status": "direction_unresolved_no_sequence",
           "confidence": "low", "warning": ""}

    decisive = {"IIIb", "IIIc"}
    if pa in decisive and pb in decisive and pa != pb:
        res.update(A_iso=pa, B_iso=pb,
                   method="translated_CDS_sequence_calibrated",
                   confidence="high")
    elif pa in decisive and pb not in decisive:
        res.update(A_iso=pa, B_iso=complement_iso(pa),
                   method="translated_CDS_sequence_calibrated", confidence="medium",
                   warning="second_exon_inferred_by_complement")
    elif pb in decisive and pa not in decisive:
        res.update(A_iso=complement_iso(pb), B_iso=pb,
                   method="translated_CDS_sequence_calibrated", confidence="medium",
                   warning="first_exon_inferred_by_complement")
    elif pa in decisive and pb in decisive and pa == pb:
        res.update(A_iso=legacy_a, B_iso=legacy_b,
                   method="order_rule_provisional",
                   status="ambiguous_sequence_direction_review", confidence="low",
                   warning="both_exons_prefer_same_isoform_order_rule_retained")
        return res
    else:
        no_seq = (metA.get("method") in ("unresolved", "order_rule_provisional")
                  and metB.get("method") in ("unresolved", "order_rule_provisional"))
        res.update(A_iso=legacy_a, B_iso=legacy_b, method="order_rule_provisional",
                   status=("direction_unresolved_no_sequence" if no_seq
                           else "ambiguous_sequence_direction_review"),
                   confidence="low",
                   warning="insufficient_sequence_direction_evidence_order_rule_retained")
        return res

    # Determine pass vs inverted relative to the legacy order rule.
    if res["A_iso"] == legacy_a and res["B_iso"] == legacy_b:
        res["status"] = "sequence_calibrated_pass"
    else:
        res["status"] = "sequence_calibrated_inverted_from_order_rule"
    return res


def classify_with_slot(
    tx_rows: List[Dict[str, str]],
    exons_by_tx: Dict[str, List[Exon]],
    slot: Dict[str, object],
    confidence_floor: str,
    confidence_reason: str,
    ctx: Optional[DirectionContext] = None,
) -> List[Dict[str, object]]:
    alt_exons: Dict[str, Exon] = slot["alt_exons"]  # type: ignore[assignment]
    tx_to_alt: Dict[str, set] = slot["tx_to_alt"]   # type: ignore[assignment]
    ordered_alts = sorted(alt_exons.values(), key=transcript_order_position)
    canonical_two_alt_slot = len(ordered_alts) == 2 and float(slot.get("mutual_exclusivity", 0)) >= 0.95
    # Order-rule (legacy) signatures: first alt exon = IIIb, second = IIIc.
    order_iiib_sig = ordered_alts[0].sig if canonical_two_alt_slot else ""
    order_iiic_sig = ordered_alts[1].sig if canonical_two_alt_slot else ""
    legacy_assignment_rule = ("conservative_exon_order_rule_applied_only_to_exact_two_alt_exon_mutually_exclusive_slot;"
                              "first=putative_FGFR2_exon8_IIIb;second=putative_FGFR2_exon9_IIIc")
    legacy_direction_validation_status = (
        "exon_structure_direction_provisional_requires_sequence_or_domain_validation"
        if canonical_two_alt_slot else "direction_not_assigned_noncanonical_slot")

    metA: Dict[str, object] = {"prefers": "unresolved", "method": "order_rule_provisional"}
    metB: Dict[str, object] = {"prefers": "unresolved", "method": "order_rule_provisional"}
    direction = {"A_iso": "IIIb", "B_iso": "IIIc", "method": "order_rule_provisional",
                 "status": "direction_not_assigned_noncanonical_slot", "confidence": "low", "warning": ""}
    if canonical_two_alt_slot:
        exonA, exonB = ordered_alts[0], ordered_alts[1]
        if ctx is not None and ctx.enabled:
            txA = [t for t, sset in tx_to_alt.items() if exonA.sig in sset]
            txB = [t for t, sset in tx_to_alt.items() if exonB.sig in sset]
            metA = calibrate_exon(ctx, exonA, txA)
            metB = calibrate_exon(ctx, exonB, txB)
            direction = resolve_pair_direction(metA, metB)
        else:
            direction = {"A_iso": "IIIb", "B_iso": "IIIc", "method": "order_rule_provisional",
                         "status": "direction_unresolved_no_sequence", "confidence": "low",
                         "warning": "sequence_calibration_disabled"}
    calibrated_iso = {order_iiib_sig: direction["A_iso"], order_iiic_sig: direction["B_iso"]}
    metrics_by_sig = {order_iiib_sig: metA, order_iiic_sig: metB}
    _legacy_iso = {order_iiib_sig: "IIIb", order_iiic_sig: "IIIc"}
    direction_is_calibrated = direction["status"] in (
        "sequence_calibrated_pass", "sequence_calibrated_inverted_from_order_rule")

    rows: List[Dict[str, object]] = []
    for tx in tx_rows:
        tx_id = transcript_internal_id(tx)
        sigs = set(tx_to_alt.get(tx_id, set()))
        sigs.update(e.sig for e in exons_by_tx.get(tx_id, []) if e.sig in alt_exons)
        matched_sig = ""
        if not canonical_two_alt_slot and sigs:
            legacy_order_iso = "ambiguous"; conf = "low"
            ev = "noncanonical_alternative_slot;direction_not_assigned_without_sequence_or_domain_validation"
        elif order_iiib_sig in sigs and order_iiic_sig not in sigs:
            legacy_order_iso = "IIIb"; conf = confidence_floor; matched_sig = order_iiib_sig
            ev = "exact_two_alt_exon_slot;contains_first_order_alt_exon;direction_sequence_calibrated"
        elif order_iiic_sig in sigs and order_iiib_sig not in sigs:
            legacy_order_iso = "IIIc"; conf = confidence_floor; matched_sig = order_iiic_sig
            ev = "exact_two_alt_exon_slot;contains_second_order_alt_exon;direction_sequence_calibrated"
        elif order_iiib_sig in sigs and order_iiic_sig in sigs:
            legacy_order_iso = "ambiguous"; conf = "low"
            ev = "contains_both_IIIb_and_IIIc_candidate_exons;possible_annotation_or_transcript_model_issue"
        elif sigs:
            legacy_order_iso = "ambiguous"; conf = "low"
            ev = "contains_noncanonical_alternative_exon_in_detected_FGFR2_III_slot"
        else:
            legacy_order_iso = "unclassified"; conf = "none"
            ev = "does_not_contain_detected_IIIb_or_IIIc_candidate_exon"

        # Sequence-calibrated + final assignment (event role is primary).
        if matched_sig and direction_is_calibrated:
            seq_iso = calibrated_iso.get(matched_sig, "")
            final_iso = seq_iso or legacy_order_iso
        else:
            seq_iso = ""
            final_iso = legacy_order_iso

        met = metrics_by_sig.get(matched_sig, {}) if matched_sig else {}
        rows.append({
            "species_input": get_any(tx, ["species_input", "species_canonical"]),
            "species_canonical": get_any(tx, ["species_canonical", "species_input"]),
            "source_db": get_any(tx, ["source_db", "source"]),
            "gene_symbol": get_any(tx, ["gene_symbol", "symbol"], "FGFR2"),
            "gene_id_internal": get_any(tx, ["gene_id_internal", "gene_id"]),
            "internal_transcript_id": tx_id,
            "transcript_id_source": transcript_source_id(tx),
            "transcript_name": get_any(tx, ["transcript_name", "name"]),
            "isoform_class": final_iso,
            "iii_isoform_assignment": final_iso,
            "confidence": conf,
            "evidence": ev,
            "slot_id": slot.get("slot_id", ""),
            "slot_confidence_reason": confidence_reason,
            "assignment_rule": legacy_assignment_rule,
            "legacy_direction_rule": LEGACY_DIRECTION_RULE,
            "legacy_order_based_isoform_assignment": legacy_order_iso,
            "sequence_calibrated_isoform_assignment": seq_iso,
            "final_iii_isoform_assignment": final_iso,
            "direction_assignment_method": direction.get("method", ""),
            "direction_validation_status": direction.get("status", ""),
            "legacy_direction_validation_status": legacy_direction_validation_status,
            "direction_confidence": direction.get("confidence", ""),
            "direction_warning": direction.get("warning", ""),
            "candidate_vs_human_IIIb_identity": met.get("id_b", ""),
            "candidate_vs_human_IIIc_identity": met.get("id_c", ""),
            "candidate_vs_human_IIIb_coverage": met.get("cov_b", ""),
            "candidate_vs_human_IIIc_coverage": met.get("cov_c", ""),
            "candidate_IIIb_direction_score": met.get("score_b", ""),
            "candidate_IIIc_direction_score": met.get("score_c", ""),
            "candidate_direction_prefers": met.get("prefers", ""),
            "left_flank_sig": slot.get("left_flank_sig", ""),
            "right_flank_sig": slot.get("right_flank_sig", ""),
            "iiib_exon_sig": order_iiib_sig,
            "iiic_exon_sig": order_iiic_sig,
            "matched_alt_exons": ";".join(sorted(sigs)),
            "alt_exon_count_in_slot": len(sigs),
            "protein_length_aa": get_any(tx, ["protein_length_aa", "protein_length"]),
            "score_from_selection": get_any(tx, ["score"]),
            "selection_role": get_any(tx, ["selection_role"]),
        })
    return rows



def make_markdown_report(summary_rows: List[Dict[str, object]], warning_rows: List[Dict[str, object]], metadata: Dict[str, object]) -> str:
    total_groups = len(summary_rows)
    slots = sum(1 for r in summary_rows if str(r.get("slot_detected")) == "1")
    both = sum(1 for r in summary_rows if r.get("status") == "both_isoforms_detected")
    lines = [
        "# FGFR2 IIIb/IIIc exon-structure classification report",
        "",
        f"Script: `{metadata.get('script_name')}` version `{metadata.get('script_version')}`",
        f"Run date: {metadata.get('run_date_utc')}",
        "",
        "## Overview",
        "",
        f"- Species/source/gene groups analysed: {total_groups}",
        f"- Groups with detected mutually exclusive exon slot: {slots}",
        f"- Groups with both IIIb and IIIc detected: {both}",
        f"- Warning rows: {len(warning_rows)}",
        "",
        "## Summary by group",
        "",
        "| Species | Source | Gene | Slot | IIIb | IIIc | Ambiguous | Unclassified | Status |",
        "|---|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in summary_rows:
        lines.append(
            f"| {r.get('species_canonical','')} | {r.get('source_db','')} | {r.get('gene_id_internal','')} | {r.get('slot_detected','')} | "
            f"{r.get('IIIb_transcripts','')} | {r.get('IIIc_transcripts','')} | {r.get('ambiguous_transcripts','')} | {r.get('unclassified_transcripts','')} | {r.get('status','')} |"
        )
    lines.extend(["", "## Interpretation note", "", "IIIb/IIIc calls are structure-derived and conservative. Missing calls can reflect true absence, incomplete transcript models, shifted exon boundaries, or insufficient annotation detail."])
    return "\n".join(lines) + "\n"


def markdown_to_simple_html(markdown_text: str) -> str:
    body = []
    for line in markdown_text.splitlines():
        if line.startswith("# "):
            body.append(f"<h1>{html.escape(line[2:])}</h1>")
        elif line.startswith("## "):
            body.append(f"<h2>{html.escape(line[3:])}</h2>")
        elif line.startswith("- "):
            body.append(f"<p>• {html.escape(line[2:])}</p>")
        elif line.startswith("|"):
            body.append(f"<pre>{html.escape(line)}</pre>")
        elif line.strip() == "":
            body.append("")
        else:
            body.append(f"<p>{html.escape(line)}</p>")
    return "<!doctype html><html><head><meta charset='utf-8'><title>FGFR2 isoform report</title></head><body>" + "\n".join(body) + "</body></html>\n"


def create_summary_plot(summary_rows: List[Dict[str, object]], outdir: Path) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None
    if not summary_rows:
        return None
    by_species: Dict[str, Counter] = defaultdict(Counter)
    for r in summary_rows:
        sp = str(r.get("species_canonical", "unknown"))
        by_species[sp]["IIIb"] += int(r.get("IIIb_transcripts", 0) or 0)
        by_species[sp]["IIIc"] += int(r.get("IIIc_transcripts", 0) or 0)
        by_species[sp]["ambiguous"] += int(r.get("ambiguous_transcripts", 0) or 0)
        by_species[sp]["unclassified"] += int(r.get("unclassified_transcripts", 0) or 0)
    species = sorted(by_species)
    labels = []
    values = []
    for sp in species:
        labels.extend([f"{sp}\nIIIb", f"{sp}\nIIIc", f"{sp}\namb", f"{sp}\nuncl"])
        values.extend([by_species[sp]["IIIb"], by_species[sp]["IIIc"], by_species[sp]["ambiguous"], by_species[sp]["unclassified"]])
    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    path = plot_dir / "isoform_summary_by_species.png"
    fig_width = max(8, min(28, len(labels) * 0.55))
    plt.figure(figsize=(fig_width, 5))
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(labels)), labels, rotation=90)
    plt.ylabel("Transcript count")
    plt.title("FGFR2 IIIb/IIIc exon-structure classification summary")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    return str(path)


def run_pipeline(args: argparse.Namespace) -> Dict[str, object]:
    args.outdir.mkdir(parents=True, exist_ok=True)
    warnings: List[Dict[str, object]] = []
    tx_rows = read_tsv(args.transcripts)
    exon_rows = read_tsv(args.exons)
    selected_rows = read_tsv(args.selected_transcripts) if args.selected_transcripts else []

    warnings.extend(validate_alias_schema(tx_rows, TRANSCRIPT_REQUIRED_ALIASES, "transcripts"))
    warnings.extend(validate_alias_schema(exon_rows, EXON_REQUIRED_ALIASES, "exons"))
    if args.strict and any(w.get("level") == "error" for w in warnings):
        write_tsv(warnings, args.outdir / "fgfr2_isoform_warnings.tsv", WARNING_FIELDS)
        raise SystemExit("Strict mode: missing required input columns or empty input")

    selected_by_internal = {transcript_internal_id(r): r for r in selected_rows if transcript_internal_id(r)}
    merged_txs: List[Dict[str, str]] = []
    for tx in tx_rows:
        tx_id = transcript_internal_id(tx)
        merged = dict(tx)
        sel = selected_by_internal.get(tx_id)
        if sel:
            for k in ["selection_role", "selection_rank", "score", "final_selection_reason"]:
                if k in sel:
                    merged[k] = sel[k]
        merged_txs.append(merged)

    exons_by_tx = load_exons(exon_rows, warnings)

    cds_rows = read_tsv(args.cds_features) if getattr(args, "cds_features", None) else []
    cds_by_tx = load_cds_by_tx(cds_rows)
    tx_by_internal = {transcript_internal_id(t): t for t in merged_txs if transcript_internal_id(t)}
    iiib_ref = read_single_fasta(args.human_iiib_segment_fasta) if getattr(args, "human_iiib_segment_fasta", None) and Path(args.human_iiib_segment_fasta).exists() else ""
    iiic_ref = read_single_fasta(args.human_iiic_segment_fasta) if getattr(args, "human_iiic_segment_fasta", None) and Path(args.human_iiic_segment_fasta).exists() else ""
    wanted_acc = {clean_id(get_any(t, ["translation_id_source"])) for t in merged_txs if get_any(t, ["translation_id_source"])}
    wanted_tx = {clean_id(get_any(t, ["transcript_id_source", "transcript_id"])) for t in merged_txs if get_any(t, ["transcript_id_source", "transcript_id"])}
    cache_index = PL.ProteinCacheIndex(args.cache, wanted_accessions=wanted_acc, wanted_transcripts=wanted_tx) if getattr(args, "cache", None) else PL.ProteinCacheIndex(None, wanted_accessions=wanted_acc, wanted_transcripts=wanted_tx)
    calibration_enabled = bool(iiib_ref and iiic_ref and cds_by_tx)
    if not calibration_enabled:
        warnings.append({"level": "warning", "warning_type": "sequence_calibration_disabled",
                         "message": "Missing cds_features or human IIIb/IIIc references; direction falls back to provisional order rule."})
    ctx = DirectionContext(
        cache_index=cache_index, protein_cache={}, cds_by_tx=cds_by_tx,
        tx_by_internal=tx_by_internal, iiib_ref=iiib_ref, iiic_ref=iiic_ref,
        args=args, enabled=calibration_enabled,
    )

    duplicate_internal = [tx for tx, c in Counter(transcript_internal_id(r) for r in merged_txs if transcript_internal_id(r)).items() if c > 1]
    for tx in duplicate_internal:
        warnings.append({"level": "warning", "warning_type": "duplicate_internal_transcript_id", "transcript_id": tx, "message": "Duplicate transcript identifier in transcript table"})

    groups: Dict[Tuple[str, str, str, str], List[Dict[str, str]]] = defaultdict(list)
    excluded_by_group: Counter = Counter()
    for tx in merged_txs:
        if args.gene_symbol and get_any(tx, ["gene_symbol", "symbol"], args.gene_symbol).upper() != args.gene_symbol.upper():
            text = " ".join(str(v) for v in tx.values()).upper()
            if args.gene_symbol.upper() not in text:
                continue
        tx_id = transcript_internal_id(tx)
        gkey = group_key_from_tx(tx)
        if not tx_id or tx_id not in exons_by_tx:
            warnings.append({"level": "warning", "warning_type": "transcript_without_exons", "species_canonical": gkey[0], "source_db": gkey[1], "gene_id_internal": gkey[2], "transcript_id": tx_id or transcript_source_id(tx), "message": "Transcript lacks linked exon rows and was not used for slot detection"})
            excluded_by_group[gkey] += 1
            continue
        ok, reasons = is_complete_candidate_tx(tx, exons_by_tx[tx_id], args.min_protein_aa, args.max_protein_aa)
        if not ok:
            for reason in reasons:
                warnings.append({"level": "warning", "warning_type": "transcript_excluded", "exclusion_reason": reason, "species_canonical": gkey[0], "source_db": gkey[1], "gene_id_internal": gkey[2], "transcript_id": tx_id, "message": "Transcript was excluded from IIIb/IIIc slot detection"})
            excluded_by_group[gkey] += 1
            continue
        groups[gkey].append(tx)

    all_evidence: List[Dict[str, object]] = []
    audit_rows: List[Dict[str, object]] = []
    alt_exon_metadata_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []

    for gkey, rows in sorted(groups.items()):
        species, source, gene, symbol = gkey
        candidates = build_slot_candidates(rows, exons_by_tx, args.min_alt_exon_len, args.max_alt_exon_len, args.min_tx_fraction_with_flanks)
        slot = choose_fgfr2_slot(candidates)
        if len(candidates) > 1:
            warnings.append({"level": "warning", "warning_type": "multiple_candidate_slots", "species_canonical": species, "source_db": source, "gene_id_internal": gene, "message": f"{len(candidates)} candidate slots detected; top-ranked slot selected"})
        for rank, cand in enumerate(candidates, start=1):
            alt_exons: Dict[str, Exon] = cand["alt_exons"]  # type: ignore[assignment]
            ordered = sorted(alt_exons.values(), key=transcript_order_position)
            if int(cand.get("alt_count", 0)) > 2:
                warnings.append({"level": "warning", "warning_type": "slot_with_more_than_two_alternative_exons", "species_canonical": species, "source_db": source, "gene_id_internal": gene, "slot_id": cand.get("slot_id", ""), "message": "Candidate slot contains more than two alternative exons; noncanonical alternatives are classified as ambiguous"})
            audit_rows.append({
                "species_canonical": species, "source_db": source, "gene_id_internal": gene, "gene_symbol": symbol,
                "candidate_rank": rank, "selected_as_fgfr2_iii_slot": "1" if slot is cand else "0",
                "slot_id": cand.get("slot_id", ""), "score": cand.get("score", ""), "alt_count": cand.get("alt_count", ""),
                "tx_with_slot": cand.get("tx_with_slot", ""), "tx_total": cand.get("tx_total", ""), "tx_fraction": cand.get("tx_fraction", ""),
                "mutual_exclusivity": cand.get("mutual_exclusivity", ""), "length_spread": cand.get("length_spread", ""),
                "alt_pair_nonoverlapping": cand.get("alt_pair_nonoverlapping", ""), "alt_pair_min_gap_bp": cand.get("alt_pair_min_gap_bp", ""), "alt_pair_order_distance_bp": cand.get("alt_pair_order_distance_bp", ""), "alt_pair_geometry_class": cand.get("alt_pair_geometry_class", ""), "fgfr2_iii_geometry_score": cand.get("fgfr2_iii_geometry_score", ""), "left_flank_sig": cand.get("left_flank_sig", ""), "right_flank_sig": cand.get("right_flank_sig", ""),
                "alternative_exon_sigs_in_transcript_order": ";".join(e.sig for e in ordered),
                "alternative_exon_lengths_in_transcript_order": ";".join(str(e.length) for e in ordered),
                "slot_confidence_reason": slot_confidence_and_reason(cand)[1],
            })

            tx_to_alt: Dict[str, set] = cand["tx_to_alt"]
            for alt_order, ex in enumerate(ordered, start=1):
                if alt_order == 1:
                    inferred_isoform = "IIIb"
                    biological_label = "FGFR2_exon8_IIIb_candidate"
                elif alt_order == 2:
                    inferred_isoform = "IIIc"
                    biological_label = "FGFR2_exon9_IIIc_candidate"
                else:
                    inferred_isoform = "noncanonical_alternative"
                    biological_label = f"additional_alternative_exon_{alt_order}"
                tx_with_this_alt = sorted(tx for tx, sigs in tx_to_alt.items() if ex.sig in sigs)
                alt_exon_metadata_rows.append({
                    "species_canonical": species,
                    "source_db": source,
                    "gene_id_internal": gene,
                    "gene_symbol": symbol,
                    "candidate_rank": rank,
                    "selected_as_fgfr2_iii_slot": "1" if slot is cand else "0",
                    "slot_id": cand.get("slot_id", ""),
                    "slot_confidence_reason": slot_confidence_and_reason(cand)[1],
                    "alternative_exon_order_in_transcript": alt_order,
                    "inferred_isoform": inferred_isoform,
                    "biological_label": biological_label,
                    "exon_sig": ex.sig,
                    "exon_id_source": ex.exon_id_source,
                    "chrom": ex.chrom,
                    "start": ex.start,
                    "end": ex.end,
                    "strand": ex.strand,
                    "length_bp": ex.length,
                    "representative_exon_rank": ex.rank,
                    "alt_pair_nonoverlapping": cand.get("alt_pair_nonoverlapping", ""),
                    "alt_pair_min_gap_bp": cand.get("alt_pair_min_gap_bp", ""),
                    "alt_pair_order_distance_bp": cand.get("alt_pair_order_distance_bp", ""),
                    "alt_pair_geometry_class": cand.get("alt_pair_geometry_class", ""),
                    "fgfr2_iii_geometry_score": cand.get("fgfr2_iii_geometry_score", ""),
                    "left_flank_sig": cand.get("left_flank_sig", ""),
                    "right_flank_sig": cand.get("right_flank_sig", ""),
                    "tx_count_with_this_alt": len(tx_with_this_alt),
                    "transcripts_with_this_alt": ";".join(tx_with_this_alt),
                    "mutual_exclusivity": cand.get("mutual_exclusivity", ""),
                    "tx_with_slot": cand.get("tx_with_slot", ""),
                    "tx_total": cand.get("tx_total", ""),
                })
        if slot is None:
            warnings.append({"level": "warning", "warning_type": "no_slot_detected", "species_canonical": species, "source_db": source, "gene_id_internal": gene, "message": "No reliable mutually exclusive FGFR2 IIIb/IIIc exon slot detected"})
            for tx in rows:
                all_evidence.append({
                    "species_input": get_any(tx, ["species_input", "species_canonical"]), "species_canonical": species, "source_db": source,
                    "gene_symbol": symbol, "gene_id_internal": gene, "internal_transcript_id": transcript_internal_id(tx),
                    "transcript_id_source": transcript_source_id(tx), "transcript_name": get_any(tx, ["transcript_name", "name"]),
                    "isoform_class": "unclassified", "iii_isoform_assignment": "unclassified", "confidence": "none",
                    "evidence": "no_mutually_exclusive_FGFR2_IIIb_IIIc_exon_slot_detected_in_this_species_source_model",
                    "slot_id": "", "slot_confidence_reason": "no_slot_detected", "assignment_rule": "not_applied",
                    "left_flank_sig": "", "right_flank_sig": "", "iiib_exon_sig": "", "iiic_exon_sig": "",
                    "matched_alt_exons": "", "alt_exon_count_in_slot": 0,
                    "protein_length_aa": get_any(tx, ["protein_length_aa", "protein_length"]), "score_from_selection": get_any(tx, ["score"]),
                    "selection_role": get_any(tx, ["selection_role"]),
                })
            summary_rows.append({
                "species_canonical": species, "source_db": source, "gene_id_internal": gene, "gene_symbol": symbol,
                "transcripts_considered": len(rows), "transcripts_excluded": excluded_by_group[gkey], "slot_detected": 0,
                "IIIb_transcripts": 0, "IIIc_transcripts": 0, "ambiguous_transcripts": 0, "unclassified_transcripts": len(rows),
                "status": "no_slot_detected", "review_note": "No reliable mutually exclusive exon slot was detected; manual review or sequence-based evidence is recommended.",
            })
            continue
        conf, conf_reason = slot_confidence_and_reason(slot)
        evidence_rows = classify_with_slot(rows, exons_by_tx, slot, confidence_floor=conf, confidence_reason=conf_reason, ctx=ctx)
        all_evidence.extend(evidence_rows)
        dstatus = evidence_rows[0].get("direction_validation_status", "") if evidence_rows else ""
        if dstatus == "sequence_calibrated_inverted_from_order_rule":
            warnings.append({"level": "warning", "warning_type": "sequence_calibrated_direction_inverted_from_order_rule",
                             "species_canonical": species, "source_db": source, "gene_id_internal": gene,
                             "message": "Sequence calibration inverted the IIIb/IIIc direction relative to the provisional exon-order rule."})
        elif dstatus == "ambiguous_sequence_direction_review":
            warnings.append({"level": "warning", "warning_type": "ambiguous_sequence_direction_review",
                             "species_canonical": species, "source_db": source, "gene_id_internal": gene,
                             "message": "Sequence calibration could not resolve IIIb/IIIc direction; provisional order rule retained for review."})
        elif dstatus == "direction_unresolved_no_sequence":
            warnings.append({"level": "warning", "warning_type": "direction_unresolved_no_sequence",
                             "species_canonical": species, "source_db": source, "gene_id_internal": gene,
                             "message": "No usable cassette sequence for direction calibration; provisional order rule retained."})
        counts = Counter(r["isoform_class"] for r in evidence_rows)
        status = "both_isoforms_detected" if counts.get("IIIb", 0) and counts.get("IIIc", 0) else "partial_or_single_isoform_detected"
        if status != "both_isoforms_detected":
            warnings.append({"level": "warning", "warning_type": "only_one_or_no_isoform_class_detected", "species_canonical": species, "source_db": source, "gene_id_internal": gene, "message": "Detected slot did not yield both IIIb and IIIc transcript classes"})
        summary_rows.append({
            "species_canonical": species, "source_db": source, "gene_id_internal": gene, "gene_symbol": symbol,
            "transcripts_considered": len(rows), "transcripts_excluded": excluded_by_group[gkey], "slot_detected": 1,
            "IIIb_transcripts": counts.get("IIIb", 0), "IIIc_transcripts": counts.get("IIIc", 0),
            "ambiguous_transcripts": counts.get("ambiguous", 0), "unclassified_transcripts": counts.get("unclassified", 0),
            "status": status,
            "review_note": "" if status == "both_isoforms_detected" else "Only one isoform class was detected; inspect annotation completeness.",
        })

    metadata: Dict[str, object] = {
        "script_name": SCRIPT_NAME,
        "script_version": SCRIPT_VERSION,
        "run_date_utc": datetime.now(timezone.utc).isoformat(),
        "input_files": {"transcripts": str(args.transcripts), "exons": str(args.exons), "selected_transcripts": str(args.selected_transcripts) if args.selected_transcripts else None},
        "parameters": {
            "gene_symbol": args.gene_symbol, "min_alt_exon_len": args.min_alt_exon_len, "max_alt_exon_len": args.max_alt_exon_len,
            "min_protein_aa": args.min_protein_aa, "max_protein_aa": args.max_protein_aa,
            "min_tx_fraction_with_flanks": args.min_tx_fraction_with_flanks, "strict": args.strict,
        },
        "input_row_counts": {"transcripts": len(tx_rows), "exons": len(exon_rows), "selected_transcripts": len(selected_rows)},
        "output_row_counts": {"evidence": len(all_evidence), "audit": len(audit_rows), "alternative_exon_metadata": len(alt_exon_metadata_rows), "summary": len(summary_rows), "warnings": len(warnings)},
        "number_groups": len(groups),
        "number_slots_detected": sum(1 for r in summary_rows if str(r.get("slot_detected")) == "1"),
    }
    # Direction-calibration statistics (per species/source/gene slot, deduplicated).
    slot_dir: Dict[Tuple[str, str, str], str] = {}
    slot_method: Dict[Tuple[str, str, str], str] = {}
    for r in all_evidence:
        if not r.get("slot_id"):
            continue
        k = (str(r.get("species_canonical", "")), str(r.get("source_db", "")), str(r.get("slot_id", "")))
        st = str(r.get("direction_validation_status", ""))
        if st:
            slot_dir[k] = st
            slot_method[k] = str(r.get("direction_assignment_method", ""))
    metadata["direction_calibration"] = {
        "calibration_enabled": calibration_enabled,
        "slots_with_direction_status": len(slot_dir),
        "direction_validation_status_counts": dict(Counter(slot_dir.values())),
        "direction_method_counts": dict(Counter(slot_method.values())),
        "species_inverted_from_order_rule": sorted({k[0] for k, v in slot_dir.items() if v == "sequence_calibrated_inverted_from_order_rule"}),
        "species_ambiguous_direction": sorted({k[0] for k, v in slot_dir.items() if v == "ambiguous_sequence_direction_review"}),
        "species_unresolved_direction": sorted({k[0] for k, v in slot_dir.items() if v == "direction_unresolved_no_sequence"}),
    }

    write_tsv(all_evidence, args.outdir / "fgfr2_isoform_evidence.tsv", EVIDENCE_FIELDS)
    write_tsv(audit_rows, args.outdir / "fgfr2_isoform_detection_audit.tsv", AUDIT_FIELDS)
    write_tsv(alt_exon_metadata_rows, args.outdir / "fgfr2_alternative_exon_metadata.tsv", ALT_EXON_METADATA_FIELDS)
    write_tsv(summary_rows, args.outdir / "fgfr2_isoform_summary.tsv", SUMMARY_FIELDS)
    write_tsv(warnings, args.outdir / "fgfr2_isoform_warnings.tsv", WARNING_FIELDS)
    (args.outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    report_md = make_markdown_report(summary_rows, warnings, metadata)
    (args.outdir / "fgfr2_isoform_report.md").write_text(report_md, encoding="utf-8")
    (args.outdir / "fgfr2_isoform_report.html").write_text(markdown_to_simple_html(report_md), encoding="utf-8")
    plot_path = create_summary_plot(summary_rows, args.outdir)
    if plot_path:
        metadata["plot_path"] = plot_path
        (args.outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")

    methods = f"""FGFR2 IIIb/IIIc isoform classification method

FGFR2 IIIb and IIIc isoforms were classified using an exon-structure-based, annotation-aware strategy. The method is based on the established biology that FGFR2 IIIb/IIIc isoforms are produced by mutually exclusive inclusion of two internal exons in the IgIII/D3 ligand-binding region. For every species/source/gene model, transcript exons were sorted by annotated exon rank. Internal exons with lengths between {args.min_alt_exon_len} and {args.max_alt_exon_len} bp were evaluated as candidate alternative IgIII exons. Candidate mutually exclusive slots were defined by identical immediate upstream and downstream flanking exons but two or more alternative internal exons. The best slot was selected conservatively by first prioritising exact two-exon, mutually exclusive, non-overlapping cassette geometry consistent with FGFR2 IIIb/IIIc, and only then considering transcript support, number of alternatives, mutual exclusivity, and alternative-exon length similarity. Overlapping splice-site-like alternatives are retained in the audit table but are not allowed to outrank separated cassette exons.

Isoform direction is assigned conservatively. The transcript-order rule is applied only when the detected slot contains exactly two mutually exclusive alternative exons with high mutual exclusivity. In that restricted case, the first alternative exon in transcript order is labelled as putative IIIb and the second as putative IIIc, reflecting the canonical FGFR2 exon 8/exon 9 arrangement. Slots with more than two alternatives, weak mutual exclusivity, or transcripts containing both alternatives are not forced into IIIb/IIIc classes and are marked ambiguous. The output explicitly records that exon-order direction remains provisional and should be checked against downstream sequence/domain evidence.

The procedure is intentionally conservative: if no mutually exclusive FGFR2 IIIb/IIIc exon slot is detected for a species/source/gene model, transcripts remain unclassified rather than being forced into an unsupported isoform class. Excluded transcripts, multiple candidate slots, noncanonical extra alternatives and missing schema fields are written to fgfr2_isoform_warnings.tsv. The main output fgfr2_isoform_evidence.tsv can be supplied to select_fgfr2_transcripts_annotation_aware_v2.py with --isoform_evidence. Detailed candidate slots, summary counts, run metadata, Markdown/HTML reports and a summary plot are written for manual review and reproducibility.
"""
    (args.outdir / "fgfr2_isoform_methods.txt").write_text(methods, encoding="utf-8")
    return metadata


EVIDENCE_FIELDS = [
    "species_input", "species_canonical", "source_db", "gene_symbol", "gene_id_internal",
    "internal_transcript_id", "transcript_id_source", "transcript_name", "isoform_class", "iii_isoform_assignment",
    "confidence", "evidence", "slot_id", "slot_confidence_reason", "assignment_rule",
    # Sequence-calibrated IIIb/IIIc direction (Step 4 hotfix). Legacy order-based
    # labels are preserved alongside the calibrated/final assignment.
    "legacy_direction_rule", "legacy_order_based_isoform_assignment",
    "sequence_calibrated_isoform_assignment", "final_iii_isoform_assignment",
    "direction_assignment_method", "direction_validation_status", "legacy_direction_validation_status",
    "direction_confidence", "direction_warning",
    "candidate_vs_human_IIIb_identity", "candidate_vs_human_IIIc_identity",
    "candidate_vs_human_IIIb_coverage", "candidate_vs_human_IIIc_coverage",
    "candidate_IIIb_direction_score", "candidate_IIIc_direction_score", "candidate_direction_prefers",
    "left_flank_sig", "right_flank_sig",
    "iiib_exon_sig", "iiic_exon_sig", "matched_alt_exons", "alt_exon_count_in_slot", "protein_length_aa", "score_from_selection", "selection_role",
]
AUDIT_FIELDS = [
    "species_canonical", "source_db", "gene_id_internal", "gene_symbol", "candidate_rank", "selected_as_fgfr2_iii_slot",
    "slot_id", "score", "alt_count", "tx_with_slot", "tx_total", "tx_fraction", "mutual_exclusivity", "length_spread",
    "alt_pair_nonoverlapping", "alt_pair_min_gap_bp", "alt_pair_order_distance_bp", "alt_pair_geometry_class", "fgfr2_iii_geometry_score",
    "left_flank_sig", "right_flank_sig", "alternative_exon_sigs_in_transcript_order", "alternative_exon_lengths_in_transcript_order", "slot_confidence_reason",
]

ALT_EXON_METADATA_FIELDS = [
    "species_canonical", "source_db", "gene_id_internal", "gene_symbol",
    "candidate_rank", "selected_as_fgfr2_iii_slot", "slot_id", "slot_confidence_reason",
    "alternative_exon_order_in_transcript", "inferred_isoform", "biological_label",
    "exon_sig", "exon_id_source", "chrom", "start", "end", "strand", "length_bp",
    "representative_exon_rank", "alt_pair_nonoverlapping", "alt_pair_min_gap_bp", "alt_pair_order_distance_bp", "alt_pair_geometry_class", "fgfr2_iii_geometry_score",
    "left_flank_sig", "right_flank_sig",
    "tx_count_with_this_alt", "transcripts_with_this_alt",
    "mutual_exclusivity", "tx_with_slot", "tx_total",
]

SUMMARY_FIELDS = [
    "species_canonical", "source_db", "gene_id_internal", "gene_symbol", "transcripts_considered", "transcripts_excluded",
    "slot_detected", "IIIb_transcripts", "IIIc_transcripts", "ambiguous_transcripts", "unclassified_transcripts", "status", "review_note",
]
WARNING_FIELDS = [
    "level", "warning_type", "file", "logical_column", "accepted_aliases", "row_number", "species_canonical", "source_db", "gene_id_internal",
    "transcript_id", "slot_id", "exclusion_reason", "message",
]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Classify FGFR2 IIIb/IIIc isoforms from exon structure.")
    parser.add_argument("--transcripts", type=Path, required=True)
    parser.add_argument("--exons", type=Path, required=True)
    parser.add_argument("--selected_transcripts", type=Path, default=None, help="Optional selected_transcripts.tsv for selection roles/scores.")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--gene_symbol", default="FGFR2")
    parser.add_argument("--min_alt_exon_len", type=int, default=120)
    parser.add_argument("--max_alt_exon_len", type=int, default=230)
    parser.add_argument("--min_protein_aa", type=int, default=550)
    parser.add_argument("--max_protein_aa", type=int, default=1000)
    parser.add_argument("--min_tx_fraction_with_flanks", type=float, default=0.10)
    parser.add_argument("--strict", action="store_true", help="Abort on missing required schema or empty inputs.")
    # Sequence-calibrated IIIb/IIIc direction (Step 4 hotfix).
    parser.add_argument("--cds_features", type=Path, default=None, help="cds_features.tsv for exon->protein AA mapping.")
    parser.add_argument("--cache", type=Path, default=None, help="NCBI datasets cache for candidate protein retrieval.")
    parser.add_argument("--human_iiib_segment_fasta", type=Path, default=None, help="Curated human IIIb cassette reference (FASTA).")
    parser.add_argument("--human_iiic_segment_fasta", type=Path, default=None, help="Curated human IIIc cassette reference (FASTA).")
    parser.add_argument("--no_ensembl_rest", action="store_true", help="Disable Ensembl REST fallback for protein retrieval.")
    parser.add_argument("--ensembl_sleep", type=float, default=0.4)
    parser.add_argument("--ensembl_timeout", type=int, default=30)
    parser.add_argument("--max_proteins_per_exon", type=int, default=6, help="Max candidate proteins evaluated per alternative exon.")
    parser.add_argument("--direction_margin", type=float, default=0.08, help="Min direction-score margin to call IIIb vs IIIc.")
    parser.add_argument("--direction_min_identity", type=float, default=0.45, help="Min direction score (identity*coverage) required to resolve direction.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    metadata = run_pipeline(args)
    print(f"Wrote {metadata['output_row_counts']['evidence']} transcript isoform evidence rows")
    print(f"Wrote {metadata['output_row_counts']['audit']} candidate slot audit rows")
    print(f"Wrote {metadata['output_row_counts'].get('alternative_exon_metadata', 0)} alternative exon metadata rows")
    print(f"Wrote {metadata['output_row_counts']['summary']} species/source/gene summary rows")
    print(f"Wrote {metadata['output_row_counts']['warnings']} warning rows")
    dc = metadata.get("direction_calibration", {})
    print(f"Direction calibration enabled: {dc.get('calibration_enabled')}")
    print(f"Direction status counts: {dc.get('direction_validation_status_counts')}")
    print(f"Species inverted from order rule: {dc.get('species_inverted_from_order_rule')}")
    print(f"Output directory: {args.outdir}")


if __name__ == "__main__":
    main()
