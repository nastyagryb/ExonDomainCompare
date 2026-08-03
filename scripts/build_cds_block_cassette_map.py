#!/usr/bin/env python3
"""
Map cassette intervals to CDS blocks.

Fixes the join bug where IIIb/IIIc cassettes were mapped onto CDS blocks by a
NON-UNIQUE NCBI/RefSeq CDS id (cds-XP_...), which collapsed many cassettes onto
the first CDS block (protein_start_aa = 1). Cassettes are now mapped by GENOMIC
(and protein) coordinate overlap against a table of UNIQUE CDS blocks.

Inputs (final tables only):
  --coordinate_audit  fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv (60 rows)
  --cds_features      02_models/cds_features.tsv (all per-transcript CDS blocks)
  --proteins          selected_fgfr2_proteins.faa (protein lengths / sequences)
  [--cds_fasta_dir]   optional dir of reconstructed CDS nucleotide FASTA (Part D patch)

Outputs (into --outdir):
  PART A  fgfr2_unique_cds_block_table.tsv
  PART B  fgfr2_cassette_cds_block_map.tsv
  PART C  fgfr2_transcript_cds_reconstruction_audit.tsv
  PART E  fgfr2_cassette_coordinate_sanity_audit.tsv
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# A full-length FGFR2 protein is ~800 aa; the IIIb/IIIc cassette lies in the IgIII/D3
# region (roughly aa 300-400). These thresholds are conservative coordinate sanity
# only (pre-InterPro; no domain calls).
FULL_LENGTH_AA = 500
MIN_PLAUSIBLE_CASSETTE_START = 150


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, rows: List[Dict[str, object]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def _int(v, default=None):
    try:
        return int(float(str(v).strip()))
    except Exception:
        return default


def norm_tx(tx: str) -> str:
    t = (tx or "").strip()
    for pref in ("rna-", "transcript:", "transcript-"):
        if t.lower().startswith(pref):
            t = t[len(pref):]
    return t


def overlap_inclusive(a0: int, a1: int, b0: int, b1: int) -> int:
    """Inclusive genomic overlap length (bp)."""
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def parse_proteins(fasta: Optional[Path]) -> Tuple[Dict[str, int], Dict[str, str]]:
    lengths: Dict[str, int] = {}
    seqs: Dict[str, str] = {}
    if not fasta or not Path(fasta).exists():
        return lengths, seqs
    pid, buf = None, []
    with open(fasta, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if pid:
                    s = "".join(buf)
                    if len(s) >= len(seqs.get(pid, "")):
                        seqs[pid] = s
                        lengths[pid] = len(s)
                pid, buf = None, []
                for tok in line[1:].strip().split("|"):
                    if tok.startswith("protein="):
                        pid = tok.split("=", 1)[1].strip()
            else:
                buf.append(line.strip())
        if pid:
            s = "".join(buf)
            if len(s) >= len(seqs.get(pid, "")):
                seqs[pid] = s
                lengths[pid] = len(s)
    return lengths, seqs


# ---------------------------------------------------------------------------
# PART A — unique CDS block table
# ---------------------------------------------------------------------------
UNIQUE_COLS = [
    "species", "transcript_id", "protein_id", "seqid", "strand",
    "genomic_start", "genomic_end", "cds_phase", "cds_rank",
    "cumulative_CDS_nt_start", "cumulative_CDS_nt_end",
    "protein_start_aa", "protein_end_aa", "source_annotation", "unique_cds_block_id",
]


def build_unique_blocks(cds_by_tx: Dict[str, List[Dict[str, str]]],
                        tx_meta: Dict[str, Dict[str, str]]) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for tx, blocks in cds_by_tx.items():
        # cds_rank already encodes biological transcript 5'->3' order (protein_start_aa
        # increases with cds_rank for both strands); sort defensively by cds_rank.
        blocks = sorted(blocks, key=lambda b: _int(b.get("cds_rank"), 0) or 0)
        meta = tx_meta.get(tx, {})
        cum = 0
        for b in blocks:
            ln = _int(b.get("cds_length_bp"), 0) or 0
            start, end = _int(b.get("start")), _int(b.get("end"))
            seqid = b.get("chrom", "")
            strand = b.get("strand", "")
            rank = b.get("cds_rank", "")
            sp = b.get("species_canonical", "")
            uid = f"{sp}|{tx}|{rank}|{seqid}:{start}-{end}:{strand}"
            out.append({
                "species": sp, "transcript_id": tx,
                "protein_id": b.get("translation_id_source", "") or meta.get("protein_id", ""),
                "seqid": seqid, "strand": strand,
                "genomic_start": start, "genomic_end": end,
                "cds_phase": b.get("phase", ""), "cds_rank": rank,
                "cumulative_CDS_nt_start": cum + 1, "cumulative_CDS_nt_end": cum + ln,
                "protein_start_aa": b.get("protein_start_aa", ""),
                "protein_end_aa": b.get("protein_end_aa", ""),
                "source_annotation": b.get("source_db", ""),
                "unique_cds_block_id": uid,
            })
            cum += ln
    return out


# ---------------------------------------------------------------------------
# PART B — cassette mapping by coordinate overlap (never by non-unique cds_id)
# ---------------------------------------------------------------------------
MAP_COLS = [
    "species", "isoform", "transcript_id", "protein_id", "resolver_source_db",
    "resolver_cds_rank", "resolver_genomic_start", "resolver_genomic_end",
    "matched_cds_rank", "matched_unique_cds_block_id",
    "matched_protein_start_aa", "matched_protein_end_aa",
    "overlap_nt_with_resolved_cassette", "overlap_aa_with_resolved_cassette",
    "overlap_fraction_of_cassette", "overlap_fraction_of_cds_block",
    "n_blocks_overlapping", "cassette_overlap_status", "mapping_method",
]


def map_cassette(coord_row: Dict[str, str], blocks: List[Dict[str, str]]):
    """Return a mapping dict for one resolved IIIb/IIIc cassette using genomic overlap
    (primary) then protein-coordinate overlap (fallback)."""
    rs, re = _int(coord_row.get("resolver_start")), _int(coord_row.get("resolver_end"))
    ns, ne = (_int(coord_row.get("native_protein_start_aa")),
              _int(coord_row.get("native_protein_end_aa")))
    cas_len_nt = (re - rs + 1) if (rs is not None and re is not None and re >= rs) else None

    best = None  # (overlap_nt, block, idx)
    n_overlap = 0
    if rs is not None and re is not None:
        for idx, b in enumerate(blocks):
            bs, be = _int(b.get("start")), _int(b.get("end"))
            if bs is None or be is None:
                continue
            ov = overlap_inclusive(rs, re, bs, be)
            if ov > 0:
                n_overlap += 1
                if best is None or ov > best[0]:
                    best = (ov, b, idx)

    method, status = "genomic_overlap", "cassette_overlap_unresolved"
    matched = {}
    if best:
        ov, b, idx = best
        bs, be = _int(b.get("start")), _int(b.get("end"))
        block_len = be - bs + 1
        frac_cas = ov / cas_len_nt if cas_len_nt else 0.0
        frac_blk = ov / block_len if block_len else 0.0
        if frac_cas >= 0.9:
            status = "cassette_block_exact_overlap"
        elif frac_cas >= 0.25 or frac_blk >= 0.25:
            status = "cassette_block_partial_overlap"
        else:
            status = "cassette_block_partial_overlap"
        if n_overlap > 1 and frac_cas < 0.9:
            status = "cassette_mapping_conflict_review"
        matched = {
            "matched_cds_rank": b.get("cds_rank", ""),
            "matched_unique_cds_block_id":
                f"{b.get('species_canonical','')}|{norm_tx(b.get('transcript_id_source'))}|"
                f"{b.get('cds_rank','')}|{b.get('chrom','')}:{bs}-{be}:{b.get('strand','')}",
            "matched_protein_start_aa": b.get("protein_start_aa", ""),
            "matched_protein_end_aa": b.get("protein_end_aa", ""),
            "overlap_nt_with_resolved_cassette": ov,
            "overlap_aa_with_resolved_cassette": ov // 3,
            "overlap_fraction_of_cassette": round(frac_cas, 3),
            "overlap_fraction_of_cds_block": round(frac_blk, 3),
        }
    else:
        # protein-coordinate fallback: best block by protein-interval overlap
        if ns is not None and ne is not None and blocks:
            bbest, bov = None, 0
            for b in blocks:
                ps, pe = _int(b.get("protein_start_aa")), _int(b.get("protein_end_aa"))
                if ps is None or pe is None:
                    continue
                ov = max(0, min(pe, ne) - max(ps, ns) + 1)
                if ov > bov:
                    bbest, bov = b, ov
            if bbest and bov > 0:
                method, status = "protein_overlap", "cassette_block_partial_overlap"
                bs, be = _int(bbest.get("start")), _int(bbest.get("end"))
                matched = {
                    "matched_cds_rank": bbest.get("cds_rank", ""),
                    "matched_unique_cds_block_id":
                        f"{bbest.get('species_canonical','')}|{norm_tx(bbest.get('transcript_id_source'))}|"
                        f"{bbest.get('cds_rank','')}|{bbest.get('chrom','')}:{bs}-{be}:{bbest.get('strand','')}",
                    "matched_protein_start_aa": bbest.get("protein_start_aa", ""),
                    "matched_protein_end_aa": bbest.get("protein_end_aa", ""),
                    "overlap_aa_with_resolved_cassette": bov,
                    "overlap_fraction_of_cassette": round(bov / (ne - ns + 1), 3) if ne > ns else "",
                }
        if not matched:
            # no CDS blocks at all -> overlay from resolved protein interval
            method, status = "resolved_protein_interval", "cassette_overlay_from_resolved_protein_interval"

    # An FGFR2 IIIb/IIIc cassette at aa 1 is a mapping artefact; prefer the validated native interval.
    if matched:
        eff_start = _int(matched.get("matched_protein_start_aa"))
        eff_end = _int(matched.get("matched_protein_end_aa"))
    else:
        eff_start, eff_end = ns, ne
    plen = _int(coord_row.get("protein_length_aa")) or _int(coord_row.get("native_protein_length_aa"))
    full_len = plen is not None and plen > FULL_LENGTH_AA
    native_plausible = ns is not None and ns >= MIN_PLAUSIBLE_CASSETTE_START
    # Reject first-block and implausibly early mappings in full-length proteins.
    block_implausible = eff_start is not None and (
        eff_start <= 1 or (full_len and eff_start < MIN_PLAUSIBLE_CASSETTE_START))
    aa1_override = False
    if block_implausible:
        if native_plausible:
            # Use the upstream-validated native cassette interval.
            eff_start, eff_end = ns, ne
            method = f"{method}+sequence_calibrated_native_override"
            status = f"{status}_native_protein_interval_override"
            aa1_override = True
        elif eff_start is not None and eff_start <= 1:
            # no plausible native alternative -> never emit the biologically-invalid aa 1
            eff_start, eff_end = None, None
            status = "cassette_unresolved_no_safe_protein_interval_aa1_suppressed"
            aa1_override = True
    if matched:
        matched["matched_protein_start_aa"] = eff_start if eff_start is not None else ""
        matched["matched_protein_end_aa"] = eff_end if eff_end is not None else ""
        if aa1_override:
            # the cds_rank-1 block match was a coordinate artifact (cassette is mid-protein);
            # clear the artifact rank so downstream sanity does not re-flag first_cds_block.
            matched["matched_cds_rank"] = ""

    row = {
        "species": coord_row.get("species_canonical", ""),
        "isoform": coord_row.get("inferred_isoform", ""),
        "transcript_id": norm_tx(coord_row.get("transcript_id_source")),
        "protein_id": coord_row.get("protein_id", ""),
        "resolver_source_db": coord_row.get("resolver_source_db", ""),
        "resolver_cds_rank": coord_row.get("resolver_cds_rank", ""),
        "resolver_genomic_start": rs if rs is not None else "",
        "resolver_genomic_end": re if re is not None else "",
        "n_blocks_overlapping": n_overlap,
        "cassette_overlap_status": status, "mapping_method": method,
        "matched_cds_rank": "", "matched_unique_cds_block_id": "",
        "matched_protein_start_aa": eff_start if eff_start is not None else "",
        "matched_protein_end_aa": eff_end if eff_end is not None else "",
        "overlap_nt_with_resolved_cassette": "", "overlap_aa_with_resolved_cassette": "",
        "overlap_fraction_of_cassette": "", "overlap_fraction_of_cds_block": "",
    }
    row.update(matched)
    return row


# ---------------------------------------------------------------------------
# PART C — transcript CDS reconstruction + translation audit
# ---------------------------------------------------------------------------
RECON_COLS = [
    "species", "isoform", "transcript_id", "protein_id", "source_annotation", "strand",
    "n_cds_blocks", "reconstructed_cds_nt_length", "reconstructed_protein_length",
    "selected_protein_length", "translation_matches_selected_protein", "translation_identity",
    "terminal_stop_codon_offset_expected", "reconstruction_status", "reconstruction_warning",
]

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def translate(nt: str) -> str:
    nt = nt.upper().replace("U", "T")
    aa = []
    for i in range(0, len(nt) - 2, 3):
        aa.append(CODON_TABLE.get(nt[i:i + 3], "X"))
    return "".join(aa)


def load_cds_fasta(cds_fasta_dir: Optional[Path]) -> Dict[str, str]:
    """Map transcript_id (normalized) -> CDS nucleotide sequence, from any FASTA in dir."""
    seqs: Dict[str, str] = {}
    if not cds_fasta_dir or not Path(cds_fasta_dir).exists():
        return seqs
    for fa in Path(cds_fasta_dir).rglob("*.f*a"):
        cur, buf = None, []
        with open(fa, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        seqs[cur] = "".join(buf)
                    head = line[1:].strip().split()[0]
                    cur, buf = norm_tx(head.split("|")[0]), []
                else:
                    buf.append(line.strip())
            if cur:
                seqs[cur] = "".join(buf)
    return seqs


def reconstruction_audit(coord_rows, cds_by_tx, prot_len, prot_seq, cds_nt):
    out = []
    for r in coord_rows:
        tx = norm_tx(r.get("transcript_id_source"))
        pid = (r.get("protein_id") or "").strip()
        blocks = cds_by_tx.get(tx, [])
        src = (r.get("resolver_source_db") or "").strip()
        n = len(blocks)
        nt_len = sum(_int(b.get("cds_length_bp"), 0) or 0 for b in blocks)
        recon_prot = nt_len // 3
        sel_len = prot_len.get(pid)
        stop_offset = ""
        identity = "not_computed_no_cds_sequence"
        status, warn = "", ""
        matches = ""

        if n == 0:
            status = "transcript_not_found_in_cds_model"
            warn = "no CDS blocks for transcript in local cds_features model"
        elif nt_len <= 0:
            status = "cds_sequence_unavailable"
        else:
            # real translation if CDS nucleotides are available
            ntseq = cds_nt.get(tx)
            if ntseq:
                prot = translate(ntseq).rstrip("*")
                ref = (prot_seq.get(pid) or "").rstrip("*")
                if ref:
                    L = min(len(prot), len(ref))
                    same = sum(1 for a, b in zip(prot[:L], ref[:L]) if a == b)
                    identity = round(same / max(1, len(ref)), 4)
                    matches = "true" if identity >= 0.99 else "false"
                    if matches == "true":
                        status = "cds_reconstruction_matches_protein"
                    else:
                        status = "cds_reconstruction_length_mismatch_review"
            if not status:
                # coordinate/length consistency proxy (no nucleotides locally)
                if sel_len is None or sel_len <= 0:
                    status = "cds_sequence_unavailable"
                    warn = "selected protein length unavailable"
                else:
                    diff = recon_prot - sel_len
                    stop_offset = "true" if diff in (1,) else "false"
                    if diff in (0, 1):
                        matches = "true"
                        status = ("cds_reconstruction_matches_with_terminal_stop_offset"
                                  if diff == 1 else "cds_reconstruction_matches_protein")
                    else:
                        matches = "false"
                        status = "cds_reconstruction_length_mismatch_review"
                        warn = f"reconstructed_protein_len({recon_prot}) vs selected({sel_len})"
        out.append({
            "species": r.get("species_canonical", ""), "isoform": r.get("inferred_isoform", ""),
            "transcript_id": tx, "protein_id": pid, "source_annotation": src,
            "strand": r.get("resolver_strand", ""), "n_cds_blocks": n,
            "reconstructed_cds_nt_length": nt_len if nt_len > 0 else "",
            "reconstructed_protein_length": recon_prot if nt_len > 0 else "",
            "selected_protein_length": sel_len if sel_len else "",
            "translation_matches_selected_protein": matches,
            "translation_identity": identity,
            "terminal_stop_codon_offset_expected": stop_offset,
            "reconstruction_status": status, "reconstruction_warning": warn,
        })
    return out


# ---------------------------------------------------------------------------
# PART E — cassette coordinate sanity audit (hard biological checks)
# ---------------------------------------------------------------------------
SANITY_COLS = [
    "species", "isoform", "protein_id", "protein_length", "cassette_start_aa",
    "cassette_end_aa", "cds_rank", "coordinate_sanity_status", "coordinate_sanity_warning",
]


def sanity_audit(coord_rows, mapping_by_key, prot_len, is_review_fn):
    out = []
    for r in coord_rows:
        sp = r.get("species_canonical", "")
        iso = r.get("inferred_isoform", "")
        pid = (r.get("protein_id") or "").strip()
        m = mapping_by_key.get((sp.lower(), iso), {})
        cas_start = _int(m.get("matched_protein_start_aa"))
        cas_end = _int(m.get("matched_protein_end_aa"))
        cds_rank = m.get("matched_cds_rank", "")
        plen = prot_len.get(pid) or _int(r.get("protein_length_aa"))
        review = is_review_fn(sp, iso)
        status, warn = "cassette_coordinate_plausible", ""

        if cas_start is None:
            status, warn = "cassette_coordinate_unresolved", "no mapped cassette start"
        else:
            full = plen is not None and plen > FULL_LENGTH_AA
            if full and cas_start < MIN_PLAUSIBLE_CASSETTE_START:
                status = "cassette_coordinate_implausible_hard_fail"
                warn = (f"cassette_start_aa={cas_start} < {MIN_PLAUSIBLE_CASSETTE_START} "
                        f"in full-length protein ({plen} aa)")
            elif str(cds_rank).strip() == "1" and full:
                status = "cassette_first_cds_block_implausible"
                warn = f"cassette mapped to CDS rank 1 in full-length protein ({plen} aa)"
            if review and status != "cassette_coordinate_plausible":
                status = status + "_review_excluded"
        out.append({
            "species": sp, "isoform": iso, "protein_id": pid,
            "protein_length": plen if plen else "",
            "cassette_start_aa": cas_start if cas_start is not None else "",
            "cassette_end_aa": cas_end if cas_end is not None else "",
            "cds_rank": cds_rank, "coordinate_sanity_status": status,
            "coordinate_sanity_warning": warn,
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Fix cassette->CDS-block mapping by coordinate overlap.")
    ap.add_argument("--coordinate_audit", type=Path, required=True)
    ap.add_argument("--cds_features", type=Path, required=True)
    ap.add_argument("--proteins", type=Path, default=None)
    ap.add_argument("--review_master", type=Path, default=None,
                    help="species_qc_master.tsv to flag review species in the sanity audit")
    ap.add_argument("--cds_fasta_dir", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    coord = read_tsv(args.coordinate_audit)
    prot_len, prot_seq = parse_proteins(args.proteins)
    cds_nt = load_cds_fasta(args.cds_fasta_dir)

    cds_by_tx: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    tx_meta: Dict[str, Dict[str, str]] = {}
    for c in read_tsv(args.cds_features):
        tx = norm_tx(c.get("transcript_id_source"))
        cds_by_tx[tx].append(c)
        tx_meta.setdefault(tx, {"protein_id": c.get("translation_id_source", "")})
    for tx in cds_by_tx:
        cds_by_tx[tx].sort(key=lambda b: _int(b.get("cds_rank"), 0) or 0)

    review_keys = set()
    if args.review_master and Path(args.review_master).exists():
        for r in read_tsv(args.review_master):
            cls = str(r.get("final_display_class", ""))
            if "review" in cls or cls.startswith("supplementary"):
                for iso in ("IIIb", "IIIc"):
                    review_keys.add((str(r.get("species", "")).lower(), iso))

    def is_review_fn(sp, iso):
        return (sp.lower(), iso) in review_keys

    args.outdir.mkdir(parents=True, exist_ok=True)

    # PART A
    uniq = build_unique_blocks(cds_by_tx, tx_meta)
    write_tsv(args.outdir / "fgfr2_unique_cds_block_table.tsv", uniq, UNIQUE_COLS)

    # PART B
    mapping = [map_cassette(r, cds_by_tx.get(norm_tx(r.get("transcript_id_source")), []))
               for r in coord]
    write_tsv(args.outdir / "fgfr2_cassette_cds_block_map.tsv", mapping, MAP_COLS)
    mapping_by_key = {(m["species"].lower(), m["isoform"]): m for m in mapping}

    # PART C
    recon = reconstruction_audit(coord, cds_by_tx, prot_len, prot_seq, cds_nt)
    write_tsv(args.outdir / "fgfr2_transcript_cds_reconstruction_audit.tsv", recon, RECON_COLS)

    # PART E
    sanity = sanity_audit(coord, mapping_by_key, prot_len, is_review_fn)
    write_tsv(args.outdir / "fgfr2_cassette_coordinate_sanity_audit.tsv", sanity, SANITY_COLS)

    # console summary
    from collections import Counter
    print(f"[OK] unique CDS blocks: {len(uniq)} rows")
    print(f"     cassette_overlap_status={dict(Counter(m['cassette_overlap_status'] for m in mapping))}")
    print(f"     mapped cassette start_aa==1: "
          f"{sum(1 for m in mapping if str(m.get('matched_protein_start_aa')) == '1')}")
    print(f"     reconstruction_status={dict(Counter(r['reconstruction_status'] for r in recon))}")
    print(f"     coordinate_sanity_status={dict(Counter(s['coordinate_sanity_status'] for s in sanity))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
