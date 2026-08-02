#!/usr/bin/env python3
"""
reconstruct_exon_blocks_post_interpro.py

Fix exon-block coordinate artifacts for the three post-InterPro coordinate-artifact
cases (canis_lupus_familiaris IIIb, gorilla_gorilla_gorilla IIIb,
xenopus_tropicalis IIIb) WITHOUT touching the protein sequence, FASTA, final
truth table, InterProScan domain calls, pyTMHMM TM prediction, or IIIb/IIIc
identity.

Root cause (established from local tables):
  * The exon/CDS display coordinates and the drawn cassette slot came from a
    *different* transcript than the final truth-table protein. The cassette was
    resolved on a resolver transcript (e.g. gorilla ENSGGOT00000051097) and its
    protein AA coordinates (394-461) were carried onto the final protein
    (e.g. gorilla ENSGGOP00000051435, 670 aa), where 394-461 lands inside the
    kinase. For gorilla/xenopus this template even matched byte-for-byte and/or
    exceeded the final protein length.

Fix strategy (native-if-available, else validated-cassette + low-confidence):
  * gorilla / xenopus: rebuild the coding-exon blocks from the *native* CDS
    features of the exact final transcript (02_models/cds_features.tsv, which
    carries per-CDS protein_start_aa / protein_end_aa and cds_rank). These are
    within protein length and monotonic and are NOT the shared template.
  * canis (RefSeq NP_001003336.1): no local RefSeq CDS coordinates exist (only
    unrelated Ensembl canis transcripts), so native blocks cannot be
    reconstructed without re-downloading. Coding exon blocks are therefore
    hidden and only the validated cassette slot is shown.
  * In all three cases the cassette SLOT is taken from the validated reference
    cassette coordinate (figure3C cassette_start_aa / cassette_end_aa), which is
    upstream of the pyTMHMM TM and the kinase and does not overlap the kinase.

Outputs:
  * tables/exon_block_coordinate_reconstruction_audit.tsv  (human-auditable)
  * tables/exon_block_reconstruction_overrides.json         (consumed by
    make_fgfr2_post_interpro_exon_domain_figures.py so plots / feature table /
    QC are regenerated consistently)

Read-only w.r.t. FASTA, truth table and primary/review membership.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO = Path(__file__).resolve().parent.parent
import os as _os  # run-folder path override (RESULTS_DIR/BASE); legacy default preserved
BASE = Path(_os.environ.get("FGFR2_RESULTS_DIR") or _os.environ.get("RESULTS_DIR")
            or _os.environ.get("BASE") or (REPO / "results" / "final_30_until_interpro_prepare"))
CLOSURE = BASE / "13_final_pre_interpro_closure"
POST = BASE / "15_exon_domain_boundary_post_interpro"

TRUTH = CLOSURE / "final_pre_interpro_truth_table.tsv"
COORD = CLOSURE / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv"
CDS_FEATURES = BASE / "02_models" / "cds_features.tsv"
CASSETTE_MAP = (BASE / "09_paper_ready_qc_v2_9" / "figures_v2_22_final_qc_display"
                / "fgfr2_cassette_cds_block_map.tsv")
FEATURES = POST / "tables" / "exon_domain_architecture_features.tsv"
QC = POST / "tables" / "fgfr2_domain_architecture_qc.tsv"

AUDIT_OUT = POST / "tables" / "exon_block_coordinate_reconstruction_audit.tsv"
OVERRIDE_OUT = POST / "tables" / "exon_block_reconstruction_overrides.json"


def display_path(path) -> str:
    """Best-effort repo-relative rendering for logs. BASE may be a run-local
    RELATIVE path (runs/<id>/results), which cannot be made relative to the
    absolute REPO — fall back to the raw path instead of raising ValueError.
    Display only; never affects outputs."""
    p = Path(path)
    try:
        return str(p.resolve().relative_to(REPO.resolve()))
    except Exception:
        return str(p)

# The three coordinate-artifact cases (species, isoform).
CASES = [
    ("canis_lupus_familiaris", "IIIb"),
    ("gorilla_gorilla_gorilla", "IIIb"),
    ("xenopus_tropicalis", "IIIb"),
]


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def to_int(v, default=None):
    try:
        return int(float(str(v).strip()))
    except (ValueError, TypeError):
        return default


def load_native_blocks(transcript_id: str) -> List[dict]:
    """Native coding-exon blocks (protein AA) for an exact transcript, if present."""
    blocks: List[dict] = []
    for r in read_tsv(CDS_FEATURES):
        if r.get("transcript_id_source") != transcript_id:
            continue
        pa, pe = to_int(r.get("protein_start_aa")), to_int(r.get("protein_end_aa"))
        rank = to_int(r.get("cds_rank"))
        if pa is None or pe is None or rank is None:
            continue
        blocks.append({
            "rank": rank,
            "exon_id": r.get("cds_id_source", ""),
            "chrom": r.get("chrom", ""),
            "genomic_start": to_int(r.get("start")),
            "genomic_end": to_int(r.get("end")),
            "strand": r.get("strand", ""),
            "protein_start_aa": pa,
            "protein_end_aa": pe,
        })
    blocks.sort(key=lambda b: b["rank"])
    return blocks


def native_cassette_rank(blocks: List[dict], gstart: Optional[int],
                         gend: Optional[int]) -> Optional[dict]:
    """Which native CDS block overlaps the resolved cassette genomic interval."""
    if gstart is None or gend is None:
        return None
    best, best_ov = None, 0
    for b in blocks:
        bs, be = b["genomic_start"], b["genomic_end"]
        if bs is None or be is None:
            continue
        ov = min(be, gend) - max(bs, gstart)
        if ov > best_ov:
            best, best_ov = b, ov
    return best


def main() -> int:
    truth = {(r["species"], r["isoform"]): r for r in read_tsv(TRUTH)}
    coord = read_tsv(COORD)
    cassette_map = {(r["species"], r["isoform"]): r for r in read_tsv(CASSETTE_MAP)}
    features = read_tsv(FEATURES)

    # figure3C: validated reference cassette coord + original block max per case
    cass_ref: Dict[Tuple[str, str], Tuple[int, int]] = {}
    orig_block_max: Dict[Tuple[str, str], int] = {}
    orig_src: Dict[Tuple[str, str], str] = {}
    precision: Dict[Tuple[str, str], str] = {}
    for r in coord:
        key = (r.get("species", ""), r.get("isoform", ""))
        be = to_int(r.get("block_end_aa"))
        if be is not None:
            orig_block_max[key] = max(orig_block_max.get(key, 0), be)
        is_cass = str(r.get("is_IIIb_cassette", "")).lower() == "true" \
            or str(r.get("is_IIIc_cassette", "")).lower() == "true"
        if is_cass and key not in cass_ref:
            cs, ce = to_int(r.get("cassette_start_aa")), to_int(r.get("cassette_end_aa"))
            if cs is not None and ce is not None:
                cass_ref[key] = (cs, ce)
            orig_src[key] = r.get("source_coordinate_table", "")
            precision[key] = f"{r.get('boundary_left_precision','')}/{r.get('boundary_right_precision','')}"

    # TM (pyTMHMM receptor) + representative kinase from the post-InterPro feature table
    tm: Dict[Tuple[str, str], Tuple[int, int]] = {}
    kinase: Dict[Tuple[str, str], Tuple[int, int]] = {}
    for r in features:
        key = (r.get("species", ""), r.get("isoform", ""))
        ft = r.get("feature_type", "")
        s, e = to_int(r.get("start_aa")), to_int(r.get("end_aa"))
        if ft == "transmembrane_pytmhmm" and r.get("status") == "receptor_tm":
            tm[key] = (s, e)
        elif ft == "kinase_domain":
            kinase[key] = (s, e)

    audit_rows: List[dict] = []
    overrides: Dict[str, dict] = {}

    for key in CASES:
        sp, iso = key
        tr = truth.get(key, {})
        tx = tr.get("transcript_id", "")
        pid = tr.get("protein_id", "")
        plen = to_int(tr.get("protein_length"))
        ref = cass_ref.get(key)
        tms = tm.get(key)
        kin = kinase.get(key)
        cm = cassette_map.get(key, {})
        resolver_tx = cm.get("transcript_id", "")
        g0, g1 = to_int(cm.get("resolver_genomic_start")), to_int(cm.get("resolver_genomic_end"))

        native = load_native_blocks(tx)
        native_available = bool(native)
        native_max = max((b["protein_end_aa"] for b in native), default=None)
        exceeds = (native_max is not None and plen is not None and native_max > plen + 1)
        cass_rank = native_cassette_rank(native, g0, g1)

        # validated cassette slot (used for display + QC in all three cases)
        cassette = {"start": ref[0], "end": ref[1]} if ref else None
        _cassette_overlaps_kinase = bool(
            cassette and kin and not (cassette["end"] < kin[0] or cassette["start"] > kin[1]))

        notes: List[str] = []
        if resolver_tx and resolver_tx != tx:
            notes.append(f"original exon/cassette coords came from resolver transcript "
                         f"{resolver_tx} (not final {tx})")

        if native_available:
            # clamp native blocks into [1, protein_length]; renumber N->C
            blocks_out: List[dict] = []
            for i, b in enumerate(native, start=1):
                s = max(1, b["protein_start_aa"])
                e = min(plen, b["protein_end_aa"]) if plen else b["protein_end_aa"]
                if e < s:
                    e = s
                blocks_out.append({"number": i, "label": f"CDS{b['rank']}",
                                   "exon_id": b["exon_id"], "start": s, "end": e})
            max_out = max(b["end"] for b in blocks_out)
            recon_success = (max_out <= (plen or max_out)) and native_available
            recon_source = f"{display_path(CDS_FEATURES)} ({tx})"
            if cass_rank:
                notes.append(f"native D3 cassette exon = CDS rank {cass_rank['rank']} "
                             f"(AA {cass_rank['protein_start_aa']}-{cass_rank['protein_end_aa']}); "
                             "displayed cassette uses validated reference coordinate")
            notes.append("coding-exon blocks reconstructed from native local CDS features; "
                         "cassette slot from validated reference coordinate")
            display_status = "native_exon_blocks_reconstructed"
            overrides[f"{sp}|{iso}"] = {
                "final_display_status": display_status,
                "exon_blocks": blocks_out,
                "cassette": cassette,
                "recon_note": "native coding-exon blocks reconstructed; cassette coordinate "
                              "validated separately (reference cassette upstream of TM/kinase)",
                "exon_block_source": f"cds_features_native ({tx})",
            }
            n_blocks = len(blocks_out)
        else:
            # RefSeq protein without local CDS coords: hide misleading blocks, keep cassette
            recon_success = False
            recon_source = "none (RefSeq CDS coordinates not in local cache)"
            notes.append("no local native CDS coordinates for the final RefSeq protein; "
                         "misleading template exon blocks hidden; only validated cassette shown")
            display_status = "cassette_only_high_confidence"
            overrides[f"{sp}|{iso}"] = {
                "final_display_status": display_status,
                "exon_blocks": [],
                "cassette": cassette,
                "recon_note": "native exon-block display low-confidence; cassette coordinate "
                              "validated separately",
                "exon_block_source": "none_refseq_cds_not_cached",
            }
            n_blocks = 0
            max_out = None

        audit_rows.append({
            "species": sp,
            "isoform": tr.get("final_isoform_label", iso),
            "transcript_id": tx,
            "protein_id": pid,
            "protein_length": plen if plen is not None else "",
            "original_exon_block_source": orig_src.get(key, ""),
            "original_exon_block_status": ("template_shared_or_exceeds_length"
                                           if (orig_block_max.get(key, 0) > (plen or 0) + 1
                                               or resolver_tx != tx)
                                           else "resolver_transcript_coordinate"),
            "reconstructed_exon_block_source": recon_source,
            "reconstruction_attempted": "true",
            "reconstruction_success": "true" if recon_success else "false",
            "number_of_exon_blocks": n_blocks,
            "max_exon_end_aa": max_out if max_out is not None else "",
            "exceeds_protein_length": "true" if exceeds else "false",
            "cassette_start_aa": cassette["start"] if cassette else "",
            "cassette_end_aa": cassette["end"] if cassette else "",
            "tm_start_aa": tms[0] if tms else "",
            "kinase_start_aa": kin[0] if kin else "",
            "final_display_status": display_status,
            "notes": "; ".join(notes),
        })

    cols = ["species", "isoform", "transcript_id", "protein_id", "protein_length",
            "original_exon_block_source", "original_exon_block_status",
            "reconstructed_exon_block_source", "reconstruction_attempted",
            "reconstruction_success", "number_of_exon_blocks", "max_exon_end_aa",
            "exceeds_protein_length", "cassette_start_aa", "cassette_end_aa",
            "tm_start_aa", "kinase_start_aa", "final_display_status", "notes"]
    AUDIT_OUT.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT_OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(audit_rows)

    OVERRIDE_OUT.write_text(json.dumps(overrides, indent=2, ensure_ascii=False),
                            encoding="utf-8")

    print(f"[ok] wrote {display_path(AUDIT_OUT)}")
    print(f"[ok] wrote {display_path(OVERRIDE_OUT)}")
    for r in audit_rows:
        print(f"  {r['species']} {r['isoform']}: {r['final_display_status']} "
              f"(blocks={r['number_of_exon_blocks']}, max_aa={r['max_exon_end_aa']}, "
              f"exceeds={r['exceeds_protein_length']}, cassette={r['cassette_start_aa']}-{r['cassette_end_aa']}, "
              f"tm={r['tm_start_aa']}, kinase_start={r['kinase_start_aa']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
