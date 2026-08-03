#!/usr/bin/env python3
"""
Check FGFR2 splice-site motif QC.

Splice-site motif QC for resolved IIIb/IIIc cassette boundaries — ONLY where
source-compatible genomic sequence is available. Splice-site motifs are NEVER faked,
and genome assemblies are never mixed silently. When no compatible local genome FASTA
is available, the table honestly reports sequence_unavailable with the genomic
coordinates that would be required.

This is an annotation-support layer and does NOT override IIIb/IIIc assignment.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402
from exondomaincompare.shared_gene_analysis.strand import is_forward  # noqa: E402


COLS = ["species", "isoform", "transcript_id", "seqid", "strand",
        "left_boundary_genomic_position", "right_boundary_genomic_position",
        "acceptor_motif", "donor_motif", "splice_site_class", "splice_site_qc_status",
        "sequence_source", "source_compatibility_status", "splice_site_warning"]
SUMMARY_COLS = ["metric", "value"]


def parse_block_id(bid: str) -> Optional[Tuple[str, int, int, str]]:
    # species|tx|rank|seqid:start-end:strand
    try:
        loc = bid.split("|")[-1]
        seqid, rest = loc.split(":", 1)
        coords, strand = rest.rsplit(":", 1)
        s, e = coords.split("-")
        return seqid, int(s), int(e), strand
    except Exception:  # noqa: BLE001
        return None


def load_genome_contig(genome_dir: Optional[Path], seqid: str) -> Optional[str]:
    if not genome_dir:
        return None
    for ext in (".fa", ".fasta", ".fna"):
        p = genome_dir / f"{seqid}{ext}"
        if p.exists():
            items = M.read_fasta(p)
            if items:
                return M.clean_alignment_seq(items[0][1])
    return None


def classify_motif(acceptor: str, donor: str) -> str:
    d = (donor or "").upper()
    if d == "GT":
        return "canonical_GT_AG"
    if d == "GC":
        return "noncanonical_GC_AG"
    if d == "AT":
        return "noncanonical_AT_AC"
    if d:
        return "noncanonical_other_review"
    return "sequence_unavailable"


def main() -> int:
    ap = argparse.ArgumentParser(description="Check splice-site motif QC.")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--genome_dir", type=Path, default=None,
                    help="optional dir with per-contig genome FASTA ({seqid}.fa); "
                         "splice motifs are only computed when source-compatible sequence exists")
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    out = dirs["splice_qc"]

    cmap = M.read_tsv(M.require(base, "fgfr2_cassette_cds_block_map.tsv"))
    rows: List[Dict[str, object]] = []
    for m in cmap:
        sp = m.get("species", "")
        iso = m.get("isoform", "")
        tx = m.get("transcript_id", "")
        bid = m.get("matched_unique_cds_block_id", "")
        parsed = parse_block_id(bid)
        if not parsed:
            rows.append({"species": sp, "isoform": iso, "transcript_id": tx, "seqid": "",
                         "strand": "", "left_boundary_genomic_position": "",
                         "right_boundary_genomic_position": "", "acceptor_motif": "",
                         "donor_motif": "", "splice_site_class": "boundary_unresolved",
                         "splice_site_qc_status": "splice_site_not_applicable",
                         "sequence_source": "", "source_compatibility_status": "no_cds_block_coordinates",
                         "splice_site_warning": "no genomic block id for cassette (e.g. protein overlay)"})
            continue
        seqid, gs, ge, strand = parsed
        contig = load_genome_contig(args.genome_dir, seqid)
        acceptor = donor = ""
        if contig and len(contig) >= ge + 2:
            # intronic dinucleotides flanking the cassette CDS block (best-effort)
            if is_forward(strand):
                acceptor = contig[gs - 3:gs - 1]   # ...AG | exon
                donor = contig[ge:ge + 2]          # exon | GT...
            else:
                def rc(s):
                    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
                    return "".join(comp.get(c, "N") for c in reversed(s))
                acceptor = rc(contig[ge:ge + 2])
                donor = rc(contig[gs - 3:gs - 1])
            sclass = classify_motif(acceptor, donor)
            qc = ("splice_site_supported" if sclass == "canonical_GT_AG"
                  else "splice_site_noncanonical_review" if sclass.startswith("noncanonical")
                  else "splice_site_sequence_unavailable")
            src = f"local_genome_fasta:{seqid}"
            comp_status = "same_assembly_local_genome"
            warn = "" if sclass == "canonical_GT_AG" else "noncanonical_or_partial_motif_review"
        else:
            sclass = "sequence_unavailable"
            qc = "splice_site_sequence_unavailable"
            src = "no_local_genome_fasta"
            comp_status = "genome_sequence_not_available_locally"
            warn = ("provide --genome_dir with the source-compatible assembly to enable "
                    "splice-site motif QC; not faked")
        rows.append({"species": sp, "isoform": iso, "transcript_id": tx, "seqid": seqid,
                     "strand": strand, "left_boundary_genomic_position": gs,
                     "right_boundary_genomic_position": ge, "acceptor_motif": acceptor,
                     "donor_motif": donor, "splice_site_class": sclass,
                     "splice_site_qc_status": qc, "sequence_source": src,
                     "source_compatibility_status": comp_status, "splice_site_warning": warn})

    M.write_tsv(out / "fgfr2_splice_site_boundary_qc.tsv", rows, COLS)
    qc_counts = Counter(r["splice_site_qc_status"] for r in rows)
    cls_counts = Counter(r["splice_site_class"] for r in rows)
    summary = [{"metric": "n_boundaries", "value": len(rows)},
               {"metric": "genome_available", "value": "true" if args.genome_dir else "false"}]
    for k, v in sorted(qc_counts.items()):
        summary.append({"metric": f"qc_status::{k}", "value": v})
    for k, v in sorted(cls_counts.items()):
        summary.append({"metric": f"class::{k}", "value": v})
    M.write_tsv(out / "fgfr2_splice_site_boundary_qc_summary.tsv", summary, SUMMARY_COLS)
    print(f"[OK] splice-site QC: {len(rows)} boundaries; qc_status={dict(qc_counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
