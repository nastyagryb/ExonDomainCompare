#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

SCRIPT_VERSION = "1.0"
FGFR_GENES = ("FGFR1", "FGFR2", "FGFR3", "FGFR4")


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    hdr = None
    seq: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(seq)))
            hdr = line[1:].strip()
            seq = []
        elif line.strip():
            seq.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(seq)))
    return out


def parse_query_header(header: str) -> Dict[str, str]:
    parts = header.split("|")
    kv: Dict[str, str] = {"fasta_id": parts[0].strip() if parts else ""}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            kv[k.strip().lower()] = v.strip()
    return {
        "fasta_id": kv.get("fasta_id", ""),
        "species": kv.get("species", ""),
        "role": kv.get("role", ""),
        "transcript_id": kv.get("transcript", kv.get("transcript_id", "")),
        "protein_id": kv.get("protein", kv.get("protein_id", "")),
        "isoform": kv.get("isoform", ""),
    }


def parse_panel_sseqid(sseqid: str) -> Dict[str, str]:
    parts = sseqid.split("|")
    kv: Dict[str, str] = {}
    gene = parts[0].strip() if parts else ""
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            kv[k.strip().lower()] = v.strip()
    return {"gene": gene, "species": kv.get("species", ""), "accession": kv.get("accession", ""),
            "source": kv.get("source", "")}


def run_blast(query: Path, panel: Path, threads: int, workdir: Path) -> List[Dict[str, str]]:
    db = workdir / "panel_db"
    subprocess.run(["makeblastdb", "-in", str(panel), "-dbtype", "prot", "-out", str(db)],
                   check=True, capture_output=True)
    out_tsv = workdir / "blast.tsv"
    fmt = "6 qseqid sseqid pident length qlen slen bitscore evalue"
    subprocess.run(["blastp", "-query", str(query), "-db", str(db), "-outfmt", fmt,
                    "-num_threads", str(threads), "-max_target_seqs", "50", "-out", str(out_tsv)],
                   check=True, capture_output=True)
    rows: List[Dict[str, str]] = []
    if out_tsv.exists():
        for line in out_tsv.read_text().splitlines():
            f = line.split("\t")
            if len(f) < 8:
                continue
            rows.append({"qseqid": f[0], "sseqid": f[1], "pident": f[2], "length": f[3],
                         "qlen": f[4], "slen": f[5], "bitscore": f[6], "evalue": f[7]})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description="Screen candidates against the multi-vertebrate FGFR panel.")
    ap.add_argument("--query_fasta", type=Path, required=True)
    ap.add_argument("--panel_fasta", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min_margin", type=float, default=0.05, help="min FGFR2 bitscore margin over next-best gene (fraction).")
    ap.add_argument("--min_coverage", type=float, default=0.6)
    args = ap.parse_args()

    if not shutil.which("makeblastdb") or not shutil.which("blastp"):
        raise SystemExit("[ERROR] blastp/makeblastdb not found on PATH.")

    args.outdir.mkdir(parents=True, exist_ok=True)
    _queries = {h.split()[0].split("|")[0]: parse_query_header(h) for h, _ in parse_fasta(args.query_fasta)}
    query_meta = {h.split("|")[0]: parse_query_header(h) for h, _ in parse_fasta(args.query_fasta)}

    with tempfile.TemporaryDirectory() as td:
        blast_rows = run_blast(args.query_fasta, args.panel_fasta, args.threads, Path(td))

    # Group hits by query.
    by_q: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for r in blast_rows:
        by_q[r["qseqid"].split("|")[0]].append(r)

    detailed: List[Dict[str, object]] = []
    warnings: List[Dict[str, str]] = []
    seen_protein: set = set()

    for fasta_id, qmeta in query_meta.items():
        prot = qmeta.get("protein_id") or fasta_id
        dedup_key = (qmeta.get("species"), qmeta.get("isoform"), prot)
        if dedup_key in seen_protein:
            continue
        seen_protein.add(dedup_key)

        hits = by_q.get(fasta_id, [])
        # best hit per gene by bitscore.
        per_gene_best: Dict[str, Dict[str, object]] = {}
        per_gene_groups: Dict[str, set] = defaultdict(set)
        for h in hits:
            pinfo = parse_panel_sseqid(h["sseqid"])
            gene = pinfo["gene"]
            if gene not in FGFR_GENES:
                continue
            try:
                bit = float(h["bitscore"]); pid = float(h["pident"]); ln = float(h["length"]); ql = float(h["qlen"])
            except ValueError:
                continue
            cov = ln / ql if ql else 0.0
            per_gene_groups[gene].add(pinfo["species"])
            cur = per_gene_best.get(gene)
            if cur is None or bit > float(cur["best_hit_score"]):
                per_gene_best[gene] = {
                    "best_hit_score": bit, "best_hit_identity": round(pid / 100.0, 4),
                    "best_hit_coverage": round(cov, 4), "best_paralog_species": pinfo["species"],
                    "best_paralog_accession": pinfo["accession"],
                }

        row: Dict[str, object] = {
            "species": qmeta.get("species", ""), "isoform": qmeta.get("isoform", ""),
            "transcript_id": qmeta.get("transcript_id", ""), "protein_id": prot,
        }

        if not per_gene_best:
            row.update({
                "best_paralog_gene": "", "best_paralog_species": "", "best_paralog_accession": "",
                "best_hit_identity": "", "best_hit_coverage": "", "best_hit_score": "",
                "next_best_paralog_gene": "", "fgfr2_margin_over_next_best": "",
                "paralog_confidence": "none", "paralog_status": "paralog_evidence_unavailable",
                "paralog_warning": "no_fgfr_panel_hits",
            })
            detailed.append(row)
            warnings.append({"species": row["species"], "protein_id": prot, "warning": "no_fgfr_panel_hits"})
            continue

        # Overall best gene + next-best gene.
        genes_sorted = sorted(per_gene_best.items(), key=lambda kv: float(kv[1]["best_hit_score"]), reverse=True)
        best_gene, best = genes_sorted[0]
        next_gene = genes_sorted[1][0] if len(genes_sorted) > 1 else ""
        _next_score = float(genes_sorted[1][1]["best_hit_score"]) if len(genes_sorted) > 1 else 0.0

        fgfr2 = per_gene_best.get("FGFR2")
        fgfr2_score = float(fgfr2["best_hit_score"]) if fgfr2 else 0.0
        # Margin of FGFR2 over the best NON-FGFR2 gene.
        non_fgfr2 = [(g, v) for g, v in genes_sorted if g != "FGFR2"]
        best_non_fgfr2 = float(non_fgfr2[0][1]["best_hit_score"]) if non_fgfr2 else 0.0
        margin = (fgfr2_score - best_non_fgfr2) / fgfr2_score if fgfr2_score > 0 else -1.0

        fgfr2_groups = per_gene_groups.get("FGFR2", set())
        n_groups = len(fgfr2_groups)
        cov_ok = bool(fgfr2) and float(fgfr2["best_hit_coverage"]) >= args.min_coverage

        warning = ""
        if best_gene != "FGFR2":
            status = "non_fgfr2_best_hit_review"
            confidence = "low"
            warning = f"best_panel_gene_is_{best_gene}_not_FGFR2"
        elif not fgfr2 or fgfr2_score <= 0:
            status = "paralog_evidence_unavailable"
            confidence = "none"
        elif margin < 0.01:
            status = "ambiguous_fgfr_paralog_review"
            confidence = "low"
            warning = f"fgfr2_margin_over_next_best_only_{round(margin,4)}"
        elif n_groups <= 1 and "homo_sapiens" in fgfr2_groups:
            status = "fgfr2_supported_human_only"
            confidence = "medium"
            warning = "fgfr2_support_from_human_reference_only"
        elif margin < args.min_margin:
            status = "fgfr2_supported_low_margin"
            confidence = "medium"
            warning = f"low_margin_{round(margin,4)}"
        elif n_groups >= 3:
            status = "fgfr2_high_confidence_multi_vertebrate"
            confidence = "high"
        else:
            status = "fgfr2_supported_low_margin"
            confidence = "medium"
            warning = f"fgfr2_support_from_{n_groups}_vertebrate_group(s)"
        if not cov_ok and status not in ("non_fgfr2_best_hit_review", "paralog_evidence_unavailable"):
            warning = (warning + ";" if warning else "") + f"low_query_coverage_{fgfr2['best_hit_coverage'] if fgfr2 else 'NA'}"

        # Report the BEST paralog hit (overall best gene) plus FGFR2 margin.
        row.update({
            "best_paralog_gene": best_gene,
            "best_paralog_species": best["best_paralog_species"],
            "best_paralog_accession": best["best_paralog_accession"],
            "best_hit_identity": best["best_hit_identity"],
            "best_hit_coverage": best["best_hit_coverage"],
            "best_hit_score": round(float(best["best_hit_score"]), 1),
            "next_best_paralog_gene": next_gene,
            "fgfr2_margin_over_next_best": round(margin, 4),
            "paralog_confidence": confidence,
            "paralog_status": status,
            "paralog_warning": warning,
        })
        detailed.append(row)
        if warning:
            warnings.append({"species": row["species"], "protein_id": prot, "warning": warning})

    detailed.sort(key=lambda r: (str(r["species"]), str(r["isoform"])))
    det_cols = ["species", "isoform", "transcript_id", "protein_id", "best_paralog_gene",
                "best_paralog_species", "best_paralog_accession", "best_hit_identity",
                "best_hit_coverage", "best_hit_score", "next_best_paralog_gene",
                "fgfr2_margin_over_next_best", "paralog_confidence", "paralog_status", "paralog_warning"]
    with open(args.outdir / "fgfr2_paralog_screen_detailed.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=det_cols)
        w.writeheader()
        for r in detailed:
            w.writerow({k: r.get(k, "") for k in det_cols})

    # Species summary.
    by_sp: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for r in detailed:
        by_sp[str(r["species"])].append(r)
    summary: List[Dict[str, object]] = []
    _RANK = {"fgfr2_high_confidence_multi_vertebrate": 5, "fgfr2_supported_low_margin": 4,
             "fgfr2_supported_human_only": 3, "ambiguous_fgfr_paralog_review": 2,
             "non_fgfr2_best_hit_review": 1, "paralog_evidence_unavailable": 0}
    for sp, rs in sorted(by_sp.items()):
        statuses = [str(r["paralog_status"]) for r in rs]
        worst = min(statuses, key=lambda s: _RANK.get(s, 0)) if statuses else "paralog_evidence_unavailable"
        margins = [float(r["fgfr2_margin_over_next_best"]) for r in rs if str(r["fgfr2_margin_over_next_best"]) not in ("", "None")]
        summary.append({
            "species": sp, "n_proteins": len(rs),
            "n_fgfr2_high_confidence": sum(1 for s in statuses if s == "fgfr2_high_confidence_multi_vertebrate"),
            "n_review": sum(1 for s in statuses if "review" in s or s == "paralog_evidence_unavailable"),
            "min_fgfr2_margin_over_next_best": round(min(margins), 4) if margins else "",
            "best_paralog_gene_set": ";".join(sorted({str(r["best_paralog_gene"]) for r in rs if r["best_paralog_gene"]})),
            "species_fgfr2_screen_status": worst,
        })
    with open(args.outdir / "fgfr2_paralog_screen_species_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        cols = ["species", "n_proteins", "n_fgfr2_high_confidence", "n_review",
                "min_fgfr2_margin_over_next_best", "best_paralog_gene_set", "species_fgfr2_screen_status"]
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        w.writerows(summary)

    with open(args.outdir / "fgfr2_paralog_screen_warnings.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=["species", "protein_id", "warning"])
        w.writeheader()
        w.writerows(warnings)

    meta = {
        "script_version": SCRIPT_VERSION,
        "panel_fasta": str(args.panel_fasta),
        "n_query_proteins_screened": len(detailed),
        "paralog_status_counts": dict(Counter(str(r["paralog_status"]) for r in detailed)),
        "species_status_counts": dict(Counter(str(s["species_fgfr2_screen_status"]) for s in summary)),
        "n_warnings": len(warnings),
    }
    (args.outdir / "fgfr2_paralog_screen_multi_vertebrate_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"[OK] multi-vertebrate paralog screen: {len(detailed)} proteins; statuses={meta['paralog_status_counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
