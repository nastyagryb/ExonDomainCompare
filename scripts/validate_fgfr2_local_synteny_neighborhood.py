#!/usr/bin/env python3
"""
Validate the FGFR2 local synteny neighbourhood.

Independent local synteny / gene-neighborhood validation around the FGFR2 locus. Validates that
the final post-rescue FGFR2 candidates lie in the expected FGFR2 genomic neighborhood (locus /
orthology context) rather than paralogous, fragmented or annotation-artifact loci.

IMPORTANT BIOLOGY:
  * Synteny validates the FGFR2 locus / orthology context only.
  * Synteny does NOT assign or relabel IIIb/IIIc; isoform labels remain sequence-calibrated.
  * Neighbor identities are normalized (symbol > curated orthology > reciprocal best hit >
    high-confidence one-way BLAST > raw ID). BLAST/RBH names are PROBABLE orthologs, not curated.
  * Uses final post-rescue transcript/protein IDs as input; upstream/legacy labels are provenance.
  * No InterProScan, no fake domain annotations.

Source annotation: local NCBI Datasets / RefSeq genomic.gff per assembly (source-compatible per
species). Neighbor protein identity is resolved with DIAMOND (preferred) or BLASTP against a human
FGFR2 5/10-neighbor reference. MCScanX block-level synteny is intentionally NOT part of this layer.

Parts implemented here: A (extraction), B (human reference), C (identity tiers), D (BLAST/RBH),
E (conservation matrix + synteny scoring/classes), I (local 5-neighbor synteny validation gate).
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402
from exondomaincompare.shared_gene_analysis.strand import is_reverse  # noqa: E402
from exondomaincompare.config import load_config  # noqa: E402

HUMAN = "homo_sapiens"
MAIN_N = 5          # main analysis window (each side)
SUPP_N = 10         # supplement window (each side)
GENERIC_PREFIXES = ("LOC", "ORF", "C1ORF", "SI:", "ZGC:", "WU:", "IM:", "BX")

# Shared curated HUMAN reference proteome (SYMBOL|protein_id headers). Broad-homology naming
# of LOC/uncharacterized neighbors needs the whole human RefSeq proteome. In the Example run
# homo_sapiens is an analysed species so its proteome is downloaded run-locally; custom runs
# usually do NOT include human as an analysed species, so we fall back to this shared HUMAN
# reference layer (the only reusable curated reference). It is NEVER a non-human Example output.
REPO_ROOT = Path(__file__).resolve().parent.parent
RUNTIME_CONFIG = load_config(repository_root=REPO_ROOT)
SHARED_HUMAN_PROTEOME = REPO_ROOT / "references" / "synteny" / "human_proteome_named.faa"
# Curated HUMAN FGFR2 gene-neighborhood reference (the 5/10 protein-coding neighbors of human
# FGFR2 + their proteins). Needed so cross-species neighbor symbols can be classified as
# ortholog_supported (matched to the human FGFR2 neighborhood) exactly like the Example run.
# Human reference layer only; never a non-human Example output.
SHARED_HUMAN_NEIGHBORHOOD = REPO_ROOT / "references" / "synteny" / "human_fgfr2_10neighbor_reference.tsv"
SHARED_HUMAN_NEIGHBOR_PROTEINS = (REPO_ROOT / "references" / "synteny" /
                                  "human_fgfr2_10neighbor_reference_proteins.faa")


def _safe_int(v: str, default: int = 0) -> int:
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return default


def load_shared_human_reference() -> "Tuple[Dict[str, object], Dict[str, str]]":
    """Reconstruct the human FGFR2 neighborhood dict + protein map from the shared curated
    HUMAN reference files. Used when the run has no run-local human genome (human is not an
    analysed species) so neighbor identity/orthology classification matches the Example."""
    if not SHARED_HUMAN_NEIGHBORHOOD.exists():
        return {}, {}
    neighbors: List[Dict[str, object]] = []
    for r in M.read_tsv(SHARED_HUMAN_NEIGHBORHOOD):
        neighbors.append({
            "rank": _safe_int(r.get("human_neighbor_rank")),
            "side": r.get("human_neighbor_side", ""),
            "gene_id": r.get("human_gene_id", ""),
            "symbol": r.get("human_gene_symbol", ""),
            "protein_id": r.get("human_protein_id", ""),
            "start": _safe_int(r.get("human_start")),
            "end": _safe_int(r.get("human_end")),
            "strand": r.get("human_strand", ""),
        })
    hproteins: Dict[str, str] = {}
    if SHARED_HUMAN_NEIGHBOR_PROTEINS.exists():
        for sid, seq in M.read_fasta(SHARED_HUMAN_NEIGHBOR_PROTEINS):
            pid = sid.split("|")[-1] if "|" in sid else sid
            hproteins[pid] = seq
    return {"neighbors": neighbors, "status": "human_reference_control",
            "source_file": str(SHARED_HUMAN_NEIGHBORHOOD)}, hproteins

# standard vs stricter (close mammal/primate) neighbor-hit thresholds
THR_STD = {"evalue": 1e-5, "pid": 30.0, "qcov": 0.50, "scov": 0.50}
THR_STRICT = {"evalue": 1e-5, "pid": 60.0, "qcov": 0.70, "scov": 0.70}
CLOSE_GROUPS = {"human_curated_positive_control", "close_primate_control",
                "known_label_risk_mammal"}


# ---------------------------------------------------------------------------
# tool discovery
# ---------------------------------------------------------------------------
def find_diamond() -> Optional[str]:
    return RUNTIME_CONFIG.executable("diamond")


def find_blastp() -> Optional[str]:
    return RUNTIME_CONFIG.executable("blastp")


def find_makeblastdb() -> Optional[str]:
    return RUNTIME_CONFIG.executable("makeblastdb")


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------
def gff_path(base: Path, taxid: str, acc: str) -> Optional[Path]:
    p = (base / "02_models" / "_ncbi_datasets_cache" / f"ncbi_{taxid}" / acc /
         "unzipped" / "ncbi_dataset" / "data" / acc / "genomic.gff")
    return p if p.exists() else None


def protein_faa_path(base: Path, taxid: str, acc: str) -> Optional[Path]:
    p = (base / "02_models" / "_ncbi_datasets_cache" / f"ncbi_{taxid}" / acc /
         "unzipped" / "ncbi_dataset" / "data" / acc / "protein.faa")
    return p if p.exists() else None


def _attrs(field9: str) -> Dict[str, str]:
    out = {}
    for tok in field9.rstrip().split(";"):
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def _geneid(attrs: Dict[str, str]) -> str:
    dbx = attrs.get("Dbxref", "")
    for part in dbx.split(","):
        if part.startswith("GeneID:"):
            return part
    return attrs.get("ID", "")


def parse_gff_for_fgfr2(gff: Path) -> Dict[str, object]:
    """Single streaming pass: collect protein-coding gene features genome-wide and a
    GeneID->protein_id map (first CDS protein per gene). Then localize the FGFR2 gene and the
    surrounding protein-coding genes on the same seqid. Returns a dict with fgfr2 + neighbors."""
    genes: List[Tuple[str, int, int, str, str, str]] = []  # seqid,start,end,strand,geneid,symbol
    pid_by_gene: Dict[str, str] = {}
    fgfr2_geneid = ""
    with open(gff, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line and line[0] == "#":
                continue
            # cheap pre-filter before splitting
            if "\tgene\t" in line:
                f = line.split("\t")
                if len(f) < 9 or f[2] != "gene":
                    continue
                a = _attrs(f[8])
                if a.get("gene_biotype") != "protein_coding":
                    continue
                gid = _geneid(a)
                sym = a.get("gene") or a.get("Name") or gid
                genes.append((f[0], int(f[3]), int(f[4]), f[6], gid, sym))
                if sym.upper() == "FGFR2" and not fgfr2_geneid:
                    fgfr2_geneid = gid
            elif "\tCDS\t" in line:
                # cheap extraction of GeneID + protein_id without full attr dict
                if "protein_id=" not in line:
                    continue
                gi = line.find("GeneID:")
                if gi < 0:
                    continue
                j = gi + 7
                k = j
                while k < len(line) and line[k].isdigit():
                    k += 1
                gid = "GeneID:" + line[j:k]
                if gid in pid_by_gene:
                    continue
                pi = line.find("protein_id=")
                ps = pi + 11
                pe = ps
                while pe < len(line) and line[pe] not in ";\n\r\t":
                    pe += 1
                pid_by_gene[gid] = line[ps:pe]
    return _localize(genes, pid_by_gene, fgfr2_geneid)


def _localize(genes, pid_by_gene, fgfr2_geneid) -> Dict[str, object]:
    out: Dict[str, object] = {"status": "fgfr2_locus_not_found", "neighbors": [], "fgfr2": None}
    frow = None
    if fgfr2_geneid:
        for g in genes:
            if g[4] == fgfr2_geneid:
                frow = g
                break
    if frow is None:
        for g in genes:
            if g[5].upper() == "FGFR2":
                frow = g
                break
    if frow is None:
        return out
    seqid, fstart, fend, fstrand, fgid, fsym = frow
    same = sorted([g for g in genes if g[0] == seqid], key=lambda g: (g[1], g[2]))
    idx = next((i for i, g in enumerate(same) if g[4] == fgid), None)
    if idx is None:
        return out
    out["fgfr2"] = {"seqid": seqid, "start": fstart, "end": fend, "strand": fstrand,
                    "gene_id": fgid, "symbol": fsym, "protein_id": pid_by_gene.get(fgid, "")}
    # genomic left (lower coord) and right (higher coord) protein-coding neighbors
    left = list(reversed(same[:idx]))          # nearest first
    right = same[idx + 1:]                      # nearest first
    plus = not is_reverse(fstrand)
    up_side = left if plus else right           # upstream = 5' of FGFR2
    down_side = right if plus else left
    neigh = []
    for side_name, lst in (("upstream", up_side), ("downstream", down_side)):
        for rank, g in enumerate(lst[:SUPP_N], start=1):
            _gseq, gs, ge, gstr, ggid, gsym = g
            dist = 0 if not (ge < fstart or gs > fend) else min(abs(fstart - ge), abs(gs - fend))
            neigh.append({"rank": rank, "side": side_name, "gene_id": ggid, "symbol": gsym,
                          "start": gs, "end": ge, "strand": gstr, "biotype": "protein_coding",
                          "protein_id": pid_by_gene.get(ggid, ""), "distance": dist})
    out["neighbors"] = neigh
    n_up = sum(1 for n in neigh if n["side"] == "upstream")
    n_down = sum(1 for n in neigh if n["side"] == "downstream")
    if n_up == 0 and n_down == 0:
        out["status"] = "fgfr2_locus_found_no_neighbors"
    elif n_up < MAIN_N or n_down < MAIN_N:
        out["status"] = "partial_neighborhood_scaffold_edge"
    else:
        out["status"] = "neighborhood_extracted"
    return out


# ---------------------------------------------------------------------------
# per-species cache (GFF parse is expensive; cache parsed neighborhood + proteins)
# ---------------------------------------------------------------------------
def cache_paths(dirs, species: str) -> Tuple[Path, Path]:
    cdir = dirs["synteny"] / "_cache"
    cdir.mkdir(parents=True, exist_ok=True)
    return cdir / f"{species}_neighborhood.json", cdir / f"{species}_neighbor_proteins.faa"


def extract_proteins(faa: Path, want: set) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not faa or not want:
        return out
    cur, seq = None, []
    with open(faa, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None and cur in want:
                    out[cur] = "".join(seq)
                cur = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if cur is not None and cur in want:
        out[cur] = "".join(seq)
    return out


def get_neighborhood(base: Path, dirs, species: str, taxid: str, acc: str,
                     refresh: bool) -> Dict[str, object]:
    cj, cf = cache_paths(dirs, species)
    if cj.exists() and not refresh:
        data = json.loads(cj.read_text(encoding="utf-8"))
        return data
    gff = gff_path(base, taxid, acc)
    if gff is None:
        data = {"status": "annotation_unavailable", "neighbors": [], "fgfr2": None,
                "assembly_accession": acc, "taxid": taxid, "source_file": ""}
        cj.write_text(json.dumps(data), encoding="utf-8")
        M.write_fasta(cf, [])
        return data
    data = parse_gff_for_fgfr2(gff)
    data["assembly_accession"] = acc
    data["taxid"] = taxid
    data["source_file"] = str(gff)
    # extract neighbor + fgfr2 protein sequences for downstream BLAST/RBH
    want = {n["protein_id"] for n in data.get("neighbors", []) if n.get("protein_id")}
    if data.get("fgfr2", {}) and data["fgfr2"].get("protein_id"):
        want.add(data["fgfr2"]["protein_id"])
    seqs = extract_proteins(protein_faa_path(base, taxid, acc), want)
    for n in data.get("neighbors", []):
        n["protein_sequence_available"] = "true" if seqs.get(n.get("protein_id")) else "false"
    cj.write_text(json.dumps(data), encoding="utf-8")
    M.write_fasta(cf, sorted(seqs.items()))
    return data


def load_neighbor_proteins(dirs, species: str) -> Dict[str, str]:
    _, cf = cache_paths(dirs, species)
    return {sid: seq for sid, seq in M.read_fasta(cf)}


# ---------------------------------------------------------------------------
# DIAMOND / BLASTP search helpers
# ---------------------------------------------------------------------------
OUTFMT = "6 qseqid sseqid pident length qlen slen evalue bitscore"


def prebuild_db(db_faa: Path, diamond: Optional[str], makeblastdb: Optional[str],
                out_dir: Path) -> Dict[str, str]:
    """Build a reusable search DB once (persistent). Returns {'engine', 'db': stem} or {}."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if diamond:
        stem = out_dir / (db_faa.stem + "_dmnd")
        if not stem.with_suffix(".dmnd").exists():
            subprocess.run([diamond, "makedb", "--in", str(db_faa), "-d", str(stem)],
                           check=True, capture_output=True, text=True)
        return {"engine": "diamond", "db": str(stem)}
    if makeblastdb:
        stem = out_dir / (db_faa.stem + "_blast")
        if not (stem.with_suffix(".pin").exists() or stem.with_suffix(".pdb").exists()):
            subprocess.run([makeblastdb, "-in", str(db_faa), "-dbtype", "prot", "-out", str(stem)],
                           check=True, capture_output=True, text=True)
        return {"engine": "blastp", "db": str(stem)}
    return {}


def search_vs_prebuilt(query_faa: Path, prebuilt: Dict[str, str], out_tsv: Path,
                       diamond: Optional[str], blastp: Optional[str], evalue: str = "1e-3",
                       sensitivity: str = "--more-sensitive") -> str:
    """Search a query against an already-built DB (no per-call makedb)."""
    if prebuilt.get("engine") == "diamond" and diamond:
        subprocess.run([diamond, "blastp", "-q", str(query_faa), "-d", prebuilt["db"], "-o",
                        str(out_tsv), "--outfmt"] + OUTFMT.split() +
                       ["--evalue", evalue, "--max-target-seqs", "25", "--quiet", sensitivity],
                       check=True, capture_output=True, text=True)
        return f"diamond {_ver(diamond)}"
    if prebuilt.get("engine") == "blastp" and blastp:
        subprocess.run([blastp, "-query", str(query_faa), "-db", prebuilt["db"], "-outfmt", OUTFMT,
                        "-evalue", evalue, "-max_target_seqs", "25", "-out", str(out_tsv)],
                       check=True, capture_output=True, text=True)
        return f"blastp {_ver(blastp, '-version')}"
    return ""


def run_search(query_faa: Path, db_faa: Path, out_tsv: Path, diamond: Optional[str],
               blastp: Optional[str], makeblastdb: Optional[str], workdir: Path) -> str:
    """Run protein search query->db (all hits), building the DB in workdir. Returns engine or ''."""
    try:
        pre = prebuild_db(db_faa, diamond, makeblastdb, workdir)
        if pre:
            return search_vs_prebuilt(query_faa, pre, out_tsv, diamond, blastp)
    except Exception:
        pass
    return ""


def _ver(tool: str, flag: str = "version") -> str:
    try:
        out = subprocess.run([tool, flag], capture_output=True, text=True).stdout
        return out.strip().split("\n")[0]
    except Exception:
        return "?"


def parse_hits(out_tsv: Path) -> List[Dict[str, object]]:
    rows = []
    if not out_tsv.exists():
        return rows
    for ln in out_tsv.read_text(encoding="utf-8").splitlines():
        p = ln.split("\t")
        if len(p) < 8:
            continue
        q, s, pid, length, qlen, slen, ev, bs = p[:8]
        ln_, ql, sl = M.to_float(length, 0), M.to_float(qlen, 1) or 1, M.to_float(slen, 1) or 1
        rows.append({"q": q, "s": s, "pident": M.to_float(pid, 0.0),
                     "evalue": M.to_float(ev, 1.0), "bitscore": M.to_float(bs, 0.0),
                     "qcov": round(ln_ / ql, 3), "scov": round(ln_ / sl, 3)})
    return rows


def best_per_query(hits: List[Dict[str, object]]) -> Dict[str, List[Dict[str, object]]]:
    by_q: Dict[str, List[Dict[str, object]]] = {}
    for h in hits:
        by_q.setdefault(h["q"], []).append(h)
    for q in by_q:
        by_q[q].sort(key=lambda h: (-h["bitscore"], h["evalue"]))
    return by_q


def is_generic_symbol(sym: str) -> bool:
    s = (sym or "").upper()
    if not s or s.startswith("GENEID:"):
        return True
    return any(s.startswith(p) for p in GENERIC_PREFIXES) or "UNCHARACTERIZED" in s


# ---------------------------------------------------------------------------
# Part A — table rows
# ---------------------------------------------------------------------------
A_COLS = ["species", "assembly_accession", "annotation_source", "annotation_release_or_date",
          "seqid", "fgfr2_gene_id", "fgfr2_transcript_id", "fgfr2_protein_id", "fgfr2_start",
          "fgfr2_end", "fgfr2_strand", "neighbor_rank", "neighbor_side", "neighbor_gene_id",
          "neighbor_symbol_raw", "neighbor_start", "neighbor_end", "neighbor_strand",
          "neighbor_biotype", "neighbor_protein_id", "neighbor_protein_sequence_available",
          "distance_to_fgfr2", "source_file", "extraction_status", "extraction_warning"]


def neighborhood_rows(species: str, data: Dict[str, object], asm: Dict[str, str],
                      final_tx: str, max_n: int) -> List[Dict[str, object]]:
    f = data.get("fgfr2") or {}
    status = data.get("status", "annotation_unavailable")
    warn = ""
    if status == "annotation_unavailable":
        warn = "no source-compatible RefSeq genomic.gff in local NCBI Datasets cache"
    elif status == "fgfr2_locus_not_found":
        warn = "FGFR2 gene feature not found in source annotation"
    elif status == "partial_neighborhood_scaffold_edge":
        warn = "fewer than 5 protein-coding neighbors on one side (scaffold/assembly edge)"
    base = {
        "species": species, "assembly_accession": data.get("assembly_accession", ""),
        "annotation_source": "NCBI RefSeq (datasets genomic.gff)",
        "annotation_release_or_date": asm.get("assembly_name", ""),
        "seqid": f.get("seqid", ""), "fgfr2_gene_id": f.get("gene_id", ""),
        "fgfr2_transcript_id": final_tx, "fgfr2_protein_id": f.get("protein_id", ""),
        "fgfr2_start": f.get("start", ""), "fgfr2_end": f.get("end", ""),
        "fgfr2_strand": f.get("strand", ""), "source_file": data.get("source_file", ""),
        "extraction_status": status, "extraction_warning": warn,
    }
    rows = []
    neigh = [n for n in data.get("neighbors", []) if n["rank"] <= max_n]
    if not neigh:
        r = dict(base)
        r.update({"neighbor_rank": "", "neighbor_side": "", "neighbor_gene_id": "",
                  "neighbor_symbol_raw": "", "neighbor_start": "", "neighbor_end": "",
                  "neighbor_strand": "", "neighbor_biotype": "", "neighbor_protein_id": "",
                  "neighbor_protein_sequence_available": "", "distance_to_fgfr2": ""})
        rows.append(r)
        return rows
    for n in sorted(neigh, key=lambda x: (x["side"], x["rank"])):
        r = dict(base)
        r.update({"neighbor_rank": n["rank"], "neighbor_side": n["side"],
                  "neighbor_gene_id": n["gene_id"], "neighbor_symbol_raw": n["symbol"],
                  "neighbor_start": n["start"], "neighbor_end": n["end"],
                  "neighbor_strand": n["strand"], "neighbor_biotype": n["biotype"],
                  "neighbor_protein_id": n["protein_id"],
                  "neighbor_protein_sequence_available": n.get("protein_sequence_available", ""),
                  "distance_to_fgfr2": n["distance"]})
        rows.append(r)
    return rows


ID_COLS = ["species", "neighbor_gene_id", "neighbor_symbol_raw", "neighbor_protein_id",
           "neighbor_rank", "neighbor_side", "human_reference_neighbor_gene_symbol",
           "human_reference_neighbor_protein_id", "normalized_neighbor_symbol",
           "normalized_neighbor_orthology_group", "identity_resolution_method",
           "identity_resolution_confidence", "best_hit_human_neighbor", "best_hit_evalue",
           "best_hit_bitscore", "best_hit_percent_identity", "best_hit_query_coverage",
           "best_hit_subject_coverage", "reciprocal_best_hit", "competing_hits",
           "identity_resolution_status", "identity_resolution_warning",
           "broad_homology_symbol", "broad_homology_protein_id", "broad_homology_percent_identity",
           "broad_homology_query_coverage", "broad_homology_subject_coverage",
           "broad_homology_evalue", "broad_homology_method"]
BROAD_DEFAULTS = {"broad_homology_symbol": "", "broad_homology_protein_id": "",
                  "broad_homology_percent_identity": "", "broad_homology_query_coverage": "",
                  "broad_homology_subject_coverage": "", "broad_homology_evalue": "",
                  "broad_homology_method": ""}
BLAST_COLS = ["species", "direction", "neighbor_gene_id", "query_protein_id",
              "human_neighbor_symbol", "human_protein_id", "percent_identity", "evalue",
              "bitscore", "query_coverage", "subject_coverage"]
RBH_COLS = ["species", "neighbor_gene_id", "neighbor_protein_id", "human_neighbor_symbol",
            "human_protein_id", "forward_bitscore", "reverse_bitscore", "reciprocal_best_hit"]
VALID_COLS = ["species", "validation_group", "fgfr2_locus_status", "n_5neighbor_slots_available",
              "n_neighbors_symbol_supported", "n_neighbors_curated_orthology_supported",
              "n_neighbors_rbh_supported", "n_neighbors_one_way_blast_supported",
              "n_neighbors_ambiguous", "n_neighbors_unmapped", "left_neighbor_support_score",
              "right_neighbor_support_score", "total_neighbor_support_score", "synteny_order_score",
              "neighbor_identity_confidence_score", "scaffold_continuity_status",
              "rescued_candidate_locus_support", "synteny_validation_class", "synteny_warning"]


def _passes(h: Dict[str, object], thr: Dict[str, float]) -> bool:
    return (h["evalue"] <= thr["evalue"] and h["pident"] >= thr["pid"]
            and h["qcov"] >= thr["qcov"] and h["scov"] >= thr["scov"])


def run_identity(base, dirs, nbh, master, human_ref5, human_ref10, h10_fa, diamond, blastp,
                 makeblastdb, refresh=False):
    """Parts C + D: resolve every neighbor's identity via symbol > RBH > one-way BLAST > raw,
    using DIAMOND/BLASTP against the human FGFR2 10-neighbor reference."""
    syn = dirs["synteny"]
    human_db = syn / "human_fgfr2_10neighbor_reference_proteins.faa"
    hpid2sym = {n["protein_id"]: sym for sym, n in human_ref10.items() if n.get("protein_id")}

    def hdr_sym(header: str) -> Tuple[str, str]:
        # human db headers are "SYMBOL|PROTEIN_ID"
        sym, _, pid = header.partition("|")
        return sym, pid

    identity_rows, blast_rows, rbh_rows, search_manifest = [], [], [], []
    have_engine = bool(diamond or (blastp and makeblastdb))
    for sp, data in nbh.items():
        neigh = [n for n in data.get("neighbors", []) if n["rank"] <= SUPP_N]
        group = master.get(sp, {}).get("validation_group", "standard_species")
        thr = THR_STRICT if group in CLOSE_GROUPS else THR_STD
        sp_seqs = load_neighbor_proteins(dirs, sp)
        fwd_best, rev_best, fwd_all = {}, {}, {}
        engine = ""
        idcache = dirs["synteny"] / "_cache" / f"{sp}_idhits.json"
        if have_engine and h10_fa and sp != HUMAN:
            if idcache.exists() and not refresh:
                d = json.loads(idcache.read_text(encoding="utf-8"))
                engine, fhits, rhits = d.get("engine", ""), d.get("fwd", []), d.get("rev", [])
            else:
                with tempfile.TemporaryDirectory() as td:
                    tdp = Path(td)
                    qf = tdp / "sp.faa"
                    M.write_fasta(qf, [(n["protein_id"], sp_seqs[n["protein_id"]]) for n in neigh
                                       if sp_seqs.get(n["protein_id"])])
                    fo = tdp / "fwd.tsv"
                    engine = run_search(qf, human_db, fo, diamond, blastp, makeblastdb, tdp)
                    fhits = parse_hits(fo)
                    ro = tdp / "rev.tsv"
                    run_search(human_db, qf, ro, diamond, blastp, makeblastdb, tdp)
                    rhits = parse_hits(ro)
                idcache.write_text(json.dumps({"engine": engine, "fwd": fhits, "rev": rhits}),
                                   encoding="utf-8")
            fwd_all = best_per_query(fhits)
            fwd_best = {q: hs[0] for q, hs in fwd_all.items()}
            rev_best = {q: hs[0] for q, hs in best_per_query(rhits).items()}
        # search manifest row
        search_manifest.append({
            "species": sp, "engine": engine or ("none" if sp != HUMAN else "skipped_human_self"),
            "human_reference_db": human_db.name, "n_query_proteins": sum(
                1 for n in neigh if sp_seqs.get(n["protein_id"])),
            "evalue_threshold": thr["evalue"], "min_percent_identity": thr["pid"],
            "min_query_coverage": thr["qcov"], "min_subject_coverage": thr["scov"],
            "threshold_set": "strict_close_mammal_primate" if group in CLOSE_GROUPS else "standard",
            "n_forward_hits": sum(len(v) for v in fwd_all.values())})
        for n in neigh:
            sym, pid = n["symbol"], n["protein_id"]
            row = {"species": sp, "neighbor_gene_id": n["gene_id"], "neighbor_symbol_raw": sym,
                   "neighbor_protein_id": pid, "neighbor_rank": n["rank"],
                   "neighbor_side": n["side"], "reciprocal_best_hit": "false",
                   "competing_hits": "", "identity_resolution_warning": "", **dict(BROAD_DEFAULTS)}
            # forward best hit metrics
            fb = fwd_best.get(pid)
            hsym, hpid = (hdr_sym(fb["s"]) if fb else ("", ""))
            if fb:
                row.update({"best_hit_human_neighbor": hsym, "best_hit_evalue": fb["evalue"],
                            "best_hit_bitscore": fb["bitscore"],
                            "best_hit_percent_identity": fb["pident"],
                            "best_hit_query_coverage": fb["qcov"],
                            "best_hit_subject_coverage": fb["scov"]})
                for h in fwd_all.get(pid, [])[:8]:
                    s2, p2 = hdr_sym(h["s"])
                    blast_rows.append({"species": sp, "direction": "forward",
                                       "neighbor_gene_id": n["gene_id"], "query_protein_id": pid,
                                       "human_neighbor_symbol": s2, "human_protein_id": p2,
                                       "percent_identity": h["pident"], "evalue": h["evalue"],
                                       "bitscore": h["bitscore"], "query_coverage": h["qcov"],
                                       "subject_coverage": h["scov"]})
                comp = [hdr_sym(h["s"])[0] for h in fwd_all.get(pid, [])[1:]
                        if h["bitscore"] >= 0.9 * fb["bitscore"]]
                row["competing_hits"] = ",".join(dict.fromkeys(comp))
            else:
                row.update({"best_hit_human_neighbor": "", "best_hit_evalue": "",
                            "best_hit_bitscore": "", "best_hit_percent_identity": "",
                            "best_hit_query_coverage": "", "best_hit_subject_coverage": ""})
            # reciprocal best hit
            rbh = False
            if fb:
                back = rev_best.get(fb["s"])
                rbh = bool(back and back["s"] == pid)
                if back:
                    rbh_rows.append({"species": sp, "neighbor_gene_id": n["gene_id"],
                                     "neighbor_protein_id": pid, "human_neighbor_symbol": hsym,
                                     "human_protein_id": hpid, "forward_bitscore": fb["bitscore"],
                                     "reverse_bitscore": back["bitscore"],
                                     "reciprocal_best_hit": "true" if rbh else "false"})
                row["reciprocal_best_hit"] = "true" if rbh else "false"
            # resolution hierarchy
            method, status, conf, href_sym, href_pid, norm, ogroup = _resolve(
                sym, sp, human_ref5, human_ref10, hpid2sym, fb, rbh, thr, row["competing_hits"])
            row.update({"human_reference_neighbor_gene_symbol": href_sym,
                        "human_reference_neighbor_protein_id": href_pid,
                        "normalized_neighbor_symbol": norm,
                        "normalized_neighbor_orthology_group": ogroup,
                        "identity_resolution_method": method,
                        "identity_resolution_confidence": conf,
                        "identity_resolution_status": status})
            if status == "unmapped_neighbor" and fb and _passes(fb, THR_STD):
                row["identity_resolution_warning"] = "passing BLAST hit present but not adopted"
            identity_rows.append(row)
    return identity_rows, blast_rows, rbh_rows, search_manifest


def human_pid2sym(gff: Path) -> Dict[str, str]:
    """Map every human RefSeq protein_id -> its gene symbol (first occurrence) from CDS lines."""
    out: Dict[str, str] = {}
    with open(gff, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "\tCDS\t" not in line or "protein_id=" not in line:
                continue
            pi = line.find("protein_id=")
            ps = pi + 11
            pe = ps
            while pe < len(line) and line[pe] not in ";\n\r\t":
                pe += 1
            pid = line[ps:pe]
            if pid in out:
                continue
            gi = line.find(";gene=")
            sym = ""
            if gi >= 0:
                gs = gi + 6
                ge = gs
                while ge < len(line) and line[ge] not in ";\n\r\t":
                    ge += 1
                sym = line[gs:ge]
            out[pid] = sym or pid
    return out


def build_human_proteome_db(base: Path, dirs, refresh: bool) -> Optional[Path]:
    """Whole human RefSeq proteome FASTA with 'SYMBOL|protein_id' headers (cached) — the broad
    reference for naming uncharacterized (LOC...) / unresolved neighbors by best human homolog."""
    cache = dirs["synteny"] / "_cache" / "human_proteome_named.faa"
    if cache.exists() and cache.stat().st_size > 0 and not refresh:
        return cache
    asm = {(r.get("species_canonical") or "").lower(): r for r in
           M.read_tsv(base / "02_models" / "ncbi_assembly_selected.tsv")}
    h = asm.get(HUMAN, {})
    gff = gff_path(base, h.get("taxid", ""), h.get("assembly_accession", ""))
    faa = protein_faa_path(base, h.get("taxid", ""), h.get("assembly_accession", ""))
    if not gff or not faa:
        # Custom runs typically have no run-local human genome (human is not an analysed
        # species). Reuse the shared curated HUMAN reference proteome so LOC/uncharacterized
        # neighbors resolve with the SAME broad-homology logic as the Example dataset.
        if SHARED_HUMAN_PROTEOME.exists() and SHARED_HUMAN_PROTEOME.stat().st_size > 0:
            cache.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(SHARED_HUMAN_PROTEOME, cache)
            return cache
        return None
    pid2sym = human_pid2sym(gff)
    seqs, cur, buf = [], None, []
    with open(faa, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append((f"{pid2sym.get(cur, cur)}|{cur}", "".join(buf)))
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur:
        seqs.append((f"{pid2sym.get(cur, cur)}|{cur}", "".join(buf)))
    M.write_fasta(cache, seqs)
    return cache


def run_broad_homology(base, dirs, nbh, identity_rows, diamond, blastp, makeblastdb, refresh):
    """For every generic (LOC.../uncharacterized) or still-unresolved neighbor, find the best human
    proteome homolog (even loose) and record symbol + percent identity + query/subject coverage.
    Adopted only as a PROBABLE/loose name; it is never treated as an FGFR2-neighborhood ortholog
    (orthology_group stays empty) so local synteny scoring is unaffected."""
    have_engine = bool(diamond or (blastp and makeblastdb))
    if not have_engine:
        return []
    db = build_human_proteome_db(base, dirs, refresh)
    if not db:
        return []
    # build the proteome search DB ONCE and reuse it for every species (avoids per-call makedb)
    prebuilt = prebuild_db(db, diamond, makeblastdb, dirs["synteny"] / "_cache" / "_proteome_db")
    if not prebuilt:
        return []
    idx = {(r["species"], r["neighbor_gene_id"]): r for r in identity_rows}
    extra = []

    def assign(r, b, weak):
        bsym, _, bpid = b["s"].partition("|")
        if is_generic_symbol(bsym):
            return
        conf = ("very_low" if weak else
                "high" if (b["pident"] >= 60 and b["qcov"] >= 0.6)
                else "medium" if (b["pident"] >= 30 and b["qcov"] >= 0.4) else "low")
        r.update({"broad_homology_symbol": bsym, "broad_homology_protein_id": bpid,
                  "broad_homology_percent_identity": b["pident"],
                  "broad_homology_query_coverage": b["qcov"],
                  "broad_homology_subject_coverage": b["scov"],
                  "broad_homology_evalue": b["evalue"],
                  "broad_homology_method": ("diamond_weak_best_hit_vs_human_proteome" if weak
                                            else "diamond_best_hit_vs_human_proteome")})
        no_name = (is_generic_symbol(r.get("normalized_neighbor_symbol", "")) or
                   r["identity_resolution_status"] in ("raw_id_only", "unmapped_neighbor"))
        if no_name:
            r["normalized_neighbor_symbol"] = bsym
            r["identity_resolution_method"] = ("broad_proteome_weak_best_hit" if weak
                                               else "broad_proteome_best_hit")
            r["identity_resolution_status"] = "broad_homology_named"
            r["identity_resolution_confidence"] = conf
            r["identity_resolution_warning"] = (
                f"{'very weak' if weak else 'loose'} human homology name "
                f"(id {b['pident']:.0f}%, qcov {int(b['qcov']*100)}%, scov {int(b['scov']*100)}%)")

    def search(sp, need, suffix, evalue, sensitivity):
        cache = dirs["synteny"] / "_cache" / f"{sp}_broadhits{suffix}.json"
        if cache.exists() and not refresh:
            return best_per_query(json.loads(cache.read_text(encoding="utf-8")))
        sp_seqs = load_neighbor_proteins(dirs, sp)
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            qf = tdp / "q.faa"
            M.write_fasta(qf, [(n["protein_id"], sp_seqs[n["protein_id"]]) for n in need
                               if sp_seqs.get(n["protein_id"])])
            oo = tdp / "broad.tsv"
            search_vs_prebuilt(qf, prebuilt, oo, diamond, blastp, evalue, sensitivity)
            hits = parse_hits(oo)
        cache.write_text(json.dumps(hits), encoding="utf-8")
        return best_per_query(hits)

    for sp, data in nbh.items():
        sp_seqs = load_neighbor_proteins(dirs, sp)
        need = []
        for n in data.get("neighbors", []):
            if n["rank"] > SUPP_N or not sp_seqs.get(n["protein_id"]):
                continue
            r = idx.get((sp, n["gene_id"]))
            if not r:
                continue
            if (is_generic_symbol(n["symbol"]) or r["identity_resolution_status"] in
                    ("raw_id_only", "unmapped_neighbor", "ambiguous_neighbor_identity")):
                need.append(n)
        if not need:
            continue
        # pass 1: standard sensitivity, E<=1e-3
        best = search(sp, need, "", "1e-3", "--more-sensitive")
        for n in need:
            hs = best.get(n["protein_id"])
            if hs:
                assign(idx[(sp, n["gene_id"])], hs[0], weak=False)
                for h in hs[:5]:
                    s2, _, p2 = h["s"].partition("|")
                    extra.append({"species": sp, "direction": "broad_vs_human_proteome",
                                  "neighbor_gene_id": n["gene_id"], "query_protein_id": n["protein_id"],
                                  "human_neighbor_symbol": s2, "human_protein_id": p2,
                                  "percent_identity": h["pident"], "evalue": h["evalue"],
                                  "bitscore": h["bitscore"], "query_coverage": h["qcov"],
                                  "subject_coverage": h["scov"]})
        # pass 2: very permissive fallback for neighbors that still have no name
        left = [n for n in need if not idx[(sp, n["gene_id"])].get("broad_homology_symbol")]
        if left:
            wbest = search(sp, left, "_weak", "1000", "--very-sensitive")
            for n in left:
                hs = wbest.get(n["protein_id"])
                if hs:
                    assign(idx[(sp, n["gene_id"])], hs[0], weak=True)
                    h = hs[0]
                    s2, _, p2 = h["s"].partition("|")
                    extra.append({"species": sp, "direction": "broad_weak_vs_human_proteome",
                                  "neighbor_gene_id": n["gene_id"], "query_protein_id": n["protein_id"],
                                  "human_neighbor_symbol": s2, "human_protein_id": p2,
                                  "percent_identity": h["pident"], "evalue": h["evalue"],
                                  "bitscore": h["bitscore"], "query_coverage": h["qcov"],
                                  "subject_coverage": h["scov"]})
    return extra


def _resolve(sym, sp, human_ref5, human_ref10, hpid2sym, fb, rbh, thr, competing):
    """Return (method, status, confidence, href_sym, href_pid, normalized, orthology_group)."""
    # 1. exact / high-confidence gene symbol match to human FGFR2 neighborhood (case-insensitive)
    if not is_generic_symbol(sym) and sym.upper() in human_ref10:
        hn = human_ref10[sym.upper()]
        canon = hn.get("symbol", sym)
        return ("exact_symbol_match", "ortholog_supported", "high", canon,
                hn.get("protein_id", ""), canon, canon.upper())
    # BLAST-based tiers (only if a forward hit exists)
    if fb is not None:
        hsym, _, hpid = fb["s"].partition("|")
        ambiguous = bool(competing)
        if ambiguous:
            return ("unresolved", "ambiguous_neighbor_identity", "low", hsym, hpid,
                    sym if not is_generic_symbol(sym) else fb["s"].split("|")[0], "")
        if rbh and _passes(fb, thr):
            return ("reciprocal_best_hit", "probable_ortholog_supported", "high", hsym, hpid,
                    hsym, hsym)
        if _passes(fb, thr):
            return ("high_confidence_one_way_blast", "probable_ortholog_supported", "medium",
                    hsym, hpid, hsym, hsym)
        if _passes(fb, THR_STD):
            return ("high_confidence_one_way_blast", "probable_ortholog_supported", "low",
                    hsym, hpid, hsym, hsym)
    # symbol present but not a human FGFR2-neighborhood gene (a real, different gene)
    if not is_generic_symbol(sym):
        return ("exact_symbol_match", "symbol_supported_only", "medium", "", "", sym, "")
    # only a raw / generic identifier and no usable hit
    if fb is None:
        return ("raw_annotation_only", "raw_id_only", "low", "", "", sym, "")
    return ("unresolved", "unmapped_neighbor", "low", "", "", sym, "")


def score_synteny(nbh, master, truth_by_sp, human_ref5, human_ref10, id_by):
    cols_order = sorted(human_ref10.values(), key=lambda n: (n["side"], n["rank"]))
    col_syms = [n["symbol"] for n in cols_order]
    matrix_rows, valid_rows = [], []
    for sp, data in nbh.items():
        mr = master.get(sp, {})
        status = data.get("status", "annotation_unavailable")
        neigh = data.get("neighbors", [])
        by_sym = {}
        for n in neigh:
            by_sym.setdefault(n["symbol"].upper(), n)
        # identity-based orthology group per species (orthology group is the canonical upper symbol)
        og = {}
        for n in neigh:
            r = id_by.get((sp, n["gene_id"]), {})
            grp = (r.get("normalized_neighbor_orthology_group", "") or "").upper()
            if grp:
                og.setdefault(grp, (n, r))
        cells = {}
        present_same_order = _present_side_ok = 0
        for hn in cols_order:
            hs, hside, hrank = hn["symbol"], hn["side"], hn["rank"]
            hu = hs.upper()
            cell = "missing_or_unmapped"
            if status in ("annotation_unavailable", "fgfr2_locus_not_found"):
                cell = "scaffold_unavailable"
            elif hu in by_sym:
                n = by_sym[hu]
                if n["side"] == hside and n["rank"] == hrank:
                    cell = "present_same_side_and_order"
                elif n["side"] == hside and n["rank"] <= MAIN_N and hrank <= MAIN_N:
                    cell = "present_same_side_reordered"
                elif n["side"] != hside:
                    cell = "present_opposite_side"
                elif n["rank"] <= MAIN_N:
                    cell = "present_elsewhere_in_5neighbor_window"
                else:
                    cell = "present_only_in_10neighbor_supplement"
                if hrank <= MAIN_N and hside == n["side"]:
                    _present_side_ok += 1
                    if n["rank"] == hrank:
                        present_same_order += 1
            elif hu in og:
                rr = og[hu][1]
                if rr.get("identity_resolution_status") == "ambiguous_neighbor_identity":
                    cell = "ambiguous_identity"
                else:
                    cell = "probable_by_blast_rbh"
                    if hrank <= MAIN_N:
                        _present_side_ok += 1
            cells[hs] = cell
        # per-species identity status counts
        sp_ids = [r for (s, _), r in id_by.items() if s == sp]
        c_sym = sum(1 for r in sp_ids if r["identity_resolution_status"]
                    in ("ortholog_supported", "symbol_supported_only"))
        c_curated = sum(1 for r in sp_ids if r["identity_resolution_method"]
                        == "curated_orthology_mapping")
        c_rbh = sum(1 for r in sp_ids if r["identity_resolution_method"] == "reciprocal_best_hit")
        c_one = sum(1 for r in sp_ids if r["identity_resolution_method"]
                    == "high_confidence_one_way_blast")
        c_amb = sum(1 for r in sp_ids if r["identity_resolution_status"]
                    == "ambiguous_neighbor_identity")
        c_unm = sum(1 for r in sp_ids if r["identity_resolution_status"]
                    in ("unmapped_neighbor", "raw_id_only"))
        n_slots = sum(1 for hn in cols_order if hn["rank"] <= MAIN_N)
        up_ref = [hn for hn in cols_order if hn["side"] == "upstream" and hn["rank"] <= MAIN_N]
        down_ref = [hn for hn in cols_order if hn["side"] == "downstream" and hn["rank"] <= MAIN_N]

        def side_score(refs):
            if not refs:
                return 0.0
            ok = sum(1 for hn in refs if cells.get(hn["symbol"]) in
                     ("present_same_side_and_order", "present_same_side_reordered",
                      "probable_by_blast_rbh"))
            return round(ok / len(refs), 3)
        left = side_score(up_ref)
        right = side_score(down_ref)
        total = round((left + right) / 2, 3)
        order = round(present_same_order / n_slots, 3) if n_slots else 0.0
        conf = round((c_sym + c_rbh + c_one) / max(len(sp_ids), 1), 3)
        # rescue locus support
        rescued = any((r.get("rescue_decision") or "").startswith("rescued")
                      for r in truth_by_sp.get(sp, []))
        if status in ("annotation_unavailable", "fgfr2_locus_not_found"):
            locus_support = "sequence_only_no_locus_coordinates" if rescued else "locus_not_found"
            scaffold = "unavailable"
        elif status == "partial_neighborhood_scaffold_edge":
            locus_support = "locus_supported_partial" if rescued else "locus_supported_partial"
            scaffold = "scaffold_edge"
        else:
            locus_support = "locus_supported" if rescued else "locus_present"
            scaffold = "continuous"
        group = mr.get("validation_group", "standard_species")
        cls, warn = _synteny_class(status, total, order, left, right, group, c_amb, c_unm)
        valid_rows.append({
            "species": sp, "validation_group": group, "fgfr2_locus_status": status,
            "n_5neighbor_slots_available": n_slots, "n_neighbors_symbol_supported": c_sym,
            "n_neighbors_curated_orthology_supported": c_curated,
            "n_neighbors_rbh_supported": c_rbh, "n_neighbors_one_way_blast_supported": c_one,
            "n_neighbors_ambiguous": c_amb, "n_neighbors_unmapped": c_unm,
            "left_neighbor_support_score": left, "right_neighbor_support_score": right,
            "total_neighbor_support_score": total, "synteny_order_score": order,
            "neighbor_identity_confidence_score": conf, "scaffold_continuity_status": scaffold,
            "rescued_candidate_locus_support": locus_support, "synteny_validation_class": cls,
            "synteny_warning": warn})
        claim = ";".join(sorted({r.get("final_claim_status_after_rescue", "")
                                 for r in truth_by_sp.get(sp, [])} - {""}))
        dec = ";".join(sorted({r.get("rescue_decision", "")
                               for r in truth_by_sp.get(sp, [])} - {""}))
        mrow = {"species": sp, "taxon_group": mr.get("taxon_group_display", mr.get("taxon_group", "")),
                "validation_group": group, "final_claim_status_after_rescue": claim,
                "rescue_decision": dec, "synteny_validation_class": cls}
        for hs in col_syms:
            mrow[hs] = cells.get(hs, "missing_or_unmapped")
        matrix_rows.append(mrow)
    return matrix_rows, valid_rows


def _synteny_class(status, total, order, left, right, group, n_amb, n_unm):
    if status in ("annotation_unavailable",):
        return "synteny_unavailable", "no source-compatible annotation available"
    if status == "fgfr2_locus_not_found":
        return "synteny_sequence_only_support", "FGFR2 locus not found in source annotation"
    close = group in CLOSE_GROUPS
    if total >= 0.8 and order >= 0.6:
        return "synteny_strong", ""
    if total >= 0.6:
        return "synteny_supported_with_minor_rearrangement", \
            ("close group with rearrangement" if close else "")
    if status == "partial_neighborhood_scaffold_edge" and total >= 0.3:
        return "synteny_partial_scaffold_limit", "fewer neighbors available (scaffold/assembly edge)"
    if total >= 0.4:
        return "synteny_partial_blast_supported", "partial neighbor support (incl. BLAST/RBH)"
    if total < 0.3 and close:
        return "synteny_conflict_review", \
            "close group expected strong local synteny but neighbor support is low"
    if total < 0.3:
        return "synteny_partial_blast_supported", \
            "low neighbor order conservation (distant taxon; not over-penalized)"
    return "synteny_partial_blast_supported", ""


SYN_INTEGRATE_COLS = ["synteny_validation_class", "combined_synteny_validation_class",
                      "total_neighbor_support_score", "neighbor_identity_confidence_score",
                      "n_neighbors_rbh_supported", "n_neighbors_one_way_blast_supported",
                      "n_neighbors_ambiguous", "n_neighbors_unmapped", "scaffold_continuity_status",
                      "rescued_candidate_locus_support", "synteny_warning"]


def integrate_synteny(base, dirs, valid_rows) -> None:
    """Part G: add per-species synteny evidence columns to the master evidence tables. No MCScanX
    columns are added (the optional MCScanX block-level layer is intentionally omitted). Synteny
    never assigns/relabels IIIb/IIIc; it only annotates locus/orthology context."""
    by_sp = {}
    for r in valid_rows:
        payload = {c: r.get(c, "") for c in SYN_INTEGRATE_COLS if c in r}
        # no MCScanX in this build: combined synteny == local synteny
        payload["combined_synteny_validation_class"] = r.get("synteny_validation_class", "")
        by_sp[(r["species"] or "").lower()] = payload

    def patch(path):
        if not path or not Path(path).exists():
            return
        rows = M.read_tsv(path)
        if not rows:
            return
        fields = list(rows[0].keys()) + [c for c in SYN_INTEGRATE_COLS if c not in rows[0]]
        for row in rows:
            p = by_sp.get((row.get("species") or "").lower(), {})
            for c in SYN_INTEGRATE_COLS:
                row[c] = p.get(c, row.get(c, ""))
        M.write_tsv(path, rows, fields)

    patch(M.locate(base, "species_qc_master.tsv", "11_pre_interpro_master"))
    patch(dirs["maps"] / "fgfr2_post_rescue_final_truth_table.tsv")
    patch(dirs["maps"] / "fgfr2_exon_type_label_reconciliation.tsv")
    patch(M.locate(base, "fgfr2_orthology_evidence.tsv"))

    # extend the corrected dataset manifest with the synteny artifacts
    man = dirs["maps"].parent / "final_corrected_pre_interpro_dataset_manifest.tsv"
    if man.exists():
        rows = M.read_tsv(man)
        have = {r.get("artifact") for r in rows}
        syn = dirs["synteny"]
        for art, p in (("synteny_local_5neighbor_table",
                        syn / "fgfr2_local_gene_neighborhood_5neighbors.tsv"),
                       ("synteny_5neighbor_validation",
                        syn / "fgfr2_5neighbor_synteny_validation.tsv"),
                       ("synteny_neighbor_identity_resolution",
                        syn / "fgfr2_neighbor_identity_resolution.tsv"),
                       ("synteny_5neighbor_validation_gate",
                        syn / "fgfr2_5neighbor_synteny_validation_gate.tsv")):
            if art in have:
                continue
            rows.append({"artifact": art, "path": str(p), "exists": "true" if p.exists() else "false",
                         "rows": str(len(M.read_tsv(p))) if p.exists() else "0",
                         "sha256": M.sha256_file(p) if p.exists() else "", "role": "synteny"})
        M.write_tsv(man, rows, list(rows[0].keys()) if rows else ["artifact"])


def write_gate(syn, human, h5_rows, identity_rows, rows5, valid_rows, have_engine,
               human_in_panel=True) -> bool:
    checks = []

    def add(check, scope, ok, detail=""):
        checks.append({"check": check, "scope": scope,
                       "status": "pass" if ok else "FAIL", "detail": detail})
    # The human synteny neighborhood is the reference layer for locus validation. When
    # homo_sapiens is NOT part of the run panel (custom run) it is legitimately absent, so
    # the two human-reference checks are not applicable (recorded as pass with an explicit
    # detail) rather than hard failures. If human IS in the panel these checks are enforced.
    if not human_in_panel:
        add("human_5neighbor_reference_exists", "human", True,
            "not_applicable: homo_sapiens not in run panel; human used as reference control only")
        add("human_fgfr2_neighborhood_extracted", "human", True,
            "not_applicable: homo_sapiens not in run panel (custom run)")
    else:
        h5_ok = bool(h5_rows) and any(r.get("human_protein_sequence_available") == "true"
                                      for r in h5_rows)
        add("human_5neighbor_reference_exists", "human", h5_ok, f"{len(h5_rows)} reference neighbors")
        hstat = human.get("status", "")
        add("human_fgfr2_neighborhood_extracted", "human",
            hstat == "neighborhood_extracted", f"status={hstat}")
    miss_status = [r for r in identity_rows if not r.get("identity_resolution_status")]
    add("plotted_labels_have_resolution_status", "identity", not miss_status,
        f"{len(miss_status)} rows missing status")
    bad_unres = [f"{r['species']}/{r['neighbor_gene_id']}" for r in identity_rows
                 if r.get("identity_resolution_status") == "unmapped_neighbor"
                 and "not adopted" in (r.get("identity_resolution_warning") or "")]
    add("no_unresolved_where_blast_could_resolve", "identity", not bad_unres,
        "; ".join(bad_unres[:6]) or "ok")
    inferred = [r for r in identity_rows if r["identity_resolution_method"]
                in ("reciprocal_best_hit", "high_confidence_one_way_blast")]
    add("blast_rbh_manifest_present_if_inferred", "identity",
        (not inferred) or have_engine, f"{len(inferred)} inferred labels")
    mis_curated = [f"{r['species']}/{r['neighbor_gene_id']}" for r in identity_rows
                   if r["identity_resolution_method"] in
                   ("reciprocal_best_hit", "high_confidence_one_way_blast")
                   and r["identity_resolution_status"] == "ortholog_supported"]
    add("no_blast_inferred_mislabeled_curated", "identity", not mis_curated,
        "; ".join(mis_curated[:6]) or "ok")
    forced = [f"{r['species']}/{r['neighbor_gene_id']}" for r in identity_rows
              if r["identity_resolution_status"] == "ambiguous_neighbor_identity"
              and r["identity_resolution_confidence"] == "high"]
    add("ambiguous_hits_not_forced", "identity", not forced, "; ".join(forced[:6]) or "ok")
    add("synteny_validates_locus_only_not_isoform", "global", True,
        "no IIIb/IIIc label is derived from synteny")
    post = bool([r for r in rows5 if r.get("fgfr2_transcript_id")])
    add("figure_tables_based_on_post_rescue_candidates", "global", post,
        "neighborhood rows carry final post-rescue transcript ids")
    silent_conf = [r["species"] for r in valid_rows
                   if r["synteny_validation_class"] == "synteny_conflict_review"
                   and not r["synteny_warning"]]
    add("synteny_conflict_not_silent", "validation", not silent_conf,
        "; ".join(silent_conf) or "ok")

    M.write_tsv(syn / "fgfr2_5neighbor_synteny_validation_gate.tsv", checks,
                ["check", "scope", "status", "detail"])
    hard = any(c["status"] != "pass" for c in checks)
    (syn / "fgfr2_5neighbor_synteny_validation_gate.json").write_text(
        json.dumps({"checks": checks, "hard_fail": hard, "timestamp": M.now_iso()}, indent=2),
        encoding="utf-8")
    return not hard


def main() -> int:
    ap = argparse.ArgumentParser(description="FGFR2 local synteny / gene-neighborhood validation.")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--refresh_cache", action="store_true",
                    help="force re-parse of genomic.gff (ignore per-species cache)")
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    syn = dirs["synteny"]

    diamond, blastp, makeblastdb = find_diamond(), find_blastp(), find_makeblastdb()
    engine_note = diamond and f"diamond:{diamond}" or (blastp and "blastp" or "none")

    asm_map = {(r.get("species_canonical") or "").lower(): r for r in
               M.read_tsv(base / "02_models" / "ncbi_assembly_selected.tsv")}
    master = {(r.get("species") or "").lower(): r for r in
              M.read_tsv(M.require(base, "species_qc_master.tsv", "11_pre_interpro_master"))}
    truth = M.read_tsv(syn.parent / "maps" / "fgfr2_post_rescue_final_truth_table.tsv")
    truth_by_sp: Dict[str, List[Dict[str, str]]] = {}
    for r in truth:
        truth_by_sp.setdefault((r.get("species") or "").lower(), []).append(r)

    species_list = sorted(master.keys(),
                          key=lambda s: M.to_int(master[s].get("phylo_order"), 999) or 999)

    # ---- Part A: extract neighborhoods (cached) ----
    nbh: Dict[str, Dict[str, object]] = {}
    rows5, rows10, manifest = [], [], []
    for sp in species_list:
        a = asm_map.get(sp, {})
        taxid, acc = a.get("taxid", ""), a.get("assembly_accession", "")
        data = get_neighborhood(base, dirs, sp, taxid, acc, args.refresh_cache)
        nbh[sp] = data
        final_tx = ";".join(sorted({r.get("transcript_id", "") for r in truth_by_sp.get(sp, [])} - {""}))
        rows5 += neighborhood_rows(sp, data, a, final_tx, MAIN_N)
        rows10 += neighborhood_rows(sp, data, a, final_tx, SUPP_N)
        f = data.get("fgfr2") or {}
        manifest.append({
            "species": sp, "assembly_accession": data.get("assembly_accession", ""),
            "taxid": taxid, "annotation_source": "NCBI RefSeq datasets genomic.gff",
            "annotation_release_or_date": a.get("assembly_name", ""),
            "seqid": f.get("seqid", ""), "fgfr2_found": "true" if f else "false",
            "n_upstream_neighbors": sum(1 for n in data.get("neighbors", [])
                                        if n["side"] == "upstream"),
            "n_downstream_neighbors": sum(1 for n in data.get("neighbors", [])
                                          if n["side"] == "downstream"),
            "extraction_status": data.get("status", ""),
            "source_file": data.get("source_file", "")})
    M.write_tsv(syn / "fgfr2_local_gene_neighborhood_5neighbors.tsv", rows5, A_COLS)
    M.write_tsv(syn / "fgfr2_local_gene_neighborhood_10neighbors_supplement.tsv", rows10, A_COLS)
    M.write_tsv(syn / "fgfr2_local_gene_neighborhood_manifest.tsv", manifest,
                list(manifest[0].keys()) if manifest else ["species"])

    # ---- Part B: human reference (5 + 10 neighbors) ----
    human = nbh.get(HUMAN, {})
    hproteins = load_neighbor_proteins(dirs, HUMAN)
    if not human.get("neighbors"):
        # No run-local human genome (human is not an analysed species in this custom run).
        # Fall back to the shared curated HUMAN FGFR2 neighborhood reference so neighbor
        # identity/orthology is classified with the SAME logic as the Example dataset.
        shared_human, shared_hprot = load_shared_human_reference()
        if shared_human.get("neighbors"):
            human = shared_human
            hproteins = shared_hprot
            print(f"[synteny] using shared HUMAN FGFR2 neighborhood reference "
                  f"({len(shared_human['neighbors'])} neighbors) — human_reference_control")
    ref_cols = ["human_neighbor_rank", "human_neighbor_side", "human_gene_id", "human_gene_symbol",
                "human_protein_id", "human_start", "human_end", "human_strand",
                "human_protein_sequence_available", "reference_warning"]

    def human_ref(max_n: int) -> Tuple[List[Dict[str, object]], List[Tuple[str, str]]]:
        rows, fasta = [], []
        for n in sorted(human.get("neighbors", []), key=lambda x: (x["side"], x["rank"])):
            if n["rank"] > max_n:
                continue
            pid = n.get("protein_id", "")
            seq = hproteins.get(pid, "")
            rows.append({"human_neighbor_rank": n["rank"], "human_neighbor_side": n["side"],
                         "human_gene_id": n["gene_id"], "human_gene_symbol": n["symbol"],
                         "human_protein_id": pid, "human_start": n["start"], "human_end": n["end"],
                         "human_strand": n["strand"],
                         "human_protein_sequence_available": "true" if seq else "false",
                         "reference_warning": "" if seq else "no reference protein sequence"})
            if seq:
                fasta.append((f"{n['symbol']}|{pid}", seq))
        return rows, fasta

    h5_rows, h5_fa = human_ref(MAIN_N)
    h10_rows, h10_fa = human_ref(SUPP_N)
    M.write_tsv(syn / "human_fgfr2_5neighbor_reference.tsv", h5_rows, ref_cols)
    M.write_tsv(syn / "human_fgfr2_10neighbor_reference.tsv", h10_rows, ref_cols)
    M.write_fasta(syn / "human_fgfr2_5neighbor_reference_proteins.faa", h5_fa)
    M.write_fasta(syn / "human_fgfr2_10neighbor_reference_proteins.faa", h10_fa)
    # human symbol dictionaries keyed by UPPERCASE symbol (cross-species case-insensitive match)
    human_ref10 = {n["symbol"].upper(): n for n in human.get("neighbors", []) if n["rank"] <= SUPP_N}
    human_ref5 = {s: n for s, n in human_ref10.items() if n["rank"] <= MAIN_N}

    # ---- Part C + D: neighbor identity resolution with BLAST/RBH ----
    identity_rows, blast_rows, rbh_rows, search_manifest = run_identity(
        base, dirs, nbh, master, human_ref5, human_ref10, h10_fa, diamond, blastp, makeblastdb,
        args.refresh_cache)
    # broad human-proteome homology naming for LOC.../unresolved neighbors (with %id/coverage)
    blast_rows += run_broad_homology(base, dirs, nbh, identity_rows, diamond, blastp, makeblastdb,
                                     args.refresh_cache)
    M.write_tsv(syn / "fgfr2_neighbor_identity_resolution.tsv", identity_rows, ID_COLS)
    M.write_tsv(syn / "neighbor_identity_blast_hits.tsv", blast_rows, BLAST_COLS)
    M.write_tsv(syn / "neighbor_identity_rbh_pairs.tsv", rbh_rows, RBH_COLS)
    M.write_tsv(syn / "neighbor_identity_search_manifest.tsv", search_manifest,
                list(search_manifest[0].keys()) if search_manifest else ["species"])

    # ---- Part E: conservation matrix + synteny validation ----
    id_by = {}
    for r in identity_rows:
        id_by[(r["species"], r["neighbor_gene_id"])] = r
    matrix_rows, valid_rows = score_synteny(nbh, master, truth_by_sp, human_ref5, human_ref10, id_by)
    M.write_tsv(syn / "fgfr2_5neighbor_conservation_matrix.tsv", matrix_rows,
                list(matrix_rows[0].keys()) if matrix_rows else ["species"])
    M.write_tsv(syn / "fgfr2_5neighbor_synteny_validation.tsv", valid_rows, VALID_COLS)

    # ---- Part G: integrate synteny evidence into master / truth / orthology / reconciliation ----
    integrate_synteny(base, dirs, valid_rows)

    # ---- Part I: local 5-neighbor synteny validation gate ----
    ok = write_gate(syn, human, h5_rows, identity_rows, rows5, valid_rows,
                    bool(diamond or (blastp and makeblastdb)),
                    human_in_panel=(HUMAN in species_list))

    print(f"[OK] synteny: species={len(species_list)} engine={engine_note} "
          f"extracted={sum(1 for s in species_list if nbh[s].get('status')=='neighborhood_extracted')} "
          f"gate={'PASS' if ok else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
