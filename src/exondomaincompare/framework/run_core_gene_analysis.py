#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import unquote

from exondomaincompare.config import discover_repository_root, load_config

PROJECT_ROOT = discover_repository_root(__file__)
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from scripts import create_new_run as cnr  # noqa: E402
from exondomaincompare.framework.gene_config import GeneConfig, GeneConfigError, load_gene_config, load_gene_config_lenient, build_generic_gene_config, GENE_DRAFT_DIR, GENE_CONFIG_DIR
from exondomaincompare.framework.build_core_gene_indices import CoreSource, build_core_indices  # noqa: E402
from exondomaincompare.framework import production_contract  # noqa: E402
from exondomaincompare.framework import run_labels  # noqa: E402
from exondomaincompare.framework.coordinate_evidence_register import (  # noqa: E402
    build_coordinate_evidence_register,
)
from exondomaincompare.contracts import (  # noqa: E402
    portable_path_reference, portable_runtime_record, stamp_payload,
)
from exondomaincompare.runs.legacy import LegacyRunAdapter  # noqa: E402
from exondomaincompare.runs.registry import resolve_run_record  # noqa: E402
from exondomaincompare.runs.layout import RunLayout, RunLayoutVersion  # noqa: E402
from scripts.make_fgfr2_post_interpro_exon_domain_figures import (  # noqa: E402
    _parse_topology_line as parse_pytmhmm_topology_line,
)
from exondomaincompare.framework.interpro_annotations import (  # noqa: E402
    load_normalized_annotations, representative_domains,
    family_annotations, feature_annotations,
)
from analyze_exon_domain_boundary_consistency import classify as classify_exon_domain_boundary  # noqa: E402
from collect_fgfr2_models_dual_source_v3 import build_cds_features_from_parts  # noqa: E402
from exondomaincompare.shared_gene_analysis import gene_locus_resolution as glr  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)
PROJECT_ROOT = RUNTIME_CONFIG.repository_root
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
RUNS_ROOT = RUNTIME_CONFIG.runs_root
NEAR_BOUNDARY_AA = 5
SYNTENY_PER_SIDE = 5


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(p: Path, data: Any) -> None:
    if isinstance(data, dict):
        run_id = str(data.get("run_id") or (
            p.parent.name if p.name == "status.json" else ""))
        data = stamp_payload(
            data,
            payload_type=("run_config" if p.name in {"run_config.json", "run.json"}
                          else "status" if p.name == "status.json" else p.stem),
            run_id=run_id,
            dataset_id=run_id,
            profile=RUNTIME_CONFIG.public_identity(),
            generator="src/exondomaincompare/framework/run_core_gene_analysis.py",
        )
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def read_json(p: Path, default: Any = None) -> Any:
    try:
        return json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:
        return default


def write_tsv(path: Path, columns: List[str], rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in columns})


# --------------------------------------------------------------------------- #
# gene config resolution
# --------------------------------------------------------------------------- #
def resolve_gene_config(gene_config_arg: Optional[str], gene_arg: Optional[str]) -> GeneConfig:
    if gene_config_arg:
        try:
            return load_gene_config(gene_config_arg)
        except GeneConfigError:
            return load_gene_config_lenient(gene_config_arg)
    if gene_arg:
        # Prefer an explicit hand-authored config if one exists (drafts/pilots).
        for rel in (f"{GENE_DRAFT_DIR}/{gene_arg}_core_only_pilot.yaml",
                    f"{GENE_DRAFT_DIR}/{gene_arg}_draft.yaml",
                    f"{GENE_CONFIG_DIR}/{gene_arg}.yaml"):
            p = PROJECT_ROOT / rel
            if p.is_file():
                return load_gene_config_lenient(p)
        # Truly generic: no file needed. Synthesize a core-only config in memory
        # so ANY protein-coding gene symbol can start an exploratory run. A run-
        # local gene_config.yaml is written by the create phase for provenance.
        print(f"[gene-config] No YAML for '{gene_arg}'; generating a generic core-only "
              "config automatically (no manual configuration required).")
        return build_generic_gene_config(gene_arg, generated_by="run_core_gene_analysis.py")
    raise SystemExit("Provide --gene-config <path> or --gene <SYMBOL>.")


# --------------------------------------------------------------------------- #
# species registry -> taxid
# --------------------------------------------------------------------------- #
def build_species_registry(run_dir: Path, species_id: str) -> Dict[str, Any]:
    outdir = run_dir / "results" / "01_species_registry"
    outdir.mkdir(parents=True, exist_ok=True)
    species_list = run_dir / "species_list.txt"
    cmd = [sys.executable, str(SCRIPTS_DIR / "build_species_registry_improved.py"),
           "--species_list", str(species_list), "--outdir", str(outdir)]
    subprocess.run(cmd, check=True, cwd=str(PROJECT_ROOT),
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    reg = outdir / "species_registry.tsv"
    rows = list(csv.DictReader(open(reg, encoding="utf-8"), delimiter="\t")) if reg.is_file() else []
    for r in rows:
        if r.get("species_id") == species_id:
            return r
    return rows[0] if rows else {"species_id": species_id, "taxid": ""}


# --------------------------------------------------------------------------- #
# locate a cached NCBI genome annotation for a taxid (offline reuse of raw data)
# --------------------------------------------------------------------------- #
def find_cached_annotation(taxid: str) -> Optional[Dict[str, str]]:
    if not taxid:
        return None
    rel = f"_ncbi_datasets_cache/ncbi_{taxid}/*/unzipped/ncbi_dataset/data/*/genomic.gff"
    # Search run caches AND the shared results/freeze cache (all RAW source data).
    search_roots = [
        (RUNS_ROOT, f"*/results/02_models/{rel}"),
        (PROJECT_ROOT / "results", f"**/02_models/{rel}"),
    ]
    for root, pattern in search_roots:
        if not root.exists():
            continue
        for gff in sorted(root.glob(pattern)):
            faa = gff.parent / "protein.faa"
            acc = gff.parent.name
            try:
                source = str(gff.relative_to(PROJECT_ROOT)).split("/")[1]
            except Exception:
                source = "cache"
            return {"gff": str(gff), "protein_faa": str(faa) if faa.is_file() else "",
                    "assembly_accession": acc, "source_run": source}
    return None


# --------------------------------------------------------------------------- #
# on-demand genome annotation retrieval (NCBI Datasets)
# --------------------------------------------------------------------------- #
def _datasets_bin() -> Optional[str]:
    return RUNTIME_CONFIG.executable("datasets")


def species_scientific_name(species_id: str, hint: str = "") -> str:
    h = (hint or "").strip()
    if h and " " in h and h.lower() != species_id.lower():
        return h
    parts = [p for p in species_id.split("_") if p]
    if len(parts) >= 2:
        return parts[0].capitalize() + " " + " ".join(parts[1:])
    return species_id.capitalize()


def download_annotation(run_dir: Path, species_id: str, taxid: str,
                        scientific_name: str, logline, timeout_s: int = 1800
                        ) -> Optional[Dict[str, str]]:
    dbin = _datasets_bin()
    if not dbin:
        logline(f"[{species_id}] NCBI Datasets CLI not found; cannot auto-retrieve genome annotation.")
        return None
    taxon = species_scientific_name(species_id, scientific_name)
    cache_root = run_dir / "results" / "02_models" / "_ncbi_datasets_cache"
    stage = cache_root / f"_download_{species_id}"
    stage.mkdir(parents=True, exist_ok=True)
    zip_path = stage / "ncbi_dataset.zip"
    dl_log = run_dir / "logs" / f"datasets_download_{species_id}.log"
    dl_log.parent.mkdir(parents=True, exist_ok=True)
    cmd = [dbin, "download", "genome", "taxon", taxon, "--reference",
           "--include", "gff3,protein", "--filename", str(zip_path)]
    logline(f"[{species_id}] Retrieving genome annotation from NCBI Datasets (taxon '{taxon}')…")
    try:
        with open(dl_log, "w", encoding="utf-8") as lf:
            lf.write("$ " + " ".join(cmd) + "\n")
            lf.flush()
            subprocess.run(cmd, check=True, stdout=lf, stderr=subprocess.STDOUT, timeout=timeout_s)
    except subprocess.TimeoutExpired:
        logline(f"[{species_id}] NCBI Datasets download timed out after {timeout_s}s "
                f"(see {dl_log.name}).")
        return None
    except subprocess.CalledProcessError as exc:
        logline(f"[{species_id}] NCBI Datasets download failed (exit {exc.returncode}; "
                f"see {dl_log.name}). The species name may be unrecognised or unavailable.")
        return None
    except Exception as exc:  # noqa: BLE001 - never crash the run on retrieval issues
        logline(f"[{species_id}] NCBI Datasets download error: {exc}")
        return None
    if not zip_path.is_file():
        logline(f"[{species_id}] NCBI Datasets produced no package.")
        return None

    resolved_taxid = (taxid or "").strip()
    acc = ""
    try:
        with zipfile.ZipFile(zip_path) as z:
            names = z.namelist()
            gff_names = [n for n in names if n.endswith("/genomic.gff")]
            if not gff_names:
                logline(f"[{species_id}] Downloaded package has no genomic.gff.")
                return None
            acc = gff_names[0].split("/data/")[1].split("/")[0]
            if not resolved_taxid:
                rep_name = next((n for n in names if n.endswith("assembly_data_report.jsonl")), "")
                if rep_name:
                    try:
                        for line in z.read(rep_name).decode("utf-8", "replace").splitlines():
                            obj = json.loads(line)
                            tx = str((obj.get("organism") or {}).get("taxId")
                                     or obj.get("taxId") or "").strip()
                            if tx and tx != "0":
                                resolved_taxid = tx
                                break
                    except Exception:
                        pass
            base = cache_root / (f"ncbi_{resolved_taxid}" if resolved_taxid
                                 else f"_download_{species_id}") / (acc or "assembly") / "unzipped"
            base.mkdir(parents=True, exist_ok=True)
            z.extractall(base)
    except Exception as exc:  # noqa: BLE001
        logline(f"[{species_id}] Could not unpack datasets package: {exc}")
        return None
    finally:
        try:
            zip_path.unlink()
            shutil.rmtree(stage, ignore_errors=True)
        except Exception:
            pass

    gff = base / "ncbi_dataset" / "data" / acc / "genomic.gff"
    if not gff.is_file():
        cand = sorted(base.glob("ncbi_dataset/data/*/genomic.gff"))
        gff = cand[0] if cand else gff
    if not gff.is_file():
        logline(f"[{species_id}] Extracted package missing genomic.gff.")
        return None
    faa = gff.parent / "protein.faa"
    size_mb = gff.stat().st_size // 1_000_000
    logline(f"[{species_id}] Retrieved assembly {acc} (taxid {resolved_taxid or 'unknown'}); "
            f"annotation {size_mb} MB cached for reuse.")
    return {"gff": str(gff), "protein_faa": str(faa) if faa.is_file() else "",
            "assembly_accession": acc, "source_run": "ncbi_datasets_download",
            "taxid": resolved_taxid}


# --------------------------------------------------------------------------- #
# gene-agnostic GFF parsing
# --------------------------------------------------------------------------- #
def _attrs(field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in field.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out




def resolve_gene_locus(gff_path: Path, gene_symbol: str, *,
                       scientific_name: str = "", taxid: str = "",
                       assembly_accession: str = "",
                       gene_lookup: Optional[Any] = None,
                       allow_network: bool = True) -> Dict[str, Any]:
    result = glr.resolve_gene_locus(
        Path(gff_path), gene_symbol, scientific_name=scientific_name, taxid=taxid,
        gene_lookup=gene_lookup, allow_network=allow_network,
        assembly_accession=assembly_accession)

    gene: Optional[Dict[str, Any]] = None
    if result.resolved and result.locus is not None:
        locus = result.locus
        gene = {
            "gene_id": locus.gene_id,
            # The annotation's own symbol, so downstream tables stay verifiable against
            # the source; the user-facing name is carried separately in `identity`.
            "gene_symbol": locus.symbol or gene_symbol,
            "dbxref": locus.dbxref,
            "biotype": locus.biotype,
            "seqid": locus.seqid, "start": locus.start, "end": locus.end,
            "strand": locus.strand,
        }

    matched_by = {
        glr.RESOLVED: "symbol",
        glr.AMBIGUOUS_FAMILY: "ambiguous",
    }.get(result.status, "none")
    if result.resolved:
        method = result.identity.resolution_method
        if method == glr.ROUTE_ANNOTATION_ALIAS:
            matched_by = "synonym"
        elif method in (glr.ROUTE_NCBI_ALIAS_GENEID, glr.ROUTE_NCBI_ALIAS_SYMBOL,
                        glr.ROUTE_NCBI_SYMBOL):
            matched_by = "ncbi_alias"

    return {
        "gene": gene,
        "matched_by": matched_by,
        "candidates": [c.locus.symbol for c in result.candidates
                       if c.decision in ("candidate", "accepted")],
        # The full locus record, including the coordinates and child counts that the
        # trimmed `gene` dict drops, so the written report is self-contained.
        "locus": result.locus.as_dict() if result.locus else None,
        "status": result.status,
        "message": result.message(),
        "detail": result.detail,
        "identity": result.identity.as_dict(),
        "candidate_inventory": [c.as_dict() for c in result.candidates],
        "ncbi_records": [r.as_dict() for r in result.ncbi_records],
        "ncbi_lookup_status": result.ncbi_lookup_status,
        "routes_attempted": list(result.routes_attempted),
    }


def parse_gene_models(gff_path: Path, gene_symbol: str, *,
                      scientific_name: str = "", taxid: str = "",
                      assembly_accession: str = "",
                      gene_lookup: Optional[Any] = None,
                      allow_network: bool = True) -> Dict[str, Any]:
    resolution = resolve_gene_locus(
        gff_path, gene_symbol, scientific_name=scientific_name, taxid=taxid,
        assembly_accession=assembly_accession, gene_lookup=gene_lookup,
        allow_network=allow_network)
    gene = resolution["gene"]
    transcripts: Dict[str, Dict[str, Any]] = {}
    rna_id_to_key: Dict[str, str] = {}
    if gene is None:
        return {"gene": None, "transcripts": [], "resolution": resolution}

    with open(gff_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            a = _attrs(f[8])
            if gene is not None and ftype in ("mRNA", "transcript"):
                if a.get("Parent") == gene["gene_id"]:
                    tid = a.get("transcript_id") or a.get("ID", "")
                    key = a.get("ID", tid)
                    tag = str(a.get("tag", "") or "")
                    # Gene-agnostic canonical heuristic: an explicit MANE / RefSeq
                    # Select tag, else a curated RefSeq transcript (NM_/NR_ prefix)
                    # over predicted (XM_/XR_). No gene-specific assumptions.
                    tl = tag.lower()
                    canonical = ("mane select" in tl or "refseq select" in tl
                                 or tid.startswith(("NM_", "NR_")))
                    transcripts[key] = {
                        "transcript_id": tid, "rna_key": key,
                        "biotype": a.get("gbkey", "") or a.get("transcript_biotype", ""),
                        "seqid": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                        "tag": tag, "canonical": canonical,
                        "exons": [], "cds": [], "protein_id": "",
                    }
                    rna_id_to_key[key] = key
            elif gene is not None and ftype == "exon":
                key = a.get("Parent", "")
                if key in transcripts:
                    transcripts[key]["exons"].append({
                        "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                        "exon_id": a.get("ID", ""),
                    })
            elif gene is not None and ftype == "CDS":
                key = a.get("Parent", "")
                if key in transcripts:
                    transcripts[key]["cds"].append({
                        "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                        "phase": f[7], "protein_id": a.get("protein_id", ""),
                    })
                    if a.get("protein_id"):
                        transcripts[key]["protein_id"] = a.get("protein_id")
    return {"gene": gene, "transcripts": list(transcripts.values()),
            "resolution": resolution}


def _exon_protein_map_for_transcript(tx: Dict[str, Any],
                                     protein_length_aa: Optional[int] = None
                                     ) -> List[Dict[str, Any]]:
    cds = [c for c in tx["cds"] if c.get("start") and c.get("end")]
    if not cds:
        return []
    tid = tx.get("transcript_id", "")
    parts = [{**c, "start": str(c["start"]), "end": str(c["end"]),
              "cds_id_source": c.get("cds_id", "")} for c in cds]
    features = build_cds_features_from_parts(
        source_db="NCBI", species_input="", species_canonical="",
        transcript_id_internal=tid, transcript_id_source=tid,
        translation_id_source=tx.get("protein_id", ""), parts=parts,
        coordinate_source="NCBI genomic GFF3",
        protein_length_aa=protein_length_aa)
    return [{
        "exon_number": int(f.coding_exon_order), "cds_start": int(f.start), "cds_end": int(f.end),
        "protein_start_aa": int(f.protein_start_aa),
        "protein_end_aa": int(f.protein_end_aa), "phase": f.phase,
        "normalized_strand": f.normalized_strand, "genomic_order": f.genomic_order,
        "transcript_order": f.transcript_order,
        "source_ordering_method": f.source_ordering_method,
    } for f in features if f.protein_start_aa and f.protein_end_aa]


# --------------------------------------------------------------------------- #
# protein FASTA extraction (from cached whole-proteome protein.faa)
# --------------------------------------------------------------------------- #
def load_fasta(protein_faa: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    if not protein_faa or not Path(protein_faa).is_file():
        return seqs
    cur = None
    buf: List[str] = []
    with open(protein_faa, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0].strip()
                buf = []
            else:
                buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


# --------------------------------------------------------------------------- #
# synteny neighbours (gene-agnostic)
# --------------------------------------------------------------------------- #
def extract_synteny_neighbors(gff_path: Path, gene: Dict[str, Any],
                              per_side: int = SYNTENY_PER_SIDE) -> List[Dict[str, Any]]:
    seqid = gene["seqid"]
    genes: List[Dict[str, Any]] = []
    with open(gff_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene" or f[0] != seqid:
                continue
            a = _attrs(f[8])
            if a.get("gene_biotype") not in (None, "protein_coding"):
                continue
            genes.append({
                "gene_id": a.get("ID", ""),
                "symbol": a.get("gene") or a.get("Name") or "",
                "start": int(f[3]), "end": int(f[4]), "strand": f[6],
            })
    genes.sort(key=lambda g: g["start"])
    idx = next((i for i, g in enumerate(genes) if g["gene_id"] == gene["gene_id"]), None)
    if idx is None:
        return []
    plus = gene["strand"] != "-"
    rows: List[Dict[str, Any]] = []
    # neighbours by genomic order; label side relative to transcription direction
    for offset in range(1, per_side + 1):
        for direction in (-1, +1):
            j = idx + direction * offset
            if 0 <= j < len(genes) and j != idx:
                g = genes[j]
                downstream = (j > idx) if plus else (j < idx)
                rows.append({
                    "neighbor_symbol": g["symbol"] or g["gene_id"],
                    "side": "downstream" if downstream else "upstream",
                    "order": offset,
                    "orientation": g["strand"],
                    "classification": "same_scaffold_ordered",
                    "status": "resolved" if g["symbol"] else "unresolved",
                    "seqid": seqid, "genomic_start": g["start"], "genomic_end": g["end"],
                    "distance_to_target": max(
                        0, gene["start"] - g["end"] if g["end"] < gene["start"]
                        else g["start"] - gene["end"]),
                })
    return rows


# --------------------------------------------------------------------------- #
# Run creation
# --------------------------------------------------------------------------- #
def _select_primary_pids(coding: List[Dict[str, Any]], lengths: Dict[str, int],
                         selection_method: str) -> Tuple[List[str], str, List[str]]:
    warnings: List[str] = []
    if not coding:
        return [], selection_method, warnings
    if selection_method == "all_protein_coding_transcripts":
        warnings.append("selection_method=all_protein_coding_transcripts: all coding "
                        "transcripts marked primary; no single canonical isoform chosen.")
        return [t["protein_id"] for t in coding], selection_method, warnings
    if selection_method == "canonical_if_available":
        canon = [t for t in coding if t.get("canonical")]
        if canon:
            return ([max(canon, key=lambda t: lengths.get(t["protein_id"], 0))["protein_id"]],
                    "canonical", warnings)
        warnings.append("canonical_if_available: no MANE/RefSeq-Select/curated transcript "
                        "found; fell back to the longest protein.")
        return ([max(coding, key=lambda t: lengths.get(t["protein_id"], 0))["protein_id"]],
                "longest_protein_fallback", warnings)
    return ([max(coding, key=lambda t: lengths.get(t["protein_id"], 0))["protein_id"]],
            "longest_protein", warnings)


def _collect_species_rows(cfg: GeneConfig, species_id: str, models: Dict[str, Any],
                          seqs: Dict[str, str], selection_method: str) -> Dict[str, Any]:
    gene = models["gene"]
    txs = models["transcripts"]

    def prot_len(pid: str) -> Optional[int]:
        s = seqs.get(pid)
        return len(s) if s else None

    coding = [t for t in txs if t.get("protein_id") and t.get("cds")]
    warnings: List[str] = []
    if not coding:
        warnings.append(f"[{species_id}] No protein-coding transcripts with CDS found for this gene.")
    lengths = {t["protein_id"]: (prot_len(t["protein_id"]) or 0) for t in coding}
    primary_pids, applied_rule, sel_warn = _select_primary_pids(coding, lengths, selection_method)
    warnings.extend(f"[{species_id}] {w}" for w in sel_warn)

    gene_models: List[Dict[str, Any]] = []
    isoforms: List[Dict[str, Any]] = []
    exon_rows: List[Dict[str, Any]] = []
    for t in txs:
        pid = t.get("protein_id", "")
        is_coding = bool(pid and t.get("cds"))
        gene_models.append({
            "analysis_id": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
            "species_id": species_id, "gene_id": gene["gene_id"],
            "transcript_id": t.get("transcript_id", ""), "protein_id": pid,
            "source": "ncbi_gff", "protein_length": prot_len(pid),
            "model_status": "protein_coding" if is_coding else (t.get("biotype") or "non_coding"),
            "notes": "",
        })
        if is_coding:
            isoforms.append({
                "species_id": species_id, "protein_id": pid,
                "transcript_id": t.get("transcript_id", ""), "isoform_label": "",
                "protein_length": prot_len(pid), "sequence_path": "results/core_gene_analysis/proteins_primary.faa",
                "primary_status": "primary" if pid in primary_pids else "alternative",
                "notes": "",
            })
            for r in _exon_protein_map_for_transcript(t, prot_len(pid)):
                exon_rows.append({
                    "species_id": species_id, "protein_id": pid,
                    "transcript_id": t.get("transcript_id", ""),
                    "exon_id": f"{pid}:cds{r['exon_number']}", "exon_number": r["exon_number"],
                    "cds_start": r["cds_start"], "cds_end": r["cds_end"],
                    "protein_start_aa": r["protein_start_aa"], "protein_end_aa": r["protein_end_aa"],
                    "phase": r["phase"], "confidence": "gff_cds_derived", "source": "ncbi_gff",
                    "normalized_strand": r["normalized_strand"],
                    "genomic_order": r["genomic_order"],
                    "transcript_order": r["transcript_order"],
                    "source_ordering_method": r["source_ordering_method"],
                })

    synteny_rows: List[Dict[str, Any]] = []
    synteny_reason = ""
    try:
        synteny_rows = extract_synteny_neighbors(Path(models["gff_path"]), gene)
        if not synteny_rows:
            synteny_reason = f"[{species_id}] No neighbouring protein-coding genes found on the gene's scaffold."
    except Exception as exc:  # noqa: BLE001
        synteny_reason = f"[{species_id}] Synteny extraction failed: {exc}"
    synteny_rows = [{**r, "species_id": species_id, "gene_symbol": cfg.gene_symbol, "source": "ncbi_gff"}
                    for r in synteny_rows]
    # The target locus is recorded explicitly, never inferred from a neighbour's
    # array position, and it is not one of the flanking neighbours.
    synteny_target_rows = [{
        "species_id": species_id,
        "gene_symbol": cfg.gene_symbol,
        "target_gene_id": gene.get("gene_id", ""),
        "target_symbol": gene.get("gene_symbol") or cfg.gene_symbol,
        "seqid": gene.get("seqid", ""),
        "genomic_start": gene.get("start", ""),
        "genomic_end": gene.get("end", ""),
        "strand": gene.get("strand", ""),
        "upstream_count_available": sum(1 for r in synteny_rows if r["side"] == "upstream"),
        "downstream_count_available": sum(1 for r in synteny_rows if r["side"] == "downstream"),
        "requested_neighbour_count": SYNTENY_PER_SIDE,
        "source": "ncbi_gff",
    }] if synteny_rows else []

    primary_seqs = [(pid, seqs.get(pid, ""), species_id) for pid in primary_pids if seqs.get(pid)]
    all_seqs = [(t["protein_id"], seqs.get(t["protein_id"], ""), species_id)
                for t in coding if seqs.get(t["protein_id"])]
    return {
        "species_id": species_id, "gene": gene, "coding": coding, "txs": txs,
        "gene_models": gene_models, "isoforms": isoforms, "exon_rows": exon_rows,
        "synteny_rows": synteny_rows, "synteny_reason": synteny_reason,
        "synteny_target_rows": synteny_target_rows,
        "primary_pids": primary_pids, "primary_seqs": primary_seqs, "all_seqs": all_seqs,
        "applied_rule": applied_rule, "warnings": warnings,
        "provenance": models.get("provenance", {}), "input_mode": models.get("input_mode", "local_cache"),
    }


def build_core_contract(run_dir: Path, cfg: GeneConfig, per_species: List[Dict[str, Any]],
                        selection_method: str) -> Dict[str, Any]:
    core_dir = run_dir / "results" / "core_gene_analysis"
    core_dir.mkdir(parents=True, exist_ok=True)

    collected = [_collect_species_rows(cfg, s["species_id"], s["models"], s["seqs"], selection_method)
                 for s in per_species]

    gene_models = [r for c in collected for r in c["gene_models"]]
    isoforms = [r for c in collected for r in c["isoforms"]]
    exon_rows = [r for c in collected for r in c["exon_rows"]]
    synteny_rows = [r for c in collected for r in c["synteny_rows"]]
    synteny_target_rows = [r for c in collected for r in c["synteny_target_rows"]]
    warnings = [w for c in collected for w in c["warnings"]]
    coding = [t for c in collected for t in c["coding"]]
    txs = [t for c in collected for t in c["txs"]]
    primary_pids = [p for c in collected for p in c["primary_pids"]]
    primary_seqs = [ps for c in collected for ps in c["primary_seqs"]]
    all_seqs = [s for c in collected for s in c["all_seqs"]]
    species_ids = [c["species_id"] for c in collected]
    n_species = len(collected)
    # representative single-species view for reports / single-species callers
    first = collected[0]
    gene = first["gene"]
    species_id = first["species_id"]
    applied_rule = first["applied_rule"]
    synteny_reason = "; ".join(c["synteny_reason"] for c in collected if c["synteny_reason"])

    write_tsv(core_dir / "gene_model_index.tsv",
              ["analysis_id", "gene_symbol", "species_id", "gene_id", "transcript_id",
               "protein_id", "source", "protein_length", "model_status", "notes"], gene_models)
    write_tsv(core_dir / "protein_isoform_index.tsv",
              ["species_id", "protein_id", "transcript_id", "isoform_label", "protein_length",
               "sequence_path", "primary_status", "notes"], isoforms)
    # The ordering columns are written beside the coordinates so a reader can see how the
    # transcript order was established, rather than having to trust that it was.
    write_tsv(core_dir / "exon_protein_map.tsv",
              ["species_id", "protein_id", "transcript_id", "exon_id", "exon_number",
               "cds_start", "cds_end", "protein_start_aa", "protein_end_aa", "phase",
               "confidence", "source", "normalized_strand", "genomic_order",
               "transcript_order", "source_ordering_method"], exon_rows)
    # domain/tm/boundary features are produced only after cluster annotation.
    write_tsv(core_dir / "domain_features.tsv",
              ["species_id", "protein_id", "domain_source", "domain_id", "domain_name",
               "start_aa", "end_aa", "score"], [])
    write_tsv(core_dir / "tm_features.tsv",
              ["species_id", "protein_id", "start_aa", "end_aa", "source"], [])
    write_tsv(core_dir / "exon_domain_boundary_distances.tsv",
              ["analysis_id", "gene_symbol", "species_id", "protein_id", "transcript_id",
               "exon_boundary_id", "boundary_position_aa", "nearest_domain_id",
               "nearest_domain_name", "nearest_domain_boundary_type", "distance_aa",
               "category", "source"], [])

    # synteny neighbours (best-effort, gene-agnostic; one block per species).
    # Flanking loci only — the target locus has its own table so it can never be
    # counted as one of its own neighbours.
    write_tsv(core_dir / "synteny_neighbors.tsv",
              ["species_id", "gene_symbol", "neighbor_symbol", "side", "order", "orientation",
               "classification", "source", "status", "seqid", "genomic_start",
               "genomic_end", "distance_to_target"], synteny_rows)
    write_tsv(core_dir / "synteny_target.tsv",
              ["species_id", "gene_symbol", "target_gene_id", "target_symbol", "seqid",
               "genomic_start", "genomic_end", "strand", "upstream_count_available",
               "downstream_count_available", "requested_neighbour_count", "source"],
              synteny_target_rows)

    # primary FASTA (one primary per species -> also the cross-species MSA input)
    fasta_reason = ""
    if not primary_seqs:
        if not coding:
            fasta_reason = "no_protein_coding_transcript"
        elif not all_seqs:
            fasta_reason = "protein_sequences_unavailable_in_source"
        else:
            fasta_reason = "primary_protein_id_not_found_in_proteome"
        warnings.append(f"No primary protein sequences resolved (reason={fasta_reason}).")
    faa = core_dir / "proteins_primary.faa"
    with open(faa, "w", encoding="utf-8") as fh:
        for pid, seq, sp in primary_seqs:
            fh.write(f">{pid} {cfg.gene_symbol}|{sp}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    # all-isoform FASTA (enables the optional exploratory event-candidate scan;
    # NOT sent to the cluster). Only isoforms whose sequence is available.
    all_faa = core_dir / "proteins_all_isoforms.faa"
    _n_all = 0
    with open(all_faa, "w", encoding="utf-8") as fh:
        for pid, seq, sp in all_seqs:
            _n_all += 1
            fh.write(f">{pid} {cfg.gene_symbol}|{sp}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    exon_map_reason = "" if exon_rows else "exon_map_unavailable"

    # mirror to the canonical closure/freeze path so existing status + cluster
    # plumbing (submit/fetch) detects the primary FASTA without special-casing.
    if primary_seqs:
        freeze_faa = run_dir / "results" / "13_final_pre_interpro_closure" / "freeze" / "final_pre_interpro_proteins_primary.faa"
        freeze_faa.parent.mkdir(parents=True, exist_ok=True)
        freeze_faa.write_text(faa.read_text(encoding="utf-8"), encoding="utf-8")

    # A dedicated, gene-agnostic collection report (per Part 3).
    collection_report = {
        "analysis_id": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
        "species_id": species_id, "species_ids": species_ids, "n_species": n_species,
        "generated_at": now_iso(),
        "input_mode": first.get("input_mode", "local_cache"),
        "source": first.get("provenance", {}),
        "per_species": [{
            "species_id": c["species_id"],
            "gene_locus": {
                "gene_id": c["gene"].get("gene_id"), "seqid": c["gene"].get("seqid"),
                "start": c["gene"].get("start"), "end": c["gene"].get("end"),
                "strand": c["gene"].get("strand"),
            },
            "selection_rule_applied": c["applied_rule"],
            "n_transcripts": len(c["txs"]), "n_protein_coding": len(c["coding"]),
            "n_primary_proteins": len(c["primary_seqs"]),
            "primary_protein_ids": [pid for pid, _, _ in c["primary_seqs"]],
            "source": c["provenance"],
        } for c in collected],
        "gene_locus": {
            "gene_id": gene.get("gene_id"), "seqid": gene.get("seqid"),
            "start": gene.get("start"), "end": gene.get("end"), "strand": gene.get("strand"),
        },
        "selection_method_requested": selection_method,
        "selection_rule_applied": applied_rule,
        "n_transcripts": len(txs), "n_protein_coding": len(coding),
        "n_primary_proteins": len(primary_seqs),
        "primary_protein_ids": [pid for pid, _, _ in primary_seqs],
        "exon_map_available": bool(exon_rows), "exon_map_reason": exon_map_reason,
        "synteny_available": bool(synteny_rows), "synteny_reason": synteny_reason,
        "fasta_reason": fasta_reason,
        "warnings": warnings,
    }
    write_json(core_dir / "core_model_collection_report.json", collection_report)

    report = {
        "analysis_id": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
        "dataset_id": f"run:{run_dir.name}", "adapter": "core_gene_runner",
        "contract_version": 1, "event_analysis_mode": cfg.event_analysis_mode,
        "has_event": cfg.has_event, "generated_at": now_iso(),
        "selection_method": selection_method,
        "summary": {
            "n_species": n_species, "n_transcripts": len(txs),
            "n_protein_coding": len(coding), "n_primary_proteins": len(primary_seqs),
            "n_exon_boundaries": 0, "n_domain_features": 0, "n_tm_features": 0,
            "n_synteny_neighbors": len(synteny_rows),
        },
        "species_ids": species_ids,
        "inputs": [c["provenance"] for c in collected],
        "outputs": ["gene_model_index.tsv", "protein_isoform_index.tsv", "exon_protein_map.tsv",
                    "domain_features.tsv", "tm_features.tsv", "synteny_neighbors.tsv",
                    "synteny_target.tsv",
                    "exon_domain_boundary_distances.tsv", "core_gene_report.json",
                    "evidence_register/coordinate_evidence_register.tsv",
                    "evidence_register/coordinate_evidence_register.json"],
        "synteny_reason": synteny_reason,
        "exon_map_reason": exon_map_reason,
        "fasta_reason": fasta_reason,
        "selection_rule_applied": applied_rule,
        "domain_status": "pending_cluster",
        "warnings": warnings, "failures": [],
        "note": "Core-only experimental pilot. Domain/TM/boundary features require the "
                "cluster InterProScan/pyTMHMM step. No event-specific analysis.",
    }
    write_json(core_dir / "core_gene_report.json", report)
    return {"core_dir": core_dir, "primary_fasta": faa, "n_primary": len(primary_seqs),
            "n_coding": len(coding), "n_synteny": len(synteny_rows),
            "exon_map_available": bool(exon_rows), "report": report,
            "primary_pids": primary_pids, "n_species": n_species,
            "species_ids": species_ids}


def phase_create(args: argparse.Namespace) -> int:
    cfg = resolve_gene_config(args.gene_config, args.gene)
    if cfg.has_event:
        raise SystemExit(f"Config '{cfg.analysis_id}' has a configured event region; the core-only "
                         "runner is for no-event configs. Use the FGFR2 pipeline for event analyses.")
    # One OR MORE species (generic). Accepts repeated --species and/or a single
    # comma/space-separated token; order is preserved, duplicates dropped. The
    # first species is the reference used for single-species report fields.
    raw_species: List[str] = []
    for tok in (args.species if isinstance(args.species, list) else [args.species]):
        # Each shell argument is one species name (may contain spaces, e.g.
        # "Gallus gallus"). Only a comma separates multiple species inside one token.
        for part in str(tok).split(","):
            if part.strip():
                raw_species.append(part.strip())
    species_ids: List[str] = []
    for tok in raw_species:
        sid = cnr.normalize_species_token(tok)
        if sid and sid not in species_ids:
            species_ids.append(sid)
    if not species_ids:
        raise SystemExit("No species provided. Use --species <name> [<name> ...].")
    species_id = species_ids[0]  # reference species (single-species report fields)
    multi = len(species_ids) > 1
    # The visible name and the directory name are separate: run_name is what the
    # user typed and may be empty, while the run_id slug always describes the run.
    run_name = run_labels.clean_run_name(args.run_name)
    reuse = (getattr(args, "reuse_run_id", "") or "").strip()
    if reuse:
        record = resolve_run_record(RUNTIME_CONFIG, reuse)
        run_id, run_dir = reuse, (
            record.path if record else RUNS_ROOT / reuse)
        if not run_dir.is_dir():
            raise SystemExit(f"--reuse-run-id: run '{reuse}' does not exist.")
        if record and record.read_only:
            raise SystemExit(
                "--reuse-run-id: registered legacy run is read-only; copy it first.")
        adapter = LegacyRunAdapter(run_dir, expected_run_id=reuse)
        # Keep the name the user gave the run when the retry does not supply one.
        if not run_name:
            existing = adapter.config()
            run_name = run_labels.clean_run_name(existing.get("run_name") or "")
        invalidate_derived_stages(run_dir)
    else:
        slug = run_labels.run_id_slug(run_name, gene_symbol=cfg.gene_symbol,
                                      species=species_ids)
        run_id, run_dir = cnr.unique_run_dir(cnr.generate_run_id(slug))
        initial_record = {
            "run_id": run_id, "dataset_id": run_id, "run_name": run_name,
            "case_study": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
            "analysis_id": cfg.analysis_id, "event_id": cfg.event_id,
            "event_type": cfg.event_type, "has_event": cfg.has_event,
            "support_level": cfg.support_level, "experimental": cfg.experimental,
            "workflow": "shared_exploratory", "created_at": now_iso(),
        }
        initial_status = {
            "run_id": run_id, "status": "created",
            "current_step": "run_created", "next_action": "automatic_precluster",
            "species_count": len(species_ids), "cluster_jobs": {},
        }
        gene_record = dict(cfg.raw)
        gene_record.update({
            "gene_symbol": cfg.gene_symbol, "analysis_id": cfg.analysis_id,
            "event_id": cfg.event_id, "event_type": cfg.event_type,
            "source": (
                f"repo:{cfg.source_path}" if cfg.source_path
                else "generated:generic-gene-config"),
        })
        RunLayout(run_dir, RunLayoutVersion.CANONICAL_V2).initialize(
            run_record=initial_record, status=initial_status,
            gene=gene_record, species=species_ids)
    adapter = LegacyRunAdapter(run_dir, expected_run_id=run_id)
    adapter.materialize_legacy_compatibility()

    log = run_dir / "logs" / "core_runner.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    def logline(m: str) -> None:
        with open(log, "a", encoding="utf-8") as fh:
            fh.write(f"[{now_iso()}] {m}\n")
        print(m)

    logline(f"Core-only pilot: analysis={cfg.analysis_id} gene={cfg.gene_symbol} "
            f"species={','.join(species_ids)}")

    # Explicit 'running' status the moment collection starts, so the UI shows a
    # live pipeline (never 'created — not started') while species are processed.
    _write_status_running(run_dir, run_id, cfg, "collecting_models",
                          species_count=len(species_ids),
                          detail=f"Collecting gene / protein models for {len(species_ids)} species…")

    # Per-species collection: registry -> cached annotation -> gene models.
    # --gff/--protein-faa user override only applies to a single-species run.
    per_species: List[Dict[str, Any]] = []
    taxids: Dict[str, str] = {}
    sci_names: Dict[str, str] = {}
    for sid in species_ids:
        reg = build_species_registry(run_dir, sid)
        taxid = str(reg.get("taxid", "") or "").strip()
        taxids[sid] = taxid
        sci_names[sid] = str(reg.get("scientific_name") or reg.get("ncbi_species") or "").strip()
        logline(f"[{sid}] Species registry: taxid={taxid or '(unknown)'}")

        ann = None
        if args.gff and not multi:
            ann = {"gff": args.gff, "protein_faa": args.protein_faa or "",
                   "assembly_accession": "user_supplied", "source_run": "user_supplied"}
        else:
            if taxid:
                ann = find_cached_annotation(taxid)
            # Auto mode: no local cache -> retrieve the annotation from NCBI Datasets
            # (gff3 + protein only) so ANY resolvable species works without a manual
            # pre-download. local_cache mode never downloads (stays fully offline).
            if ann is None and args.input_mode == "auto":
                _write_status_running(run_dir, run_id, cfg, "retrieving_models",
                                      species_count=len(species_ids),
                                      detail=f"Retrieving genome annotation for {sid} from NCBI Datasets…")
                ann = download_annotation(run_dir, sid, taxid, sci_names.get(sid, ""), logline)
                if ann and ann.get("taxid") and not taxid:
                    taxid = ann["taxid"]
                    taxids[sid] = taxid
        if ann is None:
            if taxid:
                reason = (f"[{sid}] No genome annotation available for taxid {taxid}: not in the local "
                          "cache and the NCBI Datasets retrieval did not succeed. Provide --gff/--protein-faa "
                          "or check connectivity.")
            elif args.input_mode == "auto":
                reason = (f"[{sid}] Could not resolve or retrieve a genome annotation for '{sid}'. "
                          "The species name may be unrecognised by NCBI, or the assembly is unavailable. "
                          "Check the scientific name (e.g. 'Salmo salar') or provide --gff/--protein-faa.")
            else:
                reason = (f"[{sid}] No cached genome annotation for '{sid}' (taxid unresolved) and "
                          "auto-retrieval is disabled (--input-mode local_cache). Provide --gff/--protein-faa.")
            _write_status_blocked(run_dir, cfg, sid, taxid, reason,
                                  species_count=len(species_ids))
            logline(f"[{sid}] BLOCKED: {reason}")
            print(f"\nRun folder created (blocked): {run_dir.relative_to(PROJECT_ROOT)}")
            return 3
        logline(f"[{sid}] Annotation source: {ann['gff']} (assembly {ann['assembly_accession']}, "
                f"from {ann['source_run']})")

        _write_status_running(run_dir, run_id, cfg, "resolving_gene_locus",
                              species_count=len(species_ids),
                              detail=f"Resolving the {cfg.gene_symbol} locus for {sid}…")
        models = parse_gene_models(Path(ann["gff"]), cfg.gene_symbol,
                                   scientific_name=sci_names.get(sid, ""), taxid=taxid,
                                   assembly_accession=ann["assembly_accession"])
        resolution = models.get("resolution", {}) or {}
        write_gene_resolution_report(run_dir, sid, resolution)
        if models["gene"] is None:
            # The cascade knows which of several distinguishable things went wrong; report
            # that instead of the one sentence that used to cover all of them.
            reason = f"[{sid}] {resolution.get('message') or 'Gene locus resolution failed.'}"
            _write_status_blocked(run_dir, cfg, sid, taxid, reason,
                                  species_count=len(species_ids),
                                  failed_step="gene_locus_resolution",
                                  resolution=resolution)
            logline(f"[{sid}] BLOCKED ({resolution.get('status')}): {reason}")
            return 3
        identity = resolution.get("identity", {}) or {}
        if identity.get("symbol_differs_from_source"):
            logline(f"[{sid}] Resolved '{cfg.gene_symbol}' through "
                    f"{identity.get('resolution_method')} to annotation locus "
                    f"'{identity.get('resolved_official_symbol')}' "
                    f"(GeneID {identity.get('resolved_gene_id') or 'unknown'}, "
                    f"{identity.get('source_description') or 'no source description'}).")
        models["gff_path"] = ann["gff"]
        models["gene_identity"] = identity
        models["provenance"] = {"species_id": sid, "assembly_accession": ann["assembly_accession"],
                                "genomic_gff": ann["gff"], "protein_faa": ann["protein_faa"],
                                "source_run": ann["source_run"], "taxid": taxid}
        models["input_mode"] = args.input_mode
        logline(f"[{sid}] Gene {cfg.gene_symbol}: {models['gene']['seqid']}:{models['gene']['start']}-"
                f"{models['gene']['end']} ({models['gene']['strand']}); "
                f"{len(models['transcripts'])} transcripts")
        seqs = load_fasta(Path(ann["protein_faa"])) if ann["protein_faa"] else {}
        per_species.append({"species_id": sid, "models": models, "seqs": seqs})

    taxid = taxids[species_id]
    ref_models = per_species[0]["models"]
    built = build_core_contract(run_dir, cfg, per_species, args.selection_method)
    logline(f"Core contract: {built['n_species']} species, {built['n_coding']} protein-coding, "
            f"{built['n_primary']} primary, synteny={built['n_synteny']}, "
            f"exon_map={'yes' if built['exon_map_available'] else 'no'}")

    # run_config first so the evidence layer + indices can read gene/analysis meta.
    _write_run_config(run_dir, run_id, run_name, cfg, species_ids, taxids, ref_models, built,
                      sci_names=sci_names, input_mode=args.input_mode,
                      gene_identities={s["species_id"]: s["models"].get("gene_identity", {})
                                       for s in per_species})

    # exploratory isoform-difference evidence layer (optional, non-blocking; never
    # a validated event region and never activates event-specific views).
    _build_exploratory_evidence(run_dir, cfg, logline)

    # generic core indices (picks up evidence/cluster counts for the capability report)
    core_src = CoreSource(built["core_dir"], f"run:{run_id}", cfg)
    build_core_indices(core_src, run_dir / "website_indices" / "generic")
    logline("Built generic core indices under website_indices/generic/")

    # Honest, milestone-derived status: a run must never look analysis-ready if
    # required core outputs are missing. If collection produced no gene models or
    # no primary FASTA, the run is failed/incomplete — never cluster_required.
    if built["n_coding"] == 0 or built["report"].get("summary", {}).get("n_protein_coding", 0) == 0:
        _write_status_collection_failed(
            run_dir, run_id, cfg,
            "core_model_collection produced no protein-coding transcripts with CDS.")
        logline("FAILED: no protein-coding transcripts collected.")
        return 4
    if built["n_primary"] == 0:
        _write_status_incomplete(
            run_dir, run_id, cfg,
            f"No primary protein FASTA could be produced "
            f"(reason={built['report'].get('fasta_reason') or 'unknown'}).")
        logline("INCOMPLETE: no primary protein FASTA.")
        return 4
    _write_status_cluster_required(run_dir, run_id, cfg, species_id, built["n_primary"],
                                   built["n_species"])

    # Materialize the shared generic analysis and website layers.
    _run_generic_pipeline(run_id, logline)

    # Build the shared coordinate model after the generic layer refreshes its indices.
    _build_coordinate_model_and_figures(run_dir, logline)

    evidence = build_coordinate_evidence_register(run_dir)
    logline(f"Coordinate evidence register: "
            f"{evidence['counts']['total_records']} records ({evidence['register_phase']}).")

    # An isolated test root lives outside the project, so report the path as-is there.
    shown = (run_dir.relative_to(PROJECT_ROOT)
             if run_dir.is_relative_to(PROJECT_ROOT) else run_dir)
    print(f"\nOK  core-only pilot run created: {shown}")
    print(f"    analysis={cfg.analysis_id}  support_level={cfg.support_level}  experimental={cfg.experimental}")
    print(f"    primary proteins={built['n_primary']}  synteny neighbours={built['n_synteny']}")
    print(f"    exon->protein map: {'available' if built['exon_map_available'] else 'unavailable'}")
    print(f"    synteny: {'available' if built['n_synteny'] else 'not_computed'}")
    print("\n    Next (cluster required):")
    print(f"      .venv/bin/edc cluster roundtrip --run-id {run_id}")
    return 0


def _run_generic_pipeline(run_id: str, logline) -> None:
    try:
        import sys as _sys
        scripts_dir = str(PROJECT_ROOT / "scripts")
        if scripts_dir not in _sys.path:
            _sys.path.insert(0, scripts_dir)
        from generic_gene.run_generic_gene_analysis import run as _run_generic
        res = _run_generic(run_id)
        logline(f"Shared generic pipeline: event_layer={res.get('event_layer_type')} "
                f"clusters={res.get('n_clusters')} figures={res.get('figures_generated')} "
                f"msa={res.get('msa_status')}")
    except Exception as exc:  # pragma: no cover - non-blocking
        logline(f"WARN shared generic pipeline skipped: {exc}")


def _build_exploratory_evidence(run_dir: Path, cfg: GeneConfig, logline) -> None:
    try:
        from exondomaincompare.framework.scan_isoform_event_candidates import scan as _scan, COLUMNS as _CAND_COLS
        from exondomaincompare.framework.build_event_region_evidence import build_evidence as _build_ev
        from exondomaincompare.framework.cluster_event_region_evidence import build_clusters as _build_cl
        core_dir = run_dir / "results" / "core_gene_analysis"
        cand_rows = _scan(run_dir)
        write_tsv(core_dir / "event_candidate_regions.tsv", _CAND_COLS, cand_rows)
        ev = _build_ev(core_dir, cfg.analysis_id, cfg.gene_symbol)
        cl = _build_cl(core_dir, gap=5)
        logline(f"Exploratory evidence: {ev['n_evidence']} row(s), "
                f"{cl['n_clusters']} candidate cluster(s) (not validated events).")
    except Exception as exc:  # noqa: BLE001 - evidence layer is optional
        logline(f"WARN exploratory evidence layer failed (non-blocking): {exc}")


def _merged_gene_identity(gene_symbol: str,
                          identities: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    source_symbols = {sid: (rec.get("resolved_official_symbol") or "")
                      for sid, rec in identities.items()
                      if rec.get("symbol_differs_from_source")}
    methods = sorted({rec.get("resolution_method", "") for rec in identities.values()
                      if rec.get("resolution_method")})
    descriptions = sorted({rec.get("source_description", "") for rec in identities.values()
                           if rec.get("source_description")})
    return {
        "requested_gene_symbol": gene_symbol,
        "display_symbol": gene_symbol,
        "source_symbols_by_species": source_symbols,
        "any_symbol_differs_from_source": bool(source_symbols),
        "resolution_methods": methods,
        "source_descriptions": descriptions,
    }


def _write_run_config(run_dir: Path, run_id: str, run_name: str, cfg: GeneConfig,
                      species_ids: List[str], taxids: Dict[str, str], models: Dict[str, Any],
                      built: Dict[str, Any], sci_names: Optional[Dict[str, str]] = None,
                      input_mode: str = "auto",
                      gene_identities: Optional[Dict[str, Dict[str, Any]]] = None) -> None:
    core_faa_rel = "run:results/core_gene_analysis/proteins_primary.faa"
    # For a generated (generic) config there is no repo YAML; the authoritative
    # config is the run-local copy written during scaffolding.
    generated = cfg.source_path is None
    cfg_rel = (
        portable_path_reference(
            cfg.source_path, repository_root=PROJECT_ROOT, run_root=run_dir
        )
        if cfg.source_path else "run:gene_config.yaml"
    )
    selection_policy = built["report"]["selection_method"]
    data = {
        "run_id": run_id, "run_name": run_name,
        "case_study": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
        "analysis_id": cfg.analysis_id, "event_id": cfg.event_id, "event_type": cfg.event_type,
        "event_status": cfg.event_status, "has_event": cfg.has_event,
        "support_level": cfg.support_level, "experimental": cfg.experimental,
        "gene_config": cfg_rel, "species_count": len(species_ids),
        "species_ids": species_ids, "species_taxids": taxids,
        "species_scientific_names": dict(sci_names or {}),
        # What the user asked for versus what the assembly calls it. Kept side by side so
        # the UI can show 'HBA' with 'NCBI annotation symbol: LOC122209636' underneath
        # rather than replacing one with the other.
        "gene_identity_by_species": dict(gene_identities or {}),
        "gene_identity": _merged_gene_identity(cfg.gene_symbol, gene_identities or {}),
        "species_list_path": "run:species_list.txt",
        "run_dir": "run:.", "results_dir": "run:results",
        "website_indices_dir": "run:website_indices",
        "created_at": now_iso(), "created_by": f"local_profile:{RUNTIME_CONFIG.local_profile_name}",
        "configuration_profile": RUNTIME_CONFIG.public_identity(),
        "run_mode": "core_only_pilot",
        # Workflow routing metadata (mirrors the central analysis router). This is
        # always the shared exploratory workflow; FGFR2 uses its own creator.
        "workflow": "shared_exploratory",
        "workflow_mode": "shared_exploratory",
        "event_layer": "exploratory",
        "generated_config": generated,
        "config_provenance": portable_runtime_record(
            cfg.raw.get("provenance", {}) if generated else {},
            repository_root=PROJECT_ROOT,
            run_root=run_dir,
        ),
        # Scientific selection policies (generic defaults from the working pipeline).
        "primary_model_selection_policy": selection_policy,
        "cluster_annotation_scope": "primary_isoforms",
        "cluster_mode": "external_manual",
        "cluster_input_fasta": core_faa_rel,
        "primary_fasta_path": core_faa_rel,
        "primary_fasta_count": built["n_primary"],
        "selection_method": selection_policy,
        # Input-mode metadata (Part 8): recorded so user_files can be added later
        # without a schema change. auto | local_cache | user_files.
        "input_mode": models.get("input_mode", input_mode),
        "user_files": {
            "gff": (portable_path_reference(
                models.get("provenance", {}).get("genomic_gff"),
                repository_root=PROJECT_ROOT, run_root=run_dir)
                    if models.get("provenance", {}).get("genomic_gff") else ""),
            "protein_faa": (portable_path_reference(
                models.get("provenance", {}).get("protein_faa"),
                repository_root=PROJECT_ROOT, run_root=run_dir)
                            if models.get("provenance", {}).get("protein_faa") else ""),
            "genome_fasta": "",
        },
        "annotation_provenance": portable_runtime_record(
            models.get("provenance", {}),
            repository_root=PROJECT_ROOT,
            run_root=run_dir,
        ),
        "human_reference": {},
        "notes": "Experimental core-only pilot (no event region). Not a validated analysis.",
    }
    # Same stamp as the FGFR2 creator, resolved from the gene symbol alone, so the
    # Explorer and Gallery architecture is fixed at creation for every gene.
    production_contract.stamp(data)
    target = run_dir / (
        "run.json" if (run_dir / "run.json").is_file() else "run_config.json")
    if target.name == "run.json":
        current = read_json(target, {}) or {}
        current.update(data)
        current["layout_version"] = "canonical-2.0"
        current["species_list_path"] = "run:config/species.tsv"
        current["gene_config"] = "run:config/gene.json"
        data = current
    write_json(target, data)


def _write_status_running(run_dir: Path, run_id: str, cfg: GeneConfig, step: str,
                          species_count: int = 1, detail: str = "") -> None:
    write_json(run_dir / "status.json", {
        "run_id": run_id, "status": "running", "current_step": step,
        "run_mode": "core_only_pilot", "experimental": True,
        "support_level": cfg.support_level, "has_event": cfg.has_event,
        "event_status": cfg.event_status, "species_count": species_count,
        "pre_interpro_status": "running",
        "primary_fasta_status": "not_available", "primary_fasta_count": 0,
        "cluster_analysis_status": "not_started", "post_interpro_status": "not_started",
        "next_action": "wait_pre_interpro", "detail": detail,
        "human_reference": {}, "cluster_jobs": {}, "last_updated": now_iso(),
    })


def _write_status_cluster_required(run_dir: Path, run_id: str, cfg: GeneConfig,
                                   species_id: str, primary_count: int,
                                   species_count: int = 1) -> None:
    status = {
        "run_id": run_id, "status": "cluster_required",
        "current_step": "core_fasta_ready",
        "run_mode": "core_only_pilot", "experimental": True,
        "support_level": cfg.support_level,
        "has_event": cfg.has_event, "event_status": cfg.event_status,
        "species_count": species_count,
        "pre_interpro_status": "complete",
        "primary_fasta_status": "available", "primary_fasta_count": primary_count,
        "review_fasta_status": "not_available", "review_fasta_count": 0,
        "cluster_analysis_status": "not_started", "cluster_fetch_status": "not_started",
        "post_interpro_status": "not_started", "website_indices_status": "core_only",
        "next_action": "run_cluster_roundtrip_command",
        "cluster_profile": RUNTIME_CONFIG.lrz_profile_name,
        "cluster_command": RUNTIME_CONFIG.command([
            ".venv/bin/python", "scripts/edc.py", "cluster", "roundtrip",
            "--run-id", run_id,
            "--local-profile", RUNTIME_CONFIG.local_profile_name,
            "--lrz-profile", RUNTIME_CONFIG.lrz_profile_name,
        ]),
        "human_reference": {}, "cluster_jobs": {}, "last_updated": now_iso(),
    }
    write_json(run_dir / "status.json", status)


def _write_status_collection_failed(run_dir: Path, run_id: str, cfg: GeneConfig,
                                    reason: str) -> None:
    write_json(run_dir / "status.json", {
        "run_id": run_id, "status": "core_model_collection_failed",
        "current_step": "core_model_collection_failed",
        "run_mode": "core_only_pilot", "experimental": True,
        "support_level": cfg.support_level, "has_event": cfg.has_event,
        "event_status": cfg.event_status, "species_count": 1,
        "pre_interpro_status": "failed",
        "primary_fasta_status": "not_available", "primary_fasta_count": 0,
        "failed_step": "core_model_collection", "failed_reason": reason,
        "error": reason, "next_action": "inspect_logs",
        "human_reference": {}, "cluster_jobs": {}, "last_updated": now_iso(),
    })


def _write_status_incomplete(run_dir: Path, run_id: str, cfg: GeneConfig,
                             reason: str) -> None:
    write_json(run_dir / "status.json", {
        "run_id": run_id, "status": "incomplete",
        "current_step": "core_outputs_incomplete",
        "run_mode": "core_only_pilot", "experimental": True,
        "support_level": cfg.support_level, "has_event": cfg.has_event,
        "event_status": cfg.event_status, "species_count": 1,
        "pre_interpro_status": "incomplete",
        "primary_fasta_status": "not_available", "primary_fasta_count": 0,
        "failed_step": "core_primary_fasta", "failed_reason": reason,
        "error": reason, "next_action": "inspect_logs",
        "human_reference": {}, "cluster_jobs": {}, "last_updated": now_iso(),
    })


#: Kept across an in-place retry, because none of it is derived from the failed stage and
#: all of it is expensive to obtain again: the taxon resolution already succeeded, the
#: cached annotation is hundreds of megabytes, and the InterProScan / pyTMHMM directories
#: hold results fetched from the cluster. Whether those cluster results still match the
#: rebuilt FASTA is a staleness question for the post-cluster stage, not a reason to throw
#: them away before it can look.
_RETRY_PRESERVED = ("01_species_registry", "_ncbi_datasets_cache",
                    "14_interproscan", "15_exon_domain_boundary_post_interpro")


def invalidate_derived_stages(run_dir: Path) -> List[str]:
    removed: List[str] = []
    results = run_dir / "results"
    if results.is_dir():
        for child in sorted(results.iterdir()):
            if child.name in _RETRY_PRESERVED:
                continue
            if child.name == "02_models":
                # Keep the cached archive; drop the tables derived from it.
                for inner in sorted(child.iterdir()):
                    if inner.name in _RETRY_PRESERVED:
                        continue
                    shutil.rmtree(inner, ignore_errors=True) if inner.is_dir() \
                        else inner.unlink(missing_ok=True)
                    removed.append(str(inner.relative_to(run_dir)))
                continue
            shutil.rmtree(child, ignore_errors=True) if child.is_dir() \
                else child.unlink(missing_ok=True)
            removed.append(str(child.relative_to(run_dir)))
    for name in ("website_indices", "figures"):
        target = run_dir / name
        if target.is_dir():
            shutil.rmtree(target, ignore_errors=True)
            removed.append(name)
    return removed


def write_gene_resolution_report(run_dir: Path, species_id: str,
                                 resolution: Dict[str, Any]) -> Path:
    out_dir = run_dir / "results" / "02_models"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"gene_resolution_{species_id}.json"
    write_json(path, {"species_id": species_id, "generated_at": now_iso(), **resolution})

    inventory = resolution.get("candidate_inventory") or []
    if inventory:
        fields = ["symbol", "source_gene_id", "seqid", "start", "end", "strand",
                  "biotype", "feature_type", "is_loc_labelled", "is_pseudogene",
                  "transcript_count", "protein_count", "routes", "decisive_routes",
                  "ncbi_gene_id", "ncbi_official_symbol", "ncbi_aliases",
                  "ncbi_description", "orthology_evidence", "family_evidence",
                  "decision", "reason"]
        tsv = out_dir / f"gene_candidates_{species_id}.tsv"
        with open(tsv, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            for row in inventory:
                writer.writerow({k: (", ".join(v) if isinstance(v, list) else v)
                                 for k, v in row.items() if k in fields})
    return path


def _write_status_blocked(run_dir: Path, cfg: GeneConfig, species_id: str, taxid: str,
                          reason: str, species_count: int = 1,
                          failed_step: str = "core_model_collection",
                          resolution: Optional[Dict[str, Any]] = None) -> None:
    status = {
        "run_id": run_dir.name, "status": "core_model_collection_failed",
        "current_step": "core_collection_blocked",
        "run_mode": "core_only_pilot", "experimental": True,
        "support_level": cfg.support_level, "has_event": cfg.has_event,
        "species_count": species_count, "pre_interpro_status": "failed",
        "primary_fasta_status": "not_available", "primary_fasta_count": 0,
        "failed_step": failed_step, "failed_species": species_id,
        "failed_reason": reason, "error": reason, "next_action": "retry_or_inspect",
        "human_reference": {}, "last_updated": now_iso(),
    }
    if resolution:
        # The distinguishable cause, so the UI can separate "this species has no such
        # gene" from "the alias resolved but the assembly locus could not be mapped".
        status["failure_class"] = resolution.get("status", "")
        status["failure_detail"] = resolution.get("detail", "")
        status["gene_identity"] = resolution.get("identity", {})
        status["routes_attempted"] = resolution.get("routes_attempted", [])
        status["candidate_count"] = len(resolution.get("candidate_inventory") or [])
        status["retry_supported"] = True
    write_json(run_dir / "status.json", status)
    # Preserve a richer config; otherwise persist enough detail for an honest UI status.
    cfg_path = run_dir / "run_config.json"
    if cfg_path.is_file():
        try:
            existing = json.loads(cfg_path.read_text(encoding="utf-8"))
            existing["notes"] = f"Blocked: {reason}"
            write_json(cfg_path, existing)
            return
        except Exception:
            pass
    write_json(cfg_path, {
        "run_id": run_dir.name, "case_study": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
        "analysis_id": cfg.analysis_id, "gene_config": cfg.source_path,
        "support_level": cfg.support_level, "experimental": cfg.experimental,
        "has_event": cfg.has_event, "run_mode": "core_only_pilot",
        "species_count": species_count,
        "species_list_path": "run:species_list.txt",
        "created_at": now_iso(), "notes": f"Blocked: {reason}",
    })


# --------------------------------------------------------------------------- #
# POST phase — build domain/TM/boundary features from cluster outputs
# --------------------------------------------------------------------------- #
# Generic, gene-independent boundary classifications (no cassette / event terms).
# A boundary whose protein carries no domain annotation is classified "unknown"
# (uncertain), NOT "outside_domain": absence of a domain model is not evidence
# that the boundary lies outside annotated domain space.
BOUNDARY_CATEGORIES = {
    "aligned_to_domain_boundary": "exact_edge",
    "near_domain_boundary": "near_edge",
    "within_domain": "inside_domain",
    "between_domains": "outside_domain",
    "review_or_missing": "unknown",
}


def _classify_boundary(pos: int, domains: List[Dict[str, Any]], threshold: int):
    if not domains:
        return None, "", "unknown", None
    normalized = [{
        "label": d.get("domain_name") or d.get("domain_id", ""),
        "dclass": d.get("domain_class_simplified") or d.get("feature_type", ""),
        "start": d.get("start_aa"), "end": d.get("end_aa"),
    } for d in domains]
    result = classify_exon_domain_boundary(
        pos, normalized, aligned_max=0, near_max=threshold)
    best = next((d for d in domains
                 if d.get("start_aa") == result.get("start")
                 and d.get("end_aa") == result.get("end")
                 and (d.get("domain_name") or d.get("domain_id", "")) == result.get("label")), None)
    edge = f"domain_{result.get('edge')}" if result.get("edge") else ""
    return best, edge, BOUNDARY_CATEGORIES.get(result["boundary_class"], "unknown"), result.get("dist")


def _norm_acc(acc: str) -> str:
    a = (acc or "").strip()
    return re.sub(r"\.\d+$", "", a)


def _protein_species_map(core_dir: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for name in ("protein_isoform_index.tsv", "gene_model_index.tsv", "exon_protein_map.tsv"):
        path = core_dir / name
        if not path.is_file():
            continue
        for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
            pid = (r.get("protein_id") or "").strip()
            sp = (r.get("species_id") or "").strip()
            if pid and sp:
                mapping.setdefault(pid, sp)
                mapping.setdefault(_norm_acc(pid), sp)
    return mapping


def _build_coordinate_model_and_figures(run_dir: Path, logline) -> bool:
    try:
        from exondomaincompare.shared_gene_analysis.protein_coordinate_model import build_models_for_run
        from exondomaincompare.shared_gene_analysis.validate_protein_coordinate_model import validate_index
        model_out = run_dir / "website_indices" / "generic" / "protein_coordinate_model.json"
        index = build_models_for_run(run_dir)
        errors = validate_index(index, core_dir=run_dir / "results" / "core_gene_analysis")
        if errors:
            # A model that fails validation is a technical failure of this run, not a
            # cue to serve an older coordinate track: falling back silently used to
            # leave the modern Domain Architecture and Boundary views behind while the
            # rest of the run reported success. Record what failed and stay failed.
            model_out.parent.mkdir(parents=True, exist_ok=True)
            (model_out.parent / "protein_coordinate_model_errors.json").write_text(
                json.dumps({"status": "technically_missing",
                            "reason": "coordinate_model_validation_failed",
                            "n_errors": len(errors), "errors": errors[:50]}, indent=2))
            model_out.unlink(missing_ok=True)
            logline(f"FAIL: coordinate model validation reported {len(errors)} issue(s); "
                    "Domain Architecture and Boundary stay technically_missing "
                    "(no legacy fallback). First: " + str(errors[0]))
            return False
        (model_out.parent / "protein_coordinate_model_errors.json").unlink(missing_ok=True)
        model_out.write_text(json.dumps(index, indent=2))
        logline(f"Built validated protein-coordinate model "
                f"({index.get('n_models', 0)} model(s)).")
        # One shared stage list for every run, so a species Scope of a
        # multi-species run is built by exactly the stages that build a
        # standalone single-species Gallery.
        from plotting.figure_sequence import run_figure_stages
        run_figure_stages(run_dir, model_out, logline)
        return True
    except Exception as cm_err:  # pragma: no cover - model optional
        logline(f"WARN: protein-coordinate model build skipped: {cm_err}")
        return False


def phase_post(args: argparse.Namespace) -> int:
    record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    run_dir = (
        record.path if record else (RUNS_ROOT / args.run_id).resolve())
    if not run_dir.is_dir():
        raise SystemExit(f"Run not found: {args.run_id}")
    if record and record.read_only:
        raise SystemExit("Registered legacy run is read-only; copy it before postcluster.")
    LegacyRunAdapter(
        run_dir, expected_run_id=args.run_id).materialize_legacy_compatibility()
    rc = read_json(run_dir / "run_config.json", {}) or {}
    try:
        cfg = load_gene_config_lenient(run_dir / "gene_config.yaml")
    except GeneConfigError:
        cfg = load_gene_config_lenient(rc.get("gene_config") or "configs/genes/drafts/FGFR1_core_only_pilot.yaml")

    core_dir = run_dir / "results" / "core_gene_analysis"
    ips_out = run_dir / "results" / "14_interproscan" / "primary" / "output"
    tm_out = run_dir / "results" / "15_exon_domain_boundary_post_interpro" / "pytmhmm_primary" / "output"

    # map the primary protein ids so we can recover species/protein from any
    # pipe-delimited FASTA header the cluster echoed back (real InterProScan uses
    # the FASTA header verbatim as the sequence accession in column 1).
    primary_ids = _primary_protein_ids(core_dir)

    interpro_qc = _interpro_parse_qc(ips_out, primary_ids)
    pytmhmm_qc = _pytmhmm_parse_qc(tm_out, primary_ids)
    # Returned-sequence inventory (Part 1C). A species whose primary protein came back
    # with no features is indistinguishable from a species that was never analysed
    # unless the coverage is recorded explicitly, per species, against what was
    # submitted. This is what makes a silently dropped species visible.
    inventory = _returned_sequence_inventory(core_dir, interpro_qc, pytmhmm_qc)

    qc_path = run_dir / "results" / "15_domain_architecture" / "post_cluster_qc.json"
    write_json(qc_path, {
        "interproscan": {k: v for k, v in interpro_qc.items() if k != "rows"},
        "pytmhmm": {k: v for k, v in pytmhmm_qc.items() if k != "rows"},
        "returned_sequence_inventory": inventory,
        "synthetic_outputs_accepted": False, "generated_at": now_iso(),
    })
    missing_species = [r["species_id"] for r in inventory["species"]
                       if r["interproscan_status"] == "no_returned_sequence"]
    if missing_species:
        reason = (f"cluster returned no InterProScan sequence for {missing_species}; "
                  f"submitted primaries were "
                  f"{[r['protein_id'] for r in inventory['species']]}")
        st = read_json(run_dir / "status.json", {}) or {}
        st.update({"status": "failed", "post_interpro_status": "failed",
                   "failed_reason": reason, "error": reason,
                   "next_action": "run_cluster_roundtrip_command",
                   "current_step": "failed", "last_updated": now_iso()})
        write_json(run_dir / "status.json", st)
        print(f"FAILED  returned-sequence coverage: {reason}")
        return 6
    if interpro_qc["status"] != "valid":
        st = read_json(run_dir / "status.json", {}) or {}
        reason = interpro_qc["reason"]
        st.update({"status": interpro_qc["status"],
                   "post_interpro_status": interpro_qc["status"],
                   "failed_reason": reason, "error": reason,
                   "next_action": "run_cluster_roundtrip_command",
                   "current_step": interpro_qc["status"],
                   "last_updated": now_iso()})
        write_json(run_dir / "status.json", st)
        print(f"INCOMPLETE  {interpro_qc['status']}: {reason}")
        print("            Not writing fabricated domain/TM/boundary data.")
        return 4
    domain_rows = interpro_qc["rows"]
    tm_rows = pytmhmm_qc["rows"] if pytmhmm_qc["status"] == "valid" else []
    # Attribute every domain/TM row to the species that actually owns its protein.
    # Cluster tools echo only the first FASTA header token (the bare accession, e.g.
    # 'NP_990841.2'), so the parser cannot recover the species from a pipe; we resolve
    # it here from the per-species core tables via an authoritative protein->species map.
    prot_species = _protein_species_map(core_dir)
    # No first-species fallback. Attributing an unresolvable protein to whichever
    # species happens to come first produces a complete-looking table in which one
    # species carries another's domains — the same class of silent, plausible-looking
    # wrong answer as the alphabetical primary guess. An unresolvable accession is a
    # real defect in the identifier chain and must stop the phase.
    unresolved: List[str] = []
    for row in domain_rows + tm_rows:
        pid = row.get("protein_id", "")
        resolved = (row.get("species_id")
                    or prot_species.get(pid)
                    or prot_species.get(_norm_acc(pid)))
        if not resolved:
            unresolved.append(pid)
            continue
        row["species_id"] = resolved
    if unresolved:
        reason = (f"{len(unresolved)} cluster result row(s) could not be attributed to a "
                  f"species: {sorted(set(unresolved))[:8]}. Known proteins: "
                  f"{sorted(set(prot_species.values()))}. Refusing to assign them to an "
                  f"arbitrary species.")
        st = read_json(run_dir / "status.json", {}) or {}
        st.update({"status": "failed", "post_interpro_status": "failed",
                   "failed_reason": reason, "error": reason,
                   "current_step": "failed", "last_updated": now_iso()})
        write_json(run_dir / "status.json", st)
        print(f"FAILED  identifier normalisation: {reason}")
        return 5

    # Raw-signature provenance layer: every normalised member-database hit.
    write_tsv(core_dir / "interpro_annotations.tsv",
              ["species_id", "protein_id", "layer", "member_database",
               "signature_accession", "signature_name", "interpro_accession",
               "interpro_name", "interpro_type", "start_aa", "end_aa",
               "score_or_evalue", "is_integrated"], domain_rows)

    # Curated display layers (representative domains, families, features). This is
    # the scientifically interpreted layer consumed by Domain Architecture and the
    # boundary analysis; raw signatures stay in interpro_annotations.tsv above.
    curated_rows = _curated_annotation_rows(domain_rows)
    write_tsv(core_dir / "domain_features.tsv",
              ["species_id", "protein_id", "layer", "interpro_accession", "interpro_name",
               "interpro_type", "start_aa", "end_aa", "member_databases",
               "supporting_interpro", "n_signatures", "representative_signature",
               "score_or_evalue", "domain_source", "domain_id", "domain_name",
               "feature_type", "interpro_description"], curated_rows)
    write_tsv(core_dir / "tm_features.tsv",
              ["species_id", "protein_id", "start_aa", "end_aa", "source", "topology"], tm_rows)

    # Internal coding-exon boundary distances, scoped to the primary proteins that
    # were actually submitted to InterProScan/pyTMHMM (one per species). Boundaries
    # are measured ONLY against the representative structural domain layer (type
    # DOMAIN/REPEAT); families, superfamilies, sites and disorder are never used as
    # domain edges. The terminal boundary (C-terminus = end of the last coding exon)
    # is excluded because it is not an internal junction.
    exon_map = list(csv.DictReader(open(core_dir / "exon_protein_map.tsv", encoding="utf-8"),
                                   delimiter="\t")) if (core_dir / "exon_protein_map.tsv").is_file() else []
    primary_set = set(primary_ids) | {_norm_acc(p) for p in primary_ids}
    representative_domain_rows = [d for d in curated_rows if d.get("layer") == "domain"]
    dom_by_prot: Dict[str, List[Dict[str, Any]]] = {}
    for d in representative_domain_rows:
        dpid = d.get("protein_id", "")
        dom_by_prot.setdefault(dpid, []).append(d)
        dom_by_prot.setdefault(_norm_acc(dpid), []).append(d)
    exons_by_prot: Dict[str, List[Dict[str, Any]]] = {}
    for e in exon_map:
        pid = e.get("protein_id", "")
        if pid in primary_set or _norm_acc(pid) in primary_set:
            exons_by_prot.setdefault(pid, []).append(e)
    boundary_rows: List[Dict[str, Any]] = []
    for pid, exons in exons_by_prot.items():
        exons_sorted = sorted(exons, key=lambda e: (_to_int(e.get("protein_end_aa")) or 0))
        # internal junctions only: drop the final exon (protein C-terminus)
        for e in exons_sorted[:-1]:
            end = _to_int(e.get("protein_end_aa"))
            if end is None:
                continue
            domains = dom_by_prot.get(pid) or dom_by_prot.get(_norm_acc(pid), [])
            nearest, btype, cat, dist = _classify_boundary(end, domains, NEAR_BOUNDARY_AA)
            edge_coord = (nearest or {}).get("start_aa" if btype == "domain_start" else "end_aa")
            signed_dist = end - edge_coord if edge_coord is not None else None
            nearest_acc = (nearest or {}).get("interpro_accession") or (nearest or {}).get("domain_id", "")
            # An InterPro accession cannot identify a feature when the entry is
            # repeated (FGFR1 carries three IPR007110 Ig-like domains), so persist the
            # concrete instance — accession plus its own coordinates — alongside it.
            nearest_start = (nearest or {}).get("start_aa")
            nearest_end = (nearest or {}).get("end_aa")
            instance_id = (f"{nearest_acc}:{nearest_start}-{nearest_end}"
                           if nearest and nearest_start is not None else "")
            boundary_rows.append({
                "analysis_id": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
                "species_id": e.get("species_id", ""), "protein_id": pid,
                "transcript_id": e.get("transcript_id", ""),
                "exon_boundary_id": f"{pid}:cds{e.get('exon_number','')}_end",
                "boundary_position_aa": end,
                "nearest_domain_accession": nearest_acc,
                "nearest_domain_instance_id": instance_id,
                "nearest_domain_start_aa": nearest_start,
                "nearest_domain_end_aa": nearest_end,
                "nearest_domain_id": nearest_acc,
                "nearest_domain_name": (nearest or {}).get("domain_name")
                    or (nearest or {}).get("interpro_name", ""),
                "nearest_domain_type": (nearest or {}).get("interpro_type", ""),
                "nearest_edge": btype.replace("domain_", "") if nearest else "",
                "nearest_domain_boundary_type": btype,
                "domain_edge_type": btype.replace("domain_", ""),
                "signed_distance_aa": signed_dist, "absolute_distance_aa": dist,
                "distance_aa": dist,
                "classification": cat, "category": cat,
                "domain_layer": "representative_domain",
                "source": "core_post_interpro",
            })
    write_tsv(core_dir / "exon_domain_boundary_distances.tsv",
              ["analysis_id", "gene_symbol", "species_id", "protein_id", "transcript_id",
               "exon_boundary_id", "boundary_position_aa", "nearest_domain_accession",
               "nearest_domain_instance_id", "nearest_domain_start_aa", "nearest_domain_end_aa",
               "nearest_domain_id", "nearest_domain_name", "nearest_domain_type",
               "nearest_edge", "nearest_domain_boundary_type", "domain_edge_type",
               "signed_distance_aa", "absolute_distance_aa", "distance_aa",
               "classification", "category", "domain_layer", "source"], boundary_rows)
    summary_counts: Dict[str, int] = {c: 0 for c in
                                      ("exact_edge", "near_edge", "inside_domain",
                                       "outside_domain", "unknown")}
    for row in boundary_rows:
        summary_counts[row["category"]] = summary_counts.get(row["category"], 0) + 1
    summary_rows = [{"category": key, "count": value, "near_edge_threshold_aa": NEAR_BOUNDARY_AA}
                    for key, value in sorted(summary_counts.items())]
    write_tsv(core_dir / "exon_domain_boundary_summary.tsv",
              ["category", "count", "near_edge_threshold_aa"], summary_rows)
    boundary_species = sorted({r["species_id"] for r in boundary_rows if r.get("species_id")})
    boundary_proteins = sorted({r["protein_id"] for r in boundary_rows if r.get("protein_id")})
    report_path = run_dir / "results" / "16_final_analyses" / "boundary_analysis_report.json"
    write_json(report_path, {
        "status": "available", "mode": "generic_internal_exon_boundaries",
        "gene_symbol": cfg.gene_symbol,
        "protein_scope": "primary_only",
        "isoform_scope": "primary_only",
        "boundary_scope": "internal_coding_exon_boundaries",
        "selected_primary_proteins": primary_ids,
        "species_scope": boundary_species,
        "n_species": len(boundary_species),
        "n_proteins": len(boundary_proteins),
        "domain_layer": "representative_domain",
        "near_edge_threshold_aa": NEAR_BOUNDARY_AA,
        "n_boundaries": len(boundary_rows), "category_counts": summary_counts,
        "source": "real fetched InterProScan coordinates", "generated_at": now_iso(),
    })
    candidate_rows = []
    clusters_path = run_dir / "results" / "generic_gene_analysis" / "event_region_candidate_clusters.tsv"
    clusters = list(csv.DictReader(open(clusters_path, encoding="utf-8"), delimiter="\t")) \
        if clusters_path.is_file() else []
    primary_domains = [d for d in representative_domain_rows if d.get("protein_id") in primary_ids]
    for c in clusters:
        cid = c.get("candidate_cluster_id", "")
        start = _to_int(c.get("representative_start_aa"))
        end_pos = _to_int(c.get("representative_end_aa"))
        if start is None or end_pos is None:
            continue
        overlaps = [d for d in primary_domains
                    if d.get("start_aa") is not None and d.get("end_aa") is not None
                    and start <= d["end_aa"] and end_pos >= d["start_aa"]]
        edges = [(abs(start - d[edge]), d, edge) for d in primary_domains
                 for edge in ("start_aa", "end_aa") if d.get(edge) is not None]
        nearest_dist, nearest, edge = min(edges, default=(None, {}, ""), key=lambda x: x[0])
        relation = "inside_or_overlapping" if overlaps else (
            "near_boundary" if nearest_dist is not None and nearest_dist <= NEAR_BOUNDARY_AA
            else "outside")
        candidate_rows.append({
            "candidate_id": cid, "reference_protein": primary_ids[0] if primary_ids else "",
            "aa_start": start, "aa_end": end_pos, "status": "real_domain_context",
            "overlapping_domain_ids": ";".join(d.get("domain_id", "") for d in overlaps),
            "overlapping_domain_names": ";".join(d.get("domain_name", "") for d in overlaps),
            "nearest_domain_id": nearest.get("domain_id", ""),
            "nearest_domain_edge": edge.replace("_aa", ""),
            "distance_to_nearest_edge_aa": nearest_dist,
            "relation": relation, "domain_source": nearest.get("domain_source", ""),
            "signature_accession": nearest.get("domain_id", ""),
            "evidence_confidence": "high" if overlaps else "medium",
            "domain_context_score": 10 if overlaps else 5 if relation == "near_boundary" else 0,
        })
    generic_dir = run_dir / "results" / "generic_gene_analysis"
    candidate_cols = list(candidate_rows[0]) if candidate_rows else [
        "candidate_id", "reference_protein", "aa_start", "aa_end", "status"]
    write_tsv(generic_dir / "candidate_domain_context.tsv", candidate_cols, candidate_rows)

    # refresh the core report counters
    report = read_json(core_dir / "core_gene_report.json", {}) or {}
    report.setdefault("summary", {})
    report["summary"]["n_domain_features"] = len(representative_domain_rows)
    report["summary"]["n_representative_domains"] = len(representative_domain_rows)
    report["summary"]["n_interpro_annotations_raw"] = len(domain_rows)
    report["summary"]["n_tm_features"] = len(tm_rows)
    report["summary"]["n_exon_boundaries"] = len(boundary_rows)
    report["domain_status"] = "complete"
    report["tm_status"] = pytmhmm_qc["status"]
    report["post_cluster_qc"] = str(qc_path.relative_to(run_dir))
    report["post_interpro_generated_at"] = now_iso()
    write_json(core_dir / "core_gene_report.json", report)

    core_src = CoreSource(core_dir, f"run:{run_dir.name}", cfg)
    build_core_indices(core_src, run_dir / "website_indices" / "generic")

    # refresh the shared generic layer now that domain/boundary are available
    _run_generic_pipeline(run_dir.name, lambda m: print(f"  {m}"))

    # Rebuild the coordinate model from the final domain, TM and boundary tables.
    _build_coordinate_model_and_figures(run_dir, lambda m: print(f"  {m}"))

    partial_reasons = []
    if pytmhmm_qc["status"] != "valid":
        partial_reasons.append(pytmhmm_qc["reason"])

    # Completeness is a per-species question. Deriving it from stage-file presence let
    # the two-species FGFR1 run advertise "Results ready" while one species had no
    # domain and no boundary results: the file existed and covered both species, but one
    # of them was empty. The dataset selector reads this status, so it must never claim
    # more than every species actually reached.
    species_status: Dict[str, Any] = {}
    try:
        from exondomaincompare.framework.species_completion import (
            aggregate_run_status, build_species_completion, species_status_summary,
        )
        completion = build_species_completion(run_dir)
        agg_status, agg_reasons = aggregate_run_status(
            completion, cluster_complete=True, failed_reasons=list(partial_reasons))
        species_status = species_status_summary(completion)
        final_status = agg_status
        partial_reasons = agg_reasons or partial_reasons
    except Exception as comp_err:  # pragma: no cover - aggregation is additive
        print(f"  WARN: per-species completion aggregation skipped: {comp_err}")
        final_status = "results_partial" if partial_reasons else "results_ready"

    st = read_json(run_dir / "status.json", {}) or {}
    st.update({"status": final_status, "post_interpro_status": final_status,
               "cluster_fetch_status": "complete", "next_action": "open_results",
               "failed_reason": "; ".join(partial_reasons), "error": "",
               "current_step": final_status, "last_updated": now_iso()})
    if species_status:
        st["species_status"] = species_status
        st["run_status"] = final_status
    write_json(run_dir / "status.json", st)
    evidence = build_coordinate_evidence_register(run_dir)
    print(f"OK  coordinate evidence register: "
          f"{evidence['counts']['total_records']} records ({evidence['register_phase']})")
    print(f"OK  core post-analysis: domains={len(domain_rows)} tm={len(tm_rows)} "
          f"boundaries={len(boundary_rows)}")
    return 0


# --------------------------------------------------------------------------- #
# Real cluster-output parsers (gene-agnostic; standard InterProScan / pyTMHMM)
# --------------------------------------------------------------------------- #
def _returned_sequence_inventory(core_dir: Path, interpro_qc: Dict[str, Any],
                                 pytmhmm_qc: Dict[str, Any]) -> Dict[str, Any]:
    try:
        sys.path.insert(0, str(SCRIPTS_DIR))
        from exondomaincompare.framework.primary_resolution import resolve_primaries
        primaries = resolve_primaries(core_dir)
    except Exception:
        primaries = {}

    ips_rows = interpro_qc.get("rows") or []
    tm_rows = pytmhmm_qc.get("rows") or []

    def _counts(rows: List[Dict[str, Any]], pid: str) -> int:
        norm = _norm_acc(pid)
        return sum(1 for r in rows
                   if r.get("protein_id") == pid or _norm_acc(r.get("protein_id", "")) == norm)

    species_rows: List[Dict[str, Any]] = []
    for sid, info in sorted(primaries.items()):
        pid = info.get("protein_id") or ""
        n_ips, n_tm = _counts(ips_rows, pid), _counts(tm_rows, pid)
        species_rows.append({
            "species_id": sid,
            "protein_id": pid,
            "transcript_id": info.get("transcript_id") or "",
            "resolved_from": info.get("resolved_from") or "",
            "interproscan_features": n_ips,
            "pytmhmm_features": n_tm,
            "interproscan_status": ("available" if n_ips else
                                    "returned_no_features" if ips_rows else
                                    "no_returned_sequence"),
            "pytmhmm_status": "available" if n_tm else "no_transmembrane_predicted",
        })
    return {
        "n_species_submitted": len(species_rows),
        "n_species_with_features": sum(1 for r in species_rows
                                       if r["interproscan_status"] == "available"),
        "species": species_rows,
    }


def _primary_protein_ids(core_dir: Path) -> List[str]:
    ids: List[str] = []
    iso = core_dir / "protein_isoform_index.tsv"
    if iso.is_file():
        for r in csv.DictReader(open(iso, encoding="utf-8"), delimiter="\t"):
            if r.get("protein_id") and str(r.get("primary_status", "")).lower() == "primary":
                ids.append(r["protein_id"])
    return ids


def _resolve_seq_identity(raw_acc: str, primary_ids: List[str]) -> Tuple[str, str]:
    import re
    acc = (raw_acc or "").strip()
    if "|" in acc:
        parts = acc.split("|")
        species = parts[0]
        # pick the token that matches a known primary protein id, else a token that
        # looks like a protein accession (RefSeq XP_/NP_ or Ensembl ...P...).
        pid = next((p for p in parts if p in primary_ids), "")
        if not pid:
            pid = next((p for p in parts
                        if re.match(r"^(XP_|NP_|AP_|YP_)\d", p) or re.match(r"^ENS[A-Z]*P\d", p)), "")
        if not pid:
            pid = parts[3] if len(parts) > 3 else parts[-1]
        return species, pid
    return "", acc


def _parse_interproscan(ips_out: Path, primary_ids: List[str]) -> List[Dict[str, Any]]:
    if not ips_out.is_dir():
        return []
    hits = load_normalized_annotations(ips_out)
    rows: List[Dict[str, Any]] = []
    for h in hits:
        rows.append({
            "species_id": "", "protein_id": h["protein_accession"],
            "member_database": h["member_database"],
            "signature_accession": h["signature_accession"],
            "signature_name": h["signature_name"],
            "interpro_accession": h["interpro_accession"],
            "interpro_name": h["interpro_name"],
            "interpro_type": h["interpro_type"],
            "start_aa": h["start"], "end_aa": h["end"],
            "score_or_evalue": h["score_or_evalue"],
            "is_integrated": "1" if h["is_integrated"] else "0",
            "layer": h["layer"],
        })
    return rows


def _curated_annotation_rows(raw_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_prot: Dict[str, List[Dict[str, Any]]] = {}
    for r in raw_rows:
        by_prot.setdefault(r["protein_id"], []).append(r)

    def _hit(r: Dict[str, Any]) -> Dict[str, Any]:
        return {"start": r["start_aa"], "end": r["end_aa"],
                "layer": r["layer"], "interpro_accession": r["interpro_accession"],
                "interpro_name": r["interpro_name"], "interpro_type": r["interpro_type"],
                "member_database": r["member_database"],
                "signature_accession": r["signature_accession"],
                "signature_name": r["signature_name"],
                "score_or_evalue": r.get("score_or_evalue", "")}

    curated: List[Dict[str, Any]] = []
    for pid, rows in by_prot.items():
        species = next((r["species_id"] for r in rows if r.get("species_id")), "")
        hits = [_hit(r) for r in rows]

        def _row(layer: str, ipr_acc: str, ipr_name: str, ipr_type: str,
                 start: int, end: int, member_dbs: str, supporting: str,
                 n_sig: int, rep_sig: str, score: str) -> Dict[str, Any]:
            return {
                "species_id": species, "protein_id": pid, "layer": layer,
                "interpro_accession": ipr_acc, "interpro_name": ipr_name,
                "interpro_type": ipr_type, "start_aa": start, "end_aa": end,
                "member_databases": member_dbs, "supporting_interpro": supporting,
                "n_signatures": n_sig, "representative_signature": rep_sig,
                "score_or_evalue": score,
                # legacy columns for backward compatibility
                "domain_source": member_dbs, "domain_id": ipr_acc or rep_sig,
                "domain_name": ipr_name, "feature_type": layer,
                "interpro_description": ipr_name,
            }

        for d in representative_domains(hits):
            supporting = ";".join(f"{s['interpro_accession']}:{s['interpro_name']}"
                                  for s in d["supporting_interpro"])
            curated.append(_row("domain", d["interpro_accession"], d["interpro_name"],
                                d["interpro_type"], d["start_aa"], d["end_aa"],
                                ",".join(d["member_databases"]), supporting,
                                d["n_signatures"], d["representative_signature"],
                                d["score_or_evalue"]))
        for f in family_annotations(hits):
            curated.append(_row("family", f["interpro_accession"], f["interpro_name"],
                                f["interpro_type"], f["start_aa"], f["end_aa"],
                                ",".join(f["member_databases"]), "", 0, "", ""))
        for ft in feature_annotations(hits):
            curated.append(_row("feature", ft["interpro_accession"], ft["interpro_name"],
                                ft["interpro_type"], ft["start_aa"], ft["end_aa"],
                                ft["member_database"], "", 0, "", ""))
    curated.sort(key=lambda r: (r["species_id"], r["protein_id"],
                                {"domain": 0, "family": 1, "feature": 2}.get(r["layer"], 3),
                                r["start_aa"]))
    return curated


def _parse_pytmhmm(tm_out: Path, primary_ids: List[str]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    seen: set = set()
    if not tm_out.is_dir():
        return rows

    def emit(seq_id: str, line: str) -> None:
        species, pid = _resolve_seq_identity(seq_id, primary_ids)
        for s, e in parse_pytmhmm_topology_line(line):
            # pyTMHMM positions are 0-based residue indices; convert to 1-based aa.
            key = (pid or seq_id, s + 1, e + 1)
            if key in seen:
                continue  # both combined tables list the same TM hit; keep one row
            seen.add(key)
            rows.append({"species_id": species, "protein_id": pid or seq_id,
                         "start_aa": s + 1, "end_aa": e + 1, "source": "pytmhmm",
                         "topology": "transmembrane"})

    # Prefer a single combined table (the two contain the same TM hits); only fall
    # back to per-sequence .summary files when no combined table is present.
    combined = tm_out / "pytmhmm_transmembrane_hits.tsv"
    summary_all = tm_out / "pytmhmm_summary_all.tsv"
    source = next((f for f in (combined, summary_all) if f.is_file()), None)
    if source is not None:
        for r in csv.DictReader(open(source, encoding="utf-8"), delimiter="\t"):
            sid, line = r.get("sequence_id", ""), r.get("line", "")
            if sid and line:
                emit(sid, line)
    else:
        for summ in tm_out.glob("*.summary"):
            seq_id = summ.name.rsplit(".summary", 1)[0]
            for line in summ.read_text(encoding="utf-8", errors="replace").splitlines():
                emit(seq_id, line)
    return rows


def _interpro_parse_qc(ips_out: Path, primary_ids: List[str]) -> Dict[str, Any]:
    files = [p for p in ips_out.rglob("*.tsv") if p.is_file() and p.stat().st_size > 0] \
        if ips_out.is_dir() else []
    if not files:
        return {"status": "cluster_outputs_missing", "rows": [], "files": [],
                "reason": f"No non-empty InterProScan TSV found under {ips_out}"}
    valid_lines = malformed_lines = 0
    for path in files:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 13 or _to_int(cols[6]) is None or _to_int(cols[7]) is None:
                malformed_lines += 1
            else:
                valid_lines += 1
    if not valid_lines:
        return {"status": "interpro_parse_failed", "rows": [], "files": [str(p) for p in files],
                "reason": f"No valid InterProScan rows; malformed rows={malformed_lines}"}
    rows = _parse_interproscan(ips_out, primary_ids)
    if not rows:
        return {"status": "interpro_parse_failed", "rows": [], "files": [str(p) for p in files],
                "reason": "Valid-looking TSV rows did not map to submitted protein accessions."}
    return {"status": "valid", "rows": rows, "files": [str(p) for p in files],
            "valid_lines": valid_lines, "malformed_lines": malformed_lines, "reason": ""}


def _pytmhmm_parse_qc(tm_out: Path, primary_ids: List[str]) -> Dict[str, Any]:
    files = []
    if tm_out.is_dir():
        files = [p for p in tm_out.rglob("*") if p.is_file() and p.stat().st_size > 0
                 and (p.suffix in (".tsv", ".summary") or "pytmhmm" in p.name.lower())]
    if not files:
        return {"status": "cluster_outputs_missing", "rows": [], "files": [],
                "reason": f"No non-empty pyTMHMM output found under {tm_out}"}
    try:
        rows = _parse_pytmhmm(tm_out, primary_ids)
    except Exception as exc:
        return {"status": "pytmhmm_parse_failed", "rows": [], "files": [str(p) for p in files],
                "reason": str(exc)}
    # A syntactically valid topology output may contain no transmembrane segment.
    return {"status": "valid", "rows": rows, "files": [str(p) for p in files],
            "reason": "" if rows else "Valid pyTMHMM output contains no TM segment."}


def _to_int(v: Any) -> Optional[int]:
    try:
        if v is None or v == "":
            return None
        return int(float(v))
    except (TypeError, ValueError):
        return None


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Experimental core-only live runner (no event region).")
    ap.add_argument("--gene-config", help="Path to a no-event gene config (core_only_pilot).")
    ap.add_argument("--gene", help="Gene symbol; resolves a drafts/<GENE>_core_only_pilot.yaml.")
    ap.add_argument("--species", nargs="+",
                    help="One or more species (free text or canonical id), e.g. "
                         "--species 'Gallus gallus' 'Mus musculus'. Multiple species build a "
                         "cross-species (one primary per species) analysis.")
    ap.add_argument("--run-name",
                    help="Optional human-readable run name. The run directory is "
                         "named from the gene and species regardless.")
    ap.add_argument("--support-level", default="core_only_pilot",
                    help="Informational; must be core_only_pilot.")
    ap.add_argument("--selection-method", default="canonical_if_available",
                    choices=["longest_protein", "canonical_if_available", "all_protein_coding_transcripts"],
                    help="Primary protein selection rule (documented in the collection report).")
    ap.add_argument("--input-mode", default="auto",
                    choices=["auto", "local_cache", "user_files"],
                    help="auto: local cache now, download later; local_cache: cache only; "
                         "user_files: explicit --gff/--protein-faa (future primary path).")
    ap.add_argument("--gff", help="Override: explicit genomic.gff path (whole genome; user_files mode).")
    ap.add_argument("--protein-faa", help="Override: explicit protein.faa path (user_files mode).")
    ap.add_argument("--post", action="store_true", help="Post phase: build domains/boundaries from cluster outputs.")
    ap.add_argument("--run-id", help="Run id (required with --post).")
    ap.add_argument("--reuse-run-id",
                    help="Retry this existing run in place instead of scaffolding a new "
                         "one. Keeps the run id, name, gene and species, and rebuilds the "
                         "derived stages that failed. Use after an upstream repair so the "
                         "user is not left with two runs for one request.")
    ap.add_argument("--config")
    ap.add_argument("--local-profile")
    ap.add_argument("--lrz-profile")
    args = ap.parse_args(argv)

    global RUNTIME_CONFIG, PROJECT_ROOT, SCRIPTS_DIR, RUNS_ROOT, FREEZE_ROOT
    RUNTIME_CONFIG = load_config(
        config_path=args.config, repository_root=PROJECT_ROOT,
        local_profile=args.local_profile, lrz_profile=args.lrz_profile,
    )
    PROJECT_ROOT = RUNTIME_CONFIG.repository_root
    _SCRIPTS_DIR = PROJECT_ROOT / "scripts"
    _RUNS_ROOT = RUNTIME_CONFIG.runs_root
    _FREEZE_ROOT = PROJECT_ROOT / "results" / "final_30_until_interpro_prepare"

    if args.post:
        if not args.run_id:
            raise SystemExit("--post requires --run-id")
        return phase_post(args)
    if not args.species:
        raise SystemExit("--species is required for the create phase.")
    return phase_create(args)


if __name__ == "__main__":
    raise SystemExit(main())
