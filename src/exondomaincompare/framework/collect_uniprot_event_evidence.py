#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from exondomaincompare.config import discover_repository_root, load_config
from exondomaincompare.runs.registry import RegistryError, resolve_run_record

PROJECT_ROOT = discover_repository_root(__file__)
RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)

EVIDENCE_FILE = "event_region_evidence.tsv"
REPORT_FILE = "uniprot_event_evidence_report.json"

# Must match build_event_region_evidence.py.
EVIDENCE_COLUMNS = [
    "analysis_id", "gene_symbol", "species_id", "event_candidate_id",
    "evidence_source", "evidence_status",
    "transcript_a", "transcript_b", "protein_a", "protein_b",
    "region_start_aa", "region_end_aa", "region_length_aa",
    "event_type_candidate", "exon_aligned",
    "confidence", "confidence_reason", "notes",
]

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(p: Path, default: Any = None) -> Any:
    try:
        return json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:
        return default


def write_json(p: Path, data: Any) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def read_tsv(p: Path) -> List[Dict[str, str]]:
    if not Path(p).is_file():
        return []
    with open(p, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def append_tsv(path: Path, columns: List[str], rows: List[Dict[str, Any]]) -> None:
    exists = path.is_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        if not exists:
            w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in columns})


def _http_get_json(url: str, timeout: float) -> Optional[Dict[str, Any]]:
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "ExonDomainCompare/1.0 (core-only pilot)"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310 (fixed host)
        return json.loads(resp.read().decode("utf-8"))


def query_uniprot(gene_symbol: str, taxid: str, timeout: float) -> Dict[str, Any]:
    if not gene_symbol:
        return {"ok": False, "entries": [], "error": "no_gene_symbol"}
    terms = [f'gene:{gene_symbol}']
    if taxid:
        terms.append(f'taxonomy_id:{taxid}')
    query = " AND ".join(terms)
    # Request full entries (no field restriction) so the response includes the
    # `features` (VAR_SEQ) and `comments` (alternative products) blocks we parse.
    params = {
        "query": query,
        "format": "json",
        "size": "10",
    }
    url = f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"
    try:
        data = _http_get_json(url, timeout=timeout)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
        return {"ok": False, "entries": [], "error": f"network_error:{exc}"}
    except Exception as exc:  # noqa: BLE001 - never fail the run
        return {"ok": False, "entries": [], "error": f"unexpected_error:{exc}"}
    return {"ok": True, "entries": data.get("results", []) or [], "error": ""}


def _extract_var_seq_features(entry: Dict[str, Any]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for feat in entry.get("features", []) or []:
        ftype = str(feat.get("type", "")).lower()
        if "alternative sequence" not in ftype and ftype != "var_seq":
            continue
        loc = feat.get("location", {}) or {}
        start = ((loc.get("start") or {}).get("value"))
        end = ((loc.get("end") or {}).get("value"))
        out.append({
            "start_aa": start, "end_aa": end,
            "description": feat.get("description", "") or "alternative sequence",
        })
    return out


def _has_alternative_products(entry: Dict[str, Any]) -> bool:
    for c in entry.get("comments", []) or []:
        if str(c.get("commentType", "")).lower() == "alternative products":
            return True
    return False


def collect(run_dir: Path, timeout: float, offline: bool) -> Dict[str, Any]:
    core_dir = run_dir / "results" / "core_gene_analysis"
    rc = read_json(run_dir / "run_config.json", {}) or {}
    report = read_json(core_dir / "core_gene_report.json", {}) or {}
    gene_symbol = rc.get("gene_symbol", "") or report.get("gene_symbol", "")
    analysis_id = rc.get("analysis_id", "") or report.get("analysis_id", "")
    prov = rc.get("annotation_provenance", {}) or {}
    taxid = str(prov.get("taxid", "") or "").strip()
    iso = read_tsv(core_dir / "protein_isoform_index.tsv")
    species_id = iso[0].get("species_id", "") if iso else ""
    accessions = [r.get("protein_id", "") for r in iso if r.get("protein_id")]

    base_report: Dict[str, Any] = {
        "generated_at": now_iso(),
        "gene_symbol": gene_symbol,
        "analysis_id": analysis_id,
        "taxid": taxid,
        "species_id": species_id,
        "n_local_accessions": len(accessions),
        "source": "uniprotkb_rest",
    }

    if offline:
        base_report.update({"status": "uniprot_evidence_unavailable",
                            "reason": "offline_mode_requested",
                            "n_curated_rows_appended": 0})
        write_json(core_dir / REPORT_FILE, base_report)
        return base_report

    q = query_uniprot(gene_symbol, taxid, timeout)
    if not q["ok"]:
        base_report.update({"status": "uniprot_evidence_unavailable",
                            "reason": q["error"] or "uniprot_unavailable",
                            "n_curated_rows_appended": 0})
        write_json(core_dir / REPORT_FILE, base_report)
        return base_report

    appended: List[Dict[str, Any]] = []
    entries_seen = 0
    for entry in q["entries"]:
        entries_seen += 1
        acc = entry.get("primaryAccession", "")
        var_seqs = _extract_var_seq_features(entry)
        has_alt = _has_alternative_products(entry)
        if not var_seqs and not has_alt:
            continue
        if not var_seqs:
            # Alternative products documented but no residue-level feature.
            appended.append({
                "analysis_id": analysis_id, "gene_symbol": gene_symbol,
                "species_id": species_id,
                "event_candidate_id": f"uniprot:{acc}:alt_products",
                "evidence_source": "uniprot_alternative_sequence",
                "evidence_status": "curated_annotation",
                "transcript_a": "", "transcript_b": "",
                "protein_a": acc, "protein_b": "",
                "region_start_aa": "", "region_end_aa": "", "region_length_aa": "",
                "event_type_candidate": "alternative_products",
                "exon_aligned": "no",
                "confidence": "curated",
                "confidence_reason": "UniProt documents alternative products for this entry.",
                "notes": "Curated UniProt annotation (no residue-level VAR_SEQ). "
                         "Evidence only; not a validated event region.",
            })
            continue
        for vs in var_seqs:
            s, e = vs.get("start_aa"), vs.get("end_aa")
            length = (e - s + 1) if (isinstance(s, int) and isinstance(e, int)) else ""
            appended.append({
                "analysis_id": analysis_id, "gene_symbol": gene_symbol,
                "species_id": species_id,
                "event_candidate_id": f"uniprot:{acc}:{s}-{e}",
                "evidence_source": "uniprot_alternative_sequence",
                "evidence_status": "curated_annotation",
                "transcript_a": "", "transcript_b": "",
                "protein_a": acc, "protein_b": "",
                "region_start_aa": s if s is not None else "",
                "region_end_aa": e if e is not None else "",
                "region_length_aa": length,
                "event_type_candidate": "alternative_sequence",
                "exon_aligned": "no",
                "confidence": "curated",
                "confidence_reason": f"UniProt VAR_SEQ feature: {vs.get('description','')}.",
                "notes": "Curated UniProt alternative-sequence feature. "
                         "Evidence only; not a validated event region.",
            })

    if appended:
        append_tsv(core_dir / EVIDENCE_FILE, EVIDENCE_COLUMNS, appended)
        base_report.update({"status": "uniprot_evidence_appended",
                            "entries_matched": entries_seen,
                            "n_curated_rows_appended": len(appended)})
    else:
        base_report.update({"status": "uniprot_evidence_unavailable",
                            "reason": "no_alternative_sequence_evidence_found",
                            "entries_matched": entries_seen,
                            "n_curated_rows_appended": 0})
    write_json(core_dir / REPORT_FILE, base_report)
    return base_report


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Optional, non-blocking UniProt curated event-evidence collector.")
    ap.add_argument("--run-id", required=True, help="Run id under runs/.")
    ap.add_argument("--timeout", type=float, default=8.0, help="HTTP timeout seconds (default 8).")
    ap.add_argument("--offline", action="store_true",
                    help="Skip all network access; write uniprot_evidence_unavailable.")
    args = ap.parse_args(argv)

    try:
        record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    except RegistryError as exc:
        print(f"WARN {exc} (uniprot step skipped)")
        return 0
    if record is None:
        # Even a bad run id must not raise for the caller; report and succeed.
        print(f"WARN run not found: {args.run_id} (uniprot step skipped)")
        return 0
    if record.read_only:
        print("WARN run is registered read-only (uniprot step skipped)")
        return 0
    run_dir = record.path

    try:
        rep = collect(run_dir, args.timeout, args.offline)
    except Exception as exc:  # noqa: BLE001 - never fail the run
        print(f"WARN UniProt collector error (non-blocking): {exc}")
        try:
            write_json(run_dir / "results" / "core_gene_analysis" / REPORT_FILE, {
                "generated_at": now_iso(), "status": "uniprot_evidence_unavailable",
                "reason": f"collector_error:{exc}", "n_curated_rows_appended": 0,
            })
        except Exception:
            pass
        return 0

    print(f"UniProt evidence: status={rep.get('status')} "
          f"appended={rep.get('n_curated_rows_appended', 0)}")
    if rep.get("status") == "uniprot_evidence_unavailable":
        print(f"    reason: {rep.get('reason', 'unknown')} (non-blocking; core run unaffected)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
