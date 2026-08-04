#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

# InterPro entry types -> display layer (the only place this mapping lives).
DOMAIN_ENTRY_TYPES = {"DOMAIN", "REPEAT"}
FAMILY_ENTRY_TYPES = {"FAMILY", "HOMOLOGOUS_SUPERFAMILY"}
FEATURE_ENTRY_TYPES = {"ACTIVE_SITE", "BINDING_SITE", "CONSERVED_SITE", "SITE", "PTM"}
# member databases whose (usually unintegrated) hits are short features, not domains.
FEATURE_MEMBER_DBS = {"MOBIDB_LITE", "MOBIDBLITE", "COILS", "PHOBIUS", "SIGNALP", "TMHMM"}


def layer_for(interpro_type: str, is_integrated: bool, member_database: str = "") -> str:
    t = (interpro_type or "").strip().upper()
    if is_integrated and t in DOMAIN_ENTRY_TYPES:
        return "domain"
    if is_integrated and t in FAMILY_ENTRY_TYPES:
        return "family"
    if is_integrated and t in FEATURE_ENTRY_TYPES:
        return "feature"
    if (member_database or "").strip().upper() in FEATURE_MEMBER_DBS:
        return "feature"
    return "raw"


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #
def _num_str(v: Any) -> str:
    if v is None:
        return ""
    return str(v)


def parse_interproscan_json(path: Path,
                            species_resolver: Optional[Callable[[str], str]] = None
                            ) -> List[Dict[str, Any]]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    hits: List[Dict[str, Any]] = []
    for result in data.get("results", []):
        xrefs = result.get("xref") or [{}]
        protein_acc = (xrefs[0].get("id") or "").strip()
        species = species_resolver(protein_acc) if species_resolver else ""
        for match in result.get("matches", []):
            sig = match.get("signature") or {}
            lib = ((sig.get("signatureLibraryRelease") or {}).get("library") or "").strip()
            entry = sig.get("entry") or None
            ipr_acc = (entry or {}).get("accession", "") if entry else ""
            ipr_name = (entry or {}).get("name", "") if entry else ""
            ipr_type = (entry or {}).get("type", "") if entry else ""
            is_integrated = bool(ipr_acc)
            score = match.get("evalue", match.get("score"))
            for loc in match.get("locations", []) or []:
                start, end = loc.get("start"), loc.get("end")
                if start is None or end is None:
                    continue
                hits.append({
                    "protein_accession": protein_acc,
                    "species": species,
                    "signature_accession": (sig.get("accession") or "").strip(),
                    "signature_name": (sig.get("name") or sig.get("description") or "").strip(),
                    "member_database": lib,
                    "interpro_accession": ipr_acc,
                    "interpro_name": ipr_name,
                    "interpro_type": ipr_type,
                    "start": int(start),
                    "end": int(end),
                    "score_or_evalue": _num_str(loc.get("evalue", score)),
                    "is_integrated": is_integrated,
                    "layer": layer_for(ipr_type, is_integrated, lib),
                })
    return hits


def parse_interproscan_tsv(path: Path,
                           species_resolver: Optional[Callable[[str], str]] = None,
                           type_lookup: Optional[Dict[str, str]] = None
                           ) -> List[Dict[str, Any]]:
    hits: List[Dict[str, Any]] = []
    type_lookup = type_lookup or {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            try:
                start, end = int(float(f[6])), int(float(f[7]))
            except (ValueError, TypeError):
                continue
            protein_acc = f[0].strip()
            ipr_acc = f[11] if len(f) > 11 and f[11] != "-" else ""
            ipr_name = f[12] if len(f) > 12 and f[12] != "-" else ""
            ipr_type = type_lookup.get(ipr_acc, "")
            is_integrated = bool(ipr_acc)
            hits.append({
                "protein_accession": protein_acc,
                "species": species_resolver(protein_acc) if species_resolver else "",
                "signature_accession": f[4],
                "signature_name": f[5],
                "member_database": f[3],
                "interpro_accession": ipr_acc,
                "interpro_name": ipr_name,
                "interpro_type": ipr_type,
                "start": start, "end": end,
                "score_or_evalue": f[8],
                "is_integrated": is_integrated,
                "layer": layer_for(ipr_type, is_integrated, f[3]),
            })
    return hits


def load_normalized_annotations(ips_out: Path,
                                species_resolver: Optional[Callable[[str], str]] = None
                                ) -> List[Dict[str, Any]]:
    ips_out = Path(ips_out)
    jsons = sorted(p for p in ips_out.rglob("*.json")
                   if p.is_file() and p.stat().st_size > 0)
    if jsons:
        hits: List[Dict[str, Any]] = []
        for j in jsons:
            try:
                hits.extend(parse_interproscan_json(j, species_resolver))
            except (json.JSONDecodeError, KeyError):
                continue
        if hits:
            return hits
    # TSV fallback (harvest entry types from JSON if we have any)
    type_lookup: Dict[str, str] = {}
    for j in jsons:
        try:
            for h in parse_interproscan_json(j):
                if h["interpro_accession"]:
                    type_lookup[h["interpro_accession"]] = h["interpro_type"]
        except (json.JSONDecodeError, KeyError):
            continue
    hits = []
    for tsv in sorted(p for p in ips_out.rglob("*.tsv")
                      if p.is_file() and p.stat().st_size > 0 and "transmembrane" not in p.name):
        hits.extend(parse_interproscan_tsv(tsv, species_resolver, type_lookup))
    return hits


# --------------------------------------------------------------------------- #
# representative domain layer
# --------------------------------------------------------------------------- #
def _cluster_overlaps(intervals: List[Tuple[int, int, Dict[str, Any]]]
                      ) -> List[List[Tuple[int, int, Dict[str, Any]]]]:
    clusters: List[List[Tuple[int, int, Dict[str, Any]]]] = []
    cur_end = None
    for s, e, h in sorted(intervals, key=lambda x: (x[0], x[1])):
        if clusters and cur_end is not None and s <= cur_end:
            clusters[-1].append((s, e, h))
            cur_end = max(cur_end, e)
        else:
            clusters.append([(s, e, h)])
            cur_end = e
    return clusters


def representative_domains(hits: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    dom = [(h["start"], h["end"], h) for h in hits if h.get("layer") == "domain"]
    reps: List[Dict[str, Any]] = []
    for cluster in _cluster_overlaps(dom):
        by_ipr: Dict[str, Dict[str, Any]] = {}
        for s, e, h in cluster:
            acc = h["interpro_accession"]
            node = by_ipr.setdefault(acc, {
                "interpro_accession": acc,
                "interpro_name": h["interpro_name"],
                "interpro_type": h["interpro_type"],
                "member_dbs": set(),
                "signatures": set(),
                "start": s, "end": e,
                "evalues": [],
            })
            node["member_dbs"].add(h["member_database"])
            node["signatures"].add(h["signature_accession"])
            node["start"] = min(node["start"], s)
            node["end"] = max(node["end"], e)
            ev = _to_float(h.get("score_or_evalue"))
            if ev is not None:
                node["evalues"].append(ev)

        def _rank(n: Dict[str, Any]) -> Tuple[int, int, float]:
            best_ev = min(n["evalues"]) if n["evalues"] else float("inf")
            return (len(n["member_dbs"]), n["end"] - n["start"], -best_ev)

        best = max(by_ipr.values(), key=_rank)
        contributing = sorted(
            ({"interpro_accession": n["interpro_accession"],
              "interpro_name": n["interpro_name"],
              "interpro_type": n["interpro_type"]}
             for n in by_ipr.values()),
            key=lambda x: x["interpro_accession"])
        all_member_dbs = sorted({db for n in by_ipr.values() for db in n["member_dbs"] if db})
        all_sigs = sorted({s for n in by_ipr.values() for s in n["signatures"] if s})
        reps.append({
            "interpro_accession": best["interpro_accession"],
            "interpro_name": best["interpro_name"],
            "interpro_type": best["interpro_type"],
            "start_aa": best["start"],
            "end_aa": best["end"],
            "member_databases": all_member_dbs,
            "supporting_interpro": contributing,
            "n_signatures": len(all_sigs),
            "representative_signature": (sorted(best["signatures"])[0]
                                         if best["signatures"] else ""),
            "score_or_evalue": (str(min(best["evalues"])) if best["evalues"] else ""),
        })
    reps.sort(key=lambda d: (d["start_aa"], d["end_aa"]))
    return reps


def family_annotations(hits: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_ipr: Dict[str, Dict[str, Any]] = {}
    for h in hits:
        if h.get("layer") != "family":
            continue
        acc = h["interpro_accession"]
        node = by_ipr.setdefault(acc, {
            "interpro_accession": acc, "interpro_name": h["interpro_name"],
            "interpro_type": h["interpro_type"], "start_aa": h["start"], "end_aa": h["end"],
            "member_databases": set(),
        })
        node["start_aa"] = min(node["start_aa"], h["start"])
        node["end_aa"] = max(node["end_aa"], h["end"])
        node["member_databases"].add(h["member_database"])
    out = []
    for n in by_ipr.values():
        n["member_databases"] = sorted(db for db in n["member_databases"] if db)
        out.append(n)
    out.sort(key=lambda d: (d["start_aa"], d["end_aa"]))
    return out


def feature_annotations(hits: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    seen = set()
    out: List[Dict[str, Any]] = []
    for h in sorted((h for h in hits if h.get("layer") == "feature"),
                    key=lambda x: (x["start"], x["end"])):
        key = (h["interpro_accession"] or h["signature_accession"], h["start"], h["end"])
        if key in seen:
            continue
        seen.add(key)
        out.append({
            "interpro_accession": h["interpro_accession"],
            "interpro_name": h["interpro_name"] or h["signature_name"],
            "interpro_type": h["interpro_type"] or ("DISORDER"
                                                     if h["member_database"].upper() in FEATURE_MEMBER_DBS
                                                     else ""),
            "member_database": h["member_database"],
            "start_aa": h["start"], "end_aa": h["end"],
        })
    return out


def _to_float(v: Any) -> Optional[float]:
    try:
        if v is None or v == "":
            return None
        return float(v)
    except (ValueError, TypeError):
        return None
