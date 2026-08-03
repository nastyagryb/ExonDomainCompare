from __future__ import annotations

import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import model_roles  # noqa: E402
from .boundary_classification import (  # noqa: E402
    DEFAULT_NEAR_EDGE_THRESHOLD_AA,
    canonical_class,
    classify_boundary,
    domain_instance_id,
)

SCHEMA_VERSION = 1
COORDINATE_SYSTEM = "protein_1_based_inclusive"
CORE_SUBDIR = "results/core_gene_analysis"


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_json(path: Path, default: Any = None) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def _int(v: Any, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def _rel(path: Path, root: Path) -> str:
    """Project-root-relative POSIX path; never an absolute personal path."""
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.name


def _fasta_lengths(path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    if not path.is_file():
        return lengths
    pid: Optional[str] = None
    seq_len = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if pid is not None:
                lengths[pid] = seq_len
            pid = line[1:].split()[0].strip()
            seq_len = 0
        else:
            seq_len += len(line.strip())
    if pid is not None:
        lengths[pid] = seq_len
    return lengths




def _resolve_primary_ids(core: Path) -> Dict[str, str]:
    """Map species → primary protein through the shared resolver.

    The model must be built for the same protein the cluster analysed, otherwise every
    domain, TM and boundary lookup — all of which filter on species *and* protein —
    returns nothing and the species ships as ``pending_cluster`` with an empty
    architecture. Resolution therefore goes through one shared module rather than a
    local rule, and a species with no resolvable primary raises instead of receiving an
    alphabetically chosen stand-in.
    """
    from exondomaincompare.framework.primary_resolution import resolve_primaries  # local: optional dep

    return {sid: v["protein_id"] for sid, v in resolve_primaries(core).items()}


# Clade -> display grouping, kept identical to the webapp's species badges so the
# comparative filter and the species list cannot disagree about a species' group.
_CLADE_TO_TAXON_GROUP = {
    "mammal": "Other mammals",
    "primate": "Primates",
    "bird": "Birds",
    "reptile": "Reptiles",
    "amphibian": "Amphibians",
    "fish": "Teleost fish",
}


def _species_metadata(run_dir: Path) -> Dict[str, Dict[str, str]]:
    """Clade and taxonomic group per species, from the run's own species registry.

    Carried into the model provenance so the comparative contract is self-contained:
    the taxonomic-group filter and the exported tables read it from the same place
    rather than re-joining species metadata downstream.
    """
    reg = run_dir / "results" / "01_species_registry" / "species_registry.tsv"
    meta: Dict[str, Dict[str, str]] = {}
    for row in _read_tsv(reg):
        sid = (row.get("species_id") or "").strip()
        if not sid:
            continue
        clade = (row.get("clade") or "").strip().lower()
        meta[sid] = {
            "clade": clade,
            "taxonomic_group": _CLADE_TO_TAXON_GROUP.get(clade, "Analysed species"),
        }
    return meta


def build_models_for_run(run_dir: Path, project_root: Optional[Path] = None) -> Dict[str, Any]:
    """Build one coordinate model per (species, selected primary protein).

    Returns a versioned index dict with a ``models`` list. Missing post-cluster
    inputs (pre-cluster runs) yield models with empty domain/boundary/tm tracks
    and an honest ``status`` — never fabricated features.
    """
    run_dir = Path(run_dir)
    project_root = Path(project_root) if project_root else run_dir.resolve().parents[1]
    core = run_dir / CORE_SUBDIR
    config = _read_json(run_dir / "run_config.json", {}) or {}
    status_doc = _read_json(run_dir / "status.json", {}) or {}

    exon_rows = _read_tsv(core / "exon_protein_map.tsv")
    domain_rows = _read_tsv(core / "domain_features.tsv")
    interpro_rows = _read_tsv(core / "interpro_annotations.tsv")
    tm_rows = _read_tsv(core / "tm_features.tsv")
    boundary_rows = _read_tsv(core / "exon_domain_boundary_distances.tsv")
    candidate_rows = _read_tsv(core / "event_candidate_regions.tsv")
    primary_rows = _read_tsv(core / "primary_selection_evidence.tsv")
    fasta_len = _fasta_lengths(core / "proteins_primary.faa")
    species_meta = _species_metadata(run_dir)

    gene_symbol = str(config.get("gene_symbol") or "").upper()
    sci_names = config.get("species_scientific_names") or {}
    tm_present = (core / "tm_features.tsv").is_file()

    primary_ids = {
        r.get("protein_id")
        for r in primary_rows
        if str(r.get("selected_primary", "")).lower() in ("true", "1", "yes")
    }
    # Species-scoped primaries. The run-level ``primary_ids`` set above is kept only as
    # a last resort for runs whose tables predate species-aware selection.
    primary_by_species = _resolve_primary_ids(core)

    # group exon rows by (species, protein)
    by_sp_protein: Dict[str, Dict[str, List[Dict[str, str]]]] = {}
    for r in exon_rows:
        sp = r.get("species_id") or ""
        pid = r.get("protein_id") or ""
        by_sp_protein.setdefault(sp, {}).setdefault(pid, []).append(r)

    src_files = {
        "exon_protein_map": _rel(core / "exon_protein_map.tsv", project_root),
        "domain_features": _rel(core / "domain_features.tsv", project_root),
        "interpro_annotations": _rel(core / "interpro_annotations.tsv", project_root),
        "tm_features": _rel(core / "tm_features.tsv", project_root),
        "exon_domain_boundary_distances": _rel(core / "exon_domain_boundary_distances.tsv", project_root),
    }

    models: List[Dict[str, Any]] = []
    for sp, proteins in by_sp_protein.items():
        # The primary of THIS species, resolved from the species-aware evidence. The
        # protein must be one this species actually has an exon projection for,
        # otherwise the model would describe a protein it cannot draw.
        pid = primary_by_species.get(sp)
        if pid not in proteins:
            pid = next((p for p in proteins if p in primary_ids), None) or sorted(proteins)[0]
        rows = sorted(proteins[pid], key=lambda x: _int(x.get("exon_number"), 0) or 0)

        proj_len = max((_int(x.get("protein_end_aa"), 0) or 0 for x in rows), default=0)
        protein_length = fasta_len.get(pid, proj_len) or proj_len

        transcript_id = rows[0].get("transcript_id", "") if rows else ""

        exons: List[Dict[str, Any]] = []
        for r in rows:
            s = _int(r.get("protein_start_aa"))
            e = _int(r.get("protein_end_aa"))
            n = _int(r.get("exon_number"))
            # The terminal coding exon projects one residue past the protein when
            # the last codon carries the stop; the stop is not part of the protein
            # sequence, so clamp the aa projection into [1, protein_length].
            if protein_length:
                if s is not None:
                    s = max(1, min(s, protein_length))
                if e is not None:
                    e = max(1, min(e, protein_length))
            exons.append({
                "id": r.get("exon_id") or f"{pid}:exon{n}",
                "label": f"E{n}" if n is not None else (r.get("exon_id") or ""),
                "start": s,
                "end": e,
                "source": r.get("source") or "ncbi_gff",
                "source_file": src_files["exon_protein_map"],
                "status": r.get("confidence") or "gff_cds_derived",
                "tooltip": {
                    "exon_number": n,
                    "exon_id": r.get("exon_id"),
                    "transcript_id": r.get("transcript_id"),
                    "cds_start": _int(r.get("cds_start")),
                    "cds_end": _int(r.get("cds_end")),
                    "phase": _int(r.get("phase")),
                    "protein_aa": [s, e],
                    "shared_exon_group": r.get("shared_exon_group_id") or None,
                },
            })

        rep_domains = _features_from_domain_rows(
            domain_rows, sp, pid, src_files["domain_features"], interpro_rows=interpro_rows)
        families, member_sigs, sites, disorder = _classify_interpro_layers(
            interpro_rows, sp, pid, src_files["interpro_annotations"])
        tms = _tm_features(tm_rows, sp, pid, src_files["tm_features"])
        candidates = _candidate_features(candidate_rows, sp, pid, protein_length)

        _has_domains = bool(rep_domains)
        has_interpro = bool(rep_domains or families or member_sigs or sites or disorder)
        # Post-cluster InterProScan (and its pyTMHMM step) is only complete once
        # integrated InterPro annotation exists. Pre-cluster runs must report the
        # domain / family / site / TM layers as *pending*, never as confirmed absence.
        cluster_complete = has_interpro or str(status_doc.get("post_interpro_status", "")).lower() == "complete"
        status = "available" if cluster_complete else "pending_cluster"

        boundaries = _boundary_features(
            boundary_rows, exons, rep_domains, sp, pid,
            src_files["exon_domain_boundary_distances"],
            cluster_complete=cluster_complete,
            threshold=DEFAULT_NEAR_EDGE_THRESHOLD_AA)

        transcript_models = _transcript_models(run_dir / "website_indices", sp)
        # Enrich primary exon tooltips (genomic coords / strand) from the primary
        # transcript model's rich blocks, keyed by exon id.
        primary_tm = next((t for t in transcript_models if t.get("protein_id") == pid), None)
        if primary_tm:
            by_id = {b.get("id"): b for b in primary_tm.get("blocks", [])}
            for ex in exons:
                rb = by_id.get(ex["id"])
                if rb:
                    ex["tooltip"].update({
                        "genomic_start": rb.get("genomic_start"),
                        "genomic_end": rb.get("genomic_end"),
                        "strand": rb.get("strand"),
                        "shared_exon_group": rb.get("shared_exon_group_id") or ex["tooltip"].get("shared_exon_group"),
                    })

        models.append({
            "schema_version": SCHEMA_VERSION,
            "species_id": sp,
            "scientific_name": sci_names.get(sp) or sp.replace("_", " ").title(),
            "gene_symbol": gene_symbol,
            # A generic run selects one primary protein per species, so that model is
            # this species' primary reference and carries no isoform-level claim of
            # its own. Stating it explicitly changes no behaviour here; it removes
            # the need for a renderer to infer identity from list position.
            "model_id": model_roles.model_id(gene_symbol, sp, "primary"),
            "model_role": model_roles.PRIMARY_REFERENCE,
            "is_primary_reference": True,
            "protein_id": pid,
            "transcript_id": transcript_id,
            "protein_length": protein_length,
            "coordinate_system": COORDINATE_SYSTEM,
            "status": status,
            "pending_info": (None if cluster_complete else {
                "reason": "pending_cluster",
                "cluster_status": status_doc.get("cluster_analysis_status") or status_doc.get("cluster_status"),
                "cluster_command": status_doc.get("cluster_command"),
                "pending_layers": ["representative_domains", "families_superfamilies",
                                   "member_signatures", "functional_sites",
                                   "disorder_regions", "tm_regions"],
                "message": "Domain, family, site and transmembrane layers are pending "
                           "post-cluster InterProScan for this protein.",
            }),
            "tm_analysis": {
                "performed": bool(tm_present and cluster_complete),
                "tm_region_count": len(tms),
                "pending": (not cluster_complete),
                "message": (
                    None if tms
                    else ("No transmembrane region predicted by pyTMHMM" if (tm_present and cluster_complete)
                          else "pyTMHMM pending post-cluster InterProScan")
                ),
            },
            "near_edge_threshold_aa": DEFAULT_NEAR_EDGE_THRESHOLD_AA,
            "exons": exons,
            "exon_boundaries": boundaries,
            "representative_domains": rep_domains,
            "families_superfamilies": families,
            "member_signatures": member_sigs,
            "functional_sites": sites,
            "disorder_regions": disorder,
            "tm_regions": tms,
            "candidate_regions": candidates,
            "transcript_models": transcript_models,
            "n_transcript_models": len(transcript_models),
            "alignment_mapping": {},  # native only; populated for real multi-species runs
            "provenance": {
                "run_id": config.get("run_id") or run_dir.name,
                "gene_symbol": gene_symbol,
                "species_id": sp,
                "coordinate_system": COORDINATE_SYSTEM,
                "generated_by": "src/exondomaincompare/shared_gene_analysis/protein_coordinate_model.py",
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "source_files": src_files,
                "clade": (species_meta.get(sp) or {}).get("clade") or "",
                "taxonomic_group": (species_meta.get(sp) or {}).get("taxonomic_group") or "",
                "primary_protein_source": "primary_selection_evidence.tsv",
                "protein_length_source": "proteins_primary.faa" if pid in fasta_len else "exon_projection",
            },
        })

    index = {
        "schema_version": SCHEMA_VERSION,
        "model_type": "ProteinCoordinateModelIndex",
        "run_id": config.get("run_id") or run_dir.name,
        "gene_symbol": gene_symbol,
        "coordinate_system": COORDINATE_SYSTEM,
        "species_scope": sorted(by_sp_protein.keys()),
        "n_models": len(models),
        "models": models,
    }

    # Cross-species alignment columns for every boundary. Native amino-acid positions
    # are not comparable between species of different lengths, so without this the
    # comparative layer has no evidence to group boundaries by and stays empty.
    # Single-species runs get an explicit "not applicable" record rather than nothing.
    try:
        from .msa_coordinates import (
            annotate_boundaries_with_columns, build_msa_coordinate_map,
        )
        alignment_path = run_dir / "results/generic_gene_analysis/msa/primaries_msa.aln.faa"
        coord_map = build_msa_coordinate_map(alignment_path)
        mapping_report = annotate_boundaries_with_columns(models, coord_map)
        index["msa_coordinate_map"] = {
            k: v for k, v in coord_map.items() if k != "species"
        }
        # Repository-relative: an absolute path would pin the published index to one
        # machine, which is what the no-personal-paths guard exists to prevent.
        if coord_map.get("alignment_file"):
            index["msa_coordinate_map"]["alignment_file"] = _rel(alignment_path, project_root)
        index["msa_coordinate_map"]["species"] = [
            {k: v for k, v in s.items() if k != "native_to_column"}
            for s in coord_map.get("species") or []
        ]
        index["msa_boundary_mapping"] = mapping_report
        for m in models:
            m["alignment_mapping"] = mapping_report.get(m.get("species_id") or "", {})
    except Exception as exc:  # pragma: no cover - mapping is additive
        index["msa_coordinate_map"] = {"available": False, "reason": str(exc)}

    # Derived global "Exon–Domain Boundaries" dashboard contract (page mode,
    # single-species summary + inspection cases + caption, and the multi-species
    # comparative structure). Additive and best-effort — never blocks the model.
    try:
        from .boundary_dashboard import build_boundary_dashboard
        index["boundary_dashboard"] = build_boundary_dashboard(index)
    except Exception:  # pragma: no cover - dashboard is a derived convenience
        pass
    return index


# Curated display forms for InterPro short names whose abbreviations do not read
# as prose. The InterPro accession and raw short name are always kept alongside the
# display label, so nothing is lost — this only controls what a figure or table
# shows a reader.
_DOMAIN_DISPLAY_NAMES = {
    "Ig-like_dom": "Ig-like domain",
    "Ig-like_fold": "Ig-like fold",
    "Ig-like_dom_sf": "Ig-like domain superfamily",
    "Ser-Thr/Tyr_kinase_cat_dom": "Ser-Thr/Tyr kinase domain",
    "Kinase-like_dom_sf": "Kinase-like domain superfamily",
    "Protein_kinase_ATP_BS": "Protein kinase ATP binding site",
    "Tyr_kinase_AS": "Tyrosine kinase active site",
    "FGF_rcpt_fam": "FGF receptor family",
    "RTK": "Receptor tyrosine kinase family",
    "disorder_prediction": "Predicted disorder",
}


def _pretty_domain_name(label: Any) -> str:
    """Readable form of the *real* InterPro short name — never a new name.

    Mirrors ``prettyDomainName`` in ``webapp/frontend/src/pages/viewers/common.js``
    so the numbered instance labels are identical in the model and in the UI.
    """
    raw = str(label or "domain")
    if raw in _DOMAIN_DISPLAY_NAMES:
        return _DOMAIN_DISPLAY_NAMES[raw]
    s = re.sub(r"_", " ", raw)
    return re.sub(r"\bdom\b", "domain", s).strip()


def _supporting_accessions(raw: Any) -> List[str]:
    """Parse the real ``supporting_interpro`` column (``ACC:name;ACC:name``)."""
    out: List[str] = []
    for part in str(raw or "").split(";"):
        acc = part.split(":", 1)[0].strip()
        if acc.startswith("IPR"):
            out.append(acc)
    return out


def _instance_member_signatures(interpro_rows, sp, pid, acc, supporting, start, end):
    """Contributing member-database signatures of ONE domain feature instance.

    A signature belongs to an instance when its integrated InterPro accession is
    part of that instance's real ``supporting_interpro`` set and its own interval
    lies mostly (>50 %) inside the instance span. This is what separates the three
    FGFR1 Ig-like instances: each carries its own ``PS50835`` / ``SM00409`` hit at
    its own coordinates, so the join is instance-resolving, not accession-wide.
    """
    if start is None or end is None:
        return []
    accepted = {a for a in supporting if a} | ({acc} if acc else set())
    out: List[Dict[str, Any]] = []
    for r in interpro_rows or []:
        if r.get("species_id") != sp or r.get("protein_id") != pid:
            continue
        if (r.get("interpro_accession") or "") not in accepted:
            continue
        rs, re_ = _int(r.get("start_aa")), _int(r.get("end_aa"))
        if rs is None or re_ is None:
            continue
        overlap = min(end, re_) - max(start, rs) + 1
        if overlap <= 0 or overlap * 2 < (re_ - rs + 1):
            continue
        out.append({
            "signature_accession": r.get("signature_accession") or None,
            "signature_name": r.get("signature_name") or None,
            "member_database": r.get("member_database") or None,
            "interpro_accession": r.get("interpro_accession") or None,
            "start": rs,
            "end": re_,
            "score_or_evalue": r.get("score_or_evalue") or None,
        })
    return sorted(out, key=lambda x: (x["start"], x["end"], x["signature_accession"] or ""))


def _features_from_domain_rows(rows, sp, pid, source_file,
                               interpro_rows=None) -> List[Dict[str, Any]]:
    """Representative InterPro domains as individually identified feature instances.

    Instance numbering is derived from the sorted start coordinates of the
    instances that share an InterPro accession, so ``IPR007110`` at 33–118 / 145–244
    / 253–355 becomes Ig-like domain 1 / 2 / 3 deterministically.
    """
    raw: List[Dict[str, str]] = []
    for r in rows:
        if r.get("species_id") != sp or r.get("protein_id") != pid:
            continue
        if str(r.get("layer", "domain")).lower() not in ("domain", "representative_domain"):
            continue
        raw.append(r)

    parsed: List[Dict[str, Any]] = []
    for r in raw:
        s, e = _int(r.get("start_aa")), _int(r.get("end_aa"))
        acc = r.get("interpro_accession") or r.get("domain_id") or ""
        parsed.append({"row": r, "start": s, "end": e, "acc": acc})
    parsed.sort(key=lambda p: (p["start"] if p["start"] is not None else 0,
                              p["end"] if p["end"] is not None else 0, p["acc"]))

    # 1-based instance number per accession, ordered by start coordinate.
    per_acc: Dict[str, List[Dict[str, Any]]] = {}
    for p in parsed:
        per_acc.setdefault(p["acc"], []).append(p)
    for items in per_acc.values():
        for n, p in enumerate(items, start=1):
            p["instance_number"] = n
            p["n_instances"] = len(items)

    out: List[Dict[str, Any]] = []
    for order, p in enumerate(parsed, start=1):
        r, s, e, acc = p["row"], p["start"], p["end"], p["acc"]
        label = r.get("interpro_name") or r.get("domain_name") or acc
        base = _pretty_domain_name(label)
        n, many = p["instance_number"], p["n_instances"] > 1
        short_label = f"{base} {n}" if many else base
        supporting = _supporting_accessions(r.get("supporting_interpro"))
        out.append({
            "id": f"{pid}:{acc}:{s}-{e}",
            "domain_instance_id": domain_instance_id(acc, s, e),
            "label": label,
            "short_label": short_label,
            "full_label": f"{short_label} · aa {s}–{e}",
            "instance_number": n,
            "n_instances_of_accession": p["n_instances"],
            "display_order": order,
            "start": s,
            "end": e,
            "feature_type": "representative_domain",
            "interpro_accession": acc or None,
            "supporting_interpro": supporting,
            "member_signatures": _instance_member_signatures(
                interpro_rows, sp, pid, acc, supporting, s, e),
            "source": "InterProScan",
            "source_file": source_file,
            "status": "representative_domain",
            "tooltip": {
                "interpro_accession": acc,
                "domain_instance_id": domain_instance_id(acc, s, e),
                "instance_number": n,
                "interpro_type": r.get("interpro_type") or r.get("feature_type"),
                "member_databases": r.get("member_databases"),
                "n_signatures": _int(r.get("n_signatures")),
                "representative_signature": r.get("representative_signature"),
                "score_or_evalue": r.get("score_or_evalue"),
            },
        })
    return out


_SITE_DBS = {"prosite_patterns", "prosite_profiles", "pirsr", "prints_site"}
_DISORDER_DBS = {"mobidb_lite", "mobidb-lite", "mobidblite"}


def _interpro_feature(r, kind, source_file) -> Dict[str, Any]:
    s, e = _int(r.get("start_aa")), _int(r.get("end_aa"))
    ipr = r.get("interpro_accession") or ""
    sig = r.get("signature_accession") or ""
    acc = ipr or sig
    return {
        "id": f"{r.get('protein_id')}:{kind}:{acc or 'na'}:{s}-{e}",
        "label": (r.get("interpro_name") or r.get("signature_name") or acc or kind),
        "start": s,
        "end": e,
        "feature_type": kind,
        "source": r.get("member_database") or "InterProScan",
        "interpro_accession": ipr or None,
        "signature_accession": sig or None,
        "source_file": source_file,
        "status": (r.get("layer") or kind),
        "tooltip": {
            "interpro_accession": ipr or None,
            "interpro_name": r.get("interpro_name") or None,
            "interpro_type": r.get("interpro_type") or None,
            "signature_accession": sig or None,
            "signature_name": r.get("signature_name") or None,
            "member_database": r.get("member_database") or None,
            "score_or_evalue": r.get("score_or_evalue") or None,
            "is_integrated": str(r.get("is_integrated", "")) in ("1", "true", "True"),
        },
    }


def _classify_interpro_layers(rows, sp, pid, source_file):
    """Bucket real InterProScan annotation rows into scientific feature layers.

    Purely rule-driven over the actual ``interpro_annotations.tsv`` columns
    (``interpro_type`` / ``member_database`` / ``layer``); nothing is fabricated.
    """
    families: List[Dict[str, Any]] = []
    member_sigs: List[Dict[str, Any]] = []
    sites: List[Dict[str, Any]] = []
    disorder: List[Dict[str, Any]] = []
    for r in rows:
        if r.get("species_id") != sp or r.get("protein_id") != pid:
            continue
        itype = str(r.get("interpro_type", "")).strip().lower()
        mdb = str(r.get("member_database", "")).strip().lower()
        layer = str(r.get("layer", "")).strip().lower()
        if mdb in _DISORDER_DBS or "disorder" in str(r.get("signature_name", "")).lower():
            disorder.append(_interpro_feature(r, "disorder_region", source_file))
        elif itype.endswith("site") or mdb in _SITE_DBS:
            sites.append(_interpro_feature(r, "functional_site", source_file))
        elif itype in ("family", "homologous_superfamily") or layer == "family":
            families.append(_interpro_feature(r, "family_superfamily", source_file))
        else:
            # Contributing member-database signature (PFAM/GENE3D/PRINTS/CDD/…),
            # grouped in the UI under its integrated InterPro entry when present.
            member_sigs.append(_interpro_feature(r, "member_signature", source_file))
    key = lambda d: (d["start"] if d["start"] is not None else 0)
    return (sorted(families, key=key), sorted(member_sigs, key=key),
            sorted(sites, key=key), sorted(disorder, key=key))


def _tm_features(rows, sp, pid, source_file) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for i, r in enumerate(rows):
        if r.get("species_id") != sp or r.get("protein_id") != pid:
            continue
        s, e = _int(r.get("start_aa")), _int(r.get("end_aa"))
        if s is None or e is None:
            continue
        out.append({
            "id": f"{pid}:TM{i + 1}",
            "label": f"TM{i + 1}",
            "start": s,
            "end": e,
            "source": r.get("source") or "pyTMHMM",
            "source_file": source_file,
            "status": "tm_region",
            "tooltip": {"topology": r.get("topology")},
        })
    return out


_CDS_RE = re.compile(r"cds(\d+)")


def _exon_by_number(exons: List[Dict[str, Any]], n: Optional[int]) -> Optional[Dict[str, Any]]:
    if n is None:
        return None
    return next((e for e in exons if e.get("tooltip", {}).get("exon_number") == n), None)


def _left_right_exon(boundary_id: str, position: Optional[int],
                     exons: List[Dict[str, Any]]):
    """Resolve the (left, right) coding exon flanking an internal boundary.

    Boundaries are named ``…:cds{k}_end`` (the E{k}→E{k+1} junction) and are
    positioned at the *left* exon's C-terminal aa. We map by exon number first,
    then fall back to position ordering — never recomputing any coordinate.
    """
    m = _CDS_RE.search(boundary_id or "")
    if m:
        k = int(m.group(1))
        left = _exon_by_number(exons, k)
        right = _exon_by_number(exons, k + 1)
        if left is not None and right is not None:
            return left, right
    ordered = [e for e in exons if e.get("end") is not None]
    left = None
    for i, e in enumerate(ordered):
        if position is not None and e["end"] <= position:
            left = e
            right = ordered[i + 1] if i + 1 < len(ordered) else None
    if left is not None:
        return left, right
    return None, None


def _emit_boundary(position, left, right, source_file, *, threshold, mapping_status,
                   boundary_class, signed, absd, nearest, provenance):
    """Assemble one normalized boundary object (Part 1 data contract).

    Carries the canonical field names plus legacy ``*_aa`` / ``category`` aliases
    so the shared SignedDistancePlot and evidence table consume it unchanged.
    """
    left_lbl = left.get("label") if left else None
    right_lbl = right.get("label") if right else None
    label = f"{left_lbl} → {right_lbl}" if (left_lbl and right_lbl) else (
        left_lbl or right_lbl or f"boundary@{position}")
    bid = f"{provenance.get('protein_id')}:{left_lbl}_{right_lbl}" if (left_lbl and right_lbl) \
        else f"{provenance.get('protein_id')}:b{position}"
    nd_id = nearest.get("id")
    return {
        "id": provenance.get("boundary_id") or bid,
        "boundary_id": provenance.get("boundary_id") or bid,
        "nearest_domain_instance_id": nearest.get("domain_instance_id"),
        "nearest_domain_instance_number": nearest.get("instance_number"),
        "label": label,
        "start": position,
        "end": position,
        "protein_position": position,
        "boundary_position_aa": position,           # legacy alias (SignedDistancePlot)
        "exon_boundary_id": provenance.get("boundary_id") or bid,  # legacy alias
        "left_exon_id": left.get("id") if left else None,
        "left_exon_label": left_lbl,
        "right_exon_id": right.get("id") if right else None,
        "right_exon_label": right_lbl,
        "nearest_domain_id": nd_id,
        "nearest_domain_accession": nearest.get("accession"),
        "nearest_domain_label": nearest.get("label"),
        "nearest_domain_short_label": nearest.get("short_label"),
        "nearest_domain_full_label": nearest.get("full_label"),
        "nearest_domain_name": nearest.get("label"),        # legacy alias
        "nearest_domain_start": nearest.get("start"),
        "nearest_domain_end": nearest.get("end"),
        "nearest_edge_type": nearest.get("edge_type"),
        "nearest_edge": nearest.get("edge_type"),           # legacy alias
        "nearest_edge_position": nearest.get("edge_position"),
        "signed_distance": signed,
        "signed_distance_aa": signed,                       # legacy alias
        "absolute_distance": absd,
        "absolute_distance_aa": absd,                       # legacy alias
        "boundary_class": boundary_class,
        "class": boundary_class,
        "category": boundary_class,                         # legacy alias
        "near_threshold": threshold,
        "near_edge_threshold_aa": threshold,
        "mapping_status": mapping_status,
        "domain_layer": "representative_domain",
        "source": provenance.get("source") or "core_post_interpro",
        "source_file": source_file,
        "status": boundary_class if boundary_class else mapping_status,
        "tooltip": {
            "adjacent_exon_transition": label,
            "protein_position": position,
            "nearest_domain": nearest.get("full_label") or nearest.get("label"),
            "nearest_domain_instance_id": nearest.get("domain_instance_id"),
            "nearest_edge": nearest.get("edge_type"),
            "signed_distance": signed,
            "absolute_distance": absd,
            "class": boundary_class,
            "mapping_status": mapping_status,
            "transcript_id": provenance.get("transcript_id"),
            "source": provenance.get("source") or "core_post_interpro",
        },
    }


def _boundary_features(boundary_rows, exons, rep_domains, sp, pid, source_file,
                       *, cluster_complete: bool,
                       threshold: int = DEFAULT_NEAR_EDGE_THRESHOLD_AA) -> List[Dict[str, Any]]:
    """Build the normalized internal coding-exon boundary contract.

    Three honest modes, all driven by the single coordinate model — no React-side
    coordinate maths:

    * **available + precomputed** — read the real classified rows from
      ``exon_domain_boundary_distances.tsv`` (preserving the existing signed
      distances) and enrich each with flanking exons, nearest-domain span, edge
      position, near-edge threshold and mapping status.
    * **available + no table** — recompute from exon junctions using the shared
      classifier (robust fallback; identical vocabulary).
    * **pending cluster** — emit exon-junction *positions only* (no signed
      distance, no classification) so the axis is honest but no domain relation
      is fabricated before InterProScan.

    Sign convention (identical in all three modes and in
    ``boundary_classification.classify_boundary``):
    ``signed_distance = boundary_position - nearest_edge_position``. A boundary
    N-terminal of a domain start is therefore negative, a boundary C-terminal of a
    domain end is positive, and a boundary inside a domain measured against that
    domain's end is negative.
    """
    dom_by_instance = {d["domain_instance_id"]: d for d in rep_domains
                       if d.get("domain_instance_id")}
    dom_by_id = {d.get("id"): d for d in rep_domains if d.get("id")}
    dom_by_acc: Dict[str, List[Dict[str, Any]]] = {}
    for d in rep_domains:
        if d.get("interpro_accession"):
            dom_by_acc.setdefault(d["interpro_accession"], []).append(d)
    ordered = sorted(exons, key=lambda e: (e.get("tooltip", {}).get("exon_number") or 0))
    rows = [r for r in boundary_rows if r.get("species_id") == sp and r.get("protein_id") == pid]

    def _nearest_from_domain(dom, edge_type):
        if not dom:
            return {}
        return {
            "id": dom.get("id"),
            "domain_instance_id": dom.get("domain_instance_id"),
            "instance_number": dom.get("instance_number"),
            "accession": dom.get("interpro_accession"),
            "label": dom.get("label"),
            "short_label": dom.get("short_label"),
            "full_label": dom.get("full_label"),
            "start": dom.get("start"),
            "end": dom.get("end"),
            "edge_type": edge_type,
            "edge_position": (dom.get("start") if edge_type == "start"
                              else dom.get("end") if edge_type == "end" else None),
        }

    def _instance_for_row(r, pos):
        """Resolve the domain feature INSTANCE a stored boundary row was measured against.

        The Core boundary table only persists the InterPro *accession*, which is
        ambiguous whenever an entry is repeated (FGFR1 carries three IPR007110
        Ig-like instances). The instance is therefore recovered from the geometry the
        row does carry: the stored edge type plus signed distance pin the exact edge
        coordinate ``edge = boundary_position - signed_distance``. Only when that
        geometry is absent or matches no real instance do we fall back to a full
        instance-aware recomputation. An accession is never used on its own.
        """
        explicit = r.get("nearest_domain_instance_id")
        if explicit and explicit in dom_by_instance:
            return dom_by_instance[explicit], False
        by_model_id = dom_by_id.get(r.get("nearest_domain_id"))
        if by_model_id is not None:
            return by_model_id, False
        acc = r.get("nearest_domain_accession") or r.get("nearest_domain_id")
        candidates = dom_by_acc.get(acc) or []
        edge_type = r.get("nearest_edge") or r.get("domain_edge_type")
        signed = _int(r.get("signed_distance_aa"))
        if candidates and pos is not None and signed is not None and edge_type in ("start", "end"):
            edge = pos - signed
            exact = [d for d in candidates
                     if (d.get("start") if edge_type == "start" else d.get("end")) == edge]
            if len(exact) == 1:
                return exact[0], False
            if exact:
                # Degenerate: two instances share this edge coordinate — prefer the
                # one that actually contains the boundary, else the first by start.
                inside = [d for d in exact
                          if d.get("start") is not None and d["start"] <= pos <= d["end"]]
                return (inside or exact)[0], False
        # No reconcilable geometry: re-derive the nearest instance from coordinates.
        pool = candidates or rep_domains
        rec = classify_boundary(pos, pool, threshold=threshold)
        return dom_by_instance.get(rec.get("nearest_domain_instance_id")), True

    out: List[Dict[str, Any]] = []

    if rows:
        for r in rows:
            pos = _int(r.get("boundary_position_aa"))
            bid = r.get("exon_boundary_id") or f"{pid}:b{pos}"
            left, right = _left_right_exon(bid, pos, ordered)
            acc = r.get("nearest_domain_accession") or r.get("nearest_domain_id")
            dom, rederived = _instance_for_row(r, pos)
            edge_type = r.get("nearest_edge") or r.get("domain_edge_type")
            signed = _int(r.get("signed_distance_aa"))
            absd = _int(r.get("absolute_distance_aa"))
            cls = canonical_class(r.get("category") or r.get("classification"))
            src = r.get("source")
            if rederived and dom is not None:
                # The instance (and therefore the edge) came from the shared
                # classifier, so the distances and class must come from it too.
                rec = classify_boundary(pos, [dom], threshold=threshold)
                edge_type = rec.get("nearest_edge_type")
                signed = rec.get("signed_distance")
                absd = rec.get("absolute_distance")
                cls = canonical_class(classify_boundary(
                    pos, rep_domains, threshold=threshold).get("class"))
                src = f"{src or 'core_post_interpro'}+instance_reresolved"
            nearest = _nearest_from_domain(dom, edge_type)
            if not nearest and acc:
                nearest = {"id": acc, "accession": acc, "label": r.get("nearest_domain_name"),
                           "domain_instance_id": None, "instance_number": None,
                           "start": None, "end": None, "edge_type": edge_type,
                           "edge_position": None}
            if absd is None and signed is not None:
                absd = abs(signed)
            mapping = "mapped" if nearest.get("domain_instance_id") else (
                "unavailable" if not rep_domains else "unmapped")
            out.append(_emit_boundary(
                pos, left, right, source_file, threshold=threshold, mapping_status=mapping,
                boundary_class=cls, signed=signed, absd=absd, nearest=nearest,
                provenance={"protein_id": pid, "boundary_id": bid,
                            "transcript_id": r.get("transcript_id"), "source": src}))
        return out

    if not cluster_complete:
        # pending: exon-junction positions only (no domain relation before cluster)
        for i in range(len(ordered) - 1):
            left, right = ordered[i], ordered[i + 1]
            pos = left.get("end")
            if pos is None:
                continue
            k = left.get("tooltip", {}).get("exon_number")
            bid = f"{pid}:cds{k}_end" if k is not None else f"{pid}:b{pos}"
            out.append(_emit_boundary(
                pos, left, right, source_file, threshold=threshold,
                mapping_status="pending_cluster", boundary_class=None,
                signed=None, absd=None, nearest={},
                provenance={"protein_id": pid, "boundary_id": bid, "source": "coordinate_model_pending"}))
        return out

    # available but no precomputed table: recompute from exon junctions
    for i in range(len(ordered) - 1):
        left, right = ordered[i], ordered[i + 1]
        pos = left.get("end")
        if pos is None:
            continue
        rec = classify_boundary(pos, rep_domains, threshold=threshold)
        dom = (dom_by_instance.get(rec.get("nearest_domain_instance_id"))
               or dom_by_id.get(rec.get("nearest_domain_id")))
        nearest = _nearest_from_domain(dom, rec.get("nearest_edge_type"))
        if not nearest and rec.get("nearest_domain_id"):
            nearest = {"id": rec.get("nearest_domain_id"),
                       "domain_instance_id": rec.get("nearest_domain_instance_id"),
                       "instance_number": rec.get("nearest_domain_instance_number"),
                       "accession": rec.get("nearest_domain_accession"),
                       "label": rec.get("nearest_domain_label"), "start": None, "end": None,
                       "edge_type": rec.get("nearest_edge_type"),
                       "edge_position": rec.get("nearest_edge_position")}
        k = left.get("tooltip", {}).get("exon_number")
        bid = f"{pid}:cds{k}_end" if k is not None else f"{pid}:b{pos}"
        mapping = "mapped" if rec.get("nearest_domain_id") else ("unavailable" if not rep_domains else "unmapped")
        out.append(_emit_boundary(
            pos, left, right, source_file, threshold=threshold, mapping_status=mapping,
            boundary_class=rec.get("class"), signed=rec.get("signed_distance"),
            absd=rec.get("absolute_distance"), nearest=nearest,
            provenance={"protein_id": pid, "boundary_id": bid, "source": "coordinate_model_reclassified"}))
    return out


# How many packed lanes the architecture track shows before "Show all candidates".
# Packing runs in rank order, so lane 0 already holds the highest-ranked maximal set
# of non-overlapping clusters and each further lane adds the next-best ones.
DEFAULT_CANDIDATE_LANES = 3
_CONFIDENCE_ORDER = {"high": 0, "medium": 1, "low": 2}


def _raw_candidate_id(sp: str, row: Dict[str, str], s: int, e: int) -> str:
    """Stable identifier for ONE raw ``event_candidate_regions.tsv`` row.

    The table has no id column, so the id is built from the fields that identify the
    row: the compared protein pair, the interval and the difference kind. It exists so
    a display cluster can name every raw candidate it stands for without the raw table
    being rewritten.
    """
    a = (row.get("protein_a") or "").strip()
    b = (row.get("protein_b") or "").strip()
    kind = (row.get("candidate_type") or "").strip() or "candidate"
    return f"{sp}:{a}|{b}:{s}-{e}:{kind}"


def _candidate_features(rows, sp, pid, protein_length=None) -> List[Dict[str, Any]]:
    """Exploratory candidate regions as deterministic *display clusters*.

    Two separate concerns, deliberately kept apart:

    * **Coordinate validity.** Each candidate is a difference between *two* isoforms,
      measured in that pair's own coordinate frame. Only a pair that includes the
      primary protein has coordinates on the primary protein's axis; a difference
      between two alternative isoforms is measured on an isoform that may be much
      longer. Placing those on the primary axis anyway produced features past the end
      of the protein, which failed coordinate validation and silently dropped the whole
      model — for every species, not just the one with the long isoform.

    * **Legibility.** PKM produces enough overlapping candidates that a single track
      row drew boxes and labels on top of each other. The repair is a display
      representation, never a change to the evidence: raw rows stay complete in the
      TSV, and rows are grouped only when they describe the *same primary-protein
      interval* in the *same alignment block* (difference kind + evidence source).
      An insertion and a deletion reported over the same interval are distinct blocks
      and stay distinct clusters. Every cluster then receives a deterministic rank
      (from its support), a greedy lane assignment in rank order so no two boxes ever
      share a lane, and the member ids and counts needed to state what it stands for.
    """
    grouped: Dict[tuple, Dict[str, Any]] = {}
    off_axis = 0
    for r in rows:
        sp_row = r.get("species_id")
        if sp_row and sp_row != sp:
            continue
        s = _int(r.get("candidate_start_aa") or r.get("aa_start") or r.get("start_aa"))
        e = _int(r.get("candidate_end_aa") or r.get("aa_end") or r.get("end_aa"))
        if s is None or e is None:
            continue
        # The scan states which protein a row's coordinates belong to: an insertion
        # exists only in the second protein and is measured there, a deletion or
        # substitution block is measured on the first. Where that column is present it
        # is the answer; membership of the pair is not, because the other protein's
        # coordinates can fall inside the primary's length and silently land on it.
        reference = (r.get("coordinate_reference_protein") or "").strip()
        if reference:
            if reference != pid:
                off_axis += 1
                continue
        else:
            pair = {r.get("protein_a"), r.get("protein_b")} - {None, ""}
            if pair and pid not in pair:
                off_axis += 1
                continue
        if protein_length and (s < 1 or e > protein_length):
            off_axis += 1
            continue
        ctype = (r.get("candidate_type") or "").strip()
        evidence = (r.get("evidence") or "").strip()
        key = (s, e, ctype, evidence)
        g = grouped.setdefault(key, {
            "start": s, "end": e, "candidate_type": ctype or None, "evidence": evidence or None,
            "confidences": set(), "affected": set(), "members": [], "comparisons": 0,
        })
        g["comparisons"] += 1
        g["members"].append(_raw_candidate_id(sp, r, s, e))
        if r.get("confidence"):
            g["confidences"].add(r.get("confidence"))
        for k in ("protein_a", "protein_b"):
            if r.get(k):
                g["affected"].add(r.get(k))

    def best_conf(cs):
        return sorted(cs, key=lambda c: _CONFIDENCE_ORDER.get(c, 9))[0] if cs else None

    # Positional order gives the stable C1..Cn identity; the rank below is a separate,
    # support-driven order and does not renumber anything.
    positional = sorted(grouped.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3]))
    clusters: List[Dict[str, Any]] = []
    for i, (key, g) in enumerate(positional):
        s, e = g["start"], g["end"]
        cid = f"C{i + 1}"
        supporters = sorted(g["affected"] - {pid})
        clusters.append({
            "id": cid,
            "display_cluster_id": cid,
            "label": f"{cid} · aa {s}–{e}",
            "start": s,
            "end": e,
            "length_aa": e - s + 1,
            "source": "isoform_difference_scan",
            "source_file": "results/core_gene_analysis/event_candidate_regions.tsv",
            "status": "exploratory",
            "candidate_type": g["candidate_type"],
            "evidence": g["evidence"],
            "alignment_block": f"{g['candidate_type'] or 'candidate'}:{g['evidence'] or 'na'}",
            "confidence": best_conf(g["confidences"]),
            "member_candidate_ids": sorted(g["members"]),
            "n_member_candidates": len(g["members"]),
            "n_comparisons": g["comparisons"],
            "n_supporting_isoforms": len(supporters),
            "supporting_isoforms": supporters,
            "affected_proteins": sorted(g["affected"]),
        })

    # Rank: strongest evidence first, ties broken deterministically so two runs over
    # the same table always produce the same lanes.
    ranked = sorted(clusters, key=lambda c: (
        _CONFIDENCE_ORDER.get(c["confidence"], 9),
        -c["n_comparisons"],
        -c["n_supporting_isoforms"],
        -c["length_aa"],
        c["start"],
        c["id"],
    ))
    lane_ends: List[int] = []
    for rank, c in enumerate(ranked, start=1):
        lane = next((i for i, end in enumerate(lane_ends) if c["start"] > end), len(lane_ends))
        if lane == len(lane_ends):
            lane_ends.append(c["end"])
        else:
            lane_ends[lane] = c["end"]
        c["display_rank"] = rank
        c["display_lane"] = lane
        c["is_top_ranked"] = lane < DEFAULT_CANDIDATE_LANES

    for c in clusters:
        c["tooltip"] = {
            "display_cluster_id": c["id"],
            "candidate_type": c["candidate_type"],
            "evidence": c["evidence"],
            "alignment_block": c["alignment_block"],
            "confidence": c["confidence"],
            "n_member_candidates": c["n_member_candidates"],
            "n_comparisons": c["n_comparisons"],
            "n_supporting_isoforms": c["n_supporting_isoforms"],
            "member_candidate_ids": c["member_candidate_ids"],
            "affected_proteins": c["affected_proteins"],
            "display_rank": c["display_rank"],
            "display_lane": c["display_lane"],
        }
        c["display_lane_count"] = len(lane_ends)
        c["default_display_lanes"] = DEFAULT_CANDIDATE_LANES
        if off_axis:
            # Reported, not hidden: the count says how much exploratory evidence exists
            # for this species that simply has no primary-protein coordinate.
            c["off_primary_axis_candidates"] = off_axis
    return clusters


def _transcript_models(index_dir: Path, sp: str) -> List[Dict[str, Any]]:
    """Reuse the already-built rich coordinate_track_index for per-transcript models.

    This keeps the coordinate model the single served source of truth for the Exon
    Map without re-deriving genomic/strand/shared-exon-group fields.
    """
    for name in ("coordinate_track_index.json", "generic/coordinate_track_index.json"):
        idx = _read_json(index_dir / name, None)
        if idx:
            break
    else:
        return []
    species = idx.get("species") or []
    sp_entry = next((s for s in species if s.get("species") == sp), None) or (species[0] if species else None)
    if not sp_entry:
        return []
    models = sp_entry.get("models") or []
    out: List[Dict[str, Any]] = []
    for m in models:
        blocks = [
            {
                "id": b.get("id"),
                "label": b.get("label"),
                "start": b.get("start"),
                "end": b.get("end"),
                "exon_number": b.get("exon_number"),
                "transcript_exon_number": b.get("transcript_exon_number"),
                "transcript_id": b.get("transcript_id"),
                "shared_exon_group_id": b.get("shared_exon_group_id"),
                "genomic_start": b.get("genomic_start"),
                "genomic_end": b.get("genomic_end"),
                "cds_start": b.get("cds_start"),
                "cds_end": b.get("cds_end"),
                "phase": b.get("phase"),
                "strand": b.get("strand"),
                "source": b.get("source"),
            }
            for b in (m.get("blocks") or [])
            if b.get("feature_type") in (None, "coding_exon")
        ]
        out.append({
            "protein_id": m.get("protein_id"),
            "transcript_id": m.get("transcript_id"),
            "protein_length": m.get("protein_length"),
            "curation_status": m.get("curation_status") or "predicted",
            "is_primary": bool(m.get("is_primary")),
            "role": m.get("role"),
            "exon_count": len(blocks),
            "blocks": blocks,
        })
    return out


if __name__ == "__main__":  # pragma: no cover - manual inspection helper
    import argparse

    ap = argparse.ArgumentParser(description="Build the protein coordinate model for a run.")
    ap.add_argument("run_dir")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    idx = build_models_for_run(Path(args.run_dir))
    text = json.dumps(idx, indent=2)
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"wrote {args.out} ({idx['n_models']} models)")
    else:
        print(text)
