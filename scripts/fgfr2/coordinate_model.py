"""The validated FGFR2 dataset, expressed in the shared coordinate-model contract.

The FGFR2 analysis was written before the shared, gene-agnostic figure system
existed, so it grew its own plotting code and its own card catalogue. That is why
its Figure Gallery looks nothing like the one a FGFR1 or TP53 run gets, and why a
species inside the 30-species panel has no gallery of its own.

Nothing here re-analyses FGFR2. It reads the frozen result tables and restates
them in the structure ``src/exondomaincompare/shared_gene_analysis/protein_coordinate_model.py``
produces for every other gene, so the accepted shared renderer — the same one the
Gene Explorer exports through — can draw the validated FGFR2 features. One
renderer, one set of semantic styles, one export path.

What is read, and what is deliberately not derived:

* Coding exons, the IIIb/IIIc cassette slot, the representative Ig and kinase
  domains and the pyTMHMM helix come from
  ``15_exon_domain_boundary_post_interpro/tables/exon_domain_architecture_features.tsv``.
  Those intervals are the validated selection and are copied, never recomputed.
* Member-database signatures come from ``interpro_domain_features_normalized.tsv``
  and are attached to the representative instance they support.
* Internal exon boundaries are *not* in the freeze, because the FGFR2 analysis
  classified cassette boundaries only. They are derived here with the shared
  generic classifier (``boundary_classification.py``), which is a different
  vocabulary from the frozen FGFR2 Boundary Consistency classes and is kept
  separate from them: the frozen cassette values are never touched, recomputed or
  overwritten.

One model per (species, isoform), because that is the unit the FGFR2 analysis
works in: IIIb and IIIc are two real proteins, not two views of one. Where a
comparative figure needs one row per species, :func:`comparative_primary_ids`
states the rule it uses instead of leaving the choice to iteration order.
"""
from __future__ import annotations

import csv
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from exondomaincompare.shared_gene_analysis import boundary_classification as bc  # noqa: E402
from exondomaincompare.shared_gene_analysis import model_roles as mr  # noqa: E402
from exondomaincompare.shared_gene_analysis import species_order as so  # noqa: E402

SCHEMA_VERSION = 1
COORDINATE_SYSTEM = "protein_1_based_inclusive"
GENE_SYMBOL = "FGFR2"

FREEZE = ROOT / "results" / "final_30_until_interpro_prepare"
CLOSURE = FREEZE / "13_final_pre_interpro_closure"
ARCHITECTURE = FREEZE / "15_exon_domain_boundary_post_interpro"

FEATURES_TSV = ARCHITECTURE / "tables" / "exon_domain_architecture_features.tsv"
INTERPRO_TSV = ARCHITECTURE / "tables" / "interpro_domain_features_normalized.tsv"
TRUTH_TSV = CLOSURE / "final_pre_interpro_truth_table.tsv"
FULL_MSA = CLOSURE / "MSA" / "final_fgfr2_full_length_protein_msa.aln.faa"

#: Which model of a species is its primary reference — the one a comparative figure
#: puts in that species' row. FGFR2's reference-database canonical protein is the
#: IIIc form (RefSeq NM_000141 / Ensembl ENST00000358487), so IIIc is the stated
#: default rather than whichever isoform a dictionary happened to yield first. Where
#: a species has no IIIc model the IIIb model is the reference and says so.
PRIMARY_REFERENCE_PREFERENCE = ("IIIc", "IIIb")

#: The isoforms the FGFR2 analysis validates. Used to enumerate what a species
#: *should* have, so a combination that is absent is reported as absent with its
#: reason rather than silently missing from the inventory.
VALIDATED_ISOFORMS = ("IIIb", "IIIc")

#: Curated display names for the FGFR2 representative domain labels. The stored
#: label stays alongside, so nothing is renamed — this only controls what a reader
#: sees on a figure.
_DOMAIN_DISPLAY = {
    "Ig1": "Ig-like domain 1",
    "Ig2": "Ig-like domain 2",
    "Ig3": "Ig-like domain 3",
    "kinase": "Ser-Thr/Tyr kinase domain",
}

#: Member-database preference, identical to the one the freeze used to pick a
#: representative interval, so the accession this model reports for an instance is
#: the one the validated selection was based on.
_DB_PRIORITY = {
    "CDD": 1, "Pfam": 2, "ProSiteProfiles": 3, "SMART": 4, "PRINTS": 5,
    "Gene3D": 6, "SUPERFAMILY": 7, "FunFam": 8, "PANTHER": 9, "PIRSF": 10,
    "NCBIfam": 11, "SFLD": 12, "Hamap": 12, "PIRSR": 12, "ProSitePatterns": 13,
}


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _int(value: Any) -> Optional[int]:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def _rel(path: Path) -> str:
    resolved = Path(path).resolve()
    for parent in resolved.parents:
        if (parent / "run_config.json").is_file():
            try:
                return resolved.relative_to(parent).as_posix()
            except ValueError:
                break
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return Path(path).name


# --------------------------------------------------------------------------- #
# domains
# --------------------------------------------------------------------------- #
def _signature_support(interpro_rows: Sequence[Dict[str, str]], pid: str,
                       start: int, end: int) -> Tuple[str, List[str], List[Dict[str, Any]]]:
    """The signatures that support one validated representative interval.

    An interval the freeze produced by taking the union of a cluster — the kinase
    region is one — matches no single signature exactly, so an accession cannot be
    read off by equality. Every signature of the protein that overlaps the interval
    is reported instead, and the accession named as *the* accession is the one from
    the highest-priority member database, which is the rule the freeze itself used
    to choose the interval. Where no overlapping signature carries an InterPro
    accession the field stays empty rather than being filled with a guess.
    """
    overlapping = []
    for r in interpro_rows:
        if r.get("protein_id") != pid:
            continue
        s, e = _int(r.get("domain_start_aa")), _int(r.get("domain_end_aa"))
        if s is None or e is None or e < start or s > end:
            continue
        overlapping.append(r)

    signatures = [{
        "member_database": r.get("member_database") or "",
        "signature_accession": r.get("signature_accession") or "",
        "signature_description": r.get("signature_description") or "",
        "interpro_accession": r.get("interpro_accession") or "",
        "interpro_description": r.get("interpro_description") or "",
        "start": _int(r.get("domain_start_aa")),
        "end": _int(r.get("domain_end_aa")),
    } for r in sorted(overlapping, key=lambda r: (
        _DB_PRIORITY.get(r.get("member_database") or "", 99),
        _int(r.get("domain_start_aa")) or 0))]

    accessions = list(dict.fromkeys(s["interpro_accession"] for s in signatures
                                    if s["interpro_accession"]))
    primary = next((s["interpro_accession"] for s in signatures
                    if s["interpro_accession"]), "")
    return primary, accessions, signatures


def _representative_domains(rows: Sequence[Dict[str, str]],
                            interpro_rows: Sequence[Dict[str, str]],
                            pid: str) -> List[Dict[str, Any]]:
    """The validated representative domains, as identified feature instances.

    The identity of an instance is its label and its interval, not its InterPro
    accession: three Ig-like domains of one protein share an accession, so an
    accession alone can never name one of them.

    Two numberings coexist here on purpose, and conflating them would misstate the
    freeze. ``positional_label`` is the validated FGFR2 Ig1/Ig2/Ig3 numbering,
    which the freeze assigns N→C and explicitly documents as positional only.
    ``instance_number`` is the shared contract's per-accession index. They differ
    because Ig1 and Ig3 share ``IPR007110`` while Ig2 is ``IPR013098``, so the
    second Ig-like domain by position is the first instance of its accession.
    """
    parsed = []
    for r in rows:
        if r.get("status") != "representative":
            continue
        start, end = _int(r.get("start_aa")), _int(r.get("end_aa"))
        if start is None or end is None:
            continue
        parsed.append((start, end, r))
    parsed.sort(key=lambda x: (x[0], x[1]))

    accessions = [_signature_support(interpro_rows, pid, s, e)[0] for s, e, _ in parsed]
    per_accession: Dict[str, int] = {}
    totals: Dict[str, int] = {}
    for acc in accessions:
        totals[acc] = totals.get(acc, 0) + 1

    out: List[Dict[str, Any]] = []
    for order, ((start, end, r), acc) in enumerate(zip(parsed, accessions), start=1):
        raw_label = r.get("feature_label") or ""
        display = _DOMAIN_DISPLAY.get(raw_label, raw_label)
        per_accession[acc] = per_accession.get(acc, 0) + 1
        n, total = per_accession[acc], totals[acc]
        _, supporting, signatures = _signature_support(interpro_rows, pid, start, end)
        instance_id = bc.domain_instance_id(acc or raw_label, start, end)
        out.append({
            "id": f"{pid}:{raw_label}:{start}-{end}",
            "domain_instance_id": instance_id,
            "label": raw_label,
            "positional_label": display,
            "short_label": display,
            "full_label": f"{display} · aa {start}–{end}",
            "instance_number": n,
            "n_instances_of_accession": total,
            "display_order": order,
            "start": start,
            "end": end,
            "feature_type": "representative_domain",
            "interpro_accession": acc or None,
            "supporting_interpro": supporting,
            "member_signatures": signatures,
            "source": r.get("source") or "interproscan",
            "source_file": _rel(FEATURES_TSV),
            "status": "representative_domain",
            "tooltip": {
                "interpro_accession": acc,
                "domain_instance_id": instance_id,
                "instance_number": n,
                "validated_label": raw_label,
                "n_signatures": len(signatures),
                "interval_selection": (
                    "Validated representative interval from the FGFR2 freeze. The "
                    "accession is that of the highest-priority overlapping member "
                    "signature; the kinase interval is a cluster union, so no single "
                    "signature matches it exactly."),
                "numbering": (
                    f"{display} is the validated FGFR2 positional numbering (N→C), "
                    "which the freeze documents as positional only. Instance "
                    f"{n} of {total} refers to this InterPro accession."),
            },
        })
    return out


# --------------------------------------------------------------------------- #
# exons, cassette and transmembrane helix
# --------------------------------------------------------------------------- #
def _exon_number(label: str) -> Optional[int]:
    """`exon 9 (IIIb cassette)` -> 9."""
    parts = str(label or "").split()
    return _int(parts[1]) if len(parts) > 1 else None


def _exons_and_cassette(rows: Sequence[Dict[str, str]], pid: str,
                        ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """The coding exon series, with the cassette slot in its place within it.

    The cassette is exon 9 of the series, not a feature beside it. Keeping it in
    the exon track is what lets a reader see IIIb or IIIc in the context of all
    coding exons, which is the whole point of the FGFR2 architecture figure.
    """
    collected: List[Tuple[int, Dict[str, Any]]] = []
    cassettes: List[Dict[str, Any]] = []
    for r in rows:
        ftype = r.get("feature_type") or ""
        if ftype not in ("coding_exon", "IIIb_slot", "IIIc_slot"):
            continue
        start, end = _int(r.get("start_aa")), _int(r.get("end_aa"))
        if start is None or end is None:
            continue
        label = r.get("feature_label") or ""
        number = _exon_number(label)
        is_cassette = ftype.endswith("_slot")
        cassette_id = ftype[:-5] if is_cassette else ""
        exon = {
            "id": f"{pid}:exon{number}" if number else f"{pid}:{start}-{end}",
            "label": f"E{number}" if number else label,
            "start": start,
            "end": end,
            "source": r.get("source") or "figure3C",
            "source_file": _rel(FEATURES_TSV),
            "status": r.get("status") or "coordinate_mapped",
            "is_cassette_exon": is_cassette,
            "cassette_id": cassette_id,
            "tooltip": {
                "exon_number": number,
                "validated_label": label,
                "transcript_id": r.get("transcript_id"),
                "protein_aa": [start, end],
                "cassette": cassette_id or None,
            },
        }
        collected.append((number if number is not None else start, exon))
        if is_cassette:
            cassettes.append({
                "id": f"{pid}:{cassette_id}_cassette",
                "cassette_id": cassette_id,
                "label": f"{cassette_id} cassette",
                "exon_label": f"E{number}" if number else label,
                "start": start,
                "end": end,
                "feature_type": "validated_event",
                "source": r.get("source") or "figure3C",
                "source_file": _rel(FEATURES_TSV),
                "status": r.get("status") or "coordinate_mapped",
                "tooltip": {
                    "cassette": cassette_id,
                    "exon_number": number,
                    "claim": ("Validated IIIb/IIIc cassette slot of the FGFR2 "
                              "freeze. Coordinates are the validated projection."),
                },
            })
    collected.sort(key=lambda x: x[0])
    return [e for _, e in collected], cassettes


def _tm_regions(rows: Sequence[Dict[str, str]], pid: str) -> List[Dict[str, Any]]:
    out = []
    for i, r in enumerate(rows, start=1):
        if r.get("feature_type") != "transmembrane_pytmhmm":
            continue
        start, end = _int(r.get("start_aa")), _int(r.get("end_aa"))
        if start is None or end is None:
            continue
        out.append({
            "id": f"{pid}:TM{i}",
            "label": r.get("feature_label") or "TM",
            "start": start,
            "end": end,
            "source": "pytmhmm",
            "source_file": _rel(FEATURES_TSV),
            "status": r.get("status") or "tm_region",
            "tooltip": {"topology": "transmembrane",
                        "validated_status": r.get("status") or ""},
        })
    return out


def _families(rows: Sequence[Dict[str, str]], pid: str,
              interpro_rows: Sequence[Dict[str, str]]) -> List[Dict[str, Any]]:
    """Family-level fingerprints, kept in their own neutral layer.

    A family fingerprint spans most of the protein and is not a structural domain,
    so it must not sit in the representative track where a reader would read it as
    one.
    """
    out = []
    for r in rows:
        if r.get("status") != "family_level":
            continue
        start, end = _int(r.get("start_aa")), _int(r.get("end_aa"))
        if start is None or end is None:
            continue
        acc, _, _ = _signature_support(interpro_rows, pid, start, end)
        out.append({
            "id": f"{pid}:family:{start}-{end}",
            "label": r.get("feature_label") or "family fingerprint",
            "start": start,
            "end": end,
            "feature_type": "family_superfamily",
            "interpro_accession": acc or None,
            "source": r.get("source") or "interproscan",
            "source_file": _rel(FEATURES_TSV),
            "status": "family_level",
            "tooltip": {"layer": "family / superfamily, not a structural domain"},
        })
    return out


# --------------------------------------------------------------------------- #
# internal exon boundaries
# --------------------------------------------------------------------------- #
def _boundaries(exons: Sequence[Dict[str, Any]], domains: Sequence[Dict[str, Any]],
                pid: str, transcript_id: str,
                threshold: int = bc.DEFAULT_NEAR_EDGE_THRESHOLD_AA,
                ) -> List[Dict[str, Any]]:
    """Internal coding-exon boundaries, classified with the shared generic rule.

    The FGFR2 freeze classified the two cassette boundaries only, in its own
    vocabulary and against its own thresholds. Those values are validated and are
    not read, recomputed or contradicted here. This is the generic exon-boundary
    analysis every other gene gets, applied to the FGFR2 proteins, and it is
    labelled as such in every figure and table it reaches.
    """
    out: List[Dict[str, Any]] = []
    for left, right in zip(exons, exons[1:]):
        position = left.get("end")
        if position is None:
            continue
        rec = bc.classify_boundary(position, domains, threshold=threshold)
        label = f"{left.get('label')} → {right.get('label')}"
        bid = f"{pid}:{left.get('label')}_{right.get('label')}"
        cassette = left.get("cassette_id") or right.get("cassette_id") or ""
        out.append({
            "id": bid,
            "boundary_id": bid,
            "exon_boundary_id": bid,
            "label": label,
            "start": position,
            "end": position,
            "protein_position": position,
            "boundary_position_aa": position,
            "left_exon_id": left.get("id"),
            "left_exon_label": left.get("label"),
            "right_exon_id": right.get("id"),
            "right_exon_label": right.get("label"),
            "is_cassette_boundary": bool(cassette),
            "cassette_id": cassette,
            "nearest_domain_id": rec["nearest_domain_id"],
            "nearest_domain_instance_id": rec["nearest_domain_instance_id"],
            "nearest_domain_instance_number": rec["nearest_domain_instance_number"],
            "nearest_domain_accession": rec["nearest_domain_accession"],
            "nearest_domain_label": rec["nearest_domain_label"],
            "nearest_domain_name": rec["nearest_domain_label"],
            # The names the shared comparative boundary contract reads. Emitting them
            # here is what lets the generic cross-species boundary analysis run on
            # FGFR2 unchanged, instead of an FGFR2-specific copy of it.
            "nearest_domain_short_label": rec["nearest_domain_label"],
            "nearest_domain_full_label": rec.get("nearest_domain_full_label")
                                         or rec["nearest_domain_label"],
            "nearest_domain_start": rec["nearest_domain_start"],
            "nearest_domain_end": rec["nearest_domain_end"],
            "nearest_edge_type": rec["nearest_edge_type"],
            "nearest_edge": rec["nearest_edge_type"],
            "nearest_edge_position": rec["nearest_edge_position"],
            "signed_distance": rec["signed_distance"],
            "signed_distance_aa": rec["signed_distance"],
            "absolute_distance": rec["absolute_distance"],
            "absolute_distance_aa": rec["absolute_distance"],
            "boundary_class": rec["class"],
            "class": rec["class"],
            "category": rec["class"],
            "near_threshold": threshold,
            "near_edge_threshold_aa": threshold,
            "mapping_status": "derived_from_validated_exon_projection",
            "domain_layer": "representative_domain",
            "source": "shared_generic_boundary_classification",
            "source_file": _rel(FEATURES_TSV),
            "status": rec["class"],
            "tooltip": {
                "adjacent_exon_transition": label,
                "protein_position": position,
                "nearest_domain": rec["nearest_domain_label"],
                "nearest_edge": rec["nearest_edge_type"],
                "signed_distance": rec["signed_distance"],
                "class": rec["class"],
                "transcript_id": transcript_id,
                "vocabulary": (
                    "Generic exon-boundary classes (exact / near / inside / outside). "
                    "This is not the frozen FGFR2 Boundary Consistency vocabulary, "
                    "whose validated cassette values are reported unchanged in the "
                    "Boundary Consistency views."),
            },
        })
    return out


# --------------------------------------------------------------------------- #
# the index
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Sources:
    """The four tables this index is derived from.

    Held in one object rather than read from module constants, so the same builder can
    describe the validated freeze and an individual FGFR2 run. The pre-InterPro pipeline
    writes the same table names for both, which is what makes one builder sufficient.
    """

    features: Path
    interpro: Path
    truth: Path
    msa: Path
    run_id: str = "example"
    dataset_id: str = "example"
    dataset_label: str = "validated FGFR2 30-species freeze"

    @property
    def post_cluster_available(self) -> bool:
        """Whether the domain layer exists.

        Without it a model has no domains, no exon series and no boundaries — but it still
        has an identity, a protein and a validated isoform label. That is the difference
        between a dataset with no models and a dataset whose models are awaiting a layer.
        """
        return self.features.is_file()


def freeze_sources() -> Sources:
    """The validated 30-species dataset — this module's original and default subject."""
    return Sources(features=FEATURES_TSV, interpro=INTERPRO_TSV,
                   truth=TRUTH_TSV, msa=FULL_MSA)


def run_sources(run_dir: Path) -> Sources:
    """One FGFR2 run's own tables, in the layout its pipeline writes."""
    run_dir = Path(run_dir)
    results = run_dir / "results"
    closure = results / "13_final_pre_interpro_closure"
    architecture = results / "15_exon_domain_boundary_post_interpro"
    return Sources(
        features=architecture / "tables" / "exon_domain_architecture_features.tsv",
        interpro=architecture / "tables" / "interpro_domain_features_normalized.tsv",
        truth=closure / "final_pre_interpro_truth_table.tsv",
        msa=closure / "MSA" / "final_fgfr2_full_length_protein_msa.aln.faa",
        run_id=run_dir.name,
        dataset_id=f"run:{run_dir.name}",
        dataset_label=f"FGFR2 run {run_dir.name}",
    )


def _truth_by_protein(truth: Path) -> Dict[str, Dict[str, str]]:
    return {r.get("protein_id", ""): r for r in _read_tsv(truth) if r.get("protein_id")}


def _truth_by_combination(truth: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    """``{(species, isoform): truth row}`` for every validated combination."""
    return {(r.get("species", ""), r.get("isoform", "")): r
            for r in _read_tsv(truth) if r.get("species")}


def _review_status(row: Dict[str, str]) -> str:
    """The review state of one combination, from the freeze's own status columns.

    Reported as the freeze recorded it. A flag that the freeze raised is not
    downgraded here, and none is invented.
    """
    flags = []
    for column in ("coordinate_validation_status", "protein_integrity_status",
                   "MSA_full_length_status", "MSA_cassette_status",
                   "boundary_robustness_class"):
        value = (row.get(column) or "").strip()
        if "review" in value.lower() or "outlier" in value.lower():
            flags.append(f"{column}={value}")
    readiness = (row.get("pre_interpro_readiness_class") or "").strip()
    if "supplement_review_only" in readiness:
        return "supplement_review_only"
    if flags:
        return "flagged_for_review: " + "; ".join(flags)
    return "no_review_flag"


def _exon_block_status(rows: Sequence[Dict[str, str]]) -> str:
    """How the architecture's exon blocks were obtained, from the freeze's own status.

    The freeze distinguishes three outcomes, and they decide what can honestly be
    drawn for a protein: ``coordinate_mapped`` (the projection is direct),
    ``native_exon_blocks_reconstructed`` (the blocks were reconstructed) and
    ``cassette_only_high_confidence`` (only the cassette coordinate was validated
    and there is no coding-exon series at all).
    """
    statuses = {(r.get("status") or "").strip() for r in rows
                if r.get("track") == "cassette"}
    for candidate in ("cassette_only_high_confidence",
                      "native_exon_blocks_reconstructed", "coordinate_mapped"):
        if candidate in statuses:
            return candidate
    return "not_recorded"


def _msa_columns(msa: Path) -> Dict[str, Dict[int, int]]:
    """Native position → alignment column, per protein.

    The shared MSA helper keys on species, because a generic run has one primary
    per species. FGFR2 has two, so the mapping is built per protein here; the
    alignment file and the column arithmetic are the shared ones.
    """
    from exondomaincompare.shared_gene_analysis.msa_coordinates import column_map, read_aligned_fasta
    if not msa.is_file():
        return {}
    out: Dict[str, Dict[int, int]] = {}
    for header, seq in read_aligned_fasta(msa):
        parts = [p for p in header.split("|") if p]
        pid = parts[2] if len(parts) > 2 else ""
        if pid:
            out[pid] = column_map(seq)
    return out


def build_index(threshold: int = bc.DEFAULT_NEAR_EDGE_THRESHOLD_AA,
                sources: Optional[Sources] = None) -> Dict[str, Any]:
    """One FGFR2 dataset as a shared coordinate-model index.

    ``sources`` defaults to the validated 30-species freeze, which is what every existing
    caller means.
    """
    src = sources or freeze_sources()
    features = _read_tsv(src.features)
    if not features:
        # No domain layer yet. The proteins and their validated labels do exist, and the
        # Gallery needs them to lay out its scopes, so the models are built from the truth
        # table with the domain-derived layers stated as pending rather than as empty.
        return _pre_cluster_index(src)
    interpro = _read_tsv(src.interpro)
    truth = _truth_by_protein(src.truth)
    combinations = _truth_by_combination(src.truth)
    columns = _msa_columns(src.msa)

    by_protein: Dict[str, List[Dict[str, str]]] = {}
    for r in features:
        pid = r.get("protein_id") or ""
        if pid:
            by_protein.setdefault(pid, []).append(r)

    # Which model is each species' primary reference, decided before the models are
    # built so the choice is one stated rule rather than a per-model accident.
    isoform_by_species: Dict[str, Dict[str, str]] = {}
    for pid, rows in by_protein.items():
        head = rows[0]
        isoform_by_species.setdefault(head.get("species") or "", {})[
            head.get("isoform") or ""] = pid
    reference_protein: Dict[str, str] = {}
    for species, available in isoform_by_species.items():
        for wanted in PRIMARY_REFERENCE_PREFERENCE:
            if available.get(wanted):
                reference_protein[species] = available[wanted]
                break
        else:
            reference_protein[species] = sorted(available.values())[0]

    models: List[Dict[str, Any]] = []
    for pid, rows in by_protein.items():
        head = rows[0]
        species = head.get("species") or ""
        isoform = head.get("isoform") or ""
        tr = truth.get(pid, {})
        length = _int(head.get("protein_length")) or 0
        domains = _representative_domains(rows, interpro, pid)
        exons, cassettes = _exons_and_cassette(rows, pid)
        coding_exons = [e for e in exons if not e["is_cassette_exon"]]
        boundaries = _boundaries(exons, domains, pid,
                                 head.get("transcript_id") or "", threshold)
        warnings = [r.get("feature_label") or "" for r in rows
                    if r.get("feature_type") == "warning"]
        column_of = columns.get(pid) or {}
        for b in boundaries:
            column = column_of.get(b["protein_position"])
            b["msa_column"] = column
            # Whether this boundary could be placed in the shared alignment at all.
            # A boundary that could not is reported as unmapped rather than being
            # quietly compared against boundaries that were mapped.
            b["msa_mapping_status"] = ("mapped" if column is not None
                                       else "unmapped_in_alignment")
        models.append({
            "schema_version": SCHEMA_VERSION,
            "species_id": species,
            "scientific_name": so.scientific_name(species),
            "gene_symbol": GENE_SYMBOL,
            # The model's own identity and role. A renderer is handed these; it never
            # infers which protein it is drawing from a file name or an array index.
            "model_id": mr.model_id(GENE_SYMBOL, species, isoform),
            "model_role": mr.validated_isoform_role(
                tr.get("final_isoform_label") or isoform),
            "is_primary_reference": reference_protein.get(species) == pid,
            "primary_reference_rule": (
                "Reference-database canonical isoform (IIIc) where the species has a "
                "IIIc model; otherwise its IIIb model."),
            "protein_id": pid,
            "transcript_id": head.get("transcript_id") or "",
            "protein_length": length,
            "coordinate_system": COORDINATE_SYSTEM,
            "status": "available",
            "pending_info": None,
            # FGFR2's unit of analysis is the protein: IIIb and IIIc are two real
            # proteins of one species, so the isoform is part of the model's identity.
            "isoform": isoform,
            "final_isoform_label": tr.get("final_isoform_label") or isoform,
            # What this model can actually support. A protein whose coding-exon
            # series could not be reconstructed still has a validated cassette and a
            # validated domain architecture, so it stays a real model — but the exon
            # and boundary figures do not exist for it, and saying so is the whole
            # point of carrying the status.
            "availability_status": ("cassette_only_no_exon_series"
                                    if not coding_exons else "model_available"),
            "unavailable_layers": ([] if coding_exons else
                                   ["exon_structure", "exon_domain_boundaries"]),
            "unavailable_reason": ("" if coding_exons else
                                   "no coding-exon series in the validated "
                                   "architecture; the cassette coordinate was "
                                   "validated separately"),
            "reconstruction_status": _exon_block_status(rows),
            "cds_reconstruction_status": (tr.get("CDS_reconstruction_status") or ""),
            "review_status": _review_status(tr),
            "validated_event": {
                "gene": GENE_SYMBOL,
                "event": "IIIb/IIIc mutually exclusive cassette exon",
                "isoform": tr.get("final_isoform_label") or isoform,
                "claim_status": tr.get("final_claim_status_after_rescue") or "",
                "readiness_class": tr.get("pre_interpro_readiness_class") or "",
                "rescue_decision": tr.get("rescue_decision") or "",
                "label_source": tr.get("final_label_source") or "",
                "validated": True,
            },
            "cassette_regions": cassettes,
            "tm_analysis": {
                "performed": True,
                "tm_region_count": sum(1 for r in rows
                                       if r.get("feature_type") == "transmembrane_pytmhmm"),
                "pending": False,
                "message": None,
            },
            "near_edge_threshold_aa": threshold,
            "exons": exons,
            "exon_boundaries": boundaries,
            "representative_domains": domains,
            "families_superfamilies": _families(rows, pid, interpro),
            "member_signatures": [],
            "functional_sites": [],
            "disorder_regions": [],
            "tm_regions": _tm_regions(rows, pid),
            "candidate_regions": [],
            "transcript_models": [],
            "n_transcript_models": 0,
            "alignment_mapping": {
                "available": bool(column_of),
                "n_columns": max(column_of.values()) if column_of else 0,
                "alignment_file": FULL_MSA.name if column_of else "",
            },
            "qc_warnings": warnings,
            "provenance": {
                "run_id": "example",
                "dataset": src.dataset_label,
                "gene_symbol": GENE_SYMBOL,
                "species_id": species,
                "isoform": isoform,
                "coordinate_system": COORDINATE_SYSTEM,
                "generated_by": "scripts/fgfr2/coordinate_model.py",
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "clade": so.clade_of(species),
                "taxonomic_group": so.taxon_group(species),
                "source_files": {
                    "architecture_features": _rel(src.features),
                    "interpro_signatures": _rel(src.interpro),
                    "truth_table": _rel(src.truth),
                    "full_length_msa": _rel(src.msa),
                },
                "derived_layers": {
                    "exon_boundaries": (
                        "Classified here with the shared generic classifier. The "
                        "frozen FGFR2 Boundary Consistency values are a separate "
                        "vocabulary and are not modified."),
                },
            },
        })

    order = {sid: i for i, sid in enumerate(
        so.order_species(m["species_id"] for m in models))}
    models.sort(key=lambda m: (order.get(m["species_id"], 999),
                               m.get("isoform") or "", m["protein_id"]))

    n_columns = max((m["alignment_mapping"]["n_columns"] for m in models), default=0)
    species_ids = list(dict.fromkeys(m["species_id"] for m in models))
    availability = _availability(models, species_ids, combinations)
    index = {
        "schema_version": SCHEMA_VERSION,
        "model_type": "ProteinCoordinateModelIndex",
        "run_id": src.run_id,
        "dataset": src.dataset_id,
        "gene_symbol": GENE_SYMBOL,
        "coordinate_system": COORDINATE_SYSTEM,
        "species_scope": species_ids,
        "n_models": len(models),
        "n_species": len(species_ids),
        "models": models,
        "availability": availability,
        "species_order": so.build_species_order(species_ids),
        "msa_coordinate_map": {
            "available": bool(n_columns),
            "alignment_file": src.msa.name,
            "n_columns": n_columns,
            "keyed_by": "protein_id",
            "reason": ("Full-length FGFR2 primary-protein alignment. FGFR2 has two "
                       "primary proteins per species (IIIb and IIIc), so columns are "
                       "keyed by protein rather than by species. A shared column means "
                       "the residues were aligned, not that they are equivalent."),
        },
        "comparative_primary": {
            "rule": ("One model per species: the species' primary reference, which is "
                     "its IIIc model where it has one (the reference-database "
                     "canonical isoform; RefSeq NM_000141, Ensembl ENST00000358487) "
                     "and otherwise its IIIb model. The figure states which."),
            "isoform_preference": list(PRIMARY_REFERENCE_PREFERENCE),
            "model_ids": {m["species_id"]: m["model_id"] for m in models
                          if m["is_primary_reference"]},
            "protein_ids": {m["species_id"]: m["protein_id"] for m in models
                            if m["is_primary_reference"]},
        },
        "provenance": {
            "generated_by": "scripts/fgfr2/coordinate_model.py",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "freeze_is_read_only": True,
            "note": ("Restates the frozen FGFR2 result tables in the shared contract. "
                     "No FGFR2 analysis is re-run and no validated value is changed."),
        },
    }
    index["boundary_dashboard"] = _boundary_dashboard(index)
    return index


def _boundary_dashboard(index: Dict[str, Any]) -> Dict[str, Any]:
    """The shared exon–domain boundary dashboard, over the primary reference panel.

    This is the *generic whole-protein* boundary analysis: every supported internal
    coding-exon boundary of every species, measured against the nearest
    representative domain edge. It is a different question from the validated
    IIIb/IIIc cassette-boundary analysis the freeze answers, and it makes no claim
    about the cassette; the two are reported side by side, never merged.

    The cross-species part is built from one model per species — each species'
    primary reference — because a comparable-boundary group counts species, and
    letting a species contribute both its isoforms would count it twice. Attaching
    the result to the canonical model index is what makes the Figure Gallery and the
    interactive Comparative Boundary Explorer read the same group ids, distances,
    classes and species order rather than two independently derived versions.
    """
    from exondomaincompare.shared_gene_analysis.boundary_dashboard import build_boundary_dashboard

    reduced = dict(index)
    reduced["models"] = [m for m in index["models"] if m["is_primary_reference"]]
    dashboard = build_boundary_dashboard(reduced)
    dashboard["comparative_panel"] = {
        "rule": (index.get("comparative_primary") or {}).get("rule", ""),
        "n_models_in_dataset": len(index["models"]),
        "n_models_compared": len(reduced["models"]),
        "model_ids": [m["model_id"] for m in reduced["models"]],
    }
    dashboard["scope_note"] = (
        "Supported internal coding-exon boundaries across the complete protein "
        "architecture. Separate from the validated FGFR2 IIIb/IIIc cassette-boundary "
        "analysis.")
    return dashboard


#: The layers a model only has once the cluster round-trip has returned annotation.
_CLUSTER_DERIVED_LAYERS = ("exon_structure", "domain_architecture",
                           "exon_domain_boundaries")

_PENDING_LAYER_REASON = (
    "InterProScan and pyTMHMM annotation has not been returned from the cluster yet, so "
    "this protein has no domain layer, no projected coding-exon series and no exon–domain "
    "boundaries.")


def _pre_cluster_index(src: Sources) -> Dict[str, Any]:
    """The models of a dataset whose cluster round-trip has not run.

    A protein selected by the pre-InterPro pipeline is a real model: it has an identity, a
    transcript, a length and a validated isoform label, all recorded in the closure's truth
    table. What it does not yet have are the three layers derived from cluster annotation.

    Building this instead of returning an empty index is what lets the Gallery show its
    scopes and its post-cluster sections as pending. Nothing is invented: every field comes
    from the truth table, and the missing layers are named as missing.
    """
    rows = [r for r in _read_tsv(src.truth) if r.get("protein_id") and r.get("species")]
    combinations = _truth_by_combination(src.truth)

    # The primary reference per species, by the same stated rule the full builder uses.
    by_species_isoform: Dict[str, Dict[str, str]] = {}
    for r in rows:
        by_species_isoform.setdefault(r["species"], {})[
            r.get("final_isoform_label") or r.get("isoform") or ""] = r["protein_id"]
    reference: Dict[str, str] = {}
    for species, available in by_species_isoform.items():
        for wanted in PRIMARY_REFERENCE_PREFERENCE:
            if available.get(wanted):
                reference[species] = available[wanted]
                break
        else:
            reference[species] = sorted(available.values())[0]

    models: List[Dict[str, Any]] = []
    for r in rows:
        species = r["species"]
        isoform = r.get("isoform") or ""
        label = r.get("final_isoform_label") or isoform
        models.append({
            "schema_version": SCHEMA_VERSION,
            "species_id": species,
            "scientific_name": so.scientific_name(species),
            "gene_symbol": GENE_SYMBOL,
            "model_id": mr.model_id(GENE_SYMBOL, species, isoform),
            "model_role": mr.validated_isoform_role(label),
            "is_primary_reference": reference.get(species) == r["protein_id"],
            "primary_reference_rule": (
                "Reference-database canonical isoform (IIIc) where the species has a "
                "IIIc model; otherwise its IIIb model."),
            "protein_id": r["protein_id"],
            "transcript_id": r.get("transcript_id") or "",
            "protein_length": _int(r.get("protein_length")) or 0,
            "coordinate_system": COORDINATE_SYSTEM,
            "status": "pending_cluster_annotation",
            "pending_info": {
                "waiting_for": ["interproscan", "pytmhmm"],
                "layers": list(_CLUSTER_DERIVED_LAYERS),
                "message": _PENDING_LAYER_REASON,
            },
            "isoform": isoform,
            "final_isoform_label": label,
            "availability_status": "pending_cluster_annotation",
            "unavailable_layers": list(_CLUSTER_DERIVED_LAYERS),
            "unavailable_reason": _PENDING_LAYER_REASON,
            "reconstruction_status": "pending_cluster_annotation",
            "cds_reconstruction_status": r.get("CDS_reconstruction_status") or "",
            "review_status": _review_status(r),
            "validated_event": {
                "gene": GENE_SYMBOL,
                "event": "IIIb/IIIc mutually exclusive cassette exon",
                "isoform": label,
                "claim_status": r.get("final_claim_status_after_rescue") or "",
                "readiness_class": r.get("pre_interpro_readiness_class") or "",
                "rescue_decision": r.get("rescue_decision") or "",
                "label_source": r.get("final_label_source") or "",
                "validated": True,
            },
            # Empty rather than absent: the analysis produced no rows because it has not
            # been performed. No coordinate, boundary or domain is fabricated here.
            "cassette_regions": [],
            "exons": [],
            "exon_boundaries": [],
            "representative_domains": [],
            "families_superfamilies": [],
            "member_signatures": [],
            "functional_sites": [],
            "disorder_regions": [],
            "tm_regions": [],
            "tm_analysis": {
                "performed": False, "tm_region_count": 0, "pending": True,
                "message": "pyTMHMM has not been run for this dataset yet.",
            },
            "alignment_mapping": {"n_columns": 0, "available": False},
            "provenance": {
                "dataset": src.dataset_label,
                "gene_symbol": GENE_SYMBOL,
                "species_id": species,
                "isoform": isoform,
                "coordinate_system": COORDINATE_SYSTEM,
                "generated_by": "scripts/fgfr2/coordinate_model.py",
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "clade": so.clade_of(species),
                "taxonomic_group": so.taxon_group(species),
                "source_files": {"truth_table": _rel(src.truth)},
                "derived_layers": {
                    "exon_boundaries": (
                        "Not derived: requires the cluster domain layer."),
                },
            },
        })

    order = {sid: i for i, sid in enumerate(
        so.order_species(m["species_id"] for m in models))}
    models.sort(key=lambda m: (order.get(m["species_id"], 999),
                               m.get("isoform") or "", m["protein_id"]))
    species_ids = list(dict.fromkeys(m["species_id"] for m in models))
    return {
        "schema_version": SCHEMA_VERSION,
        "model_type": "ProteinCoordinateModelIndex",
        "run_id": src.run_id,
        "dataset": src.dataset_id,
        "gene_symbol": GENE_SYMBOL,
        "coordinate_system": COORDINATE_SYSTEM,
        "cluster_annotation": "pending",
        "species_scope": species_ids,
        "n_models": len(models),
        "n_species": len(species_ids),
        "models": models,
        "availability": _availability(models, species_ids, combinations),
        "species_order": so.build_species_order(species_ids),
        "msa_coordinate_map": {
            "available": False, "alignment_file": src.msa.name, "n_columns": 0,
            "keyed_by": "protein_id",
            "reason": ("Alignment columns are mapped once the domain layer exists, so the "
                       "same builder produces both views from one pass."),
        },
        "comparative_primary": {
            "rule": ("One model per species: the species' primary reference, which is its "
                     "IIIc model where it has one and otherwise its IIIb model."),
            "isoform_preference": list(PRIMARY_REFERENCE_PREFERENCE),
            "model_ids": {m["species_id"]: m["model_id"] for m in models
                          if m["is_primary_reference"]},
            "protein_ids": {m["species_id"]: m["protein_id"] for m in models
                            if m["is_primary_reference"]},
        },
        "boundary_dashboard": {
            "available": False,
            "reason": _PENDING_LAYER_REASON,
            "scope_note": ("Supported internal coding-exon boundaries across the complete "
                           "protein architecture. Separate from the validated FGFR2 "
                           "IIIb/IIIc cassette-boundary analysis."),
        },
        "provenance": {
            "generated_by": "scripts/fgfr2/coordinate_model.py",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "freeze_is_read_only": True,
            "note": ("Pre-cluster models: identities and validated labels from the "
                     "closure truth table, with the cluster-derived layers stated as "
                     "pending. No domain, exon or boundary value is invented."),
        },
    }


def _availability(models: Sequence[Dict[str, Any]], species_ids: Sequence[str],
                  combinations: Dict[Tuple[str, str], Dict[str, str]],
                  ) -> Dict[str, Any]:
    """Which species–isoform combinations have a model, and why the others do not.

    The FGFR2 panel is 30 species and two validated isoforms, so 60 combinations
    are *expected*. 58 have an architecture model. The other two are not a gap in
    the data to be papered over: the freeze kept them for review only, with a
    stated reason, and a species that is missing one isoform model still belongs in
    the panel with the model it does have.
    """
    has_model = {(m["species_id"], m["isoform"]) for m in models}
    protein_of_model = {(m["species_id"], m["isoform"]): m["protein_id"]
                        for m in models}
    unavailable: List[Dict[str, Any]] = []
    for species in species_ids:
        for isoform in VALIDATED_ISOFORMS:
            if (species, isoform) in has_model:
                continue
            row = combinations.get((species, isoform), {})
            # Where the absent slot names the same protein as the species' available
            # slot, the isoform could not be told apart from sequence at all. That is
            # a different situation from a protein simply being missing, and a reader
            # deciding whether to trust the species needs to see which one it is.
            shared_with = next(
                (iso for iso in VALIDATED_ISOFORMS
                 if protein_of_model.get((species, iso)) == (row.get("protein_id") or "")
                 and row.get("protein_id")), "")
            unavailable.append({
                "shares_protein_with_isoform": shared_with,
                "species_id": species,
                "scientific_name": so.scientific_name(species),
                "isoform": isoform,
                "availability_status": "no_architecture_model",
                "protein_id": row.get("protein_id") or "",
                "claim_status": row.get("final_claim_status_after_rescue") or "",
                "readiness_class": row.get("pre_interpro_readiness_class") or "",
                "review_status": _review_status(row) if row else "not_in_truth_table",
                "omission_reason": (
                    (row.get("unresolved_reason_if_any")
                     or "not present in the validated architecture tables")
                    + (f"; the same protein is filed under this species' {shared_with} "
                       f"slot, so the two isoforms are not distinguishable by sequence "
                       f"here" if shared_with else "")),
                "species_still_represented": True,
            })

    by_species: Dict[str, List[Dict[str, Any]]] = {}
    for m in models:
        by_species.setdefault(m["species_id"], []).append(m)

    per_species = {}
    for species in species_ids:
        group = by_species.get(species) or []
        present = sorted(m["isoform"] for m in group)
        per_species[species] = {
            "n_models": len(present),
            "isoforms_available": present,
            "isoforms_unavailable": [i for i in VALIDATED_ISOFORMS if i not in present],
            "models": [{
                "model_id": m["model_id"],
                "isoform": m["isoform"],
                "model_role": m["model_role"],
                "is_primary_reference": m["is_primary_reference"],
                "availability_status": m["availability_status"],
                "unavailable_layers": m["unavailable_layers"],
                "reconstruction_status": m["reconstruction_status"],
                "review_status": m["review_status"],
            } for m in group],
        }

    reconstruction = {}
    for m in models:
        key = m["reconstruction_status"]
        reconstruction[key] = reconstruction.get(key, 0) + 1

    cassette_only = [m["model_id"] for m in models
                     if m["availability_status"] == "cassette_only_no_exon_series"]
    review_flagged = [m["model_id"] for m in models
                      if m["review_status"] != "no_review_flag"]

    return {
        "n_species": len(species_ids),
        "validated_isoforms": list(VALIDATED_ISOFORMS),
        "n_expected_combinations": len(species_ids) * len(VALIDATED_ISOFORMS),
        "n_models": len(models),
        "n_models_per_isoform": {
            isoform: sum(1 for m in models if m["isoform"] == isoform)
            for isoform in VALIDATED_ISOFORMS
        },
        "n_species_with_both_models": sum(
            1 for v in per_species.values() if v["n_models"] == 2),
        "n_species_with_one_model": sum(
            1 for v in per_species.values() if v["n_models"] == 1),
        "n_models_by_reconstruction_status": reconstruction,
        "cassette_only_models": cassette_only,
        "review_flagged_models": review_flagged,
        "unavailable_combinations": unavailable,
        "per_species": per_species,
        "explanation": (
            f"{len(species_ids)} species and {len(VALIDATED_ISOFORMS)} validated "
            f"isoforms give {len(species_ids) * len(VALIDATED_ISOFORMS)} expected "
            f"species-isoform combinations, of which {len(models)} have an "
            f"architecture model. The remaining {len(unavailable)} are combinations "
            f"the freeze kept for review only; every species stays in the panel with "
            f"the model or models it does have."),
    }


def primary_reference_models(index: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The one-model-per-species subset a comparative figure draws."""
    return [m for m in index.get("models") or [] if m.get("is_primary_reference")]


def models_by_species(index: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    """``{species_id: [model, …]}`` in canonical species order."""
    out: Dict[str, List[Dict[str, Any]]] = {}
    for m in index.get("models") or []:
        out.setdefault(m["species_id"], []).append(m)
    return out


INVENTORY_COLUMNS = [
    "species_id", "scientific_name", "model_id", "protein_id", "transcript_id",
    "isoform_label", "model_role", "is_primary_reference", "availability_status",
    "reconstruction_status", "review_status", "omission_reason",
]


def inventory_rows(index: Dict[str, Any]) -> List[Dict[str, str]]:
    """One row per expected species–isoform combination, present or not.

    A combination without a model gets a row too. A table that simply omitted it
    would make a species look like a one-isoform species, which is not what the
    freeze says about it.
    """
    rows: List[Dict[str, str]] = []
    for m in index.get("models") or []:
        rows.append({
            "species_id": m["species_id"],
            "scientific_name": m["scientific_name"],
            "model_id": m["model_id"],
            "protein_id": m["protein_id"],
            "transcript_id": m.get("transcript_id") or "",
            "isoform_label": m.get("final_isoform_label") or m.get("isoform") or "",
            "model_role": m["model_role"],
            "is_primary_reference": "true" if m["is_primary_reference"] else "false",
            "availability_status": m.get("availability_status") or "model_available",
            "reconstruction_status": m.get("reconstruction_status") or "",
            "review_status": m.get("review_status") or "",
            # Empty for a complete model. For a model whose exon series is missing
            # this states which layers it cannot support and why, so a reader of the
            # table never has to guess whether an absent figure is a bug.
            "omission_reason": m.get("unavailable_reason") or "",
        })
    for u in (index.get("availability") or {}).get("unavailable_combinations") or []:
        rows.append({
            "species_id": u["species_id"],
            "scientific_name": u["scientific_name"],
            "model_id": "",
            "protein_id": u.get("protein_id") or "",
            "transcript_id": "",
            "isoform_label": u["isoform"],
            "model_role": "",
            "is_primary_reference": "false",
            "availability_status": u["availability_status"],
            "reconstruction_status": "",
            "review_status": u.get("review_status") or "",
            "omission_reason": u["omission_reason"],
        })

    order = {sid: i for i, sid in enumerate(index.get("species_scope") or [])}
    rows.sort(key=lambda r: (order.get(r["species_id"], 999), r["isoform_label"]))
    return rows


def write_inventory(index: Dict[str, Any], outdir: Path) -> Path:
    """Write ``fgfr2_model_inventory.tsv``."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "fgfr2_model_inventory.tsv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=INVENTORY_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(inventory_rows(index))
    return path


def write_index(outdir: Path, sources: Optional[Sources] = None) -> Tuple[Path, Path]:
    """Write the coordinate model and the model inventory into ``outdir``."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    index = build_index(sources=sources)
    path = outdir / "protein_coordinate_model.json"
    path.write_text(json.dumps(index, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    return path, write_inventory(index, outdir / "tables")


def main(argv: Optional[List[str]] = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(ROOT / "results" / "derived" / "example"
                                        / "website_indices"),
                    help="Directory to write protein_coordinate_model.json into")
    args = ap.parse_args(argv)
    path, inventory = write_index(Path(args.out))
    index = json.loads(path.read_text(encoding="utf-8"))
    print(json.dumps({
        "model_index": _rel(path),
        "inventory": _rel(inventory),
        "msa_columns": index["msa_coordinate_map"]["n_columns"],
        "primary_reference_rows": len(index["comparative_primary"]["model_ids"]),
        "availability": {k: v for k, v in index["availability"].items()
                         if k not in ("unavailable_combinations", "per_species")},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
