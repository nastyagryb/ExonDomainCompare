#!/usr/bin/env python3
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
from exondomaincompare.contracts import freshness_verdict

# --------------------------------------------------------------------------- #
# Availability states
# --------------------------------------------------------------------------- #
#: The artefact is present, current and populated.
AVAILABLE = "available"
#: This run cannot have the artefact by design — a one-species run has no
#: cross-species comparison. Not a defect and not worth a retry.
NOT_APPLICABLE = "not_applicable"
#: The stage ran and honestly found nothing: no second isoform in this species, no
#: supported boundary. The data, not the software, is the limit.
SCIENTIFICALLY_UNAVAILABLE = "scientifically_unavailable"
#: The stage is marked complete but its output is not on disk. Always a defect, and
#: always worth a retry. This is the state the *Equus quagga* views should have shown.
TECHNICALLY_MISSING = "technically_missing"
#: The producing stage reported an error.
FAILED = "failed"
#: The producing stage has not run yet or is running.
PENDING = "pending"
#: The artefact exists but is older than an input it derives from, so what it says
#: describes a state the run has left behind.
STALE = "stale"

#: Every state, in the order of decreasing usability.
STATES = (AVAILABLE, SCIENTIFICALLY_UNAVAILABLE, NOT_APPLICABLE, STALE,
          TECHNICALLY_MISSING, FAILED, PENDING)

#: States that must never let a run count as finished.
BLOCKING_STATES = (TECHNICALLY_MISSING, FAILED, STALE)


# --------------------------------------------------------------------------- #
# Stage dependency order
# --------------------------------------------------------------------------- #
#: The FGFR2 pre-InterPro chain, upstream first. Each entry names the stage and the
#: directory under ``results/`` whose newest file dates the stage. Changing a stage
#: invalidates every stage after it — that is the whole point of keeping this ordered
#: rather than storing a per-stage "complete" flag that survives its own inputs.
FGFR2_STAGE_ORDER: Tuple[Tuple[str, str], ...] = (
    ("species_registry", "01_species_registry"),
    ("gene_models", "02_models"),
    ("transcript_selection", "03_selection_initial"),
    ("isoform_evidence", "04_isoform_evidence_v2_3_human_calibrated"),
    ("label_reconciliation", "05b_selection_with_isoforms_v2_7_marker_validated"),
    ("protein_export", "06_protein_export_v2_7_marker_validated"),
    ("cluster_input", "07_interpro_prepare_v2_7_marker_validated"),
    ("msa_boundary", "12_msa_boundary_robustness_pre_interpro"),
    ("closure", "13_final_pre_interpro_closure"),
    ("interproscan", "14_interproscan"),
    ("domain_architecture", "15_exon_domain_boundary_post_interpro"),
    ("final_analyses", "16_final_thesis_analyses"),
)

#: Derived layer, rebuilt from the closure and the post-cluster stages. Not a stage
#: folder: it lives beside ``results/``.
WEBSITE_INDICES = "website_indices"


@dataclass
class ViewRequirement:
    view: str
    label: str
    index: str
    #: Paths relative to the closure directory. Their absence is what separates
    #: "never written" from "written and empty".
    closure_inputs: Tuple[str, ...] = ()
    #: True when a run may legitimately lack this view. A one-species run has no
    #: cross-species views; a pre-cluster run has no domain views.
    needs_two_species: bool = False
    needs_cluster_outputs: bool = False
    #: Only FGFR2-style runs with a configured event region have cassette views.
    needs_event: bool = False


#: The views whose emptiness prompted this module. Each names the index the backend
#: serves and the closure artefacts the index is built from, so that a missing view can
#: say which file is missing rather than blaming the species.
FGFR2_VIEWS: Tuple[ViewRequirement, ...] = (
    ViewRequirement(
        "overview", "Overview", "run_index.json",
        closure_inputs=("final_pre_interpro_truth_table.tsv",)),
    ViewRequirement(
        "gene_explorer", "Gene explorer", "species_index.json",
        closure_inputs=("final_pre_interpro_truth_table.tsv",)),
    ViewRequirement(
        "exon_map", "Exon map", "coordinate_track_index.json",
        closure_inputs=("tables/figure3C_exon_to_protein_cassette_coordinate_map.tsv",)),
    ViewRequirement(
        "event_region", "Cassette", "cassette_residue_index.json",
        closure_inputs=("tables/figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv",
                        "tables/figure3B_IIIb_IIIc_cassette_amino_acid_motif_map.tsv"),
        needs_event=True),
    ViewRequirement(
        "msa", "MSA", "msa_index.json",
        closure_inputs=("MSA/final_fgfr2_full_length_protein_msa.aln.faa",)),
    ViewRequirement(
        "synteny", "Synteny", "synteny_locus_index.json"),
    ViewRequirement(
        "figure_gallery", "Figure gallery", "figure_index.json"),
    ViewRequirement(
        "domain_architecture", "Domain architecture", "species_domain_architecture.json",
        needs_cluster_outputs=True),
    ViewRequirement(
        "exon_domain_boundaries", "Exon–domain boundaries", "domain_architecture_index.json",
        needs_cluster_outputs=True),
    ViewRequirement(
        "boundary_consistency", "Boundary consistency", "boundary_consistency_summary.json",
        needs_cluster_outputs=True, needs_event=True),
    ViewRequirement(
        "downloads", "Data & downloads", "download_index.json"),
)

#: Views that must be available before a run may be called finished. Deliberately the
#: scientific core rather than every view: a run whose synteny is honestly unavailable is
#: still a finished run, but one without its coordinate map or alignment is not.
REQUIRED_FOR_READY: Tuple[str, ...] = (
    "overview", "gene_explorer", "exon_map", "msa", "downloads",
)
#: The views the cluster round-trip produces. Required because that round-trip is required
#: work — not because its output happens to be on disk.
#:
#: This was previously appended only ``if _cluster_outputs_present(run_dir)``, which asked
#: for the evidence exclusively in the case where it was already there. A pre-cluster FGFR2
#: run therefore had to satisfy only the pre-cluster views and was published as finished
#: while its domain and boundary layers were correctly reported as pending.
REQUIRED_FOR_READY_WITH_DOMAINS: Tuple[str, ...] = (
    "domain_architecture", "exon_domain_boundaries",
)


def _newest_mtime(path: Path, ignore: Sequence[str] = ()) -> Optional[float]:
    if path.is_file():
        return path.stat().st_mtime
    if not path.is_dir():
        return None
    newest: Optional[float] = None
    for child in path.rglob("*"):
        if not child.is_file():
            continue
        if any(part in ignore for part in child.parts):
            continue
        mtime = child.stat().st_mtime
        if newest is None or mtime > newest:
            newest = mtime
    return newest




def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _has_payload(data: Any) -> bool:
    if isinstance(data, list):
        return len(data) > 0
    if not isinstance(data, dict):
        return False
    if data.get("available") is False:
        return False
    # Sequence evidence can establish a cassette result without residue rows.
    collections = ("species", "rows", "alignments", "items", "figures", "cards",
                   "groups", "entries", "blocks", "sequence_evidence")
    present = [k for k in collections if k in data]
    for key in present:
        value = data.get(key)
        if isinstance(value, list) and value:
            return True
        if isinstance(value, dict):
            if key == "alignments":
                if any(bool(v.get("available")) for v in value.values()
                       if isinstance(v, dict)):
                    return True
            elif value:
                return True
    if present:
        # The index is a collection and every collection in it is empty.
        return False
    # A summary index: gate status, counters, KPIs: carries its content in its own
    # keys. Requiring an ``available`` flag it never had reported the run overview as
    # an absent scientific result.
    if data.get("available") is True:
        return True
    metadata = {"generated_at", "created_at", "run_id", "case_study", "title", "step",
                "note", "availability"}
    return any(k not in metadata for k in data)


# --------------------------------------------------------------------------- #
# Staleness
# --------------------------------------------------------------------------- #
@dataclass
class StageState:
    stage: str
    directory: str
    exists: bool
    newest_mtime: Optional[float]


def stage_states(run_dir: Path,
                 order: Sequence[Tuple[str, str]] = FGFR2_STAGE_ORDER,
                 ) -> List[StageState]:
    results = Path(run_dir) / "results"
    out: List[StageState] = []
    for stage, folder in order:
        path = results / folder
        newest = _newest_mtime(path)
        out.append(StageState(stage=stage, directory=folder,
                              exists=newest is not None, newest_mtime=newest))
    return out


def indices_are_stale(run_dir: Path, tolerance_s: float = 1.0,
                      ) -> Tuple[bool, str]:
    run_dir = Path(run_dir)
    indices = run_dir / WEBSITE_INDICES
    freshness = indices / "_freshness.json"
    if freshness.is_file():
        record = _read_json(freshness)
        if isinstance(record, dict):
            current, reason = freshness_verdict(run_dir, record)
            return not current, reason
    built = _newest_mtime(indices)
    if built is None:
        return False, "no website indices have been built"

    # Only stages the indices are actually derived from can make them stale. The
    # species registry is upstream of everything, so a registry rebuild that produced
    # identical models must not by itself condemn the indices.
    derived_from = {"closure", "interproscan", "domain_architecture", "final_analyses"}
    for state in stage_states(run_dir):
        if state.stage not in derived_from or state.newest_mtime is None:
            continue
        if state.newest_mtime > built + tolerance_s:
            return True, (f"stage {state.directory} was rewritten after the indices "
                          f"were built")
    return False, "legacy unversioned indices pass the mtime compatibility adapter"



@dataclass
class ViewState:
    view: str
    label: str
    state: str
    reason: str
    index: str
    missing_inputs: List[str] = field(default_factory=list)

    @property
    def available(self) -> bool:
        return self.state == AVAILABLE

    def as_dict(self) -> Dict[str, Any]:
        return {"view": self.view, "label": self.label, "state": self.state,
                "reason": self.reason, "index": self.index,
                "available": self.available,
                "missing_inputs": list(self.missing_inputs)}


def _cluster_outputs_present(run_dir: Path) -> bool:
    ips = (run_dir / "results" / "14_interproscan" / "primary" / "output")
    tm = (run_dir / "results" / "15_exon_domain_boundary_post_interpro"
          / "pytmhmm_primary" / "output")
    has_ips = ips.is_dir() and any(p.suffix == ".tsv" for p in ips.rglob("*"))
    has_tm = tm.is_dir() and any(p.is_file() for p in tm.rglob("*"))
    if not (has_ips and has_tm):
        return False
    # Only a positive mismatch removes the outputs. A run whose payload carries no
    # sequence digest cannot be checked either way, and treating "unverifiable" as
    # "superseded" would withhold results that are in fact current.
    try:
        from exondomaincompare.shared_gene_analysis.cluster_output_freshness import STALE, evaluate
        return evaluate(run_dir)["status"] != STALE
    except Exception:
        return True


def _tsv_has_rows(path: Path) -> bool:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            next(handle, None)            # header
            return next(handle, None) is not None
    except OSError:
        return False


def _fasta_has_records(path: Path) -> bool:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            return any(line.startswith(">") for line in handle)
    except OSError:
        return False


def _input_state(path: Path) -> str:
    if not path.is_file():
        return "missing"
    if path.suffix in (".faa", ".fasta", ".fa"):
        return "available" if _fasta_has_records(path) else "empty"
    if path.suffix == ".tsv":
        return "available" if _tsv_has_rows(path) else "empty"
    return "available" if path.stat().st_size > 0 else "empty"


def view_states(run_dir: Path, n_species: int = 0, has_event: bool = True,
                views: Sequence[ViewRequirement] = FGFR2_VIEWS,
                pre_interpro_complete: bool = True,
                ) -> List[ViewState]:
    run_dir = Path(run_dir)
    closure = run_dir / "results" / "13_final_pre_interpro_closure"
    indices = run_dir / WEBSITE_INDICES
    stale, stale_reason = indices_are_stale(run_dir)
    have_cluster = _cluster_outputs_present(run_dir)

    out: List[ViewState] = []
    for req in views:
        def state(s: str, reason: str, missing: Sequence[str] = ()) -> ViewState:
            return ViewState(req.view, req.label, s, reason, req.index, list(missing))

        if req.needs_event and not has_event:
            out.append(state(NOT_APPLICABLE,
                             "This gene has no configured event region, so there is no "
                             "cassette layer to show."))
            continue
        if req.needs_two_species and n_species < 2:
            out.append(state(NOT_APPLICABLE,
                             "A single-species run has nothing to compare across species."))
            continue
        if req.needs_cluster_outputs and not have_cluster:
            out.append(state(PENDING,
                             "Domain and transmembrane annotation has not been returned "
                             "from the cluster yet."))
            continue
        if not pre_interpro_complete:
            out.append(state(PENDING, "The local analysis has not finished yet."))
            continue

        # What the index is built from decides the diagnosis. An input that is missing
        # or empty explains the empty view; an input that is fine but an index that is
        # not points at the index build.
        missing = [rel for rel in req.closure_inputs
                   if _input_state(closure / rel) != "available"]
        empty_inputs = [rel for rel in req.closure_inputs
                        if _input_state(closure / rel) == "empty"]
        index_path = indices / req.index
        data = _read_json(index_path) if index_path.exists() else None

        if missing and len(missing) == len(req.closure_inputs) and req.closure_inputs:
            if empty_inputs and not [m for m in missing if m not in empty_inputs]:
                out.append(state(SCIENTIFICALLY_UNAVAILABLE,
                                 f"{req.label} produced no rows for this run's models.",
                                 missing))
            else:
                out.append(state(TECHNICALLY_MISSING,
                                 f"Expected {req.label.lower()} outputs were not "
                                 "generated. Retry local analysis.", missing))
            continue

        if data is None:
            out.append(state(TECHNICALLY_MISSING,
                             f"The {req.label.lower()} index was not written even though "
                             "its source data exist. Retry local analysis.",
                             [req.index]))
            continue

        if stale:
            out.append(state(STALE,
                             f"The {req.label.lower()} index is older than the analysis "
                             f"outputs it summarizes: {stale_reason}. Rebuild the "
                             "indices."))
            continue

        if not _has_payload(data):
            if req.needs_cluster_outputs:
                # The cluster returned annotation for these proteins, so a domain view
                # with nothing in it means the local post-cluster step did not derive
                # from it — a defect, not a species without domains.
                out.append(state(TECHNICALLY_MISSING,
                                 f"Cluster annotation was returned but the "
                                 f"{req.label.lower()} layer was not derived from it. "
                                 "Retry local analysis.", [req.index]))
                continue
            # Source data present, index current, still nothing in it. That is a real
            # scientific absence and may be stated as one.
            out.append(state(SCIENTIFICALLY_UNAVAILABLE,
                             f"No supported {req.label.lower()} result was recovered "
                             "for this run's models."))
            continue

        out.append(state(AVAILABLE, f"{req.label} is available."))
    return out


# --------------------------------------------------------------------------- #
# Readiness
# --------------------------------------------------------------------------- #
@dataclass
class Readiness:
    ready: bool
    reason: str
    blocking: List[ViewState] = field(default_factory=list)
    views: List[ViewState] = field(default_factory=list)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "ready": self.ready,
            "reason": self.reason,
            "blocking_views": [v.as_dict() for v in self.blocking],
            "views": {v.view: v.as_dict() for v in self.views},
        }


def models_run(run_dir: Path) -> bool:
    closure = Path(run_dir) / "results" / "13_final_pre_interpro_closure"
    return (closure / "final_pre_interpro_truth_table.tsv").is_file()


def readiness(run_dir: Path, n_species: int = 0, has_event: bool = True,
              pre_interpro_complete: bool = True) -> Readiness:
    run_dir = Path(run_dir)
    states = view_states(run_dir, n_species=n_species, has_event=has_event,
                         pre_interpro_complete=pre_interpro_complete)
    by_view = {s.view: s for s in states}

    # Applicable pending views block readiness; not-applicable views do not.
    required = list(REQUIRED_FOR_READY) + list(REQUIRED_FOR_READY_WITH_DOMAINS)

    blocking = [by_view[v] for v in required
                if v in by_view and by_view[v].state in BLOCKING_STATES]
    if blocking:
        first = blocking[0]
        return Readiness(False,
                         f"{first.label}: {first.reason}", blocking, states)

    pending = [by_view[v] for v in required
               if v in by_view and by_view[v].state == PENDING]
    if pending:
        return Readiness(False, f"{pending[0].label}: {pending[0].reason}",
                         pending, states)

    return Readiness(True, "Every required view has current data.", [], states)


__all__ = [
    "AVAILABLE", "NOT_APPLICABLE", "SCIENTIFICALLY_UNAVAILABLE", "TECHNICALLY_MISSING",
    "FAILED", "PENDING", "STALE", "STATES", "BLOCKING_STATES",
    "FGFR2_STAGE_ORDER", "FGFR2_VIEWS", "REQUIRED_FOR_READY",
    "REQUIRED_FOR_READY_WITH_DOMAINS", "WEBSITE_INDICES",
    "ViewRequirement", "ViewState", "StageState", "Readiness",
    "stage_states", "indices_are_stale", "view_states", "readiness", "models_run",
]
