#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

TWO_SPECIES_RUN = ROOT / "runs" / "2026-07-26_2157_fgfr1_gallus_mus_core_pilot"
GALLUS_RUN = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"
TP53_RUN = ROOT / "runs" / "2026-07-21_1436_custom_run"

pytestmark = pytest.mark.skipif(
    not TWO_SPECIES_RUN.is_dir(),
    reason=f"reference two-species run not present: {TWO_SPECIES_RUN.name}")


def _tsv(path: Path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def core() -> Path:
    return TWO_SPECIES_RUN / "results" / "core_gene_analysis"


@pytest.fixture(scope="module")
def models():
    from exondomaincompare.shared_gene_analysis.protein_coordinate_model import build_models_for_run
    idx = build_models_for_run(TWO_SPECIES_RUN)
    return {m["species_id"]: m for m in idx["models"]}, idx


# --------------------------------------------------------------------------- #
# 1-2. cluster input and returned output cover every species
# --------------------------------------------------------------------------- #
def test_every_selected_species_appears_once_in_the_cluster_input():
    fasta = (TWO_SPECIES_RUN / "results/14_interproscan/primary/input"
             / "final_pre_interpro_proteins_primary.faa")
    headers = [l.strip()[1:] for l in fasta.read_text().splitlines()
               if l.startswith(">")]
    species = [h.rpartition("|")[2] for h in headers]
    wanted = [s.strip() for s in
              (TWO_SPECIES_RUN / "species_list.txt").read_text().split() if s.strip()]

    assert sorted(species) == sorted(wanted), (
        f"cluster input must hold exactly one primary per selected species; "
        f"selected {wanted}, submitted {species}")
    assert len(set(headers)) == len(headers), f"duplicate FASTA headers: {headers}"


def test_every_returned_sequence_maps_to_exactly_one_species(core: Path):
    ips = _tsv(TWO_SPECIES_RUN
               / "results/14_interproscan/primary/output/input.fasta.tsv")
    returned = {r[list(r.keys())[0]] for r in ips} if ips else set()
    if not returned:  # header-less TSV: read the first column directly
        raw = (TWO_SPECIES_RUN
               / "results/14_interproscan/primary/output/input.fasta.tsv").read_text()
        returned = {l.split("\t")[0] for l in raw.splitlines() if l.strip()}

    from exondomaincompare.framework.primary_resolution import resolve_primaries
    primaries = resolve_primaries(core)
    by_protein = {v["protein_id"]: sid for sid, v in primaries.items()}

    assert len(by_protein) == len(primaries), (
        "two species must not resolve to the same protein: "
        f"{ {s: v['protein_id'] for s, v in primaries.items()} }")
    for pid in returned:
        assert pid in by_protein, (
            f"returned sequence {pid} maps to no species; known primaries "
            f"{sorted(by_protein)}")


# --------------------------------------------------------------------------- #
# 3-4. both species processed, and the primary agrees with what was analysed
# --------------------------------------------------------------------------- #
def test_the_resolved_primary_is_the_protein_the_cluster_analysed(core: Path):
    from exondomaincompare.framework.primary_resolution import resolve_primaries
    fasta = (TWO_SPECIES_RUN / "results/14_interproscan/primary/input"
             / "final_pre_interpro_proteins_primary.faa")
    submitted = {}
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            head = line[1:].strip()
            submitted[head.rpartition("|")[2]] = head.split()[0]

    for sid, v in resolve_primaries(core).items():
        assert v["protein_id"] == submitted.get(sid), (
            f"{sid}: the coordinate model would be built for {v['protein_id']} but the "
            f"cluster analysed {submitted.get(sid)}. Every domain and boundary lookup "
            f"filters on species AND protein, so a mismatch silently empties the species.")


def test_both_species_have_real_domain_and_boundary_results(models):
    by_species, _ = models
    assert set(by_species) == {"gallus_gallus", "mus_musculus"}
    for sid, m in by_species.items():
        assert m["status"] == "available", f"{sid} is {m['status']}"
        assert m["representative_domains"], f"{sid} has no representative domains"
        assert m["exon_boundaries"], f"{sid} has no exon boundaries"
        assert m["tm_regions"], f"{sid} has no TM region (FGFR1 is a receptor)"


def test_mus_musculus_reproduces_its_real_architecture(models):
    by_species, _ = models
    mus, gallus = by_species["mus_musculus"], by_species["gallus_gallus"]

    assert mus["protein_id"] != gallus["protein_id"]
    assert mus["protein_length"] != gallus["protein_length"], (
        "identical lengths would suggest one species' data was reused for the other")
    assert len(mus["representative_domains"]) == len(gallus["representative_domains"])

    mus_spans = {(d["start"], d["end"]) for d in mus["representative_domains"]}
    gallus_spans = {(d["start"], d["end"]) for d in gallus["representative_domains"]}
    assert not (mus_spans & gallus_spans), (
        f"Mus and Gallus share domain spans {mus_spans & gallus_spans}, which means one "
        f"species' annotation was applied to the other")


# --------------------------------------------------------------------------- #
# 5. no first-species-only logic in the post-cluster builders
# --------------------------------------------------------------------------- #
def test_the_coordinate_model_does_not_guess_a_primary_alphabetically():
    src = (ROOT / "src/exondomaincompare/shared_gene_analysis"
           / "protein_coordinate_model.py").read_text()
    assert "or sorted(proteins)[0]" not in src.split("primary_by_species.get(sp)")[0], (
        "an alphabetical fallback must not precede the species-scoped resolution")
    assert "primary_by_species" in src, (
        "the model must resolve the primary per species, not from a run-level set")


def test_a_species_without_a_primary_raises_instead_of_being_guessed(tmp_path: Path):
    from exondomaincompare.framework.primary_resolution import PrimaryResolutionError, resolve_primaries
    core = tmp_path / "core"
    core.mkdir()
    (core / "protein_isoform_index.tsv").write_text(
        "species_id\tprotein_id\ttranscript_id\tprotein_length\tprimary_status\n"
        "gallus_gallus\tNP_1.1\tNM_1.1\t100\tprimary\n"
        "mus_musculus\tNP_2.1\tNM_2.1\t110\talternative\n"
        "mus_musculus\tNP_3.1\tNM_3.1\t120\talternative\n")
    with pytest.raises(PrimaryResolutionError) as err:
        resolve_primaries(core)
    assert "mus_musculus" in str(err.value)
    assert "Refusing to guess" in str(err.value)


# --------------------------------------------------------------------------- #
# 6. the selection evidence table is species-aware
# --------------------------------------------------------------------------- #
def test_the_selection_evidence_table_carries_a_species_column():
    from exondomaincompare.framework.primary_selection import _TSV_FIELDS
    assert "species_id" in _TSV_FIELDS
    assert "species_primary" in _TSV_FIELDS


def test_the_writer_marks_one_primary_per_species(tmp_path: Path):
    from exondomaincompare.framework.primary_selection import write_selection_evidence
    report = {
        "primary_protein_id": "NP_g.1",
        "species_primaries": [
            {"species_id": "gallus_gallus", "primary_protein_id": "NP_g.1"},
            {"species_id": "mus_musculus", "primary_protein_id": "NP_m.1"},
        ],
        "proteins": [
            {"protein_id": "NP_g.1", "species_id": "gallus_gallus", "transcript_id": "a"},
            {"protein_id": "NP_m.1", "species_id": "mus_musculus", "transcript_id": "b"},
            {"protein_id": "NP_m.2", "species_id": "mus_musculus", "transcript_id": "c"},
        ],
    }
    tsv, js = tmp_path / "e.tsv", tmp_path / "e.json"
    write_selection_evidence(report, tsv, js)
    rows = _tsv(tsv)

    marked = {r["species_id"] for r in rows if r["species_primary"] == "true"}
    assert marked == {"gallus_gallus", "mus_musculus"}, (
        "each species needs its own primary; marking only the run-level primary is what "
        "left every Mus row false")
    assert sum(1 for r in rows if r["selected_primary"] == "true") == 1, (
        "the run-level primary stays single, for single-species consumers")


# --------------------------------------------------------------------------- #
# 7-8. per-species completion contract and run-level aggregation
# --------------------------------------------------------------------------- #
def test_every_species_gets_a_completion_object():
    from exondomaincompare.framework.species_completion import REQUIRED_ANALYSES, build_species_completion
    completion = build_species_completion(TWO_SPECIES_RUN)
    assert set(completion) == {"gallus_gallus", "mus_musculus"}
    for sid, rec in completion.items():
        assert rec["primary_protein"], f"{sid} has no primary protein recorded"
        for key in REQUIRED_ANALYSES:
            assert key in rec, f"{sid} is missing the {key} state"


def test_a_complete_two_species_run_is_results_ready():
    from exondomaincompare.framework.species_completion import aggregate_run_status, build_species_completion
    status, reasons = aggregate_run_status(build_species_completion(TWO_SPECIES_RUN))
    assert status == "results_ready", f"{status}: {reasons}"


def test_one_incomplete_species_makes_the_run_partial(tmp_path: Path):
    from exondomaincompare.framework.species_completion import aggregate_run_status, build_species_completion

    sandbox = tmp_path / "run"
    shutil.copytree(TWO_SPECIES_RUN, sandbox)
    model_path = sandbox / "website_indices/generic/protein_coordinate_model.json"
    model = json.loads(model_path.read_text())
    for m in model["models"]:
        if m["species_id"] == "mus_musculus":
            m["status"] = "pending_cluster"
            m["representative_domains"] = []
            m["tm_regions"] = []
            m["tm_analysis"] = {"performed": False, "pending": True}
    model_path.write_text(json.dumps(model))

    status, reasons = aggregate_run_status(build_species_completion(sandbox))
    assert status == "post_cluster_partial", (
        "a run with one species lacking domains must not advertise Results ready")
    assert any("mus_musculus" in r for r in reasons), reasons


def test_an_analysis_that_ran_and_found_nothing_is_not_incomplete():
    if not TP53_RUN.is_dir():
        pytest.skip("TP53 regression run not present")
    from exondomaincompare.framework.species_completion import aggregate_run_status, build_species_completion
    completion = build_species_completion(TP53_RUN)
    for sid, rec in completion.items():
        assert rec["pytmhmm"] == "available", (
            f"{sid}: pyTMHMM ran and predicted no TM region, which must not be "
            f"reported as missing")
    status, _ = aggregate_run_status(completion)
    assert status == "results_ready"


# --------------------------------------------------------------------------- #
# 9-10. single-species regressions still hold
# --------------------------------------------------------------------------- #
def test_the_gallus_single_species_reference_is_unchanged():
    if not GALLUS_RUN.is_dir():
        pytest.skip("Gallus reference run not present")
    from exondomaincompare.shared_gene_analysis.protein_coordinate_model import build_models_for_run
    idx = build_models_for_run(GALLUS_RUN)
    assert idx["n_models"] == 1
    m = idx["models"][0]
    assert m["protein_id"] == "NP_990841.2"
    assert m["transcript_id"] == "NM_205510.2"
    assert m["protein_length"] == 817
    assert len(m["exons"]) == 17
    assert len(m["exon_boundaries"]) == 16
    assert len(m["representative_domains"]) == 4
    assert len(m["tm_regions"]) == 1

    from collections import Counter
    classes = Counter(b["boundary_class"] for b in m["exon_boundaries"])
    assert classes["near_domain_edge"] == 6
    assert classes["inside_domain"] == 8
    assert classes["outside_annotated_domains"] == 2
