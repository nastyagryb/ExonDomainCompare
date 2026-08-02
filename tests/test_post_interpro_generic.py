"""Correctness tests for the generic post-InterPro data flow.

Covers the bugs fixed in the Phase-1 post-InterPro audit:
  * protein-ID normalisation between FASTA / InterProScan / pyTMHMM
  * domain / TM attribution to the correct species and protein
  * generic exon-domain boundary classification (exact/near/inside/outside/unknown)
  * generic Domain Architecture index shape (no FGFR2 IIIb/IIIc panels)
  * preservation of the FGFR2 specialization (validated example keeps IIIb/IIIc)

The real two-species FGFR1 run (Gallus + Mus, real InterProScan/pyTMHMM) is used
as an end-to-end fixture when present.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from framework import run_core_gene_analysis as rcga  # noqa: E402

FGFR1_RUN = PROJECT_ROOT / "runs" / "2026-07-14_2153_fgfr1_gallus_mus_twospecies"
FGFR1_INDICES = FGFR1_RUN / "website_indices"
FGFR2_EXAMPLE = (PROJECT_ROOT / "results" / "final_30_until_interpro_prepare"
                 / "13_final_pre_interpro_closure" / "website_indices")

GENERIC_CATEGORIES = {"exact_edge", "near_edge", "inside_domain", "outside_domain", "unknown"}


# --------------------------------------------------------------------------- #
# protein-ID normalisation
# --------------------------------------------------------------------------- #
def test_norm_acc_strips_refseq_version():
    assert rcga._norm_acc("NP_990841.2") == "NP_990841"
    assert rcga._norm_acc("XP_024998538.1") == "XP_024998538"
    assert rcga._norm_acc("NP_990841") == "NP_990841"
    assert rcga._norm_acc("  NP_034336.2  ") == "NP_034336"


def test_protein_species_map(tmp_path):
    core = tmp_path / "core"
    core.mkdir()
    with open(core / "protein_isoform_index.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species_id", "protein_id", "primary_status"])
        w.writerow(["gallus_gallus", "NP_990841.2", "primary"])
        w.writerow(["mus_musculus", "NP_034336.2", "primary"])
        w.writerow(["gallus_gallus", "XP_024998538.1", "alternative"])
    m = rcga._protein_species_map(core)
    assert m["NP_990841.2"] == "gallus_gallus"
    assert m["NP_990841"] == "gallus_gallus"          # version-stripped key
    assert m["NP_034336.2"] == "mus_musculus"
    assert m["XP_024998538.1"] == "gallus_gallus"


def test_species_assignment_does_not_collapse_to_one_species(tmp_path):
    """Regression: a domain row for gallus's protein must not be labelled mus."""
    core = tmp_path / "core"
    core.mkdir()
    with open(core / "protein_isoform_index.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species_id", "protein_id", "primary_status"])
        w.writerow(["gallus_gallus", "NP_990841.2", "primary"])
        w.writerow(["mus_musculus", "NP_034336.2", "primary"])
    prot_species = rcga._protein_species_map(core)
    # cluster tools echo only the bare accession -> parser leaves species empty
    rows = [
        {"protein_id": "NP_990841.2", "species_id": ""},
        {"protein_id": "NP_034336.2", "species_id": ""},
    ]
    default_species = next(iter(dict.fromkeys(prot_species.values())), "")
    for r in rows:
        r["species_id"] = (r.get("species_id")
                           or prot_species.get(r["protein_id"])
                           or prot_species.get(rcga._norm_acc(r["protein_id"]))
                           or default_species)
    assert rows[0]["species_id"] == "gallus_gallus"
    assert rows[1]["species_id"] == "mus_musculus"


# --------------------------------------------------------------------------- #
# boundary classification
# --------------------------------------------------------------------------- #
def test_classify_boundary_unknown_when_no_domains():
    nearest, edge, cat, dist = rcga._classify_boundary(120, [], 5)
    assert cat == "unknown"
    assert dist is None


def test_classify_boundary_exact_near_inside_outside():
    dom = [{"domain_name": "IgI", "domain_id": "d1", "start_aa": 100, "end_aa": 200}]
    # exact edge (boundary lands on the domain end)
    assert rcga._classify_boundary(200, dom, 5)[2] == "exact_edge"
    # near edge (within threshold)
    assert rcga._classify_boundary(203, dom, 5)[2] == "near_edge"
    # inside the domain
    assert rcga._classify_boundary(150, dom, 5)[2] == "inside_domain"
    # clearly outside annotated domain space
    assert rcga._classify_boundary(500, dom, 5)[2] == "outside_domain"


# --------------------------------------------------------------------------- #
# end-to-end: real FGFR1 two-species run
# --------------------------------------------------------------------------- #
requires_run = pytest.mark.skipif(
    not (FGFR1_INDICES / "domain_architecture_index.json").is_file(),
    reason="FGFR1 two-species post-InterPro run not present",
)


@requires_run
def test_domain_index_is_generic_and_species_correct():
    d = json.loads((FGFR1_INDICES / "domain_architecture_index.json").read_text())
    assert d.get("mode") == "generic"
    assert d.get("available") is True
    species = {s["species"]: s for s in d["species"]}
    assert "gallus_gallus" in species and "mus_musculus" in species
    # No FGFR2 IIIb/IIIc cassette panels anywhere in a generic index.
    for s in d["species"]:
        assert "panels" not in s
    # gallus's primary protein carries its own domains (not orphaned onto mus).
    gallus = {p["protein_id"]: p for p in species["gallus_gallus"]["proteins"]}
    assert "NP_990841.2" in gallus
    assert len(gallus["NP_990841.2"]["domains"]) > 0
    assert gallus["NP_990841.2"]["role"] == "primary"
    # NP_990841.2 must appear ONLY under gallus (no species mixing).
    mus = {p["protein_id"] for p in species["mus_musculus"]["proteins"]}
    assert "NP_990841.2" not in mus


@requires_run
def test_tm_features_attributed_per_species():
    rows = list(csv.DictReader(
        open(FGFR1_RUN / "results" / "core_gene_analysis" / "tm_features.tsv",
             encoding="utf-8"), delimiter="\t"))
    by_pid = {r["protein_id"]: r["species_id"] for r in rows}
    assert by_pid.get("NP_990841.2") == "gallus_gallus"
    assert by_pid.get("NP_034336.2") == "mus_musculus"


@requires_run
def test_boundary_index_numeric_and_primary_scoped():
    d = json.loads((FGFR1_INDICES / "exon_domain_boundaries_index.json").read_text())
    assert d.get("protein_scope") == "primary_only"
    assert d.get("near_edge_threshold_aa") == 5
    assert set(d.get("species_scope", [])) == {"gallus_gallus", "mus_musculus"}
    for p in d["proteins"]:
        for b in p["boundaries"]:
            # every boundary distance is numeric (no "— aa")
            assert isinstance(b["absolute_distance_aa"], int)
            assert b["category"] in GENERIC_CATEGORIES
    # counts use the explicit generic vocabulary only
    assert set(d["category_counts"]).issuperset({"exact_edge", "near_edge"})
    assert set(d["category_counts"]).issubset(GENERIC_CATEGORIES)


@requires_run
def test_no_species_mixing_in_boundaries():
    d = json.loads((FGFR1_INDICES / "exon_domain_boundaries_index.json").read_text())
    for p in d["proteins"]:
        if p["protein_id"] == "NP_990841.2":
            assert p["species_id"] == "gallus_gallus"
        if p["protein_id"] == "NP_034336.2":
            assert p["species_id"] == "mus_musculus"


# --------------------------------------------------------------------------- #
# FGFR2 specialization preserved
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not (FGFR2_EXAMPLE / "species_domain_architecture.json").is_file(),
                    reason="FGFR2 example freeze not present")
def test_fgfr2_example_keeps_cassette_panels():
    # The validated FGFR2 Domain Architecture is driven by species_domain_architecture.json
    # (served to the component via legacy_fgfr2_indices). It must still expose the IIIb/IIIc
    # cassette panels, i.e. the frontend cassette-detection must classify it as FGFR2, not generic.
    d = json.loads((FGFR2_EXAMPLE / "species_domain_architecture.json").read_text())
    has_panels = any(
        (s.get("panels") or {}).get("IIIb") or (s.get("panels") or {}).get("IIIc")
        for s in d.get("species", []))
    assert has_panels
    assert d.get("mode") != "generic"
