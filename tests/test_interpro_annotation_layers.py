"""Phase 1B tests: generic InterPro annotation semantics and domain layering.

Covers:
  * InterPro entry-type -> layer mapping (domain / family / feature / raw)
  * integrated vs unintegrated hits
  * duplicate / overlap collapse into representative domains
  * representative-domain selection (member-DB consensus rule)
  * separation of domains, families, features and (external) topology
  * TP53 is rendered from its real p53 entry types, with NO FGFR keywords
  * no alternative-protein selector when only the primary is annotated
  * Gene-Explorer summary and global Boundary page use identical boundary rows
  * FGFR2 specialization (IIIb/IIIc panels) is untouched by the generic layer
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from framework import interpro_annotations as ia  # noqa: E402

FGFR1_RUN = PROJECT_ROOT / "runs" / "2026-07-14_2153_fgfr1_gallus_mus_twospecies"
TP53_RUN = PROJECT_ROOT / "runs" / "2026-07-16_1642_tp53_human_core_pilot"
FGFR2_EXAMPLE = (PROJECT_ROOT / "results" / "final_30_until_interpro_prepare"
                 / "13_final_pre_interpro_closure" / "website_indices")

DOMAIN_INDEX = "generic/domain_architecture_index.json"
BOUNDARY_INDEX = "generic/exon_domain_boundaries_index.json"

FGFR_KEYWORDS = ("ig-like", "ig_like", "kinase", "fgfr", "fgf_rcpt", "tyr_kinase")


def _load(run: Path, rel: str):
    p = run / "website_indices" / rel
    if not p.is_file():
        pytest.skip(f"index not present: {p}")
    return json.loads(p.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
# entry-type -> layer mapping
# --------------------------------------------------------------------------- #
def test_layer_for_entry_types():
    assert ia.layer_for("DOMAIN", True) == "domain"
    assert ia.layer_for("REPEAT", True) == "domain"
    assert ia.layer_for("FAMILY", True) == "family"
    assert ia.layer_for("HOMOLOGOUS_SUPERFAMILY", True) == "family"
    assert ia.layer_for("ACTIVE_SITE", True) == "feature"
    assert ia.layer_for("BINDING_SITE", True) == "feature"
    assert ia.layer_for("CONSERVED_SITE", True) == "feature"


def test_layer_for_integration_and_disorder():
    # a DOMAIN-typed hit that is NOT integrated must not be promoted to a domain
    assert ia.layer_for("DOMAIN", False) == "raw"
    # unintegrated hit with no type is raw ...
    assert ia.layer_for("", False, "PFAM") == "raw"
    # ... unless it is a recognised feature predictor (MobiDB-lite disorder)
    assert ia.layer_for("", False, "MOBIDB_LITE") == "feature"


# --------------------------------------------------------------------------- #
# duplicate / overlap collapse + representative selection
# --------------------------------------------------------------------------- #
def _hit(acc, name, typ, s, e, db, layer="domain"):
    return {"interpro_accession": acc, "interpro_name": name, "interpro_type": typ,
            "start": s, "end": e, "member_database": db, "signature_accession": f"{db}:{acc}",
            "signature_name": name, "score_or_evalue": "", "layer": layer}


def test_representative_domains_collapse_same_region():
    # four different InterPro DOMAIN entries all describe ONE region -> 1 domain
    hits = [
        _hit("IPR007110", "Ig-like_dom", "DOMAIN", 33, 118, "SMART"),
        _hit("IPR003599", "Ig_sub", "DOMAIN", 39, 118, "SMART"),
        _hit("IPR013151", "Immunoglobulin_dom", "DOMAIN", 41, 112, "PFAM"),
        _hit("IPR003598", "Ig_sub2", "DOMAIN", 45, 107, "SMART"),
    ]
    reps = ia.representative_domains(hits)
    assert len(reps) == 1
    assert reps[0]["start_aa"] == 33 and reps[0]["end_aa"] == 118
    # the collapsed region keeps every contributing InterPro entry for provenance
    accs = {s["interpro_accession"] for s in reps[0]["supporting_interpro"]}
    assert accs == {"IPR007110", "IPR003599", "IPR013151", "IPR003598"}


def test_representative_selection_prefers_member_db_consensus():
    # same region: entry A supported by 3 member DBs, entry B by 1 -> A wins
    hits = [
        _hit("IPR-A", "Consensus", "DOMAIN", 100, 200, "PFAM"),
        _hit("IPR-A", "Consensus", "DOMAIN", 100, 200, "SMART"),
        _hit("IPR-A", "Consensus", "DOMAIN", 100, 200, "CDD"),
        _hit("IPR-B", "Single", "DOMAIN", 95, 205, "PANTHER"),
    ]
    reps = ia.representative_domains(hits)
    assert len(reps) == 1
    assert reps[0]["interpro_accession"] == "IPR-A"


def test_layers_separate_domain_family_feature():
    hits = [
        _hit("IPR007110", "Ig-like_dom", "DOMAIN", 33, 118, "SMART", "domain"),
        _hit("IPR016248", "FGF_rcpt_fam", "FAMILY", 6, 811, "PIRSF", "family"),
        _hit("IPR008266", "Tyr_kinase_AS", "ACTIVE_SITE", 615, 627, "PROSITE", "feature"),
    ]
    assert len(ia.representative_domains(hits)) == 1
    assert len(ia.family_annotations(hits)) == 1
    assert len(ia.feature_annotations(hits)) == 1


# --------------------------------------------------------------------------- #
# real FGFR1 run: non-redundant representative domains
# --------------------------------------------------------------------------- #
def test_fgfr1_representative_domain_layer_is_nonredundant():
    ips = FGFR1_RUN / "results" / "14_interproscan" / "primary" / "output"
    if not ips.is_dir():
        pytest.skip("FGFR1 InterProScan output not present")
    hits = ia.load_normalized_annotations(ips)
    by_prot = {}
    for h in hits:
        by_prot.setdefault(h["protein_accession"], []).append(h)
    for pid, ph in by_prot.items():
        reps = ia.representative_domains(ph)
        # FGFR1 has far fewer representative domains than raw hits (3 Ig + 1 kinase)
        assert 3 <= len(reps) <= 6, (pid, len(reps))
        assert len(reps) < len([h for h in ph if h["layer"] == "domain"])


def test_fgfr1_domain_index_counts_are_representative():
    idx = _load(FGFR1_RUN, DOMAIN_INDEX)
    assert idx["mode"] == "generic"
    for sp in idx["species"]:
        prim = next(p for p in sp["proteins"] if p["role"] == "primary")
        assert 3 <= len(prim["domains"]) <= 6
        assert prim["families"], "family layer must be populated separately"
        assert prim["annotated"] is True
        # domains carry an InterPro DOMAIN/REPEAT type only
        assert all(d["interpro_type"] in ("DOMAIN", "REPEAT") for d in prim["domains"])


# --------------------------------------------------------------------------- #
# TP53: real p53 entry types, NO FGFR keywords
# --------------------------------------------------------------------------- #
def test_tp53_domains_are_p53_not_fgfr():
    idx = _load(TP53_RUN, DOMAIN_INDEX)
    sp = idx["species"][0]
    prim = next(p for p in sp["proteins"] if p["role"] == "primary")
    names = " ".join((d["interpro_name"] or "").lower() for d in prim["domains"])
    assert "p53" in names
    for kw in FGFR_KEYWORDS:
        assert kw not in names, f"generic TP53 domains must not contain FGFR keyword {kw!r}"
    assert prim["domains"], "TP53 must show real representative domains"


def test_tp53_only_primary_is_annotated_selector():
    idx = _load(TP53_RUN, DOMAIN_INDEX)
    sp = idx["species"][0]
    # exactly one annotated protein -> UI renders a static label, not a dropdown
    assert sp["n_annotated_proteins"] == 1
    annotated = [p for p in sp["proteins"] if p["annotated"]]
    assert len(annotated) == 1 and annotated[0]["role"] == "primary"
    # the many unsubmitted isoforms must be present but not annotated
    assert any(p["role"] == "alternative" and not p["annotated"] for p in sp["proteins"])


def test_tp53_boundary_summary_populated_and_representative_layer():
    b = _load(TP53_RUN, BOUNDARY_INDEX)
    assert b["available"] is True
    assert b["domain_layer"] == "representative_domain"
    proteins = b["proteins"]
    assert proteins, "TP53 must expose a species-specific boundary protein"
    # every boundary distance is numeric and classified generically
    for p in proteins:
        for bd in p["boundaries"]:
            assert isinstance(bd["boundary_position_aa"], int)
            assert bd["classification"] in ("exact_edge", "near_edge", "inside_domain",
                                            "outside_domain", "unknown")


# --------------------------------------------------------------------------- #
# unify: identical rows feed local summary and global page
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("run", [FGFR1_RUN, TP53_RUN])
def test_boundary_rows_identical_scope_local_and_global(run):
    b = _load(run, BOUNDARY_INDEX)
    # the same `proteins` array is what both the Gene-Explorer summary and the
    # global Boundary page group by species; counts must be derivable from it.
    total = sum(len(p["boundaries"]) for p in b["proteins"])
    assert total == b["n_boundaries"]
    recomputed = {}
    for p in b["proteins"]:
        for bd in p["boundaries"]:
            recomputed[bd["classification"]] = recomputed.get(bd["classification"], 0) + 1
    for cat, n in recomputed.items():
        assert b["category_counts"].get(cat, 0) == n


# --------------------------------------------------------------------------- #
# FGFR2 specialization preserved
# --------------------------------------------------------------------------- #
def test_fgfr2_specialization_untouched():
    panels = FGFR2_EXAMPLE / "species_domain_architecture.json"
    if not panels.is_file():
        pytest.skip("FGFR2 validated example not present")
    data = json.loads(panels.read_text(encoding="utf-8"))
    has_cassette = any(s.get("panels", {}).get("IIIb") or s.get("panels", {}).get("IIIc")
                       for s in data.get("species", []))
    assert has_cassette, "FGFR2 IIIb/IIIc cassette panels must remain present"
