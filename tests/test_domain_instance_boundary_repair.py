"""Regression tests for the repeated-domain-instance bug.

FGFR1 (Gallus gallus, NP_990841.2) carries **three** Ig-like domains that all
share the InterPro accession ``IPR007110`` but are distinct feature instances at
aa 33–118, 145–244 and 253–355. Resolving a domain by accession alone collapsed
them onto one another, so the nearest-domain coordinates persisted on a boundary
belonged to a different instance than the one the distance was measured against,
and the Boundary "Advanced analysis" domain filter could return zero rows for two
of the three instances while the third absorbed all Ig-like boundaries.

Everything asserted here is checked against the real reference runs — the FGFR1
Gallus core pilot and the TP53 Danio run — never against fixtures. The frontend
assertions execute the real ``common.js`` filter module through node, so the
option list and the filter resolution are tested as shipped, not as prose.
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SGA = ROOT / "scripts" / "shared_gene_analysis"
if str(SGA) not in sys.path:
    sys.path.insert(0, str(SGA))

from exondomaincompare.shared_gene_analysis import boundary_classification as bc  # noqa: E402
from exondomaincompare.shared_gene_analysis import protein_coordinate_model as pcm  # noqa: E402
from exondomaincompare.shared_gene_analysis import validate_protein_coordinate_model as vc  # noqa: E402

FGFR1_RUN = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"
TP53_RUN = ROOT / "runs" / "2026-07-21_1436_custom_run"
FRONTEND = ROOT / "webapp" / "frontend"

NEAR_EDGE_THRESHOLD_AA = 5

IG1 = "IPR007110:33-118"
IG2 = "IPR007110:145-244"
IG3 = "IPR007110:253-355"
KINASE = "IPR001245:476-750"

# The real FGFR1 Ig-like region: boundary label -> (position, instance id, edge
# type, edge position, signed distance, canonical class).
FGFR1_IG_BOUNDARIES = [
    ("E1 → E2", 31, IG1, "start", 33, -2, bc.NEAR),
    ("E2 → E3", 119, IG1, "end", 118, 1, bc.NEAR),
    ("E3 → E4", 147, IG2, "start", 145, 2, bc.NEAR),
    ("E4 → E5", 205, IG2, "end", 244, -39, bc.INSIDE),
    ("E5 → E6", 247, IG2, "end", 244, 3, bc.NEAR),
    ("E6 → E7", 310, IG3, "end", 355, -45, bc.INSIDE),
    ("E7 → E8", 359, IG3, "end", 355, 4, bc.NEAR),
]

pytestmark = [
    pytest.mark.skipif(not FGFR1_RUN.is_dir(), reason="FGFR1 Gallus pilot run not present"),
]


# --------------------------------------------------------------------------- #
# helpers — real data only
# --------------------------------------------------------------------------- #
def build(run: Path) -> dict:
    return pcm.build_models_for_run(run)["models"][0]


def served(run: Path) -> dict:
    """The coordinate model actually served to the frontend for this run."""
    path = run / "website_indices" / "generic" / "protein_coordinate_model.json"
    if not path.is_file():
        pytest.skip(f"served coordinate model missing for {run.name}")
    return json.loads(path.read_text(encoding="utf-8"))["models"][0]


@pytest.fixture(scope="module")
def fgfr1() -> dict:
    return build(FGFR1_RUN)


@pytest.fixture(scope="module")
def tp53() -> dict:
    if not TP53_RUN.is_dir():
        pytest.skip("TP53 Danio run not present")
    return build(TP53_RUN)


def by_label(model: dict) -> dict:
    return {b["label"]: b for b in model["exon_boundaries"]}


# --------------------------------------------------------------------------- #
# Task A — domain-instance identity in the coordinate model
# --------------------------------------------------------------------------- #
def test_every_domain_instance_has_a_unique_instance_id(fgfr1):
    domains = fgfr1["representative_domains"]
    ids = [d["domain_instance_id"] for d in domains]
    assert len(ids) == 4
    assert len(set(ids)) == len(ids), f"instance ids collapse: {ids}"
    assert ids == [IG1, IG2, IG3, KINASE]


def test_instance_id_is_derived_from_accession_and_coordinates(fgfr1):
    for d in fgfr1["representative_domains"]:
        assert d["domain_instance_id"] == f"{d['interpro_accession']}:{d['start']}-{d['end']}"


def test_instance_numbering_follows_sorted_start_coordinates(fgfr1):
    ig = [d for d in fgfr1["representative_domains"] if d["interpro_accession"] == "IPR007110"]
    assert [(d["instance_number"], d["start"], d["end"]) for d in ig] == [
        (1, 33, 118), (2, 145, 244), (3, 253, 355)]
    kinase = next(d for d in fgfr1["representative_domains"]
                  if d["interpro_accession"] == "IPR001245")
    assert kinase["instance_number"] == 1


def test_domain_instances_carry_the_full_identity_contract(fgfr1):
    required = ("domain_instance_id", "interpro_accession", "short_label", "full_label",
                "instance_number", "start", "end", "feature_type", "source",
                "member_signatures", "display_order")
    for d in fgfr1["representative_domains"]:
        for field in required:
            assert field in d, f"{d.get('id')} missing {field}"
        assert d["feature_type"] == "representative_domain"
        assert d["source"] == "InterProScan"
        assert isinstance(d["member_signatures"], list)
    assert [d["display_order"] for d in fgfr1["representative_domains"]] == [1, 2, 3, 4]


def test_repeated_entry_gets_numbered_labels_and_coordinate_suffixed_full_labels(fgfr1):
    by_id = {d["domain_instance_id"]: d for d in fgfr1["representative_domains"]}
    assert by_id[IG1]["short_label"] == "Ig-like domain 1"
    assert by_id[IG2]["short_label"] == "Ig-like domain 2"
    assert by_id[IG3]["short_label"] == "Ig-like domain 3"
    assert by_id[IG1]["full_label"] == "Ig-like domain 1 · aa 33–118"
    assert by_id[IG2]["full_label"] == "Ig-like domain 2 · aa 145–244"
    assert by_id[IG3]["full_label"] == "Ig-like domain 3 · aa 253–355"
    # A single-instance entry is NOT numbered, and keeps its real InterPro identity.
    assert by_id[KINASE]["short_label"] == "Ser-Thr/Tyr kinase domain"
    assert by_id[KINASE]["interpro_accession"] == "IPR001245"
    assert by_id[KINASE]["full_label"].endswith("· aa 476–750")


def test_member_signatures_are_resolved_per_instance_not_per_accession(fgfr1):
    """Each Ig-like instance owns its own real PS50835 hit at its own coordinates."""
    by_id = {d["domain_instance_id"]: d for d in fgfr1["representative_domains"]}
    for iid, (start, end) in ((IG1, (33, 118)), (IG2, (145, 244)), (IG3, (253, 355))):
        sigs = by_id[iid]["member_signatures"]
        assert sigs, f"{iid} has no contributing member signature"
        ig_like = [s for s in sigs if s["signature_accession"] == "PS50835"]
        assert [(s["start"], s["end"]) for s in ig_like] == [(start, end)], \
            f"{iid} picked up a PS50835 hit of another instance: {ig_like}"


# --------------------------------------------------------------------------- #
# Task B — boundaries persist the instance actually used
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("label,pos,instance,edge,edge_pos,signed,cls", FGFR1_IG_BOUNDARIES)
def test_fgfr1_ig_boundary_matches_the_real_reference_values(
        fgfr1, label, pos, instance, edge, edge_pos, signed, cls):
    b = by_label(fgfr1)[label]
    assert b["protein_position"] == pos
    assert b["nearest_domain_instance_id"] == instance
    assert b["nearest_edge_type"] == edge
    assert b["nearest_edge_position"] == edge_pos
    assert b["signed_distance"] == signed
    assert b["absolute_distance"] == abs(signed)
    assert b["boundary_class"] == cls


@pytest.mark.parametrize("label,pos,instance,edge,edge_pos,signed,cls", FGFR1_IG_BOUNDARIES)
def test_fgfr1_ig_boundary_carries_the_coordinates_of_the_instance_used(
        fgfr1, label, pos, instance, edge, edge_pos, signed, cls):
    """The bug: the persisted span belonged to a different IPR007110 instance."""
    b = by_label(fgfr1)[label]
    dom = next(d for d in fgfr1["representative_domains"]
               if d["domain_instance_id"] == instance)
    assert (b["nearest_domain_start"], b["nearest_domain_end"]) == (dom["start"], dom["end"])
    assert b["nearest_domain_accession"] == dom["interpro_accession"]
    assert b["nearest_domain_instance_number"] == dom["instance_number"]
    assert b["nearest_domain_label"] == dom["label"]
    assert b["nearest_domain_full_label"] == dom["full_label"]
    assert b["nearest_domain_id"] == dom["id"]


def test_all_three_ig_instances_own_boundaries(fgfr1):
    """No instance may absorb the boundaries of its siblings."""
    per_instance: dict[str, list[str]] = {}
    for b in fgfr1["exon_boundaries"]:
        per_instance.setdefault(b["nearest_domain_instance_id"], []).append(b["label"])
    assert per_instance[IG1] == ["E1 → E2", "E2 → E3"]
    assert per_instance[IG2] == ["E3 → E4", "E4 → E5", "E5 → E6"]
    assert per_instance[IG3] == ["E6 → E7", "E7 → E8"]
    assert len(per_instance[KINASE]) == 9


def test_fgfr1_class_distribution_over_16_internal_boundaries(fgfr1):
    bounds = fgfr1["exon_boundaries"]
    assert len(fgfr1["exons"]) == 17
    assert len(bounds) == 16
    counts = {c: 0 for c in bc.CANONICAL_CLASSES}
    for b in bounds:
        counts[b["boundary_class"]] += 1
    assert counts[bc.NEAR] == 6
    assert counts[bc.INSIDE] == 8
    assert counts[bc.OUTSIDE] == 2
    assert counts[bc.EXACT] == 0 and counts[bc.UNAVAILABLE] == 0
    assert fgfr1["near_edge_threshold_aa"] == NEAR_EDGE_THRESHOLD_AA


def _assert_invariants(model: dict) -> None:
    by_instance = {d["domain_instance_id"]: d for d in model["representative_domains"]}
    thr = model["near_edge_threshold_aa"]
    for b in model["exon_boundaries"]:
        iid = b["nearest_domain_instance_id"]
        if iid is None:
            assert b["signed_distance"] is None
            continue
        dom = by_instance[iid]
        pos, signed = b["protein_position"], b["signed_distance"]
        assert (b["nearest_domain_start"], b["nearest_domain_end"]) == (dom["start"], dom["end"])
        if b["nearest_edge_type"] == "start":
            assert b["nearest_edge_position"] == dom["start"]
            assert signed == pos - dom["start"], f"{b['label']}: {signed} != {pos} - {dom['start']}"
        else:
            assert b["nearest_edge_type"] == "end"
            assert b["nearest_edge_position"] == dom["end"]
            assert signed == pos - dom["end"], f"{b['label']}: {signed} != {pos} - {dom['end']}"
        assert b["absolute_distance"] == abs(signed)
        if b["boundary_class"] == bc.INSIDE:
            assert dom["start"] <= pos <= dom["end"], \
                f"{b['label']} is inside_domain but aa {pos} is outside {iid}"
        if b["boundary_class"] == bc.NEAR:
            assert abs(signed) <= thr


def test_fgfr1_boundary_invariants_hold(fgfr1):
    _assert_invariants(fgfr1)


def test_tp53_boundary_invariants_hold(tp53):
    _assert_invariants(tp53)


def test_tp53_domain_instances_are_unique_and_numbered(tp53):
    ids = [d["domain_instance_id"] for d in tp53["representative_domains"]]
    assert ids == ["IPR011615:69-258", "IPR010991:294-334"]
    assert len(set(ids)) == len(ids)
    assert all(d["instance_number"] == 1 for d in tp53["representative_domains"])


def test_sign_convention_is_documented_and_applied(fgfr1):
    """Negative before a start, positive after an end, negative inside vs the end."""
    assert "signed_distance = boundary_position - nearest_edge_position" in (
        bc.__doc__ or "")
    b = by_label(fgfr1)
    assert b["E1 → E2"]["signed_distance"] < 0    # before the Ig-like 1 start
    assert b["E5 → E6"]["signed_distance"] > 0    # after the Ig-like 2 end
    assert b["E4 → E5"]["signed_distance"] < 0    # inside, measured to the end


@pytest.mark.parametrize("run", [FGFR1_RUN, TP53_RUN])
def test_validator_enforces_instance_identity_and_invariants(run):
    if not run.is_dir():
        pytest.skip(f"{run.name} not present")
    index = pcm.build_models_for_run(run)
    assert vc.validate_index(index, core_dir=run / "results" / "core_gene_analysis") == []
    for model in index["models"]:
        assert vc.domain_instance_errors(model) == []
        assert vc.boundary_instance_errors(model) == []


def test_validator_rejects_a_boundary_pointing_at_the_wrong_instance(fgfr1):
    """The guard must actually fire — this is the shape of the original bug."""
    broken = json.loads(json.dumps(fgfr1))
    b = next(x for x in broken["exon_boundaries"] if x["label"] == "E1 → E2")
    b["nearest_domain_instance_id"] = IG3
    b["nearest_domain_start"], b["nearest_domain_end"] = 253, 355
    b["nearest_edge_position"] = 253
    errors = vc.boundary_instance_errors(broken)
    assert errors, "collapsing an Ig-like instance must be reported"
    assert any("signed_distance" in e for e in errors)


def _fgfr1_core_rows() -> tuple[list[dict], list[dict], list[dict]]:
    core = FGFR1_RUN / "results" / "core_gene_analysis"
    read = lambda name: list(csv.DictReader(  # noqa: E731
        (core / name).open(newline="", encoding="utf-8"), delimiter="\t"))
    domains = pcm._features_from_domain_rows(
        read("domain_features.tsv"), "gallus_gallus", "NP_990841.2",
        "results/core_gene_analysis/domain_features.tsv")
    exons = [{"id": f"NP_990841.2:exon{i}", "label": f"E{i}", "start": 1, "end": 1,
              "tooltip": {"exon_number": i}} for i in range(1, 18)]
    return read("exon_domain_boundary_distances.tsv"), domains, exons


def _resolved(rows: list[dict], domains: list[dict], exons: list[dict]) -> list[tuple]:
    out = pcm._boundary_features(
        rows, exons, domains, "gallus_gallus", "NP_990841.2",
        "results/core_gene_analysis/exon_domain_boundary_distances.tsv",
        cluster_complete=True)
    return [(b["label"], b["nearest_domain_instance_id"], b["nearest_edge_type"],
             b["nearest_edge_position"], b["signed_distance"], b["boundary_class"])
            for b in out]


def test_explicit_instance_id_column_is_used_when_present():
    """New Core runs persist nearest_domain_instance_id directly; honour it."""
    rows, domains, exons = _fgfr1_core_rows()
    by_edge = {(d["interpro_accession"], "start", d["start"]): d for d in domains}
    by_edge.update({(d["interpro_accession"], "end", d["end"]): d for d in domains})
    enriched = []
    for r in rows:
        r = dict(r)
        edge = r["nearest_edge"]
        edge_pos = int(r["boundary_position_aa"]) - int(r["signed_distance_aa"])
        r["nearest_domain_instance_id"] = by_edge[
            (r["nearest_domain_accession"], edge, edge_pos)]["domain_instance_id"]
        enriched.append(r)
    assert _resolved(enriched, domains, exons)[:7] == [
        (label, inst, edge, edge_pos, signed, cls)
        for label, _pos, inst, edge, edge_pos, signed, cls in FGFR1_IG_BOUNDARIES]


def test_rows_without_any_geometry_are_re_resolved_instance_aware():
    """A legacy row carrying only an accession must never collapse onto one instance."""
    rows, domains, exons = _fgfr1_core_rows()
    stripped = [{k: v for k, v in r.items()
                 if k not in ("signed_distance_aa", "absolute_distance_aa", "distance_aa")}
                for r in rows]
    assert _resolved(stripped, domains, exons)[:7] == [
        (label, inst, edge, edge_pos, signed, cls)
        for label, _pos, inst, edge, edge_pos, signed, cls in FGFR1_IG_BOUNDARIES]


def test_served_coordinate_model_carries_the_repaired_instances():
    """The JSON the webapp reads must agree with the freshly built model."""
    for run in (FGFR1_RUN, TP53_RUN):
        if not run.is_dir():
            continue
        fresh, disk = build(run), served(run)
        for a, b in zip(fresh["exon_boundaries"], disk["exon_boundaries"]):
            assert a["label"] == b["label"]
            assert a["nearest_domain_instance_id"] == b["nearest_domain_instance_id"]
            assert (a["nearest_domain_start"], a["nearest_domain_end"]) \
                == (b["nearest_domain_start"], b["nearest_domain_end"])
            assert a["signed_distance"] == b["signed_distance"]
        assert vc.boundary_instance_errors(disk) == []
        assert vc.domain_instance_errors(disk) == []


# --------------------------------------------------------------------------- #
# Task D — the Boundary "Advanced analysis" filters (real frontend module)
# --------------------------------------------------------------------------- #
def _node(script: str) -> dict:
    """Run a snippet against the real frontend modules and return its JSON result."""
    if not (FRONTEND / "node_modules").is_dir():
        pytest.skip("frontend node_modules not installed")
    proc = subprocess.run([  # noqa: S603 - fixed argv, no shell
        "node", "--input-type=module", "-e", script,
    ], cwd=FRONTEND, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


def _model_literal(run: Path) -> str:
    model = served(run)
    return json.dumps({"representative_domains": model["representative_domains"],
                       "exon_boundaries": model["exon_boundaries"]})


@pytest.fixture(scope="module")
def fgfr1_js() -> dict:
    """Options + per-option filter results produced by the shipped common.js."""
    return _node(f"""
      const {{ domainFilterOptions, filterBoundaries }} =
        await import('./src/pages/viewers/common.js');
      const model = {_model_literal(FGFR1_RUN)};
      const options = domainFilterOptions(model.representative_domains);
      const results = {{}};
      for (const o of options) {{
        results[o.value] = filterBoundaries(model.exon_boundaries, {{
          domainFilter: o.value, domainOptions: options,
        }}).map((b) => b.label);
      }}
      console.log(JSON.stringify({{ options, results }}));
    """)


def test_domain_filter_option_list_comes_from_the_real_data(fgfr1_js):
    assert [o["label"] for o in fgfr1_js["options"]] == [
        "All domains",
        "All Ig-like domains (3)",
        "Ig-like domain 1 · aa 33–118",
        "Ig-like domain 2 · aa 145–244",
        "Ig-like domain 3 · aa 253–355",
        "Ser-Thr/Tyr kinase domain · aa 476–750",
    ]
    # A single-instance entry gets no "all instances" group option.
    assert [o["value"] for o in fgfr1_js["options"] if o["kind"] == "group"] \
        == ["grp:IPR007110"]


def test_every_domain_filter_option_resolves_to_instance_ids_only(fgfr1_js):
    for o in fgfr1_js["options"]:
        if o["value"] == "all":
            assert o["instanceIds"] == []
            continue
        assert o["instanceIds"], f"{o['value']} resolves to no instance"
        for iid in o["instanceIds"]:
            assert re.fullmatch(r"IPR\d{6}:\d+-\d+", iid), \
                f"{o['value']} resolves by accession alone ({iid})"


def test_each_ig_like_instance_returns_its_own_non_empty_row_set(fgfr1_js):
    res = fgfr1_js["results"]
    assert res[f"inst:{IG1}"] == ["E1 → E2", "E2 → E3"]
    assert res[f"inst:{IG2}"] == ["E3 → E4", "E4 → E5", "E5 → E6"]
    assert res[f"inst:{IG3}"] == ["E6 → E7", "E7 → E8"]
    for iid in (IG1, IG2, IG3):
        assert res[f"inst:{iid}"], f"instance {iid} wrongly returns zero rows"
    # No instance absorbs another instance's boundaries.
    assert not (set(res[f"inst:{IG1}"]) & set(res[f"inst:{IG2}"]))
    assert not (set(res[f"inst:{IG2}"]) & set(res[f"inst:{IG3}"]))


def test_all_ig_like_domains_option_returns_the_union_of_the_three_instances(fgfr1_js):
    res = fgfr1_js["results"]
    union = res[f"inst:{IG1}"] + res[f"inst:{IG2}"] + res[f"inst:{IG3}"]
    assert res["grp:IPR007110"] == union
    assert len(res["grp:IPR007110"]) == 7
    assert res["all"] == union + res[f"inst:{KINASE}"]
    assert len(res["all"]) == 16


def test_filters_never_resolve_by_accession_alone(fgfr1_js):
    """An accession-shaped filter value must select nothing, not every instance."""
    out = _node(f"""
      const {{ domainFilterOptions, filterBoundaries }} =
        await import('./src/pages/viewers/common.js');
      const model = {_model_literal(FGFR1_RUN)};
      const options = domainFilterOptions(model.representative_domains);
      const bare = filterBoundaries(model.exon_boundaries,
        {{ domainFilter: 'IPR007110', domainOptions: options }}).length;
      const legacyId = filterBoundaries(model.exon_boundaries,
        {{ domainFilter: 'NP_990841.2:IPR007110:253-355', domainOptions: options }}).length;
      console.log(JSON.stringify({{ bare, legacyId }}));
    """)
    assert out["bare"] == 0
    assert out["legacyId"] == 0


def test_all_linked_filters_run_through_the_one_central_rule(fgfr1_js):
    """mapping / boundary / |dist| bounds / candidate / class / sort all apply."""
    out = _node(f"""
      const {{ domainFilterOptions, filterBoundaries }} =
        await import('./src/pages/viewers/common.js');
      const model = {_model_literal(FGFR1_RUN)};
      const options = domainFilterOptions(model.representative_domains);
      const rows = model.exon_boundaries;
      const f = (o) => filterBoundaries(rows, {{ domainOptions: options, ...o }})
        .map((b) => b.label);
      console.log(JSON.stringify({{
        mapped: f({{ mappingFilter: 'mapped' }}).length,
        unmapped: f({{ mappingFilter: 'unmapped' }}).length,
        oneBoundary: f({{ exonFilter: 'NP_990841.2:cds4_end' }}),
        nearOnly: f({{ distMax: 5 }}),
        farOnly: f({{ distMin: 100 }}),
        classNear: f({{ classFilter: new Set(['near_domain_edge']) }}).length,
        candidate: f({{ candOnly: true, candidate: {{ start: 31, end: 118 }} }}),
        byDistance: f({{ sort: 'distance' }}).slice(0, 3),
        combined: f({{ domainFilter: 'inst:{IG2}', classFilter: new Set(['near_domain_edge']) }}),
      }}));
    """)
    assert out["mapped"] == 16 and out["unmapped"] == 0
    assert out["oneBoundary"] == ["E4 → E5"]
    assert out["nearOnly"] == ["E1 → E2", "E2 → E3", "E3 → E4", "E5 → E6",
                               "E7 → E8", "E9 → E10"]
    assert out["farOnly"] == ["E12 → E13"]
    assert out["classNear"] == 6
    # Candidate C1 spans aa 31–118, so only the aa-31 junction falls inside it.
    assert out["candidate"] == ["E1 → E2"]
    assert out["byDistance"] == ["E2 → E3", "E1 → E2", "E3 → E4"]
    assert out["combined"] == ["E3 → E4", "E5 → E6"]


def test_tp53_domain_filter_options_have_no_spurious_group_option():
    if not TP53_RUN.is_dir():
        pytest.skip("TP53 Danio run not present")
    out = _node(f"""
      const {{ domainFilterOptions, filterBoundaries }} =
        await import('./src/pages/viewers/common.js');
      const model = {_model_literal(TP53_RUN)};
      const options = domainFilterOptions(model.representative_domains);
      console.log(JSON.stringify({{
        values: options.map((o) => o.value),
        dnaBd: filterBoundaries(model.exon_boundaries,
          {{ domainFilter: 'inst:IPR011615:69-258', domainOptions: options }}).length,
      }}));
    """)
    assert out["values"] == ["all", "inst:IPR011615:69-258", "inst:IPR010991:294-334"]
    assert out["dnaBd"] == 6


# --------------------------------------------------------------------------- #
# Task D — the central filtered array really feeds every linked view
# --------------------------------------------------------------------------- #
def _explorer() -> str:
    return (FRONTEND / "src" / "pages" / "viewers" / "BoundaryExplorer.jsx").read_text(
        encoding="utf-8")


def test_boundary_explorer_filters_through_the_shared_instance_aware_rule():
    text = _explorer()
    assert "filterBoundaries" in text and "domainFilterOptions" in text
    assert "const filteredBoundaries = useMemo" in text
    assert "const filtered = filteredBoundaries;" in text
    # The old accession-only comparison must be gone.
    assert "b.nearest_domain_accession !== domainFilter" not in text
    assert "b.nearest_domain_id !== domainFilter" not in text


def test_every_linked_view_reads_the_central_filtered_array():
    text = _explorer()
    views = {
        "summary counts": r"const counts = useMemo[\s\S]{0,400}?filteredBoundaries",
        "header badge": r"\$\{filteredBoundaries\.length\} of \$\{boundaries\.length\} shown",
        "architecture figure": r"tracks\.boundaries && filtered\.map",
        "signed-distance plot": r"<SignedDistancePlot rows=\{filtered\}",
        "evidence table": r"<tbody>\{filtered\.map",
        "visible TSV export": r"const rows = filtered\.map",
    }
    for name, pattern in views.items():
        assert re.search(pattern, text), f"{name} does not read filteredBoundaries"


def test_domain_filter_options_are_derived_not_hard_coded():
    text = _explorer()
    assert "domainFilterOptions(speciesModel?.representative_domains" in text
    for literal in ("IPR007110", "IPR001245", "Ig-like domain 1", "aa 33–118"):
        assert literal not in text, f"hard-coded domain literal {literal!r}"


def test_all_advanced_filters_are_present_resettable_and_counted():
    text = _explorer()
    for state in ("domainFilter", "mappingFilter", "exonFilter", "candOnly",
                  "distMin", "distMax", "sort"):
        assert f"const [{state}" in text, f"missing filter state {state}"
    reset = text[text.index("const resetFilters"):]
    reset = reset[:reset.index("};")]
    for setter in ("setClassFilter", "setDomainFilter", "setMappingFilter", "setExonFilter",
                   "setCandOnly", "setDistMin", "setDistMax", "setSort"):
        assert setter in reset, f"Reset all filters does not reset {setter}"
    assert "activeFilterCount" in text and "Reset all filters" in text


def test_boundary_transition_filter_is_labelled_with_real_exon_transitions(fgfr1):
    labels = [b["label"] for b in fgfr1["exon_boundaries"]]
    assert labels[:3] == ["E1 → E2", "E2 → E3", "E3 → E4"]
    assert labels[-1] == "E16 → E17"
    text = _explorer()
    assert re.search(r"boundaries\.map\(\(b\) => \([\s\S]{0,200}?\{b\.label\}", text), \
        "the boundary filter must be labelled with the real E1 → E2 transitions"


def test_no_dead_controls_left_in_the_advanced_filter_block():
    block = _explorer()
    block = block[block.index('className="bnd-filters"'):]
    block = block[:block.index('className="table-scroll"')]
    assert "disabled" not in block or "disabled={!domains.length}" in block
    # Every rendered control is wired to a state setter.
    for setter in ("setDomainFilter", "setMappingFilter", "setExonFilter",
                   "setDistMin", "setDistMax", "setCandOnly", "setSort", "resetFilters"):
        assert setter in block, f"control for {setter} is missing from the filter block"
