"""Regressions for the cross-pipeline release repair.

Four defects are covered here, each of which reached a reader as a wrong or
missing scientific statement rather than as an error:

1. a strand spelling decided the biological exon order, so the same gene was
   assembled 5'→3' from one source and 3'→5' from another;
2. a stop codon was counted as a residue, so a coding exon projected one amino
   acid past the end of its protein and the coordinate model failed validation
   silently, which is what left a completed run showing an obsolete gallery;
3. availability was maintained separately by six figure stages, so the record
   could describe a card set the figures had moved on from;
4. run readiness was derived from a persisted pre-cluster field instead of the
   run's artifacts, so a repaired stage kept a finished run out of results_ready.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from exondomaincompare.shared_gene_analysis.strand import (  # noqa: E402
    MINUS, PLUS, is_forward, is_reverse, normalize_strand, same_strand,
    strand_sign, strand_symbol,
)

FGFR2_RUN = "2026-07-29_1634_fgfr2_homo_sapiens_felis_catus"
BCL2L1_RUN = "2026-07-29_1646_bcl2l1_homo_sapiens_mus_musculus"


def run_dir(run_id: str) -> Path:
    return REPO_ROOT / "runs" / run_id


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path):
    import csv
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --------------------------------------------------------------------------- #
# Strand normalisation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("value", ["-", "-1", -1, -1.0, "reverse", "rev", "R",
                                   " -1 ", "MINUS", "antisense"])
def test_every_reverse_spelling_normalises_to_minus(value):
    assert normalize_strand(value) == MINUS
    assert is_reverse(value)
    assert not is_forward(value)


@pytest.mark.parametrize("value", ["+", "+1", 1, 1.0, "1", "forward", "fwd",
                                   "F", " +1 ", "PLUS", "sense"])
def test_every_forward_spelling_normalises_to_plus(value):
    assert normalize_strand(value) == PLUS
    assert is_forward(value)
    assert not is_reverse(value)


@pytest.mark.parametrize("value", [None, "", ".", "?", "na", "unknown", 0, "0"])
def test_unavailable_strand_is_unknown_not_forward(value):
    assert normalize_strand(value) is None


def test_ensembl_and_refseq_spellings_are_the_same_strand():
    """The defect in one line: Ensembl says -1 where RefSeq says -."""
    assert normalize_strand("-1") == normalize_strand("-")
    assert same_strand("-1", "-", -1, "reverse")
    assert same_strand("+1", "+", 1, "forward")
    assert not same_strand("-1", "+1")


def test_two_unknown_strands_do_not_agree():
    # Two absences of information are not an agreement about orientation.
    assert not same_strand("", "")


def test_symbol_and_sign_round_trip():
    assert strand_symbol("-1") == "-"
    assert strand_symbol(1) == "+"
    assert strand_symbol("unknown") == ""
    assert strand_sign("-1") == -1
    assert strand_sign("unknown", default=PLUS) == PLUS


def test_no_production_code_compares_a_raw_strand_value():
    """The normaliser is the only place that may decide what a strand means."""
    import re
    pattern = r"strand[a-z_]*\s*[=!]=\s*[\"']?-1?[\"']?"
    offenders = []
    for path in SCRIPTS.rglob("*.py"):
        if path.name == "strand.py":
            continue
        for line_number, line in enumerate(
                path.read_text(encoding="utf-8").splitlines(), start=1):
            if re.search(pattern, line) and not line.lstrip().startswith("#"):
                offenders.append(f"{path}:{line_number}:{line}")
    assert not offenders, "raw strand comparison outside the normaliser:\n" + \
        "\n".join(offenders)


# --------------------------------------------------------------------------- #
# Transcript ordering and protein projection
# --------------------------------------------------------------------------- #
def _cds_parts_minus_strand(strand):
    """Two CDS parts of a minus-strand transcript, in genomic ascending order."""
    return [
        {"seqid": "1", "start": "1000", "end": "1059", "strand": strand,
         "phase": "0", "exon_id": "exonB"},
        {"seqid": "1", "start": "2000", "end": "2029", "strand": strand,
         "phase": "0", "exon_id": "exonA"},
    ]


def _build(strand, protein_length_aa=30):
    from collect_fgfr2_models_dual_source_v3 import build_cds_features_from_parts
    return build_cds_features_from_parts(
        source_db="ensembl", species_input="test", species_canonical="test",
        transcript_id_internal="T1", transcript_id_source="T1",
        translation_id_source="P1",
        parts=_cds_parts_minus_strand(strand),
        coordinate_source="test",
        protein_length_aa=protein_length_aa,
    )


def _model(features):
    return [(f.cds_id_source, f.protein_start_aa, f.protein_end_aa)
            for f in features]


@pytest.mark.parametrize("spelling", ["-", -1, "reverse"])
def test_minus_strand_ordering_is_identical_across_source_spellings(spelling):
    """A source spelling must never change the biological exon order."""
    assert _model(_build(spelling)) == _model(_build("-1"))


def test_minus_strand_transcript_starts_at_the_highest_coordinate():
    """5'→3' on the minus strand runs from the highest genomic coordinate down."""
    features = _build("-1")
    assert features[0].cds_id_source == "exonA"
    assert features[0].protein_start_aa == "1"


def test_genomic_order_is_preserved_separately_from_transcript_order():
    """Source exon identity must survive independently of display order."""
    features = _build("-1")
    assert [f.transcript_order for f in features] == ["1", "2"]
    assert [f.genomic_order for f in features] == ["2", "1"]


def test_an_already_ordered_table_is_not_reversed_twice():
    """A transcript given in 5'→3' order with explicit ranks stays in it."""
    from collect_fgfr2_models_dual_source_v3 import build_cds_features_from_parts
    parts = [
        {"seqid": "1", "start": "2000", "end": "2029", "strand": "-1",
         "phase": "0", "exon_id": "exonA", "rank": "1"},
        {"seqid": "1", "start": "1000", "end": "1059", "strand": "-1",
         "phase": "0", "exon_id": "exonB", "rank": "2"},
    ]
    features = build_cds_features_from_parts(
        source_db="ensembl", species_input="t", species_canonical="t",
        transcript_id_internal="T1", transcript_id_source="T1",
        translation_id_source="P1", parts=parts, coordinate_source="test",
        protein_length_aa=30)
    assert [f.cds_id_source for f in features] == ["exonA", "exonB"]
    assert features[0].protein_start_aa == "1"


def test_protein_coordinates_are_monotonic_and_do_not_overlap():
    features = _build("-1")
    spans = [(int(f.protein_start_aa), int(f.protein_end_aa)) for f in features
             if f.protein_start_aa and f.protein_end_aa]
    assert spans == sorted(spans)
    for (_, prev_end), (next_start, _) in zip(spans, spans[1:]):
        assert next_start >= prev_end


def test_projection_never_exceeds_the_protein_length():
    """A stop codon is not a residue: BCL2L1 mouse is 235 aa, not 236."""
    features = _build("-1", protein_length_aa=25)
    ends = [int(f.protein_end_aa) for f in features if f.protein_end_aa]
    assert max(ends) <= 25


def test_ordering_fields_are_recorded_explicitly():
    features = _build("-1")
    first = features[0]
    assert first.normalized_strand == "-1"
    assert first.transcript_order == "1"
    assert first.coding_exon_order == "1"
    assert first.source_ordering_method


# --------------------------------------------------------------------------- #
# Comparative boundary groups (FGFR2 human/cat)
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def fgfr2_boundary_rows():
    path = (run_dir(FGFR2_RUN) / "results" / "generic_gene_analysis"
            / "comparative" / "comparable_boundary_groups.tsv")
    if not path.is_file():
        pytest.skip(f"{FGFR2_RUN} comparable boundary table not built")
    return read_tsv(path)


def test_two_species_comparison_yields_real_groups(fgfr2_boundary_rows):
    groups = {r["comparable_boundary_group_id"] for r in fgfr2_boundary_rows}
    assert len(groups) >= 2, "a two-species dataset must still produce groups"


def test_both_dataset_species_are_represented(fgfr2_boundary_rows):
    species = {r["species_id"] for r in fgfr2_boundary_rows}
    assert species == {"homo_sapiens", "felis_catus"}


def test_every_group_carries_both_species(fgfr2_boundary_rows):
    by_group = {}
    for row in fgfr2_boundary_rows:
        by_group.setdefault(row["comparable_boundary_group_id"], set()).add(
            row["species_id"])
    incomplete = {g: s for g, s in by_group.items() if len(s) < 2}
    assert not incomplete, f"single-species comparable groups: {incomplete}"


def test_groups_are_not_keyed_on_exon_rank_alone(fgfr2_boundary_rows):
    """Grouping evidence must be an alignment position, not an exon label."""
    methods = {r["mapping_method"] for r in fgfr2_boundary_rows}
    assert methods, "no mapping method recorded"
    assert not (methods & {"exon_rank", "exon_label", "genomic_order"})
    assert all(r["msa_column"] for r in fgfr2_boundary_rows)


def test_each_observation_carries_its_full_identity(fgfr2_boundary_rows):
    required = ("run_id", "species_id", "protein_id", "transcript_id",
                "boundary_id", "native_protein_position", "msa_column",
                "nearest_domain_instance_id", "nearest_edge", "signed_distance",
                "mapping_method", "mapping_confidence", "boundary_class")
    for row in fgfr2_boundary_rows:
        missing = [k for k in required if not row.get(k)]
        assert not missing, f"{row.get('boundary_id')} missing {missing}"


def test_domain_instance_ids_are_real_interpro_instances(fgfr2_boundary_rows):
    for row in fgfr2_boundary_rows:
        assert row["nearest_domain_instance_id"].startswith("IPR")


def test_signed_distances_agree_between_orthologous_boundaries(fgfr2_boundary_rows):
    """Human and cat FGFR2 place the same boundary at the same domain offset."""
    by_group = {}
    for row in fgfr2_boundary_rows:
        by_group.setdefault(row["comparable_boundary_group_id"], {})[
            row["species_id"]] = int(row["signed_distance"])
    spreads = [max(v.values()) - min(v.values()) for v in by_group.values()
               if len(v) == 2]
    assert spreads and max(spreads) <= 5, f"distance spreads: {spreads}"


# --------------------------------------------------------------------------- #
# Comparative synteny
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def fgfr2_synteny_svg():
    path = (run_dir(FGFR2_RUN) / "figures" / "comparative"
            / "cmp_comparative_synteny.svg")
    if not path.is_file():
        pytest.skip("comparative synteny figure not built")
    return path.read_text(encoding="utf-8")


def test_comparative_synteny_shows_both_species(fgfr2_synteny_svg):
    assert "H. sapiens" in fgfr2_synteny_svg or "Homo sapiens" in fgfr2_synteny_svg
    assert "F. catus" in fgfr2_synteny_svg or "Felis catus" in fgfr2_synteny_svg


def test_comparative_synteny_discloses_its_species_coverage(fgfr2_synteny_svg):
    """The caption has to name what the figure covers and what it does not."""
    assert "Species shown" in fgfr2_synteny_svg
    for word in ("complete", "partial", "no neighbourhood available"):
        assert word in fgfr2_synteny_svg


def test_comparative_synteny_names_a_target_locus(fgfr2_synteny_svg):
    assert "target locus" in fgfr2_synteny_svg


# --------------------------------------------------------------------------- #
# Figure registration
# --------------------------------------------------------------------------- #
def test_registration_rejects_a_card_from_another_run(tmp_path):
    from plotting.figure_registration import normalise_index
    run = tmp_path / "2026-01-01_0000_gene_species"
    (run / "figures").mkdir(parents=True)
    (run / "figures" / "own.png").write_bytes(b"x")
    doc = {"figures": [
        {"figure_id": "own", "status": "available",
         "png_url": f"runs/{run.name}/figures/own.png"},
        {"figure_id": "foreign", "status": "available",
         "png_url": "runs/2025-01-01_0000_other_run/figures/foreign.png"},
    ]}
    report = normalise_index(doc, run)
    assert [c["figure_id"] for c in doc["figures"]] == ["own"]
    assert report["n_rejected"] == 1
    assert "another run" in report["rejected"][0]["reasons"][0]


def test_a_card_whose_output_is_gone_is_reported_not_hidden(tmp_path):
    """An absent expected output is a gap the run must show, not one to delete."""
    from plotting.figure_registration import normalise_index
    run = tmp_path / "2026-01-01_0000_gene_species"
    (run / "figures").mkdir(parents=True)
    doc = {"figures": [{"figure_id": "ghost", "status": "available",
                        "png_url": f"runs/{run.name}/figures/ghost.png"}]}
    report = normalise_index(doc, run)
    assert doc["figures"][0]["status"] == "technically_missing"
    assert doc["available"] == []
    assert report["n_technically_missing"] == 1


def test_registration_rejects_a_superseded_card(tmp_path):
    from plotting.figure_registration import normalise_index
    run = tmp_path / "2026-01-01_0000_gene_species"
    (run / "figures").mkdir(parents=True)
    for name in ("old.png", "new.png"):
        (run / "figures" / name).write_bytes(b"x")
    doc = {"figures": [
        {"figure_id": "old", "status": "available",
         "png_url": f"runs/{run.name}/figures/old.png"},
        {"figure_id": "new", "status": "available", "supersedes": ["old"],
         "png_url": f"runs/{run.name}/figures/new.png"},
    ]}
    normalise_index(doc, run)
    assert [c["figure_id"] for c in doc["figures"]] == ["new"]


def test_registration_derives_availability_from_the_cards(tmp_path):
    """The availability record can no longer describe a superseded card set."""
    from plotting.figure_registration import normalise_index
    run = tmp_path / "2026-01-01_0000_gene_species"
    (run / "figures").mkdir(parents=True)
    for name in ("a.png", "b.png"):
        (run / "figures" / name).write_bytes(b"x")
    doc = {
        "figures": [
            {"figure_id": "a", "status": "available",
             "png_url": f"runs/{run.name}/figures/a.png"},
            {"figure_id": "b", "status": "pending_cluster",
             "png_url": f"runs/{run.name}/figures/b.png"},
        ],
        # A legacy entry no stage owns any more, and a card that no longer exists.
        "available": [{"id": "legacy_card"}, "removed_card"],
    }
    normalise_index(doc, run)
    assert doc["available"] == ["a"]
    assert doc["pending"] == ["b"]


def test_registration_records_the_scope_of_every_card(tmp_path):
    from plotting.figure_registration import normalise_index
    run = tmp_path / "2026-01-01_0000_gene_species"
    (run / "figures").mkdir(parents=True)
    for name in ("s.png", "c.png"):
        (run / "figures" / name).write_bytes(b"x")
    doc = {"figures": [
        {"figure_id": "s", "status": "available", "species_id": "homo_sapiens",
         "png_url": f"runs/{run.name}/figures/s.png"},
        {"figure_id": "c", "status": "available",
         "png_url": f"runs/{run.name}/figures/c.png"},
    ]}
    normalise_index(doc, run)
    scopes = {c["figure_id"]: c["scope"] for c in doc["figures"]}
    assert scopes == {"s": "species", "c": "comparative"}


def test_the_registration_stage_is_last_in_the_one_build_sequence():
    from plotting.figure_sequence import FIGURE_STAGES
    assert FIGURE_STAGES[-1][1] == "plotting.figure_registration"


# --------------------------------------------------------------------------- #
# Gallery inventories of the repaired runs
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def bcl2l1_index():
    path = run_dir(BCL2L1_RUN) / "website_indices" / "figures_index.json"
    if not path.is_file():
        pytest.skip(f"{BCL2L1_RUN} figure index not built")
    return read_json(path)


def test_bcl2l1_comparative_scope_is_populated(bcl2l1_index):
    comparative = [f for f in bcl2l1_index["figures"]
                   if f.get("scope") == "comparative"]
    assert len(comparative) >= 10


def test_bcl2l1_comparative_scope_carries_the_expected_analyses(bcl2l1_index):
    ids = {f["figure_id"] for f in bcl2l1_index["figures"]}
    for expected in ("cmp_primary_msa_overview", "cmp_msa_aligned_exon_architecture",
                     "cmp_domain_architecture_native", "cmp_domain_annotation_matrix",
                     "cmp_boundary_matrix", "cmp_paired_signed_distance",
                     "cmp_boundary_position_consistency", "cmp_comparative_synteny",
                     "cmp_isoform_diversity"):
        assert expected in ids, f"missing comparative figure {expected}"


def test_bcl2l1_species_scopes_reuse_the_modern_single_species_figures(bcl2l1_index):
    for species in ("homo_sapiens", "mus_musculus"):
        ids = {f["figure_id"] for f in bcl2l1_index["figures"]
               if f.get("species_id") == species}
        for suffix in ("primary_exon_projection", "integrated_domain_architecture",
                       "boundary_on_architecture", "signed_boundary_distances",
                       "boundary_class_summary", "full_isoform_alignment"):
            assert any(suffix in i for i in ids), f"{species} lacks {suffix}"
        assert any("local_gene_neighbourhood" in i for i in ids)


def test_bcl2l1_availability_record_matches_the_registered_cards(bcl2l1_index):
    ids = [f["figure_id"] for f in bcl2l1_index["figures"]]
    assert all(isinstance(x, str) for x in bcl2l1_index["available"])
    assert set(bcl2l1_index["available"]) <= set(ids)
    assert bcl2l1_index["run_id"] == BCL2L1_RUN


@pytest.fixture(scope="module")
def fgfr2_index():
    path = run_dir(FGFR2_RUN) / "website_indices" / "figure_index.json"
    if not path.is_file():
        pytest.skip(f"{FGFR2_RUN} figure index not built")
    return read_json(path)


def test_fgfr2_run_uses_the_modern_catalogue(fgfr2_index):
    categories = set(fgfr2_index["categories"])
    for expected in ("Comparative exon structure", "FGFR2 cassette evidence",
                     "Comparative sequence analysis", "Comparative domain architecture",
                     "FGFR2 IIIb/IIIc Boundary Consistency",
                     "Comparative exon–domain boundaries",
                     "Comparative genomic context"):
        assert expected in categories
    assert fgfr2_index["default_scope"] == "comparative"


def test_fgfr2_run_has_no_old_numbered_main_cards(fgfr2_index):
    """Figure 1–12 numbering belongs to supplements, never to a main card."""
    import re
    numbered = [f["figure_id"] for f in fgfr2_index["figures"]
                if f.get("kind") == "main"
                and re.search(r"figure_\d", f["figure_id"], re.I)]
    assert not numbered, numbered


def test_fgfr2_comparative_boundary_cards_are_available(fgfr2_index):
    cards = [f for f in fgfr2_index["figures"]
             if f["figure_id"].startswith("fgfr2_cmp_all_exon")]
    assert len(cards) == 3
    for card in cards:
        assert card["status"] == "available", card["figure_id"]
        assert card["scope"] == "comparative"
        formats = (card.get("modes") or [{}])[0].get("formats") or {}
        for fmt in ("png", "svg", "pdf", "tsv"):
            assert formats.get(fmt), f"{card['figure_id']} has no {fmt}"


def test_fgfr2_comparative_cards_reference_only_this_run(fgfr2_index):
    for card in fgfr2_index["figures"]:
        for mode in card.get("modes") or []:
            for path in (mode.get("formats") or {}).values():
                if "runs/" in path:
                    assert FGFR2_RUN in path, f"{card['figure_id']}: {path}"


def test_fgfr2_species_scopes_exist_for_every_dataset_species(fgfr2_index):
    scopes = {f.get("species_id") for f in fgfr2_index["figures"]
              if f.get("scope") == "species"}
    assert scopes == {"homo_sapiens", "felis_catus"}


# --------------------------------------------------------------------------- #
# Comparative availability wording
# --------------------------------------------------------------------------- #
def test_a_comparative_card_never_claims_a_protein_lacks_the_analysis():
    """`Not available for this protein` is a species-scope statement only."""
    source = (REPO_ROOT / "webapp" / "frontend" / "src" / "pages"
              / "FigureGallery.jsx").read_text(encoding="utf-8")
    assert "unavailableText(f)" in source
    marker = source.index("function unavailableText")
    body = source[marker:marker + 900]
    assert 'scope === "comparative"' in body
    assert "Not available for this dataset" in body


# --------------------------------------------------------------------------- #
# Cluster-output freshness
# --------------------------------------------------------------------------- #
def _fake_run(tmp_path, fasta_seq, scored_seq):
    import hashlib
    run = tmp_path / "2026-01-01_0000_gene_species"
    ip_in = run / "results" / "14_interproscan" / "primary" / "input"
    ip_out = run / "results" / "14_interproscan" / "primary" / "output"
    tm_out = (run / "results" / "15_exon_domain_boundary_post_interpro"
              / "pytmhmm_primary" / "output")
    for d in (ip_in, ip_out, tm_out):
        d.mkdir(parents=True, exist_ok=True)
    (ip_in / "final_pre_interpro_proteins_primary.faa").write_text(
        f">P1\n{fasta_seq}\n")
    (ip_out / "input.fasta.json").write_text(json.dumps({"results": [{
        "md5": hashlib.md5(scored_seq.upper().encode()).hexdigest(),
        "sequence": scored_seq, "xref": [{"id": "P1"}]}]}))
    (tm_out / "pytmhmm_summary_all.tsv").write_text("protein_id\tn_tm\nP1\t0\n")
    return run


def test_matching_cluster_outputs_are_reusable(tmp_path):
    from exondomaincompare.shared_gene_analysis.cluster_output_freshness import evaluate
    report = evaluate(_fake_run(tmp_path, "MKV", "MKV"))
    assert report["status"] == "fresh"
    assert report["usable"] is True


def test_cluster_outputs_for_a_superseded_sequence_are_stale(tmp_path):
    """A repaired coordinate model must not reuse annotations by filename."""
    from exondomaincompare.shared_gene_analysis.cluster_output_freshness import evaluate
    report = evaluate(_fake_run(tmp_path, "MKVA", "MKV"))
    assert report["status"] == "stale"
    assert report["interproscan"]["mismatched"] == ["P1"]


def test_stale_cluster_outputs_block_results_ready(tmp_path):
    from exondomaincompare.framework.species_completion import (
        STATE_STALE, aggregate_run_status,
    )
    completion = {"homo_sapiens": {
        "species_id": "homo_sapiens", "cluster_outputs": STATE_STALE,
        "complete": False, "blocking_analyses": ["cluster_outputs"]}}
    state, reasons = aggregate_run_status(completion, cluster_complete=True)
    assert state != "results_ready"
    assert reasons


def test_not_applicable_does_not_block_results_ready():
    from exondomaincompare.framework.species_completion import (
        REQUIRED_ANALYSES, STATE_AVAILABLE, STATE_NOT_APPLICABLE,
        aggregate_run_status,
    )
    record = {k: STATE_AVAILABLE for k in REQUIRED_ANALYSES}
    record["boundary_analysis"] = STATE_NOT_APPLICABLE
    record["species_id"] = "gallus_gallus"
    record["blocking_analyses"] = []
    record["complete"] = True
    state, reasons = aggregate_run_status({"gallus_gallus": record})
    assert state == "results_ready"
    assert reasons == []


# --------------------------------------------------------------------------- #
# Run status of the repaired and reference runs
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("run_id", [FGFR2_RUN, BCL2L1_RUN])
def test_repaired_runs_are_consistently_results_ready(run_id):
    path = run_dir(run_id) / "status.json"
    if not path.is_file():
        pytest.skip(f"{run_id} not present")
    status = read_json(path)
    assert status["status"] == "results_ready"
    assert status.get("run_status") == "results_ready"
    assert status.get("cluster_output_status") == "fresh"
    assert status.get("explorable") is True
    for stale in ("failed_step", "failed_reason"):
        assert not status.get(stale), f"{run_id} kept {stale}"
    for stage in ("pre_interpro_status", "post_interpro_status",
                  "website_indices_status"):
        assert status.get(stage) == "complete", f"{run_id}: {stage}"


@pytest.mark.parametrize("run_id", [
    "2026-07-29_1306_mc1r_gallus_gallus",
    "2026-07-29_1347_hba_panthera_leo",
    "2026-07-29_1526_akt1_mus_musculus",
    "2026-07-26_2157_fgfr1_gallus_mus_core_pilot",
    "2026-07-23_1100_fgfr1_gallus_core_pilot",
])
def test_accepted_reference_runs_keep_their_ready_status(run_id):
    """The new gates must not demote a run that was already correct."""
    directory = run_dir(run_id)
    if not (directory / "status.json").is_file():
        pytest.skip(f"{run_id} not present")
    from exondomaincompare.shared_gene_analysis.finalize_run_status import evaluate_run
    report = evaluate_run(directory)
    assert report["status"] == "results_ready", report["reason"]
    assert report["cluster_outputs"] == "fresh"


def test_readiness_is_not_derived_from_the_persisted_precluster_field():
    """The regression that kept a repaired FGFR2 run out of results_ready."""
    source = (SCRIPTS / "run_post_interpro_for_run.py").read_text(encoding="utf-8")
    assert '"complete" if prev.get("pre_interpro_status")' not in source
    assert "_finalize_run_status(rp.run_dir)" in source


# --------------------------------------------------------------------------- #
# Freeze protection
# --------------------------------------------------------------------------- #
def test_the_validated_freeze_is_bytewise_unchanged():
    import subprocess
    freeze = "results/final_30_until_interpro_prepare"
    if not (REPO_ROOT / freeze).is_dir():
        pytest.skip("validated freeze not present")
    out = subprocess.run(["git", "status", "--porcelain", "--", freeze],
                         cwd=REPO_ROOT, capture_output=True, text=True)
    assert out.stdout.strip() == "", f"freeze modified:\n{out.stdout}"
