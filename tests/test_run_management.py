from __future__ import annotations

import json
import re
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exondomaincompare.framework import run_labels as rl  # noqa: E402

FRONTEND = ROOT / "webapp" / "frontend" / "src"


def src(rel: str) -> str:
    return (FRONTEND / rel).read_text(encoding="utf-8")


# --------------------------------------------------------------------------- #
# Part 11 — the New Run form asks for biology, not for internals
# --------------------------------------------------------------------------- #
def test_new_run_form_has_no_advanced_section_or_species_preset():
    text = src("pages/runworkflow/CreateRunPanel.jsx")
    assert "rs-advanced" not in text
    assert "Advanced" not in text
    assert "Species panel preset" not in text
    assert "Custom species (default)" not in text
    assert "full30" not in text and '"pilot"' not in text


def test_new_run_form_does_not_explain_the_workflow_router():
    text = src("pages/runworkflow/CreateRunPanel.jsx")
    rendered = text.split("return (", 1)[1]
    for phrase in ("validated IIIb/IIIc pipeline", "exploratory workflow",
                   "Validated workflow", "Generic exploratory workflow",
                   "run configuration is generated automatically"):
        assert phrase not in rendered, f"the form still explains internals: {phrase}"


def test_new_run_form_keeps_the_four_real_inputs():
    text = src("pages/runworkflow/CreateRunPanel.jsx")
    for field in ("Gene symbol", "Run name", "Species", "Upload .txt",
                  "Create and start run"):
        assert field in text
    # Species parsing behaviour that must survive the simplification.
    assert 'replace(/[;,]/g, "\\n")' in text        # commas and semicolons
    assert 'startsWith("#")' in text                # comment lines
    assert "new Set()" in text                      # duplicate removal


def test_an_empty_run_name_is_submitted_as_empty():
    text = src("pages/runworkflow/CreateRunPanel.jsx")
    assert "custom_run" not in text
    assert "run_name: runName.trim()," in text


# --------------------------------------------------------------------------- #
# Part 12 — run_id and run_name are different things
# --------------------------------------------------------------------------- #
def test_a_visible_run_name_keeps_spaces_capitals_and_punctuation():
    assert rl.clean_run_name("  FGFR1 chicken and mouse (validation)  ") == \
        "FGFR1 chicken and mouse (validation)"
    assert rl.clean_run_name("TP53 — mammals, v2") == "TP53 mammals, v2"
    assert rl.clean_run_name("a" * 400) == "a" * rl.MAX_RUN_NAME


def test_placeholder_names_count_as_no_name():
    for placeholder in ("custom_run", "Custom Run", "run", "full30_run", "  "):
        assert rl.clean_run_name(placeholder) == ""


def test_a_pipeline_generated_default_is_not_treated_as_a_user_name():
    assert rl.clean_run_name("fgfr1_gallus_mus_core_pilot", "FGFR1") == ""
    assert rl.clean_run_name("tp53_human_core_pilot", "TP53") == ""
    assert rl.clean_run_name("tpm1_human_mouse_twospecies", "TPM1") == ""
    # Anything a person plausibly typed keeps its exact wording.
    assert rl.clean_run_name("FGFR1 gallus mus", "FGFR1") == "FGFR1 gallus mus"
    assert rl.clean_run_name("Salmo salar", "FGFR2") == "Salmo salar"
    assert rl.clean_run_name("mammals_only", "FGFR1") == "mammals_only"


def test_the_run_id_slug_is_filesystem_safe_and_never_the_raw_name():
    slug = rl.run_id_slug("FGFR1 chicken & mouse!")
    assert re.fullmatch(r"[a-z0-9_]+", slug), slug
    assert " " not in slug and "&" not in slug
    assert rl.run_id_slug("", gene_symbol="FGFR1",
                          species=["gallus_gallus", "mus_musculus"]) == \
        "fgfr1_gallus_gallus_mus_musculus"
    assert rl.run_id_slug("", gene_symbol="TP53", species=["danio_rerio"]) == "tp53_danio_rerio"
    assert rl.run_id_slug("", gene_symbol="FGFR2",
                          species=[f"sp_{i}" for i in range(30)]) == "fgfr2_30species"


def test_duplicate_visible_run_names_are_allowed():
    name = "FGFR1 validation"
    assert rl.clean_run_name(name) == rl.clean_run_name(name)
    # The uniqueness rule lives in the directory allocator, not in the label.
    import create_new_run as cnr
    assert "while path.exists()" in \
        Path(cnr.__file__).read_text(encoding="utf-8")


def test_a_run_without_a_name_is_titled_from_gene_and_species():
    assert rl.display_label("", gene_symbol="FGFR1",
                            species=["gallus_gallus", "mus_musculus"]) == \
        "FGFR1 · Gallus gallus + Mus musculus"
    assert rl.display_label("custom_run", gene_symbol="TP53", species=["danio_rerio"]) == \
        "TP53 · Danio rerio"
    assert rl.display_label("", gene_symbol="FGFR2",
                            species=[f"sp_{i}" for i in range(30)]).startswith("FGFR2 · Sp 0 + 29")


def test_a_named_run_shows_its_name():
    assert rl.display_label("My FGFR1 run", gene_symbol="FGFR1",
                            species=["gallus_gallus"]) == "My FGFR1 run"


def test_the_species_binomial_capitalises_only_the_genus():
    assert rl.short_species("gallus_gallus") == "Gallus gallus"
    assert rl.short_species("canis_lupus_familiaris") == "Canis lupus familiaris"


def test_the_backend_derives_the_slug_and_stores_the_clean_name():
    main = (ROOT / "webapp" / "backend" / "main.py").read_text(encoding="utf-8")
    assert "run_labels.clean_run_name(req.run_name)" in main
    assert "run_labels.run_id_slug(run_name" in main
    assert "record = cnr.build_run_config(" in main
    assert "RunLayout(run_dir, RunLayoutVersion.CANONICAL_V2).initialize" in main


def test_the_core_runner_separates_the_label_from_the_directory_name():
    runner = (ROOT / "src" / "exondomaincompare" / "framework"
              / "run_core_gene_analysis.py").read_text(encoding="utf-8")
    assert "run_labels.clean_run_name(args.run_name)" in runner
    assert "run_labels.run_id_slug(run_name" in runner
    assert "cnr.generate_run_id(slug)" in runner


# --------------------------------------------------------------------------- #
# Part 13 — newest first
# --------------------------------------------------------------------------- #
def _run(run_id: str, created_at: str = "") -> dict:
    return {"run_id": run_id, "created_at": created_at}


def test_runs_are_ordered_by_created_at_descending():
    now = datetime(2026, 7, 29, 12, tzinfo=timezone.utc)
    runs = [
        _run("a", (now - timedelta(days=2)).isoformat()),
        _run("b", now.isoformat()),
        _run("c", (now - timedelta(hours=1)).isoformat()),
    ]
    assert [r["run_id"] for r in rl.sort_runs(runs)] == ["b", "c", "a"]


def test_a_legacy_run_falls_back_to_the_timestamp_in_its_run_id():
    runs = [
        _run("2026-07-21_1436_custom_run"),
        _run("2026-07-29_0900_fgfr1_gallus"),
        _run("2026-07-23_1100_fgfr1_gallus_core_pilot"),
    ]
    assert [r["run_id"] for r in rl.sort_runs(runs)] == [
        "2026-07-29_0900_fgfr1_gallus",
        "2026-07-23_1100_fgfr1_gallus_core_pilot",
        "2026-07-21_1436_custom_run",
    ]


def test_created_at_wins_over_the_run_id_timestamp():
    runs = [
        _run("2026-01-01_0000_old_id", "2026-07-29T10:00:00+00:00"),
        _run("2026-07-29_1200_new_id", "2026-02-01T10:00:00+00:00"),
    ]
    assert rl.sort_runs(runs)[0]["run_id"] == "2026-01-01_0000_old_id"


def test_the_order_is_stable_for_runs_with_no_time_information():
    runs = [_run("zulu"), _run("alpha"), _run("mike")]
    once = [r["run_id"] for r in rl.sort_runs(runs)]
    assert once == [r["run_id"] for r in rl.sort_runs(list(reversed(runs)))]


def test_the_backend_never_orders_runs_by_filesystem_time():
    main = (ROOT / "webapp" / "backend" / "main.py").read_text(encoding="utf-8")
    assert "run_labels.sort_runs(" in main
    listing = main.split("def local_runs(", 1)[1].split("\n@app", 1)[0]
    for forbidden in ("st_mtime", "getmtime", ".stat()"):
        assert forbidden not in listing


# --------------------------------------------------------------------------- #
# Part 14 / 15 — the card shows the run, not the machinery
# --------------------------------------------------------------------------- #
def test_technical_details_and_logs_are_gone_from_the_run_page():
    page = src("pages/RunWorkflowPage.jsx")
    assert "Technical details" not in page
    assert "rw-tech" not in page
    assert "Show technical log" not in page
    assert "rw-log-pre" not in page
    assert "PreInterproPanel" not in page
    assert not (FRONTEND / "pages/runworkflow/PreInterproPanel.jsx").exists()


def test_a_finished_run_shows_no_stage_checklist():
    page = src("pages/RunWorkflowPage.jsx")
    assert 'status !== "results_ready" && (' in page
    assert "<RunStatusStepper" in page  # still used while a run is working


def test_the_portable_roundtrip_command_outranks_legacy_status_text():
    page = src("pages/RunWorkflowPage.jsx")
    portable = page.index("commands?.cluster_roundtrip?.portable_command")
    persisted = page.index("model.cluster_command", portable)
    assert portable < persisted


def test_a_failed_run_stays_actionable_without_a_terminal_dump():
    page = src("pages/RunWorkflowPage.jsx")
    panel = page.split("function RunFailurePanel", 1)[1].split("\nfunction ", 1)[0]
    assert "Retry with the same input" in panel
    assert "Download diagnostics" in panel
    assert "api.runDiagnosticsUrl" in panel
    assert "<pre" not in panel


def test_the_card_shows_the_display_name_with_the_run_id_as_metadata():
    card = src("pages/runworkflow/RunCard.jsx")
    assert "run.display_name" in card
    assert "Run ID:" in card
    assert "run.completion_summary" in card
    assert "run.failure_summary" in card


def test_runs_are_grouped_and_each_group_keeps_the_newest_first():
    page = src("pages/RunWorkflowPage.jsx")
    assert '["active", "Active"]' in page
    assert '["attention", "Attention required"]' in page
    assert '["completed", "Completed"]' in page
    # Grouping filters the already-ordered list rather than re-sorting it.
    assert "runs.filter((r) => (r.group" in page


@pytest.mark.parametrize("status,group", [
    ("results_ready", "completed"),
    ("running", "active"),
    ("pre_interpro_running", "active"),
    ("cluster_required", "active"),
    ("failed", "attention"),
    ("core_model_collection_failed", "attention"),
    ("stopped", "attention"),
])
def test_status_grouping(status, group):
    assert rl.run_group(status) == group


def test_the_completion_summary_replaces_the_stage_checklist():
    text = rl.completion_summary(primary_fasta_count=2,
                                 available_views={"domain_architecture": True,
                                                  "boundary": True})
    assert text == ("2 primary proteins analysed. Domain architecture and "
                    "Boundary analysis are available.")
    assert rl.completion_summary(primary_fasta_count=1,
                                 available_views={}) == "1 primary protein analysed."


def test_the_failure_summary_names_the_species_and_the_stage():
    assert rl.describe_failure(failed_stage="domain_architecture",
                               failed_species="mus_musculus") == \
        "Mus musculus: domain architecture could not be generated."


def test_logs_are_preserved_and_downloadable():
    main = (ROOT / "webapp" / "backend" / "main.py").read_text(encoding="utf-8")
    assert "/api/local-runs/{run_id}/diagnostics" in main
    bundle = main.split("def local_run_diagnostics(", 1)[1].split("\n@app", 1)[0]
    assert 'run_dir / "logs"' in bundle and "zipfile.ZipFile" in bundle
    # The endpoints that expose logs for diagnostics are still there; removing
    # the accordion hid the logs from the interface, it did not delete them.
    assert "/api/local-runs/{run_id}/logs/core" in main
    assert "/api/local-runs/{run_id}/logs/preinterpro" in main
    assert "def append_local_run_log" in main


# --------------------------------------------------------------------------- #
# Existing run creation still works
# --------------------------------------------------------------------------- #
def test_existing_runs_keep_their_stored_names():
    runs_root = ROOT / "runs"
    if not runs_root.is_dir():
        pytest.skip("no runs directory")
    for cfg_path in sorted(runs_root.glob("*/run_config.json"))[:20]:
        cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
        label = rl.display_label(cfg.get("run_name"),
                                 gene_symbol=cfg.get("gene_symbol", ""),
                                 species=cfg.get("species_ids") or [],
                                 species_count=cfg.get("species_count", 0),
                                 run_id=cfg.get("run_id", cfg_path.parent.name))
        assert label and label != "custom_run"
        assert "custom_run" not in label
