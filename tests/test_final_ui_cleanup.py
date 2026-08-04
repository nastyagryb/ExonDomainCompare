from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FRONTEND = ROOT / "webapp" / "frontend" / "src"
VIEWERS = FRONTEND / "pages" / "viewers"


def source(relative: str) -> str:
    return (FRONTEND / relative).read_text(encoding="utf-8")


def test_evidence_stack_has_no_provenance_markers_or_source_link():
    overview = source("pages/Overview.jsx")
    css = source("App.css")

    assert 'className={`cell st-${d.class}`}' in overview
    assert "tone-corrected::after" not in css
    assert "tone-rescued::after" not in css
    assert "tone-offset::after" not in css
    assert 'label="Source table"' not in overview


def test_requested_explanatory_ui_is_removed_from_shared_views():
    checks = {
        "pages/Overview.jsx": "Cassette-end boundaries sit near",
        "pages/DataDownloads.jsx": "Run-level artefacts containing",
        "pages/viewers/CoordinateTrack.jsx": "Domain layer pending InterProScan",
        "pages/viewers/CassetteExplorer.jsx": "Gold markers indicate",
        "pages/viewers/DomainArchitecture.jsx": "Annotation flags",
        "pages/viewers/BoundaryDetailTrack.jsx": "Display note",
        "pages/viewers/MsaExplorer.jsx": "Cross-species alignment of one primary protein",
    }
    for relative, phrase in checks.items():
        assert phrase not in source(relative)

    msa = source("pages/viewers/MsaExplorer.jsx")
    assert "Column annotations are generic" not in msa


def test_story_view_is_not_exposed_by_the_frontend():
    assert not (VIEWERS / "SpeciesStory.jsx").exists()
    assert "SpeciesStory" not in source("pages/Overview.jsx")
    assert "SpeciesStory" not in source("pages/GeneExplorer.jsx")
    assert '["story", "Story"]' not in source("pages/GeneExplorer.jsx")


def test_page_navigation_is_restored_from_browser_history():
    app = source("App.jsx")
    assert 'searchParams.set("page", target)' in app
    assert 'window.history[method]' in app
    assert 'window.addEventListener("popstate", restorePage)' in app


def test_gallery_downloads_keep_the_application_loaded():
    gallery = source("pages/FigureGallery.jsx")

    assert gallery.count('target="_blank" rel="noreferrer"') >= 7
    for label in ("SVG", "PDF", "PNG", "TSV"):
        matching_lines = [line for line in gallery.splitlines() if f">{label}</a>" in line]
        assert matching_lines
        assert all('target="_blank"' in line for line in matching_lines)


def test_public_docs_explain_windows_setup_and_local_run_storage():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    windows = (ROOT / "docs" / "WINDOWS.md").read_text(encoding="utf-8")

    assert "docs/WINDOWS.md" in readme
    assert "~/.local/share/ExonDomainCompare/runs" in readme
    assert "~/Library/Application Support/ExonDomainCompare/runs" in readme
    assert "WSL2" in windows
    assert ".venv/bin/edc doctor" in windows
