"""Part 3 / Part 4 / Part 6 — Comparative Figure Gallery, package builder, linked selection."""
from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "2026-07-26_2157_fgfr1_gallus_mus_core_pilot"
MODEL = RUN / "website_indices" / "generic" / "protein_coordinate_model.json"
CMP = RUN / "website_indices" / "generic" / "comparative_dataset_index.json"
FIG_IDX = RUN / "website_indices" / "figures_index.json"

pytestmark = pytest.mark.skipif(
    not MODEL.is_file(), reason="FGFR1 Gallus+Mus reference run not present")


@pytest.fixture(scope="module")
def comparative_index(tmp_path_factory):
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    sys.path.insert(0, str(ROOT / "scripts" / "shared_gene_analysis"))
    from comparative_dataset import build_comparative_dataset
    return build_comparative_dataset(RUN)


def test_comparative_dataset_index_is_available_for_two_species(comparative_index):
    assert comparative_index["available"] is True
    assert comparative_index["n_species"] == 2
    assert comparative_index["msa"]["available"] is True
    assert comparative_index["msa"]["n_columns"] > 0
    assert len(comparative_index["msa_aligned_exons"]) > 0
    assert len(comparative_index["msa_aligned_domains"]) > 0
    assert len(comparative_index["comparable_domain_groups"]) >= 1
    # Domain matrix never uses the forbidden "absent" vocabulary.
    for row in comparative_index["domain_annotation_matrix"]:
        assert row["state"] != "absent"
        assert row["state"] in {
            "detected", "not detected", "pending", "unavailable", "uncertain mapping",
        }


def test_comparable_domains_are_not_grouped_by_accession_alone(comparative_index):
    # FGFR1 has repeated Ig-like domains; groups must carry instance numbers.
    ig = [g for g in comparative_index["comparable_domain_groups"]
          if g.get("interpro_accession") == "IPR007110"]
    if len(ig) < 2:
        pytest.skip("this run does not expose repeated Ig-like instances")
    nums = {g["instance_number"] for g in ig}
    assert len(nums) == len(ig)


def test_comparative_gallery_cards_are_registered():
    doc = json.loads(FIG_IDX.read_text())
    cmp_cards = [f for f in doc.get("figures") or [] if f.get("scope") == "comparative"]
    assert len(cmp_cards) >= 8
    cats = {f.get("category") for f in cmp_cards}
    assert "Comparative exon structure" in cats
    assert "Comparative exon–domain boundaries" in cats
    for f in cmp_cards:
        assert f.get("png_url") or f.get("png_path")
        assert f.get("svg_url") or f.get("svg_path")
        assert f.get("pdf_url") or f.get("pdf_path")
        assert f.get("scientific_question")
        assert f.get("interpretation")
        # One card per figure — formats are not separate cards.
        assert f.get("kind") != "format"


def test_species_specific_gallery_cards_remain_for_both_species():
    doc = json.loads(FIG_IDX.read_text())
    by_sp = {}
    for f in doc.get("figures") or []:
        sid = f.get("species_id")
        if sid:
            by_sp.setdefault(sid, []).append(f)
    assert "gallus_gallus" in by_sp
    assert "mus_musculus" in by_sp
    # Multi-species must not empty a species' own Gallery structure.
    assert len(by_sp["gallus_gallus"]) >= 3
    assert len(by_sp["mus_musculus"]) >= 3


def test_package_capabilities_are_availability_aware():
    import sys
    sys.path.insert(0, str(ROOT / "scripts" / "shared_gene_analysis"))
    from package_builder import capabilities, resolve_dependencies
    caps = capabilities(RUN, "comparative")
    assert caps["multi_species"] is True
    assert [s["id"] for s in caps["scopes"]][:2] == ["comparative", "all"]
    assert any(s["id"] == "gallus_gallus" for s in caps["scopes"])
    assert "recommended" in caps["presets"]
    resolved = resolve_dependencies(caps["presets"]["recommended"]["items"],
                                    "comparative")
    assert "primary_proteins_msa" in resolved
    assert "boundary_long_table" in resolved
    # Unavailable items carry an exact reason, never a silent omit flag alone.
    for item in caps["items"].values():
        if not item["available"]:
            assert item["reason"]


def test_recommended_package_builds_valid_zip_with_manifest():
    import sys
    sys.path.insert(0, str(ROOT / "scripts" / "shared_gene_analysis"))
    from package_builder import PACKAGES_ROOT, build_package
    job = build_package(RUN, {"preset": "recommended", "scope": "comparative"})
    assert job.status == "ready", job.error
    assert job.zip_path.startswith("package:")
    zip_path = PACKAGES_ROOT / job.zip_path[len("package:"):]
    assert zip_path.is_file()
    assert zip_path.stat().st_size > 1000
    with zipfile.ZipFile(zip_path) as zf:
        assert zf.testzip() is None
        names = zf.namelist()
        assert any(n.endswith("README.md") for n in names)
        assert any(n.endswith("manifest.json") for n in names)
        assert any(n.endswith("package_selection.json") for n in names)
        assert any("comparative_results.xlsx" in n for n in names)
        assert any("species_inventory.tsv" in n for n in names)
        # A comparative-only package carries no per-species tree.
        assert not any("/species/" in n for n in names)


def test_dependency_resolution_pulls_msa_for_aligned_exons():
    import sys
    sys.path.insert(0, str(ROOT / "scripts" / "shared_gene_analysis"))
    from package_builder import resolve_dependencies
    resolved = resolve_dependencies(["msa_aligned_exons"], "comparative")
    assert resolved.index("primary_proteins_msa") < resolved.index("msa_aligned_exons")


def test_fgfr2_freeze_untouched_by_package_builder():
    # The package builder must not rewrite the validated FGFR2 example indices.
    example = (ROOT / "results" / "final_30_until_interpro_prepare"
               / "13_final_pre_interpro_closure" / "website_indices")
    if not example.is_dir():
        pytest.skip("FGFR2 example missing")
    before = {p.name: p.stat().st_mtime for p in example.glob("*.json")}
    import sys
    sys.path.insert(0, str(ROOT / "scripts" / "shared_gene_analysis"))
    from package_builder import capabilities
    # Capabilities against the multi-species run only.
    capabilities(RUN, "comparative")
    after = {p.name: p.stat().st_mtime for p in example.glob("*.json")}
    assert before == after


def test_single_species_run_produces_no_comparative_gallery_cards():
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from plotting.generate_comparative_gallery_figures import generate
    # Prefer a known single-species FGFR1 run if present.
    candidates = sorted((ROOT / "runs").glob("*fgfr1*gallus*"))
    single = None
    for c in candidates:
        pcm = c / "website_indices" / "generic" / "protein_coordinate_model.json"
        if not pcm.is_file():
            continue
        models = json.loads(pcm.read_text()).get("models") or []
        if len(models) == 1:
            single = (c, pcm)
            break
    if not single:
        pytest.skip("no single-species FGFR1 run available")
    run_dir, pcm = single
    res = generate(run_dir, pcm)
    assert res.get("skipped") == "single_species" or res.get("cards", 0) == 0
