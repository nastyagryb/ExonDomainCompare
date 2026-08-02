#!/usr/bin/env python3
"""Build the two read-only datasets distributed with ExonDomainCompare."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import zipfile
from pathlib import Path
from typing import Any, Iterable


FGFR2_SOURCE_REL = Path("results/final_30_until_interpro_prepare")
FGFR2_DERIVED_REL = Path("results/derived/example/website_indices")
BCL2L1_RUN_ID = "2026-07-29_1646_bcl2l1_homo_sapiens_mus_musculus"
BCL2L1_SOURCE_REL = Path("runs") / BCL2L1_RUN_ID
TEXT_SUFFIXES = {
    ".csv", ".faa", ".fasta", ".fa", ".gff", ".gff3", ".json", ".md",
    ".svg", ".tsv", ".txt", ".yaml", ".yml",
}
PUBLIC_FILE_SUFFIXES = {
    ".csv", ".faa", ".fasta", ".fa", ".json", ".md", ".pdf", ".png",
    ".svg", ".tsv", ".txt", ".yaml", ".yml", ".zip",
}
FORBIDDEN_TEXT = (
    "/Users/", "/home/", "/dss/", "/gpfs/",
)
IGNORED_PARTS = {".DS_Store", "__pycache__"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_file(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def copy_tree(
    source: Path, target: Path, *, excluded: tuple[Path, ...] = ()
) -> None:
    if not source.is_dir():
        raise FileNotFoundError(source)
    for path in sorted(source.rglob("*")):
        relative = path.relative_to(source)
        if any(relative == item or item in relative.parents for item in excluded):
            continue
        if any(part in IGNORED_PARTS for part in relative.parts):
            continue
        if path.is_symlink():
            raise RuntimeError(f"Release datasets must not contain symlinks: {path}")
        if path.is_file() and path.suffix != ".pyc":
            copy_file(path, target / relative)


def json_strings(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for child in value.values():
            yield from json_strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from json_strings(child)
    elif isinstance(value, str):
        yield value


def external_reference(value: str) -> str:
    assembly = re.search(r"(GCF_\d+\.\d+)", value)
    name = Path(value).name or "source_file"
    if assembly:
        return f"external:NCBI/{assembly.group(1)}/{name}"
    return f"external_input/{name}"


def sanitize_json(value: Any, source_repo: Path, key: str = "") -> Any:
    if isinstance(value, dict):
        return {name: sanitize_json(child, source_repo, name) for name, child in value.items()}
    if isinstance(value, list):
        return [sanitize_json(child, source_repo, key) for child in value]
    if not isinstance(value, str):
        return value
    normalized = value.replace("\\", "/")
    if "_ncbi_datasets_cache" in normalized:
        return external_reference(normalized)
    source_text = source_repo.as_posix().rstrip("/") + "/"
    if source_text in normalized:
        suffix = normalized.split(source_text, 1)[1]
        if "_ncbi_datasets_cache" in suffix:
            return external_reference(suffix)
        return suffix
    if normalized.startswith(("/Users/", "/home/", "/dss/", "/gpfs/")):
        return external_reference(normalized)
    return value


def sanitize_json_files(root: Path, source_repo: Path) -> None:
    for path in sorted(root.rglob("*.json")):
        data = json.loads(path.read_text(encoding="utf-8"))
        projected = sanitize_json(data, source_repo)
        path.write_text(
            json.dumps(projected, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )


def sanitize_text_files(root: Path, source_repo: Path) -> None:
    absolute = re.compile(r"/(?:Users|home|dss|gpfs)/[^\t\r\n <>'\"]+")
    source_prefix = source_repo.as_posix().rstrip("/") + "/"

    def replace(match: re.Match[str]) -> str:
        value = match.group(0)
        if value.startswith(source_prefix) and "_ncbi_datasets_cache" not in value:
            return value[len(source_prefix):]
        return external_reference(value)

    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES \
                or path.suffix.lower() == ".json":
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        projected = absolute.sub(replace, text)
        if projected != text:
            path.write_text(projected, encoding="utf-8")


def referenced_fgfr2_files(source_repo: Path) -> set[Path]:
    roots = (
        source_repo / FGFR2_SOURCE_REL / "13_final_pre_interpro_closure" / "website_indices",
        source_repo / FGFR2_DERIVED_REL,
    )
    found: set[Path] = set()
    prefix = FGFR2_SOURCE_REL.as_posix() + "/"
    for root in roots:
        for path in sorted(root.rglob("*.json")):
            data = json.loads(path.read_text(encoding="utf-8"))
            for value in json_strings(data):
                if not value.startswith(prefix) or "_ncbi_datasets_cache" in value:
                    continue
                candidate = source_repo / value
                if candidate.is_file() and candidate.suffix.lower() != ".zip" \
                        and "archive" not in candidate.relative_to(
                            source_repo / FGFR2_SOURCE_REL).parts:
                    found.add(candidate)
    return found


def referenced_fgfr2_derived_files(source_repo: Path) -> set[Path]:
    root = source_repo / FGFR2_DERIVED_REL
    prefix = "results/derived/example/"
    found: set[Path] = set()
    for path in sorted(root.rglob("*.json")):
        data = json.loads(path.read_text(encoding="utf-8"))
        for value in json_strings(data):
            if not value.startswith(prefix):
                continue
            candidate = source_repo / value
            if candidate.is_file():
                found.add(candidate)
    return found


def project_fgfr2(source_repo: Path, datasets_root: Path) -> Path:
    source = source_repo / FGFR2_SOURCE_REL
    target = datasets_root / "fgfr2_30_species"
    closure_rel = Path("13_final_pre_interpro_closure")
    copy_tree(
        source / closure_rel,
        target / closure_rel,
        excluded=(Path("archive"), Path("final_pre_interpro_run_log.txt")),
    )
    copy_tree(source_repo / FGFR2_DERIVED_REL, target / "derived" / "website_indices")
    for path in sorted(referenced_fgfr2_files(source_repo)):
        relative = path.relative_to(source)
        destination = target / relative
        if not destination.exists():
            copy_file(path, destination)
    derived_source = source_repo / "results" / "derived" / "example"
    for path in sorted(referenced_fgfr2_derived_files(source_repo)):
        relative = path.relative_to(derived_source)
        destination = target / "derived" / relative
        if not destination.exists():
            copy_file(path, destination)
    download_index = target / closure_rel / "website_indices" / "download_index.json"
    downloads = json.loads(download_index.read_text(encoding="utf-8"))
    downloads = [item for item in downloads if item.get("format") != "zip"]
    download_index.write_text(
        json.dumps(downloads, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    sanitize_json_files(target, source_repo)
    sanitize_text_files(target, source_repo)
    (target / "DATASET.md").write_text(
        "# FGFR2 IIIb/IIIc — 30 species\n\n"
        "This is the validated, read-only thesis dataset distributed with "
        "ExonDomainCompare. It is a compact website projection of the accepted "
        "FGFR2 freeze; the raw NCBI caches and the full working tree are not included.\n\n"
        "The projection does not change validated biology, coordinates, labels or "
        "accepted figures. Source paths that identify local machines or omitted raw "
        "caches are represented as external provenance references.\n",
        encoding="utf-8",
    )
    metadata = {
        "schema_version": "1.0",
        "dataset_id": "example",
        "release_id": "fgfr2_30_species",
        "gene_symbol": "FGFR2",
        "species_count": 30,
        "pipeline_type": "validated_event_pipeline",
        "support_level": "validated_event_analysis",
        "event": "FGFR2 IIIb/IIIc mutually exclusive cassette",
        "read_only": True,
        "source_record": FGFR2_SOURCE_REL.as_posix(),
        "distribution": "compact website projection",
        "omitted": ["raw NCBI caches", "mutable logs", "temporary files", "full working tree"],
    }
    (target / "dataset.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return target


BCL2L1_RESULT_DIRS = (
    "01_species_registry",
    "07_msa",
    "13_final_pre_interpro_closure",
    "14_interproscan/primary/output",
    "15_domain_architecture",
    "15_exon_domain_boundary_post_interpro/pytmhmm_primary/output",
    "16_final_analyses",
    "core_gene_analysis",
    "generic_gene_analysis",
)


def repair_bcl2l1_figure_sources(target: Path) -> None:
    for relative in (Path("website_indices/figures_index.json"),
                     Path("website_indices/generic/figures_index.json")):
        path = target / relative
        data = json.loads(path.read_text(encoding="utf-8"))
        for card in data.get("figures", []):
            values = [card.get("source_table", ""), *(card.get("source_files") or [])]
            if not any("isoform_alignment.fasta" in value for value in values):
                continue
            species_id = str(card.get("species_id") or "")
            source = f"results/generic_gene_analysis/msa/isoform_msa__{species_id}.aln.faa"
            card["source_table"] = source
            card["source_files"] = [source]
        path.write_text(
            json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    for relative in (Path("website_indices/available_views.json"),
                     Path("website_indices/generic/available_views.json")):
        path = target / relative
        data = json.loads(path.read_text(encoding="utf-8"))
        views = data.get("available_views") or data.get("views")
        if isinstance(views, dict):
            views["figure_gallery"] = True
            if "figures" in views:
                views["figures"] = True
        path.write_text(
            json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def public_bcl2l1_config(source: Path, source_repo: Path) -> dict[str, Any]:
    data = json.loads((source / "run_config.json").read_text(encoding="utf-8"))
    data = sanitize_json(data, source_repo)
    data.update({
        "gene_config": "run:gene_config.yaml",
        "species_list_path": "run:species_list.txt",
        "run_dir": "run:.",
        "results_dir": "run:results",
        "website_indices_dir": "run:website_indices",
        "cluster_input_fasta": "run:results/core_gene_analysis/proteins_primary.faa",
        "primary_fasta_path": "run:results/core_gene_analysis/proteins_primary.faa",
        "read_only": True,
        "release_projection": {
            "status": "bundled_read_only_example",
            "scientific_semantics": "exploratory_not_validated",
        },
    })
    return data


def public_bcl2l1_status(source: Path) -> dict[str, Any]:
    raw = json.loads((source / "status.json").read_text(encoding="utf-8"))
    keep = {
        "run_id", "status", "current_step", "run_mode", "experimental",
        "support_level", "has_event", "event_status", "species_count",
        "pre_interpro_status", "primary_fasta_status", "primary_fasta_count",
        "review_fasta_status", "review_fasta_count", "cluster_analysis_status",
        "cluster_fetch_status", "post_interpro_status", "website_indices_status",
        "next_action", "human_reference", "gene_symbol", "analysis_id",
        "pipeline_type", "event_layer_type", "event_analysis_enabled",
        "shared_pipeline", "stage_status", "cluster_status", "website_indices",
        "species_status", "run_status", "explorable", "readiness_reason",
        "cluster_output_status", "status_source",
    }
    projected = {key: raw[key] for key in keep if key in raw}
    projected.update({
        "read_only": True,
        "distribution_status": "bundled_read_only_example",
        "scientific_semantics": "exploratory_not_validated",
    })
    return projected


def project_bcl2l1(source_repo: Path, datasets_root: Path) -> Path:
    source = source_repo / BCL2L1_SOURCE_REL
    target = datasets_root / "runs" / BCL2L1_RUN_ID
    target.mkdir(parents=True)
    for name in ("gene_config.yaml", "species_list.txt"):
        copy_file(source / name, target / name)
    copy_tree(source / "website_indices", target / "website_indices")
    for name in BCL2L1_RESULT_DIRS:
        copy_tree(source / "results" / name, target / "results" / name)
    config = public_bcl2l1_config(source, source_repo)
    status = public_bcl2l1_status(source)
    (target / "run_config.json").write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (target / "status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    public = target / "website_indices" / "public"
    public.mkdir(parents=True, exist_ok=True)
    (public / "run_config.json").write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    repair_bcl2l1_figure_sources(target)
    sanitize_json_files(target, source_repo)
    sanitize_text_files(target, source_repo)
    (target / "DATASET.md").write_text(
        "# BCL2L1 — Homo sapiens and Mus musculus\n\n"
        "This bundled read-only dataset demonstrates the generic exploratory "
        "isoform workflow. Protein-difference Candidates are exploratory evidence, "
        "not validated splicing events. The dataset must not be presented as validated "
        "FGFR2-like event biology.\n\n"
        "Logs, SSH/Slurm state, user names, cluster paths, raw NCBI caches and private "
        "machine paths are intentionally omitted.\n",
        encoding="utf-8",
    )
    metadata = {
        "schema_version": "1.0",
        "dataset_id": f"run:{BCL2L1_RUN_ID}",
        "release_id": "bcl2l1_human_mouse",
        "run_id": BCL2L1_RUN_ID,
        "gene_symbol": "BCL2L1",
        "species": ["Homo sapiens", "Mus musculus"],
        "species_count": 2,
        "pipeline_type": "shared_gene_pipeline",
        "support_level": "generic_gene_analysis_experimental",
        "scientific_semantics": "exploratory_not_validated",
        "read_only": True,
        "source_record": BCL2L1_SOURCE_REL.as_posix(),
        "distribution": "sanitized website projection",
        "omitted": ["raw NCBI caches", "logs", "SSH and Slurm state", "private paths"],
    }
    (target / "dataset.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return target


def write_dataset_manifest(root: Path) -> None:
    excluded = {"file_manifest.tsv", "SHA256SUMS"}
    files = [
        path for path in sorted(root.rglob("*"))
        if path.is_file() and path.name not in excluded
    ]
    manifest = root / "file_manifest.tsv"
    lines = ["path\tbytes\tsha256"]
    for path in files:
        lines.append(f"{path.relative_to(root).as_posix()}\t{path.stat().st_size}\t{sha256(path)}")
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    checksum_files = files + [manifest]
    (root / "SHA256SUMS").write_text(
        "".join(f"{sha256(path)}  {path.relative_to(root).as_posix()}\n"
                for path in checksum_files),
        encoding="utf-8",
    )


def validate_text(root: Path) -> None:
    violations: list[str] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if path.stat().st_size >= 100 * 1024 * 1024:
            violations.append(f"oversized:{path.relative_to(root)}")
        if path.suffix.lower() in TEXT_SUFFIXES:
            text = path.read_text(encoding="utf-8", errors="replace")
            for token in FORBIDDEN_TEXT:
                if token in text:
                    violations.append(f"forbidden:{token}:{path.relative_to(root)}")
        if path.suffix.lower() == ".zip":
            with zipfile.ZipFile(path) as archive:
                for name in archive.namelist():
                    if any(token in name for token in FORBIDDEN_TEXT):
                        violations.append(f"forbidden-zip-name:{name}")
                    if Path(name).suffix.lower() in TEXT_SUFFIXES:
                        content = archive.read(name).decode("utf-8", errors="replace")
                        for token in FORBIDDEN_TEXT:
                            if token in content:
                                violations.append(f"forbidden-zip:{token}:{name}")
    if violations:
        raise RuntimeError("Release dataset validation failed:\n" + "\n".join(violations))


def validate_references(fgfr2: Path, bcl2l1: Path) -> None:
    missing: list[str] = []
    legacy_prefix = FGFR2_SOURCE_REL.as_posix() + "/"
    derived_prefix = "results/derived/example/"
    for path in sorted(fgfr2.rglob("*.json")):
        for value in json_strings(json.loads(path.read_text(encoding="utf-8"))):
            if value.startswith(legacy_prefix) and Path(value).suffix.lower() in PUBLIC_FILE_SUFFIXES:
                target = fgfr2 / Path(value).relative_to(FGFR2_SOURCE_REL)
                if not target.is_file():
                    missing.append(f"FGFR2:{value}")
            if value.startswith(derived_prefix) and Path(value).suffix.lower() in PUBLIC_FILE_SUFFIXES:
                target = fgfr2 / "derived" / Path(value).relative_to("results/derived/example")
                if not target.is_file():
                    missing.append(f"FGFR2:{value}")
    for path in sorted(bcl2l1.rglob("*.json")):
        for value in json_strings(json.loads(path.read_text(encoding="utf-8"))):
            candidate = value.split("?path=", 1)[1].split("&", 1)[0] if "?path=" in value else value
            run_prefix = f"runs/{BCL2L1_RUN_ID}/"
            if candidate.startswith(run_prefix):
                candidate = candidate[len(run_prefix):]
            if candidate.startswith("results/") and Path(candidate).suffix.lower() in PUBLIC_FILE_SUFFIXES:
                if not (bcl2l1 / candidate).is_file():
                    missing.append(f"BCL2L1:{candidate}")
    if missing:
        raise RuntimeError("Referenced release files are missing:\n" + "\n".join(sorted(set(missing))))


def build(source_repo: Path, release_repo: Path) -> None:
    source_repo = source_repo.resolve()
    release_repo = release_repo.resolve()
    if source_repo == release_repo:
        raise RuntimeError("Source and release repositories must be different.")
    if not (source_repo / FGFR2_SOURCE_REL).is_dir():
        raise FileNotFoundError("Validated FGFR2 source dataset is unavailable.")
    if not (source_repo / BCL2L1_SOURCE_REL).is_dir():
        raise FileNotFoundError("BCL2L1 source run is unavailable.")
    if not (release_repo / ".git").is_dir():
        raise RuntimeError("Release target must be the separate Git repository.")
    datasets = release_repo / "datasets"
    if datasets.exists():
        raise FileExistsError("Release datasets already exist; refusing to overwrite them.")
    datasets.mkdir()
    fgfr2 = project_fgfr2(source_repo, datasets)
    bcl2l1 = project_bcl2l1(source_repo, datasets)
    registry = {
        "schema_version": "1.0",
        "default_dataset": "example",
        "datasets": [
            {"id": "example", "path": "fgfr2_30_species", "read_only": True},
            {"id": f"run:{BCL2L1_RUN_ID}", "path": f"runs/{BCL2L1_RUN_ID}",
             "read_only": True},
        ],
    }
    (datasets / "registry.json").write_text(
        json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    (datasets / "README.md").write_text(
        "# Bundled datasets\n\n"
        "These compact, checksummed, read-only projections are displayed immediately "
        "after a clean installation. `fgfr2_30_species` is the validated thesis dataset; "
        "the BCL2L1 human/mouse run is an exploratory generic demonstration. Mutable new "
        "runs are stored separately in the user's configured application data directory.\n",
        encoding="utf-8",
    )
    validate_text(datasets)
    validate_references(fgfr2, bcl2l1)
    write_dataset_manifest(fgfr2)
    write_dataset_manifest(bcl2l1)
    root_files = [path for path in sorted(datasets.rglob("*"))
                  if path.is_file() and path != datasets / "SHA256SUMS"]
    (datasets / "SHA256SUMS").write_text(
        "".join(f"{sha256(path)}  {path.relative_to(datasets).as_posix()}\n"
                for path in root_files), encoding="utf-8")
    validate_text(datasets)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-repo", type=Path, required=True)
    parser.add_argument("--release-repo", type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.source_repo, args.release_repo)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Release datasets written to {(args.release_repo / 'datasets').resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
