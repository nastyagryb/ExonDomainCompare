#!/usr/bin/env python3
"""Create a new, self-contained FGFR2 IIIb/IIIc run folder under ``runs/``.

This is step 1 of the new run logic: SAFE local run-folder creation only. It

  * never touches the validated example freeze
    (``results/final_30_until_interpro_prepare/`` stays read-only),
  * creates ``runs/<run_id>/`` with the full step layout, config/status files,
    a normalized species list, and copy-paste command templates,
  * does NOT run any pipeline, InterProScan, pyTMHMM, SSH, or SLURM command.

Cluster jobs are always submitted explicitly from the local terminal using the
scripts referenced in ``00_README_NEXT_STEPS.md``; the webapp never collects
LRZ passwords or 2FA codes.

Examples
--------
    python scripts/create_new_run.py --preset full30 --run-name fgfr2_test
    python scripts/create_new_run.py --preset pilot  --run-name fgfr2_pilot
    python scripts/create_new_run.py --species-list path/to/species.txt --run-name custom_species
    python scripts/create_new_run.py --species "homo_sapiens,mus_musculus,gallus_gallus" --run-name custom3
"""
from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.framework import production_contract  # noqa: E402
from exondomaincompare.contracts import portable_path_reference, stamp_payload  # noqa: E402
from exondomaincompare.config import discover_repository_root, load_config  # noqa: E402
from exondomaincompare.runs.layout import RunLayout, RunLayoutVersion  # noqa: E402

RUNTIME_CONFIG = load_config(repository_root=discover_repository_root(__file__))
REPO = RUNTIME_CONFIG.repository_root
# Isolated test roots and user profiles redirect the registry through one resolver.
RUNS_ROOT = RUNTIME_CONFIG.runs_root
REFERENCE_DIR = REPO / "reference"
FULL30_LIST = REFERENCE_DIR / "Species_list_final_30.txt"
PILOT_GLOBS = [
    "Species_list_pilot.txt",
    "Species_list_pilot_*.txt",
    "pilot_species*.txt",
]

CASE_STUDY = "FGFR2_IIIb_IIIc"
GENE_SYMBOL = "FGFR2"

# Gene/event generalization layer (additive; does not change FGFR2 logic).
# New runs record which configured gene/event analysis they belong to. FGFR2
# IIIb/IIIc is the only active analysis; older runs default to it on read.
ANALYSIS_ID = "FGFR2_IIIb_IIIc"
EVENT_ID = "FGFR2_IIIb_IIIc_cassette"
EVENT_TYPE = "mutually_exclusive_cassette"
GENE_CONFIG_REL = "configs/genes/FGFR2_IIIb_IIIc.yaml"

# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _rel(path: Path) -> str:
    """Portable logical path; never serialize a personal absolute path."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO))
    except ValueError:
        try:
            return "run:" + str(resolved.relative_to(RUNS_ROOT))
        except ValueError:
            return portable_path_reference(resolved, repository_root=REPO)


def sanitize_run_name(name: str) -> str:
    name = (name or "").strip().lower().replace(" ", "_")
    name = re.sub(r"[^a-z0-9_-]+", "", name)   # keep letters, digits, _ and -
    name = re.sub(r"_{2,}", "_", name).strip("_-")
    return name or "run"


def sanitize_run_id(run_id: str) -> str:
    run_id = (run_id or "").strip().replace(" ", "_")
    run_id = re.sub(r"[^A-Za-z0-9_-]+", "", run_id)
    return run_id


def generate_run_id(run_name: str, when: Optional[datetime] = None) -> str:
    when = when or datetime.now()
    stamp = when.strftime("%Y-%m-%d_%H%M")
    return f"{stamp}_{sanitize_run_name(run_name)}"


def unique_run_dir(run_id: str) -> tuple[str, Path]:
    """Return a (run_id, path) pair that does not collide, appending _2, _3 …"""
    candidate = run_id
    path = RUNS_ROOT / candidate
    n = 2
    while path.exists():
        candidate = f"{run_id}_{n}"
        path = RUNS_ROOT / candidate
        n += 1
    return candidate, path


# Ensembl-style species identifier, e.g. gallus_gallus, canis_lupus_familiaris.
SPECIES_ID_RE = re.compile(r"^[a-z][a-z0-9]+_[a-z0-9_]+$")


def normalize_species_token(raw: str) -> str:
    """Normalize one species entry to an Ensembl-style identifier.

    'Gallus gallus' -> 'gallus_gallus'; 'homo_sapiens' -> 'homo_sapiens'.
    Lowercases, converts internal whitespace to underscores, collapses repeats
    and strips leading/trailing underscores. Does not validate.
    """
    s = (raw or "").strip().lower()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"_{2,}", "_", s).strip("_")
    return s


def human_reference_info(species: List[str]) -> dict:
    """Describe how curated human FGFR2 IIIb/IIIc is used for this run.

    Human is ALWAYS available as a curated reference/control layer taken from the
    validated example dataset (marker comparison, IIIb/IIIc orientation, calibration).
    It is only counted as an analysed species when the user explicitly selected
    homo_sapiens; otherwise the analysed species panel is left unchanged.
    """
    in_panel = any((s or "").strip().lower() == "homo_sapiens" for s in species)
    if in_panel:
        return {
            "enabled": True,
            "source": "run_panel_and_validated_example_dataset",
            "human_role": "analysed_species_plus_reference_control",
            "homo_sapiens_in_panel": True,
            "species_panel_unchanged": True,
            "note": ("homo_sapiens is part of the selected analysed panel and is also "
                     "linked to the curated human FGFR2 IIIb/IIIc reference."),
        }
    return {
        "enabled": True,
        "source": "validated_example_dataset",
        "human_role": "human_reference_control",
        "homo_sapiens_in_panel": False,
        "species_panel_unchanged": True,
        "note": ("Curated human FGFR2 IIIb/IIIc is reused only as a reference/control "
                 "layer and is not added as an analysed species."),
    }




def species_error_message(raw: str) -> str:
    """Clear, actionable error for an invalid species entry.

    Only offers "Did you mean X?" when the normalized suggestion is itself a
    valid identifier (e.g. legacy 'gallus gallus' -> 'gallus_gallus');
    otherwise explains the expected format.
    """
    suggestion = normalize_species_token(raw)
    if suggestion and SPECIES_ID_RE.match(suggestion):
        return f'Invalid species identifier: "{raw}". Did you mean "{suggestion}"?'
    return (f'Invalid species identifier: "{raw}". Expected a lowercase '
            'genus_species identifier, e.g. "homo_sapiens".')


def passthrough_species_lines(raw_lines: List[str]) -> List[str]:
    """Trim, drop empties/comments, de-duplicate — NO case/space transform.

    Used for curated preset files (e.g. the validated full30 reference list),
    which must be written exactly as maintained so preset runs stay unchanged.
    """
    seen = set()
    out: List[str] = []
    for line in raw_lines:
        s = (line or "").strip()
        if not s or s.startswith("#"):
            continue
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def normalize_species_lines(raw_lines: List[str]) -> List[str]:
    """Normalize, drop empties/comments, validate, and de-duplicate (order kept).

    Used for user-supplied custom input. Raises SystemExit with a clear,
    actionable message if any entry is not a valid lowercase-underscore
    identifier after normalization.
    """
    seen = set()
    out: List[str] = []
    for line in raw_lines:
        s = (line or "").strip()
        if not s or s.startswith("#"):
            continue
        token = normalize_species_token(s)
        if not token or not SPECIES_ID_RE.match(token):
            raise SystemExit(species_error_message(s))
        if token not in seen:
            seen.add(token)
            out.append(token)
    return out


def species_from_inline(text: str) -> List[str]:
    parts = re.split(r"[,\n;]", text or "")
    return normalize_species_lines(parts)


def species_from_file(path: Path, normalize: bool = True) -> List[str]:
    if not path.exists():
        raise SystemExit(f"ERROR: species list file not found: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    return normalize_species_lines(lines) if normalize else passthrough_species_lines(lines)


def find_pilot_list() -> Optional[Path]:
    for pattern in PILOT_GLOBS:
        matches = sorted(REFERENCE_DIR.glob(pattern))
        if matches:
            return matches[0]
    return None


def resolve_species(args: argparse.Namespace) -> tuple[List[str], str]:
    """Return (species_lines, source_label). Exactly one source is required."""
    sources = [bool(args.preset), bool(args.species_list), bool(args.species)]
    if sum(sources) != 1:
        raise SystemExit(
            "ERROR: provide exactly one species source: "
            "--preset, --species-list, or --species."
        )

    if args.preset:
        # Curated preset files are written through unchanged (no normalization),
        # so the validated full30 preset run stays byte-identical to the freeze.
        if args.preset == "full30":
            return (species_from_file(FULL30_LIST, normalize=False),
                    f"preset:full30 ({_rel(FULL30_LIST)})")
        if args.preset == "pilot":
            pilot = find_pilot_list()
            if pilot is None:
                raise SystemExit(
                    "ERROR: no pilot species file found. Looked for "
                    + ", ".join(str(REFERENCE_DIR / p) for p in PILOT_GLOBS)
                    + ".\nCreate one of these files or use --species-list / --species."
                )
            return species_from_file(pilot, normalize=False), f"preset:pilot ({_rel(pilot)})"
        raise SystemExit(f"ERROR: unknown preset '{args.preset}'.")

    # User-supplied custom sources ARE normalized + validated.
    if args.species_list:
        p = Path(args.species_list)
        return species_from_file(p, normalize=True), f"species_list:{_rel(p)}"

    return species_from_inline(args.species), "species:inline"


# --------------------------------------------------------------------------- #
# folder + file creation
# --------------------------------------------------------------------------- #




def build_run_config(run_id: str, run_name: str, run_dir: Path,
                     species: List[str], species_list_path: Path,
                     source_label: str) -> dict:
    cfg = {
        "run_id": run_id,
        "run_name": run_name,
        "case_study": CASE_STUDY,
        "gene_symbol": GENE_SYMBOL,
        # Gene/event generalization metadata (additive; FGFR2 remains canonical).
        "analysis_id": ANALYSIS_ID,
        "event_id": EVENT_ID,
        "event_type": EVENT_TYPE,
        "gene_config": "run:config/gene.json",
        "species_count": len(species),
        "species_list_path": "run:config/species.tsv",
        "run_dir": "run:.",
        "results_dir": "run:scientific",
        "website_indices_dir": "run:website/indices",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "created_by": f"local_profile:{RUNTIME_CONFIG.local_profile_name}",
        "configuration_profile": RUNTIME_CONFIG.public_identity(),
        "primary_fasta_expected_path": "run:inputs/cluster/primary.faa",
        "cluster_mode": "manual_local_scripts",
        "cluster_submission_policy": "commands_generated_for_local_terminal_execution",
        "notes": (
            "Run folder created by scripts/create_new_run.py (step 1: safe folder "
            "creation only). No pipeline, InterProScan, pyTMHMM, SSH, or SLURM "
            "command was executed. The validated example freeze at "
            "results/final_30_until_interpro_prepare/ is untouched."
        ),
    }
    # Stamp the production architecture at creation time, resolved from the gene symbol
    # alone. Every later stage reads it from here instead of re-deciding, so no run can
    # acquire a different Explorer or Gallery architecture than the one it was created
    # with, and no run needs a migration command to get the modern one.
    production_contract.stamp(cfg)
    # Record logical provenance; a user list's absolute location is only an
    # execution-time detail and must not become part of the run contract.
    portable_source_label = source_label
    if source_label.startswith("species_list:"):
        portable_source_label = f"species_list:external:{species_list_path.name}"
    if source_label.startswith("preset:"):
        cfg["source_preset"] = source_label.split(" ", 1)[0].split(":", 1)[1]
    else:
        cfg["source_species_list"] = portable_source_label
    cfg["source_label"] = portable_source_label
    cfg["human_reference"] = human_reference_info(species)
    return stamp_payload(
        cfg,
        payload_type="run_config",
        run_id=run_id,
        dataset_id=run_id,
        profile=RUNTIME_CONFIG.public_identity(),
        generator="scripts/create_new_run.py",
    )


def build_status(run_id: str, species: List[str], run_dir: Path) -> dict:
    compatibility_roundtrip = RUNTIME_CONFIG.command([
        "python", "scripts/interpro_cluster/run_cluster_roundtrip.py",
        "--run-id", run_id,
        "--local-profile", RUNTIME_CONFIG.local_profile_name,
        "--lrz-profile", RUNTIME_CONFIG.lrz_profile_name,
    ])
    roundtrip = RUNTIME_CONFIG.command([
        ".venv/bin/edc", "cluster", "roundtrip", "--run-id", run_id,
        "--local-profile", RUNTIME_CONFIG.local_profile_name,
        "--lrz-profile", RUNTIME_CONFIG.lrz_profile_name,
    ])
    return stamp_payload({
        "run_id": run_id,
        "status": "created",
        "current_step": "run_created",
        "next_action": "run_pre_interpro",
        "species_count": len(species),
        "pre_interpro_status": "not_started",
        "primary_fasta_status": "not_available",
        "primary_fasta_count": 0,
        "review_fasta_status": "not_available",
        "interproscan_status": "not_started",
        "pytmhmm_status": "not_started",
        "post_interpro_status": "not_started",
        "website_indices_status": "not_started",
        "human_reference": human_reference_info(species),
        "last_updated": datetime.now(timezone.utc).isoformat(),
        "cluster_jobs": {},
        "cluster_profile": RUNTIME_CONFIG.lrz_profile_name,
        "cluster_continuation_command": roundtrip,
        "cluster_continuation_compatibility_command": compatibility_roundtrip,
        "next_actions": [
            "Run the pre-InterPro pipeline for this run "
            "(python scripts/run_pre_interpro_for_run.py --run-id "
            f"{run_id}) [planned].",
            "After the primary FASTA exists, submit the cluster analysis using the "
            "local command "
            f"(python scripts/interpro_cluster/submit_cluster_analysis.py --run-id {run_id}).",
        ],
    }, payload_type="status", run_id=run_id, dataset_id=run_id,
       profile=RUNTIME_CONFIG.public_identity(),
       generator="scripts/create_new_run.py")


def script_status_label(rel_path: str) -> str:
    return "available" if (REPO / rel_path).exists() else "planned"




# --------------------------------------------------------------------------- #
# validation + summary
# --------------------------------------------------------------------------- #
def validate_run(run_dir: Path, species: List[str]) -> None:
    problems: List[str] = []
    if not (run_dir / "config" / "species.tsv").exists():
        problems.append("config/species.tsv missing")
    if len(species) <= 0:
        problems.append("species_count is 0")
    if not (run_dir / "run.json").exists():
        problems.append("run.json missing")
    if not (run_dir / "status.json").exists():
        problems.append("status.json missing")
    if not (run_dir / "config" / "gene.json").is_file():
        problems.append("config/gene.json missing")
    expected = {
        "run.json", "status.json", "config/gene.json", "config/species.tsv"}
    actual = {
        path.relative_to(run_dir).as_posix()
        for path in run_dir.rglob("*") if path.is_file()
    }
    if actual != expected:
        problems.append(
            "initial tree is not minimal: " + ", ".join(sorted(actual - expected)))
    if problems:
        raise SystemExit("ERROR: run validation failed:\n  - " + "\n  - ".join(problems))


def main(argv: Optional[List[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        description="Create a new FGFR2 IIIb/IIIc run folder (safe, local, no execution).")
    ap.add_argument("--preset", choices=["full30", "pilot"],
                    help="use a bundled species list preset")
    ap.add_argument("--species-list", help="path to a species list file (one per line)")
    ap.add_argument("--species", help="inline species, comma- or newline-separated")
    ap.add_argument("--run-name", default="run", help="human-readable run name")
    ap.add_argument("--run-id", help="explicit run id (skips auto-generation)")
    ap.add_argument("--force", action="store_true",
                    help="with an existing --run-id: fill missing files only (never deletes)")
    ap.add_argument("--config", help="explicit configuration TOML")
    ap.add_argument("--local-profile", help="named local profile")
    ap.add_argument("--lrz-profile", help="named LRZ profile")
    args = ap.parse_args(argv)

    global RUNTIME_CONFIG, REPO, RUNS_ROOT, REFERENCE_DIR, FULL30_LIST
    RUNTIME_CONFIG = load_config(
        config_path=args.config,
        repository_root=REPO,
        local_profile=args.local_profile,
        lrz_profile=args.lrz_profile,
    )
    REPO = RUNTIME_CONFIG.repository_root
    RUNS_ROOT = RUNTIME_CONFIG.runs_root
    REFERENCE_DIR = REPO / "reference"
    _FULL30_LIST = REFERENCE_DIR / "Species_list_final_30.txt"

    species, source_label = resolve_species(args)
    if not species:
        raise SystemExit("ERROR: the resolved species list is empty.")

    RUNS_ROOT.mkdir(parents=True, exist_ok=True)

    forced_existing = False
    if args.run_id:
        run_id = sanitize_run_id(args.run_id)
        if not run_id:
            raise SystemExit("ERROR: --run-id sanitized to an empty string.")
        run_dir = RUNS_ROOT / run_id
        if run_dir.exists():
            if not args.force:
                raise SystemExit(
                    f"ERROR: run id '{run_id}' already exists at {_rel(run_dir)}.\n"
                    "Use --force to fill missing directories/files (nothing is deleted)."
                )
            forced_existing = True
            print(f"WARNING: --force on existing run '{run_id}'. Existing files are "
                  "kept; only missing directories/files are created.")
    else:
        base_id = generate_run_id(args.run_name)
        run_id, run_dir = unique_run_dir(base_id)

    if forced_existing:
        raise SystemExit(
            "ERROR: --force is not supported for canonical runs; use retry/resume.")
    species_path = run_dir / "config" / "species.tsv"
    cfg = build_run_config(
        run_id, args.run_name, run_dir, species, species_path, source_label)
    gene = {
        "gene_symbol": GENE_SYMBOL,
        "analysis_id": ANALYSIS_ID,
        "event_id": EVENT_ID,
        "event_type": EVENT_TYPE,
        "source": f"repo:{GENE_CONFIG_REL}",
    }
    RunLayout(run_dir, RunLayoutVersion.CANONICAL_V2).initialize(
        run_record=cfg,
        status=build_status(run_id, species, run_dir),
        gene=gene,
        species=species,
    )

    # 6) validation
    validate_run(run_dir, species)

    # 7) summary
    pre_cmd = f"python scripts/run_pre_interpro_for_run.py --run-id {run_id}"
    print("\n" + "=" * 64)
    print("  NEW RUN CREATED")
    print("=" * 64)
    print(f"  run_id        : {run_id}")
    print(f"  run_dir       : {_rel(run_dir)}")
    print(f"  species_count : {len(species)}")
    print(f"  source        : {source_label}")
    print(f"  next command  : {pre_cmd}"
          + ("" if script_status_label('scripts/run_pre_interpro_for_run.py') == 'available'
             else "   [planned]"))
    print("  layout        : canonical-2.0 (lazy)")
    print("=" * 64 + "\n")


if __name__ == "__main__":
    main()
