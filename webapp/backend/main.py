from __future__ import annotations

import importlib.util
import json
import logging
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import uuid
import zipfile
from concurrent.futures import Future
from datetime import datetime, timezone
from pathlib import Path
from threading import Lock
from typing import Any, Dict, List, Mapping, Optional

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

try:
    from .canonical_dataset import build_canonical_dataset_model
except ImportError:  # uvicorn main:app from webapp/backend
    from canonical_dataset import build_canonical_dataset_model

# --- project wiring -------------------------------------------------------- #
APP_VERSION = "1.0.0-phase1"
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))
from exondomaincompare.contracts import stamp_payload, verify_payload_contract  # noqa: E402
from exondomaincompare.runs.legacy import LegacyRunAdapter, LegacyRunError  # noqa: E402
from exondomaincompare.runs.registry import (  # noqa: E402
    RegistryError, RunCollisionError, discover_runs, hide_discovered_run,
    resolve_run_record, unregister,
)
from exondomaincompare.config import (CONFIG_ENV, DATA_ENV, LOCAL_PROFILE_ENV,
                                      LRZ_PROFILE_ENV, RUNS_ENV, load_config)
from exondomaincompare.runs.layout import RunLayout, RunLayoutVersion  # noqa: E402
from exondomaincompare.shared_gene_analysis.analysis_availability import (  # noqa: E402
    index_version,
)

RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)
PROJECT_ROOT = RUNTIME_CONFIG.repository_root
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
BUNDLED_FGFR2_ROOT = PROJECT_ROOT / "datasets" / "fgfr2_30_species"
LEGACY_FGFR2_ROOT = PROJECT_ROOT / "results" / "final_30_until_interpro_prepare"
RESULTS_ROOT = (
    BUNDLED_FGFR2_ROOT
    if (BUNDLED_FGFR2_ROOT / "13_final_pre_interpro_closure" / "website_indices").is_dir()
    else LEGACY_FGFR2_ROOT
)
EXAMPLE_RUN_DIR = RESULTS_ROOT / "13_final_pre_interpro_closure"
WEB_RUNS_ROOT = PROJECT_ROOT / "results" / "web_runs"
WEB_STATE_DIR = RUNTIME_CONFIG.paths.config / "state"
LEGACY_WEB_STATE_DIR = PROJECT_ROOT / "results" / "web_state"
SHARED_NCBI_CACHE = RESULTS_ROOT / "02_models" / "_ncbi_datasets_cache"
REFERENCE_SPECIES_LIST = PROJECT_ROOT / "reference" / "Species_list_final_30.txt"
RUNNER_SCRIPT = PROJECT_ROOT / "run_fgfr2_pipeline_current_final_pre_interpro.sh"
LOCAL_PYTHON_RUNTIME = RUNTIME_CONFIG.local_python()
PYTHON = LOCAL_PYTHON_RUNTIME.selected


def _local_python_command(script: Path, *args: str) -> List[str]:
    return [PYTHON, str(script), *map(str, args)]


def _local_python_module_command(module: str, *args: str) -> List[str]:
    return [PYTHON, "-m", module, *map(str, args)]


CURRENT_RUN_PTR = WEB_STATE_DIR / "current_run.json"
LEGACY_CURRENT_RUN_PTR = LEGACY_WEB_STATE_DIR / "current_run.json"


_bwi_spec = importlib.util.spec_from_file_location(
    "build_website_indices", SCRIPTS_DIR / "build_website_indices.py"
)
if _bwi_spec is None or _bwi_spec.loader is None:
    raise ImportError(
        f"Could not load build_website_indices from {SCRIPTS_DIR}"
    )
bwi = importlib.util.module_from_spec(_bwi_spec)
sys.modules["build_website_indices"] = bwi
_bwi_spec.loader.exec_module(bwi)

_cnr_spec = importlib.util.spec_from_file_location(
    "create_new_run", SCRIPTS_DIR / "create_new_run.py"
)
if _cnr_spec is None or _cnr_spec.loader is None:
    raise ImportError(f"Could not load create_new_run from {SCRIPTS_DIR}")
cnr = importlib.util.module_from_spec(_cnr_spec)
sys.modules["create_new_run"] = cnr
_cnr_spec.loader.exec_module(cnr)

from exondomaincompare.framework import run_labels  # noqa: E402

try:
    from build_species_registry_improved import (  # type: ignore
        lookup_known_species as _lookup_known_species,
        slug_species as _slug_species,
    )
except Exception:
    _lookup_known_species = None
    _slug_species = None

try:
    from exondomaincompare.framework import gene_config as _gene_config  # type: ignore
except Exception:
    _gene_config = None

try:
    from exondomaincompare.framework.core_run_milestones import (  # type: ignore
        evaluate_core_run as _evaluate_core_run,
        is_core_only_run as _is_core_only_run,
    )
except Exception:  # pragma: no cover
    _evaluate_core_run = None
    _is_core_only_run = None

try:
    from exondomaincompare.shared_gene_analysis.run_availability import (  # type: ignore
        readiness as _readiness,
        view_states as _view_states,
    )
except Exception:  # pragma: no cover
    _readiness = None
    _view_states = None

try:
    from exondomaincompare.framework import analysis_router as _router  # type: ignore
except Exception:  # pragma: no cover
    _router = None


def _resolve_workflow(gene_symbol: Optional[str], mode: str = "auto") -> Dict[str, Any]:
    if _router is not None:
        try:
            return _router.resolve_gene_workflow(gene_symbol, mode=mode).to_dict()
        except Exception:
            pass
    sym = str(gene_symbol or "").strip().upper()
    validated = sym == "FGFR2"
    return {
        "gene_symbol": sym,
        "workflow": "validated_event_analysis" if validated else "shared_exploratory",
        "event_layer": "validated" if validated else "exploratory",
        "is_validated": validated,
        "has_event": validated,
        "support_level": "validated_event_analysis" if validated else "core_only_pilot",
        "case_study": "FGFR2_IIIb_IIIc" if validated else (
            f"{sym}_core_only_pilot" if sym else "core_only_pilot"),
        "analysis_id": "FGFR2_IIIb_IIIc" if validated else (
            f"{sym}_core_only_pilot" if sym else "core_only_pilot"),
        "creator": "run_pre_interpro_for_run.py" if validated else "run_core_gene_analysis.py",
        "mode": mode,
    }

_GENE_META_CACHE: Dict[str, Dict[str, Any]] = {}


_COLLECTION_STATUS_REL = Path("results") / "02_models" / "collection_status.json"


def _collection_status(run_dir: Optional[Path]) -> Dict[str, Any]:
    if run_dir is None:
        return {}
    data = read_json(run_dir / _COLLECTION_STATUS_REL, {}) or {}
    return data if isinstance(data, dict) else {}


def _precluster_failure(run_dir: Optional[Path], status: Dict[str, Any]
                        ) -> Dict[str, str]:
    contract = _collection_status(run_dir)
    raw = str(status.get("failed_reason") or status.get("error") or "").strip()
    message = str(contract.get("message") or "").strip()
    if message:
        return {
            "stage": "pre_cluster_data_acquisition",
            "cause": message,
            "collection_status": str(contract.get("status") or ""),
            "next_action": str(contract.get("next_action") or "retry_local_preparation"),
        }
    looks_like_traceback = raw.startswith("Traceback") or "  File \"" in raw
    return {
        "stage": "pre_cluster_data_acquisition" if raw else "",
        "cause": ("The run stopped before the gene and protein models were collected. "
                  "Download the diagnostics for the full log."
                  if looks_like_traceback or not raw else raw.splitlines()[0]),
        "collection_status": "",
        "next_action": "retry_local_preparation" if raw else "",
    }


def _gene_identity_for_run(run_dir: Optional[Path]) -> Dict[str, Any]:
    if run_dir is None:
        return {}
    rc = read_json(run_dir / "run_config.json", {}) or {}
    identity = rc.get("gene_identity")
    if isinstance(identity, dict) and identity:
        per_species = rc.get("gene_identity_by_species")
        return {**identity,
                "by_species": per_species if isinstance(per_species, dict) else {}}
    return {}


def _gene_meta_for_run(run_dir: Optional[Path]) -> Dict[str, Any]:
    fallback = {
        "analysis_id": "FGFR2_IIIb_IIIc",
        "gene_symbol": "FGFR2",
        "event_id": "FGFR2_IIIb_IIIc_cassette",
        "event_type": "mutually_exclusive_cassette",
        "event_display_name": "IIIb/IIIc cassette",
        "analysis_modes": {"core_gene_analysis": True, "event_analysis": "configured"},
        "event_analysis_mode": "configured",
        "has_event": True,
        "ui_labels": {
            "gene_explorer": "Gene Explorer",
            "event_region": "Cassette",
            "event_region_full": "IIIb/IIIc cassette",
            "event_discriminating_columns": "IIIb/IIIc-discriminating columns",
            "boundary_relation": "Boundary Consistency",
            "domain_relation_description": "Cassette-to-domain boundary consistency",
            "reference_comparison": "Human comparison",
        },
    }
    if _gene_config is None:
        return fallback
    key = str(run_dir) if run_dir is not None else "__default__"
    if key in _GENE_META_CACHE:
        return _GENE_META_CACHE[key]
    try:
        if run_dir is not None:
            rc = read_json(run_dir / "run_config.json", {}) or {}
            cfg = _gene_config.resolve_run_analysis(rc, run_dir)
        else:
            cfg = _gene_config.default_gene_config()
        meta = {
            "analysis_id": cfg.analysis_id,
            "gene_symbol": cfg.gene_symbol,
            "event_id": cfg.event_id,
            "event_type": cfg.event_type,
            "event_display_name": cfg.event_display_name,
            "analysis_modes": cfg.analysis_modes,
            "event_analysis_mode": cfg.event_analysis_mode,
            "has_event": cfg.has_event,
            "ui_labels": cfg.ui_labels,
        }
    except Exception:
        meta = fallback
    _GENE_META_CACHE[key] = meta
    return meta


app = FastAPI(title="ExonDomainCompare API", version=APP_VERSION)
app.add_middleware(
    CORSMiddleware,
    # accept the Vite dev server on any local port (5173, 5174, 4173, …)
    allow_origin_regex=r"http://(localhost|127\.0\.0\.1):\d+",
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

LOGGER = logging.getLogger("exondomaincompare.api")


def _validate_package_builder_capability() -> None:
    sys.path.insert(0, str(SCRIPTS_DIR / "shared_gene_analysis"))
    try:
        from exondomaincompare.shared_gene_analysis.package_builder import workbook_capability  # type: ignore
        ok, detail = workbook_capability()
    except Exception as exc:  # pragma: no cover - import path problem only
        ok, detail = False, str(exc)
    if ok:
        LOGGER.info("Package builder: workbook support available.")
    else:
        LOGGER.warning(
            "Package builder: workbook support unavailable — %s "
            "Install the backend requirements into this interpreter: "
            "%s -m pip install -r webapp/backend/requirements.txt",
            detail, sys.executable)


_validate_package_builder_capability()

PRESET_PILOT = ["Homo sapiens", "Mus musculus", "Gallus gallus", "Danio rerio"]
PRESET_VALIDATION = [
    "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Canis lupus familiaris",
    "Bos taurus", "Gallus gallus", "Anolis carolinensis", "Xenopus tropicalis",
    "Danio rerio", "Oryzias latipes",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def free_gb(path: Path) -> float:
    usage = shutil.disk_usage(path if path.exists() else PROJECT_ROOT)
    return round(usage.free / (1024 ** 3), 2)


def read_json(path: Path, fallback: Any = None) -> Any:
    if not Path(path).exists():
        if Path(path).name == "run_config.json" and (
                Path(path).parent / "run.json").is_file():
            try:
                return LegacyRunAdapter(Path(path).parent).config()
            except LegacyRunError:
                return fallback
        return fallback
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return fallback


def write_json(path: Path, data: Any) -> None:
    if isinstance(data, dict) and Path(path).name in {"run_config.json", "status.json"}:
        run_id = str(data.get("run_id") or Path(path).parent.name)
        data = stamp_payload(
            data,
            payload_type="run_config" if Path(path).name == "run_config.json" else "status",
            run_id=run_id, dataset_id=run_id,
            profile=RUNTIME_CONFIG.public_identity(),
            generator="webapp/backend/main.py",
        )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def _runtime_env() -> Dict[str, str]:
    env = os.environ.copy()
    env[RUNS_ENV] = str(RUNTIME_CONFIG.runs_root)
    env[LOCAL_PROFILE_ENV] = RUNTIME_CONFIG.local_profile_name
    env[LRZ_PROFILE_ENV] = RUNTIME_CONFIG.lrz_profile_name
    return env


def _public_run_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT))
    except ValueError:
        try:
            return "runs/" + str(path.resolve().relative_to(LOCAL_RUNS_ROOT.resolve()))
        except ValueError:
            return f"external:{path.name}"


def read_species_file(path: Path) -> List[str]:
    if not Path(path).exists():
        return []
    seen, out = set(), []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out




def parse_species_with_stats(text: Optional[str]) -> Dict[str, Any]:
    tokens: List[str] = []
    invalid: List[str] = []
    for raw in (text or "").replace(";", "\n").replace(",", "\n").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        token = cnr.normalize_species_token(s)
        if not token:
            continue
        tokens.append(token)
        if not cnr.SPECIES_ID_RE.match(token):
            invalid.append(s)
    counts: Dict[str, int] = {}
    for t in tokens:
        counts[t] = counts.get(t, 0) + 1
    seen, species = set(), []
    for t in tokens:
        if t not in seen:
            seen.add(t)
            species.append(t)
    duplicates = sorted({t for t, c in counts.items() if c > 1})
    return {
        "species": species,
        "count": len(species),
        "raw_count": len(tokens),
        "duplicates": duplicates,
        "invalid": sorted(set(invalid)),
        "preview": species[:60],
    }


def set_current_run(kind: str, run_dir: Path, label: str, run_id: str = "") -> Dict[str, Any]:
    root_id = ""
    if run_id and run_id != "example":
        try:
            record = resolve_run_record(RUNTIME_CONFIG, run_id)
            root_id = record.root_id if record else ""
        except RegistryError:
            root_id = ""
    ptr = {
        "schema_version": "2.0",
        "kind": kind,
        "dataset_id": "example" if run_id == "example" else f"run:{run_id}",
        "label": label,
        "run_id": run_id,
        "root_id": root_id,
        "selected_at": now_iso(),
    }
    write_json(CURRENT_RUN_PTR, ptr)
    return {**ptr, "run_dir": _public_run_path(run_dir)}


def get_current_run(*, public: bool = False) -> Optional[Dict[str, Any]]:
    ptr = read_json(CURRENT_RUN_PTR, None)
    if not ptr:
        ptr = read_json(LEGACY_CURRENT_RUN_PTR, None)
    if not ptr:
        return None
    ptr = dict(ptr)
    run_id = str(ptr.get("run_id") or "")
    if run_id == "example" or ptr.get("kind") == "example":
        run_dir = EXAMPLE_RUN_DIR
    elif run_id:
        try:
            record = resolve_run_record(RUNTIME_CONFIG, run_id)
        except RegistryError:
            return None
        if record is None:
            return None
        stored_root = str(ptr.get("root_id") or "")
        if stored_root and stored_root != record.root_id:
            return None
        run_dir = record.path
        ptr["root_id"] = record.root_id
    else:
        legacy = Path(str(ptr.get("run_dir") or ""))
        if not legacy.exists():
            return None
        run_dir = legacy
    ptr["run_dir"] = _public_run_path(run_dir) if public else str(run_dir)
    return ptr


def ensure_indices(run_dir: Path, rebuild: bool = False) -> Path:
    idx = Path(run_dir) / "website_indices"
    if rebuild or not (idx / "run_index.json").exists():
        bwi.write_all(Path(run_dir), idx)
    return idx


def current_index(name: str, rebuild: bool = False) -> Any:
    ptr = get_current_run()
    if not ptr:
        raise HTTPException(status_code=409, detail="No run loaded. Load the example freeze or start/open a run.")
    idx = ensure_indices(Path(ptr["run_dir"]), rebuild=rebuild)
    data = read_json(idx / name, None)
    if data is None:
        raise HTTPException(status_code=404, detail=f"Index not available: {name}")
    return data



DATASET_EXAMPLE = "example"

EXAMPLE_DERIVED_INDICES_DIR = (
    RESULTS_ROOT / "derived" / "website_indices"
    if RESULTS_ROOT == BUNDLED_FGFR2_ROOT
    else PROJECT_ROOT / "results" / "derived" / "example" / "website_indices"
)


def resolve_dataset(dataset: Optional[str]) -> Dict[str, Any]:
    if not dataset or dataset == DATASET_EXAMPLE:
        return {
            "kind": "example",
            "run_id": "example",
            "run_base": None,
            "closure_dir": EXAMPLE_RUN_DIR,
            "indices_dir": EXAMPLE_RUN_DIR / "website_indices",
            "derived_indices_dir": EXAMPLE_DERIVED_INDICES_DIR,
            "read_only": True,
        }
    run_id = dataset[4:] if dataset.startswith("run:") else dataset
    run_base = _safe_local_run_dir(run_id)
    if not run_base.is_dir():
        raise HTTPException(status_code=404, detail=f"Dataset run '{run_id}' not found.")
    record = resolve_run_record(RUNTIME_CONFIG, run_id)
    adapter = LegacyRunAdapter(run_base, expected_run_id=run_id)
    return {
        "kind": "run",
        "run_id": run_id,
        "run_base": run_base,
        "closure_dir": run_base / "results" / "13_final_pre_interpro_closure",
        "indices_dir": (
            run_base / "website" / "indices"
            if adapter.is_canonical and (run_base / "website" / "indices").is_dir()
            else run_base / "website_indices"),
        "derived_indices_dir": None,
        "read_only": bool(record.read_only) if record else False,
        "root_id": record.root_id if record else "",
        "layout_version": adapter.describe().layout_version,
    }


def dataset_index(name: str, dataset: Optional[str], rebuild: bool = False) -> Any:
    ds = resolve_dataset(dataset)
    idx = ds["indices_dir"]
    if ds["kind"] == "run":
        need = rebuild or not (idx / "run_index.json").exists()
        if need and (ds["closure_dir"] / "final_pre_interpro_truth_table.tsv").exists():
            try:
                bwi.write_all(ds["closure_dir"], idx)
            except Exception:
                pass  # partial run: serve whatever indices could be built
    derived = ds.get("derived_indices_dir")
    source_path = derived / name if derived else None
    data = read_json(source_path, None) if source_path else None
    if data is None:
        source_path = idx / name
        data = read_json(source_path, None)
    if data is None:
        raise HTTPException(status_code=404, detail=f"Index not available: {name}")
    if ds["kind"] == "run":
        expected = ds["run_id"]
        sidecar = read_json(source_path.parent / "_payload_contracts.json", None)
        if isinstance(sidecar, dict):
            valid, reason = verify_payload_contract(
                source_path, sidecar,
                expected_run_id=expected, expected_dataset_id=expected)
            if not valid:
                raise HTTPException(
                    status_code=409,
                    detail=f"Selected-run payload contract rejected {name}: {reason}")
        if isinstance(data, dict):
            metadata = data.get("_exondomain") if isinstance(data.get("_exondomain"), dict) else {}
            actual_run = str(data.get("run_id") or metadata.get("run_id") or "")
            actual_dataset = str(data.get("dataset_id") or metadata.get("dataset_id") or "")
            if actual_run and actual_run != expected:
                raise HTTPException(status_code=409, detail="Selected-run payload identity mismatch")
            if actual_dataset and actual_dataset not in {expected, f"run:{expected}"}:
                raise HTTPException(status_code=409, detail="Selected dataset payload identity mismatch")
            data.setdefault("run_id", expected)
            data.setdefault("dataset_id", expected)
    return _with_availability(data, name, ds)


_INDEX_ANALYSIS = {
    "msa_index.json": "protein_isoform_comparison",
    "isoform_alignment_index.json": "protein_isoform_comparison",
    "event_candidate_evidence_index.json": "protein_difference_candidate_analysis",
    "exon_domain_boundary_index.json": "boundary_analysis",
    "exon_domain_boundaries_index.json": "boundary_analysis",
    "domain_architecture_index.json": "domain_architecture",
}


def _with_availability(data: Any, name: str, ds: Mapping[str, Any]) -> Any:
    run_base = ds.get("run_base")
    if not run_base or not isinstance(data, dict):
        return data
    analysis = _INDEX_ANALYSIS.get(name)
    try:
        from exondomaincompare.shared_gene_analysis import analysis_availability as aa
        if not aa.has_core_tables(Path(run_base)):
            return _with_event_view_state(data, name, Path(run_base))
        if not analysis:
            return data
        state = aa.build_manifest(Path(run_base)).by_name().get(analysis)
    except Exception:
        return data
    if state is None:
        return data
    if state.status == aa.NOT_APPLICABLE:
        data["available"] = False
        data["status"] = aa.NOT_APPLICABLE
    elif data.get("available") is True:
        return data
    data["availability"] = aa.availability_block(state)
    return data


_EVENT_STATE_LABELS = {
    "pending": ("Pending cluster annotation", "pending"),
    "technically_missing": ("Expected output missing", "missing"),
    "stale": ("Rebuild required", "stale"),
    "failed": ("Analysis failed", "failed"),
    "scientifically_unavailable": ("No supported result", ""),
    "not_applicable": ("Not applicable", ""),
}


def _with_event_view_state(data: Dict[str, Any], name: str, run_base: Path) -> Dict[str, Any]:
    if data.get("available") is True:
        return data
    try:
        from exondomaincompare.shared_gene_analysis.run_availability import FGFR2_VIEWS, models_run
        if not models_run(run_base):
            return data
        requirement = next((v for v in FGFR2_VIEWS if v.index == name), None)
        if requirement is None:
            return data
        ready = _run_readiness(run_base, has_event=bool(
            _gene_meta_for_run(run_base).get("has_event", True)))
        state = next((s for s in (ready.views if ready else []) if s.view == requirement.view),
                     None)
    except Exception:
        return data
    if state is None or state.available:
        return data
    label, badge = _EVENT_STATE_LABELS.get(state.state, ("", ""))
    if not label:
        return data
    data["status"] = state.state
    data["availability"] = {
        "state": state.state,
        "label": label,
        "reason": state.reason,
        "reason_code": state.state,
        "badge": badge,
        "prerequisite_name": "",
        "prerequisite_count": None,
    }
    return data


def index_for_request(name: str, dataset: Optional[str], rebuild: bool = False) -> Any:
    if dataset:
        return dataset_index(name, dataset, rebuild=rebuild)
    return current_index(name, rebuild=rebuild)


def _ensure_run_indices(run_dir: Path, rebuild: bool = False) -> None:
    closure = run_dir / "results" / "13_final_pre_interpro_closure"
    idx = run_dir / "website_indices"
    if not (closure / "final_pre_interpro_truth_table.tsv").exists():
        return
    if rebuild or not (idx / "run_index.json").exists():
        try:
            bwi.write_all(closure, idx)
        except Exception:
            pass


def _count_fasta_records(faa: Path) -> int:
    if not Path(faa).is_file():
        return 0
    try:
        n = 0
        with open(faa, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
        return n
    except Exception:
        return 0


HUMAN_REFERENCE_CACHE = WEB_STATE_DIR / "human_reference_control.json"


def build_human_reference_control(rebuild: bool = False) -> Dict[str, Any]:
    if not rebuild:
        cached = read_json(HUMAN_REFERENCE_CACHE, None)
        if cached:
            return cached
    species_idx = read_json(EXAMPLE_RUN_DIR / "website_indices" / "species_index.json", []) or []
    human_row = next((r for r in species_idx
                      if isinstance(r, dict) and r.get("species", "").lower() == "homo_sapiens"), None)
    isoforms: Dict[str, Any] = {}
    if human_row:
        for iso in human_row.get("isoforms", []):
            key = (iso.get("final_isoform_label") or iso.get("isoform") or "").strip()
            if key:
                isoforms[key] = {
                    "protein_id": iso.get("protein_id", ""),
                    "transcript_id": iso.get("transcript_id", ""),
                    "protein_length": iso.get("protein_length"),
                    "sequence_md5": iso.get("sequence_md5", ""),
                    "interpro_included": iso.get("interpro_included", ""),
                }

    cassette_block = None
    ex_cass = read_json(EXAMPLE_RUN_DIR / "website_indices" / "cassette_residue_index.json", {}) or {}
    human_ref_residues = ex_cass.get("human_reference") or {}
    if human_ref_residues.get("IIIb") or human_ref_residues.get("IIIc"):
        cassette_block = {
            "panels": ["IIIb", "IIIc"],
            "human_reference": {
                "IIIb": human_ref_residues.get("IIIb", []),
                "IIIc": human_ref_residues.get("IIIc", []),
            },
            "discriminating": ex_cass.get("discriminating") or {},
            "cassette_length": ex_cass.get("cassette_length") or {},
            "source_table": (ex_cass.get("source_tables") or {}).get("residue_map", ""),
            "note": ("Human IIIb/IIIc cassette residues from the validated example dataset; "
                     "reused as reference/control only."),
        }

    data = {
        "source": "validated_example_dataset",
        "source_dataset": "FGFR2 IIIb/IIIc — 30 vertebrates (freeze)",
        "source_path": str(EXAMPLE_RUN_DIR.relative_to(PROJECT_ROOT)),
        "validation_status": "validated" if isoforms else "unavailable",
        "species": "homo_sapiens",
        "display_species_name": "Homo sapiens",
        "isoforms": isoforms,
        "cassette": cassette_block,
        "note": ("Curated human FGFR2 IIIb/IIIc from the validated example dataset. "
                 "Reused by custom runs as a fixed reference/control layer only; "
                 "not counted as an analysed species unless explicitly selected."),
        "built_at": now_iso(),
    }
    try:
        write_json(HUMAN_REFERENCE_CACHE, data)
    except Exception:
        pass
    return data


_CALM_STATUS = {
    "created": "Ready to prepare FASTA",
    "pre_interpro_running": "Preparing FASTA",
    "cluster_required": "Cluster input ready",
    "cluster_running": "Cluster annotation running",
    "cluster_fetch_complete": "Cluster results fetched · Post-InterPro required",
    "post_interpro_running": "Post-InterPro analysis running",
    "results_ready": "Results ready",
    "failed": "Needs attention",
    "stopped": "Stopped by user",
}
_NEXT_ACTION = {
    "created": ("run_pre_interpro", "Prepare FASTA (local pre-InterPro)"),
    "pre_interpro_running": ("wait_pre_interpro", "Preparing FASTA…"),
    "cluster_required": ("run_cluster_roundtrip_command", "Run cluster annotation (one terminal command)"),
    "cluster_running": ("wait_cluster", "Cluster annotation running on LRZ"),
    "cluster_fetch_complete": ("run_post_interpro", "Run Post-InterPro analysis"),
    "post_interpro_running": ("wait_post_interpro", "Post-InterPro analysis running…"),
    "results_ready": ("open_results", "Open results"),
    "failed": ("inspect_logs", "Needs attention — see technical logs"),
    "stopped": ("restart_or_delete", "Restart or delete this run"),
}


_CORE_STATUS_LABELS = {
    "created_not_started": "Created — not started",
    "running": "Running — pre-InterPro pipeline",
    "core_model_collection_failed": "Model collection failed",
    "incomplete": "Incomplete — outputs missing",
    "cluster_required": "Cluster annotation required",
    "cluster_running": "Cluster annotation running",
    "post_interpro_incomplete": "Post-InterPro incomplete",
    "results_ready": "Results ready",
}
_CORE_NEXT_LABELS = {
    "run_core_collection": "Re-run core model collection",
    "wait_pre_interpro": "Pre-InterPro pipeline running…",
    "inspect_logs": "Needs attention — see logs / missing outputs",
    "run_cluster_roundtrip_command": "Run cluster annotation (one terminal command)",
    "wait_cluster": "Cluster annotation running on LRZ",
    "run_core_post": "Build domain architecture from cluster outputs",
    "open_results": "Open results",
}


def _run_species_count(run_dir: Path) -> int:
    slist = run_dir / "species_list.txt"
    if slist.is_file():
        n = len([ln for ln in slist.read_text(encoding="utf-8").splitlines() if ln.strip()])
        if n:
            return n
    psel = read_json(run_dir / "website_indices" / "primary_selection_index.json", {}) or {}
    return len(psel.get("species_primaries") or []) or int(psel.get("n_species") or 1) or 1


def _derive_core_status_model(run_dir: Path) -> Dict[str, Any]:
    rep = _evaluate_core_run(run_dir)
    run_id = run_dir.name
    wi = run_dir / "website_indices"
    gene_meta = _gene_meta_for_run(run_dir)
    inferred = rep["inferred_status"]

    gav = read_json(wi / "generic" / "available_views.json", {}) or {}
    gviews = gav.get("available_views", {}) if isinstance(gav, dict) else {}
    results_ready = inferred == "results_ready"
    # Only expose data views when the underlying core outputs actually exist.
    have_models = "model_collection" in rep["completed_milestones"]
    available_views = {
        "overview": have_models,
        "gene_explorer": have_models and bool(gviews.get("gene_explorer") or gviews.get("gene_models")),
        "domain_architecture": results_ready and bool(gviews.get("domain_architecture")),
        "exon_domain_boundaries": results_ready and bool(gviews.get("exon_domain_boundaries")),
        "synteny": have_models and bool(gviews.get("synteny")),
        "figure_gallery": bool(gviews.get("figure_gallery") or gviews.get("figures")),
        # event-specific views are never available for a core-only (no-event) run
        "boundary_consistency": False,
        "event_region": False,
    }


    capability = read_json(wi / "generic" / "gene_capability_report.json", None)
    counts = rep["counts"]
    if not isinstance(capability, dict) or not capability:
        capability = {
            "core_model_collection": "available" if have_models else "failed",
            "protein_isoforms_count": counts.get("protein_isoforms", 0),
            "primary_proteins_count": counts.get("primary_proteins", 0),
            "gene_models_count": counts.get("gene_models", 0),
            "species_count": _run_species_count(run_dir),
            "exon_map": "available" if counts.get("exon_map_rows", 0) > 0 else "unavailable",
            "synteny": "available" if counts.get("synteny_neighbors", 0) > 0 else "unavailable",
            "synteny_neighbours_count": counts.get("synteny_neighbors", 0),
            "cluster_status": "complete" if results_ready else "pending",
            "domain_architecture": "available" if results_ready else "pending",
            "exon_domain_boundaries": "available" if results_ready else "pending",
            "event_configured": False,
            "exploratory_event_evidence": "none",
            "candidate_clusters_count": 0,
            "event_analysis_enabled": False,
            "support_level": "core_only_pilot",
            "cluster_command": rep["cluster_command"],
        }

    next_action = rep["suggested_next_action"]

    _STAGE = {
        "created_not_started": ("created", "Created"),
        "running": ("running", "Running — pre-InterPro pipeline"),
        "core_model_collection_failed": ("models_failed", "Model collection failed"),
        "incomplete": ("incomplete", "Incomplete"),
        "cluster_required": ("cluster_required", "Models ready · cluster required"),
        "cluster_running": ("cluster_submitted", "Cluster annotation running"),
        "post_interpro_incomplete": ("cluster_complete", "Cluster complete · post-domain pending"),
        "results_ready": ("results_ready", "Results ready"),
    }
    stage, stage_label = _STAGE.get(inferred, ("cluster_required", "Models ready"))
    return {
        "run_id": run_id,
        "stage": stage,
        "stage_label": stage_label,
        "analysis_id": gene_meta["analysis_id"],
        "gene_symbol": gene_meta["gene_symbol"],
        "gene_identity": _gene_identity_for_run(run_dir),
        "event_id": gene_meta["event_id"],
        "event_type": gene_meta["event_type"],
        "event_display_name": gene_meta["event_display_name"],
        "analysis_modes": gene_meta.get("analysis_modes",
                                        {"core_gene_analysis": True, "event_analysis": "optional"}),
        "event_analysis_mode": gene_meta.get("event_analysis_mode", "optional"),
        "has_event": False,
        "event_message": ("No event region is configured or detected for this gene. "
                          "Core gene-level analysis is available."),
        "ui_labels": gene_meta["ui_labels"],
        "status": inferred,
        "status_label": _CORE_STATUS_LABELS.get(inferred, inferred),
        "current_step": rep.get("detail") or rep.get("run_mode", "core_only_pilot"),
        "next_action": next_action,
        "next_action_label": _CORE_NEXT_LABELS.get(next_action, next_action),
        "last_error": rep.get("failed_reason") or "",
        "failed_stage": rep.get("failed_stage") or "",
        "failed_species": rep.get("failed_species") or "",
        "pre_interpro_status": (
            "running" if inferred == "running"
            else "failed" if inferred in ("core_model_collection_failed", "incomplete")
            else "complete" if have_models else "not_started"),
        "primary_fasta_status": "available" if rep["counts"]["primary_proteins"] > 0 else "not_available",
        "primary_fasta_count": rep["counts"]["primary_proteins"],
        "review_fasta_count": 0,
        "post_interpro_status": "complete" if results_ready else "not_started",
        "cluster_analysis_status": (read_json(run_dir / "status.json", {}) or {}).get(
            "cluster_analysis_status", "not_started"),
        "cluster_fetch_status": (read_json(run_dir / "status.json", {}) or {}).get(
            "cluster_fetch_status", "not_started"),
        "cluster_required": inferred == "cluster_required",
        "cluster_command": rep["cluster_command"],
        "explorable": available_views["overview"],
        "available_views": available_views,
        "human_reference": {},
        # Core-only extras: honest milestone detail for the dashboard.
        "core_only": True,
        "pipeline_type": "core_gene_pipeline",
        # Shared-pipeline routing: every gene runs the same conceptual
        # pipeline; only the event layer differs. A core-only gene uses the
        # exploratory event-evidence layer and has no validated event.
        "event_layer_type": "exploratory_event_evidence",
        "has_validated_event": False,
        "experimental": True,
        "support_level": capability.get("support_level", "core_only_pilot"),
        "capability": capability,
        "milestones": rep["milestones"],
        "completed_milestones": rep["completed_milestones"],
        "missing_required": rep["missing_required"],
        "missing_optional": rep["missing_optional"],
        "counts": rep["counts"],
        "logs": rep["logs"],
    }


def derive_status_model(run_dir: Path) -> Dict[str, Any]:
    if _evaluate_core_run is not None and _is_core_only_run is not None:
        try:
            if _is_core_only_run(run_dir):
                return _derive_core_status_model(run_dir)
        except Exception:
            pass  # fall back to the generic model on any error

    st = read_json(run_dir / "status.json", {}) or {}
    files = _local_run_file_checks(run_dir)
    closure = run_dir / "results" / "13_final_pre_interpro_closure"
    primary_faa = closure / "freeze" / "final_pre_interpro_proteins_primary.faa"
    review_faa = closure / "freeze" / "final_pre_interpro_proteins_all_review_included.faa"
    run_id = run_dir.name
    wi = run_dir / "website_indices"

    proc_running = _preinterpro_running(run_dir) is not None
    post_proc_running = _postinterpro_running(run_dir)
    overall = (st.get("status") or "").lower()
    pre = (st.get("pre_interpro_status") or "not_started").lower()
    post = (st.get("post_interpro_status") or "not_started").lower()
    cluster = (st.get("cluster_analysis_status") or "").lower()
    primary_ok = files["primary_fasta"]
    ips_ok = files["interproscan_output"]
    tm_ok = files["pytmhmm_output"]
    indices_ok = files["website_indices"]
    last_error = st.get("failed_reason") or st.get("error") or ""
    post_complete = post in ("complete", "completed")
    failure = _precluster_failure(run_dir, st)
    pre_ok = primary_ok or pre in ("complete", "completed", "done")
    failed = (overall in ("failed", "error")
              or (bool(last_error) and not post_complete and not pre_ok
                  and overall not in ("complete", "running")))

    def has_idx(name: str) -> bool:
        return (wi / name).exists()

    def idx_available(name: str) -> bool:
        p = wi / name
        if not p.exists():
            return False
        data = read_json(p, {}) or {}
        if isinstance(data, dict) and "available" in data:
            return bool(data["available"])
        return True

    gene_meta = _gene_meta_for_run(run_dir)
    has_event = bool(gene_meta.get("has_event", True))

    domain_ok = idx_available("species_domain_architecture.json")
    available_views = {
        # --- core (gene-agnostic) views: never require a configured event region ---
        "overview": has_idx("run_index.json") or files["truth_table"],
        "gene_explorer": has_idx("species_index.json"),
        "domain_architecture": domain_ok,
        # generic all-exon exon-domain boundary view: derivable wherever domain
        # architecture (domains + exon map) is available.
        "exon_domain_boundaries": domain_ok,
        "synteny": idx_available("synteny_locus_index.json"),
        "figure_gallery": has_idx("figure_index.json"),
        # --- event-specific views: only when an event region is configured ---
        # (FGFR2 has_event=True -> unchanged; core-only genes hide these.)
        "boundary_consistency": has_event and idx_available("boundary_consistency_summary.json"),
        "event_region": has_event and idx_available("cassette_residue_index.json"),
    }

    # Core-only runs publish availability through their generic indices.
    if not has_event:
        gav = read_json(wi / "generic" / "available_views.json", {}) or {}
        gviews = gav.get("available_views", {}) if isinstance(gav, dict) else {}
        if gviews:
            available_views["overview"] = available_views["overview"] or bool(gviews.get("overview"))
            available_views["gene_explorer"] = (available_views["gene_explorer"]
                                                or bool(gviews.get("gene_explorer"))
                                                or bool(gviews.get("gene_models")))
            available_views["domain_architecture"] = (available_views["domain_architecture"]
                                                       or bool(gviews.get("domain_architecture")))
            available_views["exon_domain_boundaries"] = (available_views["exon_domain_boundaries"]
                                                         or bool(gviews.get("exon_domain_boundaries")))
            available_views["synteny"] = available_views["synteny"] or bool(gviews.get("synteny"))

    # On-disk artefacts override stale persisted completion fields.
    ready = _run_readiness(run_dir, has_event=has_event)
    view_reasons: Dict[str, Any] = {}
    if ready is not None:
        view_reasons = ready.as_dict()["views"]
        for name, node in view_reasons.items():
            if name in available_views:
                # An index that exists but is stale or empty is not an available view.
                available_views[name] = bool(node["available"])
            else:
                available_views[name] = bool(node["available"])

    stopped = (overall == "stopped" or pre == "stopped" or post == "stopped")
    if ready is not None:
        results_ready = ready.ready
    else:
        results_ready = post_complete or (indices_ok
                                          and available_views["boundary_consistency"])
    if results_ready:
        status = "results_ready"
    elif proc_running or (pre == "running" and not stopped):
        status = "pre_interpro_running"
    elif post_proc_running or (post == "running" and not stopped):
        status = "post_interpro_running"
    elif stopped:
        status = "stopped"
    elif failed:
        status = "failed"
    elif post == "running":
        status = "post_interpro_running"
    elif ips_ok and tm_ok:
        status = "cluster_fetch_complete"
    elif cluster in ("submitted", "running", "queued", "pending"):
        status = "cluster_running"
    elif primary_ok or pre in ("complete", "completed", "done"):
        status = "cluster_required"
    else:
        status = overall or "created"

    next_action, next_label = _NEXT_ACTION.get(status, ("", ""))
    primary_count = _count_fasta_records(primary_faa) if primary_ok else int(st.get("primary_fasta_count") or 0)
    # The cluster round-trip is required whenever pre-InterPro produced the primary FASTA
    # but the post-InterPro (domain/boundary) layer is not complete yet. Independent of the
    # exact calm status label so the UI can always decide whether to surface the command.
    cluster_required = bool(primary_ok and not post_complete)

    event_message = ("" if has_event else
                     "No event region is configured or detected for this gene. "
                     "Core gene-level analysis is still available.")

    return {
        "run_id": run_id,
        "analysis_id": gene_meta["analysis_id"],
        "gene_symbol": gene_meta["gene_symbol"],
        "gene_identity": _gene_identity_for_run(run_dir),
        "event_id": gene_meta["event_id"],
        "event_type": gene_meta["event_type"],
        "event_display_name": gene_meta["event_display_name"],
        "analysis_modes": gene_meta.get("analysis_modes",
                                        {"core_gene_analysis": True, "event_analysis": "configured"}),
        "event_analysis_mode": gene_meta.get("event_analysis_mode", "configured"),
        "has_event": has_event,
        "event_message": event_message,
        "ui_labels": gene_meta["ui_labels"],
        "core_only": False,
        "pipeline_type": ("validated_event_pipeline" if has_event else "core_gene_pipeline"),
        # Shared-pipeline routing. FGFR2 (has_event) carries the
        # validated IIIb/IIIc event layer; any other gene uses the exploratory
        # event-evidence layer. The frontend routes on these, not on isCoreOnly.
        "event_layer_type": ("validated_fgfr2_iiib_iiic" if has_event else "exploratory_event_evidence"),
        "has_validated_event": bool(has_event),
        "status": status,
        "status_label": _CALM_STATUS.get(status, status),
        "current_step": st.get("current_step", ""),
        "next_action": next_action,
        "next_action_label": next_label,
        # The recorded cause, not the traceback that noticed it.
        "last_error": (failure["cause"] if failed else ""),
        "failed_stage": (failure["stage"] if failed else ""),
        "collection_status": failure["collection_status"],
        "can_retry_local_preparation": bool(
            failed and failure["next_action"] == "retry_local_preparation"),
        "primary_fasta_status": "available" if primary_ok else "not_available",
        "primary_fasta_count": primary_count,
        "review_fasta_count": (_count_fasta_records(review_faa)
                               if review_faa.exists() else int(st.get("review_fasta_count") or 0)),
        "post_interpro_status": post,
        "cluster_analysis_status": cluster or "not_started",
        "cluster_fetch_status": (st.get("cluster_fetch_status") or "not_started"),
        "cluster_required": cluster_required,
        "cluster_command": (
            _cluster_roundtrip_command(run_id)
            if primary_ok else ""),
        "explorable": available_views["overview"],
        "available_views": available_views,
        # Why each view is or is not available, so an empty page can name the cause
        # instead of implying that the species lacks the biology.
        "view_availability": view_reasons,
        "readiness": (ready.as_dict() | {"views": {}} if ready is not None else {}),
        "human_reference": st.get("human_reference") or {},
    }


def _run_readiness(run_dir: Path, has_event: bool = True):
    if _readiness is None:
        return None
    try:
        st = read_json(run_dir / "status.json", {}) or {}
        n_species = int(st.get("species_count") or 0)
        primary = (run_dir / "results" / "13_final_pre_interpro_closure" / "freeze"
                   / "final_pre_interpro_proteins_primary.faa")
        return _readiness(run_dir, n_species=n_species, has_event=has_event,
                          pre_interpro_complete=primary.exists())
    except Exception:
        return None


def closure_dir_for(run_base: Path) -> Path:
    return Path(run_base) / "13_final_pre_interpro_closure"


# --- run status (jobs) ----------------------------------------------------- #
def run_state_dir(run_id: str) -> Path:
    return WEB_RUNS_ROOT / run_id


def run_status_path(run_id: str) -> Path:
    return run_state_dir(run_id) / "web_status.json"


def run_log_path(run_id: str) -> Path:
    return run_state_dir(run_id) / "web_run.log"






def list_runs() -> List[Dict[str, Any]]:
    runs: List[Dict[str, Any]] = []
    if EXAMPLE_RUN_DIR.exists():
        rm = read_json(EXAMPLE_RUN_DIR / "final_pre_interpro_run_mode.json", {}) or {}
        runs.append({
            "kind": "example", "run_id": rm.get("run_id", "example"),
            "label": "Example · FGFR2 IIIb/IIIc — 30 vertebrates",
            "run_dir": _public_run_path(EXAMPLE_RUN_DIR), "status": "finished",
        })
    for p in sorted(WEB_RUNS_ROOT.glob("run_*"), reverse=True):
        st = read_json(p / "web_status.json", None)
        if st:
            runs.append({**st, "kind": st.get("kind", "job"),
                         "run_dir": _public_run_path(closure_dir_for(p))})
    return runs


# Request models.
class RunRequest(BaseModel):
    case_study: str = Field(default="FGFR2 IIIb/IIIc")
    species_source: str = Field(default="full_reference")  # full_reference|pilot|validation|custom
    species_list_text: Optional[str] = None
    preset_name: Optional[str] = None
    live_refresh: bool = False
    debug: bool = False


class OpenRunRequest(BaseModel):
    path: str


class ParseRequest(BaseModel):
    text: str = ""


# --------------------------------------------------------------------------- #
# meta endpoints
# --------------------------------------------------------------------------- #
@app.get("/api/health")
def health() -> Dict[str, Any]:
    cur = get_current_run(public=True)
    return {
        "status": "ok", "service": "ExonDomainCompare API", "version": APP_VERSION,
        "project_root": ".",
        "example_available": EXAMPLE_RUN_DIR.exists(),
        "shared_ncbi_cache_exists": SHARED_NCBI_CACHE.exists(),
        "runner_available": RUNNER_SCRIPT.exists(),
        "free_disk_gb": free_gb(PROJECT_ROOT),
        "current_run": cur,
    }


@app.get("/api/presets")
def presets() -> Dict[str, Any]:
    ref = read_species_file(REFERENCE_SPECIES_LIST)
    out = {
        "full_reference": {
            "id": "full_reference",
            "label": f"Full reference panel — {len(ref) or 30} vertebrates",
            "file": str(REFERENCE_SPECIES_LIST.relative_to(PROJECT_ROOT)) if REFERENCE_SPECIES_LIST.exists() else "",
            "available": REFERENCE_SPECIES_LIST.exists(),
            "species": ref,
        },
        "pilot": {"id": "pilot", "label": f"Pilot panel — {len(PRESET_PILOT)} vertebrates",
                  "available": True, "species": PRESET_PILOT},
        "validation": {"id": "validation", "label": f"Validation panel — {len(PRESET_VALIDATION)} vertebrates",
                       "available": True, "species": PRESET_VALIDATION},
    }
    return out


@app.post("/api/species/parse")
def species_parse(req: ParseRequest) -> Dict[str, Any]:
    return parse_species_with_stats(req.text)


# --------------------------------------------------------------------------- #
# run selection / lifecycle
# --------------------------------------------------------------------------- #
@app.get("/api/runs")
def runs() -> List[Dict[str, Any]]:
    return list_runs()


@app.get("/api/runs/current")
def runs_current() -> Dict[str, Any]:
    cur = get_current_run(public=True)
    if not cur:
        raise HTTPException(status_code=404, detail="No run loaded.")
    return cur


# Local run discovery never executes cluster or remote commands.
LOCAL_RUNS_ROOT = RUNTIME_CONFIG.runs_root


def _run_species_ids(run_dir: Path, cfg: Dict[str, Any]) -> List[str]:
    ids = [s for s in (cfg.get("species_ids") or []) if s]
    if ids:
        return ids
    try:
        adapted = LegacyRunAdapter(run_dir).species()
        if adapted:
            return adapted
    except LegacyRunError:
        pass
    listing = run_dir / "species_list.txt"
    if listing.is_file():
        try:
            return [ln.strip() for ln in listing.read_text(encoding="utf-8").splitlines()
                    if ln.strip() and not ln.startswith("#")]
        except OSError:
            return []
    return []


def _local_run_summary(run_dir: Path) -> Dict[str, Any]:
    adapter = LegacyRunAdapter(run_dir)
    cfg = adapter.config()
    status = adapter.status()
    description = adapter.describe()
    model = derive_status_model(run_dir)
    run_id = cfg.get("run_id", run_dir.name)
    species_ids = _run_species_ids(run_dir, cfg)
    species_count = cfg.get("species_count", status.get("species_count", 0)) or len(species_ids)
    # A legacy run has no usable run_name; its title is derived rather than
    # written back, so a validated result file is never rewritten for a label.
    labels = run_labels.run_display_fields(
        cfg, run_id=run_id, species_count=species_count,
        gene_symbol=model.get("gene_symbol") or cfg.get("gene_symbol") or "",
        species=species_ids, status=model["status"])
    return {
        "run_id": run_id,
        **labels,
        "species_ids": species_ids,
        "run_dir": _public_run_path(run_dir),
        "species_count": species_count,
        "case_study": cfg.get("case_study", ""),
        "created_at": cfg.get("created_at", ""),
        "completion_summary": run_labels.completion_summary(
            primary_fasta_count=model["primary_fasta_count"],
            available_views=model["available_views"],
            species_count=species_count),
        "failure_summary": (run_labels.describe_failure(
            failed_stage=model.get("failed_stage", ""),
            failed_species=model.get("failed_species", ""),
            last_error=model.get("last_error", ""))
            if run_labels.run_group(model["status"]) == "attention" else ""),
        "status": model["status"],
        "status_label": model["status_label"],
        "current_step": model["current_step"],
        "next_action": model["next_action"],
        "next_action_label": model["next_action_label"],
        "primary_fasta_count": model["primary_fasta_count"],
        "explorable": model["explorable"],
        "available_views": model["available_views"],
        "human_reference": model["human_reference"],
        "gene_symbol": model["gene_symbol"],
        "analysis_id": model["analysis_id"],
        "has_event": model["has_event"],
        "event_message": model["event_message"],
        "experimental": bool(cfg.get("experimental") or status.get("experimental")),
        "support_level": cfg.get("support_level") or status.get("support_level") or "",
        "run_mode": cfg.get("run_mode") or status.get("run_mode") or "",
        "core_only": bool(model.get("core_only")),
        "pipeline_type": model.get("pipeline_type", "validated_event_pipeline"),
        "event_layer_type": model.get("event_layer_type",
                                      "validated_fgfr2_iiib_iiic" if model.get("has_event", True)
                                      else "exploratory_event_evidence"),
        "has_validated_event": bool(model.get("has_validated_event", model.get("has_event", True))),
        "missing_required": model.get("missing_required", []),
        "last_error": model.get("last_error", ""),
        "failed_stage": model.get("failed_stage", ""),
        "failed_species": model.get("failed_species", ""),
        "stage": model.get("stage", ""),
        "stage_label": model.get("stage_label", ""),
        "pre_interpro_status": model.get("pre_interpro_status", ""),
        "has_config": bool(cfg),
        "has_status": (run_dir / "status.json").exists(),
        "has_readme": (run_dir / "00_README_NEXT_STEPS.md").exists(),
        "layout_version": description.layout_version,
        "legacy_adapter": not adapter.is_canonical,
        "legacy_partial": description.partial,
    }


def _local_run_summaries(*, include_bundled: bool) -> List[Dict[str, Any]]:
    records, collisions = discover_runs(RUNTIME_CONFIG)
    runs_out = []
    for record in records:
        if record.kind == "bundled_example" and not include_bundled:
            continue
        try:
            summary = _local_run_summary(record.path)
        except LegacyRunError:
            continue
        summary.update({
            "root_id": record.root_id,
            "root_kind": record.kind,
            "read_only": record.read_only,
            "explicit_binding": record.explicit,
        })
        runs_out.append(summary)
    for run_id, candidates in collisions.items():
        runs_out.append({
            "run_id": run_id,
            "run_name": run_id,
            "label": f"{run_id} · storage collision",
            "status": "collision",
            "status_label": "Needs storage selection",
            "collision": True,
            "collision_roots": [row.root_id for row in candidates],
            "explorable": False,
            "read_only": True,
        })
    return run_labels.sort_runs(runs_out)


@app.get("/api/local-runs")
def local_runs() -> List[Dict[str, Any]]:
    return _local_run_summaries(include_bundled=False)


@app.get("/api/local-runs/{run_id}")
def local_run_detail(run_id: str) -> Dict[str, Any]:
    run_dir = _safe_local_run_dir(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")
    adapter = LegacyRunAdapter(run_dir, expected_run_id=run_id)
    record = resolve_run_record(RUNTIME_CONFIG, run_id)
    return {
        "run_id": run_id,
        "run_dir": _public_run_path(run_dir),
        "config": adapter.config(),
        "status": adapter.status(),
        "readme_available": (run_dir / "00_README_NEXT_STEPS.md").exists(),
        "layout": adapter.describe().__dict__,
        "read_only": bool(record.read_only) if record else False,
        "root_id": record.root_id if record else "",
    }


# Workflow endpoints create folders; cluster commands remain terminal-only.
class CreateLocalRunRequest(BaseModel):
    # The visible label. Optional: when it is empty the run is titled from its
    # gene and species instead, and the run_id slug is derived the same way.
    run_name: str = Field(default="")
    preset: Optional[str] = None            # "full30" | "pilot"
    species_text: Optional[str] = None
    species_file_content: Optional[str] = None
    case_study: str = Field(default="FGFR2_IIIb_IIIc")
    # Gene-aware New Run form. The gene symbol drives the central analysis router:
    # FGFR2 → validated (frozen) workflow, every other gene → shared exploratory
    # workflow. Defaults to FGFR2 so existing callers keep working unchanged.
    gene_symbol: Optional[str] = Field(default=None)
    mode: str = Field(default="auto")       # "auto" (only mode for now)


def _safe_local_run_dir(run_id: str) -> Path:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", run_id or ""):
        raise HTTPException(status_code=400, detail="Invalid run id.")
    try:
        record = resolve_run_record(RUNTIME_CONFIG, run_id)
    except RunCollisionError as exc:
        raise HTTPException(
            status_code=409,
            detail=(
                f"Run '{run_id}' exists in multiple registered roots; "
                "select an explicit registry binding."),
        ) from exc
    except RegistryError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    return record.path if record else (LOCAL_RUNS_ROOT / run_id).resolve()


def _require_writable_run(run_id: str) -> Path:
    run_dir = _safe_local_run_dir(run_id)
    try:
        record = resolve_run_record(RUNTIME_CONFIG, run_id)
    except RegistryError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    if record and record.read_only:
        raise HTTPException(
            status_code=409,
            detail="This registered legacy run is read-only. Copy it before retry/resume.")
    return run_dir


def _cluster_roundtrip_command(run_id: str) -> str:
    argv = [".venv/bin/edc", "cluster", "roundtrip", "--run-id", run_id]
    assignments = []
    for name in (DATA_ENV, CONFIG_ENV, RUNS_ENV, LOCAL_PROFILE_ENV, LRZ_PROFILE_ENV):
        value = os.environ.get(name)
        if not value:
            continue
        if name in (DATA_ENV, CONFIG_ENV, RUNS_ENV):
            value = str(Path(value).expanduser().resolve())
        assignments.append(f"{name}={value}")
    if assignments:
        argv = ["env", *assignments, *argv]
    return RUNTIME_CONFIG.command(argv)


def _local_run_commands(run_id: str, run_dir: Path) -> Dict[str, Any]:
    # Workflow-aware next steps: the validated (FGFR2) path and the shared
    # exploratory path share the SAME cluster round-trip, but differ in the
    # pre/post local steps (pre-InterPro wrapper vs. the gene-agnostic core
    # runner). The gene symbol is resolved through the central router.
    rc = read_json(run_dir / "run_config.json", {}) or {}
    wf = _resolve_workflow(rc.get("gene_symbol") or "FGFR2")
    shared = not wf["is_validated"]
    gene = wf["gene_symbol"]
    def command(*parts: str, cluster: bool = False) -> str:
        # Keep the accepted copy/paste command stable. The selected profile is
        # returned separately and resolved by the same user config at execution.
        return RUNTIME_CONFIG.command(["python", *parts])

    if shared:
        commands = [
            {"id": "core", "title": "1 · Core gene analysis (pre-cluster)",
             "command": command("-m", "exondomaincompare.framework.run_core_gene_analysis",
                                "--gene", gene, "--species", "<species>",
                                "--input-mode", "auto"),
             "explanation": "Build the gene-agnostic core contract, exploratory evidence and "
                            "pre-cluster indices for this run.",
             "cluster": False},
        ]
    else:
        commands = [
            {"id": "pre_interpro", "title": "1 · Pre-InterPro pipeline",
             "command": command("scripts/run_pre_interpro_for_run.py",
                                "--run-id", run_id),
             "explanation": "Run the local pre-InterPro pipeline; produces the primary/review FASTA "
                            "in this run's freeze folder.",
             "cluster": False},
        ]
    commands += [
        {"id": "submit", "title": "2 · Submit cluster jobs",
         "command": command("scripts/interpro_cluster/submit_cluster_analysis.py",
                            "--run-id", run_id, cluster=True),
         "explanation": "Submit InterProScan + pyTMHMM on the cluster from YOUR local terminal. "
                        "Login / 2FA happen locally; the webapp never sees them.",
         "cluster": True},
        {"id": "check", "title": "3 · Check cluster jobs",
         "command": command("scripts/interpro_cluster/check_cluster_analysis.py",
                            "--run-id", run_id, cluster=True),
         "explanation": "Poll the cluster job status from your local terminal.",
         "cluster": True},
        {"id": "fetch", "title": "4 · Fetch cluster outputs",
         "command": command("scripts/interpro_cluster/fetch_cluster_analysis.py",
                            "--run-id", run_id, cluster=True),
         "explanation": "Download InterProScan + pyTMHMM outputs into this run folder.",
         "cluster": True},
    ]
    if shared:
        commands.append(
            {"id": "post_core", "title": "5 · Post-cluster analysis (domains/boundaries)",
             "command": command("-m", "exondomaincompare.framework.run_core_gene_analysis",
                                "--post", "--run-id", run_id),
             "explanation": "Parse real InterProScan/pyTMHMM outputs, rebuild domain architecture, "
                            "compute all-exon boundary analysis and regenerate indices.",
             "cluster": False})
    else:
        commands.append(
            {"id": "post_interpro", "title": "5 · Post-InterPro analysis",
             "command": command("scripts/run_post_interpro_for_run.py",
                                "--run-id", run_id),
             "explanation": "Run local domain-architecture, boundary consistency, audit and website indices.",
             "cluster": False})

    roundtrip_script = SCRIPTS_DIR / "interpro_cluster" / "run_cluster_roundtrip.py"
    roundtrip_command = _cluster_roundtrip_command(run_id)
    return {
        "run_id": run_id,
        "project_root": ".",
        "run_dir": _public_run_path(run_dir),
        "readme_path": _public_run_path(run_dir / "00_README_NEXT_STEPS.md"),
        "configuration_profile": RUNTIME_CONFIG.public_identity(),
        "workflow": wf["workflow"],
        "event_layer": wf["event_layer"],
        "commands": commands,
        "cluster_roundtrip": {
            "available": roundtrip_script.exists(),
            "title": "Cluster annotation (one command)",
            "command": roundtrip_command,
            "portable_command": roundtrip_command,
            "explanation": "Runs submit → poll → fetch → post-InterPro in one local terminal command. "
                           "Login / 2FA happen locally; never in the webapp.",
        },
    }


def _local_run_file_checks(run_dir: Path) -> Dict[str, Any]:
    closure = run_dir / "results" / "13_final_pre_interpro_closure"
    primary = closure / "freeze" / "final_pre_interpro_proteins_primary.faa"
    ips_out = run_dir / "results" / "14_interproscan" / "primary" / "output"
    tm_out = (run_dir / "results" / "15_exon_domain_boundary_post_interpro"
              / "pytmhmm_primary" / "output")
    wi = run_dir / "website_indices"

    def any_file(folder: Path, suffix: Optional[str] = None) -> bool:
        if not folder.is_dir():
            return False
        for p in folder.rglob("*"):
            if p.is_file() and (suffix is None or p.name.endswith(suffix)):
                return True
        return False

    return {
        "primary_fasta": primary.exists(),
        "truth_table": (closure / "final_pre_interpro_truth_table.tsv").exists(),
        "interproscan_output": any_file(ips_out, ".tsv"),
        "pytmhmm_output": any_file(tm_out),
        "website_indices": (wi / "run_index.json").exists(),
        "boundary_index": (wi / "boundary_consistency_index.json").exists(),
        "explorable": (wi / "run_index.json").exists()
                      or (closure / "final_pre_interpro_truth_table.tsv").exists(),
        "closure_dir": _public_run_path(closure),
    }


def _req_gene_symbol(req: CreateLocalRunRequest) -> str:
    if (req.gene_symbol or "").strip():
        return req.gene_symbol.strip().upper()
    cs = (req.case_study or "").strip().upper()
    if cs.startswith("FGFR2") or not cs:
        return "FGFR2"
    # case_study like "FGFR1_core_only_pilot" → gene "FGFR1"
    return cs.split("_", 1)[0]


def _species_from_req(req: CreateLocalRunRequest) -> List[str]:
    inline = (req.species_file_content or "").strip() or (req.species_text or "").strip()
    ns = cnr.argparse.Namespace(preset=None, species_list=None, species=None)
    if req.preset in ("full30", "pilot"):
        ns.preset = req.preset
    elif inline:
        ns.species = inline
    else:
        raise HTTPException(status_code=400,
                            detail="Provide a preset, species text, or an uploaded species list.")
    try:
        species, _ = cnr.resolve_species(ns)
    except SystemExit as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    if not species:
        raise HTTPException(status_code=400, detail="The resolved species list is empty.")
    return species


# A syntactically plausible gene symbol: starts with a letter, then letters /
# digits / . _ - (covers FGFR2, TP53, FOXP1, TPM1, C1orf112, HLA-A, …). This is
# a syntactic gate only — biological existence is verified per-species by the
# core runner against the genome annotation (honest data-source errors).
_GENE_SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._-]{0,63}$")


def _is_plausible_gene_symbol(symbol: str) -> bool:
    return bool(_GENE_SYMBOL_RE.match((symbol or "").strip()))


def _launch_shared_run(gene_symbol: str, species: List[str], run_name: str,
                       mode: str = "auto", reuse_run_id: str = "") -> Dict[str, Any]:
    runner = PROJECT_ROOT / "src" / "exondomaincompare" / "framework" / \
        "run_core_gene_analysis.py"
    if not runner.exists():
        raise HTTPException(status_code=500, detail="Core gene analysis module not found.")
    # No pre-existing gene YAML is required: the core runner generates a generic
    # core-only config automatically for any gene. We only gate on a syntactically
    # plausible symbol here; whether the gene actually exists in a given species is
    # verified by the runner against the genome annotation (real biological errors).
    if not _is_plausible_gene_symbol(gene_symbol):
        raise HTTPException(
            status_code=400,
            detail=(f"'{gene_symbol}' is not a valid gene symbol. Enter a symbol such as "
                    "FOXP1, TP53, TPM1 or FGFR1."))
    if not species:
        raise HTTPException(status_code=400,
                            detail="Add at least one species (one scientific name per line).")

    # Pass ALL species to the generic runner (species count is a dataset property,
    # not a separate code path). argparse consumes every value after --species up
    # to the next flag, so one/two/many species use the exact same launcher.
    cmd = _local_python_module_command(
        "exondomaincompare.framework.run_core_gene_analysis",
        "--gene", gene_symbol, "--species", *species,
        "--input-mode", "auto")
    label = run_labels.clean_run_name(run_name)
    if label:
        cmd += ["--run-name", label]
    # Retrying a failed run repairs that run. Creating a second one for the same request
    # leaves the user to work out which of the two is the answer.
    if reuse_run_id:
        cmd += ["--reuse-run-id", reuse_run_id]

    LOCAL_RUNS_ROOT.mkdir(parents=True, exist_ok=True)
    before = {p.name for p in LOCAL_RUNS_ROOT.iterdir() if p.is_dir()}
    launch_log = LOCAL_RUNS_ROOT / f".shared_launch_{int(time.time())}.log"
    lf = launch_log.open("w", encoding="utf-8")
    try:
        proc = subprocess.Popen(cmd, cwd=str(PROJECT_ROOT), stdout=lf,
                                stderr=subprocess.STDOUT, start_new_session=True,
                                env=_runtime_env())
    finally:
        lf.close()

    # The runner creates runs/<id>/ early in its create phase; poll briefly.
    new_id: Optional[str] = None
    for _ in range(80):  # ~8s
        now = {p.name for p in LOCAL_RUNS_ROOT.iterdir() if p.is_dir()}
        diff = sorted(now - before)
        if diff:
            new_id = diff[-1]
            break
        if proc.poll() is not None and proc.returncode != 0:
            break
        time.sleep(0.1)

    base = {"command": " ".join(cmd), "pid": proc.pid,
            "workflow": "shared_exploratory", "started": True, "gene_symbol": gene_symbol}
    if reuse_run_id:
        # No new folder appears, so the "did it start" signal is the existing run itself.
        _LOCAL_PROCS[reuse_run_id] = proc
        try:
            (LOCAL_RUNS_ROOT / reuse_run_id / "logs").mkdir(parents=True, exist_ok=True)
            launch_log.replace(LOCAL_RUNS_ROOT / reuse_run_id / "logs" / "shared_launch.log")
        except Exception:
            pass
        return {**base, "run_id": reuse_run_id, "in_place": True}
    if new_id:
        _LOCAL_PROCS[new_id] = proc
        try:
            (LOCAL_RUNS_ROOT / new_id / "logs").mkdir(parents=True, exist_ok=True)
            launch_log.replace(LOCAL_RUNS_ROOT / new_id / "logs" / "shared_launch.log")
        except Exception:
            pass
        return {**base, "run_id": new_id}

    log_tail = ""
    try:
        log_tail = launch_log.read_text(encoding="utf-8")[-2000:]
    except Exception:
        pass
    return {**base, "run_id": None, "provisioning": True, "log_tail": log_tail,
            "note": ("The shared exploratory run is provisioning. A run configuration is "
                     "generated automatically — no gene YAML is needed. If it does not appear "
                     "in the run list, the gene or a species could not be retrieved; see the "
                     "run log for the specific data-source reason.")}


def _do_create_local_run(req: CreateLocalRunRequest) -> Path:
    if "FGFR2" not in (req.case_study or "").upper().replace(" ", ""):
        raise HTTPException(status_code=400,
                            detail="This workflow supports the FGFR2 IIIb/IIIc case study.")

    inline = (req.species_file_content or "").strip() or (req.species_text or "").strip()
    ns = cnr.argparse.Namespace(preset=None, species_list=None, species=None)
    if req.preset in ("full30", "pilot"):
        ns.preset = req.preset
    elif inline:
        ns.species = inline
    else:
        raise HTTPException(status_code=400,
                            detail="Provide a preset, species text, or an uploaded species list.")

    try:
        species, source_label = cnr.resolve_species(ns)
    except SystemExit as exc:  # cnr raises SystemExit on bad input
        raise HTTPException(status_code=400, detail=str(exc))
    if not species:
        raise HTTPException(status_code=400, detail="The resolved species list is empty.")

    LOCAL_RUNS_ROOT.mkdir(parents=True, exist_ok=True)
    run_name = run_labels.clean_run_name(req.run_name)
    slug = run_labels.run_id_slug(run_name, gene_symbol=_req_gene_symbol(req),
                                  species=species)
    run_id, run_dir = cnr.unique_run_dir(cnr.generate_run_id(slug))

    species_path = run_dir / "config" / "species.tsv"
    record = cnr.build_run_config(
        run_id, run_name, run_dir, species, species_path, source_label)
    gene = {
        "gene_symbol": cnr.GENE_SYMBOL,
        "analysis_id": cnr.ANALYSIS_ID,
        "event_id": cnr.EVENT_ID,
        "event_type": cnr.EVENT_TYPE,
        "source": f"repo:{cnr.GENE_CONFIG_REL}",
    }
    RunLayout(run_dir, RunLayoutVersion.CANONICAL_V2).initialize(
        run_record=record,
        status=cnr.build_status(run_id, species, run_dir),
        gene=gene,
        species=species,
    )
    cnr.validate_run(run_dir, species)
    return run_dir


def _launch_pre_interpro(run_dir: Path, mode: str = "cached") -> Dict[str, Any]:
    run_id = run_dir.name
    _require_writable_run(run_id)
    wrapper = SCRIPTS_DIR / "run_pre_interpro_for_run.py"
    if not wrapper.exists():
        raise HTTPException(status_code=500, detail="run_pre_interpro_for_run.py not found.")

    logs_dir = run_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    launch_log = logs_dir / "pre_interpro_launch.log"
    pipeline_log = logs_dir / "pre_interpro_pipeline.log"
    pipeline_err = logs_dir / "pre_interpro_pipeline.err"

    cmd = _local_python_command(wrapper, "--run-id", run_id, "--env", mode)
    if mode == "live":
        cmd.append("--force")

    lf = launch_log.open("w", encoding="utf-8")
    try:
        proc = subprocess.Popen(cmd, cwd=str(PROJECT_ROOT), stdout=lf,
                                stderr=subprocess.STDOUT, start_new_session=True,
                                env=_runtime_env())
    finally:
        lf.close()
    _LOCAL_PROCS[run_id] = proc

    meta = {
        "pid": proc.pid,
        "command": " ".join(cmd),
        "start_time": now_iso(),
        "mode": mode,
        "log_path": _public_run_path(pipeline_log),
        "err_path": _public_run_path(pipeline_err),
        "launch_log": _public_run_path(launch_log),
    }
    _update_run_status(run_dir, pre_interpro_status="running",
                       current_step="pre_interpro_analysis", status="running",
                       pre_interpro_process=meta, error="")
    return meta


@app.get("/api/analysis-router")
def analysis_router(gene: str = Query(default="FGFR2"),
                    mode: str = Query(default="auto")) -> Dict[str, Any]:
    wf = _resolve_workflow(gene, mode)
    wf["validated_genes"] = (_router.list_validated_genes()
                             if _router is not None else ["FGFR2"])
    return wf


def _resolve_species_preview(species_text: str) -> Dict[str, Any]:
    resolved: List[Dict[str, Any]] = []
    invalid: List[Dict[str, str]] = []
    seen: set = set()
    duplicates = 0
    raw_text = species_text or ""
    for raw_line in re.split(r"[\n;,]", raw_text):
        s = raw_line.strip()
        if not s or s.startswith("#"):
            continue
        token = cnr.normalize_species_token(s)
        if not token or not cnr.SPECIES_ID_RE.match(token):
            invalid.append({"raw": s, "message": cnr.species_error_message(s)})
            continue
        if token in seen:
            duplicates += 1
            continue
        seen.add(token)
        entry = None
        if _lookup_known_species is not None:
            try:
                _sci, entry = _lookup_known_species(s)
            except Exception:
                entry = None
        if entry is not None:
            resolved.append({
                "input": s,
                "species_id": entry["ensembl_species"],
                "scientific_name": entry["ncbi_species"],
                "taxid": entry["taxid"],
                "common_name": entry.get("common_name", ""),
                "known": True,
            })
        else:
            resolved.append({
                "input": s, "species_id": token,
                "scientific_name": s, "taxid": "", "common_name": "",
                "known": False,
            })
    return {"count": len(resolved), "duplicates_removed": duplicates,
            "resolved": resolved, "invalid": invalid}


class ResolveInputsRequest(BaseModel):
    gene_symbol: Optional[str] = Field(default=None)
    species_text: Optional[str] = Field(default=None)
    species_file_content: Optional[str] = Field(default=None)
    mode: str = Field(default="auto")


@app.post("/api/runs/resolve-inputs")
def resolve_run_inputs(req: ResolveInputsRequest) -> Dict[str, Any]:
    sym = (req.gene_symbol or "").strip().upper()
    gene_valid = _is_plausible_gene_symbol(sym)
    gene_msg = "" if gene_valid or not sym else (
        f"'{req.gene_symbol}' is not a valid gene symbol. Enter a symbol such as "
        "FOXP1, TP53, TPM1 or FGFR1.")

    wf = _resolve_workflow(sym, req.mode) if sym else None
    text = (req.species_file_content or "").strip() or (req.species_text or "")
    species = _resolve_species_preview(text)

    note = None
    if wf is not None and not wf.get("is_validated"):
        note = ("Generic exploratory workflow. No gene-specific configuration is required — "
                "a run configuration will be generated automatically.")
    elif wf is not None and wf.get("is_validated"):
        note = ("Validated FGFR2 IIIb/IIIc workflow (events, cassette, human comparison). "
                "Uses its frozen, hand-authored specialization.")

    return {
        "gene": {"symbol": sym, "valid": gene_valid, "message": gene_msg},
        "workflow": wf,
        "workflow_note": note,
        "species": species,
    }


@app.post("/api/local-runs/create")
def create_local_run(req: CreateLocalRunRequest) -> Dict[str, Any]:
    wf = _resolve_workflow(_req_gene_symbol(req), req.mode)
    if not wf["is_validated"]:
        launched = _launch_shared_run(wf["gene_symbol"], _species_from_req(req),
                                      req.run_name, req.mode)
        return {**launched, "workflow": wf}
    run_dir = _do_create_local_run(req)
    run_id = run_dir.name
    return {
        "run_id": run_id,
        "run_dir": _public_run_path(run_dir),
        "species_count": _local_run_summary(run_dir)["species_count"],
        "summary": _local_run_summary(run_dir),
        "next_steps": _local_run_commands(run_id, run_dir),
        "workflow": wf,
    }


@app.post("/api/local-runs/start")
def create_and_start_local_run(req: CreateLocalRunRequest) -> Dict[str, Any]:
    other = _any_local_run_running()
    if other:
        raise HTTPException(status_code=409,
                            detail=f"Another local run is already running ({other}). "
                                   "Wait for it to finish before starting a new one.")

    wf = _resolve_workflow(_req_gene_symbol(req), req.mode)
    if not wf["is_validated"]:
        launched = _launch_shared_run(wf["gene_symbol"], _species_from_req(req),
                                      req.run_name, req.mode)
        run_id = launched.get("run_id")
        run_dir = _safe_local_run_dir(run_id) if run_id else None
        if run_dir and run_dir.is_dir():
            _ensure_run_indices(run_dir)
            launched["summary"] = _local_run_summary(run_dir)
            launched["status_model"] = derive_status_model(run_dir)
            launched["next_steps"] = _local_run_commands(run_id, run_dir)
        return {**launched, "workflow": wf}

    run_dir = _do_create_local_run(req)
    run_id = run_dir.name
    meta = _launch_pre_interpro(run_dir, mode="cached")
    _ensure_run_indices(run_dir)
    return {
        "run_id": run_id,
        "run_dir": _public_run_path(run_dir),
        "started": True,
        "process": meta,
        "summary": _local_run_summary(run_dir),
        "status_model": derive_status_model(run_dir),
        "next_steps": _local_run_commands(run_id, run_dir),
        "workflow": wf,
    }


@app.get("/api/local-runs/{run_id}/commands")
def local_run_commands(run_id: str) -> Dict[str, Any]:
    run_dir = _safe_local_run_dir(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")
    return _local_run_commands(run_id, run_dir)


@app.get("/api/local-runs/{run_id}/status")
def local_run_status(run_id: str) -> Dict[str, Any]:
    run_dir = _safe_local_run_dir(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")
    return read_json(run_dir / "status.json", {}) or {}


@app.post("/api/local-runs/{run_id}/refresh")
def local_run_refresh(run_id: str) -> Dict[str, Any]:
    run_dir = _safe_local_run_dir(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")
    # Build partial run-local indices as soon as pre-InterPro produced a closure,
    # so overview / gene explorer become explorable before cluster annotation.
    _ensure_run_indices(run_dir)
    return {
        "run_id": run_id,
        "run_dir": _public_run_path(run_dir),
        "config": read_json(run_dir / "run_config.json", {}) or {},
        "status": read_json(run_dir / "status.json", {}) or {},
        "status_model": derive_status_model(run_dir),
        "files": _local_run_file_checks(run_dir),
        "summary": _local_run_summary(run_dir),
    }


@app.get("/api/web-runs/legacy")
def legacy_web_runs() -> Dict[str, Any]:
    items: List[Dict[str, Any]] = []
    if WEB_RUNS_ROOT.exists():
        for p in sorted(WEB_RUNS_ROOT.glob("run_*"), reverse=True):
            st = read_json(p / "web_status.json", {}) or {}
            items.append({
                "run_id": p.name,
                "path": str(p.relative_to(PROJECT_ROOT)),
                "status": st.get("status", "unknown"),
                "created_at": st.get("created_at", ""),
            })
    return {
        "count": len(items),
        "runs": items,
        "note": "Legacy web-runs are not part of the new runs/<run_id>/ system.",
    }


# Local preparation never executes LRZ, SSH, SLURM or cluster annotation.
class StartPreInterproRequest(BaseModel):
    mode: str = Field(default="cached")     # cached (local-safe) | live (full refresh)
    confirm_live_run: bool = False
    allow_concurrent: bool = False


# In-process registry of launched local wrappers so we can reap them (poll()
# reaps a finished child, avoiding zombie PIDs that would look "alive").
_LOCAL_PROCS: Dict[str, "subprocess.Popen"] = {}


def _pid_alive(pid: Optional[int]) -> bool:
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def _proc_running(run_id: str, pid: Optional[int]) -> bool:
    proc = _LOCAL_PROCS.get(run_id)
    if proc is not None:
        if proc.poll() is None:
            return True
        _LOCAL_PROCS.pop(run_id, None)
        return False
    return _pid_alive(pid)


def append_local_run_log(run_dir: Path, message: str) -> None:
    try:
        logs_dir = run_dir / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        with (logs_dir / "web_actions.log").open("a", encoding="utf-8") as fh:
            fh.write(f"{now_iso()}  {message.rstrip()}\n")
    except Exception:
        pass


def _update_run_status(run_dir: Path, **fields: Any) -> Dict[str, Any]:
    p = run_dir / "status.json"
    st = read_json(p, {}) or {}
    st.update(fields)
    st["last_updated"] = now_iso()
    write_json(p, st)
    return st


def _preinterpro_running(run_dir: Path) -> Optional[Dict[str, Any]]:
    st = read_json(run_dir / "status.json", {}) or {}
    proc = st.get("pre_interpro_process") or {}
    return proc if _proc_running(run_dir.name, proc.get("pid")) else None


def _postinterpro_running(run_dir: Path) -> bool:
    run_id = run_dir.name
    proc = _LOCAL_POST_PROCS.get(run_id)
    if proc is not None:
        if proc.poll() is None:
            return True
        _LOCAL_POST_PROCS.pop(run_id, None)
    st = read_json(run_dir / "status.json", {}) or {}
    meta = st.get("post_interpro_process") or {}
    return _pid_alive(meta.get("pid"))


def _any_local_run_running(exclude: Optional[str] = None) -> Optional[str]:
    if not LOCAL_RUNS_ROOT.exists():
        return None
    for d in LOCAL_RUNS_ROOT.iterdir():
        if not d.is_dir() or d.name == exclude \
                or d.name.startswith(".") or d.name.startswith("_"):
            continue
        if _preinterpro_running(d):
            return d.name
    return None


@app.post("/api/local-runs/{run_id}/start-preinterpro")
def start_preinterpro(run_id: str, req: StartPreInterproRequest) -> Dict[str, Any]:
    run_dir = _require_writable_run(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")

    if _preinterpro_running(run_dir):
        raise HTTPException(status_code=409, detail="Pre-InterPro is already running for this run.")
    if not req.allow_concurrent:
        other = _any_local_run_running(exclude=run_id)
        if other:
            raise HTTPException(status_code=409,
                                detail=f"Another local run is already running ({other}). "
                                       "Pass allow_concurrent=true to override.")

    mode = (req.mode or "cached").lower()
    if mode not in ("cached", "live"):
        raise HTTPException(status_code=400, detail="mode must be 'cached' or 'live'.")
    if mode == "live" and not req.confirm_live_run:
        raise HTTPException(status_code=400,
                            detail="Live / full-refresh runs require confirm_live_run=true.")

    meta = _launch_pre_interpro(run_dir, mode=mode)
    return {"run_id": run_id, "status": "started", "process": meta}


@app.post("/api/local-runs/{run_id}/stop-preinterpro")
def stop_preinterpro(run_id: str) -> Dict[str, Any]:
    run_dir = _require_writable_run(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")
    st = read_json(run_dir / "status.json", {}) or {}
    proc = dict(st.get("pre_interpro_process") or {})
    pid = proc.get("pid")
    proc_obj = _LOCAL_PROCS.get(run_id)
    terminated = False
    if _proc_running(run_id, pid):
        def _sig(sig):
            try:
                os.killpg(os.getpgid(int(pid)), sig)
            except (OSError, ProcessLookupError):
                try:
                    os.kill(int(pid), sig)
                except OSError:
                    pass
        _sig(signal.SIGTERM)
        for _ in range(20):
            if not _proc_running(run_id, pid):
                break
            time.sleep(0.1)
        if _proc_running(run_id, pid):
            _sig(signal.SIGKILL)
        terminated = True
    # reap the child so it does not linger as a zombie
    if proc_obj is not None:
        try:
            proc_obj.wait(timeout=3)
        except Exception:
            pass
        _LOCAL_PROCS.pop(run_id, None)
    proc["stopped_at"] = now_iso()
    _update_run_status(run_dir, pre_interpro_status="stopped",
                       current_step="pre_interpro_stopped", pre_interpro_process=proc)
    return {"run_id": run_id, "stopped": terminated}


def _terminate_pid(pid: Optional[int]) -> bool:
    if not _pid_alive(pid):
        return False

    def _sig(sig):
        try:
            os.killpg(os.getpgid(int(pid)), sig)
        except (OSError, ProcessLookupError):
            try:
                os.kill(int(pid), sig)
            except OSError:
                pass

    _sig(signal.SIGTERM)
    for _ in range(20):
        if not _pid_alive(pid):
            break
        time.sleep(0.1)
    if _pid_alive(pid):
        _sig(signal.SIGKILL)
    return True


@app.post("/api/local-runs/{run_id}/stop")
def stop_local_run(run_id: str) -> Dict[str, Any]:
    run_dir = _require_writable_run(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")

    st = read_json(run_dir / "status.json", {}) or {}
    stopped_any = False
    stopped_steps: List[str] = []

    # Local pre-InterPro subprocess
    pre_meta = dict(st.get("pre_interpro_process") or {})
    pre_pid = pre_meta.get("pid")
    if _proc_running(run_id, pre_pid):
        _terminate_pid(pre_pid)
        obj = _LOCAL_PROCS.pop(run_id, None)
        if obj is not None:
            try:
                obj.wait(timeout=3)
            except Exception:
                pass
        pre_meta["stopped_at"] = now_iso()
        stopped_any = True
        stopped_steps.append("pre_interpro")

    # Local post-InterPro subprocess
    post_meta = dict(st.get("post_interpro_process") or {})
    post_pid = post_meta.get("pid")
    post_obj = _LOCAL_POST_PROCS.get(run_id)
    post_alive = (post_obj is not None and post_obj.poll() is None) or _pid_alive(post_pid)
    if post_alive:
        _terminate_pid(post_pid)
        if post_obj is not None:
            try:
                post_obj.wait(timeout=3)
            except Exception:
                pass
        _LOCAL_POST_PROCS.pop(run_id, None)
        post_meta["stopped_at"] = now_iso()
        stopped_any = True
        stopped_steps.append("post_interpro")

    # Only mutate the run status when we actually terminated a local process. If
    # nothing local was running (already finished, or only cluster jobs), leave
    # the status untouched and just report the cluster note.
    if stopped_any:
        updates: Dict[str, Any] = {
            "status": "stopped",
            "last_error": "stopped by user",
            "next_action": "restart_or_delete",
            "current_step": "stopped_by_user",
        }
        if "pre_interpro" in stopped_steps:
            updates["pre_interpro_status"] = "stopped"
            updates["pre_interpro_process"] = pre_meta
        if "post_interpro" in stopped_steps:
            updates["post_interpro_status"] = "stopped"
            updates["post_interpro_process"] = post_meta
        _update_run_status(run_dir, **updates)
        append_local_run_log(run_dir, "Run stopped by user via webapp "
                                       f"(steps: {', '.join(stopped_steps)}).")

    return {
        "run_id": run_id,
        "stopped": stopped_any,
        "stopped_steps": stopped_steps,
        "cluster_note": ("Cluster jobs cannot be stopped from the webapp. "
                         "Use your terminal / SLURM if cluster jobs are still running."),
        "status_model": derive_status_model(run_dir),
    }


@app.post("/api/local-runs/{run_id}/refresh-all")
def local_run_refresh_all(run_id: str) -> Dict[str, Any]:
    run_dir = _safe_local_run_dir(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")
    # Force a rebuild so freshly-fetched post-InterPro outputs (e.g. boundary
    # consistency) become visible without a browser reload.
    _ensure_run_indices(run_dir, rebuild=True)
    return {
        "run_id": run_id,
        "run_dir": _public_run_path(run_dir),
        "config": read_json(run_dir / "run_config.json", {}) or {},
        "status": read_json(run_dir / "status.json", {}) or {},
        "status_model": derive_status_model(run_dir),
        "files": _local_run_file_checks(run_dir),
        "summary": _local_run_summary(run_dir),
    }


@app.delete("/api/local-runs/{run_id}")
def delete_local_run(run_id: str) -> Dict[str, Any]:
    if not run_id or run_id in ("example", "current") or "/" in run_id or "\\" in run_id \
            or run_id.startswith(".") or run_id.startswith("_"):
        raise HTTPException(status_code=400, detail="Invalid run id.")

    try:
        record = resolve_run_record(RUNTIME_CONFIG, run_id)
    except RegistryError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    if record is None:
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")
    run_dir = record.path
    if record.kind in {"legacy_external", "repository_legacy"}:
        result = hide_discovered_run(
            RUNTIME_CONFIG, run_id=run_id, root_id=record.root_id)
        return {
            "run_id": run_id,
            "deleted": False,
            "unregistered": True,
            "source_unchanged": True,
            "registry_records_removed": result,
        }
    runs_root = LOCAL_RUNS_ROOT.resolve()
    if run_dir.parent != runs_root or run_dir == runs_root:
        raise HTTPException(status_code=400, detail="Canonical run path is outside configured storage.")
    if str(run_dir) == str(EXAMPLE_RUN_DIR.resolve()) \
            or str(run_dir).startswith(str(RESULTS_ROOT.resolve())):
        raise HTTPException(status_code=400, detail="Refusing to delete the example/freeze dataset.")
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")

    # Stop any local subprocess we started for this run (best-effort).
    try:
        stop_local_run(run_id)
    except HTTPException:
        pass
    except Exception:
        pass

    trash_root = RUNTIME_CONFIG.paths.deleted
    trash_root.mkdir(parents=True, exist_ok=True)
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    dest = trash_root / f"{run_id}_{ts}"
    try:
        shutil.move(str(run_dir), str(dest))
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Could not delete run: {exc}")
    try:
        unregister(RUNTIME_CONFIG, run_id=run_id)
    except RegistryError:
        pass

    return {
        "run_id": run_id,
        "deleted": True,
        "archived_to": _public_run_path(dest),
    }


@app.get("/api/local-runs/{run_id}/logs/preinterpro")
def local_run_preinterpro_logs(run_id: str,
                               tail: int = Query(200, ge=1, le=5000)) -> Dict[str, Any]:
    run_dir = _safe_local_run_dir(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")

    def _tail(p: Path) -> List[str]:
        if not p.exists():
            return []
        return p.read_text(encoding="utf-8", errors="replace").splitlines()[-tail:]

    logs_dir = run_dir / "logs"
    st = read_json(run_dir / "status.json", {}) or {}
    proc = st.get("pre_interpro_process") or {}
    return {
        "run_id": run_id,
        "running": _proc_running(run_id, proc.get("pid")),
        "log": _tail(logs_dir / "pre_interpro_pipeline.log"),
        "err": _tail(logs_dir / "pre_interpro_pipeline.err"),
        "launch_log": _tail(logs_dir / "pre_interpro_launch.log"),
        "status": st,
    }


@app.get("/api/local-runs/{run_id}/logs/core")
def local_run_core_logs(run_id: str,
                        tail: int = Query(400, ge=1, le=8000)) -> Dict[str, Any]:
    run_dir = _safe_local_run_dir(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")

    def _tail(p: Path) -> List[str]:
        if not p.exists():
            return []
        return p.read_text(encoding="utf-8", errors="replace").splitlines()[-tail:]

    logs_dir = run_dir / "logs"
    st = read_json(run_dir / "status.json", {}) or {}
    download_logs: Dict[str, List[str]] = {}
    if logs_dir.is_dir():
        for lp in sorted(logs_dir.glob("datasets_download_*.log")):
            download_logs[lp.name] = _tail(lp)
    return {
        "run_id": run_id,
        "status": st.get("status", ""),
        "failed_stage": st.get("failed_step", ""),
        "failed_reason": st.get("failed_reason", "") or st.get("error", ""),
        "failed_species": st.get("failed_species", ""),
        "detail": st.get("detail", ""),
        "core_runner": _tail(logs_dir / "core_runner.log"),
        "shared_launch": _tail(logs_dir / "shared_launch.log"),
        "downloads": download_logs,
    }


@app.get("/api/local-runs/{run_id}/diagnostics")
def local_run_diagnostics(run_id: str) -> FileResponse:
    run_dir = _safe_local_run_dir(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")
    out_dir = WEB_STATE_DIR / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    archive = out_dir / f"{run_id}_diagnostics.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
        for name in ("run_config.json", "status.json"):
            path = run_dir / name
            if path.is_file():
                zf.write(path, f"{run_id}/{name}")
        logs_dir = run_dir / "logs"
        if logs_dir.is_dir():
            for path in sorted(logs_dir.rglob("*")):
                if path.is_file():
                    zf.write(path, f"{run_id}/logs/{path.relative_to(logs_dir)}")
    return FileResponse(archive, media_type="application/zip",
                        filename=archive.name)


@app.post("/api/local-runs/{run_id}/retry-local-preparation")
def retry_local_preparation(run_id: str) -> Dict[str, Any]:
    run_dir = _require_writable_run(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")
    if _preinterpro_running(run_dir) is not None:
        raise HTTPException(status_code=409,
                            detail="Local preparation is already running for this run.")

    model = derive_status_model(run_dir)
    if model.get("primary_fasta_status") == "available":
        raise HTTPException(
            status_code=409,
            detail="This run already has its primary FASTA; there is no failed local "
                   "stage to rebuild.")

    # `live` rather than `cached`: the point of the retry is to redo the stages that
    # produced nothing, and a cached retry would reuse exactly those empty outputs.
    meta = _launch_pre_interpro(run_dir, mode="live")
    return {
        "run_id": run_id,
        "resumed_in_place": True,
        "launch": meta,
        "note": ("Local preparation restarted for this run. The run id, run name and "
                 "requested species are unchanged."),
    }


@app.post("/api/local-runs/{run_id}/retry-precluster")
def retry_precluster(run_id: str) -> Dict[str, Any]:
    run_dir = _require_writable_run(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")

    rc = read_json(run_dir / "run_config.json", {}) or {}
    gene_symbol = (rc.get("gene_symbol") or "").strip()
    if not gene_symbol:
        raise HTTPException(status_code=400,
                            detail="This run has no gene symbol to retry.")
    wf = _resolve_workflow(gene_symbol)
    if wf["is_validated"]:
        raise HTTPException(status_code=400,
                            detail="Retry applies to shared exploratory runs only; the "
                                   "validated FGFR2 workflow is unchanged.")

    # Only retry runs that actually stopped in a failed/incomplete pre-cluster state.
    model = derive_status_model(run_dir)
    if model.get("status") not in ("core_model_collection_failed", "incomplete", "failed"):
        raise HTTPException(status_code=409,
                            detail=f"Run is not in a failed pre-cluster state "
                                   f"(status: {model.get('status')}). Nothing to retry.")

    # Species from the authoritative run-local list (falls back to config).
    species: List[str] = LegacyRunAdapter(run_dir).species()
    if not species:
        species = [s for s in (rc.get("species_ids") or []) if s]
    if not species:
        raise HTTPException(status_code=400, detail="No species recorded on this run to retry.")

    launched = _launch_shared_run(gene_symbol, species,
                                  rc.get("run_name") or "", mode="auto",
                                  reuse_run_id=run_id)
    return {"retried_from": run_id, "run_id": run_id, "in_place": True,
            "gene_symbol": gene_symbol, "species": species, "launch": launched}


_LOCAL_POST_PROCS: Dict[str, "subprocess.Popen"] = {}


@app.post("/api/local-runs/{run_id}/start-post-interpro")
def start_post_interpro(run_id: str) -> Dict[str, Any]:
    run_dir = _require_writable_run(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Local run '{run_id}' not found.")

    existing = _LOCAL_POST_PROCS.get(run_id)
    if existing is not None and existing.poll() is None:
        raise HTTPException(status_code=409, detail="Post-InterPro is already running for this run.")

    files = _local_run_file_checks(run_dir)
    if not (files["interproscan_output"] and files["pytmhmm_output"]):
        raise HTTPException(status_code=409,
                            detail="Cluster outputs not present yet. Run the cluster round-trip first.")

    wrapper = SCRIPTS_DIR / "run_post_interpro_for_run.py"
    if not wrapper.exists():
        raise HTTPException(status_code=500, detail="run_post_interpro_for_run.py not found.")

    logs_dir = run_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    launch_log = logs_dir / "post_interpro_launch.log"
    cmd = _local_python_command(wrapper, "--run-id", run_id)
    lf = launch_log.open("w", encoding="utf-8")
    try:
        proc = subprocess.Popen(cmd, cwd=str(PROJECT_ROOT), stdout=lf,
                                stderr=subprocess.STDOUT, start_new_session=True,
                                env=_runtime_env())
    finally:
        lf.close()
    _LOCAL_POST_PROCS[run_id] = proc
    proc_meta = {"pid": proc.pid, "command": " ".join(cmd),
                 "start_time": now_iso(),
                 "launch_log": _public_run_path(launch_log)}
    _update_run_status(run_dir, post_interpro_status="running",
                       current_step="post_interpro_analysis", status="running", error="",
                       post_interpro_process=proc_meta)
    return {"run_id": run_id, "status": "started", "process": proc_meta}


# --------------------------------------------------------------------------- #
# Datasets — the single source of truth for the UI's dataset switcher.
# --------------------------------------------------------------------------- #
@app.get("/api/datasets")
def list_datasets() -> Dict[str, Any]:
    datasets: List[Dict[str, Any]] = []
    if EXAMPLE_RUN_DIR.exists():
        datasets.append({
            "id": "example",
            "kind": "example",
            "label": "Example 30-species dataset",
            "switcher_label": "Example · FGFR2 IIIb/IIIc (30 species)",
            "sublabel": "Validated 30-species thesis dataset",
            "gene_symbol": "FGFR2",
            "pipeline_type": "validated_event_pipeline",
            "support_level": "validated_event_analysis",
            "status": "results_ready",
            "status_label": "Results ready",
            "read_only": True,
            "explorable": True,
            "available_views": {
                "overview": True, "gene_explorer": True, "domain_architecture": True,
                "boundary_consistency": True, "figure_gallery": True,
            },
        })
    runs_out = []
    for s in _local_run_summaries(include_bundled=True):
            if s.get("collision"):
                continue
            core_only = bool(s.get("core_only"))
            gene = s.get("gene_symbol") or ""
            # A clear, consistent switcher label so the dropdown always matches the
            # displayed content. The visible run name wins when the user gave one;
            # otherwise the label describes the biology rather than the directory.
            label = s.get("display_name") or s["run_id"]
            if s.get("run_name"):
                sw = f"{gene} · {s['run_name']}" if gene else s["run_name"]
            elif core_only:
                nsp = s.get("species_count") or 0
                scope = "single-species" if nsp <= 1 else f"{nsp}-species"
                sw = f"{gene or s['run_id']} · {scope} gene analysis"
            else:
                sw = label
            runs_out.append({
                "id": f"run:{s['run_id']}",
                "kind": "run",
                "run_id": s["run_id"],
                "label": label,
                "display_name": label,
                "switcher_label": sw,
                "sublabel": s.get("status_label", ""),
                "gene_symbol": gene,
                "pipeline_type": s.get("pipeline_type", "validated_event_pipeline"),
                "support_level": s.get("support_level", ""),
                "core_only": core_only,
                "experimental": bool(s.get("experimental")),
                "status": s.get("status", "created"),
                "status_label": s.get("status_label", ""),
                "species_count": s.get("species_count", 0),
                "created_at": s.get("created_at", ""),
                "read_only": bool(s.get("read_only")),
                "bundled_example": s.get("root_kind") == "bundled_example",
                "root_id": s.get("root_id", ""),
                "layout_version": s.get("layout_version", ""),
                "explorable": s.get("explorable", False),
                "available_views": s.get("available_views", {}),
                "human_reference": s.get("human_reference", {}),
            })
    datasets.extend(run_labels.sort_runs(runs_out))
    return {"datasets": datasets, "active_default": "example"}


@app.get("/api/datasets/human-reference")
def dataset_human_reference(rebuild: bool = False) -> Dict[str, Any]:
    return build_human_reference_control(rebuild=rebuild)


@app.get("/api/analysis-capabilities")
def analysis_capabilities() -> Dict[str, Any]:
    if _gene_config is None:
        return {
            "supported": ["FGFR2_IIIb_IIIc"],
            "drafts": [],
            "module_capabilities": {},
            "note": "Only FGFR2 is currently runnable. (Config layer unavailable.)",
        }
    try:
        analyses = _gene_config.discover_analyses()
        caps = _gene_config.capability_summary()
    except Exception:
        return {
            "supported": ["FGFR2_IIIb_IIIc"], "drafts": [],
            "module_capabilities": {},
            "note": "Only FGFR2 is currently runnable.",
        }
    detectors: Dict[str, Any] = {"supported": [], "planned": []}
    try:
        for name, spec in _gene_config.supported_detectors().items():
            detectors["supported"].append({"detector": name,
                                            "analysis_id": spec.get("analysis_id"),
                                            "event_type": spec.get("event_type")})
        for name, spec in _gene_config.planned_detectors().items():
            detectors["planned"].append({"detector": name,
                                         "event_type": spec.get("event_type"),
                                         "status": spec.get("status")})
    except Exception:
        pass
    return {
        "supported": analyses.get("supported", []),
        "core_only_pilots": analyses.get("core_only_pilots", []),
        "drafts": analyses.get("drafts", []),
        "supported_detail": analyses.get("supported_detail", []),
        "core_only_pilots_detail": analyses.get("core_only_pilots_detail", []),
        "drafts_detail": analyses.get("drafts_detail", []),
        "detectors": detectors,
        "support_levels": {
            "validated_event_analysis": "Runnable, validated event-specific analysis (FGFR2 IIIb/IIIc).",
            "core_only_pilot": "Experimental core-only proof of concept; no event region.",
            "draft_not_runnable": "Config exists but no working runner yet.",
        },
        "module_capabilities": {
            "reusable_modules": caps.get("reusable_modules", []),
            "gene_specific_modules": caps.get("gene_specific_modules", []),
            "partial_or_event_specific_modules": caps.get("partial_or_event_specific_modules", []),
        },
        "note": "Support levels: validated_event_analysis (runnable), core_only_pilot "
                "(experimental, core-only), draft_not_runnable (config only). Normal "
                "runnable requires an active config AND a supported event detector.",
    }


@app.get("/api/datasets/{dataset_id}/status")
def dataset_status(dataset_id: str) -> Dict[str, Any]:
    ds = resolve_dataset(dataset_id)
    if ds["kind"] == "example":
        gm = _gene_meta_for_run(None)
        return {
            "id": "example", "kind": "example", "status": "results_ready",
            "status_label": "Results ready", "explorable": True, "read_only": True,
            "analysis_id": gm["analysis_id"], "gene_symbol": gm["gene_symbol"],
            "event_id": gm["event_id"], "event_type": gm["event_type"],
            "event_display_name": gm["event_display_name"], "ui_labels": gm["ui_labels"],
            "analysis_modes": gm.get("analysis_modes",
                                     {"core_gene_analysis": True, "event_analysis": "configured"}),
            "event_analysis_mode": gm.get("event_analysis_mode", "configured"),
            "has_event": bool(gm.get("has_event", True)),
            "core_only": False,
            "pipeline_type": "validated_event_pipeline",
            "event_layer_type": ("validated_fgfr2_iiib_iiic" if bool(gm.get("has_event", True))
                                 else "exploratory_event_evidence"),
            "has_validated_event": bool(gm.get("has_event", True)),
            "event_message": "",
            "available_views": {
                "overview": True, "gene_explorer": True, "domain_architecture": True,
                "exon_domain_boundaries": True, "synteny": True,
                "boundary_consistency": bool(gm.get("has_event", True)),
                "event_region": bool(gm.get("has_event", True)),
                "figure_gallery": True,
            },
        }
    run_dir = ds["run_base"]
    _ensure_run_indices(run_dir)
    model = derive_status_model(run_dir)
    model["index_version"] = index_version(run_dir)
    return {"id": f"run:{ds['run_id']}", "kind": "run", **model,
            "config": read_json(run_dir / "run_config.json", {}) or {},
            "files": _local_run_file_checks(run_dir)}


@app.post("/api/runs/example/load")
def load_example() -> Dict[str, Any]:
    if not EXAMPLE_RUN_DIR.exists():
        raise HTTPException(status_code=404,
                            detail="Example dataset not found. Please open a final run folder.")
    # Read-only: the freeze ships prebuilt website_indices; never rewrite them.
    ensure_indices(EXAMPLE_RUN_DIR, rebuild=False)
    return set_current_run("example", EXAMPLE_RUN_DIR,
                           "Example · FGFR2 IIIb/IIIc — 30 vertebrates", run_id="example")


@app.post("/api/runs/open")
def open_run(req: OpenRunRequest) -> Dict[str, Any]:
    p = (PROJECT_ROOT / req.path).resolve() if not Path(req.path).is_absolute() else Path(req.path).resolve()
    if not str(p).startswith(str(PROJECT_ROOT.resolve())):
        raise HTTPException(status_code=400, detail="Path must be inside the project.")
    # allow either the closure dir itself or a run base containing it
    if (p / "final_pre_interpro_truth_table.tsv").exists():
        run_dir = p
    elif (closure_dir_for(p) / "final_pre_interpro_truth_table.tsv").exists():
        run_dir = closure_dir_for(p)
    else:
        raise HTTPException(status_code=404, detail="No final_pre_interpro_truth_table.tsv found in that folder.")
    ensure_indices(run_dir, rebuild=True)
    return set_current_run("previous", run_dir, f"Run · {run_dir.parent.name}",
                           run_id=read_json(run_dir / "final_pre_interpro_run_mode.json", {}).get("run_id", ""))


@app.post("/api/runs")
def create_run(req: RunRequest) -> Dict[str, Any]:
    raise HTTPException(
        status_code=410,
        detail=("Automatic web runs are disabled. Create a run with "
                "POST /api/local-runs/create (runs/<run_id>/), then run "
                "`python scripts/run_pre_interpro_for_run.py --run-id <run_id>` "
                "yourself in a local terminal."),
    )


@app.get("/api/runs/{run_id}/status")
def run_status(run_id: str) -> Dict[str, Any]:
    if run_id in ("example", "current"):
        cur = get_current_run()
        if not cur:
            raise HTTPException(status_code=404, detail="No run loaded.")
        run_dir = Path(cur["run_dir"])
        steps = bwi.build_run_index(run_dir).get("steps", [])
        return {"run_id": cur.get("run_id", run_id), "kind": cur.get("kind"),
                "status": "finished", "steps": steps,
                "run_mode": read_json(run_dir / "final_pre_interpro_run_mode.json", {})}
    st = read_json(run_status_path(run_id), None)
    if st is None:
        raise HTTPException(status_code=404, detail=f"Run not found: {run_id}")
    # merge live step status from the closure dir if present
    cdir = closure_dir_for(run_state_dir(run_id))
    if (cdir / "final_pre_interpro_step_status.tsv").exists():
        try:
            st["steps"] = bwi.build_run_index(cdir).get("steps", [])
        except Exception:
            pass
        st["run_mode"] = read_json(cdir / "final_pre_interpro_run_mode.json", {})
    return st


@app.get("/api/runs/{run_id}/logs")
def run_logs(run_id: str, tail: int = Query(400, ge=1, le=8000)) -> Dict[str, Any]:
    lp = run_log_path(run_id)
    lines: List[str] = []
    if lp.exists():
        lines = lp.read_text(encoding="utf-8", errors="replace").splitlines()[-tail:]
    return {"run_id": run_id, "lines": lines}


# Dataset endpoints accept example or run:<run_id> selectors.
_DATASET_MODEL_INFLIGHT_LOCK = Lock()
_DATASET_MODEL_INFLIGHT: Dict[str, Future[Dict[str, Any]]] = {}


@app.get("/api/runs/current/dataset-model")
def current_dataset_model(dataset: Optional[str] = Query(None)) -> Dict[str, Any]:
    selected = dataset or DATASET_EXAMPLE
    with _DATASET_MODEL_INFLIGHT_LOCK:
        pending = _DATASET_MODEL_INFLIGHT.get(selected)
        owner = pending is None
        if owner:
            pending = Future()
            _DATASET_MODEL_INFLIGHT[selected] = pending
    if not owner:
        return pending.result()
    try:
        model = build_canonical_dataset_model(resolve_dataset(selected))
        result = _prune_missing_file_links(model, selected)
        pending.set_result(result)
        return result
    except Exception as exc:
        pending.set_exception(exc)
        raise
    finally:
        with _DATASET_MODEL_INFLIGHT_LOCK:
            if _DATASET_MODEL_INFLIGHT.get(selected) is pending:
                _DATASET_MODEL_INFLIGHT.pop(selected, None)


@app.get("/api/runs/current/summary")
def current_summary(rebuild: bool = False, dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("run_index.json", dataset, rebuild=rebuild)


@app.get("/api/runs/current/evidence-stack")
def current_evidence_stack(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("evidence_stack.json", dataset)


@app.get("/api/runs/current/species")
def current_species(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("species_index.json", dataset)


@app.get("/api/runs/current/species/{species}")
def current_species_one(species: str, dataset: Optional[str] = Query(None)) -> Any:
    data = index_for_request("species_index.json", dataset)
    key = species.lower()
    for row in data:
        if row.get("species", "").lower() == key or row.get("display_species_name", "").lower() == key:
            return row
    raise HTTPException(status_code=404, detail=f"Species not found: {species}")


@app.get("/api/runs/current/figures")
def current_figures(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("figure_index.json", dataset)


@app.get("/api/runs/current/downloads")
def current_downloads(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("download_index.json", dataset)


@app.get("/api/runs/current/freeze")
def current_freeze(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("freeze_index.json", dataset)


# Interactive cassette, coordinate, MSA, synteny and story indices.
def _species_slice(index_name: str, species: str, dataset: Optional[str] = None) -> Dict[str, Any]:
    data = index_for_request(index_name, dataset)
    key = species.lower()
    for row in data.get("species", []):
        if row.get("species", "").lower() == key or row.get("display_species_name", "").lower() == key:
            return row
    raise HTTPException(status_code=404, detail=f"Species not found in {index_name}: {species}")


@app.get("/api/runs/current/cassette")
def current_cassette(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("cassette_residue_index.json", dataset)


@app.get("/api/runs/current/species/{species}/cassette")
def current_species_cassette(species: str, dataset: Optional[str] = Query(None)) -> Any:
    return _species_slice("cassette_residue_index.json", species, dataset)


@app.get("/api/runs/current/coordinates")
def current_coordinates(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("coordinate_track_index.json", dataset)


@app.get("/api/runs/current/species/{species}/coordinates")
def current_species_coordinates(species: str, dataset: Optional[str] = Query(None)) -> Any:
    return _species_slice("coordinate_track_index.json", species, dataset)


@app.get("/api/runs/current/msa")
def current_msa(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("msa_index.json", dataset)


@app.get("/api/runs/current/synteny")
def current_synteny(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("synteny_locus_index.json", dataset)


@app.get("/api/runs/current/species/{species}/synteny")
def current_species_synteny(species: str, dataset: Optional[str] = Query(None)) -> Any:
    return _species_slice("synteny_locus_index.json", species, dataset)


@app.get("/api/runs/current/story")
def current_story(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("species_story_index.json", dataset)


@app.get("/api/runs/current/species/{species}/story")
def current_species_story(species: str, dataset: Optional[str] = Query(None)) -> Any:
    return _species_slice("species_story_index.json", species, dataset)


# --------------------------------------------------------------------------- #
# Post-InterPro / pyTMHMM domain-architecture indices (step 15)
# --------------------------------------------------------------------------- #
@app.get("/api/runs/current/domain-architecture")
def current_domain_architecture(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("domain_architecture_index.json", dataset)


@app.get("/api/runs/current/domain-architecture/summary")
def current_domain_architecture_summary(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("domain_architecture_summary.json", dataset)


@app.get("/api/runs/current/domain-architecture/species")
def current_domain_architecture_species(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("species_domain_architecture.json", dataset)


@app.get("/api/runs/current/domain-architecture/qc")
def current_domain_architecture_qc(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("domain_architecture_qc.json", dataset)


@app.get("/api/runs/current/species/{species}/domain-architecture")
def current_species_domain_architecture(species: str, dataset: Optional[str] = Query(None)) -> Any:
    return _species_slice("species_domain_architecture.json", species, dataset)


# Core-only datasets read their generic run-local indices.
def _resolve_shared_index(run_base: Path, name: str) -> Any:
    for p in (run_base / "website_indices" / name,
              run_base / "website_indices" / "generic" / name):
        data = read_json(p, None)
        if data is not None:
            return data
    return None


def _core_generic_index(name: str, dataset: Optional[str]) -> Any:
    ds = resolve_dataset(dataset)
    if ds["kind"] != "run":
        raise HTTPException(status_code=404, detail="Core indices are only available for runs.")
    data = _resolve_shared_index(ds["run_base"], name)
    if data is None:
        raise HTTPException(status_code=404, detail=f"Core index not available: {name}")
    return _with_availability(data, name, ds)


# Shared indices use run-local data; FGFR2 retains validated event endpoints.
_SHARED_INDEX_FILES = {
    "overview": "overview_index.json",
    "evidence-stack": "evidence_stack.json",
    "gene-explorer": "gene_explorer_index.json",
    "protein-architecture": "protein_architecture_index.json",
    "synteny": "synteny_index.json",
    "event-evidence": "event_evidence_index.json",
    "domain-architecture": "domain_architecture_index.json",
    "exon-domain-boundaries": "exon_domain_boundaries_index.json",
    "figures": "figures_index.json",
    "available-views": "available_views.json",
}


def _shared_index(name: str, dataset: Optional[str]) -> Any:
    ds = resolve_dataset(dataset)
    if ds["kind"] != "run":
        raise HTTPException(status_code=404,
                            detail="Shared indices are only available for runs.")
    fname = _SHARED_INDEX_FILES.get(name)
    if not fname:
        raise HTTPException(status_code=404, detail=f"Unknown shared index: {name}")
    data = _resolve_shared_index(ds["run_base"], fname)
    if data is None:
        raise HTTPException(status_code=404, detail=f"Shared index not available: {name}")
    return _with_availability(data, fname, ds)


@app.get("/api/runs/current/shared/{name}")
def current_shared_index(name: str, dataset: Optional[str] = Query(None)) -> Any:
    return _shared_index(name, dataset)


@app.get("/api/runs/current/core/summary")
def current_core_summary(dataset: Optional[str] = Query(None)) -> Any:
    return _core_generic_index("dataset_summary.json", dataset)


@app.get("/api/runs/current/core/evidence-stack")
def current_core_evidence_stack(dataset: Optional[str] = Query(None)) -> Any:
    return _core_generic_index("evidence_stack.json", dataset)


@app.get("/api/runs/current/core/primary-selection")
def current_core_primary_selection(dataset: Optional[str] = Query(None)) -> Any:
    return _core_generic_index("primary_selection_index.json", dataset)


@app.get("/api/runs/current/core/protein-architecture")
def current_core_protein_architecture(dataset: Optional[str] = Query(None)) -> Any:
    ds = resolve_dataset(dataset)
    if ds["kind"] == "run":
        data = _resolve_shared_index(ds["run_base"], "protein_architecture_index.json")
        if data is not None:
            return data
    return current_core_exon_protein_map(dataset)


@app.get("/api/runs/current/core/event-evidence-index")
def current_core_event_evidence_index(dataset: Optional[str] = Query(None)) -> Any:
    return _core_generic_index("event_evidence_index.json", dataset)


@app.get("/api/runs/current/core/figures")
def current_core_figures(dataset: Optional[str] = Query(None)) -> Any:
    return _core_generic_index("figures_index.json", dataset)


@app.get("/api/runs/current/core/gene-analysis")
def current_core_gene_analysis(dataset: Optional[str] = Query(None)) -> Any:
    return _core_generic_index("gene_analysis_index.json", dataset)


@app.get("/api/runs/current/core/domain-architecture")
def current_core_domain_architecture(dataset: Optional[str] = Query(None)) -> Any:
    return _core_generic_index("domain_architecture_index.json", dataset)


@app.get("/api/runs/current/core/synteny")
def current_core_synteny(dataset: Optional[str] = Query(None)) -> Any:
    return _core_generic_index("synteny_index.json", dataset)


@app.get("/api/runs/current/core/exon-domain-boundaries")
def current_core_exon_domain_boundaries(dataset: Optional[str] = Query(None)) -> Any:
    return _core_generic_index("exon_domain_boundary_index.json", dataset)


@app.get("/api/runs/current/core/event-candidates")
def current_core_event_candidates(dataset: Optional[str] = Query(None)) -> Any:
    ds = resolve_dataset(dataset)
    if ds["kind"] != "run":
        return {"available": False, "candidates": [], "reason": "not_a_run"}
    p = ds["run_base"] / "results" / "core_gene_analysis" / "event_candidate_regions.tsv"
    if not p.is_file():
        return {"available": False, "candidates": [],
                "reason": "no_event_candidate_scan",
                "message": "No event region configured or detected. Core gene-level analysis is available."}
    rows = []
    try:
        import csv as _csv
        with open(p, encoding="utf-8", newline="") as fh:
            rows = list(_csv.DictReader(fh, delimiter="\t"))
    except Exception:
        rows = []
    return {
        "available": bool(rows),
        "exploratory": True,
        "candidates": rows,
        "message": ("Potential isoform-specific regions detected. These are exploratory "
                    "and not validated event regions." if rows else
                    "No event region configured or detected. Core gene-level analysis is available."),
    }


@app.get("/api/runs/current/core/exon-protein-map")
def current_core_exon_protein_map(dataset: Optional[str] = Query(None)) -> Any:
    ds = resolve_dataset(dataset)
    if ds["kind"] != "run":
        raise HTTPException(status_code=404, detail="Exon/protein track is only available for runs.")
    import csv as _csv
    core = ds["run_base"] / "results" / "core_gene_analysis"

    def _rows(name: str) -> List[Dict[str, str]]:
        p = core / name
        if not p.is_file():
            return []
        try:
            with open(p, encoding="utf-8", newline="") as fh:
                return list(_csv.DictReader(fh, delimiter="\t"))
        except Exception:
            return []

    def _int(v: Any) -> Optional[int]:
        try:
            return int(float(v))
        except (TypeError, ValueError):
            return None

    gene_symbol = ""
    iso_rows = _rows("protein_isoform_index.tsv")
    gm_rows = _rows("gene_model_index.tsv")
    if gm_rows:
        gene_symbol = gm_rows[0].get("gene_symbol", "")
    exon_rows = _rows("exon_protein_map.tsv")
    domain_rows = _rows("domain_features.tsv")
    tm_rows = _rows("tm_features.tsv")
    cluster_rows = _rows("event_region_candidate_clusters.tsv")

    # protein length + role + transcript from the isoform index
    meta: Dict[str, Dict[str, Any]] = {}
    for r in iso_rows:
        pid = r.get("protein_id", "")
        if not pid:
            continue
        meta[pid] = {
            "transcript_id": r.get("transcript_id", ""),
            "length_aa": _int(r.get("protein_length")),
            "role": "primary" if str(r.get("primary_status", "")).lower() == "primary" else "alternative",
        }

    selection_method = ""
    rc = read_json(ds["run_base"] / "run_config.json", {}) or {}
    selection_method = rc.get("selection_method", "")

    def exons_for(pid: str) -> List[Dict[str, Any]]:
        out = []
        for e in exon_rows:
            if e.get("protein_id") != pid:
                continue
            out.append({
                "exon_id": e.get("exon_id", ""),
                "exon_number": _int(e.get("exon_number")),
                "protein_start_aa": _int(e.get("protein_start_aa")),
                "protein_end_aa": _int(e.get("protein_end_aa")),
                "confidence": e.get("confidence", ""),
            })
        out.sort(key=lambda x: (x["protein_start_aa"] or 0))
        return out

    def domains_for(pid: str) -> List[Dict[str, Any]]:
        return [{
            "domain_source": d.get("domain_source", ""), "domain_id": d.get("domain_id", ""),
            "domain_name": d.get("domain_name", ""),
            "start_aa": _int(d.get("start_aa")), "end_aa": _int(d.get("end_aa")),
        } for d in domain_rows if d.get("protein_id") == pid]

    def tms_for(pid: str) -> List[Dict[str, Any]]:
        return [{"start_aa": _int(t.get("start_aa")), "end_aa": _int(t.get("end_aa")),
                 "source": t.get("source", "")} for t in tm_rows if t.get("protein_id") == pid]

    # Candidate clusters overlaid on the PRIMARY protein (the representative
    # coordinates are relative to the reference/primary isoform). Exploratory only.
    def clusters_for(pid: str, is_primary: bool) -> List[Dict[str, Any]]:
        if not is_primary:
            return []
        out = []
        for c in cluster_rows:
            out.append({
                "candidate_cluster_id": c.get("candidate_cluster_id", ""),
                "start_aa": _int(c.get("representative_start_aa")),
                "end_aa": _int(c.get("representative_end_aa")),
                "length_aa": _int(c.get("representative_length_aa")),
                "support_count": _int(c.get("support_count")),
                "confidence": c.get("confidence", ""),
                "exon_aligned": (_int(c.get("exon_aligned_support")) or 0) > 0,
                "evidence_status": "exploratory",
            })
        return out

    # group by species
    by_species: Dict[str, Dict[str, Any]] = {}
    for r in gm_rows:
        pid = r.get("protein_id", "")
        if not pid or str(r.get("model_status", "")) not in ("protein_coding",):
            continue
        sp = r.get("species_id", "")
        node = by_species.setdefault(sp, {"species_id": sp, "proteins": []})
        if any(p["protein_id"] == pid for p in node["proteins"]):
            continue
        m = meta.get(pid, {})
        is_primary = m.get("role") == "primary"
        node["proteins"].append({
            "protein_id": pid,
            "transcript_id": m.get("transcript_id", r.get("transcript_id", "")),
            "length_aa": m.get("length_aa") or _int(r.get("protein_length")),
            "role": m.get("role", "alternative"),
            "selection_method": selection_method if is_primary else "",
            "exons": exons_for(pid),
            "candidate_regions": clusters_for(pid, is_primary),
            "domains": domains_for(pid),
            "tm_regions": tms_for(pid),
        })
    # primary first within each species
    for node in by_species.values():
        node["proteins"].sort(key=lambda p: (p["role"] != "primary", p["protein_id"]))

    domain_status = "available" if domain_rows else "pending_cluster"
    return {
        "gene_symbol": gene_symbol,
        "species": list(by_species.values()),
        "domain_status": domain_status,
        "selection_method": selection_method,
    }


def _read_core_tsv(dataset: Optional[str], filename: str) -> Dict[str, Any]:
    ds = resolve_dataset(dataset)
    if ds["kind"] != "run":
        return {"available": False, "rows": [], "reason": "not_a_run"}
    p = ds["run_base"] / "results" / "core_gene_analysis" / filename
    if not p.is_file():
        return {"available": False, "rows": [], "reason": "file_not_present"}
    try:
        import csv as _csv
        with open(p, encoding="utf-8", newline="") as fh:
            rows = list(_csv.DictReader(fh, delimiter="\t"))
    except Exception:
        rows = []
    return {"available": bool(rows), "rows": rows, "reason": "" if rows else "empty"}


@app.get("/api/runs/current/core/event-evidence")
def current_core_event_evidence(dataset: Optional[str] = Query(None)) -> Any:
    ev = _read_core_tsv(dataset, "event_region_evidence.tsv")
    clusters = _read_core_tsv(dataset, "event_region_candidate_clusters.tsv")
    # legacy raw pairwise table kept available for transparency
    raw = _read_core_tsv(dataset, "event_candidate_regions.tsv")
    uniprot_report = None
    ds = resolve_dataset(dataset)
    if ds["kind"] == "run":
        uniprot_report = read_json(
            ds["run_base"] / "results" / "core_gene_analysis"
            / "uniprot_event_evidence_report.json", None)
    has_any = ev["available"] or clusters["available"] or raw["available"]
    return {
        "available": has_any,
        "exploratory": True,
        "evidence_status": "exploratory",
        "clusters": clusters["rows"],
        "evidence": ev["rows"],
        "raw_candidates": raw["rows"],
        "n_clusters": len(clusters["rows"]),
        "n_evidence": len(ev["rows"]),
        "n_raw_candidates": len(raw["rows"]),
        "uniprot": uniprot_report,
        "explanation": ("Isoform A and Isoform B are two protein isoforms compared within the "
                        "same gene. Candidate regions are sequence differences between isoforms. "
                        "They may be biologically interesting, but they are exploratory and NOT "
                        "validated event regions."),
        "message": ("Exploratory isoform-difference evidence. These candidate regions are not "
                    "validated events and do not enable event-specific analysis." if has_any else
                    "No isoform-difference evidence for this gene. Core gene-level analysis is available."),
    }


@app.get("/api/runs/current/core/capability")
def current_core_capability(dataset: Optional[str] = Query(None)) -> Any:
    ds = resolve_dataset(dataset)
    if ds["kind"] != "run":
        raise HTTPException(status_code=404, detail="Capability report is only available for runs.")
    p = ds["run_base"] / "website_indices" / "generic" / "gene_capability_report.json"
    data = read_json(p, None)
    if data is None:
        raise HTTPException(status_code=404, detail="gene_capability_report.json not available.")
    return data


@app.get("/api/local-runs/{run_id}/core-validate")
def local_run_core_validate(run_id: str) -> Dict[str, Any]:
    run_dir = _safe_local_run_dir(run_id)
    if not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Run not found: {run_id}")
    if _evaluate_core_run is None:
        raise HTTPException(status_code=500, detail="Milestone evaluator unavailable.")
    return _evaluate_core_run(run_dir)


# --------------------------------------------------------------------------- #
# Module 1 — exon-domain boundary consistency indices (step 16)
# --------------------------------------------------------------------------- #
@app.get("/api/runs/current/boundary-consistency")
def current_boundary_consistency(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("boundary_consistency_index.json", dataset)


@app.get("/api/runs/current/boundary-consistency/summary")
def current_boundary_consistency_summary(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("boundary_consistency_summary.json", dataset)


@app.get("/api/runs/current/boundary-consistency/matrix")
def current_boundary_consistency_matrix(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("boundary_consistency_matrix.json", dataset)


@app.get("/api/runs/current/boundary-consistency/outliers")
def current_boundary_consistency_outliers(dataset: Optional[str] = Query(None)) -> Any:
    return index_for_request("boundary_consistency_outliers.json", dataset)


# --------------------------------------------------------------------------- #
# Configurable scientific package builder.
# --------------------------------------------------------------------------- #
class PackageSelection(BaseModel):
    preset: str = "recommended"
    scope: str = ""
    items: List[str] = Field(default_factory=list)
    species: List[str] = Field(default_factory=list)
    completed_species_only: bool = False
    formats: List[str] = Field(
        default_factory=lambda: ["tsv", "xlsx", "faa", "json", "svg", "pdf", "png"])


def _resolve_active_run_dir(dataset: Optional[str] = None) -> Path:
    if dataset:
        desc = resolve_dataset(dataset)
        if desc.get("kind") == "example":
            raise HTTPException(
                status_code=400,
                detail="Package builder is for analysis runs; the FGFR2 example "
                       "keeps its validated Files view.",
            )
        return Path(desc["run_base"])
    ptr = read_json(CURRENT_RUN_PTR, None) or {}
    rid = ptr.get("run_id") or ""
    if rid and rid != "example":
        p = _safe_local_run_dir(rid)
        if p.is_dir():
            return p
    # Dataset query param may also arrive as the active switcher id.
    raise HTTPException(status_code=404, detail="No active run for package builder")


@app.get("/api/runs/current/package-capabilities")
def package_capabilities(scope: Optional[str] = Query(None),
                         dataset: Optional[str] = Query(None)) -> Any:
    sys.path.insert(0, str(SCRIPTS_DIR / "shared_gene_analysis"))
    from exondomaincompare.shared_gene_analysis.package_builder import capabilities  # type: ignore
    return capabilities(_resolve_active_run_dir(dataset), scope)


@app.get("/api/runs/current/package-catalogue")
def package_catalogue(scope: Optional[str] = Query(None),
                      dataset: Optional[str] = Query(None)) -> Any:
    return package_capabilities(scope=scope, dataset=dataset)


@app.post("/api/runs/current/packages")
def create_package(body: PackageSelection, dataset: Optional[str] = Query(None)) -> Any:
    sys.path.insert(0, str(SCRIPTS_DIR / "shared_gene_analysis"))
    from exondomaincompare.shared_gene_analysis.package_builder import PackageJob, build_package  # type: ignore
    run_dir = _resolve_active_run_dir(dataset)
    job = PackageJob(
        job_id=f"pkg_{uuid.uuid4().hex[:12]}",
        run_id=run_dir.name,
        created_at=datetime.now(timezone.utc).isoformat(),
    )
    selection = body.model_dump() if hasattr(body, "model_dump") else body.dict()
    job = build_package(run_dir, selection, job=job)
    return _job_payload(job)


@app.get("/api/runs/current/packages/{job_id}")
def package_status(job_id: str) -> Any:
    sys.path.insert(0, str(SCRIPTS_DIR / "shared_gene_analysis"))
    from exondomaincompare.shared_gene_analysis.package_builder import get_job  # type: ignore
    job = get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Package job not found: {job_id}")
    return _job_payload(job)


def _job_payload(job: Any) -> Dict[str, Any]:
    return {
        "job_id": job.job_id,
        "run_id": job.run_id,
        "status": job.status,
        "progress": job.progress,
        "message": job.message,
        "error": job.error,
        "scope": job.scope,
        "preset": job.preset,
        "package_name": job.package_name,
        "zip_path": job.zip_path,
        "n_files": job.n_files,
        "estimated_bytes": job.estimated_bytes,
        "selected_items": job.selected_items,
        "skipped_items": job.skipped_items,
        "selected_species": job.selected_species,
        "warnings": job.warnings,
        "manifest": job.manifest,
        "created_at": job.created_at,
        "updated_at": job.updated_at,
    }


# --------------------------------------------------------------------------- #
# file serving (sandboxed to project root)
# --------------------------------------------------------------------------- #
def _resolve_dataset_file_path(path: str, dataset: str) -> Optional[Path]:
    logical = Path(path)
    if not dataset:
        return None
    descriptor = resolve_dataset(dataset)
    run_id = str(descriptor.get("run_id") or "")
    if descriptor.get("kind") == "example":
        roots = [
            Path(RESULTS_ROOT),
            Path(descriptor["closure_dir"]),
            Path(EXAMPLE_DERIVED_INDICES_DIR).parent,
        ]
    else:
        roots = [Path(descriptor["run_base"]), Path(descriptor["closure_dir"])]

    resolved_roots = [base.resolve() for base in roots]
    if logical.is_absolute():
        absolute = logical.resolve()
        for root in resolved_roots:
            try:
                absolute.relative_to(root)
                return absolute
            except ValueError:
                pass
        if run_id and run_id in logical.parts:
            marker = logical.parts.index(run_id)
            logical = Path(*logical.parts[marker + 1:])
        else:
            return None

    relative = logical
    if len(logical.parts) >= 2 and logical.parts[0] == "runs" \
            and logical.parts[1] == run_id:
        relative = Path(*logical.parts[2:])
    elif descriptor.get("kind") == "example":
        for prefix in (
            Path("results") / "final_30_until_interpro_prepare",
            Path("datasets") / "fgfr2_30_species",
        ):
            try:
                relative = logical.relative_to(prefix)
                break
            except ValueError:
                pass
        derived_prefix = Path("results") / "derived" / "example"
        try:
            derived_relative = logical.relative_to(derived_prefix)
            relative = Path("derived") / derived_relative
        except ValueError:
            pass

    safe_candidates: List[Path] = []
    for root in resolved_roots:
        candidate = (root / relative).resolve()
        try:
            candidate.relative_to(root)
        except ValueError:
            continue
        safe_candidates.append(candidate)
        if candidate.is_file():
            return candidate
    return safe_candidates[0] if safe_candidates else None


def _resolve_public_file_path(path: str, dataset: Optional[str] = None) -> Path:
    if path.startswith("package:"):
        relative = Path(path[len("package:"):])
        candidate = (RUNTIME_CONFIG.paths.packages / relative).resolve()
        root = RUNTIME_CONFIG.paths.packages.resolve()
        try:
            candidate.relative_to(root)
        except ValueError:
            raise HTTPException(status_code=400, detail="Invalid package path") from None
        return candidate
    dataset_candidate = _resolve_dataset_file_path(path, dataset or "")
    if dataset_candidate is not None and dataset_candidate.is_file():
        return dataset_candidate
    logical = Path(path)
    legacy_fgfr2 = Path("results") / "final_30_until_interpro_prepare"
    legacy_derived = Path("results") / "derived" / "example"
    try:
        bundled_suffix = logical.relative_to(legacy_fgfr2)
    except ValueError:
        bundled_suffix = None
    try:
        derived_suffix = logical.relative_to(legacy_derived)
    except ValueError:
        derived_suffix = None
    if derived_suffix is not None and RESULTS_ROOT == BUNDLED_FGFR2_ROOT:
        root = (RESULTS_ROOT / "derived").resolve()
        candidate = (root / derived_suffix).resolve()
    elif bundled_suffix is not None and RESULTS_ROOT == BUNDLED_FGFR2_ROOT:
        root = RESULTS_ROOT.resolve()
        candidate = (root / bundled_suffix).resolve()
    elif len(logical.parts) >= 2 and logical.parts[0] == "runs":
        root = _safe_local_run_dir(logical.parts[1]).resolve()
        candidate = root.joinpath(*logical.parts[2:]).resolve()
    else:
        candidate = (PROJECT_ROOT / logical).resolve()
        root = PROJECT_ROOT.resolve()
    try:
        candidate.relative_to(root)
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid path") from None
    return candidate


_DIRECT_FILE_LINK_KEYS = {
    "file", "report", "source_file", "source_table", "thumbnail",
    "png", "svg", "pdf",
}
_FILE_LINK_MAP_KEYS = {"formats", "source_tables"}
_DOWNLOAD_SUFFIXES = {
    ".csv", ".faa", ".fasta", ".json", ".md", ".pdf", ".png",
    ".svg", ".tsv", ".txt", ".xlsx", ".zip",
}


def _looks_like_file_link(value: str) -> bool:
    return value.startswith("package:") or Path(value.split("?", 1)[0]).suffix.lower() \
        in _DOWNLOAD_SUFFIXES


def _file_link_available(path: str, dataset: str,
                         availability: Dict[str, bool]) -> bool:
    if path not in availability:
        try:
            availability[path] = _resolve_public_file_path(
                path, dataset=dataset
            ).is_file()
        except HTTPException:
            availability[path] = False
    return availability[path]


def _prune_missing_file_links(value: Any, dataset: str,
                              availability: Optional[Dict[str, bool]] = None) -> Any:
    if availability is None:
        availability = {}
    if isinstance(value, list):
        cleaned = []
        for item in value:
            if isinstance(item, dict) and item.get("path") and item.get("format"):
                if not _file_link_available(
                    str(item["path"]), dataset, availability
                ):
                    continue
            cleaned.append(_prune_missing_file_links(item, dataset, availability))
        return cleaned
    if not isinstance(value, dict):
        return value
    cleaned: Dict[str, Any] = {}
    for key, item in value.items():
        if key in _DIRECT_FILE_LINK_KEYS and isinstance(item, str) and item \
                and _looks_like_file_link(item) \
                and not item.startswith(("/api/", "http://", "https://")):
            if not _file_link_available(item, dataset, availability):
                continue
        if key in _FILE_LINK_MAP_KEYS and isinstance(item, dict):
            mapped = {}
            for item_key, path in item.items():
                if not isinstance(path, str) or not path:
                    continue
                if path.startswith(("/api/", "http://", "https://")):
                    mapped[item_key] = path
                    continue
                if not _looks_like_file_link(path):
                    continue
                if _file_link_available(path, dataset, availability):
                    mapped[item_key] = path
            cleaned[key] = mapped
        else:
            cleaned[key] = _prune_missing_file_links(item, dataset, availability)
    return cleaned


@app.get("/api/runs/{run_id}/files")
def run_local_file(run_id: str, path: str, inline: bool = False):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", run_id or ""):
        raise HTTPException(status_code=400, detail="Invalid run id")
    base = _safe_local_run_dir(run_id).resolve()
    requested = (base / path).resolve()
    try:
        requested.relative_to(base)
    except ValueError:
        raise HTTPException(status_code=400, detail="Path must stay inside the selected run")
    if not requested.is_file():
        raise HTTPException(status_code=404, detail=f"Run-local file not found: {path}")
    return FileResponse(str(requested), filename=None if inline else requested.name)


@app.get("/api/download")
def download(path: str, inline: bool = False,
             dataset: Optional[str] = Query(None)):
    p = _resolve_public_file_path(path, dataset=dataset)
    if not p.exists() or not p.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    if inline:
        return FileResponse(str(p))
    return FileResponse(str(p), filename=p.name)


@app.get("/api/file-preview")
def file_preview(path: str, max_bytes: int = 20000) -> Dict[str, Any]:
    p = _resolve_public_file_path(path)
    if not p.exists() or not p.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    size = p.stat().st_size
    cap = max(1, min(int(max_bytes or 20000), 200000))
    try:
        with p.open("r", encoding="utf-8", errors="replace") as fh:
            text = fh.read(cap + 1)
    except Exception as exc:  # pragma: no cover
        raise HTTPException(status_code=400, detail=f"Cannot preview file: {exc}")
    return {
        "path": path, "name": p.name, "format": p.suffix.lstrip(".").lower(),
        "size_bytes": size, "truncated": len(text) > cap, "text": text[:cap],
    }
