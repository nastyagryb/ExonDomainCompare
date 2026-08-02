#!/usr/bin/env bash
###############################################################################
# FGFR2 IIIb/IIIc comparative pipeline runner — v3 (Tasks 1-6 hardened)
#
# Improvements over run_fgfr2_pipeline_current_v2_22_references_fixed.sh:
#   * Strict mode + per-step error trap (no silent failures).
#   * Directory names match the committed result layout exactly, so the runner
#     can resume from cached inputs without path guessing.
#   * CLI calls verified against each script's actual argparse (the old runner
#     mis-called Step 5b [missing --cache], Step 6 export, and Step 9 [wrong
#     script name + wrong args]).
#   * Step 5b and Step 6 now receive the NCBI datasets --cache so protein
#     sequences resolve locally instead of failing.
#   * The coordinate resolver (Tasks 3-5) runs as an explicit step that writes
#     the coordinate audit + pair-level QC summary, before the figures step.
#   * Every step verifies its key outputs and records them in a run manifest
#     (JSON), with a final human-readable summary.
#
# Usage:
#   chmod +x run_fgfr2_pipeline_current_v3.sh
#   ./run_fgfr2_pipeline_current_v3.sh
#   BASE=results/final_30_until_interpro_prepare RUN_FROM=5b RUN_TO=10 ./run_fgfr2_pipeline_current_v3.sh
#   FORCE=1 RUN_FROM=10 ./run_fgfr2_pipeline_current_v3.sh
#
# Environment overrides: BASE, SCRIPTS, PYTHON, CACHE, RUN_FROM, RUN_TO, FORCE,
#   RUN_STEP7, USE_NCBI_DATASETS, DATASETS_BIN, THREADS, MAX_MAIN_SPECIES,
#   NO_ENSEMBL_REST (1 to disable remote protein fetches), SPECIES_LIST,
#   HUMAN_IIIB_FASTA, HUMAN_IIIC_FASTA, FGFR_PARALOG_REF_FASTA.
###############################################################################
set -Eeuo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE="${BASE:-results/final_30_until_interpro_prepare}"
SCRIPTS="${SCRIPTS:-scripts}"
PYTHON="${PYTHON:-python}"
FORCE="${FORCE:-0}"
RUN_FROM="${RUN_FROM:-1}"
RUN_TO="${RUN_TO:-11}"               # full pre-InterProScan endpoint (Tasks 8-12)
RUN_STEP7="${RUN_STEP7:-1}"          # InterPro FASTA prep (on: scope ends at InterPro preparation)
USE_NCBI_DATASETS="${USE_NCBI_DATASETS:-1}"
DATASETS_BIN="${DATASETS_BIN:-datasets}"
THREADS="${THREADS:-4}"
MAX_MAIN_SPECIES="${MAX_MAIN_SPECIES:-16}"
NO_ENSEMBL_REST="${NO_ENSEMBL_REST:-0}"
REF_DIR="${REF_DIR:-references}"
FGFR2_III_SEGMENTS_DIR="${FGFR2_III_SEGMENTS_DIR:-${REF_DIR}/fgfr2_iii_segments}"
GENE_SYMBOL="${GENE_SYMBOL:-FGFR2}"

# Directory layout — matches the committed result tree so cached reuse works.
D01="${BASE}/01_species_registry"
D02="${BASE}/02_models"
D03="${BASE}/03_selection_initial"
D04="${BASE}/04_isoform_evidence_v2_3_human_calibrated"
D05="${BASE}/05_selection_with_isoforms_v2_3_human_calibrated"
D05B="${BASE}/05b_selection_with_isoforms_v2_7_marker_validated"
D06="${BASE}/06_protein_export_v2_7_marker_validated"
D06B="${BASE}/06b_paralog_screen_v2_7_marker_validated"
D06E="${BASE}/06e_III_region_anchor_v5_9_marker_validated"
D07="${BASE}/07_interpro_prepare_v2_7_marker_validated"
D06F="${BASE}/06f_orthology_evidence"          # Addendum B orthology evidence
D09="${BASE}/09_paper_ready_qc_v2_9"
D10="${D10:-${D09}/figures_v2_22_final_qc_display}"
D11="${D11:-${BASE}/11_pre_interpro_master}"  # Tasks 8/10/11/12 pre-InterPro package
# Addendum A: multi-vertebrate FGFR1/2/3/4 paralog panel (preferred; human-only kept as legacy control).
MULTI_PANEL_FASTA="${MULTI_PANEL_FASTA:-${REF_DIR}/fgfr_paralog_panel_multi_vertebrate.fasta}"
ALLOW_PANEL_BUILD="${ALLOW_PANEL_BUILD:-1}"     # build panel via NCBI datasets if missing

# NCBI datasets cache produced by Step 2 (holds protein.faa + gff per species).
CACHE="${CACHE:-${D02}/_ncbi_datasets_cache}"

# Reference inputs.
# Robust species-list resolution (final_30 closure):
#   1) $SPECIES_LIST if provided AND present; 2) reference/Species_list_final_30.txt;
#   3) case-insensitive search for Species_list_final_30.txt under project/reference dirs;
#   4) species_list.txt fallback; 5) clear error listing every checked path.
resolve_species_list() {
  local checked=()
  if [[ -n "${SPECIES_LIST:-}" ]]; then
    if [[ -s "$SPECIES_LIST" ]]; then echo "$SPECIES_LIST"; return 0; fi
    checked+=("$SPECIES_LIST (from \$SPECIES_LIST, missing)")
  fi
  local cands=(
    "reference/Species_list_final_30.txt" "reference/species_list_final_30.txt"
    "${REF_DIR}/Species_list_final_30.txt" "${REF_DIR}/species_list_final_30.txt"
    "./Species_list_final_30.txt" "./species_list_final_30.txt"
  )
  local c
  for c in "${cands[@]}"; do
    checked+=("$c")
    [[ -s "$c" ]] && { echo "$c"; return 0; }
  done
  # case-insensitive search under common roots
  local hit
  hit="$(find reference "${REF_DIR}" . -maxdepth 2 -iname 'Species_list_final_30.txt' -type f 2>/dev/null | head -n1)"
  if [[ -n "$hit" && -s "$hit" ]]; then echo "$hit"; return 0; fi
  checked+=("find -iname Species_list_final_30.txt (none)")
  # last-resort generic fallback
  for c in "${REF_DIR}/species_list.txt" "./species_list.txt" "reference/species_list.txt"; do
    checked+=("$c")
    [[ -s "$c" ]] && { echo "$c"; return 0; }
  done
  echo "[ERROR] No species list found. Checked paths:" >&2
  printf '  - %s\n' "${checked[@]}" >&2
  echo "Set SPECIES_LIST=<path> or place reference/Species_list_final_30.txt." >&2
  return 1
}
SPECIES_LIST="$(resolve_species_list)" || exit 2
echo "[CONFIG] SPECIES_LIST resolved -> ${SPECIES_LIST}"
HUMAN_IIIB_FASTA="${HUMAN_IIIB_FASTA:-${FGFR2_III_SEGMENTS_DIR}/human_FGFR2_IIIb_segment.fasta}"
HUMAN_IIIC_FASTA="${HUMAN_IIIC_FASTA:-${FGFR2_III_SEGMENTS_DIR}/human_FGFR2_IIIc_segment.fasta}"
REFERENCE_MANIFEST="${REFERENCE_MANIFEST:-reference/fgfr2_human_iiib_iiic_manifest.tsv}"
HUMAN_IIIB_PROTEIN_FASTA="${HUMAN_IIIB_PROTEIN_FASTA:-reference/human_FGFR2_IIIb_protein.faa}"
HUMAN_IIIC_PROTEIN_FASTA="${HUMAN_IIIC_PROTEIN_FASTA:-reference/human_FGFR2_IIIc_protein.faa}"
if [[ -z "${FGFR_PARALOG_REF_FASTA:-}" ]]; then
  if [[ -s "${REF_DIR}/human_FGFR1_2_3_4.fasta" ]]; then
    FGFR_PARALOG_REF_FASTA="${REF_DIR}/human_FGFR1_2_3_4.fasta"
  else
    FGFR_PARALOG_REF_FASTA="${REF_DIR}/human_FGFR1_2_3_4_raw.fasta"
  fi
fi

MANIFEST="${MANIFEST:-${BASE}/run_manifest_v3.json}"

# ---------------------------------------------------------------------------
# Logging + manifest helpers
# ---------------------------------------------------------------------------
log(){ printf '\n[RUN] %s\n' "$*"; }
warn(){ printf '[WARN] %s\n' "$*" >&2; }
fail(){ printf '[ERROR] %s\n' "$*" >&2; exit 2; }
script(){ [[ -s "${SCRIPTS}/$1" ]] || fail "Missing script: ${SCRIPTS}/$1"; }
need(){ [[ -s "$1" ]] || fail "Required input missing/empty: $1"; }

MANIFEST_TMP="$(mktemp)"
echo "[]" > "$MANIFEST_TMP"
CURRENT_STEP="init"
FAILED_STEPS=()

manifest_add(){
  # manifest_add <step> <status> <detail>
  "$PYTHON" - "$MANIFEST_TMP" "$1" "$2" "$3" <<'PY'
import json, sys, datetime
path, step, status, detail = sys.argv[1:5]
data = json.load(open(path))
data.append({"step": step, "status": status, "detail": detail,
             "timestamp": datetime.datetime.now().isoformat(timespec="seconds")})
json.dump(data, open(path, "w"), indent=2)
PY
}

on_err(){
  local rc=$?
  warn "Step '${CURRENT_STEP}' failed with exit code ${rc}"
  manifest_add "$CURRENT_STEP" "failed" "exit_code=${rc}"
  finalize_manifest
  exit "$rc"
}
trap on_err ERR

finalize_manifest(){
  mkdir -p "$(dirname "$MANIFEST")"
  cp "$MANIFEST_TMP" "$MANIFEST"
}

check_outputs(){
  # check_outputs <step> <file...>  -> fail if any missing/empty
  local step="$1"; shift
  local missing=()
  for f in "$@"; do
    [[ -s "$f" ]] || missing+=("$f")
  done
  if (( ${#missing[@]} )); then
    fail "Step ${step} did not produce expected outputs: ${missing[*]}"
  fi
  manifest_add "$step" "ok" "outputs=$*"
}

# ---------------------------------------------------------------------------
# Step ordering
# ---------------------------------------------------------------------------
step_rank(){
  case "$1" in
    1) echo 10;; 2) echo 20;; 3) echo 30;; 4) echo 40;; 5) echo 50;; 5b) echo 55;;
    6) echo 60;; 6b) echo 65;; 6e) echo 68;; 7) echo 70;; 9) echo 90;; 10) echo 100;; 11) echo 110;;
    *) fail "Unknown step label: $1";;
  esac
}
should_run(){
  local r from to
  r=$(step_rank "$1"); from=$(step_rank "$RUN_FROM"); to=$(step_rank "$RUN_TO")
  [[ "$r" -ge "$from" && "$r" -le "$to" ]]
}
run_if_needed(){
  local step="$1" sentinel="$2"; shift 2
  CURRENT_STEP="$step"
  if ! should_run "$step"; then echo "[SKIP] Step $step outside RUN_FROM/RUN_TO"; return 0; fi
  if [[ "$FORCE" != "1" && -s "$sentinel" ]]; then
    echo "[SKIP] Step $step output exists: $sentinel"
    manifest_add "$step" "cached" "$sentinel"
    return 0
  fi
  "$@"
}

mkdir -p "$BASE"
echo "[CONFIG] BASE=$BASE"
echo "[CONFIG] SCRIPTS=$SCRIPTS  PYTHON=$($PYTHON --version 2>&1)"
echo "[CONFIG] CACHE=$CACHE"
echo "[CONFIG] RUN_FROM=$RUN_FROM RUN_TO=$RUN_TO FORCE=$FORCE NO_ENSEMBL_REST=$NO_ENSEMBL_REST"
echo "[CONFIG] HUMAN_IIIB_FASTA=$HUMAN_IIIB_FASTA"
echo "[CONFIG] HUMAN_IIIC_FASTA=$HUMAN_IIIC_FASTA"
echo "[CONFIG] FGFR_PARALOG_REF_FASTA=$FGFR_PARALOG_REF_FASTA"

# Whether Step 5b must fail (vs warn) when the NCBI cache is missing.
STEP5B_REQUIRE_CACHE="${STEP5B_REQUIRE_CACHE:-1}"

ENSEMBL_ARG=(); [[ "$NO_ENSEMBL_REST" == "1" ]] && ENSEMBL_ARG=(--no_ensembl_rest)
# NOTE: do NOT compute a cache argument globally here. In a from-scratch run the
# NCBI datasets cache is only created by Step 2, so any cache test before Step 2
# would be stale. Each step that needs the cache resolves it at call time.

# ---------------------------------------------------------------------------
# 1. Species registry
# ---------------------------------------------------------------------------
step1_run(){
  log "Step 1: species registry"
  local s="build_species_registry.py"
  [[ -s "${SCRIPTS}/${s}" ]] || s="build_species_registry_improved.py"
  script "$s"; need "$SPECIES_LIST"; mkdir -p "$D01"
  "$PYTHON" "${SCRIPTS}/${s}" --species_list "$SPECIES_LIST" --outdir "$D01"
  check_outputs 1 "$D01/species_registry.tsv"
}
run_if_needed 1 "$D01/species_registry.tsv" step1_run

# ---------------------------------------------------------------------------
# 2. Dual-source FGFR2 collection (models + CDS features + datasets cache)
# ---------------------------------------------------------------------------
step2_run(){
  log "Step 2: collect FGFR2 models + CDS features"
  script "collect_fgfr2_models_dual_source_v3.py"; need "$D01/species_registry.tsv"; mkdir -p "$D02"
  local DS=(); [[ "$USE_NCBI_DATASETS" == "1" ]] && DS=(--use_ncbi_datasets --datasets_bin "$DATASETS_BIN")
  local GC=(); [[ -s "${GENE_CONFIG:-}" ]] && GC=(--gene_config "$GENE_CONFIG")
  "$PYTHON" "$SCRIPTS/collect_fgfr2_models_dual_source_v3.py" \
    --species_registry "$D01/species_registry.tsv" \
    --gene_symbol "$GENE_SYMBOL" \
    --outdir "$D02" \
    --paralog_reference_fasta "$FGFR_PARALOG_REF_FASTA" \
    ${GC[@]+"${GC[@]}"} \
    ${DS[@]+"${DS[@]}"}
  check_outputs 2 "$D02/cds_features.tsv" "$D02/transcripts.tsv" "$D02/exons.tsv"
}
run_if_needed 2 "$D02/cds_features.tsv" step2_run

# ---------------------------------------------------------------------------
# 3. Initial transcript selection
# ---------------------------------------------------------------------------
step3_run(){
  log "Step 3: initial annotation-aware transcript selection"
  script "select_fgfr2_transcripts_annotation_aware_v2.py"
  need "$D02/transcripts.tsv"; need "$D02/exons.tsv"; mkdir -p "$D03"
  "$PYTHON" "$SCRIPTS/select_fgfr2_transcripts_annotation_aware_v2.py" \
    --transcripts "$D02/transcripts.tsv" --exons "$D02/exons.tsv" \
    --outdir "$D03" --gene_symbol "$GENE_SYMBOL"
  check_outputs 3 "$D03/selected_transcripts.tsv"
}
run_if_needed 3 "$D03/selected_transcripts.tsv" step3_run

# ---------------------------------------------------------------------------
# 4. IIIb/IIIc exon-structure classification
# ---------------------------------------------------------------------------
step4_run(){
  log "Step 4: classify IIIb/IIIc by exon structure"
  script "classify_fgfr2_IIIb_IIIc_by_exon_structure_v2_3_human_calibrated.py"
  need "$D02/transcripts.tsv"; need "$D02/exons.tsv"; need "$D03/selected_transcripts.tsv"; mkdir -p "$D04"
  need "$HUMAN_IIIB_FASTA"; need "$HUMAN_IIIC_FASTA"
  # Resolve the cache argument at call time (Step 2 may have just created it).
  local CACHE_ARG=()
  if [[ -d "$CACHE" ]]; then
    CACHE_ARG=(--cache "$CACHE")
  elif [[ "$STEP5B_REQUIRE_CACHE" == "1" ]]; then
    fail "Step 4: required NCBI datasets cache missing at $CACHE (set STEP5B_REQUIRE_CACHE=0 to proceed without sequence calibration)."
  else
    warn "Step 4: NCBI cache not found at $CACHE; direction calibration may fall back to provisional order rule."
  fi
  "$PYTHON" "$SCRIPTS/classify_fgfr2_IIIb_IIIc_by_exon_structure_v2_3_human_calibrated.py" \
    --transcripts "$D02/transcripts.tsv" --exons "$D02/exons.tsv" \
    --selected_transcripts "$D03/selected_transcripts.tsv" \
    --cds_features "$D02/cds_features.tsv" \
    --human_iiib_segment_fasta "$HUMAN_IIIB_FASTA" \
    --human_iiic_segment_fasta "$HUMAN_IIIC_FASTA" \
    ${CACHE_ARG[@]+"${CACHE_ARG[@]}"} ${ENSEMBL_ARG[@]+"${ENSEMBL_ARG[@]}"} \
    --outdir "$D04" --gene_symbol "$GENE_SYMBOL"
  check_outputs 4 "$D04/fgfr2_isoform_evidence.tsv"
}
run_if_needed 4 "$D04/fgfr2_isoform_evidence.tsv" step4_run

# ---------------------------------------------------------------------------
# 5. Re-selection with isoform evidence
# ---------------------------------------------------------------------------
step5_run(){
  log "Step 5: re-select transcripts with IIIb/IIIc evidence"
  script "select_fgfr2_transcripts_annotation_aware_v2.py"
  need "$D02/transcripts.tsv"; need "$D02/exons.tsv"; need "$D04/fgfr2_isoform_evidence.tsv"; mkdir -p "$D05"
  "$PYTHON" "$SCRIPTS/select_fgfr2_transcripts_annotation_aware_v2.py" \
    --transcripts "$D02/transcripts.tsv" --exons "$D02/exons.tsv" \
    --isoform_evidence "$D04/fgfr2_isoform_evidence.tsv" \
    --outdir "$D05" --gene_symbol "$GENE_SYMBOL"
  check_outputs 5 "$D05/selected_transcripts.tsv"
}
run_if_needed 5 "$D05/selected_transcripts.tsv" step5_run

# ---------------------------------------------------------------------------
# 5b. Protein-level IIIb/IIIc validation  (Task 2 — now receives --cache)
# ---------------------------------------------------------------------------
step5b_run(){
  log "Step 5b: protein-validate IIIb/IIIc candidates"
  script "protein_validate_fgfr2_III_candidate_selection_v2_6.py"
  need "$D05/selected_transcripts.tsv"; need "$D02/transcripts.tsv"; need "$D04/fgfr2_isoform_evidence.tsv"
  need "$HUMAN_IIIB_FASTA"; need "$HUMAN_IIIC_FASTA"; mkdir -p "$D05B"
  # Resolve the cache argument at call time (Step 2 may have just created it).
  local CACHE_ARG=()
  if [[ -d "$CACHE" ]]; then
    CACHE_ARG=(--cache "$CACHE")
  elif [[ "$STEP5B_REQUIRE_CACHE" == "1" ]]; then
    fail "Step 5b: required NCBI datasets cache missing at $CACHE (set STEP5B_REQUIRE_CACHE=0 to proceed without it)."
  else
    warn "Step 5b: NCBI cache not found at $CACHE; NCBI proteins may be unavailable."
  fi
  local MANIFEST_ARG=(); [[ -s "$REFERENCE_MANIFEST" ]] && MANIFEST_ARG=(--reference_manifest "$REFERENCE_MANIFEST")
  local HUMANPROT_ARG=()
  if [[ -s "$HUMAN_IIIB_PROTEIN_FASTA" && -s "$HUMAN_IIIC_PROTEIN_FASTA" ]]; then
    HUMANPROT_ARG=(--human_iiib_protein_fasta "$HUMAN_IIIB_PROTEIN_FASTA" --human_iiic_protein_fasta "$HUMAN_IIIC_PROTEIN_FASTA")
  fi
  "$PYTHON" "$SCRIPTS/protein_validate_fgfr2_III_candidate_selection_v2_6.py" \
    --selected "$D05/selected_transcripts.tsv" \
    --transcripts "$D02/transcripts.tsv" \
    --isoform_evidence "$D04/fgfr2_isoform_evidence.tsv" \
    --human_iiib_segment_fasta "$HUMAN_IIIB_FASTA" \
    --human_iiic_segment_fasta "$HUMAN_IIIC_FASTA" \
    ${MANIFEST_ARG[@]+"${MANIFEST_ARG[@]}"} ${HUMANPROT_ARG[@]+"${HUMANPROT_ARG[@]}"} ${CACHE_ARG[@]+"${CACHE_ARG[@]}"} ${ENSEMBL_ARG[@]+"${ENSEMBL_ARG[@]}"} \
    --outdir "$D05B" \
    --output_selected "$D05B/selected_transcripts.tsv"
  check_outputs 5b "$D05B/selected_transcripts.tsv" \
    "$D05B/fgfr2_III_candidate_protein_validation.tsv" \
    "$D05B/fgfr2_III_candidate_protein_validation_summary.tsv" \
    "$D05B/fgfr2_III_candidate_protein_validation_warnings.tsv" \
    "$D05B/fgfr2_III_candidate_protein_validation_debug_controls.tsv" \
    "$D05B/fgfr2_III_close_primate_control_diagnostics.tsv" \
    "$D05B/fgfr2_III_final_selected_protein_validation_summary.tsv" \
    "$D05B/fgfr2_III_reference_manifest_check.tsv"
}
run_if_needed 5b "$D05B/selected_transcripts.tsv" step5b_run

# ---------------------------------------------------------------------------
# 6. Export selected proteins  (real CLI: --cache/--out/--report/--region_qc...)
# ---------------------------------------------------------------------------
step6_run(){
  log "Step 6: export selected FGFR2 proteins"
  script "export_selected_fgfr2_proteins_complete_v2_1_region_qc.py"
  need "$D05B/selected_transcripts.tsv"; mkdir -p "$D06"
  [[ -d "$CACHE" ]] || fail "Step 6 requires NCBI datasets cache at $CACHE"
  "$PYTHON" "$SCRIPTS/export_selected_fgfr2_proteins_complete_v2_1_region_qc.py" \
    --selected "$D05B/selected_transcripts.tsv" \
    --cache "$CACHE" \
    --out "$D06/selected_fgfr2_proteins.faa" \
    --report "$D06/protein_export_report.tsv" \
    --md_report "$D06/protein_export_report.md" \
    --html_report "$D06/protein_export_report.html" \
    --metadata "$D06/protein_export_metadata.json" \
    --region_qc "$D06/fgfr2_exported_protein_region_qc.tsv" \
    --warnings "$D06/protein_export_warnings.tsv" \
    ${ENSEMBL_ARG[@]+"${ENSEMBL_ARG[@]}"}
  check_outputs 6 "$D06/selected_fgfr2_proteins.faa" "$D06/protein_export_report.tsv"
}
run_if_needed 6 "$D06/selected_fgfr2_proteins.faa" step6_run

SELECTED_PROTEIN_FASTA="${SELECTED_PROTEIN_FASTA:-$D06/selected_fgfr2_proteins.faa}"
PROTEIN_REPORT="${PROTEIN_REPORT:-$D06/protein_export_report.tsv}"

# ---------------------------------------------------------------------------
# 6b. FGFR paralog screen (Task 7 deferred — kept as-is, cached if present)
# ---------------------------------------------------------------------------
step6b_run(){
  log "Step 6b: FGFR paralog screen (legacy human-only control + multi-vertebrate panel + orthology)"
  script "validate_fgfr2_paralog_screen_v2_cautious.py"
  script "build_fgfr_paralog_panel_multi_vertebrate.py"
  script "screen_fgfr2_paralogs_multi_vertebrate.py"
  script "build_fgfr2_orthology_evidence.py"
  need "$SELECTED_PROTEIN_FASTA"; need "$FGFR_PARALOG_REF_FASTA"; mkdir -p "$D06B" "$D06F"
  local REPORT_ARG=(); [[ -s "$PROTEIN_REPORT" ]] && REPORT_ARG=(--protein_report "$PROTEIN_REPORT")

  # 6b-1: legacy human-only FGFR1/2/3/4 panel screen (kept as legacy control).
  "$PYTHON" "$SCRIPTS/validate_fgfr2_paralog_screen_v2_cautious.py" \
    --query_fasta "$SELECTED_PROTEIN_FASTA" \
    --reference_fasta "$FGFR_PARALOG_REF_FASTA" \
    ${REPORT_ARG[@]+"${REPORT_ARG[@]}"} \
    --outdir "$D06B" --prefix fgfr2 --threads "$THREADS" --fallback_pairwise
  # Preserve legacy canonical-named outputs before the multi-vertebrate screen overwrites them.
  for n in species_summary warnings; do
    [[ -s "$D06B/fgfr2_paralog_screen_${n}.tsv" ]] && cp "$D06B/fgfr2_paralog_screen_${n}.tsv" "$D06B/legacy_human_panel_paralog_screen_${n}.tsv"
  done

  # 6b-2: build multi-vertebrate FGFR1/2/3/4 panel if missing (NCBI datasets; cached/reproducible).
  if [[ ! -s "$MULTI_PANEL_FASTA" ]]; then
    if [[ "$ALLOW_PANEL_BUILD" == "1" ]]; then
      "$PYTHON" "$SCRIPTS/build_fgfr_paralog_panel_multi_vertebrate.py" \
        --panel_fasta "$MULTI_PANEL_FASTA" --outdir "$D06B" --datasets_bin "$DATASETS_BIN" || \
        warn "Multi-vertebrate panel build failed; will fall back to human-only legacy panel"
    else
      warn "Multi-vertebrate panel missing and ALLOW_PANEL_BUILD=0"
    fi
  else
    echo "[RUN] Using existing multi-vertebrate panel: $MULTI_PANEL_FASTA"
  fi

  # 6b-3: multi-vertebrate paralog screen (writes canonical detailed/species_summary/warnings).
  local PANEL="$MULTI_PANEL_FASTA"; [[ -s "$PANEL" ]] || PANEL="$FGFR_PARALOG_REF_FASTA"
  "$PYTHON" "$SCRIPTS/screen_fgfr2_paralogs_multi_vertebrate.py" \
    --query_fasta "$SELECTED_PROTEIN_FASTA" \
    --panel_fasta "$PANEL" \
    --outdir "$D06B" --threads "$THREADS"

  # Note: orthology evidence (Addendum B) is built in Step 11, after Step 10
  # produces the pair-level QC it depends on.
  check_outputs 6b \
    "$D06B/fgfr2_paralog_screen_detailed.tsv" \
    "$D06B/fgfr2_paralog_screen_species_summary.tsv" \
    "$D06B/fgfr2_paralog_screen_warnings.tsv"
}
run_if_needed 6b "$D06B/fgfr2_paralog_screen_detailed.tsv" step6b_run

# ---------------------------------------------------------------------------
# 6e. Dynamic human-calibrated IIIb/IIIc anchoring (Task 6 — similarity classes)
# ---------------------------------------------------------------------------
step6e_run(){
  log "Step 6e: dynamic IIIb/IIIc protein-region anchors + similarity classes"
  script "map_human_IIIb_IIIc_region_to_orthologs_FINAL_v5_7_dynamic_human_qc.py"
  need "$SELECTED_PROTEIN_FASTA"; need "$HUMAN_IIIB_FASTA"; need "$HUMAN_IIIC_FASTA"; mkdir -p "$D06E"
  "$PYTHON" "$SCRIPTS/map_human_IIIb_IIIc_region_to_orthologs_FINAL_v5_7_dynamic_human_qc.py" \
    --query_fasta "$SELECTED_PROTEIN_FASTA" \
    --human_iiib_segment_fasta "$HUMAN_IIIB_FASTA" \
    --human_iiic_segment_fasta "$HUMAN_IIIC_FASTA" \
    --outdir "$D06E" --prefix fgfr2 \
    --exons_tsv "$D02/exons.tsv" \
    --isoform_evidence_tsv "$D04/fgfr2_isoform_evidence.tsv"
  check_outputs 6e "$D06E/fgfr2_III_pair_audit.tsv"
}
run_if_needed 6e "$D06E/fgfr2_III_pair_audit.tsv" step6e_run

# ---------------------------------------------------------------------------
# 7. InterProScan-ready clean FASTA (optional; off unless RUN_STEP7=1)
# ---------------------------------------------------------------------------
step7_run(){
  log "Step 7: prepare InterProScan-ready clean FASTA"
  script "prepare_interpro_clean_fasta_v2.py"; need "$SELECTED_PROTEIN_FASTA"; mkdir -p "$D07"
  local REPORT_ARG=(); [[ -s "$PROTEIN_REPORT" ]] && REPORT_ARG=(--protein_export_report "$PROTEIN_REPORT")
  "$PYTHON" "$SCRIPTS/prepare_interpro_clean_fasta_v2.py" \
    --input "$SELECTED_PROTEIN_FASTA" --outdir "$D07" --prefix FGFR2 --split_size 50 ${REPORT_ARG[@]+"${REPORT_ARG[@]}"}
  check_outputs 7 \
    "$D07/fgfr2_interpro_clean_unique.fasta" \
    "$D07/fgfr2_interpro_unique_mapping.tsv" \
    "$D07/fgfr2_interpro_id_mapping.tsv" \
    "$D07/fgfr2_interpro_prepare_summary.tsv" \
    "$D07/fgfr2_interpro_prepare_warnings.tsv" \
    "$D07/interproscan_run_instructions.md" \
    "$D07/interproscan_input_manifest.tsv"
}
if should_run 7; then
  if [[ "$RUN_STEP7" != "1" ]]; then echo "[SKIP] Step 7 disabled (RUN_STEP7=0)"; CURRENT_STEP=7; manifest_add 7 "skipped" "RUN_STEP7=0";
  else run_if_needed 7 "$D07/fgfr2_interpro_clean_unique.fasta" step7_run; fi
fi

# ---------------------------------------------------------------------------
# 9. Paper-ready QC package  (correct script name + real CLI)
# ---------------------------------------------------------------------------
step9_run(){
  log "Step 9: paper-ready QC package"
  script "make_fgfr2_paper_ready_qc_package_v2_9.py"
  need "$D06E/fgfr2_III_pair_audit.tsv"; mkdir -p "$D09"
  local ANCHOR_ARG=(); [[ -s "$D06E/fgfr2_III_region_anchor_map.tsv" ]] && ANCHOR_ARG=(--anchor_map "$D06E/fgfr2_III_region_anchor_map.tsv")
  local DIFF_ARG=(); [[ -s "$D06E/fgfr2_III_pair_difference_positions.tsv" ]] && DIFF_ARG=(--pair_difference_positions "$D06E/fgfr2_III_pair_difference_positions.tsv")
  # Step 9 is validated against the legacy human-only panel summary (expects the
  # exact status "all_high_confidence_FGFR2"); keep that input stable. The
  # multi-vertebrate screen + orthology evidence feed the species QC master instead.
  local PARA_ARG=()
  if [[ -s "$D06B/legacy_human_panel_paralog_screen_species_summary.tsv" ]]; then
    PARA_ARG=(--paralog_species_summary "$D06B/legacy_human_panel_paralog_screen_species_summary.tsv")
  elif [[ -s "$D06B/fgfr2_paralog_screen_species_summary.tsv" ]]; then
    PARA_ARG=(--paralog_species_summary "$D06B/fgfr2_paralog_screen_species_summary.tsv")
  fi
  "$PYTHON" "$SCRIPTS/make_fgfr2_paper_ready_qc_package_v2_9.py" \
    --pair_audit "$D06E/fgfr2_III_pair_audit.tsv" \
    ${ANCHOR_ARG[@]+"${ANCHOR_ARG[@]}"} ${DIFF_ARG[@]+"${DIFF_ARG[@]}"} ${PARA_ARG[@]+"${PARA_ARG[@]}"} \
    --outdir "$D09" --prefix fgfr2
  check_outputs 9 "$D09/fgfr2_paper_ready_species_qc.tsv"
}
run_if_needed 9 "$D09/fgfr2_paper_ready_species_qc.tsv" step9_run

# ---------------------------------------------------------------------------
# 10a. Coordinate resolver (Tasks 3-5): refined labels + native/normalized
#       coordinates + conservative pair-level QC + main-figure eligibility.
# ---------------------------------------------------------------------------
find_alt_exons(){
  if [[ -n "${ALT_EXONS:-}" && -s "$ALT_EXONS" ]]; then echo "$ALT_EXONS"; return; fi
  for c in "$D04/fgfr2_alternative_exon_metadata.tsv" \
           "${BASE}/04_exon_structure/fgfr2_alternative_exon_metadata.tsv"; do
    [[ -s "$c" ]] && { echo "$c"; return; }
  done
  find "$BASE" -path '*fgfr2_alternative_exon_metadata.tsv' -type f 2>/dev/null | head -n 1 || true
}

step10_run(){
  log "Step 10: coordinate resolver + final paper figures"
  script "resolve_fgfr2_IIIb_IIIc_exons_exact_v2_22.py"
  script "make_fgfr2_exact_exon_to_protein_figures_v2_22.py"
  need "$D05B/selected_transcripts.tsv"; need "$D02/cds_features.tsv"; need "$D06E/fgfr2_III_pair_audit.tsv"
  mkdir -p "$D10"
  local alt; alt="$(find_alt_exons)"
  local ALT_ARG=(); [[ -s "$alt" ]] && { ALT_ARG=(--alt_exons "$alt"); echo "[RUN] ALT_EXONS=$alt"; }

  # 10a — resolver (writes coordinate audit + pair-level QC summary).
  "$PYTHON" "$SCRIPTS/resolve_fgfr2_IIIb_IIIc_exons_exact_v2_22.py" \
    --pair_audit "$D06E/fgfr2_III_pair_audit.tsv" \
    --cds_features "$D02/cds_features.tsv" \
    --selected "$D05B/selected_transcripts.tsv" \
    ${ALT_ARG[@]+"${ALT_ARG[@]}"} \
    --outdir "$D10" --prefix fgfr2
  check_outputs 10a \
    "$D10/fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv" \
    "$D10/fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv" \
    "$D10/fgfr2_pair_level_qc_summary.tsv"

  # 10b — figures (visualises the QC fields produced above).
  local QC_ARG=(); [[ -s "$D09/fgfr2_paper_ready_species_qc.tsv" ]] && QC_ARG=(--species_qc "$D09/fgfr2_paper_ready_species_qc.tsv")
  local EXONS_ARG=(); [[ -s "$D02/exons.tsv" ]] && EXONS_ARG=(--exons "$D02/exons.tsv")
  "$PYTHON" "$SCRIPTS/make_fgfr2_exact_exon_to_protein_figures_v2_22.py" \
    --selected "$D05B/selected_transcripts.tsv" \
    --cds_features "$D02/cds_features.tsv" \
    --pair_audit "$D06E/fgfr2_III_pair_audit.tsv" \
    ${ALT_ARG[@]+"${ALT_ARG[@]}"} ${QC_ARG[@]+"${QC_ARG[@]}"} ${EXONS_ARG[@]+"${EXONS_ARG[@]}"} \
    --outdir "$D10" --prefix fgfr2 --max_main_species "$MAX_MAIN_SPECIES"
  check_outputs 10b "$D10/fgfr2_pair_level_qc_summary.tsv"
}
# Force step 10 by default so the resolver always re-derives QC from current 6e.
if should_run 10; then
  CURRENT_STEP=10
  if [[ "$FORCE" != "1" && -s "$D10/fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv" && -s "$D10/fgfr2_pair_level_qc_summary.tsv" ]]; then
    echo "[SKIP] Step 10 outputs exist (set FORCE=1 to rerun): $D10"
    manifest_add 10 "cached" "$D10"
  else
    step10_run
  fi
fi

# ---------------------------------------------------------------------------
# 11. Pre-InterProScan package: species QC master + figure tables + figures + report
#     (Tasks 8, 10, 11, 12). Reads existing outputs; no online protein export.
# ---------------------------------------------------------------------------
step11_run(){
  log "Step 11: orthology evidence + canonical QC master + reproducible phylogenetic order + publication figures"
  script "build_fgfr2_orthology_evidence.py"
  script "build_species_phylogenetic_order.py"
  script "build_species_qc_master_pre_interpro.py"
  script "build_cds_phase_boundary_audit.py"
  script "build_cds_block_cassette_map.py"
  script "build_phase_rescue.py"
  script "build_refined_uncertainty_classes.py"
  script "patch_ncbi_cds_boundaries.py"
  script "make_pre_interpro_figure_tables.py"
  script "final_pre_interpro_validation.py"
  script "make_publication_figures_pre_interpro.py"
  script "make_all_figures.py"
  need "$D04/fgfr2_isoform_evidence.tsv"
  need "$D05B/fgfr2_III_final_selected_protein_validation_summary.tsv"
  need "$D09/fgfr2_paper_ready_species_qc.tsv"
  need "$D10/fgfr2_pair_level_qc_summary.tsv"
  need "$D07/fgfr2_interpro_id_mapping.tsv"
  need "$D06B/fgfr2_paralog_screen_detailed.tsv"
  local D11PUB="$BASE/11_publication_figures_pre_interpro"
  mkdir -p "$D11" "$D06F" "$D11PUB/tables"

  # Addendum B — orthology evidence (needs Step 10 pair QC, built here).
  "$PYTHON" "$SCRIPTS/build_fgfr2_orthology_evidence.py" \
    --genes "$D02/genes.tsv" \
    --paralog_detailed "$D06B/fgfr2_paralog_screen_detailed.tsv" \
    --protein_validation_summary "$D05B/fgfr2_III_final_selected_protein_validation_summary.tsv" \
    --isoform_evidence "$D04/fgfr2_isoform_evidence.tsv" \
    --pair_qc "$D10/fgfr2_pair_level_qc_summary.tsv" \
    --outdir "$D06F"

  local PARALOG_ARG=(); [[ -s "$D06B/fgfr2_paralog_screen_species_summary.tsv" ]] && PARALOG_ARG=(--paralog_summary "$D06B/fgfr2_paralog_screen_species_summary.tsv")
  local ORTHO_ARG=(); [[ -s "$D06F/fgfr2_orthology_species_summary.tsv" ]] && ORTHO_ARG=(--orthology_summary "$D06F/fgfr2_orthology_species_summary.tsv")

  # Sprint Part 2 — reproducible phylogenetic/taxonomic species order (built BEFORE
  # the master so the order columns can be integrated into species_qc_master.tsv).
  "$PYTHON" "$SCRIPTS/build_species_phylogenetic_order.py" \
    --registry "$D01/species_registry.tsv" \
    --outdir "$D11PUB/tables"

  # CDS-boundary sprint Part A — explainable CDS phase/boundary audit (built BEFORE
  # the master so the explainability summary can be integrated). Propagates
  # explainability columns into the final coordinate outputs (no resolver re-run).
  local PROT_ARG=(); [[ -s "$D06/selected_fgfr2_proteins.faa" ]] && PROT_ARG=(--proteins "$D06/selected_fgfr2_proteins.faa")
  "$PYTHON" "$SCRIPTS/build_cds_phase_boundary_audit.py" \
    --coordinate_audit "$D10/fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv" \
    --cds_features "$D02/cds_features.tsv" \
    ${PROT_ARG[@]+"${PROT_ARG[@]}"} \
    --outdir "$D10" \
    --update_coordinate_audit \
    --update_exon_cds_mapping "$D10/fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv" \
    --update_pair_qc "$D10/fgfr2_pair_level_qc_summary.tsv"

  # Task 8 / Addendum E + Sprint Part 2 + Part A — canonical species_qc_master.tsv (+ alias).
  "$PYTHON" "$SCRIPTS/build_species_qc_master_pre_interpro.py" \
    --registry "$D01/species_registry.tsv" \
    --isoform_evidence "$D04/fgfr2_isoform_evidence.tsv" \
    --protein_validation_summary "$D05B/fgfr2_III_final_selected_protein_validation_summary.tsv" \
    --paper_ready_qc "$D09/fgfr2_paper_ready_species_qc.tsv" \
    --pair_qc "$D10/fgfr2_pair_level_qc_summary.tsv" \
    --interpro_id_mapping "$D07/fgfr2_interpro_id_mapping.tsv" \
    --phylo_order "$D11PUB/tables/species_phylogenetic_order.tsv" \
    --cds_audit "$D10/cds_phase_boundary_audit.tsv" \
    ${PARALOG_ARG[@]+"${PARALOG_ARG[@]}"} ${ORTHO_ARG[@]+"${ORTHO_ARG[@]}"} \
    --outdir "$D11" --prefix fgfr2

  # Task 10 — tidy figure input tables (read canonical master; complementary).
  "$PYTHON" "$SCRIPTS/make_pre_interpro_figure_tables.py" \
    --master "$D11/species_qc_master.tsv" \
    --pair_qc "$D10/fgfr2_pair_level_qc_summary.tsv" \
    --outdir "$D11/figure_tables"

  # Sprint Parts 1,3,4,5,6,7 — MANDATORY central entry point: validation gate +
  # publication figures (SVG/PDF/PNG) + plotting tables + captions + manifests + reports.
  "$PYTHON" "$SCRIPTS/make_all_figures.py" --base "$BASE"

  local PUBFIGS=(\
    "Figure_1_framework_pre_interpro" \
    "Figure_2_exon_to_protein_architecture_pre_interpro" \
    "Figure_3_IgIII_cassette_zoom_pre_interpro" \
    "Figure_4_species_evidence_matrix_pre_interpro" \
    "Figure_5_native_vs_normalized_coordinate_qc_pre_interpro" \
    "Supplement_Figure_1_all_species_native_tracks_pre_interpro" \
    "Supplement_Figure_2_review_cases_pre_interpro" \
    "Supplement_Figure_3_interproscan_input_readiness_pre_interpro")
  # Figures may be legitimately skipped for small/custom runs (e.g. a single-species
  # run has no review cases → Supplement Figure 2). make_publication_figures records
  # those with status=skipped_empty in the publication manifest. Such optional figures
  # must NOT be hard-required here, or a valid custom run fails at Step 11. For the
  # full-30 freeze nothing is skipped, so this block is a no-op there.
  local PUBMANIFEST="$D11PUB/metadata/publication_figure_manifest.tsv"
  local SKIPPED_IDS=""
  if [[ -s "$PUBMANIFEST" ]]; then
    SKIPPED_IDS="$(awk -F'\t' 'NR>1 && $2=="skipped_empty"{print $1}' "$PUBMANIFEST")"
  fi
  local FIG_CHECK=()
  for f in "${PUBFIGS[@]}"; do
    local _skip=0 sid
    for sid in $SKIPPED_IDS; do
      [[ -n "$sid" && "$f" == "$sid"* ]] && { _skip=1; break; }
    done
    if [[ "$_skip" == "1" ]]; then
      echo "[INFO] Step 11: figure '$f' is skipped_empty per manifest — treated as optional (small/custom run)"
      continue
    fi
    for ext in svg pdf png; do FIG_CHECK+=("$D11PUB/figures/${f}.${ext}"); done
  done

  check_outputs 11 \
    "$D11/species_qc_master.tsv" \
    "$D11/species_qc_master_pre_interpro.tsv" \
    "$D06F/fgfr2_orthology_evidence.tsv" \
    "$D06F/fgfr2_orthology_species_summary.tsv" \
    "$D10/cds_phase_boundary_audit.tsv" \
    "$D10/cds_phase_boundary_legacy_vs_refined_counts.tsv" \
    "$D10/cds_phase_boundary_explainability_summary.tsv" \
    "$D10/fgfr2_unique_cds_block_table.tsv" \
    "$D10/fgfr2_cassette_cds_block_map.tsv" \
    "$D10/fgfr2_transcript_cds_reconstruction_audit.tsv" \
    "$D10/fgfr2_cassette_coordinate_sanity_audit.tsv" \
    "$D10/fgfr2_ncbi_cds_boundary_patch_report.tsv" \
    "$D10/cds_phase_rescue_audit.tsv" \
    "$D10/fgfr2_refined_uncertainty_classes.tsv" \
    "$D11PUB/tables/species_phylogenetic_order.tsv" \
    "$D11PUB/tables/species_taxonomy_metadata.tsv" \
    "$D11PUB/tables/figure1_framework_counts_pre_interpro.tsv" \
    "$D11PUB/tables/figure2_exon_to_protein_architecture_tracks.tsv" \
    "$D11PUB/tables/figure3_igIII_cassette_zoom_tracks.tsv" \
    "$D11PUB/tables/figure4_species_evidence_matrix.tsv" \
    "$D11PUB/tables/figure5_native_vs_normalized_coordinate_qc.tsv" \
    "$D11PUB/tables/supplement_interproscan_input_readiness.tsv" \
    "${FIG_CHECK[@]}" \
    "$D11PUB/captions/figure_captions_pre_interpro.md" \
    "$D11PUB/metadata/publication_figure_manifest.tsv" \
    "$D11PUB/metadata/final_pre_interpro_validation_report.tsv" \
    "$D11PUB/metadata/final_pre_interpro_validation_report.json" \
    "$D11PUB/metadata/output_file_manifest_pre_interpro.tsv" \
    "$D11PUB/QC_migration_report_tasks_7_to_12_pre_interpro.md" \
    "$D11PUB/methods_update_pre_interpro.md" \
    "$D11PUB/results_summary_pre_interpro.md"
}
if should_run 11; then
  CURRENT_STEP=11
  if [[ "$FORCE" != "1" && -s "$D11/species_qc_master_pre_interpro.tsv" && -s "$BASE/11_publication_figures_pre_interpro/QC_migration_report_tasks_7_to_12_pre_interpro.md" ]]; then
    echo "[SKIP] Step 11 outputs exist (set FORCE=1 to rerun): $D11"
    manifest_add 11 "cached" "$D11"
  else
    step11_run
  fi
fi

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
finalize_manifest
log "Pipeline finished"
echo "[SUMMARY] Run manifest: $MANIFEST"
echo "[SUMMARY] Final QC/figure directory: $D10"
if [[ -d "$D10" ]]; then
  find "$D10" -maxdepth 1 -type f \( -name '*.tsv' -o -name '*.json' \) 2>/dev/null | sort || true
fi
if [[ -d "$D11" ]]; then
  echo "[SUMMARY] Canonical QC master: $D11/species_qc_master.tsv (+ pre_interpro alias)"
fi
if [[ -d "$BASE/11_publication_figures_pre_interpro" ]]; then
  echo "[SUMMARY] Publication package: $BASE/11_publication_figures_pre_interpro"
  echo "          - figures/ : Figures 1-5 + Supplements 1-3 in SVG/PDF/PNG"
  echo "          - tables/ : plotting tables + species_phylogenetic_order.tsv"
  echo "          - captions/figure_captions_pre_interpro.md"
  echo "          - metadata/ : publication_figure_manifest.tsv, final_pre_interpro_validation_report.{tsv,json}"
  echo "          - QC_migration_report_tasks_7_to_12_pre_interpro.md, results_summary, methods_update"
  echo "          - InterProScan-ready FASTA: $D07/fgfr2_interpro_clean_unique.fasta"
  echo "          - InterProScan NOT executed yet (domain annotation pending)"
fi
echo "[SUMMARY] Step statuses:"
"$PYTHON" - "$MANIFEST" <<'PY'
import json, sys
for e in json.load(open(sys.argv[1])):
    print(f"  - step {e['step']:>3}: {e['status']:<8} {e.get('detail','')[:90]}")
PY
