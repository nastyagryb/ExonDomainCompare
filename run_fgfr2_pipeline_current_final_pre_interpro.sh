#!/usr/bin/env bash
set -Eeuo pipefail

BASE="${BASE:-results/final_30_until_interpro_prepare}"
SCRIPTS="${SCRIPTS:-scripts}"
PYTHON="${PYTHON:-.venv/bin/python3}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
SKIP_V3="${SKIP_V3:-0}"
SKIP_MSA="${SKIP_MSA:-0}"
FORCE="${FORCE:-0}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

CLOSURE="${BASE}/13_final_pre_interpro_closure"
LOG="${CLOSURE}/final_pre_interpro_run_log.txt"
STATUS="${CLOSURE}/final_pre_interpro_step_status.tsv"
STEP=0
RC=0

mkdir -p "$CLOSURE"
: > "$LOG"
echo -e "step_id\tstep_name\tcommand\tstart_time\tend_time\truntime_seconds\treturn_code\tstatus\toutput_files\twarning_summary" > "$STATUS"

log() { echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$LOG"; }

record_step() {
  local id="$1" name="$2" cmd="$3" t0="$4" t1="$5" rc="$6" st="$7" outs="$8" warn="${9:-}"
  local dur=$(( t1 - t0 ))
  echo -e "${id}\t${name}\t${cmd}\t$(date -u -r "$t0" +%Y-%m-%dT%H:%M:%SZ)\t$(date -u -r "$t1" +%Y-%m-%dT%H:%M:%SZ)\t${dur}\t${rc}\t${st}\t${outs}\t${warn}" >> "$STATUS"
}

run_step() {
  local id="$1" name="$2" cmd="$3"
  STEP=$((STEP + 1))
  log ">>> STEP ${id}: ${name}"
  log "    CMD: ${cmd}"
  local t0 t1 rc st outs warn=""
  t0=$(date +%s)
  if eval "$cmd" >> "$LOG" 2>&1; then
    rc=0; st="pass"
  else
    rc=$?; st="failed"
    RC=$rc
    warn="step failed rc=${rc}"
    log "!!! FAILED step ${id} (rc=${rc})"
  fi
  t1=$(date +%s)
  outs="${4:-}"
  record_step "$id" "$name" "$cmd" "$t0" "$t1" "$rc" "$st" "$outs" "$warn"
  if [[ "$rc" -ne 0 ]]; then
    log "Aborting: hard gate / required step failed."
    exit "$rc"
  fi
}

run_step_optional() {
  local id="$1" name="$2" cmd="$3"
  STEP=$((STEP + 1))
  log ">>> STEP ${id}: ${name} (optional)"
  log "    CMD: ${cmd}"
  local t0 t1 rc st outs warn=""
  t0=$(date +%s)
  if eval "$cmd" >> "$LOG" 2>&1; then
    rc=0; st="pass"
  else
    rc=$?; st="skipped_optional_failed"
    warn="optional step failed rc=${rc}; continuing (not required for FASTA closure)"
    log "~~~ OPTIONAL step ${id} failed (rc=${rc}) — continuing; not required for FASTA closure"
    rc=0
  fi
  t1=$(date +%s)
  outs="${4:-}"
  record_step "$id" "$name" "$cmd" "$t0" "$t1" "$rc" "$st" "$outs" "$warn"
}

# Choose hard vs optional execution for paper figures based on PAPER_FIGURES_OPTIONAL.
PAPER_FIGURES_OPTIONAL="${PAPER_FIGURES_OPTIONAL:-0}"
run_paper_figure_step() {
  if [[ "$PAPER_FIGURES_OPTIONAL" == "1" ]]; then
    run_step_optional "$@"
  else
    run_step "$@"
  fi
}

# Run-mode transparency (Part A)
USED_CACHED_V3="false"; [[ "$SKIP_V3" == "1" ]] && USED_CACHED_V3="true"
USED_CACHED_MSA="false"; [[ "$SKIP_MSA" == "1" ]] && USED_CACHED_MSA="true"

# Robust species-list resolution (final_30): mirrors run_fgfr2_pipeline_current_v3.sh
REF_DIR="${REF_DIR:-references}"
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
  for c in "${cands[@]}"; do checked+=("$c"); [[ -s "$c" ]] && { echo "$c"; return 0; }; done
  local hit
  hit="$(find reference "${REF_DIR}" . -maxdepth 2 -iname 'Species_list_final_30.txt' -type f 2>/dev/null | head -n1)"
  [[ -n "$hit" && -s "$hit" ]] && { echo "$hit"; return 0; }
  for c in "${REF_DIR}/species_list.txt" "./species_list.txt" "reference/species_list.txt"; do
    checked+=("$c"); [[ -s "$c" ]] && { echo "$c"; return 0; }
  done
  echo "[ERROR] No species list found. Checked paths:" >&2
  printf '  - %s\n' "${checked[@]}" >&2
  return 1
}
SPECIES_LIST_RESOLVED="$(resolve_species_list)"
SPECIES_LIST_STATUS="resolved"
if [[ -z "$SPECIES_LIST_RESOLVED" ]]; then
  SPECIES_LIST_RESOLVED=""
  SPECIES_LIST_STATUS="unresolved"
else
  export SPECIES_LIST="$SPECIES_LIST_RESOLVED"
fi

log "================================================================="
log "FGFR2 final pre-InterPro closure run  RUN_ID=${RUN_ID}"
log "BASE=${BASE}"
log "species_list_resolved=${SPECIES_LIST_RESOLVED:-<none>} (${SPECIES_LIST_STATUS})"
log "used_cached_v3_outputs=${USED_CACHED_V3}  used_cached_msa_outputs=${USED_CACHED_MSA}"

# Record species-list resolution EARLY so it is always captured, even if a later
# step (e.g. A1) hard-fails before the closure finalisation block runs.
mkdir -p "$CLOSURE"
record_step "A0" "species_list_resolution" "resolve_species_list" "$(date +%s)" "$(date +%s)" \
  "$([[ "$SPECIES_LIST_STATUS" == "resolved" ]] && echo 0 || echo 1)" \
  "$([[ "$SPECIES_LIST_STATUS" == "resolved" ]] && echo pass || echo failed)" \
  "${SPECIES_LIST_RESOLVED:-<none>}" \
  "species_list_status=${SPECIES_LIST_STATUS}"
cat > "${CLOSURE}/final_pre_interpro_run_mode.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "full_clean_run_completed": false,
  "used_cached_v3_outputs": ${USED_CACHED_V3},
  "used_cached_msa_outputs": ${USED_CACHED_MSA},
  "species_list_resolved": "${SPECIES_LIST_RESOLVED}",
  "species_list_status": "${SPECIES_LIST_STATUS}"
}
EOF
if [[ "$USED_CACHED_V3" == "true" || "$USED_CACHED_MSA" == "true" ]]; then
  log "MODE: CACHED DEBUG RUN — this is NOT a full clean end-to-end run."
else
  log "MODE: FULL CLEAN END-TO-END RUN (Steps 1-11 + MSA/rescue/synteny + closure)."
fi
log "================================================================="

# Clean closure output (fresh freeze each run; upstream cache preserved)
if [[ -d "$CLOSURE" ]]; then
  rm -rf "${CLOSURE}/MSA" "${CLOSURE}/figures" "${CLOSURE}/gates" \
         "${CLOSURE}/freeze" "${CLOSURE}/archive" "${CLOSURE}/tables" 2>/dev/null || true
fi
mkdir -p "$CLOSURE"

# --- Part A: full pipeline ---
if [[ "$SKIP_V3" != "1" ]]; then
  run_step "A1" "pipeline_v3_steps_1_11" \
    "BASE=${BASE} PYTHON=${PYTHON} ./run_fgfr2_pipeline_current_v3.sh" \
    "11_pre_interpro_master/species_qc_master.tsv"
else
  record_step "A1" "pipeline_v3_steps_1_11" "skipped" "$(date +%s)" "$(date +%s)" 0 "skipped_optional" "" "SKIP_V3=1"
  log "SKIP_V3=1 — reusing existing Steps 1-11 outputs"
fi

if [[ "$SKIP_MSA" != "1" ]]; then
  MSA_EXTRA=""
  run_step "A2" "msa_rescue_synteny_module" \
    "${PYTHON} ${SCRIPTS}/run_fgfr2_msa_boundary_module.py --base ${BASE} ${MSA_EXTRA}" \
    "12_msa_boundary_robustness_pre_interpro/maps/fgfr2_post_rescue_final_truth_table.tsv"
else
  record_step "A2" "msa_rescue_synteny_module" "skipped" "$(date +%s)" "$(date +%s)" 0 "skipped_optional" "" "SKIP_MSA=1"
fi

if [[ "$SKIP_V3" != "1" ]] || [[ -f "${BASE}/11_publication_figures_pre_interpro/figures/Figure_1_framework_pre_interpro.svg" ]]; then
  run_step "A2b" "publication_figures_1_4" \
    "${PYTHON} ${SCRIPTS}/make_all_figures.py --base ${BASE}" \
    "11_publication_figures_pre_interpro/figures/Figure_1_framework_pre_interpro.svg"
else
  record_step "A2b" "publication_figures_1_4" "skipped" "$(date +%s)" "$(date +%s)" 0 "skipped_optional" "" "SKIP_V3=1 and no cached pub figures"
fi

run_paper_figure_step "A3" "synteny_paper_figures" \
  "${PYTHON} ${SCRIPTS}/make_fgfr2_synteny_figures_paper.py --base ${BASE}" \
  "12_msa_boundary_robustness_pre_interpro/figures/Figure_9A_FGFR2_local_synteny_5neighbor_paper.svg"

run_paper_figure_step "A4" "framework_evidence_figure" \
  "${PYTHON} ${SCRIPTS}/make_fgfr2_final_framework_figure.py --base ${BASE}" \
  "13_final_pre_interpro_closure/figures/Figure_Final_Framework_Evidence_Stack.svg"

# full_clean_run_completed is only true when A1 and A2 were NOT skipped and RC==0
FULL_CLEAN="false"
if [[ "$USED_CACHED_V3" == "false" && "$USED_CACHED_MSA" == "false" && "$RC" -eq 0 ]]; then
  FULL_CLEAN="true"
fi

cat > "${CLOSURE}/final_pre_interpro_run_mode.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "full_clean_run_completed": ${FULL_CLEAN},
  "used_cached_v3_outputs": ${USED_CACHED_V3},
  "used_cached_msa_outputs": ${USED_CACHED_MSA},
  "species_list_resolved": "${SPECIES_LIST_RESOLVED}",
  "species_list_status": "${SPECIES_LIST_STATUS}"
}
EOF

run_step "A5" "final_closure_sprint" \
  "${PYTHON} ${SCRIPTS}/run_fgfr2_final_pre_interpro_closure.py --base ${BASE} --run-id ${RUN_ID} --run-mode-json ${CLOSURE}/final_pre_interpro_run_mode.json" \
  "13_final_pre_interpro_closure/final_pre_interpro_truth_table.tsv"

# record run-mode flags as a meta row in the step status (Part A)
record_step "A6" "run_mode_summary" "n/a" "$(date +%s)" "$(date +%s)" 0 \
  "$([[ "$FULL_CLEAN" == "true" ]] && echo pass || echo pass_with_warnings)" \
  "final_pre_interpro_run_mode.json" \
  "full_clean_run_completed=${FULL_CLEAN};used_cached_v3_outputs=${USED_CACHED_V3};used_cached_msa_outputs=${USED_CACHED_MSA}"

log "================================================================="
log "full_clean_run_completed=${FULL_CLEAN}"
if [[ "$RC" -eq 0 ]]; then
  log "[DONE] Final pre-InterPro closure completed successfully."
  log "Truth table: ${CLOSURE}/final_pre_interpro_truth_table.tsv"
  log "Primary FASTA: ${CLOSURE}/freeze/final_pre_interpro_proteins_primary.faa"
  log "Archive: ${CLOSURE}/archive/FGFR2_final_pre_interpro_freeze_${RUN_ID}.zip"
  log "Gate: ${CLOSURE}/gates/final_pre_interpro_cross_table_consistency_gate.tsv"
else
  log "[FAIL] Closure run ended with errors (rc=${RC})."
fi
log "Run log: ${LOG}"
log "Step status: ${STATUS}"
exit "$RC"
