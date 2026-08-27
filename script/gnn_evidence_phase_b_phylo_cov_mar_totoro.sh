#!/bin/bash
# Totoro launcher: Phase B lane 3a — phylo_MAR + covariate_MAR only (NOT MNAR).
# 162 cells × 30 seeds = 4,860 fits; lambda_mode=fixed_1.
#
# G6: only covariate_MAR is MAR-labeled in output RDS (g6_mar_label field).
# Sibling lane runs MNAR via PIGAUTO_LIST_JOBS=mnar launcher.
#
# Usage:
#   cd ~/pigauto_gnn_evidence_phase_b
#   nohup bash script/gnn_evidence_phase_b_phylo_cov_mar_totoro.sh > logs/phase_b_3a.log 2>&1 &
# Poll:
#   bash script/gnn_evidence_phase_b_phylo_cov_mar_totoro.sh status
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"
mkdir -p logs results_phase_b

N_JOBS=4860
NWORKERS="${NWORKERS:-100}"
WALL_CEILING_SEC="${WALL_CEILING_SEC:-12600}"  # ~3.5 h for 2/3 of Phase B @100 workers
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PIGAUTO_PKG_PATH="${ROOT}"
export PIGAUTO_OUT_DIR="${ROOT}/results_phase_b"
export PIGAUTO_LAMBDA_MODE="${PIGAUTO_LAMBDA_MODE:-fixed_1}"
export PIGAUTO_MECHANISM_ARM=MAR
export PIGAUTO_SKIP_EXISTING="${PIGAUTO_SKIP_EXISTING:-1}"
export CUDA_VISIBLE_DEVICES=""

JOB_IDS_FILE="${ROOT}/logs/phase_b_3a_job_ids.txt"

gen_job_ids() {
  PIGAUTO_LIST_JOBS=phylo_cov Rscript "${ROOT}/script/gnn_evidence_campaign_phase_b.R" \
    > "${JOB_IDS_FILE}"
  local n
  n=$(wc -l < "${JOB_IDS_FILE}" | tr -d ' ')
  if [[ "${n}" -ne "${N_JOBS}" ]]; then
    echo "[totoro-3a] FATAL: expected ${N_JOBS} job IDs, got ${n}" >&2
    exit 1
  fi
}

status_cmd() {
  local done fail total mnar_done
  done=$(find "${ROOT}/results_phase_b" -maxdepth 1 -name 'gnn_phase_b_job_*.rds' 2>/dev/null | wc -l | tr -d ' ')
  total="${N_JOBS}"
  fail=0
  if [[ -d "${ROOT}/results_phase_b" ]]; then
    fail=$(Rscript -e '
      fs <- sort(Sys.glob(file.path(Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results_phase_b"),
                                   "gnn_phase_b_job_*.rds")))
      if (!length(fs)) { cat(0); quit(save="no") }
      rows <- lapply(fs, readRDS)
      cat(sum(!vapply(rows, function(r) isTRUE(r$fit_ok), logical(1))))
    ' 2>/dev/null || echo 0)
    mnar_done=$(Rscript -e '
      fs <- sort(Sys.glob(file.path(Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results_phase_b"),
                                   "gnn_phase_b_job_*.rds")))
      if (!length(fs)) { cat(0); quit(save="no") }
      rows <- lapply(fs, readRDS)
      cat(sum(vapply(rows, function(r) identical(r$missing_mechanism, "MNAR"), logical(1))))
    ' 2>/dev/null || echo "?")
  else
    mnar_done="?"
  fi
  echo "[totoro-3a] lane=phylo_MAR+covariate_MAR  progress: ${done}/${total} RDS (all mechanisms on disk)  failures=${fail}"
  echo "[totoro-3a] G6 check: MNAR RDS on disk=${mnar_done} (expect 0 for lane 3a-only completion)"
  if [[ -f "${ROOT}/logs/phase_b_3a.log" ]]; then
    echo "[totoro-3a] log tail:"
    tail -n 5 "${ROOT}/logs/phase_b_3a.log" || true
  fi
}

case "${1:-run}" in
  status)
    status_cmd
    exit 0
    ;;
  run)
    ;;
  *)
    echo "usage: $0 [run|status]" >&2
    exit 1
    ;;
esac

echo "[totoro-3a] ROOT=${ROOT}"
echo "[totoro-3a] lane 3a: phylo_MAR + covariate_MAR  N_JOBS=${N_JOBS}  NWORKERS=${NWORKERS}"
echo "[totoro-3a] G6: covariate_MAR only labeled MAR in RDS g6_mar_label field"
echo "[totoro-3a] G8 wall ceiling: ${WALL_CEILING_SEC}s (~$(( WALL_CEILING_SEC / 3600 ))h)"
echo "[totoro-3a] loadavg: $(cat /proc/loadavg)"
echo "[totoro-3a] git HEAD: $(git -C "${ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "[totoro-3a] started at $(date -Is)"

gen_job_ids
START_TS=$(date +%s)

cat "${JOB_IDS_FILE}" | xargs -P "${NWORKERS}" -I{} bash -c '
  jid={}
  out="'"${ROOT}"'/results_phase_b/gnn_phase_b_job_$(printf "%04d" ${jid}).rds"
  if [[ "'"${PIGAUTO_SKIP_EXISTING}"'" == "1" && -f "${out}" ]]; then exit 0; fi
  PIGAUTO_JOB_ID=${jid} Rscript "'"${ROOT}"'/script/gnn_evidence_campaign_phase_b.R" \
    > "'"${ROOT}"'/logs/gnn_phase_b_job_${jid}.log" 2>&1
'

END_TS=$(date +%s)
WALL=$((END_TS - START_TS))

DONE_3A=$(Rscript -e '
  ids <- as.integer(readLines("'"${JOB_IDS_FILE}"'"))
  out <- Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results_phase_b")
  fs <- file.path(out, sprintf("gnn_phase_b_job_%04d.rds", ids))
  cat(sum(file.exists(fs)))
' 2>/dev/null || echo "?")

FAIL=$(Rscript -e '
  ids <- as.integer(readLines("'"${JOB_IDS_FILE}"'"))
  out <- Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results_phase_b")
  fs <- file.path(out, sprintf("gnn_phase_b_job_%04d.rds", ids))
  fs <- fs[file.exists(fs)]
  if (!length(fs)) { cat(0); quit(save="no") }
  rows <- lapply(fs, readRDS)
  cat(sum(!vapply(rows, function(r) isTRUE(r$fit_ok), logical(1))))
' 2>/dev/null || echo "?")

FAIL_PCT=$(awk -v f="${FAIL}" -v n="${N_JOBS}" 'BEGIN { if (n>0) printf "%.1f", 100*f/n; else print "?" }')
echo "[totoro-3a] lane 3a complete: ${DONE_3A}/${N_JOBS} RDS  failures=${FAIL}/${N_JOBS} (${FAIL_PCT}%)  wall_sec=${WALL}"

if awk "BEGIN { exit !(${FAIL_PCT:-0} > 20) }" 2>/dev/null; then
  echo "[totoro-3a] G8 STOP: >20% failures — pause and report" >&2
  exit 2
fi
if [[ "${WALL}" -gt "${WALL_CEILING_SEC}" ]]; then
  echo "[totoro-3a] G8 STOP: wall > $(( WALL_CEILING_SEC / 3600 ))h — pause and report" >&2
  exit 3
fi

echo "[totoro-3a] done at $(date -Is)"
