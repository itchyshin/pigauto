#!/bin/bash
# Totoro CPU launcher — Phase B lane 3b: MNAR arm only (2,430 fits).
# Parallel with lane 3a (phylo_MAR + covariate_MAR).
#
# Usage on Totoro after rsync:
#   cd ~/pigauto_gnn_evidence_phase_b_mnar
#   nohup bash script/gnn_evidence_phase_b_mnar_totoro.sh > logs/phase_b_mnar.log 2>&1 &
# Poll:
#   bash script/gnn_evidence_phase_b_mnar_totoro.sh status
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"
mkdir -p logs results_phase_b_mnar

N_JOBS=2430
NWORKERS="${NWORKERS:-100}"
WALL_CEILING_SEC="${WALL_CEILING_SEC:-5400}"  # ~1.5 h G8 ceiling (Phase A anchor)
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PIGAUTO_PKG_PATH="${ROOT}"
export PIGAUTO_OUT_DIR="${ROOT}/results_phase_b_mnar"
export PIGAUTO_LAMBDA_MODE="${PIGAUTO_LAMBDA_MODE:-fixed_1}"
export PIGAUTO_SKIP_EXISTING="${PIGAUTO_SKIP_EXISTING:-1}"
export PIGAUTO_MECHANISM_ARM="MNAR"
export PIGAUTO_JOB_PREFIX="gnn_phase_b_mnar"
export CUDA_VISIBLE_DEVICES=""

status_cmd() {
  local done fail total
  done=$(find "${ROOT}/results_phase_b_mnar" -maxdepth 1 -name 'gnn_phase_b_mnar_job_*.rds' 2>/dev/null | wc -l | tr -d ' ')
  total="${N_JOBS}"
  fail=0
  if [[ -d "${ROOT}/results_phase_b_mnar" ]]; then
    fail=$(Rscript -e '
      fs <- sort(Sys.glob(file.path(Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results_phase_b_mnar"),
                                   "gnn_phase_b_mnar_job_*.rds")))
      if (!length(fs)) { cat(0); quit(save="no") }
      rows <- lapply(fs, readRDS)
      cat(sum(!vapply(rows, function(r) isTRUE(r$fit_ok), logical(1))))
    ' 2>/dev/null || echo 0)
  fi
  echo "[totoro] Phase B MNAR progress: ${done}/${total} RDS  failures=${fail}"
  if [[ -f "${ROOT}/logs/phase_b_mnar.log" ]]; then
    echo "[totoro] phase_b_mnar.log tail:"
    tail -n 5 "${ROOT}/logs/phase_b_mnar.log" || true
  fi
  if [[ -f "${ROOT}/results_phase_b_mnar/gnn_phase_b_mnar_summary.csv" ]]; then
    echo "[totoro] summary present: results_phase_b_mnar/gnn_phase_b_mnar_summary.csv"
  fi
}

collect_cmd() {
  echo "[totoro] collecting Phase B MNAR summary"
  PIGAUTO_RESULTS_DIR="${ROOT}/results_phase_b_mnar" \
    Rscript "${ROOT}/script/collect_gnn_evidence_phase_b_mnar.R"
}

case "${1:-run}" in
  status)
    status_cmd
    exit 0
    ;;
  collect)
    collect_cmd
    exit 0
    ;;
  run)
    ;;
  *)
    echo "usage: $0 [run|status|collect]" >&2
    exit 1
    ;;
esac

echo "[totoro] ROOT=${ROOT}"
echo "[totoro] Phase B lane 3b (MNAR): N_JOBS=${N_JOBS}  NWORKERS=${NWORKERS}"
echo "[totoro] G6: MNAR is not MAR (g6_mar_label=FALSE for all cells)"
echo "[totoro] G8 wall ceiling: ${WALL_CEILING_SEC}s (~$(( WALL_CEILING_SEC / 3600 ))h)"
echo "[totoro] loadavg: $(cat /proc/loadavg)"
echo "[totoro] git HEAD: $(git -C "${ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "[totoro] started at $(date -Is)"

START_TS=$(date +%s)

seq 0 $((N_JOBS - 1)) | xargs -P "${NWORKERS}" -I{} bash -c '
  echo "[totoro] starting Phase B MNAR job {} at $(date -Is)"
  PIGAUTO_JOB_ID={} Rscript "'"${ROOT}"'/script/gnn_evidence_campaign_phase_b.R" \
    > "'"${ROOT}"'/logs/gnn_phase_b_mnar_job_{}.log" 2>&1
  echo "[totoro] finished Phase B MNAR job {} at $(date -Is)"
'

END_TS=$(date +%s)
WALL=$((END_TS - START_TS))
echo "[totoro] all Phase B MNAR jobs dispatched; wall_sec=${WALL}"

collect_cmd

FAIL=$(Rscript -e '
  df <- read.csv(file.path(Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results_phase_b_mnar"),
                           "gnn_phase_b_mnar_summary.csv"))
  cat(sum(!df$fit_ok))
' 2>/dev/null || echo "?")
FAIL_PCT=$(awk -v f="${FAIL}" -v n="${N_JOBS}" 'BEGIN { if (n>0) printf "%.1f", 100*f/n; else print "?" }')
echo "[totoro] G8 check: failures=${FAIL}/${N_JOBS} (${FAIL_PCT}%)  wall_sec=${WALL}"
if awk "BEGIN { exit !(${FAIL_PCT:-0} > 20) }" 2>/dev/null; then
  echo "[totoro] G8 STOP: >20% failures — pause and report" >&2
  exit 2
fi
if [[ "${WALL}" -gt "${WALL_CEILING_SEC}" ]]; then
  echo "[totoro] G8 STOP: wall > $(( WALL_CEILING_SEC / 3600 ))h — pause and report" >&2
  exit 3
fi

echo "[totoro] Phase B MNAR done at $(date -Is)"
