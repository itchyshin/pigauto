#!/bin/bash
# Totoro launcher: F2 @ lambda_DGP=1 confirmatory arm (prereg §3.4).
# 5 cells × 60 seeds = 300 fits; lambda_mode=fixed_1.
#
# Usage:
#   nohup bash script/gnn_evidence_f2_confirm_totoro.sh > logs/f2_confirm.log 2>&1 &
#   bash script/gnn_evidence_f2_confirm_totoro.sh status
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"
mkdir -p logs results_confirm

N_JOBS=300
NWORKERS="${NWORKERS:-100}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PIGAUTO_PKG_PATH="${ROOT}"
export PIGAUTO_OUT_DIR="${ROOT}/results_confirm"
export PIGAUTO_SKIP_EXISTING="${PIGAUTO_SKIP_EXISTING:-1}"
export CUDA_VISIBLE_DEVICES=""

status_cmd() {
  local done fail
  done=$(find "${ROOT}/results_confirm" -maxdepth 1 -name 'gnn_confirm_job_*.rds' 2>/dev/null | wc -l | tr -d ' ')
  fail=0
  if [[ "${done}" -gt 0 ]]; then
    fail=$(Rscript -e '
      fs <- sort(Sys.glob(file.path(Sys.getenv("PIGAUTO_OUT_DIR"), "gnn_confirm_job_*.rds")))
      rows <- lapply(fs, readRDS)
      cat(sum(!vapply(rows, function(r) isTRUE(r$fit_ok), logical(1))))
    ' 2>/dev/null || echo 0)
  fi
  echo "[totoro-f2-confirm] progress: ${done}/${N_JOBS} RDS  failures=${fail}  lambda_mode=fixed_1"
  if [[ -f "${ROOT}/logs/f2_confirm.log" ]]; then
    echo "[totoro-f2-confirm] log tail:"
    tail -n 3 "${ROOT}/logs/f2_confirm.log" || true
  fi
}

case "${1:-run}" in
  status) status_cmd; exit 0 ;;
  run) ;;
  *)
    echo "usage: $0 [run|status]" >&2
    exit 1
    ;;
esac

echo "[totoro-f2-confirm] ROOT=${ROOT}  jobs=${N_JOBS}  workers=${NWORKERS}"
echo "[totoro-f2-confirm] cells: F2 lambda=1 n=300 miss 10/30/50 + n=1000 miss 10/30"
echo "[totoro-f2-confirm] started at $(date -Is)"
START_TS=$(date +%s)

seq 0 $((N_JOBS - 1)) | xargs -P "${NWORKERS}" -I{} bash -c '
  jid={}
  out="'"${ROOT}"'/results_confirm/gnn_confirm_job_$(printf "%04d" ${jid}).rds"
  if [[ "'"${PIGAUTO_SKIP_EXISTING}"'" == "1" && -f "${out}" ]]; then exit 0; fi
  PIGAUTO_JOB_ID=${jid} Rscript "'"${ROOT}"'/script/gnn_evidence_f2_confirm.R" \
    > "'"${ROOT}"'/logs/gnn_confirm_job_${jid}.log" 2>&1
'

END_TS=$(date +%s)
WALL=$((END_TS - START_TS))
echo "[totoro-f2-confirm] done wall_sec=${WALL} at $(date -Is)"
