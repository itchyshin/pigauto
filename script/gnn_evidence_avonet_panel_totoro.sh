#!/bin/bash
# Totoro CPU launcher for AVONET300 multi-seed validation panel.
# 15 jobs = 5 seeds × 3 methods (pigauto_fixed1, pigauto_bayes, bace).
# Usage on Totoro after rsync:
#   cd ~/pigauto_gnn_evidence_campaign
#   nohup bash script/gnn_evidence_avonet_panel_totoro.sh > logs/avonet_panel.log 2>&1 &
# Poll:
#   bash script/gnn_evidence_avonet_panel_totoro.sh status
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"
mkdir -p logs results_avonet_panel

N_JOBS=15
NWORKERS="${NWORKERS:-15}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PIGAUTO_PKG_PATH="${ROOT}"
export PIGAUTO_OUT_DIR="${ROOT}/results_avonet_panel"
export PIGAUTO_SKIP_EXISTING="${PIGAUTO_SKIP_EXISTING:-1}"
export CUDA_VISIBLE_DEVICES=""

status_cmd() {
  local done fail
  done=$(find "${ROOT}/results_avonet_panel" -maxdepth 1 -name 'avonet_panel_job_*.rds' 2>/dev/null | wc -l | tr -d ' ')
  fail=0
  if [[ "${done}" -gt 0 ]]; then
    fail=$(PIGAUTO_OUT_DIR="${ROOT}/results_avonet_panel" Rscript -e '
      fs <- sort(Sys.glob(file.path(Sys.getenv("PIGAUTO_OUT_DIR"), "avonet_panel_job_*.rds")))
      rows <- lapply(fs, readRDS)
      cat(sum(!vapply(rows, function(r) isTRUE(r$fit_ok), logical(1))))
    ' 2>/dev/null || echo 0)
  fi
  echo "[totoro-avonet-panel] progress: ${done}/${N_JOBS} RDS  failures=${fail}"
  echo "[totoro-avonet-panel] seeds=2026-2030  miss=30%  methods=fixed_1,bayes,bace"
  Rscript -e 'cat("Rphylopars:", requireNamespace("Rphylopars", quietly=TRUE),
                " BACE:", requireNamespace("BACE", quietly=TRUE), "\n")' 2>/dev/null || true
  if [[ -f "${ROOT}/logs/avonet_panel.log" ]]; then
    echo "[totoro-avonet-panel] log tail:"
    tail -n 5 "${ROOT}/logs/avonet_panel.log" || true
  fi
}

collect_cmd() {
  echo "[totoro] collecting AVONET panel summary"
  PIGAUTO_AVONET_DIR="${ROOT}/results_avonet_panel" \
    PIGAUTO_OUT_MD="${ROOT}/script/returned_gnn_campaign/AVONET_PANEL_SUMMARY.md" \
    Rscript "${ROOT}/script/collect_gnn_evidence_avonet_panel.R"
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
echo "[totoro] N_JOBS=${N_JOBS}  NWORKERS=${NWORKERS}"
echo "[totoro] loadavg: $(cat /proc/loadavg)"
echo "[totoro] git HEAD: $(git -C "${ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "[totoro] started at $(date -Is)"

START_TS=$(date +%s)

seq 0 $((N_JOBS - 1)) | xargs -P "${NWORKERS}" -I{} bash -c '
  echo "[totoro] starting AVONET job {} at $(date -Is)"
  PIGAUTO_JOB_ID={} Rscript "'"${ROOT}"'/script/gnn_evidence_avonet_panel.R" \
    > "'"${ROOT}"'/logs/avonet_panel_job_{}.log" 2>&1
  echo "[totoro] finished AVONET job {} at $(date -Is)"
'

END_TS=$(date +%s)
WALL=$((END_TS - START_TS))
echo "[totoro] all AVONET jobs dispatched; wall_sec=${WALL}"

collect_cmd

echo "[totoro] done at $(date -Is)"
