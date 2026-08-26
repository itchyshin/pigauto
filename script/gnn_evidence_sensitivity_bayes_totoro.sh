#!/bin/bash
# Totoro launcher: bayes sensitivity arm (λ_DGP ∈ {0.2, 0.5} only).
# Prereg §4: separate sensitivity arm, not primary estimand.
# Output: results_bayes/ (does not overwrite primary results/)
#
# Usage:
#   nohup bash script/gnn_evidence_sensitivity_bayes_totoro.sh > logs/sensitivity_bayes.log 2>&1 &
#   bash script/gnn_evidence_sensitivity_bayes_totoro.sh status
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"
mkdir -p logs results_bayes

NWORKERS="${NWORKERS:-100}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PIGAUTO_PKG_PATH="${ROOT}"
export PIGAUTO_OUT_DIR="${ROOT}/results_bayes"
export PIGAUTO_LAMBDA_MODE="bayes"
export PIGAUTO_SKIP_EXISTING="${PIGAUTO_SKIP_EXISTING:-1}"
export CUDA_VISIBLE_DEVICES=""

JOB_IDS_FILE="${ROOT}/logs/bayes_low_lambda_job_ids.txt"
if [[ ! -f "${JOB_IDS_FILE}" ]]; then
  Rscript -e '
    cells <- expand.grid(
      family = c("F1", "F2", "F3"),
      n_species = c(100L, 300L, 1000L),
      phylo_lambda = c(0.2, 0.5, 1.0),
      miss_frac = c(0.10, 0.30, 0.50),
      stringsAsFactors = FALSE
    )
    cells$cell_id <- seq_len(nrow(cells)) - 1L
    reps <- data.frame(rep = seq_len(30L))
    jobs <- merge(cells, reps, by = NULL)
    jobs$job_id <- seq_len(nrow(jobs)) - 1L
    low <- jobs[jobs$phylo_lambda %in% c(0.2, 0.5), "job_id"]
    writeLines(as.character(low), "'"${JOB_IDS_FILE}"'")
    cat("wrote", length(low), "job ids\n")
  '
fi

N_JOBS=$(wc -l < "${JOB_IDS_FILE}" | tr -d ' ')

status_cmd() {
  local done fail
  done=$(find "${ROOT}/results_bayes" -maxdepth 1 -name 'gnn_campaign_job_*.rds' 2>/dev/null | wc -l | tr -d ' ')
  fail=0
  if [[ "${done}" -gt 0 ]]; then
    fail=$(Rscript -e '
      fs <- sort(Sys.glob(file.path(Sys.getenv("PIGAUTO_OUT_DIR"), "gnn_campaign_job_*.rds")))
      rows <- lapply(fs, readRDS)
      cat(sum(!vapply(rows, function(r) isTRUE(r$fit_ok), logical(1))))
    ' 2>/dev/null || echo 0)
  fi
  echo "[totoro-bayes] progress: ${done}/${N_JOBS} RDS  failures=${fail}  lambda_mode=bayes"
  if [[ -f "${ROOT}/logs/sensitivity_bayes.log" ]]; then
    echo "[totoro-bayes] log tail:"
    tail -n 3 "${ROOT}/logs/sensitivity_bayes.log" || true
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

echo "[totoro-bayes] ROOT=${ROOT}  jobs=${N_JOBS}  workers=${NWORKERS}"
echo "[totoro-bayes] started at $(date -Is)"
START_TS=$(date +%s)

cat "${JOB_IDS_FILE}" | xargs -P "${NWORKERS}" -I{} bash -c '
  jid={}
  out="'"${ROOT}"'/results_bayes/gnn_campaign_job_$(printf "%04d" ${jid}).rds"
  if [[ "'"${PIGAUTO_SKIP_EXISTING}"'" == "1" && -f "${out}" ]]; then exit 0; fi
  PIGAUTO_JOB_ID=${jid} Rscript "'"${ROOT}"'/script/gnn_evidence_campaign.R" \
    > "'"${ROOT}"'/logs/gnn_campaign_bayes_job_${jid}.log" 2>&1
'

END_TS=$(date +%s)
WALL=$((END_TS - START_TS))
echo "[totoro-bayes] done wall_sec=${WALL} at $(date -Is)"
