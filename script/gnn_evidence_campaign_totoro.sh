#!/bin/bash
# Totoro CPU launcher for GNN evidence Phase A campaign (2,430 fits).
# Usage on Totoro after rsync:
#   cd ~/pigauto_gnn_evidence_campaign
#   nohup bash script/gnn_evidence_campaign_totoro.sh > logs/campaign.log 2>&1 &
# Poll:
#   tail -f logs/campaign.log
#   bash script/gnn_evidence_campaign_totoro.sh status
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"
mkdir -p logs results

N_JOBS=2430
NWORKERS="${NWORKERS:-100}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PIGAUTO_PKG_PATH="${ROOT}"
export PIGAUTO_OUT_DIR="${ROOT}/results"
export PIGAUTO_LAMBDA_MODE="${PIGAUTO_LAMBDA_MODE:-fixed_1}"
export PIGAUTO_SKIP_EXISTING="${PIGAUTO_SKIP_EXISTING:-1}"
export CUDA_VISIBLE_DEVICES=""

status_cmd() {
  local done fail total
  done=$(find "${ROOT}/results" -maxdepth 1 -name 'gnn_campaign_job_*.rds' 2>/dev/null | wc -l | tr -d ' ')
  total="${N_JOBS}"
  fail=0
  if [[ -d "${ROOT}/results" ]]; then
    fail=$(Rscript -e '
      fs <- sort(Sys.glob(file.path(Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results"),
                                   "gnn_campaign_job_*.rds")))
      if (!length(fs)) { cat(0); quit(save="no") }
      rows <- lapply(fs, readRDS)
      cat(sum(!vapply(rows, function(r) isTRUE(r$fit_ok), logical(1))))
    ' 2>/dev/null || echo 0)
  fi
  echo "[totoro] progress: ${done}/${total} RDS  failures=${fail}"
  if [[ -f "${ROOT}/logs/campaign.log" ]]; then
    echo "[totoro] campaign.log tail:"
    tail -n 5 "${ROOT}/logs/campaign.log" || true
  fi
  if [[ -f "${ROOT}/results/gnn_campaign_summary.csv" ]]; then
    echo "[totoro] summary present: results/gnn_campaign_summary.csv"
  fi
}

collect_cmd() {
  echo "[totoro] collecting summary"
  Rscript -e '
    out <- Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results")
    fs <- sort(Sys.glob(file.path(out, "gnn_campaign_job_*.rds")))
    stopifnot(length(fs) > 0)
    rows <- lapply(fs, readRDS)
    df <- do.call(rbind, lapply(rows, function(r) {
      data.frame(
        job_id=r$job_id, cell_id=r$cell_id, rep=r$rep,
        family=r$family, n_species=r$n_species,
        phylo_lambda=r$phylo_lambda, miss_frac=r$miss_frac,
        seed=r$seed, lambda_mode=r$lambda_mode,
        fit_ok=r$fit_ok, fit_sec=r$fit_sec,
        blend_loss=r$blend_loss, baseline_loss=r$baseline_loss,
        paired_delta=r$paired_delta,
        r_cal_gnn_mean=r$r_cal_gnn_mean,
        r_cal_bm_mean=r$r_cal_bm_mean,
        r_cal_mean_mean=r$r_cal_mean_mean,
        floor_fired=r$floor_fired,
        error=if (is.null(r$error)) NA else r$error,
        stringsAsFactors=FALSE
      )
    }))
    write.csv(df, file.path(out, "gnn_campaign_summary.csv"), row.names=FALSE)
    saveRDS(df, file.path(out, "gnn_campaign_summary.rds"))
    cat("\n=== campaign summary ===\n")
    cat("fits:", nrow(df), "\n")
    cat("failures:", sum(!df$fit_ok), "\n")
    cat("max fit_sec:", max(df$fit_sec, na.rm=TRUE), "\n")
    by_cell <- aggregate(
      cbind(fit_ok = as.integer(df$fit_ok), paired_delta = df$paired_delta) ~
        family + n_species + phylo_lambda + miss_frac,
      df, function(x) c(n=length(x), fail=sum(!df$fit_ok[seq_along(x)]))
    )
    print(head(by_cell, 3))
  '
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
echo "[totoro] N_JOBS=${N_JOBS}  NWORKERS=${NWORKERS}  lambda_mode=${PIGAUTO_LAMBDA_MODE}"
echo "[totoro] loadavg: $(cat /proc/loadavg)"
echo "[totoro] git HEAD: $(git -C "${ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "[totoro] started at $(date -Is)"

START_TS=$(date +%s)

seq 0 $((N_JOBS - 1)) | xargs -P "${NWORKERS}" -I{} bash -c '
  echo "[totoro] starting job {} at $(date -Is)"
  PIGAUTO_JOB_ID={} Rscript "'"${ROOT}"'/script/gnn_evidence_campaign.R" \
    > "'"${ROOT}"'/logs/gnn_campaign_job_{}.log" 2>&1
  echo "[totoro] finished job {} at $(date -Is)"
'

END_TS=$(date +%s)
WALL=$((END_TS - START_TS))
echo "[totoro] all jobs dispatched; wall_sec=${WALL}"

collect_cmd

FAIL=$(Rscript -e '
  df <- read.csv(file.path(Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results"),
                           "gnn_campaign_summary.csv"))
  cat(sum(!df$fit_ok))
' 2>/dev/null || echo "?")
FAIL_PCT=$(awk -v f="${FAIL}" -v n="${N_JOBS}" 'BEGIN { if (n>0) printf "%.1f", 100*f/n; else print "?" }')
echo "[totoro] G8 check: failures=${FAIL}/${N_JOBS} (${FAIL_PCT}%)  wall_sec=${WALL}"
if awk "BEGIN { exit !(${FAIL_PCT:-0} > 20) }" 2>/dev/null; then
  echo "[totoro] G8 STOP: >20% failures — pause and report" >&2
  exit 2
fi
if [[ "${WALL}" -gt 7200 ]]; then
  echo "[totoro] G8 STOP: wall > 2h — pause and report" >&2
  exit 3
fi

echo "[totoro] done at $(date -Is)"
