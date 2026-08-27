#!/bin/bash
# Totoro CPU launcher for GNN evidence sentinel pre-run (12 fits).
# Usage on Totoro after rsync:
#   cd ~/pigauto_gnn_sentinel_prerun
#   bash script/run_gnn_evidence_sentinel_prerun_totoro.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"
mkdir -p logs results

NWORKERS="${NWORKERS:-12}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PIGAUTO_PKG_PATH="${ROOT}"
export PIGAUTO_OUT_DIR="${ROOT}/results"
export CUDA_VISIBLE_DEVICES=""

echo "[totoro] ROOT=${ROOT}"
echo "[totoro] NWORKERS=${NWORKERS}"
echo "[totoro] loadavg: $(cat /proc/loadavg)"
echo "[totoro] git HEAD: $(git -C "${ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"

seq 0 11 | xargs -P "${NWORKERS}" -I{} bash -c '
  echo "[totoro] starting job {} at $(date -Is)"
  PIGAUTO_JOB_ID={} Rscript "'"${ROOT}"'/script/gnn_evidence_sentinel_prerun.R" \
    > "'"${ROOT}"'/logs/gnn_sentinel_job_{}.log" 2>&1
  echo "[totoro] finished job {} at $(date -Is)"
'

echo "[totoro] all jobs dispatched; collecting summary"
Rscript -e '
  out <- Sys.getenv("PIGAUTO_OUT_DIR", "'"${ROOT}"'/results")
  fs <- sort(Sys.glob(file.path(out, "gnn_sentinel_job_*.rds")))
  stopifnot(length(fs) > 0)
  rows <- lapply(fs, readRDS)
  df <- do.call(rbind, lapply(rows, function(r) {
    data.frame(
      job_id=r$job_id, family=r$family, n_species=r$n_species, seed=r$seed,
      fit_ok=r$fit_ok, fit_sec=r$fit_sec,
      blend_loss=r$blend_loss, baseline_loss=r$baseline_loss,
      paired_delta=r$paired_delta,
      r_cal_gnn_mean=if (is.null(r$r_cal_gnn_mean)) NA else r$r_cal_gnn_mean,
      r_cal_bm_mean=if (is.null(r$r_cal_bm_mean)) NA else r$r_cal_bm_mean,
      r_cal_mean_mean=if (is.null(r$r_cal_mean_mean)) NA else r$r_cal_mean_mean,
      floor_fired=if (is.null(r$floor_fired)) NA else r$floor_fired,
      error=if (is.null(r$error)) NA else r$error,
      stringsAsFactors=FALSE
    )
  }))
  write.csv(df, file.path(out, "gnn_sentinel_prerun_summary.csv"), row.names=FALSE)
  cat("\n=== timing summary ===\n")
  print(aggregate(fit_sec ~ family + n_species, df, mean))
  cat("\ntotal wall (max fit_sec):", max(df$fit_sec, na.rm=TRUE), "sec\n")
  cat("failures:", sum(!df$fit_ok), "/", nrow(df), "\n")
'

echo "[totoro] done at $(date -Is)"
