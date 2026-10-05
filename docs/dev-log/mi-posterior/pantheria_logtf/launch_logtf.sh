#!/bin/bash
# PanTHERIA log_transform = FALSE sensitivity (Shinichi 2026-10-05): 6 cells, pigauto dbe3304 (main 45011fc + opt-in MI_REALDATA_LOG_TRANSFORM), default residual_prior "sep".
# Mirrors script/mi_realdata/12_totoro_run.sh (same env, same 01_run.R call); only the cell list differs.
set -euo pipefail
CODE_DIR=$HOME/hsq_work/mi_real_sep/code_dbe3304
OUT_DIR=$HOME/hsq_work/mi_real_sep/out_logtf_dbe3304
mkdir -p "$OUT_DIR/logs"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=""
export MI_POST_SHA="$(cat "$CODE_DIR/SHA")" MI_REALDATA_OFFLINE=1 PIGAUTO_PKG_PATH="$CODE_DIR"
export R_LIBS="${MI_REALDATA_RLIB:-/home/snakagaw/R/lib}${R_LIBS:+:${R_LIBS}}"
unset MI_POST_NITER MI_POST_BURNIN; export MI_REALDATA_LOG_TRANSFORM=FALSE
export HARNESS="$CODE_DIR/script/mi_realdata" OUT_DIR
cd "$HARNESS"
Rscript -e '
suppressPackageStartupMessages(devtools::load_all(Sys.getenv("PIGAUTO_PKG_PATH"), quiet = TRUE))
for (p in c("phylolm", "ape", "Matrix")) if (!requireNamespace(p, quietly = TRUE)) stop("not installed: ", p)
stopifnot(identical(pigauto:::.mip_resolve_control(list(), m = 5L)$residual_prior, "sep"))
cat("PREFLIGHT_OK sep default\n")'
[ "${CONFIRM:-no}" = yes ] || { echo "dry run: set CONFIRM=yes to launch"; exit 0; }
Rscript -e 'source("lib.R"); p <- planned_cells(); p <- p[p$dataset == "pantheria", ]; cat(sprintf("%s:%s:%d\n", p$dataset, p$arm, p$seed), sep = "")' > "$OUT_DIR/logs/cells.txt"
run_cell() {
  local d="$1" a="$2" s="$3"; local name="${d}-${a}-m${s}"; local log="${OUT_DIR}/logs/${name}.log"
  echo "CELL_START ${name} $(date '+%F %T')"
  if nice -n 10 Rscript 01_run.R "${d}" "${a}" "${s}" "${OUT_DIR}" > "${log}" 2>&1 && grep -q '^MI_REALDATA_CELL_OK$' "${log}"; then
    echo "CELL_OK ${name} $(date '+%F %T')"; else echo "CELL_FAILED ${name} $(date '+%F %T') (see ${log})"; fi
}
export -f run_cell
setsid nohup bash -c '
  xargs -P 6 -a "$OUT_DIR/logs/cells.txt" -I{} bash -c "IFS=: read -r d a s <<< \"{}\"; run_cell \"\$d\" \"\$a\" \"\$s\"" >> "$OUT_DIR/logs/run.log" 2>&1
  echo "CELLS_DONE ok=$(grep -c ^CELL_OK "$OUT_DIR/logs/run.log" || true) failed=$(grep -c ^CELL_FAILED "$OUT_DIR/logs/run.log" || true) of 6 $(date "+%F %T")" >> "$OUT_DIR/logs/run.log"
' > /dev/null 2>&1 &
echo "launched pgid $(ps -o pgid= $! | tr -d ' ')"
