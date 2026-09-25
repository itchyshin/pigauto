#!/bin/bash
# mi-posterior real-data harness (G8, .unlazy/mi-posterior/GATES.md) on
# Totoro. PRIMARY route (2026-09-24); 10_fir.sbatch is the DRAC fallback.
#
# Usage (on Totoro):
#   bash 12_totoro_run.sh <code_dir> <out_dir> <n_jobs>
#
#   code_dir  a copy of the arc/mi-posterior worktree at a COMMITTED state,
#             containing:
#               SHA   one line, `git rev-parse HEAD` of that commit
#               script/mi_realdata/inputs/<cell>/mask_receipt.rds
#               script/mi_realdata/inputs/<cell>/conformal_metrics.rds
#             for all 10 planned cells. Totoro has no git history for
#             arc/mondrian-realdata, so 00_fetch_masks.R cannot run there,
#             and 01_run.R runs with MI_REALDATA_OFFLINE=1 (never calls git).
#   out_dir   receipts go to <out_dir>/<cell>/mi_posterior.rds, one log per
#             cell to <out_dir>/logs/<cell>.log, the run log to
#             <out_dir>/logs/run.log.
#   n_jobs    cells run at once, one core each (xargs -P). Must match
#             ^[1-9][0-9]*$ (no sign, no leading zero, so bash never reads
#             it as octal) and is refused above 150 (our Totoro cap, D-143).
#             Only 10 cells exist, so more than 10 adds nothing.
#
# Why one core per cell: multi_impute(draws_method = "posterior") has no
# n_cores control; one cell runs its chains one after another in one R
# process. BLAS/OpenMP threads are pinned to 1, and every R process runs
# under nice -n 10 (Totoro is shared).
#
# On the Mac, before copying (inputs/ is git-ignored, so rsync it
# explicitly; `git archive` would drop it):
#   cd <worktree>/script/mi_realdata && Rscript 00_fetch_masks.R --all --with-conformal
#   cd <worktree> && git rev-parse HEAD > SHA     # commit first: the SHA must name the code
#   rsync -a --exclude .git <worktree>/ snakagaw@totoro.biology.ualberta.ca:<code_dir>/
# Then, on Totoro:
#   bash <code_dir>/script/mi_realdata/12_totoro_run.sh <code_dir> <out_dir> 10
# Back on the Mac:
#   rsync -a snakagaw@totoro.biology.ualberta.ca:<out_dir>/ <worktree>/script/mi_realdata/returned/
#   Rscript script/mi_realdata/02_summarise.R && Rscript script/mi_realdata/03_acceptance.R
#
# pigauto is loaded with devtools::load_all(<code_dir>) (PIGAUTO_PKG_PATH;
# pigauto has no compiled code, so parallel load_all() has no build race).
# Its Imports (torch, ggplot2, ...) and phylolm/ape must be installed;
# MI_REALDATA_RLIB (default /home/snakagaw/R/lib) is put in front of the
# library path via R_LIBS (R skips it if the directory does not exist, and
# the user's own library stays on the path). A preflight checks all of this
# before any cell starts.
#
# Completed cells are skipped: a receipt with status "ok", the current
# schema, the same SHA and a sampler record with no override. Anything else
# (including the status "running" placeholder a killed run leaves) is
# rerun and overwritten. Paths may contain spaces: every Rscript call runs
# a script relative to the current directory or via -e, never by an
# absolute --file= path (Rscript would pass spaces there as `~+~`).
#
# Output lines: CELL_START / CELL_OK / CELL_SKIPPED / CELL_FAILED <cell>,
# then one final line:
#   CELLS_DONE ok=<n> skipped=<n> failed=<n> of 10
# The exit status is non-zero if any cell failed.

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: bash 12_totoro_run.sh <code_dir> <out_dir> <n_jobs>" >&2
  exit 2
fi
N_JOBS="$3"
if ! [[ "${N_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "REFUSED: n_jobs must be a positive integer without a sign or leading zero, got '${N_JOBS}'" >&2
  exit 2
fi
if (( N_JOBS > 150 )); then
  echo "REFUSED: n_jobs=${N_JOBS} exceeds the 150-core Totoro cap (D-143)" >&2
  exit 2
fi
CODE_DIR="$(cd "$1" && pwd)"
HARNESS="${CODE_DIR}/script/mi_realdata"
if [[ ! -s "${CODE_DIR}/SHA" ]]; then
  echo "REFUSED: ${CODE_DIR}/SHA missing or empty (write it on the Mac: git rev-parse HEAD > SHA)" >&2
  exit 2
fi
mkdir -p "$2/logs"
OUT_DIR="$(cd "$2" && pwd)"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MI_POST_SHA="$(cat "${CODE_DIR}/SHA")"
export MI_REALDATA_OFFLINE=1
export PIGAUTO_PKG_PATH="${CODE_DIR}"
export R_LIBS="${MI_REALDATA_RLIB:-/home/snakagaw/R/lib}${R_LIBS:+:${R_LIBS}}"
export CUDA_VISIBLE_DEVICES=""
unset MI_POST_NITER MI_POST_BURNIN   # the campaign runs at pigauto's sampler defaults
export HARNESS OUT_DIR

echo "[totoro] code_dir=${CODE_DIR} sha=${MI_POST_SHA}"
echo "[totoro] out_dir=${OUT_DIR} n_jobs=${N_JOBS} loadavg: $(cat /proc/loadavg 2>/dev/null || uptime)"

# ---- cell list, taken from lib.R::planned_cells() (one source of truth) ----
CELLS="$(Rscript -e 'source(file.path(Sys.getenv("HARNESS"), "lib.R")); p <- planned_cells(); cat(sprintf("%s:%s:%d\n", p$dataset, p$arm, p$seed), sep = "")')"
N_CELLS="$(printf '%s\n' "${CELLS}" | grep -c .)"

# ---- preflight 1: offline inputs for every cell ----
missing=0
while IFS=: read -r d a s; do
  for f in mask_receipt.rds conformal_metrics.rds; do
    if [[ ! -s "${HARNESS}/inputs/${d}-${a}-m${s}/${f}" ]]; then
      echo "MISSING_INPUT ${d}-${a}-m${s}/${f}"
      missing=1
    fi
  done
done <<< "${CELLS}"
if (( missing )); then
  echo "REFUSED: inputs missing. On the Mac run 'Rscript 00_fetch_masks.R --all --with-conformal' and copy script/mi_realdata/inputs/." >&2
  exit 2
fi

# ---- preflight 2: pigauto from code_dir, with the posterior method, and phylolm ----
Rscript -e '
suppressPackageStartupMessages(devtools::load_all(Sys.getenv("PIGAUTO_PKG_PATH"), quiet = TRUE))
for (p in c("phylolm", "ape", "Matrix")) if (!requireNamespace(p, quietly = TRUE)) stop("not installed: ", p)
stopifnot("posterior" %in% eval(formals(pigauto::multi_impute)$draws_method))
cat("PREFLIGHT_OK\n")'

# ---- one cell ----
run_cell() {
  local d="$1" a="$2" s="$3"
  local name="${d}-${a}-m${s}"
  local receipt="${OUT_DIR}/${name}/mi_posterior.rds"
  local log="${OUT_DIR}/logs/${name}.log"
  if [[ -s "${receipt}" ]] && Rscript -e '
      source(file.path(Sys.getenv("HARNESS"), "lib.R"))
      r <- readRDS(commandArgs(TRUE)[1])
      done <- identical(r$status, "ok") && identical(r$receipt_schema, receipt_schema_current) &&
        identical(r$code_sha, Sys.getenv("MI_POST_SHA")) &&
        !is.null(r$sampler) && "overrides" %in% names(r$sampler) && length(r$sampler$overrides) == 0L
      quit(status = if (done) 0L else 1L)' "${receipt}" 2>/dev/null; then
    echo "CELL_SKIPPED ${name} (ok receipt for this SHA exists)"
    return 0
  fi
  echo "CELL_START ${name} $(date '+%F %T')"
  # Relative script path on purpose (cwd is ${HARNESS}, set before xargs):
  # an absolute path with a space would reach R as `~+~` and break here().
  if nice -n 10 Rscript 01_run.R "${d}" "${a}" "${s}" "${OUT_DIR}" > "${log}" 2>&1 &&
      grep -q '^MI_REALDATA_CELL_OK$' "${log}"; then
    echo "CELL_OK ${name} $(date '+%F %T')"
  else
    echo "CELL_FAILED ${name} $(date '+%F %T') (see ${log})"
  fi
}
export -f run_cell

cd "${HARNESS}"
printf '%s\n' "${CELLS}" |
  xargs -P "${N_JOBS}" -I{} bash -c 'IFS=: read -r d a s <<< "{}"; run_cell "$d" "$a" "$s"' |
  tee "${OUT_DIR}/logs/run.log"

n_ok="$(grep -c '^CELL_OK ' "${OUT_DIR}/logs/run.log" || true)"
n_skip="$(grep -c '^CELL_SKIPPED ' "${OUT_DIR}/logs/run.log" || true)"
n_fail="$(grep -c '^CELL_FAILED ' "${OUT_DIR}/logs/run.log" || true)"
echo "CELLS_DONE ok=${n_ok} skipped=${n_skip} failed=${n_fail} of ${N_CELLS}" | tee -a "${OUT_DIR}/logs/run.log"
(( n_fail == 0 && n_ok + n_skip == N_CELLS ))
