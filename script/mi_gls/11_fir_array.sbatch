#!/bin/bash
# PRIMARY ROUTE IS TOTORO: script/mi_gls/12_totoro_campaign.sh. This fir (DRAC) array is the fallback.
#SBATCH --job-name=pigauto_mi_posterior
#SBATCH --account=def-snakagaw
#SBATCH --time=SET_ME           # REQUIRED; sbatch rejects this token on purpose. See "Sizing" below.
#SBATCH --mem=4G                # finite starting value; refine with `seff <jobid>` after a one-task pre-run
#SBATCH --cpus-per-task=1
#SBATCH --array=0-47
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=itchyshin@gmail.com
#SBATCH --output=logs/mi_posterior-%A_%a.out

# pigauto: posterior-MI validation sweep (arc/mi-posterior), fir CPU array.
# 24 regimes (script/mi_gls/regimes.R) x 200 reps = 4800 cells, split into
# 48 tasks of 100 reps each. Task t (0-indexed) handles global cell indices
# [t*100, t*100+99] over a row-major (regime_id 1..24) x (rep 1..200) grid,
# so every task holds 100 reps of ONE regime.
#
# Sizing (review cluster#3, cluster#4, estimands#4). Size --time and --mem
# from a timed REAL cell, not from gate_convergence.R (one posterior_full
# fit at m = 5 on n = 300 only). The n = 1000 regimes (3, 4, 7, 8, 11, 12,
# 15, 16, 18, 20, 22, 24) dominate, and the both-missing ones cost most:
# posterior_none (both-missing regimes only) comes from the same chain run
# (param_uncertainty = "both") but adds m more completed datasets and
# their downstream fits, and the conformal arm adds one impute(gnn = FALSE)
# fit per cell.
# On a compute node:
#   salloc --account=def-snakagaw --time=1:00:00 --mem=4G --cpus-per-task=1
#   cd "$CODE_DIR" && OPENBLAS_NUM_THREADS=1 /usr/bin/time -v \
#     Rscript script/mi_gls/01_cell_v2.R 24 1 "$SLURM_TMPDIR/sz"
# Then --time = that wall time x 100 reps x 1.5, and --mem = 1.5 x its
# "Maximum resident set size". Or submit --array=0 first, read `seff`, and
# only then submit 1-47. Edit the --time line in the copy you submit
# (bootstrap step 4): sbatch reads #SBATCH lines before the command line,
# so an invalid SET_ME here fails the submission even with --time given.
#
# One-time bootstrap on the fir login node:
#   1. Freeze the code (one git SHA; /project quota is full as of 2026-09,
#      so $HOME, never /scratch):
#        SHA=$(git -C <worktree> rev-parse HEAD)      # commit first
#        CODE_DIR=$HOME/mi_post_code/$SHA
#        mkdir -p "$CODE_DIR" && git -C <worktree> archive "$SHA" | tar -x -C "$CODE_DIR"
#        echo "$SHA" > "$CODE_DIR/SHA"
#   2. Check the R deps on the login node (this array never installs):
#        module load StdEnv/2023 gcc/12.3 r/4.4.0
#        cd "$CODE_DIR" && Rscript -e 'suppressMessages(devtools::load_all(quiet = TRUE)); for (p in c("phylolm", "ape", "nlme")) if (!requireNamespace(p, quietly = TRUE)) message("MISSING ", p)'
#      Install anything missing HERE (login node), into R_LIBS_USER.
#   3. Make the output and log directory BEFORE sbatch (slurm opens
#      --output=logs/... relative to the submit directory before any line
#      of this script runs, and does not create directories):
#        OUT_DIR=$SCRATCH/mi_posterior/$SHA
#        mkdir -p "$OUT_DIR/logs"
#   4. Submit FROM $OUT_DIR, so logs/ resolves to $OUT_DIR/logs and
#      results land in $OUT_DIR:
#        cp "$CODE_DIR/script/mi_gls/11_fir_array.sbatch" "$OUT_DIR/" && (edit --time, --mem there)
#        cd "$OUT_DIR" && sbatch --export=ALL,CODE_DIR="$CODE_DIR" 11_fir_array.sbatch
#   Resubmitting the same command resumes: cells whose .rds exists are skipped.
#   5. Summarise into the paths the GATES.md G6/G7 CHECK lines read, copy
#      them back into the worktree, and gate there (as in 12_totoro_campaign.sh):
#        cd "$CODE_DIR" && D=docs/dev-log/mi-posterior && Rscript script/mi_gls/03_summarise_v2.R \
#          "$OUT_DIR" $D/sim_summary.csv $D/sim_summary.md $D/cell_coverage.csv
#        (Mac, worktree root) rsync the three files into docs/dev-log/mi-posterior/,
#        then run the G6 and G7 CHECK lines.

set -euo pipefail
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-}"
echo "SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID:-}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-}"

module load StdEnv/2023 gcc/12.3 cuda/12.6 r/4.4.0
# LD_LIBRARY_PATH is unset on fir after these module loads; guard it for set -u.
export LD_LIBRARY_PATH="$EBROOTCUDA/targets/x86_64-linux/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

: "${CODE_DIR:?set CODE_DIR to the frozen code directory (bootstrap step 1), e.g. sbatch --export=ALL,CODE_DIR=...}"
if [[ ! -f "${CODE_DIR}/SHA" || ! -f "${CODE_DIR}/script/mi_gls/01_cell_v2.R" ]]; then
  echo "CODE_DIR=${CODE_DIR} is not a frozen code directory (needs SHA and script/mi_gls/01_cell_v2.R)" >&2
  exit 1
fi
# A git checkout (.git directory or worktree file) is refused: load_all()
# would also source untracked or ignored files that no SHA records.
if [[ -e "${CODE_DIR}/.git" ]]; then
  echo "CODE_DIR=${CODE_DIR} is a git checkout; unpack a git archive instead (bootstrap step 1)" >&2
  exit 1
fi
export MI_POST_SHA="$(tr -d '[:space:]' < "${CODE_DIR}/SHA")"

OUT_DIR="${OUT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
mkdir -p "${OUT_DIR}"
# One code version per output directory; a resumed run must use the same SHA.
[[ -f "${OUT_DIR}/CODE_SHA" ]] || echo "${MI_POST_SHA}" > "${OUT_DIR}/CODE_SHA"
if [[ "$(tr -d '[:space:]' < "${OUT_DIR}/CODE_SHA")" != "${MI_POST_SHA}" ]]; then
  echo "OUT_DIR ${OUT_DIR} holds results from code $(cat "${OUT_DIR}/CODE_SHA"), not ${MI_POST_SHA}; use a fresh OUT_DIR" >&2
  exit 1
fi

# ---- R deps: checked across ALL library paths, never installed here -------
# (review cluster#1: nlme lives only in the r/4.4.0 system library, so a
# home-library-only check always failed; a compute node must never install.)
# (No backslashes in this R code: Rscript -e rewrites them.)
if ! (cd "${CODE_DIR}" && Rscript -e '
    imp <- read.dcf("DESCRIPTION", fields = "Imports")[1, 1]
    imp <- trimws(sub("[(].*$", "", strsplit(gsub("[[:space:]]+", " ", imp), ",")[[1]]))
    need <- unique(c("phylolm", "ape", "nlme", "devtools", imp[nzchar(imp)]))
    miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss)) {
      message("MISSING R packages: ", paste(miss, collapse = " "))
      message("libPaths: ", paste(.libPaths(), collapse = " "))
      quit(status = 1)
    }'); then
  echo "missing R deps: install them on the LOGIN node (bootstrap step 2), then resubmit" >&2
  exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
N_REGIMES=24
N_REPS=200
REPS_PER_TASK=100
START=$(( TASK_ID * REPS_PER_TASK ))
END=$(( START + REPS_PER_TASK - 1 ))

echo "Task ${TASK_ID}: global cell indices ${START}-${END} (of $((N_REGIMES * N_REPS)))"
echo "Code: ${CODE_DIR} (SHA ${MI_POST_SHA}); output: ${OUT_DIR}"
echo "R version: $(R --version | head -1)"

# A failing cell prints CELL_FAILED and the loop moves on (review
# cluster#6); the task exits non-zero at the end if any cell failed, so
# slurm still reports it.
RAN=0; SKIPPED=0; FAILED=0
for GLOBAL in $(seq "${START}" "${END}"); do
  REGIME=$(( GLOBAL / N_REPS + 1 ))
  REP=$(( GLOBAL % N_REPS + 1 ))
  if [[ -f "${OUT_DIR}/regime_${REGIME}_rep_${REP}.rds" ]]; then
    SKIPPED=$(( SKIPPED + 1 ))
    continue
  fi
  echo "--- regime=${REGIME} rep=${REP} ---"
  rc=0
  (cd "${CODE_DIR}" && Rscript script/mi_gls/01_cell_v2.R "${REGIME}" "${REP}" "${OUT_DIR}") || rc=$?
  RAN=$(( RAN + 1 ))
  if (( rc != 0 )); then
    echo "CELL_FAILED regime=${REGIME} rep=${REP} status=${rc}"
    FAILED=$(( FAILED + 1 ))
  fi
done

echo "TASK_DONE task=${TASK_ID} ran=${RAN} skipped_existing=${SKIPPED} failed=${FAILED} at $(date)"
if (( FAILED > 0 )); then exit 1; fi
