#!/bin/bash
#SBATCH --job-name=pigauto_mi_posterior
#SBATCH --account=def-snakagaw
#SBATCH --time=00:00:00        # PLACEHOLDER -- set from script/mi_gls/gate_convergence.R's
                                # per-cell wall time x 100 reps/task before submitting
#SBATCH --mem=0G                # PLACEHOLDER -- set from `seff <jobid>` after the first real run
#SBATCH --cpus-per-task=1
#SBATCH --array=0-47
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=itchyshin@gmail.com
#SBATCH --output=logs/mi_posterior-%A_%a.out

# pigauto -- posterior-MI validation sweep (arc/mi-posterior), Fir CPU array.
# 24 regimes (script/mi_gls/regimes.R) x 200 reps = 4800 cells, split into
# 48 tasks of 100 reps each. Task t (0-indexed) handles global cell indices
# [t*100, t*100+99] over a row-major (regime_id 1..24) x (rep 1..200) grid.
#
# --time and --mem are PLACEHOLDERS (see comments above): right-size them
# from script/mi_gls/gate_convergence.R's wall time for ONE posterior_full
# fit on regime 1 and regime 23 before the first real submission, then
# refine with `seff <jobid>` after it (see ~/shinichi-brain/tools/
# drac-setup.md "golden rule" 5).
#
# One-time bootstrap on the login node before the first submission:
#   1. Build the source tarball this array installs per-task:
#        cd <worktree root> && R CMD build --no-build-vignettes .
#        mv pigauto_*.tar.gz script/mi_gls/pigauto_mi_posterior_src.tar.gz
#   2. sbatch script/mi_gls/11_fir_array.sbatch
#
# The persistent-CRAN-deps install (phylolm et al., below) is serialised
# across the 48 concurrent array tasks with flock, so no separate install
# job is required -- but it only needs outbound internet once; if the
# compute node has none, `module load httpproxy` (below) is the DRAC
# escape hatch (see ~/shinichi-brain/tools/drac-setup.md).

set -euo pipefail
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-}"
echo "SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID:-}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-}"

module load StdEnv/2023 gcc/12.3 cuda/12.6 r/4.4.0
export LD_LIBRARY_PATH="$EBROOTCUDA/targets/x86_64-linux/lib:$LD_LIBRARY_PATH"

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${ROOT}"
mkdir -p logs

# ---- persistent CRAN deps (phylolm, ape, nlme, devtools), installed ONCE --
# into a shared, PERSISTENT library on /project (never /scratch -- purged
# ~60 d, no backup; ~/shinichi-brain/tools/drac-setup.md rule 2). 48 array
# tasks would otherwise race on `install.packages()`/`R CMD INSTALL` to the
# same directory; flock serialises them onto one winner, everyone else
# blocks on the lock and then finds the deps already present.
# fir: /project quota is full (2026-09); the home library already holds phylolm 2.6.5.
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
mkdir -p "${R_LIBS_USER}"
LOCKFILE="${R_LIBS_USER}/.pigauto_mi_posterior_install.lock"
(
  flock -x 200
  if Rscript --vanilla -e '
      lib <- Sys.getenv("R_LIBS_USER")
      need <- c("phylolm", "ape", "nlme", "devtools")
      have <- rownames(installed.packages(lib.loc = lib))
      quit(status = as.integer(length(setdiff(need, have)) > 0L))
    ' 2>/dev/null; then
    echo "[install] all persistent CRAN deps already present in ${R_LIBS_USER}"
  else
    echo "[install] installing missing CRAN deps into ${R_LIBS_USER}"
    module load httpproxy || true   # DRAC compute nodes may have no direct internet
    Rscript --vanilla -e '
      options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
      lib <- Sys.getenv("R_LIBS_USER")
      dir.create(lib, recursive = TRUE, showWarnings = FALSE)
      .libPaths(c(lib, .libPaths()))
      need <- c("phylolm", "ape", "nlme", "devtools")
      have <- rownames(installed.packages(lib.loc = lib))
      miss <- setdiff(need, have)
      if (length(miss)) install.packages(miss, lib = lib, Ncpus = 4)
    '
  fi
) 200>"${LOCKFILE}"

# ---- pigauto itself, from a frozen source tarball, into THIS TASK's ------
# node-local $SLURM_TMPDIR. Never a shared directory: 48 concurrent
# `R CMD INSTALL`s of the package under test would race; $SLURM_TMPDIR is
# private to this job step and wiped at job end, which is fine since every
# task rebuilds from the identical tarball.
PKG_TARBALL="${ROOT}/script/mi_gls/pigauto_mi_posterior_src.tar.gz"
if [[ ! -f "${PKG_TARBALL}" ]]; then
  echo "missing ${PKG_TARBALL} -- build it on the login node first:" >&2
  echo "  R CMD build --no-build-vignettes ${ROOT}" >&2
  echo "  mv pigauto_*.tar.gz ${PKG_TARBALL}" >&2
  exit 1
fi
TASK_LIB="${SLURM_TMPDIR}/Rlib"
mkdir -p "${TASK_LIB}"
R CMD INSTALL --no-docs --no-byte-compile --library="${TASK_LIB}" "${PKG_TARBALL}"
export R_LIBS="${TASK_LIB}:${R_LIBS_USER}"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

OUTDIR="${SCRATCH:-${ROOT}/results}/mi_posterior/regime_reps"   # results to /scratch
mkdir -p "${OUTDIR}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
N_REGIMES=24
N_REPS=200
REPS_PER_TASK=100
START=$(( TASK_ID * REPS_PER_TASK ))
END=$(( START + REPS_PER_TASK - 1 ))

echo "Task ${TASK_ID}: global cell indices ${START}-${END} (of $((N_REGIMES * N_REPS)))"
echo "R version: $(R --version | head -1)"

for GLOBAL in $(seq "${START}" "${END}"); do
  REGIME=$(( GLOBAL / N_REPS + 1 ))
  REP=$(( GLOBAL % N_REPS + 1 ))
  echo "--- regime=${REGIME} rep=${REP} ---"
  Rscript "${ROOT}/script/mi_gls/01_cell_v2.R" "${REGIME}" "${REP}" "${OUTDIR}"
done

echo "Task ${TASK_ID} finished at: $(date)"
