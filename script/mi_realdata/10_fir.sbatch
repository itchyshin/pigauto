#!/bin/bash
#SBATCH --job-name=pigauto_mi_realdata
#SBATCH --account=def-snakagaw               # TODO: confirm/replace account
#SBATCH --time=02:00:00                      # PRE-RUN DEFAULT only; replace from the pre-run's seff
#SBATCH --mem=16G                            # PRE-RUN DEFAULT only; replace from the pre-run's seff
#SBATCH --cpus-per-task=1
#SBATCH --array=0-9
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/mi_realdata-%A_%a.out
#
# mi-posterior real-data harness (G8, .unlazy/mi-posterior/GATES.md) on
# Fir. DRAC FALLBACK ONLY: Totoro is the primary route
# (script/mi_realdata/12_totoro_run.sh, 2026-09-24). Use this file only if
# Totoro is unavailable.
#
# One array task per planned (dataset, arm, seed) cell, 10 total:
#   PanTHERIA mcar        x seeds 20260818, 20260819, 20260820  (idx 0-2)
#   PanTHERIA structured  x seeds 20260818, 20260819, 20260820  (idx 3-5)
#   AVONET    mcar        x seeds 20260818, 20260819, 20260820  (idx 6-8)
#   FishBase  structured  x seed  20260818                      (idx 9)
# (must match script/mi_realdata/lib.R::planned_cells() exactly)
#
# One CPU per task: multi_impute(draws_method = "posterior") has no n_cores
# control and runs its chains one after another in one R process, so extra
# cores would sit idle. Parallelism is across cells (the array).
#
# Prerequisites, run ONCE on the login node FROM THE REPO ROOT before
# `sbatch` (compute nodes are not relied on for git history or network):
#   module load StdEnv/2023 gcc/12.3 cuda/12.6 r/4.4.0
#   (cd script/mi_realdata && Rscript 00_fetch_masks.R --all --with-conformal)
#   git rev-parse HEAD > SHA                           # recorded in every receipt
#   R CMD build . --no-build-vignettes                 # -> pigauto_<ver>.tar.gz at the repo root
#   mkdir -p logs                                      # Slurm opens --output before the script runs
#   (00_fetch_masks.R needs origin/arc/mondrian-realdata; if this clone lacks
#   it, run that step on the Mac and copy script/mi_realdata/inputs/ over.)
#
# Per D-139 (estimate before you run): --time=02:00:00 and --mem=16G above
# are PRE-RUN DEFAULTS, sized for the smallest cell only, not the campaign.
# (A deliberately invalid --time token was tried and dropped on
# 2026-09-24: sbatch parses the #SBATCH lines before the command line and
# exits on an invalid --time, so a command-line override never applies.)
# First check the script on the Fir login node, then pre-run the smallest
# cell, from the repo root:
#   sbatch --test-only --array=6 script/mi_realdata/10_fir.sbatch
#   sbatch --array=6 script/mi_realdata/10_fir.sbatch    # avonet mcar 20260818
# then `seff <jobid>`, and edit --time/--mem here (pad about 2-3x, and size
# for FishBase, idx 9, the largest cell; its time is NOT known from the
# pre-run) before submitting --array=0-9. Do not submit --array=0-9 with
# the pre-run defaults.

set -euo pipefail
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-}"

module load StdEnv/2023 gcc/12.3 cuda/12.6 r/4.4.0

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${ROOT}"
mkdir -p logs

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
# No git on compute nodes: read conformal metrics from inputs/ only.
export MI_REALDATA_OFFLINE=1
if [[ -s "${ROOT}/SHA" ]]; then
  export MI_POST_SHA="$(cat "${ROOT}/SHA")"
else
  echo "No ${ROOT}/SHA; run 'git rev-parse HEAD > SHA' on the login node first." >&2
  exit 1
fi
# torch/lantern needs libcudart from the CUDA module (see script/submit_gnn_attribution_low_lambda_fir.sh
# for the precedent on this cluster).
export LD_LIBRARY_PATH="${EBROOTCUDA:+$EBROOTCUDA/targets/x86_64-linux/lib:}${LD_LIBRARY_PATH:-}"

# ---- fresh per-job R library on node-local $SLURM_TMPDIR, install pigauto
# from a pre-built tarball (never rely on a shared/persistent library that
# an array job would race to write). See ~/shinichi-brain/tools/drac-setup.md.
TARBALL="$(ls -t "${ROOT}"/pigauto_*.tar.gz 2>/dev/null | head -1 || true)"
if [[ -z "${TARBALL}" ]]; then
  echo "No pigauto_*.tar.gz at ${ROOT}. Run 'R CMD build .' on the login node first." >&2
  exit 1
fi
export R_LIBS_JOB="${SLURM_TMPDIR}/Rlib"
mkdir -p "${R_LIBS_JOB}"
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_JOB"), .libPaths())); install.packages(commandArgs(TRUE), repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_JOB"))' "${TARBALL}"
export R_LIBS="${R_LIBS_JOB}:${R_LIBS:-}"

# ---- array index -> (dataset, arm, seed); must match planned_cells() ----
DATASETS=(pantheria pantheria pantheria pantheria pantheria pantheria avonet avonet avonet fishbase)
ARMS=(mcar mcar mcar structured structured structured mcar mcar mcar structured)
SEEDS=(20260818 20260819 20260820 20260818 20260819 20260820 20260818 20260819 20260820 20260818)

i="${SLURM_ARRAY_TASK_ID}"
DATASET="${DATASETS[$i]}"
ARM="${ARMS[$i]}"
SEED="${SEEDS[$i]}"
echo "Cell: dataset=${DATASET} arm=${ARM} seed=${SEED} sha=${MI_POST_SHA}"

# /scratch, not /project or /scratch's parent -- purged ~60d, fine for a
# rerunnable receipt; copy anything worth keeping to /project afterwards.
OUTDIR="/scratch/${USER}/pigauto_mi_realdata/returned"
mkdir -p "${OUTDIR}"

cd "${ROOT}/script/mi_realdata"
Rscript 01_run.R "${DATASET}" "${ARM}" "${SEED}" "${OUTDIR}" \
  2>&1 | tee "${ROOT}/logs/mi_realdata-${SLURM_ARRAY_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID}.log"

echo "Rscript finished at: $(date)"
