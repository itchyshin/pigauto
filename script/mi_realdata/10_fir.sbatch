#!/bin/bash
#SBATCH --job-name=pigauto_mi_realdata
#SBATCH --account=def-snakagaw               # TODO: confirm/replace account
#SBATCH --time=00:00:00                      # TODO: set from a pre-run test (design.md/G8); no estimate exists yet
#SBATCH --mem=0G                             # TODO: set from a pre-run test
#SBATCH --cpus-per-task=4
#SBATCH --array=0-9
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/mi_realdata-%A_%a.out
#
# mi-posterior real-data harness (G8, .unlazy/mi-posterior/GATES.md) on
# Fir. One array task per planned (dataset, arm, seed) cell -- 10 total:
#   PanTHERIA mcar        x seeds 20260818, 20260819, 20260820  (idx 0-2)
#   PanTHERIA structured  x seeds 20260818, 20260819, 20260820  (idx 3-5)
#   AVONET    mcar        x seeds 20260818, 20260819, 20260820  (idx 6-8)
#   FishBase  structured  x seed  20260818                      (idx 9)
# (must match script/mi_realdata/lib.R::planned_cells() exactly)
#
# Prerequisites (run ONCE on the login node before `sbatch`, from
# $HOME/pigauto or wherever this worktree lives on Fir -- compute nodes
# should not be relied on for outbound git/network access):
#   module load StdEnv/2023 gcc/12.3 cuda/12.6 r/4.4.0
#   cd script/mi_realdata && Rscript 00_fetch_masks.R --all
#   cd ../.. && R CMD build . --no-build-vignettes   # -> pigauto_<ver>.tar.gz
# The tarball path below assumes that file sits at the repo root.
#
# Per D-139 (estimate before you run): --time/--mem above are PLACEHOLDERS.
# Run ONE array task by hand first (or `--array=0`), check `seff <jobid>`,
# and fill in real numbers before submitting the full --array=0-9.

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
# torch/lantern needs libcudart from the CUDA module (see script/submit_gnn_attribution_low_lambda_fir.sh
# for the precedent on this cluster).
export LD_LIBRARY_PATH="${EBROOTCUDA:+$EBROOTCUDA/targets/x86_64-linux/lib:}${LD_LIBRARY_PATH:-}"

# ---- fresh per-job R library on node-local $SLURM_TMPDIR, install pigauto
# from a pre-built tarball (never rely on a shared/persistent library that
# an array job would race to write). See ~/shinichi-brain/tools/drac-setup.md.
TARBALL="$(ls -t "${ROOT}"/pigauto_*.tar.gz | head -1)"
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
echo "Cell: dataset=${DATASET} arm=${ARM} seed=${SEED}"

# /scratch, not /project or /scratch's parent -- purged ~60d, fine for a
# rerunnable receipt; copy anything worth keeping to /project afterwards.
OUTDIR="/scratch/${USER}/pigauto_mi_realdata/returned"
mkdir -p "${OUTDIR}"

cd "${ROOT}/script/mi_realdata"
Rscript 01_run.R "${DATASET}" "${ARM}" "${SEED}" "${OUTDIR}" \
  2>&1 | tee "${ROOT}/logs/mi_realdata-${SLURM_ARRAY_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID}.log"

echo "Rscript finished at: $(date)"
