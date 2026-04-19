#!/bin/bash
#SBATCH --job-name=pigauto_calibration_grid
#SBATCH --account=aip-<PI>
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-60
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=itchyshin@gmail.com
#SBATCH --output=pigauto-cal-%A_%a.out

# pigauto v0.9.1 calibration grid on Vulcan (PAICE).
# Each array task = one (scenario, mechanism, trait_type) cell x N_REPS reps.
#
# BEFORE SUBMITTING:
#   1. Replace aip-<PI> above with aip-snakagaw (or your PAICE account).
#   2. Verify pigauto and torch load: R -e 'library(pigauto); library(torch)'.
#   3. Verify user library: R -e '.libPaths()'.
#   4. Optionally edit N_REPS / CELL_GRID in Rscripts/run_calibration_cell.R.

set -euo pipefail

module load r/4.4.0
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
export CUDA_VISIBLE_DEVICES=""   # CPU-only; remove line if you want GPU

mkdir -p results logs

echo "=== Array task $SLURM_ARRAY_TASK_ID starting at $(date) ==="
echo "=== Host: $(hostname)   Account: $SLURM_JOB_ACCOUNT   Mem: ${SLURM_MEM_PER_NODE}M ==="

Rscript Rscripts/run_calibration_cell.R $SLURM_ARRAY_TASK_ID

echo "=== Array task $SLURM_ARRAY_TASK_ID finished at $(date) ==="
