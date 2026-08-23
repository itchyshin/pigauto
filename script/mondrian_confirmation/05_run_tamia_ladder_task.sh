#!/bin/bash
# One exclusive GPU task from 04_run_tamia_gpu_smoke.sbatch.
set -euo pipefail
PIGAUTO_STAGEC_ROOT="$1"
n_tips="$2"
attempt_id="${3:?expected Tamia Slurm job id}"
export R_LIBS_USER="$PIGAUTO_STAGEC_ROOT/Rlib"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="${EBROOTCUDA:?cuda module not loaded}/targets/x86_64-linux/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Attempts are immutable: a repaired ladder must never reuse a prior
# split.rds and thereby falsely appear to have executed.
result_dir="$PIGAUTO_STAGEC_ROOT/results/attempt-${attempt_id}/fishbase-${n_tips}-split-gpu"
mkdir -p "$result_dir"
nvidia-smi --query-gpu=index,uuid,pci.bus_id,name,driver_version,memory.total --format=csv,noheader > "$result_dir/gpu.txt"
{
  printf 'CUDA_VISIBLE_DEVICES=%s\n' "${CUDA_VISIBLE_DEVICES:-UNSET}"
  printf 'SLURM_PROCID=%s\n' "${SLURM_PROCID:-UNSET}"
  printf 'SLURM_LOCALID=%s\n' "${SLURM_LOCALID:-UNSET}"
  Rscript --vanilla -e 'stopifnot(requireNamespace("torch", quietly = TRUE)); cat("cuda_available=", torch::cuda_is_available(), "\\n", sep = ""); if (!torch::cuda_is_available()) stop("CUDA is unavailable inside this allocation")'
} > "$result_dir/device.txt"
Rscript --vanilla "$PIGAUTO_STAGEC_ROOT/source/pigauto/script/mondrian_confirmation/01_run_masked_confirmation.R" \
  "$PIGAUTO_STAGEC_ROOT/input/fishbase-${n_tips}-input.rds" \
  "$result_dir" \
  20260820 10 split > "$result_dir/runner.log" 2>&1
