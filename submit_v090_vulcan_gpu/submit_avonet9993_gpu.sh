#!/bin/bash
#SBATCH --job-name=pigauto_avonet9993_gpu
#SBATCH --account=aip-snakagaw
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=itchyshin@gmail.com
#SBATCH --output=pigauto-avonet9993-%j.out

# pigauto × BACE head-to-head on AVONET 9,993 species (GPU).
#
# Pigauto's GNN trains on GPU (the O(n²) attention at n=10k fills the
# device and benefits from parallel kernels).  BACE stays on CPU
# (MCMCglmm has no GPU support).  One SLURM job, one GPU allocated,
# BACE running on the 8 CPU cores in parallel while pigauto runs on GPU.
#
# BEFORE SUBMITTING:
#   1. Verify pigauto + BACE + torch all load:
#        R -e 'library(pigauto); library(BACE); library(torch); torch::torch_randn(3L)'
#      (If torch complains about libtorch: torch::install_torch() on the
#      login node first.  It has internet; compute nodes don't.)
#   2. Verify GPU visible from R:
#        R -e 'library(torch); cat("cuda_is_available:", torch::cuda_is_available(), "\n")'
#      Must print TRUE.
#
# Submit:
#   cd ~/pigauto_vulcan_gpu
#   sbatch submit_avonet9993_gpu.sh
#
# Expected wall:
#   - preprocess + graph       : ~5 min  (Rphylopars + cophenetic)
#   - fit_baseline              : ~5–15 min (Rphylopars at n=10k)
#   - pigauto GNN train on GPU  : ~30–60 min (500 epochs, attention at n=10k)
#   - pigauto predict           : ~5 min
#   - BACE MCMCglmm on 7 traits : ~4–8 hr (dominates wall-clock)
#   - Total                     : ~6–10 hr
# 12-hr wall is comfortably above that ceiling.

set -euo pipefail

echo "=== Job started at $(date) ==="
echo "=== Host: $(hostname) ==="
echo "=== SLURM_JOB_ID: ${SLURM_JOB_ID} ==="
echo "=== GPU(s) allocated: ${CUDA_VISIBLE_DEVICES:-none} ==="

module load r/4.4.0 cuda/12 cudnn/9 2>/dev/null || module load r/4.4.0
export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.4"

# Check torch sees the GPU (fatal if not)
R -e 'suppressPackageStartupMessages(library(torch));
      if (!torch::cuda_is_available()) {
        stop("CUDA not available inside the job. Check that --gres=gpu:1 landed and libtorch sees the card.")
      }
      cat("cuda_is_available:", torch::cuda_is_available(), "\n")
      cat("cuda_device_count:", torch::cuda_device_count(), "\n")'

echo "=== Running driver at $(date) ==="
Rscript Rscripts/run_avonet9993_bace.R \
  2>&1 | tee run_avonet9993_bace_${SLURM_JOB_ID}.log

echo "=== Job finished at $(date) ==="

# Package results for rsync_results_back.sh
tar czf results-avonet9993-${SLURM_JOB_ID}.tar.gz \
  bench_avonet9993_bace.rds \
  bench_avonet9993_bace.md \
  run_avonet9993_bace_${SLURM_JOB_ID}.log 2>/dev/null || true

echo "=== Results archived to results-avonet9993-${SLURM_JOB_ID}.tar.gz ==="
