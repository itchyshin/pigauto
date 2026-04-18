#!/bin/bash
#SBATCH --job-name=pigauto_avonet_missingness
#SBATCH --account=def-<PI>
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=itchyshin@gmail.com
#SBATCH --output=%x-%j.out

# pigauto v0.9.0 — AVONET full-scale missingness sweep
# n = 9,993 species; 3 missingness levels x 3 methods x 500 epochs each
# Estimated runtime: 10-20 hr; wall limit set to 24 hr for safety.
#
# BEFORE SUBMITTING:
#   1. Replace def-<PI> above with your actual CC account code (e.g. def-senak).
#   2. Verify pigauto is installed: R -e 'library(pigauto)' (no error = good).
#   3. Verify torch is installed: R -e 'library(torch); torch::torch_randn(3)'.
#   4. Submit from the directory containing Rscripts/ and submit_avonet.sh:
#        sbatch submit_avonet.sh

set -euo pipefail
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "Working directory: $(pwd)"

# Load R module (check available versions with: module spider r)
module load r/4.4.0

# torch uses ~/.local/share/torch by default; tell it explicitly if needed.
# Uncomment and adjust if torch complains about missing libtorch:
# export TORCH_HOME="${HOME}/.local/share/torch"

# CC CPU-only nodes have no GPU; ensure torch uses CPU backend.
export CUDA_VISIBLE_DEVICES=""

echo "R version: $(R --version | head -1)"

# Run the benchmark. Rscript is called from the job submission directory
# (SLURM_SUBMIT_DIR), which is where results will land.
cd "${SLURM_SUBMIT_DIR}"

Rscript Rscripts/bench_avonet_missingness.R \
  2>&1 | tee bench_avonet_missingness-${SLURM_JOB_ID}.log

echo "Rscript finished at: $(date)"

# Package results: .rds + .md + the captured log
tar czf results-avonet-${SLURM_JOB_ID}.tar.gz \
  bench_avonet_missingness.rds \
  bench_avonet_missingness.md \
  bench_avonet_missingness-${SLURM_JOB_ID}.log 2>/dev/null || true

echo "Results archived to: results-avonet-${SLURM_JOB_ID}.tar.gz"
echo "Job finished at: $(date)"
