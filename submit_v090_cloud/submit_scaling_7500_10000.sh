#!/bin/bash
#SBATCH --job-name=pigauto_scaling_ext
#SBATCH --account=def-<PI>
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=itchyshin@gmail.com
#SBATCH --output=%x-%j.out

# pigauto v0.9.0 — Scaling curve extension: n in {5000, 7500, 10000}
# Mixed types (2 continuous + 1 binary + 1 K=4 categorical), 25% MCAR.
# n=5000 re-run for continuity with local results; 7500 and 10000 are new.
# Estimated runtime: 6-10 hr; wall limit set to 12 hr for safety.
#
# BEFORE SUBMITTING:
#   1. Replace def-<PI> above with your actual CC account code (e.g. def-senak).
#   2. Verify pigauto is installed: R -e 'library(pigauto)' (no error = good).
#   3. Verify torch is installed: R -e 'library(torch); torch::torch_randn(3)'.
#   4. Submit from the directory containing Rscripts/ and submit_scaling_7500_10000.sh:
#        sbatch submit_scaling_7500_10000.sh

set -euo pipefail
echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "Working directory: $(pwd)"

# Load R module (check available versions with: module spider r)
module load r/4.4.0

# CC CPU-only nodes have no GPU; ensure torch uses CPU backend.
export CUDA_VISIBLE_DEVICES=""

echo "R version: $(R --version | head -1)"

# Run the benchmark. Results land in SLURM_SUBMIT_DIR.
cd "${SLURM_SUBMIT_DIR}"

Rscript Rscripts/bench_scaling_v090_extended.R \
  2>&1 | tee bench_scaling_v090_extended-${SLURM_JOB_ID}.log

echo "Rscript finished at: $(date)"

# Package results: .rds + .md + .png (if produced) + the captured log
tar czf results-scaling-${SLURM_JOB_ID}.tar.gz \
  bench_scaling_v090_extended.rds \
  bench_scaling_v090_extended.md \
  bench_scaling_v090_extended-${SLURM_JOB_ID}.log \
  $(ls bench_scaling_v090_extended.png 2>/dev/null || true) 2>/dev/null || true

echo "Results archived to: results-scaling-${SLURM_JOB_ID}.tar.gz"
echo "Job finished at: $(date)"
