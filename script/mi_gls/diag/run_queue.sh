#!/usr/bin/env bash
# script/mi_gls/diag/run_queue.sh <joblist> <slots>
# Runs each line of <joblist> as a shell command, <slots> at a time.
# BLAS/OpenMP threads are pinned to 1 HERE, in the launching shell, so every
# R child inherits them before BLAS initialises.
set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export R_LIBS_USER=/home/snakagaw/R/lib
xargs -P "$2" -d '\n' -I{} bash -c '{}' < "$1"
echo "QUEUE_DONE $1 $(date)"
