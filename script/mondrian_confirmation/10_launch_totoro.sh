#!/bin/bash
# Launch one (dataset, arm, seed) cell on Totoro: both methods, sequential,
# one compute thread, detached. Usage:
#   10_launch_totoro.sh <root> <dataset> <arm> <seed> [epochs] [methods]
# <root> holds lib/, input/<dataset>-full-input.rds, source/ (this repo) and results/.
set -euo pipefail
ROOT="$1"; DS="$2"; ARM="$3"; SEED="$4"; EP="${5:-500}"; METHODS="${6:-mondrian,split}"  # mondrian first: split reuses its thresholds for stratum labels
OUT="$ROOT/results/$DS-$ARM-m$SEED"
mkdir -p "$OUT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export R_LIBS="$ROOT/lib:$HOME/R/lib"
nohup Rscript --vanilla "$ROOT/source/script/mondrian_confirmation/01_run_masked_confirmation.R" \
  "$ROOT/input/$DS-full-input.rds" "$OUT" "$ARM" "$SEED" "$EP" "$METHODS" \
  > "$OUT/run.log" 2>&1 &
echo $! > "$OUT/pid"
echo "launched $DS $ARM $SEED pid $(cat "$OUT/pid")"
