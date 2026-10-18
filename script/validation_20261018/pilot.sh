#!/usr/bin/env bash
set -euo pipefail
# Serial pilot: one CPU, 30 minutes overall, ten minutes per job.
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 PIG_TORCH_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
OUT="${1:-script/validation_20261018/results/pilot}"
mkdir -p "$OUT"
command -v timeout >/dev/null 2>&1 || { echo "GNU timeout is required; no unbounded fallback" >&2; exit 2; }
[[ "$(uname -s)" == Linux ]] || { echo "Use this bounded resource-measured launcher on Linux" >&2; exit 2; }
/usr/bin/time -v timeout --signal=TERM --kill-after=20s 1800s bash -c '
  set -euo pipefail
  Rscript script/validation_20261018/run.R --mode pilot --out "$1" --budget-sec 1700 --job-timeout 600
  Rscript script/validation_20261018/summarise.R --pilot --out "$1"
' validation-pilot "$OUT"
