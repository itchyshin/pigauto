#!/usr/bin/env bash
set -euo pipefail
: "${VALIDATION_APPROVAL:?Exact measured-pilot approval receipt required}"
: "${VALIDATION_WORKERS:?Approved worker count required}"
: "${VALIDATION_CONTEXT:?Choose totoro or drac}"
: "${VALIDATION_JOB_TIMEOUT_SEC:?Approved per-job timeout required}"
: "${VALIDATION_BUDGET_SEC:?Approved per-worker wall budget required}"
[[ "$VALIDATION_WORKERS" =~ ^[0-9]+$ ]] && (( VALIDATION_WORKERS >= 1 && VALIDATION_WORKERS <= 48 )) || exit 2
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PIG_TORCH_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
command -v timeout >/dev/null || exit 2
export VALIDATION_APPROVAL VALIDATION_WORKERS VALIDATION_BUDGET_SEC VALIDATION_JOB_TIMEOUT_SEC VALIDATION_CONTEXT
Rscript script/validation_20261018/run.R --mode approval-check --approval "$VALIDATION_APPROVAL" --workers "$VALIDATION_WORKERS" --budget-sec "$VALIDATION_BUDGET_SEC" --job-timeout "$VALIDATION_JOB_TIMEOUT_SEC"
run_all() {
  pids=()
  cleanup() { for p in "${pids[@]}"; do kill "$p" 2>/dev/null || true; done; }
  trap cleanup EXIT
  for ((i=1;i<=VALIDATION_WORKERS;i++)); do
    Rscript script/validation_20261018/run.R --mode production --approval "$VALIDATION_APPROVAL" --workers "$VALIDATION_WORKERS" --worker-index "$i" --budget-sec "$VALIDATION_BUDGET_SEC" --job-timeout "$VALIDATION_JOB_TIMEOUT_SEC" &
    pids+=("$!")
  done
  failed=0
  for p in "${pids[@]}"; do wait "$p" || failed=1; done
  return "$failed"
}
export -f run_all
/usr/bin/time -v timeout --signal=TERM --kill-after=20s "$VALIDATION_BUDGET_SEC" bash -c run_all
Rscript script/validation_20261018/summarise.R
