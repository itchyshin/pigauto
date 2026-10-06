#!/usr/bin/env bash
set -euo pipefail

audit_dir=/home/snakagaw/pigauto_cran_011_audit
out_dir="$audit_dir/recovery-parallel"
mkdir -p "$out_dir"
export R_LIBS="$audit_dir/library:/home/snakagaw/R/lib"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

run_seed() {
  local seed="$1"
  timeout --signal=TERM --kill-after=30s 7200 Rscript --vanilla \
    "$audit_dir/run.R" gllvm "$seed" "$out_dir" \
    >"$out_dir/gllvm-$seed.log" 2>&1
}

start=$(date +%s)
echo "PARALLEL_FEASIBILITY_START $(date --iso-8601=seconds)"
run_seed 2026100602 &
pid2=$!
run_seed 2026100603 &
pid3=$!
cleanup() {
  kill "$pid2" "$pid3" 2>/dev/null || true
  wait "$pid2" "$pid3" 2>/dev/null || true
}
trap cleanup EXIT INT TERM
set +e
wait "$pid2"
status2=$?
wait "$pid3"
status3=$?
set -e
elapsed=$(($(date +%s) - start))
echo "PARALLEL_FEASIBILITY_EXIT seed2=$status2 seed3=$status3 elapsed_s=$elapsed"
test "$status2" -eq 0
test "$status3" -eq 0
