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
pids=()
cleanup() {
  if ((${#pids[@]})); then
    kill "${pids[@]}" 2>/dev/null || true
    wait "${pids[@]}" 2>/dev/null || true
  fi
}
trap cleanup EXIT INT TERM
echo "GLLVM_REMAINING_START $(date --iso-8601=seconds)"

waves=("2026100604 2026100605" "2026100606 2026100607"
       "2026100608 2026100609" "2026100610")
for wave_index in "${!waves[@]}"; do
  read -r -a seeds <<<"${waves[$wave_index]}"
  wave_start=$(date +%s)
  pids=()
  for seed in "${seeds[@]}"; do
    run_seed "$seed" &
    pids+=("$!")
  done
  statuses=()
  for pid in "${pids[@]}"; do
    set +e
    wait "$pid"
    statuses+=("$?")
    set -e
  done
  pids=()
  wave_elapsed=$(($(date +%s) - wave_start))
  elapsed=$(($(date +%s) - start))
  echo "GLLVM_WAVE_END index=$((wave_index + 1)) seeds=${waves[$wave_index]} statuses=${statuses[*]} wave_s=$wave_elapsed total_s=$elapsed"
  for status in "${statuses[@]}"; do
    if ((status != 0)); then
      echo "GLLVM_STOP failed_seed_wave=$((wave_index + 1))"
      exit 1
    fi
  done
  remaining_waves=$((${#waves[@]} - wave_index - 1))
  projected_total=$((elapsed + remaining_waves * wave_elapsed))
  if ((projected_total > 10800)); then
    echo "GLLVM_STOP projected_total_s=$projected_total exceeds_3h=10800"
    exit 2
  fi
done
echo "GLLVM_REMAINING_OK elapsed_s=$(($(date +%s) - start))"
