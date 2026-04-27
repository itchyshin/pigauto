#!/bin/bash
# script/run_all_per_type_benches.sh
#
# Sequentially re-run all 9 bundled per-type benches with the row-fix
# (commit a3e6d39) applied.  Caches the pre-fix RDS files in
# script/_pre_row_fix_cache/ via prep step (already done before this
# script runs).
#
# Each bench is sequential to itself (uses internal mc.cores=16 PSOCK
# parallelism), so we run them one at a time to avoid oversubscription.
#
# Usage:
#   bash script/run_all_per_type_benches.sh

set -e
cd "/Users/z3437171/Dropbox/Github Local/pigauto"

LOG_DIR="script"
RSCRIPT=/usr/local/bin/Rscript

# Throttle PSOCK parallelism: 16 workers swamped this 18-core machine
# (each ~1.2 GB resident -> ~19 GB total -> swapped).  8 workers at
# ~10 GB total leaves headroom for OS and the parent R process.
export MC_CORES=8

# bench_continuous already started by hand; if it's still running we
# wait for it. Otherwise start it.
benches=(
  bench_continuous
  bench_binary
  bench_count
  bench_ordinal
  bench_categorical
  bench_proportion
  bench_zi_count
  bench_multi_proportion
  bench_missingness_mechanism
)

for b in "${benches[@]}"; do
  log_file="${LOG_DIR}/${b}_postfix.log"
  if [ -f "${LOG_DIR}/${b}.rds" ] && grep -q "All cells already done\|=== DONE\|done$" "${log_file}" 2>/dev/null; then
    echo "[$(date +%H:%M:%S)] ${b} already done -- skipping"
    continue
  fi
  echo "[$(date +%H:%M:%S)] === Starting ${b} ==="
  "$RSCRIPT" "${LOG_DIR}/${b}.R" > "${log_file}" 2>&1
  rc=$?
  echo "[$(date +%H:%M:%S)] === Finished ${b} (rc=${rc}) ==="
  if [ "${rc}" != "0" ]; then
    echo "[$(date +%H:%M:%S)] ${b} returned non-zero exit; continuing anyway"
  fi
done

echo "[$(date +%H:%M:%S)] === All 9 per-type benches done ==="
