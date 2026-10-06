#!/usr/bin/env bash
# One registered standard-control gllvm seed. Estimated ceiling: 120 minutes.
set -uo pipefail
audit_dir=/home/snakagaw/pigauto_cran_011_audit
export R_LIBS="$audit_dir/library:/home/snakagaw/R/lib"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
printf 'SUPERVISOR_START %s\n' "$(date -Is)"
timeout --signal=TERM --kill-after=30s 7200 Rscript --vanilla \
  "$audit_dir/run.R" gllvm 2026100601 "$audit_dir/recovery"
status=$?
printf 'SUPERVISOR_EXIT %s %s\n' "$status" "$(date -Is)"
exit "$status"
