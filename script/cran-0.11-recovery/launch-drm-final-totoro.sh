#!/usr/bin/env bash
# Replay the registered drmTMB recovery seeds against the final 0.11.0 adapter.
set -uo pipefail

audit_dir=/home/snakagaw/pigauto_cran_011_audit
out_dir="$audit_dir/recovery-final"
mkdir -p "$out_dir"
export R_LIBS="$audit_dir/library:/home/snakagaw/R/lib"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

case "${1:-}" in
  smoke) seeds=(2026100601) ;;
  remaining) seeds=(2026100602 2026100603 2026100604 2026100605 2026100606
                    2026100607 2026100608 2026100609 2026100610) ;;
  *) echo 'Usage: launch-drm-final-totoro.sh smoke|remaining' >&2; exit 2 ;;
esac

printf 'SUPERVISOR_START drm %s %s\n' "$1" "$(date -Is)"
for seed in "${seeds[@]}"; do
  printf 'SEED_START drm %s %s\n' "$seed" "$(date -Is)"
  timeout --signal=TERM --kill-after=30s 900 Rscript --vanilla \
    "$audit_dir/run.R" drm "$seed" "$out_dir" \
    > "$audit_dir/drm-final-$seed.log" 2>&1
  rc=$?
  printf 'SEED_EXIT drm %s %s %s\n' "$seed" "$rc" "$(date -Is)"
  if [ "$rc" -ne 0 ]; then exit "$rc"; fi
done
printf 'SUPERVISOR_EXIT drm %s 0 %s\n' "$1" "$(date -Is)"
