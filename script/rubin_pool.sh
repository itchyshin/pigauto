#!/usr/bin/env bash
# Pool the rubin-freq-bace campaign results from every host into per-host folders on the Mac (not in Dropbox:
# thousands of rds files). Results of record: BACE from results/bace on nibi, fir and rorqual; frequentist from
# fir ~/pigauto_rubin_f3/results/freq (the fixed rerun). Also copies the as-shipped BACE failure records.
#   bash script/rubin_pool.sh [POOL_DIR]      (default ~/pigauto_rubin_pool)
# Never aggregate while this runs (v1 lesson).
set -euo pipefail
POOL="${1:-$HOME/pigauto_rubin_pool}"
mkdir -p "$POOL"/{bace,bace_asshipped_failed,freq}/{nibi,fir,rorqual}
pull() { rsync -a --include='*.rds' --exclude='*' "$1:$2/" "$POOL/$3/" 2>/dev/null || echo "note: nothing at $1:$2"; }
pull nibi    projects/def-snakagaw/snakagaw/pigauto_rubin/results/bace                   bace/nibi
pull fir     pigauto_rubin/results/bace                                                  bace/fir
pull rorqual projects/def-snakagaw/snakagaw/pigauto_rubin/results/bace                   bace/rorqual
pull nibi    projects/def-snakagaw/snakagaw/pigauto_rubin/results/bace_asshipped_failed  bace_asshipped_failed/nibi
pull fir     pigauto_rubin/results/bace_asshipped_failed                                 bace_asshipped_failed/fir
pull rorqual projects/def-snakagaw/snakagaw/pigauto_rubin/results/bace_asshipped_failed  bace_asshipped_failed/rorqual
pull fir     pigauto_rubin_f3/results/freq                                               freq/fir
for d in bace bace_asshipped_failed freq; do echo "$d: $(find "$POOL/$d" -name '*.rds' | wc -l) files"; done
