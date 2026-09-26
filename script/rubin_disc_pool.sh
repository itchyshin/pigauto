#!/usr/bin/env bash
# Pool the discrete re-run of the Rubin study (results_disc/, written with --save_imp --discrete) from every host into
# per-host folders on the Mac, separate from the continuous pool (~/pigauto_rubin_pool, never touched here).
# Launch record: BACE n = 300 seeds 1-100 on fir and 101-200 on rorqual; BACE n = 100 started on nibi, and its 880 queued
# seeds were moved to rorqual (400) and fir (480) at 10:12; freq + castor on Totoro.
#   bash script/rubin_disc_pool.sh [POOL_DIR]      (default ~/pigauto_rubin_disc_pool)
# Never aggregate while this runs.
set -euo pipefail
POOL="${1:-$HOME/pigauto_rubin_disc_pool}"
mkdir -p "$POOL"/{bace,freq}/{nibi,fir,rorqual,totoro}
pull() { rsync -a --include='*.rds' --exclude='*' "$1:$2/" "$POOL/$3/" 2>/dev/null || echo "note: nothing at $1:$2"; }
for set in bace freq; do
  pull nibi    projects/def-snakagaw/snakagaw/pigauto_rubin/results_disc/$set  $set/nibi
  pull fir     pigauto_rubin/results_disc/$set                                 $set/fir
  pull rorqual projects/def-snakagaw/snakagaw/pigauto_rubin/results_disc/$set  $set/rorqual
  pull totoro  pigauto_rubin/results_disc/$set                                 $set/totoro
done
for d in bace freq; do echo "$d: $(find "$POOL/$d" -name '*.rds' | wc -l | tr -d ' ') files"; done
