#!/usr/bin/env bash
# script/mi_gls/diag/make_jobs.sh
#
# Writes the job lists for the posterior-MI campaign diagnosis
# (docs/dev-log/mi-posterior/diagnosis.md). One shell command per line;
# run_queue.sh executes a list with a fixed number of parallel slots.
#
# Usage: make_jobs.sh <base_dir> <diag_code_dir> <frozen_pkg_dir> <prior_pkg_dir>
#   base_dir        e.g. /home/snakagaw/pigauto_mi_posterior/diag_<sha>
#   diag_code_dir   git archive of the diagnosis commit (holds these scripts)
#   frozen_pkg_dir  frozen campaign code (69670d44f9/code): package + DGP
#   prior_pkg_dir   scratch copy of frozen_pkg_dir with ONE prior change
set -euo pipefail
BASE=$1; DIAG=$2; PKG=$3; PRIOR=$4
PKG_SHA=$(cat "$PKG/../CODE_SHA" 2>/dev/null || echo 69670d44f9)
DIAG_SHA=$(cat "$DIAG/SHA")
OUT=$BASE/out; LOG=$BASE/logs
mkdir -p "$OUT" "$LOG"
R="nice -n 10 Rscript"
rerun() { # regime rep tag niter burnin twin pkgdir
  local env="MI_POST_SHA=$PKG_SHA DIAG_SHA=$DIAG_SHA DIAG_TWIN=$6"
  [ -n "$4" ] && env="$env MI_POST_NITER=$4 MI_POST_BURNIN=$5"
  echo "cd $7 && env $env $R $DIAG/script/mi_gls/diag/rerun_cell.R $1 $2 $OUT/rerun $3 > $LOG/$3_r$1_$2.log 2>&1"
}
oracle() { # regime rep
  echo "cd $PKG && env MI_POST_SHA=$PKG_SHA DIAG_SHA=$DIAG_SHA $R $DIAG/script/mi_gls/diag/oracle_cell.R $1 $2 $OUT/oracle > $LOG/oracle_r$1_$2.log 2>&1"
}
NC="5:59 5:65 5:120 5:182 5:183 21:15 21:20 21:28 21:68 21:73 21:146 21:157 21:175 21:182 21:188 21:190 23:5 23:6 23:37 23:59 23:76 23:79 23:87 23:88 23:111 23:144 23:153 23:170 23:174 23:186 23:196 23:198"
# ---- sampler queue, longest jobs first ----
{
  for c in $NC; do rerun ${c%%:*} ${c##*:} f2_x4 20000 4000 0 "$PKG"; done
  for c in $NC; do rerun ${c%%:*} ${c##*:} f2_x2 10000 2000 0 "$PKG"; done
  for k in $(seq 1 30); do rerun 1 $k twin "" "" 1 "$PKG"; done
  for k in $(seq 1 15); do rerun 1 $k h1 "" "" 0 "$PKG"; done
  for k in $(seq 1 40); do rerun 3 $k twin "" "" 1 "$PKG"; done
  for k in $(seq 1 15); do rerun 3 $k h1 "" "" 0 "$PKG"; done
  for c in $NC; do rerun ${c%%:*} ${c##*:} f2_x1 "" "" 0 "$PKG"; done
} > "$BASE/jobs_sampler.txt"
# ---- prior-sensitivity queue (launched only if H1 looks likely) ----
{
  for k in $(seq 1 15); do rerun 1 $k prior "" "" 0 "$PRIOR"; done
  for k in $(seq 1 15); do rerun 3 $k prior "" "" 0 "$PRIOR"; done
} > "$BASE/jobs_prior.txt"
# ---- oracle queue: n = 1000 regimes first (slower) ----
{
  for k in $(seq 1 60); do oracle 3 $k; oracle 11 $k; done
  for k in $(seq 1 40); do oracle 7 $k; oracle 15 $k; oracle 4 $k; done
  for r in 1 5 9 13 2; do for k in $(seq 1 200); do oracle $r $k; done; done
} > "$BASE/jobs_oracle.txt"
wc -l "$BASE"/jobs_*.txt
