#!/usr/bin/env bash
# Totoro driver for the pig_post CAMPAIGN of the Rubin study (pigauto posterior MI, PR #189).
# Approved by Shinichi 2026-10-01 ("go Totoro, option b"): 3,600 fits, <= 110 cores, option (b) = n = 100 fits may
# extend up to 6 times (posterior_control$max_extend = 6) instead of pigauto's default 3. Pre-run and budget:
# docs/dev-log/arc/2026-10-01-rubin-pigauto-prerun.md (about 980 core-hours measured, 1,225 with margin, ~11-12 h).
#
#   ssh totoro; cd ~/pigauto_rubin
#   PIGAUTO_SHA=<merged main sha> CONFIRM=yes bash script/rubin_pig_campaign_totoro.sh 110
#
# Grid: n {100, 300, 1000} x lambda {0.3, 0.7, 1} x rho {0, 0.5} x seeds 1..200, arm pig_post, M = 20, MCAR 30%: the
# same datasets as the stored freq and BACE campaign. Longest cells first. Each cell is capped at 3 h. Output goes to its
# own directory; resume skips existing rds, so a stopped run restarts with the same command.
set -euo pipefail

PAR="${1:-110}"
ROOT="${RUBIN_ROOT:-$HOME/pigauto_rubin}"
OUT="$ROOT/results_pig"; LOG="$ROOT/logs"; mkdir -p "$OUT" "$LOG"
[ "$PAR" -le 110 ] || { echo "PAR=$PAR exceeds the 110 cores approved for this campaign"; exit 1; }

free -g
BUSY=$(ps -u "$USER" -o pcpu= | awk '{s+=$1} END {printf "%d", s/100 + 0.5}')
echo "cores already busy for $USER: $BUSY"
[ $(( BUSY + PAR )) -le 150 ] || { echo "PAR=$PAR + $BUSY busy exceeds the 150-core cap (D-143); lower PAR to $(( 150 - BUSY ))"; exit 1; }

export R_LIBS="$ROOT/rlib:${R_LIBS:-}"
: "${PIGAUTO_SHA:?set PIGAUTO_SHA to the merged main sha}"
Rscript -e "d <- utils::packageDescription('pigauto'); s <- if (is.null(d\$RemoteSha)) '' else d\$RemoteSha;
  stopifnot('installed pigauto is not the expected sha' = identical(s, '$PIGAUTO_SHA'),
            'installed pigauto has no posterior draws' = 'posterior' %in% eval(formals(pigauto::multi_impute)\$draws_method));
  cat('pigauto', s, 'OK\n')"
for f in campaign_gnn_off_lib.R rubin_lib.R rubin_freq.R rubin_bace.R rubin_cell.R rubin_pigauto.R; do
  [ -s "$ROOT/script/$f" ] || { echo "missing $ROOT/script/$f (rsync the lane's script/ first)"; exit 1; }
done
grep -q 'pig_max_extend' "$ROOT/script/rubin_cell.R" || { echo "deployed rubin_cell.R lacks --pig_max_extend"; exit 1; }

[ "${CONFIRM:-no}" = yes ] || { echo "dry run: set CONFIRM=yes to launch"; exit 0; }

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

JOBS="$LOG/campaign_pig_jobs.txt"; : > "$JOBS"
for n in 1000 300 100; do
  extra=""; [ "$n" = 100 ] && extra="--pig_max_extend 6"
  for lambda in 0.3 0.7 1; do for rho in 0 0.5; do for seed in $(seq 1 200); do
    echo "--n $n --seed $seed --lambda $lambda --rho $rho --M 20 --arms pig_post $extra --out $OUT" >> "$JOBS"
  done; done; done
done
echo "[$(date +%FT%T)] campaign_pig jobs=$(wc -l < "$JOBS") parallel=$PAR out=$OUT sha=$PIGAUTO_SHA" | tee -a "$LOG/campaign_pig.log"

setsid nohup bash -c "
  xargs -P $PAR -L 1 -a '$JOBS' bash -c 'timeout 10800 Rscript \"$ROOT/script/rubin_cell.R\" \"\$@\" >> \"$LOG/campaign_pig_cells.log\" 2>&1' _
  echo \"[\$(date +%FT%T)] campaign done\" >> '$LOG/campaign_pig.log'
" > /dev/null 2>&1 &
echo "launched; pgid $(ps -o pgid= $! | tr -d ' '); stop with: kill -- -<pgid>"
