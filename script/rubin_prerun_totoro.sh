#!/usr/bin/env bash
# Totoro driver for the BACE settings pre-run of the rubin-freq-bace lane
# (docs/dev-log/arc/2026-09-24-rubin-prerun-plan.md). NOT launched without Shinichi's approval (D-139).
#
#   ssh totoro            # through the ~/.ssh/cm-* ControlMaster socket
#   cd ~/pigauto_rubin
#   CONFIRM=yes bash script/rubin_prerun_totoro.sh 100      # 100 concurrent single-threaded cells
#
# Grid: BACE runs {5, 10, 15} x nitt {50k, 100k} (burnin 20%, thin keeping 1,600 samples), n_final = M = 20,
# cells lambda {0.3, 0.7} x rho {0, 0.5} x n {100, 300}, MCAR 30%, seeds 1..5 -> 6 x 8 x 5 = 240 BACE fits.
# Each setting writes to its own directory because the cell tag does not carry runs/nitt.
# Resume: rubin_cell.R skips a cell whose rds exists. Detached with setsid; stop with kill -- -<pgid>.
set -euo pipefail

PAR="${1:-100}"
ROOT="${RUBIN_ROOT:-$HOME/pigauto_rubin}"
OUT="$ROOT/prerun"; LOG="$ROOT/logs"; mkdir -p "$OUT" "$LOG"
[ "$PAR" -le 150 ] || { echo "PAR=$PAR exceeds the 150-core Totoro cap (D-143)"; exit 1; }

# Shared-machine preflight: print free memory and the heaviest users before anything starts.
free -g
ps -eo user,rss --no-headers | awk '{r[$1]+=$2} END {for (u in r) printf "%-12s %8.1f GB\n", u, r[u]/1048576}' \
  | sort -k2 -nr | head -8
# D-143 caps THIS USER at 150 Totoro cores across all lanes (the lambda-default lane runs as the same user).
BUSY=$(ps -u "$USER" -o pcpu= | awk '{s+=$1} END {printf "%d", s/100 + 0.5}')
echo "cores already busy for $USER: $BUSY"
[ $(( BUSY + PAR )) -le 150 ] || { echo "PAR=$PAR + $BUSY busy exceeds the 150-core cap; lower PAR to $(( 150 - BUSY ))"; exit 1; }
AVAIL_GB=$(free -g | awk '/^Mem:/ {print $7}')
NEED_GB=$(( PAR * 2 ))            # BACE at n <= 300 stays under ~2 GB RSS per process (v1)
[ "$AVAIL_GB" -ge "$NEED_GB" ] || { echo "only ${AVAIL_GB} GB available, need ~${NEED_GB} GB"; exit 1; }

# Same BACE build as the Mac smoke: the installed bace_final_imp() must draw posterior predictive values
# (.predict_bace(..., sample = TRUE)). The stale 2026-04-01 source takes posterior means instead.
Rscript -e 'stopifnot("BACE build lacks sample = TRUE in bace_final_imp" =
  any(grepl("sample = TRUE", deparse(BACE:::bace_final_imp), fixed = TRUE))); cat("BACE build OK\n")'
# Code deployed and loadable before any process starts.
for f in campaign_gnn_off_lib.R rubin_lib.R rubin_freq.R rubin_bace.R rubin_cell.R; do
  [ -s "$ROOT/script/$f" ] || { echo "missing $ROOT/script/$f (rsync the lane's script/ first)"; exit 1; }
done

[ "${CONFIRM:-no}" = yes ] || { echo "dry run: set CONFIRM=yes to launch (D-139: needs Shinichi's approval)"; exit 0; }

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

JOBS="$LOG/prerun_jobs.txt"; : > "$JOBS"
for runs in 5 10 15; do
  for nitt in 50000 100000; do
    burnin=$(( nitt / 5 )); thin=$(( (nitt - burnin) / 1600 ))
    for lambda in 0.3 0.7; do for rho in 0 0.5; do for n in 100 300; do for seed in 1 2 3 4 5; do
      echo "--n $n --seed $seed --lambda $lambda --rho $rho --M 20 --arms bace,bace_chain,bace_resid" \
           "--bace_nitt $nitt --bace_burnin $burnin --bace_thin $thin --bace_runs $runs" \
           "--out $OUT/runs${runs}_nitt${nitt}" >> "$JOBS"
    done; done; done; done
  done
done
# Longest fits first so the tail of the run is short.
sort -t' ' -k2,2nr "$JOBS" -o "$JOBS"
echo "[$(date +%FT%T)] prerun jobs=$(wc -l < "$JOBS") parallel=$PAR out=$OUT" | tee -a "$LOG/prerun.log"

setsid nohup bash -c "
  cd '$ROOT' && xargs -P '$PAR' -L 1 -a '$JOBS' -I{} sh -c \
    'Rscript script/rubin_cell.R {} >> \"$LOG/prerun.cells.log\" 2>&1' \
  ; echo \"[\$(date +%FT%T)] prerun DONE\" >> '$LOG/prerun.log'
" > "$LOG/prerun.nohup.log" 2>&1 &
echo "pgid=$! (kill -- -$! to stop)" | tee -a "$LOG/prerun.log"
