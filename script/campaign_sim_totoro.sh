#!/usr/bin/env bash
# Totoro driver for the four-arm imputation simulation (arc/imputation-sim).
#
#   ssh totoro
#   cd ~/pigauto_sim
#   bash script/campaign_sim_totoro.sh prerun 16         # 16 cells x 1 rep, all arms
#   bash script/campaign_sim_totoro.sh core   36         # 18 cells x 200 reps, BACE on seeds 1..100
#   bash script/campaign_sim_totoro.sh avonet 20
#
# Second argument = number of concurrent cell processes. Each process uses PIG_TORCH_THREADS
# (default 4) torch threads and 1 BLAS thread, so 36 processes = 144 threads (Totoro cap 150, D-143).
# Resume: a (cell, seed) whose rds exists is skipped by campaign_sim_cell.R. Detached with setsid so a
# dropped ssh cannot kill or orphan the run. Progress: tail -f logs/<stage>.log ; stop: kill -- -<pgid>
set -euo pipefail

STAGE="${1:?stage: prerun|core|factorial|avonet}"
PAR="${2:-36}"
ROOT="${PIG_SIM_ROOT:-$HOME/pigauto_sim}"
OUT="$ROOT/results/$STAGE"; [ "$STAGE" = prerun ] && OUT="$ROOT/prerun"
LOG="$ROOT/logs"; mkdir -p "$OUT" "$LOG"

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PIG_TORCH_THREADS="${PIG_TORCH_THREADS:-4}"
export NOT_CRAN=true

ARMS_ALL="gnn_on,gnn_off,gnn_off_rphylopars,freq,bace,floor"
ARMS_NOBACE="gnn_on,gnn_off,gnn_off_rphylopars,freq,floor"

# Expand the design into one line per (cell, seed): "<args for campaign_sim_cell.R>"
JOBS="$LOG/${STAGE}_jobs.txt"
Rscript "$ROOT/script/campaign_sim_design.R" --stage "$STAGE" | awk -F, -v out="$OUT" -v all="$ARMS_ALL" -v nob="$ARMS_NOBACE" '
NR > 1 {
  for (s = 1; s <= $10; s++) {
    arms = (s <= $11) ? all : nob
    printf "--dgp %s --evo %s --lambda %s --rho %s --miss %s --frac %s --n %s --seed %d --arms %s --driver --thresholds fixed --out %s\n",
           $3, $4, $5, $6, $7, $8, $9, s, arms, out
  }
}' > "$JOBS"
echo "[$(date +%FT%T)] stage=$STAGE jobs=$(wc -l < "$JOBS") parallel=$PAR out=$OUT" | tee -a "$LOG/$STAGE.log"

# One process per (cell, seed). xargs -P bounds concurrency; setsid detaches the whole group.
setsid nohup bash -c "
  cd '$ROOT' && xargs -P '$PAR' -L 1 -a '$JOBS' -I{} sh -c \
    'Rscript script/campaign_sim_cell.R {} >> \"$LOG/$STAGE.cells.log\" 2>&1' \
  ; echo \"[\$(date +%FT%T)] stage=$STAGE DONE\" >> '$LOG/$STAGE.log'
" > "$LOG/$STAGE.nohup.log" 2>&1 &
echo "pgid=$! (kill -- -$! to stop)" | tee -a "$LOG/$STAGE.log"
