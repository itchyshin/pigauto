#!/usr/bin/env bash
# Totoro driver for part of the rubin-freq-bace campaign (no queue): runs one rubin_cell.R per line of a task file,
# PAR at a time, detached with setsid. Refuses if this user would exceed D-143's 150 cores or free memory is short.
#   ssh totoro; cd ~/pigauto_rubin
#   bash script/rubin_totoro.sh <tasks_file> <PAR> <GB_per_task>
# Task file: one line of rubin_cell.R arguments per fit. Resume is free (rubin_cell.R skips an existing rds).
# Stop: kill -- -<pgid> (printed at launch).
set -euo pipefail
TASKS="${1:?tasks file}"; PAR="${2:?parallel}"; GB="${3:?GB per task}"
ROOT="${RUBIN_ROOT:-$HOME/pigauto_rubin}"; LOG="$ROOT/logs"; mkdir -p "$LOG"
BUSY=$(ps -u "$USER" -o pcpu= | awk '{s+=$1} END {printf "%d", s/100 + 0.5}')
AVAIL=$(free -g | awk '/^Mem:/ {print $7}')
echo "busy cores for $USER: $BUSY; available memory: ${AVAIL} GB; tasks: $(wc -l < "$TASKS")"
[ $(( BUSY + PAR )) -le 150 ] || { echo "PAR=$PAR + $BUSY busy exceeds the 150-core cap (D-143)"; exit 1; }
[ "$AVAIL" -ge $(( PAR * GB + 100 )) ] || { echo "need ~$(( PAR * GB + 100 )) GB free, have $AVAIL"; exit 1; }
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
TAG=$(basename "$TASKS" .txt)
setsid nohup bash -c "
  cd '$ROOT' && xargs -P '$PAR' -L 1 -a '$TASKS' -I{} sh -c 'Rscript script/rubin_cell.R {} >> \"$LOG/$TAG.cells.log\" 2>&1' \
  ; echo \"[\$(date +%FT%T)] $TAG DONE\" >> '$LOG/$TAG.log'
" > "$LOG/$TAG.nohup.log" 2>&1 &
echo "[$(date +%FT%T)] $TAG started, pgid=$! (kill -- -$! to stop)" | tee -a "$LOG/$TAG.log"
