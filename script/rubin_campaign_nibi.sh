#!/usr/bin/env bash
# nibi job-array driver for the rubin-freq-bace CAMPAIGN (core cells: lambda {0.3, 0.7, 1} x rho {0, 0.5} per n).
# NOT to be launched without Shinichi's approval of the settings and budget
# (docs/dev-log/arc/2026-09-25-rubin-campaign-plan.md, D-139).
#
#   ssh nibi; cd ~/projects/def-snakagaw/snakagaw/pigauto_rubin
#   ARMSET=freq  N=1000 TIME=01:30:00 BLOCK=10 SEEDS=1-200 bash script/rubin_campaign_nibi.sh
#   ARMSET=bace  N=300  TIME=04:30:00 BLOCK=1  SEEDS=1-200 RUNS=5 NITT=50000 MEM=16G bash script/rubin_campaign_nibi.sh
#
# ARMSET=freq runs freqA,freqB; ARMSET=bace runs bace,bace_chain,bace_resid on one BACE fit. Each task runs
# BLOCK consecutive seeds of one cell. Seeds are shared across arm sets, so a (cell, seed) is the same simulated
# dataset for every arm: freq and BACE results pair by seed. nibi caps a user at 1,000 submitted array tasks:
# the script refuses a submission above THROTTLE_CAP (default 900) and prints the block size that fits.
# Resume is free: rubin_cell.R skips an existing rds, so re-submitting runs only what is missing. Size --time
# from seff on the first finished tasks, not from a guess.
set -euo pipefail

ARMSET="${ARMSET:?freq|bace}"; N="${N:?n}"; TIME="${TIME:?--time HH:MM:SS}"
BLOCK="${BLOCK:-1}"; SEEDS="${SEEDS:-1-200}"; MEM="${MEM:-8G}"; THROTTLE="${THROTTLE:-300}"
RUNS="${RUNS:-5}"; NITT="${NITT:-50000}"; CAP="${THROTTLE_CAP:-900}"
ROOT="${RUBIN_ROOT:-$HOME/projects/def-snakagaw/snakagaw/pigauto_rubin}"
ENV_SH="${RUBIN_ENV:-$HOME/projects/def-snakagaw/snakagaw/pigauto_sim/env.sh}"
# RUBIN_OUT: results parent folder (default results/); CELL_FLAGS: extra rubin_cell.R flags, e.g. the discrete
# re-run: RUBIN_OUT=$ROOT/results_disc CELL_FLAGS="--save_imp --discrete" TAG=disc
OUT="${RUBIN_OUT:-$ROOT/results}/$ARMSET"; LOG="$ROOT/logs"; mkdir -p "$OUT" "$LOG"
S0="${SEEDS%-*}"; S1="${SEEDS#*-}"

case "$ARMSET" in
  freq) ARMS="freqA,freqB"; EXTRA="" ;;
  bace) ARMS="bace,bace_chain,bace_resid"
        BURNIN=$(( NITT / 5 )); THIN=$(( (NITT - BURNIN) / 1600 ))
        EXTRA="--bace_nitt $NITT --bace_burnin $BURNIN --bace_thin $THIN --bace_runs $RUNS" ;;
  *) echo "ARMSET must be freq or bace"; exit 1 ;;
esac
EXTRA="$EXTRA ${CELL_FLAGS:-}"

TASKS="$LOG/campaign_${ARMSET}_n${N}_s${S0}-${S1}${TAG:+_$TAG}_tasks.txt"; : > "$TASKS"
# RETRY_FILE: explicit (lambda rho seed) triples, one per line, instead of the full grid (reruns of named fits).
if [ -n "${RETRY_FILE:-}" ]; then
  while read -r lambda rho seed; do
    printf '%s\t%d\t%d\n' "--n $N --lambda $lambda --rho $rho --M 20 --arms $ARMS $EXTRA --out $OUT" "$seed" "$seed" >> "$TASKS"
  done < "$RETRY_FILE"
else
for lambda in ${LAMBDAS:-0.3 0.7 1}; do for rho in 0 0.5; do
  for (( s = S0; s <= S1; s += BLOCK )); do
    e=$(( s + BLOCK - 1 )); [ "$e" -gt "$S1" ] && e=$S1
    printf '%s\t%d\t%d\n' "--n $N --lambda $lambda --rho $rho --M 20 --arms $ARMS $EXTRA --out $OUT" "$s" "$e" >> "$TASKS"
  done
done; done
fi
NT=$(wc -l < "$TASKS")
if [ "$NT" -gt "$CAP" ]; then
  echo "$NT tasks exceeds the $CAP cap; use BLOCK >= $(( ( (S1 - S0 + 1) * 6 + CAP - 1 ) / CAP )) or split SEEDS"; exit 1
fi
echo "armset=$ARMSET n=$N seeds=$S0..$S1 block=$BLOCK tasks=$NT time=$TIME mem=$MEM"
[ "${CONFIRM:-no}" = yes ] || { echo "dry run: set CONFIRM=yes to submit (needs Shinichi's approval)"; exit 0; }

SB="$LOG/campaign_${ARMSET}_n${N}_s${S0}-${S1}${TAG:+_$TAG}.sbatch"
cat > "$SB" <<SBEOF
#!/bin/bash
#SBATCH --account=def-snakagaw_cpu
#SBATCH --job-name=rubin_${ARMSET}_n${N}${TAG:+_$TAG}
#SBATCH --time=$TIME
#SBATCH --cpus-per-task=1
#SBATCH --mem=$MEM
#SBATCH --array=1-${NT}%${THROTTLE}
#SBATCH --output=$LOG/%x-%A_%a.out
set -euo pipefail
source "$ENV_SH"
cd "$ROOT"
line=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$TASKS")
args=\$(printf '%s' "\$line" | cut -f1); s0=\$(printf '%s' "\$line" | cut -f2); s1=\$(printf '%s' "\$line" | cut -f3)
echo "host=\$(hostname) task=\$SLURM_ARRAY_TASK_ID seeds=\$s0..\$s1"
for seed in \$(seq "\$s0" "\$s1"); do Rscript script/rubin_cell.R \$args --seed "\$seed"; done
SBEOF
sbatch "$SB"
