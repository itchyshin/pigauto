#!/usr/bin/env bash
# nibi job-array driver for the BACE settings pre-run of the rubin-freq-bace lane
# (docs/dev-log/arc/2026-09-24-rubin-prerun-plan.md). Used instead of Totoro on 2026-09-24 because this
# user's other lanes already held ~200 Totoro cores, over D-143's 150.
#
#   ssh nibi; cd ~/projects/def-snakagaw/snakagaw/pigauto_rubin
#   bash script/rubin_prerun_nibi.sh check        # one short validation job (gate tests + a BACE smoke)
#   bash script/rubin_prerun_nibi.sh submit       # the four pre-run arrays
#
# One array task = one BACE fit (arms bace, bace_chain, bace_resid share it). Arrays are split by
# (n, nitt) so each --time is sized to its own slowest setting: cost = nitt x (runs + 40) units at
# 0.93 ms/unit (n = 100) and 2.44 ms/unit (n = 300), rates from v1's nibi BACE fits, plus ~70% margin.
# Resume is free: rubin_cell.R skips a cell whose rds exists. This script only writes and submits sbatch
# files; it computes nothing on the login node.
set -euo pipefail

MODE="${1:?check|submit}"
ROOT="${RUBIN_ROOT:-$HOME/projects/def-snakagaw/snakagaw/pigauto_rubin}"
ENV_SH="${RUBIN_ENV:-$HOME/projects/def-snakagaw/snakagaw/pigauto_sim/env.sh}"
OUT="$ROOT/prerun"; LOG="$ROOT/logs"; mkdir -p "$OUT" "$LOG"
MEM="${MEM:-16G}"; THROTTLE="${THROTTLE:-60}"

for f in campaign_gnn_off_lib.R rubin_lib.R rubin_freq.R rubin_bace.R rubin_cell.R; do
  [ -s "$ROOT/script/$f" ] || { echo "missing $ROOT/script/$f (rsync the lane's script/ first)"; exit 1; }
done

if [ "$MODE" = check ]; then
  SB="$LOG/prerun_check.sbatch"
  cat > "$SB" <<EOF
#!/bin/bash
#SBATCH --account=def-snakagaw_cpu
#SBATCH --job-name=rubin_check
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=$LOG/%x-%j.out
set -euo pipefail
source "$ENV_SH"
cd "$ROOT"
Rscript -e 'stopifnot(any(grepl("sample = TRUE", deparse(BACE:::bace_final_imp), fixed = TRUE))); cat("BACE build OK\n")'
# testthat is not in the shared nibi library (v1's env); the unit tests already pass on the Mac against
# byte-identical BACE internals, so they run here only when testthat is present.
if Rscript -e 'quit(status = !requireNamespace("testthat", quietly = TRUE))'; then
for t in test-bace.R test-bace-chain.R test-lib-pool.R test-freq-condmean.R; do
  Rscript -e "r <- as.data.frame(testthat::test_file('script/tests-rubin/\$t', reporter = 'silent')); stopifnot(nrow(r) > 0, sum(r\\\$failed) == 0, !any(r\\\$error)); cat('\$t PASS\n')"
done
else echo "testthat absent: unit tests skipped (they pass on the Mac)"; fi
Rscript script/rubin_cell.R --n 60 --seed 2 --lambda 0.7 --rho 0.5 --M 20 --smoke --out "$LOG/check_smoke"
Rscript script/rubin_checks.R --gate G-S4a --dir "$LOG/check_smoke"
EOF
  sbatch "$SB"; exit 0
fi

[ "$MODE" = submit ] || { echo "mode must be check or submit"; exit 1; }
# group: n nitt --time
while read -r N NITT TIME; do
  TASKS="$LOG/prerun_n${N}_nitt${NITT}_tasks.txt"; : > "$TASKS"
  burnin=$(( NITT / 5 )); thin=$(( (NITT - burnin) / 1600 ))
  for runs in 15 10 5; do for lambda in 0.3 0.7; do for rho in 0 0.5; do for seed in 1 2 3 4 5; do
    echo "--n $N --seed $seed --lambda $lambda --rho $rho --M 20 --arms bace,bace_chain,bace_resid --bace_nitt $NITT --bace_burnin $burnin --bace_thin $thin --bace_runs $runs --out $OUT/runs${runs}_nitt${NITT}" >> "$TASKS"
  done; done; done; done
  NT=$(wc -l < "$TASKS")
  SB="$LOG/prerun_n${N}_nitt${NITT}.sbatch"
  cat > "$SB" <<EOF
#!/bin/bash
#SBATCH --account=def-snakagaw_cpu
#SBATCH --job-name=rubin_pre_n${N}_${NITT}
#SBATCH --time=$TIME
#SBATCH --cpus-per-task=1
#SBATCH --mem=$MEM
#SBATCH --array=1-${NT}%${THROTTLE}
#SBATCH --output=$LOG/%x-%A_%a.out
set -euo pipefail
source "$ENV_SH"
cd "$ROOT"
args=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$TASKS")
echo "host=\$(hostname) task=\$SLURM_ARRAY_TASK_ID args=\$args"
Rscript script/rubin_cell.R \$args
EOF
  echo "n=$N nitt=$NITT tasks=$NT time=$TIME"; sbatch "$SB"
done <<'GROUPS'
100 50000 01:15:00
100 100000 02:30:00
300 50000 03:15:00
300 100000 06:30:00
GROUPS
