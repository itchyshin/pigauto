#!/usr/bin/env bash
# DRAC job-array driver for the factorial stage of the imputation simulation (nibi, rorqual, fir: same script).
#
#   ssh nibi; cd ~/projects/def-snakagaw/snakagaw/pigauto_sim
#   bash script/campaign_sim_nibi_array.sh factorial 100  00:50:00     # n = 100 cells
#   bash script/campaign_sim_nibi_array.sh factorial 1000 04:30:00     # n = 1000 cells
#   HALF=A|B splits the design cells alternately across two clusters (A = odd rows, B = even rows);
#   unset = every cell. Finished (cell, seed) rds are skipped, so any cluster can take over another's remainder.
#
# One array task = one cell x BLOCK consecutive seeds (default 5). Arrays are per n so --time can be
# sized from the pre-run walls (seff on the first finished task, then re-submit the rest with the
# measured time + 30%). %THROTTLE caps concurrent tasks. Resume is free: finished (cell, seed) rds are
# skipped by campaign_sim_cell.R, so re-submitting the same array only runs what is missing.
# Never run this on the login node's shell for compute: it only writes and submits the sbatch script.
set -euo pipefail

STAGE="${1:?stage}"; N="${2:?n}"; TIME="${3:?--time HH:MM:SS}"
BLOCK="${BLOCK:-5}"; HALF="${HALF:-}"; THROTTLE="${THROTTLE:-400}"
# Measured on fir 2026-09-20 (seff over a completed BACE n = 1000 task): CPU efficiency 24.9% of a
# 4-core allocation, memory efficiency 92% of 16 GB with 200 of 600 tasks killed OUT_OF_MEMORY.
# BACE is single-threaded and memory-hungry, so a Bayesian-only array asks for 1 core and 32 GB;
# the GNN arm is the only one that uses several threads. Override with CPUS= / MEM=.
CPUS="${CPUS:-4}"; MEM="${MEM:-16G}"
ROOT="${PIG_SIM_ROOT:-$HOME/projects/def-snakagaw/snakagaw/pigauto_sim}"
OUT="$ROOT/results/$STAGE"; LOG="$ROOT/logs"; mkdir -p "$OUT" "$LOG"

# ARMS / ARMS_BACE override the arm sets, as on Totoro. SEEDS=bace caps the seed range at the
# design's bace_reps so a Bayesian-only wave does not run seeds the other arms own.
ARMS_ALL="${ARMS_BACE:-gnn_on,gnn_off,gnn_off_rphylopars,freq,bace,floor}"
ARMS_NOBACE="${ARMS:-gnn_on,gnn_off,gnn_off_rphylopars,freq,floor}"

# Task table: one line per (cell, seed block). Column 1 = the cell args, column 2 = first seed, 3 = last.
TASKS="$LOG/${STAGE}_n${N}${HALF:+_half$HALF}_tasks.txt"
source "$ROOT/env.sh"
Rscript "$ROOT/script/campaign_sim_design.R" --stage "$STAGE" | awk -F, -v n="$N" -v blk="$BLOCK" -v half="$HALF" -v seeds="${SEEDS:-all}" '
NR > 1 && $9 == n && (half == "" || (half == "A" && (NR % 2) == 0) || (half == "B" && (NR % 2) == 1)) {
  lim = (seeds == "bace") ? $11 : $10
  for (s = 1; s <= lim; s += blk) {
    e = s + blk - 1; if (e > lim) e = lim
    printf "--dgp %s --evo %s --lambda %s --rho %s --miss %s --frac %s --n %s --driver --thresholds fixed\t%d\t%d\t%d\n", $3, $4, $5, $6, $7, $8, $9, s, e, $11
  }
}' > "$TASKS"
NT=$(wc -l < "$TASKS")
echo "stage=$STAGE n=$N half=${HALF:-all} tasks=$NT block=$BLOCK time=$TIME throttle=$THROTTLE"

SB="$LOG/${STAGE}_n${N}${HALF:+_half$HALF}.sbatch"
cat > "$SB" <<EOF
#!/bin/bash
#SBATCH --account=def-snakagaw_cpu
#SBATCH --job-name=pig_${STAGE}_n${N}${HALF}
#SBATCH --time=$TIME
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --array=1-${NT}%${THROTTLE}
#SBATCH --output=$LOG/%x-%A_%a.out
set -euo pipefail
source "$ROOT/env.sh"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PIG_TORCH_THREADS=$CPUS NOT_CRAN=true
export PIG_BACE_RUNS=${PIG_BACE_RUNS:-5} PIG_BACE_NFINAL=${PIG_BACE_NFINAL:-20}
cd "$ROOT"
line=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$TASKS")
cellargs=\$(printf '%s' "\$line" | cut -f1); s0=\$(printf '%s' "\$line" | cut -f2); s1=\$(printf '%s' "\$line" | cut -f3); nb=\$(printf '%s' "\$line" | cut -f4)
echo "host=\$(hostname) task=\$SLURM_ARRAY_TASK_ID seeds=\$s0..\$s1"
for seed in \$(seq "\$s0" "\$s1"); do
  if [ "\$seed" -le "\$nb" ]; then arms="$ARMS_ALL"; else arms="$ARMS_NOBACE"; fi
  Rscript script/campaign_sim_cell.R \$cellargs --seed "\$seed" --arms "\$arms" --out "$OUT"
done
EOF
sbatch "$SB"
