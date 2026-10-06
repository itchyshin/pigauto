#!/bin/bash
# Mondrian real-data campaign on DRAC fir: one (dataset, arm, seed) cell per
# array task, Mondrian then split, one compute thread. Cells per the
# pre-registration and Amendment 1. Set --time from the measured pre-run
# (docs/dev-log/mondrian-realdata/RUNLOG.md) at submission: sbatch --time=HH:MM:SS
#SBATCH --account=def-snakagaw
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=1-9
#SBATCH --output=/scratch/snakagaw/pigauto-mondrian-rd/logs/cell-%A_%a.out
set -uo pipefail
module load StdEnv/2023 gcc/12.3 cuda/12.6 r/4.4.0
export LD_LIBRARY_PATH="$EBROOTCUDA/targets/x86_64-linux/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
P=/scratch/snakagaw/pigauto-mondrian-rd
export R_LIBS=$P/lib:$HOME/R/x86_64-pc-linux-gnu-library/4.4
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CELLS=(
  "pantheria mcar 20260818" "pantheria mcar 20260819" "pantheria mcar 20260820"
  "pantheria structured 20260818" "pantheria structured 20260819" "pantheria structured 20260820"
  "avonet mcar 20260818" "avonet mcar 20260819" "avonet mcar 20260820"
)
set -- ${CELLS[$((SLURM_ARRAY_TASK_ID - 1))]}
DS=$1; ARM=$2; SEED=$3
OUT=$P/results/$DS-$ARM-m$SEED
mkdir -p "$OUT"
tar xzf "$P/pigauto_7af133a.tar.gz" -C "$SLURM_TMPDIR"
cd "$P"
Rscript --vanilla "$SLURM_TMPDIR"/pigauto/script/mondrian_confirmation/01_run_masked_confirmation.R \
  "input/$DS-full-input.rds" "$OUT" "$ARM" "$SEED" 500 mondrian,split
