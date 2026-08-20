#!/bin/bash
# Bounded Stage-C feasibility ladder: FishBase 300/600/1000/1500 tips,
# split only. Tamia allocates H100s by whole four-GPU nodes, so this uses all
# four GPUs for independent receipts rather than reserving idle devices. It
# measures device availability and wall-time scaling; it is not calibration
# evidence.
#SBATCH --account=aip-snakagaw
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=h100:4
#SBATCH --mem=128G
#SBATCH --time=01:00:00

set -euo pipefail
PIGAUTO_STAGEC_ROOT="${PIGAUTO_STAGEC_ROOT:-/project/aip-snakagaw/pigauto-stagec-gpu}"
module load r/4.5.0 cuda/12.6
export LD_LIBRARY_PATH="$EBROOTCUDA/targets/x86_64-linux/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="$PIGAUTO_STAGEC_ROOT/Rlib"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mkdir -p "$PIGAUTO_STAGEC_ROOT/results"

status=0
for n_tips in 300 600 1000 1500; do
  # The outer allocation is typed (`h100`); child requests must preserve that
  # type.  An untyped --gpus-per-task request is rejected by Tamia Slurm.
  srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=12 \
    --gres=gpu:h100:1 --gpu-bind=single:1 \
    "$PIGAUTO_STAGEC_ROOT/source/pigauto/script/mondrian_confirmation/05_run_tamia_ladder_task.sh" \
    "$PIGAUTO_STAGEC_ROOT" "$n_tips" "$SLURM_JOB_ID" &
done
for pid in $(jobs -p); do wait "$pid" || status=1; done
exit "$status"
