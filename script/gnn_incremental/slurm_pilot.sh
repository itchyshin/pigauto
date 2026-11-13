#!/usr/bin/env bash
# Submit through sbatch on Tamia. This script itself runs only in a job.
# Usage: sbatch --export=ALL,PIGAUTO_SOURCE_BUNDLE=sha256:... slurm_pilot.sh
#        RUN_DIR
# The pilot is intentionally five datasets (one per frozen scenario), with
# both gate settings per dataset. It is not the production campaign.

#SBATCH --account=aip-snakagaw
#SBATCH --partition=cpubase_bycore_b1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --job-name=pigauto-gnn-pilot
#SBATCH --output=pilot-%j.out
#SBATCH --error=pilot-%j.err

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: slurm_pilot.sh RUN_DIR" >&2
  exit 2
fi

source_dir=/project/aip-snakagaw/pigauto-gnn-validation/source
library_dir=/project/aip-snakagaw/pigauto-gnn-validation/Rlib
base_library=/project/aip-snakagaw/pigauto-stagec-gpu/Rlib
torch_home=/project/aip-snakagaw/pigauto-gnn-validation/TORCH_HOME
run_dir=$1

module load StdEnv/2023 gcc/12.3 r/4.5.0
export R_LIBS_USER="${library_dir}:${base_library}"
export TORCH_HOME="$torch_home" CUDA=cpu TORCH_VERIFY_LOAD=FALSE
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PIGAUTO_SOURCE_BUNDLE="${PIGAUTO_SOURCE_BUNDLE:?source bundle identity required}"
export OMP_THREAD_LIMIT=1

cd "$source_dir"
R CMD INSTALL --library="$library_dir" .
Rscript --vanilla script/gnn_incremental/campaign.R init "$run_dir" 1 20261019 2000
bash script/gnn_incremental/run_tasks.sh "$run_dir" 4
