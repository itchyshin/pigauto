#!/usr/bin/env bash
# Run the bounded five-dataset pilot on Totoro after a live capacity check.
# Usage: PIGAUTO_SOURCE_BUNDLE=sha256:... bash totoro_pilot.sh RUN_DIR

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: totoro_pilot.sh RUN_DIR" >&2
  exit 2
fi

source_dir=${PIGAUTO_CAMPAIGN_SOURCE_DIR:-/home/snakagaw/hsq_work/pigauto-gnn-incremental-20261019-corrected/source}
library_dir=${PIGAUTO_R_LIBRARY:-/home/snakagaw/R/pigauto-gnn-validation/Rlib}
base_library=${PIGAUTO_BASE_R_LIBRARY:-/home/snakagaw/R/lib}
run_dir=$1

export R_LIBS_USER="${library_dir}:${base_library}"
export TORCH_VERIFY_LOAD=FALSE
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 OMP_THREAD_LIMIT=1
export PIGAUTO_SOURCE_BUNDLE="${PIGAUTO_SOURCE_BUNDLE:?source bundle identity required}"

cd "$source_dir"
R CMD INSTALL --library="$library_dir" .
Rscript --vanilla script/gnn_incremental/campaign.R init "$run_dir" 1 20261019 2000
bash script/gnn_incremental/run_tasks.sh "$run_dir" 4
