#!/usr/bin/env bash
# Launch the frozen 2,500-dataset production campaign on Totoro.
# Requires explicit user approval and a freshly checked live capacity.

set -euo pipefail

if [[ "${PIGAUTO_PRODUCTION_APPROVED:-}" != "YES" ]]; then
  echo "Production is gated. Set PIGAUTO_PRODUCTION_APPROVED=YES only after explicit user approval." >&2
  exit 3
fi

if [[ $# -ne 1 ]]; then
  echo "Usage: PIGAUTO_PRODUCTION_APPROVED=YES PIGAUTO_SOURCE_BUNDLE=sha256:... PIGAUTO_CAMPAIGN_SOURCE_DIR=/path/to/source $0 RUN_DIR" >&2
  exit 2
fi
if [[ "${PIGAUTO_TOTORO_CAPACITY_CHECKED:-}" != "YES" ]]; then
  echo "Totoro production is gated until a fresh live capacity and CPU-owner check is recorded." >&2
  exit 3
fi

source_dir=${PIGAUTO_CAMPAIGN_SOURCE_DIR:?exact production source directory is required}
source_bundle=${PIGAUTO_SOURCE_BUNDLE:?exact source-bundle hash is required}
library_dir=${PIGAUTO_R_LIBRARY:-/home/snakagaw/R/pigauto-gnn-validation/Rlib}
base_library=${PIGAUTO_BASE_R_LIBRARY:-/home/snakagaw/R/lib}
run_dir=$1
workers=${PIGAUTO_WORKERS:-48}

if [[ "$(hostname -s)" != "totoro" ]]; then
  echo "This launcher is for Totoro only; use the approved DRAC route if Totoro is unavailable." >&2
  exit 2
fi
if [[ ! "$workers" =~ ^[1-9][0-9]*$ || "$workers" -gt 48 ]]; then
  echo "PIGAUTO_WORKERS must be an integer from 1 to 48." >&2
  exit 2
fi
if [[ ! "$run_dir" =~ ^/home/snakagaw/hsq_work/ ]]; then
  echo "RUN_DIR must be under /home/snakagaw/hsq_work (persistent storage)." >&2
  exit 2
fi
if [[ -e "$run_dir" ]]; then
  echo "RUN_DIR already exists; refusing to overwrite or reuse a run." >&2
  exit 2
fi
if [[ ! "$source_bundle" =~ ^sha256:[0-9a-f]{64}$ ]]; then
  echo "PIGAUTO_SOURCE_BUNDLE must be sha256:<64 lowercase hex characters>." >&2
  exit 2
fi
if [[ ! -f "$source_dir/script/gnn_incremental/campaign.R" ||
      ! -f "$source_dir/script/gnn_incremental/scenario_manifest.csv" ]]; then
  echo "The source directory lacks the frozen campaign runner or scenario manifest." >&2
  exit 2
fi

export R_LIBS_USER="${library_dir}:${base_library}"
export TORCH_VERIFY_LOAD=FALSE
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 OMP_THREAD_LIMIT=1
export PIGAUTO_SOURCE_BUNDLE="$source_bundle"

cd "$source_dir"
R CMD INSTALL --library="$library_dir" .
Rscript --vanilla script/gnn_incremental/campaign.R init "$run_dir" 500 20261019 2000
bash script/gnn_incremental/run_tasks.sh "$run_dir" "$workers"
