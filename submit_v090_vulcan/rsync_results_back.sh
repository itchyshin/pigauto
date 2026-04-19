#!/bin/bash
# rsync_results_back.sh
#
# Pull the per-cell RDS files (results/cell_*.rds) from Vulcan (PAICE)
# back to this machine, plus the SLURM stdout logs.
#
# Usage:
#   bash rsync_results_back.sh <CC_LOGIN>
#
# Example:
#   bash rsync_results_back.sh snakagaw@vulcan.alliancecan.ca
#
# Results land in submit_v090_vulcan/returned/ (created if absent).
# Each cell_<TASK_ID>.rds is a tidy data.frame; the local merge step
# rbind()s them into the grid-level RDS used by the pkgdown report.

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <CC_LOGIN>"
  echo "  e.g.: $0 snakagaw@vulcan.alliancecan.ca"
  exit 1
fi

CC_LOGIN="$1"
REMOTE_DIR="~/pigauto_vulcan"
LOCAL_DIR="$(dirname "$0")/returned"

mkdir -p "${LOCAL_DIR}"
mkdir -p "${LOCAL_DIR}/logs"

echo "Syncing per-cell results from ${CC_LOGIN}:${REMOTE_DIR}/results/"
echo "  -> ${LOCAL_DIR}/"

# Per-cell tidy RDS (one per array task, idempotent re-pull is fine)
rsync -avz --progress \
  "${CC_LOGIN}:${REMOTE_DIR}/results/cell_*.rds" \
  "${LOCAL_DIR}/"

# SLURM stdout logs (--output=pigauto-cal-%A_%a.out from submit script).
# Optional but handy when a task died mid-rep.
rsync -avz --progress \
  "${CC_LOGIN}:${REMOTE_DIR}/pigauto-cal-*.out" \
  "${LOCAL_DIR}/logs/" || true

echo ""
echo "Done. Per-cell RDS in ${LOCAL_DIR}/, logs in ${LOCAL_DIR}/logs/"
echo ""
echo "To merge into a grid-level data.frame:"
echo "  R -e 'files <- list.files(\"${LOCAL_DIR}\", pattern = \"^cell_.*\\\\.rds$\", full.names = TRUE);"
echo "        grid  <- do.call(rbind, lapply(files, readRDS));"
echo "        saveRDS(grid, \"${LOCAL_DIR}/calibration_grid.rds\");"
echo "        cat(sprintf(\"merged %d cells into calibration_grid.rds\\n\", length(files)))'"
