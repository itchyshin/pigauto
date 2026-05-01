#!/bin/bash
# rsync_results_back.sh
#
# Pull the AVONET 9,993 × BACE GPU bench results + log back from Vulcan.
#
# Usage:
#   bash rsync_results_back.sh <CC_LOGIN>
#
# Example:
#   bash rsync_results_back.sh snakagaw@vulcan.alliancecan.ca
#
# Pulls:
#   - bench_avonet9993_bace.rds / .md             (tidy + human-readable)
#   - results-avonet9993-*.tar.gz                 (archived log + output)
#   - run_avonet9993_bace_*.log                   (SLURM stdout)
#
# Lands in submit_v090_vulcan_gpu/returned/.  Also copies the .rds
# to script/ so the local HTML generator can pick it up.

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <CC_LOGIN>"
  echo "  e.g.: $0 snakagaw@vulcan.alliancecan.ca"
  exit 1
fi

CC_LOGIN="$1"
REMOTE_DIR="~/pigauto_vulcan_gpu"
BUNDLE_DIR="$(dirname "$0")"
LOCAL_DIR="${BUNDLE_DIR}/returned"
SCRIPT_DIR="$(dirname "$(dirname "$(realpath "$0")")")/script"

mkdir -p "${LOCAL_DIR}"

echo "Syncing bench output from ${CC_LOGIN}:${REMOTE_DIR}/"
echo "  -> ${LOCAL_DIR}/"

rsync -avz --progress \
  "${CC_LOGIN}:${REMOTE_DIR}/bench_avonet9993_bace.rds" \
  "${CC_LOGIN}:${REMOTE_DIR}/bench_avonet9993_bace.md" \
  "${CC_LOGIN}:${REMOTE_DIR}/results-avonet9993-*.tar.gz" \
  "${CC_LOGIN}:${REMOTE_DIR}/run_avonet9993_bace_*.log" \
  "${LOCAL_DIR}/" 2>/dev/null || true

# Mirror the RDS into script/ so make_bench_avonet9993_bace_html.R can
# find it at the conventional location (script/bench_*.rds).
if [[ -f "${LOCAL_DIR}/bench_avonet9993_bace.rds" ]]; then
  cp "${LOCAL_DIR}/bench_avonet9993_bace.rds" "${SCRIPT_DIR}/"
  echo "Copied bench_avonet9993_bace.rds -> ${SCRIPT_DIR}/"
fi

echo ""
echo "Done. Files in ${LOCAL_DIR}/"
echo ""
echo "To render the HTML report:"
echo "  Rscript script/make_bench_avonet9993_bace_html.R"
