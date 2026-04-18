#!/bin/bash
# rsync_results_back.sh
#
# Pull results archives from Compute Canada back to this machine.
# Usage:
#   bash rsync_results_back.sh <CC_LOGIN>
#
# Example:
#   bash rsync_results_back.sh snakagawa@fir.alliancecan.ca
#
# The archives land in submit_v090_cloud/returned/ (created if absent).
# Unpack them there:
#   cd returned && for f in *.tar.gz; do tar xzf "$f"; done

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <CC_LOGIN>"
  echo "  e.g.: $0 snakagawa@fir.alliancecan.ca"
  exit 1
fi

CC_LOGIN="$1"
REMOTE_DIR="~/pigauto_cloud"
LOCAL_DIR="$(dirname "$0")/returned"

mkdir -p "${LOCAL_DIR}"

echo "Syncing results from ${CC_LOGIN}:${REMOTE_DIR}/results-*.tar.gz"
echo "  -> ${LOCAL_DIR}/"

rsync -avz --progress \
  "${CC_LOGIN}:${REMOTE_DIR}/results-*.tar.gz" \
  "${LOCAL_DIR}/"

echo ""
echo "Done. Archives in ${LOCAL_DIR}/"
echo "To unpack:"
echo "  cd ${LOCAL_DIR} && for f in *.tar.gz; do tar xzf \"\$f\"; done"
