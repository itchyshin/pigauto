#!/bin/bash
# Push evidence lane to Totoro and optionally pull results back.
# Usage:
#   bash script/rsync_gnn_evidence_sentinel_prerun.sh push
#   bash script/rsync_gnn_evidence_sentinel_prerun.sh pull
set -euo pipefail

MODE="${1:-push}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REMOTE="totoro:~/pigauto_gnn_sentinel_prerun"

RSYNC_EXCLUDES=(
  --exclude '.git'
  --exclude 'BACE'
  --exclude 'dev/gnn_attribution_*'
  --exclude 'script/returned_gnn_attr'
  --exclude '.unlazy'
  --exclude 'checkpoints*'
  --exclude 'useful'
  --exclude 'avonet'
)

case "${MODE}" in
  push)
    echo "Syncing ${ROOT} -> ${REMOTE}"
    rsync -avz "${RSYNC_EXCLUDES[@]}" \
      --exclude 'results' --exclude 'logs' \
      "${ROOT}/" "${REMOTE}/"
    ;;
  pull)
    LOCAL="${ROOT}/script/returned_gnn_sentinel/"
    mkdir -p "${LOCAL}"
    echo "Syncing ${REMOTE}/results + logs -> ${LOCAL}"
    rsync -avz "${REMOTE}/results/" "${LOCAL}results/" || true
    rsync -avz "${REMOTE}/logs/" "${LOCAL}logs/" || true
    ;;
  *)
    echo "usage: $0 push|pull" >&2
    exit 1
    ;;
esac
