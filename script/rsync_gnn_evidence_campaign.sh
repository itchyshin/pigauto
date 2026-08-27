#!/bin/bash
# Push/pull GNN evidence Phase A campaign to Totoro.
# Usage:
#   bash script/rsync_gnn_evidence_campaign.sh push
#   bash script/rsync_gnn_evidence_campaign.sh pull
set -euo pipefail

MODE="${1:-push}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REMOTE="totoro:~/pigauto_gnn_evidence_campaign"

RSYNC_EXCLUDES=(
  --exclude '.git'
  --exclude 'BACE'
  --exclude 'dev/gnn_attribution_*'
  --exclude 'script/returned_gnn_attr'
  --exclude 'script/returned_gnn_sentinel'
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
      --exclude 'results_bayes' --exclude 'results_confirm' \
      --exclude 'results_avonet_panel' \
      "${ROOT}/" "${REMOTE}/"
    if [[ -d "${ROOT}/BACE" ]]; then
      echo "Syncing BACE/ for panel comparator"
      rsync -avz "${ROOT}/BACE/" "${REMOTE}/BACE/"
    fi
    ;;
  pull)
    LOCAL="${ROOT}/script/returned_gnn_campaign/"
    mkdir -p "${LOCAL}"
    echo "Syncing ${REMOTE}/results + logs -> ${LOCAL}"
    rsync -avz "${REMOTE}/results/" "${LOCAL}results/" || true
    rsync -avz "${REMOTE}/logs/" "${LOCAL}logs/" || true
    rsync -avz "${REMOTE}/results_bayes/" "${LOCAL}results_bayes/" 2>/dev/null || true
    rsync -avz "${REMOTE}/results_confirm/" "${LOCAL}results_confirm/" 2>/dev/null || true
    rsync -avz "${REMOTE}/results_avonet_panel/" "${LOCAL}results_avonet_panel/" 2>/dev/null || true
    rsync -avz "${REMOTE}/results_avonet_panel/AVONET_PANEL_SUMMARY.md" \
      "${LOCAL}AVONET_PANEL_SUMMARY.md" 2>/dev/null || true
    ;;
  *)
    echo "usage: $0 push|pull" >&2
    exit 1
    ;;
esac
