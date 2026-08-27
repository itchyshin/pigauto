#!/bin/bash
# Push/pull GNN evidence Phase B lane 3b (MNAR arm) to Totoro.
# Usage:
#   bash script/rsync_gnn_evidence_phase_b_mnar.sh push
#   bash script/rsync_gnn_evidence_phase_b_mnar.sh pull
set -euo pipefail

MODE="${1:-push}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REMOTE="totoro:~/pigauto_gnn_evidence_phase_b_mnar"

RSYNC_EXCLUDES=(
  --exclude '.git'
  --exclude 'BACE'
  --exclude 'dev/gnn_attribution_*'
  --exclude 'script/returned_gnn_attr'
  --exclude 'script/returned_gnn_sentinel'
  --exclude 'script/returned_gnn_campaign'
  --exclude '.unlazy'
  --exclude 'checkpoints*'
  --exclude 'useful'
  --exclude 'avonet'
)

case "${MODE}" in
  push)
    echo "Syncing ${ROOT} -> ${REMOTE}"
    rsync -avz "${RSYNC_EXCLUDES[@]}" \
      --exclude 'results_phase_b_mnar' --exclude 'logs' \
      "${ROOT}/" "${REMOTE}/"
    ;;
  pull)
    LOCAL="${ROOT}/script/returned_gnn_campaign/"
    mkdir -p "${LOCAL}"
    echo "Syncing ${REMOTE}/results_phase_b_mnar + logs -> ${LOCAL}"
    rsync -avz "${REMOTE}/results_phase_b_mnar/" "${LOCAL}results_phase_b_mnar/" || true
    rsync -avz "${REMOTE}/logs/gnn_phase_b_mnar_job_*.log" "${LOCAL}logs/" 2>/dev/null || true
    rsync -avz "${REMOTE}/logs/phase_b_mnar.log" "${LOCAL}logs/" 2>/dev/null || true
    ;;
  *)
    echo "usage: $0 push|pull" >&2
    exit 1
    ;;
esac
