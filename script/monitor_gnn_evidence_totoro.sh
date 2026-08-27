#!/bin/bash
# Passive monitor: poll Phase B lanes 3a+3b on Totoro every 10 min.
# On completion: pull, collect, PHASE_B_SUMMARY.md, checkpoint, commit (no push).
#
# Usage:
#   nohup bash script/monitor_gnn_evidence_totoro.sh >> script/returned_gnn_campaign/monitor.log 2>&1 &
#   bash script/monitor_gnn_evidence_totoro.sh once
#   bash script/monitor_gnn_evidence_totoro.sh finalize   # force collect if already complete
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOG="${ROOT}/script/returned_gnn_campaign/monitor.log"
CHECKPOINT="${ROOT}/LOOP/checkpoint.md"
INTERVAL="${MONITOR_INTERVAL_SEC:-600}"
MAX_POLLS="${MONITOR_MAX_POLLS:-0}"
TARGET_3A=4860
TARGET_3B=2430

poll_once() {
  local wall done3a done3b phase_b_done
  wall=$(ssh totoro '
    echo "=== $(date -Is) load=$(cut -d" " -f1-3 /proc/loadavg) ==="
    bash ~/pigauto_gnn_evidence_phase_b/script/gnn_evidence_phase_b_phylo_cov_mar_totoro.sh status 2>/dev/null | head -3
    bash ~/pigauto_gnn_evidence_phase_b_mnar/script/gnn_evidence_phase_b_mnar_totoro.sh status 2>/dev/null | head -2
  ' 2>&1)
  echo "${wall}"
  echo "${wall}" >> "${LOG}"

  done3a=$(echo "${wall}" | grep -oE '[0-9]+/4860' | head -1 | cut -d/ -f1)
  done3b=$(echo "${wall}" | grep -oE '[0-9]+/2430' | head -1 | cut -d/ -f1)
  done3a=${done3a:-0}; done3b=${done3b:-0}

  if [[ -f "${CHECKPOINT}" ]]; then
    ts=$(date '+%Y-%m-%d %H:%M %Z')
    sed -i.bak \
      -e "s/| Phase B lane 3b (MNAR) | \*\*[^|]*\*\* — [0-9]*\/2430 RDS @ .*/| Phase B lane 3b (MNAR) | **RUNNING** — ${done3b}\/2430 RDS @ ${ts} |/" \
      -e "s/| Phase B lane 3a (phylo\/cov MAR) | \*\*[^|]*\*\* — [0-9]*\/4860 RDS @ .*/| Phase B lane 3a (phylo\/cov MAR) | **RUNNING** — ${done3a}\/4860 RDS @ ${ts} |/" \
      "${CHECKPOINT}" 2>/dev/null || true
    rm -f "${CHECKPOINT}.bak"
  fi

  phase_b_done=0
  [[ "${done3a}" -ge "${TARGET_3A}" && "${done3b}" -ge "${TARGET_3B}" ]] && phase_b_done=1
  echo "[monitor] Phase B: 3a=${done3a}/${TARGET_3A} 3b=${done3b}/${TARGET_3B} done=${phase_b_done}" | tee -a "${LOG}"
  return $((1 - phase_b_done))
}

on_complete() {
  echo "[monitor] Phase B complete — running finalize at $(date -Is)" | tee -a "${LOG}"
  bash "${ROOT}/script/finalize_gnn_evidence_phase_b.sh" 2>&1 | tee -a "${LOG}"
}

case "${1:-loop}" in
  once)
    if poll_once; then
      on_complete
      exit 0
    fi
    exit 1
    ;;
  finalize)
    on_complete
    exit 0
    ;;
  loop)
    n=0
    while true; do
      n=$((n + 1))
      if poll_once; then
        echo "monitor exit — Phase B complete at $(date -Is)" | tee -a "${LOG}"
        on_complete
        exit 0
      fi
      if [[ "${MAX_POLLS}" -gt 0 && "${n}" -ge "${MAX_POLLS}" ]]; then
        echo "monitor exit — max polls ${MAX_POLLS} at $(date -Is)" | tee -a "${LOG}"
        exit 0
      fi
      sleep "${INTERVAL}"
    done
    ;;
  *)
    echo "usage: $0 [loop|once|finalize]" >&2
    exit 1
    ;;
esac
