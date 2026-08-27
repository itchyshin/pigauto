#!/bin/bash
# Passive monitor: poll Phase B 3a, 3b, and AVONET panel on Totoro every 10 min.
# Usage:
#   nohup bash script/monitor_gnn_evidence_totoro.sh >> script/returned_gnn_campaign/monitor.log 2>&1 &
#   bash script/monitor_gnn_evidence_totoro.sh once
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOG="${ROOT}/script/returned_gnn_campaign/monitor.log"
CHECKPOINT="${ROOT}/LOOP/checkpoint.md"
INTERVAL="${MONITOR_INTERVAL_SEC:-600}"
MAX_POLLS="${MONITOR_MAX_POLLS:-0}"

poll_once() {
  local wall done3a done3b doneav all_done
  wall=$(ssh totoro '
    echo "=== $(date -Is) load=$(cut -d" " -f1-3 /proc/loadavg) ==="
    bash ~/pigauto_gnn_evidence_phase_b/script/gnn_evidence_phase_b_phylo_cov_mar_totoro.sh status 2>/dev/null | head -3
    bash ~/pigauto_gnn_evidence_phase_b_mnar/script/gnn_evidence_phase_b_mnar_totoro.sh status 2>/dev/null | head -2
    bash ~/pigauto_gnn_evidence_campaign/script/gnn_evidence_avonet_panel_totoro.sh status 2>/dev/null | head -4
  ' 2>&1)
  echo "${wall}"
  echo "${wall}" >> "${LOG}"

  done3a=$(echo "${wall}" | grep -oE '[0-9]+/4860' | head -1 | cut -d/ -f1)
  done3b=$(echo "${wall}" | grep -oE '[0-9]+/2430' | head -1 | cut -d/ -f1)
  doneav=$(echo "${wall}" | grep -oE '[0-9]+/15 RDS' | head -1 | cut -d/ -f1)
  done3a=${done3a:-0}; done3b=${done3b:-0}; doneav=${doneav:-0}

  if [[ -f "${CHECKPOINT}" ]]; then
    ts=$(date '+%Y-%m-%d %H:%M %Z')
    sed -i.bak \
      -e "s/| Phase B lane 3b (MNAR) | \*\*[^|]*\*\* — [0-9]*\/2430 RDS @ .*/| Phase B lane 3b (MNAR) | **RUNNING** — ${done3b}\/2430 RDS @ ${ts} |/" \
      -e "s/| Phase B lane 3a (phylo\/cov MAR) | \*\*[^|]*\*\* — [0-9]*\/4860 RDS @ .*/| Phase B lane 3a (phylo\/cov MAR) | **RUNNING** — ${done3a}\/4860 RDS @ ${ts} |/" \
      -e "s/| AVONET panel | \*\*[^|]*\*\* — .*/| AVONET panel | **RUNNING** — ${doneav}\/15 RDS @ ${ts} |/" \
      "${CHECKPOINT}" 2>/dev/null || true
    rm -f "${CHECKPOINT}.bak"
  fi

  all_done=0
  [[ "${done3a}" -ge 4860 && "${done3b}" -ge 2430 && "${doneav}" -ge 15 ]] && all_done=1
  return $((1 - all_done))
}

case "${1:-loop}" in
  once) poll_once; exit $? ;;
  loop)
    n=0
    while true; do
      n=$((n + 1))
      if poll_once; then
        echo "monitor exit — all jobs complete at $(date -Is)" | tee -a "${LOG}"
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
    echo "usage: $0 [loop|once]" >&2
    exit 1
    ;;
esac
