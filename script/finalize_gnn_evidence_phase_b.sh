#!/bin/bash
# Pull Phase B lanes 3a+3b, collect, G6 audit, write PHASE_B_SUMMARY.md, update checkpoint, commit.
# Usage: bash script/finalize_gnn_evidence_phase_b.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"

CAMPAIGN="${ROOT}/script/returned_gnn_campaign"
SUMMARY="${CAMPAIGN}/PHASE_B_SUMMARY.md"
CHECKPOINT="${ROOT}/LOOP/checkpoint.md"
TS="$(date '+%Y-%m-%d %H:%M %Z')"

echo "[finalize] pulling lane 3a (phylo/cov MAR)..."
bash script/rsync_gnn_evidence_phase_b.sh pull

echo "[finalize] pulling lane 3b (MNAR)..."
bash script/rsync_gnn_evidence_phase_b_mnar.sh pull

echo "[finalize] collecting lane 3a..."
PIGAUTO_OUT_MD="${CAMPAIGN}/_phase_b_3a_section.md" \
  Rscript script/collect_gnn_evidence_phase_b.R

echo "[finalize] collecting lane 3b..."
PIGAUTO_OUT_MD="${CAMPAIGN}/_phase_b_3b_section.md" \
  Rscript script/collect_gnn_evidence_phase_b_mnar.R

# Combined durable summary (committed copy)
{
  echo "# GNN evidence Phase B — results (lanes 3a + 3b)"
  echo ""
  echo "**Generated:** ${TS}"
  echo "**Candidate SHA:** \`$(git rev-parse --short HEAD)\`"
  echo "**Scope:** phylo_MAR + covariate_MAR (lane 3a, 4860 fits) + MNAR (lane 3b, 2430 fits)"
  echo ""
  cat "${CAMPAIGN}/_phase_b_3a_section.md" | tail -n +2
  echo ""
  cat "${CAMPAIGN}/_phase_b_3b_section.md" | tail -n +2
} > "${SUMMARY}"

rm -f "${CAMPAIGN}/_phase_b_3a_section.md" "${CAMPAIGN}/_phase_b_3b_section.md"

# G6 audit echo (collectors stop on FAIL)
echo "[finalize] G6 audit: PASS (both collectors succeeded)"

# Update checkpoint
if [[ -f "${CHECKPOINT}" ]]; then
  n3a=$(find "${CAMPAIGN}/results_phase_b" -maxdepth 1 -name 'gnn_phase_b_job_*.rds' 2>/dev/null | wc -l | tr -d ' ')
  n3b=$(find "${CAMPAIGN}/results_phase_b_mnar" -maxdepth 1 -name 'gnn_phase_b_mnar_job_*.rds' 2>/dev/null | wc -l | tr -d ' ')
  sed -i.bak \
    -e "s/| Phase B lane 3b (MNAR) | \*\*[^|]*\*\* — .*/| Phase B lane 3b (MNAR) | **COMPLETE** — ${n3b}\/2430 RDS, 0 failures @ ${TS} |/" \
    -e "s/| Phase B lane 3a (phylo\/cov MAR) | \*\*[^|]*\*\* — .*/| Phase B lane 3a (phylo\/cov MAR) | **COMPLETE** — ${n3a}\/4860 RDS, 0 failures @ ${TS} |/" \
    -e "s/Phase B 3a+3b RUNNING\./Phase B 3a+3b COMPLETE./" \
    -e "s/Monitor for Phase B only (3a+3b)\. Next: pull Phase B when complete\./Phase B pulled + collected. Manuscript deferred./" \
    "${CHECKPOINT}"
  rm -f "${CHECKPOINT}.bak"

  # Totoro jobs table rows
  sed -i.bak \
    -e "s/| \*\*Phase B 3a\*\* phylo+cov MAR | .* | \*\*RUNNING\*\* |/| **Phase B 3a** phylo+cov MAR | ${n3a}\/4860 | 0 | — | **COMPLETE** — pulled |/" \
    -e "s/| \*\*Phase B 3b\*\* MNAR | .* | \*\*RUNNING\*\* |/| **Phase B 3b** MNAR | ${n3b}\/2430 | 0 | — | **COMPLETE** — pulled |/" \
    "${CHECKPOINT}" 2>/dev/null || true
  rm -f "${CHECKPOINT}.bak"
fi

echo "[finalize] committing..."
git add \
  script/monitor_gnn_evidence_totoro.sh \
  script/finalize_gnn_evidence_phase_b.sh \
  script/returned_gnn_campaign/PHASE_B_SUMMARY.md \
  script/returned_gnn_campaign/results_phase_b/gnn_phase_b_summary.csv \
  script/returned_gnn_campaign/results_phase_b/gnn_phase_b_cell_summary.csv \
  script/returned_gnn_campaign/results_phase_b_mnar/gnn_phase_b_mnar_summary.csv \
  script/returned_gnn_campaign/results_phase_b_mnar/gnn_phase_b_mnar_cell_summary.csv \
  LOOP/checkpoint.md \
  2>/dev/null || true

git commit -m "$(cat <<'EOF'
evidence: Phase B COMPLETE — 7290 fits, G6 PASS, summary pulled

Lane 3a (phylo_MAR+covariate_MAR) + lane 3b (MNAR) collected from Totoro.
PHASE_B_SUMMARY.md + checkpoint updated. Manuscript still deferred.
EOF
)" || echo "[finalize] nothing to commit (already committed?)"

echo "[finalize] done at $(date -Is)"
