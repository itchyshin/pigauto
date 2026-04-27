#!/bin/bash
# script/run_phase23.sh
#
# Phase 2: bench_transformer_ablation_nosf smoke (architectures with the
#   safety floor OFF -- resolves transformer vs legacy GNN vs plain GCN
#   that the safety_floor=TRUE smoke could not).
# Phase 3: bench_ou_regime smoke (non-BM DGPs).
#
# Both at PIGAUTO_TIER=smoke and PIGAUTO_BENCH_EPOCHS=200 to fit the
# ~1pm deadline.  Each ~2 hours budget.
#
# Launched by a wakeup AFTER run_all_per_type_benches.sh finishes, so we
# don't oversubscribe the machine.

set -e
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
RSCRIPT=/usr/local/bin/Rscript
LOG_DIR="script"

# -----------------------------------------------------------------------
echo "[$(date +%H:%M:%S)] === Starting bench_transformer_ablation_nosf (smoke, EPOCHS=200) ==="
PIGAUTO_TIER=smoke PIGAUTO_BENCH_EPOCHS=200 \
  "$RSCRIPT" "${LOG_DIR}/bench_transformer_ablation_nosf.R" \
  > "${LOG_DIR}/bench_transformer_ablation_nosf_smoke.log" 2>&1
rc=$?
echo "[$(date +%H:%M:%S)] === Finished bench_transformer_ablation_nosf (rc=${rc}) ==="

# -----------------------------------------------------------------------
echo "[$(date +%H:%M:%S)] === Starting bench_ou_regime (smoke, EPOCHS=200) ==="
PIGAUTO_TIER=smoke PIGAUTO_BENCH_EPOCHS=200 \
  "$RSCRIPT" "${LOG_DIR}/bench_ou_regime.R" \
  > "${LOG_DIR}/bench_ou_regime_smoke.log" 2>&1
rc=$?
echo "[$(date +%H:%M:%S)] === Finished bench_ou_regime (rc=${rc}) ==="

echo "[$(date +%H:%M:%S)] === Phase 2+3 DONE ==="
