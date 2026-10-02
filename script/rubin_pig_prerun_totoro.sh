#!/usr/bin/env bash
# Totoro driver for the pig_post pre-run (pigauto posterior MI, PR #189) of the Rubin study.
# docs/dev-log/arc/2026-10-01-rubin-pigauto-prerun.md. Approved as part of the 2026-10-01 plan (36 fits, <= 40 cores,
# <= 1 h); the full campaign is NOT launched by this script (D-287: needs Shinichi's approval).
#
#   ssh totoro            # through the ~/.ssh/cm-* ControlMaster socket
#   cd ~/pigauto_rubin
#   PIGAUTO_SHA=<merged main sha> CONFIRM=yes bash script/rubin_pig_prerun_totoro.sh 36
#
# Grid: n {100, 300, 1000} x lambda {0.3, 0.7, 1} x rho {0, 0.5} x seeds {1, 2} = 36 cells, arm pig_post only,
# M = 20, MCAR 30%. Each cell is capped at 1 h (timeout 3600): a capped cell is a lower bound on fit time, and the
# pre-run cannot outlast the approved hour. Same datasets as the stored campaign (rubin_cell.R regenerates them from (n, lambda, rho, seed)).
# Output goes to its own directory (the cell tag does not name the arms). Resume: existing rds are skipped.
set -euo pipefail

PAR="${1:-36}"
ROOT="${RUBIN_ROOT:-$HOME/pigauto_rubin}"
ARMS="${ARMS:-pig_post}"   # pig_conf = the opt-in conformal negative control
OUT="$ROOT/prerun_${ARMS//,/_}"; [ "$ARMS" = pig_post ] && OUT="$ROOT/prerun_pig"; LOG="$ROOT/logs"; mkdir -p "$OUT" "$LOG"
[ "$PAR" -le 40 ] || { echo "PAR=$PAR exceeds the 40 cores approved for this pre-run"; exit 1; }

free -g
BUSY=$(ps -u "$USER" -o pcpu= | awk '{s+=$1} END {printf "%d", s/100 + 0.5}')
echo "cores already busy for $USER: $BUSY"
[ $(( BUSY + PAR )) -le 150 ] || { echo "PAR=$PAR + $BUSY busy exceeds the 150-core cap (D-143)"; exit 1; }

# pigauto comes from a private library ($ROOT/rlib) so other Totoro lanes keep their installed build.
export R_LIBS="$ROOT/rlib:${R_LIBS:-}"
# pig_conf fits pigauto's GNN, which needs the torch runtime: a private copy in $ROOT/torch_home (torch::install_torch()
# with TORCH_HOME set), so the shared ~/R/lib torch package is untouched. Checked before launch: the 2026-10-01 first
# pig_conf launch failed in every cell on a missing runtime.
export TORCH_HOME="$ROOT/torch_home"
case "$ARMS" in *pig_conf*) Rscript -e 'stopifnot("torch runtime missing (set TORCH_HOME, run torch::install_torch())" = torch::torch_is_installed()); cat("torch runtime OK\n")' ;; esac
# The installed pigauto must be the merged-main build, installed from GitHub so DESCRIPTION carries RemoteSha.
: "${PIGAUTO_SHA:?set PIGAUTO_SHA to the merged main sha}"
Rscript -e "d <- utils::packageDescription('pigauto'); s <- if (is.null(d\$RemoteSha)) '' else d\$RemoteSha;
  stopifnot('installed pigauto is not the expected sha' = identical(s, '$PIGAUTO_SHA'),
            'installed pigauto has no posterior draws' = 'posterior' %in% eval(formals(pigauto::multi_impute)\$draws_method));
  cat('pigauto', s, 'OK\n')"
for f in campaign_gnn_off_lib.R rubin_lib.R rubin_freq.R rubin_bace.R rubin_cell.R rubin_pigauto.R; do
  [ -s "$ROOT/script/$f" ] || { echo "missing $ROOT/script/$f (rsync the lane's script/ first)"; exit 1; }
done

[ "${CONFIRM:-no}" = yes ] || { echo "dry run: set CONFIRM=yes to launch"; exit 0; }

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

JOBS="$LOG/prerun_${ARMS//,/_}_jobs.txt"; : > "$JOBS"
for n in 1000 300 100; do for lambda in 0.3 0.7 1; do for rho in 0 0.5; do for seed in 1 2; do
  echo "--n $n --seed $seed --lambda $lambda --rho $rho --M 20 --arms $ARMS --out $OUT" >> "$JOBS"
done; done; done; done
echo "[$(date +%FT%T)] prerun_${ARMS} jobs=$(wc -l < "$JOBS") parallel=$PAR out=$OUT sha=$PIGAUTO_SHA" | tee -a "$LOG/prerun_${ARMS//,/_}.log"

setsid nohup bash -c "
  xargs -P $PAR -L 1 -a '$JOBS' bash -c 'timeout 3600 Rscript \"$ROOT/script/rubin_cell.R\" \"\$@\" >> \"$LOG/prerun_${ARMS//,/_}_cells.log\" 2>&1' _
  echo \"[\$(date +%FT%T)] prerun done\" >> '$LOG/prerun_${ARMS//,/_}.log'
" > /dev/null 2>&1 &
echo "launched; pgid $(ps -o pgid= $! | tr -d ' '); stop with: kill -- -<pgid>"
