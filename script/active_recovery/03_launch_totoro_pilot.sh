#!/usr/bin/env bash
# Launch the pre-registered Stage-A 100-replicate-per-cell Totoro pilot.
# Run from the pigauto repository root after checking out the recorded SHA.
set -euo pipefail

readonly TOTORO_CORE_CAP=150
readonly NWORKERS="${NWORKERS:-96}"
readonly REPLICATES="${REPLICATES:-100}"
readonly EPOCHS="${EPOCHS:-100}"
readonly CAMPAIGN_DIR="${1:?usage: $0 CAMPAIGN_DIR}"

if (( NWORKERS < 1 || NWORKERS > TOTORO_CORE_CAP )); then
  echo "NWORKERS=${NWORKERS} violates the ${TOTORO_CORE_CAP}-core Totoro cap" >&2
  exit 2
fi
if (( REPLICATES != 100 )); then
  echo "This launcher is registered only for the 100-replicate pilot" >&2
  exit 2
fi
if [[ ! -f DESCRIPTION || ! -f script/active_recovery/01_run_active_recovery_cell.R ]]; then
  echo "run from the active-recovery pigauto checkout" >&2
  exit 2
fi
if [[ -e "$CAMPAIGN_DIR" ]]; then
  echo "campaign directory already exists: $CAMPAIGN_DIR" >&2
  exit 2
fi

mkdir -p "$CAMPAIGN_DIR/results"
git rev-parse HEAD > "$CAMPAIGN_DIR/source_sha.txt"
{
  printf 'campaign=stage-a-active-recovery-pilot\n'
  printf 'replicates_per_cell=%s\n' "$REPLICATES"
  printf 'n_values=100,300\n'
  printf 'lambda_values=1,0.2\n'
  printf 'families=continuous,binary\n'
  printf 'workers=%s\n' "$NWORKERS"
  printf 'openblas_threads=1\n'
  printf 'epochs=%s\n' "$EPOCHS"
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$CAMPAIGN_DIR/config.txt"
Rscript -e 'writeLines(capture.output(sessionInfo()), commandArgs(TRUE)[1])' \
  "$CAMPAIGN_DIR/sessionInfo.txt"

for n in 100 300; do
  for lambda in 1 0.2; do
    for family in continuous binary; do
      for replicate in $(seq 1 "$REPLICATES"); do
        printf '%s\t%s\t%s\t%s\n' "$n" "$lambda" "$family" "$replicate"
      done
    done
  done
done > "$CAMPAIGN_DIR/tasks.tsv"

if [[ $(wc -l < "$CAMPAIGN_DIR/tasks.tsv" | tr -d ' ') != 800 ]]; then
  echo "pilot task manifest must contain exactly 800 receipts" >&2
  exit 2
fi

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export CAMPAIGN_DIR EPOCHS

run_one() {
  local n="$1" lambda="$2" family="$3" replicate="$4"
  Rscript script/active_recovery/01_run_active_recovery_cell.R \
    "$n" "$lambda" "$family" "$replicate" "$CAMPAIGN_DIR/results" "$EPOCHS"
}
export -f run_one

printf 'launching 800 receipts with %s workers\n' "$NWORKERS" | tee "$CAMPAIGN_DIR/launcher.log"
< "$CAMPAIGN_DIR/tasks.tsv" xargs -P "$NWORKERS" -n 4 bash -c 'run_one "$@"' _ \
  >> "$CAMPAIGN_DIR/launcher.log" 2>&1
printf 'completed_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$CAMPAIGN_DIR/config.txt"
Rscript script/active_recovery/02_summarise_active_recovery.R \
  "$CAMPAIGN_DIR/results" "$CAMPAIGN_DIR/summary.rds" \
  > "$CAMPAIGN_DIR/summary.txt" 2>&1
