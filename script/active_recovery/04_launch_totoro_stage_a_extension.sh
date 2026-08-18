#!/usr/bin/env bash
# Extend the audited Stage-A pilot from 100 to 1,000 total replicates per cell.
# Run from the pigauto repository root at the recorded extension SHA.
set -euo pipefail

readonly TOTORO_CORE_CAP=150
readonly NWORKERS="${NWORKERS:-96}"
readonly REPLICATE_START=101
readonly REPLICATE_END=1000
readonly EXPECTED_TASKS=7200
readonly EXTENSION_DIR="${1:?usage: $0 EXTENSION_DIR PILOT_DIR}"
readonly PILOT_DIR="${2:?usage: $0 EXTENSION_DIR PILOT_DIR}"
readonly PILOT_SHA="beee5df360759f53cd22b595366436bca50845a5"

if (( NWORKERS < 1 || NWORKERS > TOTORO_CORE_CAP )); then
  echo "NWORKERS=${NWORKERS} violates the ${TOTORO_CORE_CAP}-core Totoro cap" >&2
  exit 2
fi
if [[ ! -f DESCRIPTION || ! -f script/active_recovery/01_run_active_recovery_cell.R ]]; then
  echo "run from the active-recovery pigauto checkout" >&2
  exit 2
fi
if [[ -e "$EXTENSION_DIR" || ! -d "$PILOT_DIR/results" ]]; then
  echo "extension directory must be new and pilot results must exist" >&2
  exit 2
fi
if [[ $(find "$PILOT_DIR/results" -type f -name '*.rds' | wc -l | tr -d ' ') != 800 ]]; then
  echo "expected exactly 800 completed pilot receipts" >&2
  exit 2
fi
if ! Rscript -e 'fs <- list.files(commandArgs(TRUE)[1L], pattern = "\\.rds$", full.names = TRUE); x <- do.call(rbind, lapply(fs, readRDS)); stopifnot(all(x$source_sha == "beee5df360759f53cd22b595366436bca50845a5"))' \
  "$PILOT_DIR/results"; then
  echo "pilot source SHA receipt check failed" >&2
  exit 2
fi

mkdir -p "$EXTENSION_DIR/results"
git rev-parse HEAD > "$EXTENSION_DIR/source_sha.txt"
sha256sum script/active_recovery/00_prepare_active_recovery.R \
  script/active_recovery/01_run_active_recovery_cell.R \
  script/active_recovery/02_summarise_active_recovery.R > "$EXTENSION_DIR/runner_sha256.txt"
{
  printf 'campaign=stage-a-active-recovery-extension\n'
  printf 'pilot_source_sha=%s\n' "$PILOT_SHA"
  printf 'replicate_range=%s-%s\n' "$REPLICATE_START" "$REPLICATE_END"
  printf 'n_values=100,300\n'
  printf 'lambda_values=1,0.2\n'
  printf 'families=continuous,binary\n'
  printf 'workers=%s\n' "$NWORKERS"
  printf 'openblas_threads=1\n'
  printf 'epochs=100\n'
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$EXTENSION_DIR/config.txt"
Rscript -e 'writeLines(capture.output(sessionInfo()), commandArgs(TRUE)[1])' \
  "$EXTENSION_DIR/sessionInfo.txt"

for n in 100 300; do
  for lambda in 1 0.2; do
    for family in continuous binary; do
      for replicate in $(seq "$REPLICATE_START" "$REPLICATE_END"); do
        printf '%s\t%s\t%s\t%s\n' "$n" "$lambda" "$family" "$replicate"
      done
    done
  done
done > "$EXTENSION_DIR/tasks.tsv"
if [[ $(wc -l < "$EXTENSION_DIR/tasks.tsv" | tr -d ' ') != "$EXPECTED_TASKS" ]]; then
  echo "extension task manifest must contain ${EXPECTED_TASKS} receipts" >&2
  exit 2
fi

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export EXTENSION_DIR

run_one() {
  local n="$1" lambda="$2" family="$3" replicate="$4"
  Rscript script/active_recovery/01_run_active_recovery_cell.R \
    "$n" "$lambda" "$family" "$replicate" "$EXTENSION_DIR/results" 100
}
export -f run_one

printf 'launching 7200 extension receipts with %s workers\n' "$NWORKERS" | tee "$EXTENSION_DIR/launcher.log"
< "$EXTENSION_DIR/tasks.tsv" xargs -P "$NWORKERS" -n 4 bash -c 'run_one "$@"' _ \
  >> "$EXTENSION_DIR/launcher.log" 2>&1
printf 'completed_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$EXTENSION_DIR/config.txt"
