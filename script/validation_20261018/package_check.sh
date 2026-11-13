#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
: "${VALIDATION_CHECK_DIR:?Set an isolated check output directory outside frozen source}"
mkdir -p "$VALIDATION_CHECK_DIR"
cd "$VALIDATION_CHECK_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PIG_TORCH_THREADS=1
/usr/bin/time -v timeout --signal=TERM --kill-after=20s 1200s bash -c '
  set -euo pipefail
  R CMD build "$1"
  pkg=$(find . -maxdepth 1 -name "pigauto_*.tar.gz" -print)
  [[ -n "$pkg" ]] || exit 2
  R CMD check --no-manual "$pkg"
' validation-check "$ROOT"
