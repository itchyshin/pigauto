#!/usr/bin/env bash
set -euo pipefail

runner="$(cd "$(dirname "$0")" && pwd)/totoro_production.sh"

set +e
output=$(env -u PIGAUTO_PRODUCTION_APPROVED bash "$runner" 2>&1)
status=$?
set -e
if [[ "$status" -ne 3 || "$output" != *"Production is gated"* ]]; then
  echo "FAIL: launcher did not refuse production without approval." >&2
  exit 1
fi

set +e
output=$(PIGAUTO_PRODUCTION_APPROVED=YES bash "$runner" 2>&1)
status=$?
set -e
if [[ "$status" -ne 2 || "$output" != *"Usage:"* ]]; then
  echo "FAIL: approved launcher did not stop safely when RUN_DIR was absent." >&2
  exit 1
fi

set +e
output=$(PIGAUTO_PRODUCTION_APPROVED=YES bash "$runner" \
  /home/snakagaw/hsq_work/should-not-create 2>&1)
status=$?
set -e
if [[ "$status" -ne 3 || "$output" != *"fresh live capacity"* ]]; then
  echo "FAIL: launcher did not require a fresh Totoro capacity check." >&2
  exit 1
fi

echo "PRODUCTION_LAUNCHER_GUARD_PASS"
