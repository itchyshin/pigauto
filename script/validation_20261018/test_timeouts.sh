#!/usr/bin/env bash
set -euo pipefail
[[ "$(uname -s)" == Linux ]] || exit 77
command -v timeout >/dev/null || exit 77
work="$(mktemp -d)"
export VALIDATION_TIMEOUT_TEST_PID="$work/child.pid"
# The R subprocess inherits its parent's group when its own timeout is zero.
# GNU timeout owns and kills the complete outer group.
set +e
timeout --signal=KILL 2s Rscript -e 'system2("sh",c("-c",shQuote("echo $$ > \"$VALIDATION_TIMEOUT_TEST_PID\"; exec sleep 30")),timeout=0)' >/dev/null 2>&1
status=$?
set -e
[[ "$status" == 137 || "$status" == 124 ]] || exit 1
[[ -s "$VALIDATION_TIMEOUT_TEST_PID" ]] || exit 1
pid="$(cat "$VALIDATION_TIMEOUT_TEST_PID")"
state="$(ps -o stat= -p "$pid" || true)"
if [[ -n "$state" && "$state" != Z* ]]; then
  kill -KILL "$pid" 2>/dev/null || true
  echo "TIMEOUT_FAIL: child remained live" >&2
  exit 1
fi
rm -f "$VALIDATION_TIMEOUT_TEST_PID"
rmdir "$work"
echo TIMEOUT_GROUP_PASS
