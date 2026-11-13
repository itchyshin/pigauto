#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: run_tasks.sh RUN_DIR [WORKERS]" >&2
  exit 2
fi
run_dir=$1
workers=${2:-2}
if [[ ! -f "$run_dir/run_config.rds" || ! -f "$run_dir/tasks.csv" ]]; then
  echo "Run directory must be initialized before dispatch." >&2
  exit 2
fi
if [[ ! "$workers" =~ ^[1-9][0-9]*$ || "$workers" -gt 48 ]]; then
  echo "WORKERS must be an integer from 1 to 48." >&2
  exit 2
fi
mkdir -p "$run_dir/logs"
lock_dir="$run_dir/logs/driver.lock"
if ! mkdir "$lock_dir" 2>/dev/null; then
  echo "Another launcher owns this run directory (or a stale lock needs inspection): $lock_dir" >&2
  exit 2
fi
trap 'rmdir "$lock_dir" 2>/dev/null || true' EXIT
stop_owned_workers() {
  if [[ -f "$run_dir/logs/xargs.pid" ]]; then
    kill -TERM "$(cat "$run_dir/logs/xargs.pid")" 2>/dev/null || true
  fi
  for pid_file in "$run_dir"/logs/task_process_*.pid; do
    [[ -f "$pid_file" ]] || continue
    task_pid=$(cat "$pid_file")
    kill -TERM -- "-$task_pid" 2>/dev/null || kill -TERM "$task_pid" 2>/dev/null || true
    rm -f "$pid_file"
  done
  rm -f "$run_dir/logs/xargs.pid"
  exit "$1"
}
trap 'stop_owned_workers 130' INT
trap 'stop_owned_workers 143' TERM
printf '%s\n' "$$" > "$run_dir/logs/driver.pid"
rm -f "$run_dir/logs/failures.tsv"
n_tasks=$(awk 'END { print NR - 1 }' "$run_dir/tasks.csv")
if [[ "$n_tasks" -lt 1 ]]; then
  echo "No campaign tasks were found." >&2
  exit 2
fi

set +e
seq 1 "$n_tasks" | xargs -P "$workers" -I '{}' bash -c '
  idx=$1
  run_dir=$2
  printf "%s\n" "$$" > "$run_dir/logs/worker_${idx}.pid"
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
  setsid /usr/bin/time -v Rscript --vanilla script/gnn_incremental/campaign.R \
    task "$run_dir" "$idx" > "$run_dir/logs/task_${idx}.log" 2>&1 &
  task_pid=$!
  printf "%s\n" "$task_pid" > "$run_dir/logs/task_process_${idx}.pid"
  wait "$task_pid"
  rc=$?
  rm -f "$run_dir/logs/task_process_${idx}.pid"
  if [[ "$rc" -ne 0 ]]; then
    printf "%s\t%s\n" "$idx" "$rc" >> "$run_dir/logs/failures.tsv"
  fi
  exit "$rc"
' _ '{}' "$run_dir" &
xargs_pid=$!
printf '%s\n' "$xargs_pid" > "$run_dir/logs/xargs.pid"
wait "$xargs_pid"
task_status=$?
rm -f "$run_dir/logs/xargs.pid"
set -e

Rscript --vanilla script/gnn_incremental/audit_run.R --run-dir "$run_dir" \
  > "$run_dir/logs/audit.log" 2>&1
Rscript --vanilla script/gnn_incremental/summarise.R --verify --run-dir "$run_dir" \
  > "$run_dir/logs/summary.log" 2>&1
if [[ "$task_status" -ne 0 ]]; then
  echo "One or more task processes failed; see logs/failures.tsv and per-task logs." >&2
  exit "$task_status"
fi
echo "CAMPAIGN_TASKS_AND_SUMMARY_PASS"
