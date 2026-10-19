#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
receipt_dir="${1:-$repo_root/docs/dev-log/validation-20261018/receipts}"
socket="${HOME}/.ssh/cm-snakagaw@totoro.biology.ualberta.ca:22"
host="snakagaw@totoro.biology.ualberta.ca"
remote_root="/home/snakagaw/hsq_work/pigauto_validation_20261018/source-5e25a197"
remote_traits="/home/snakagaw/pigauto_sim/dlam/bien/bien_trait_means.rds"
remote_tree="/home/snakagaw/pigauto_sim/dlam/bien/bien_tree.rds"

[[ -S "$socket" ]] || { echo "TOTORO_EXISTING_SOCKET_MISSING" >&2; exit 2; }
receipt="$receipt_dir/production-bien-source-roster.csv"
[[ -f "$receipt" ]] || { echo "BIEN_SOURCE_ROSTER_RECEIPT_MISSING" >&2; exit 2; }

remote_output="$(ssh -S "$socket" -o ControlMaster=no -o BatchMode=yes "$host" \
  "cd '$remote_root' && VALIDATION_BIEN_TRAITS='$remote_traits' VALIDATION_BIEN_TREE='$remote_tree' Rscript -" \
  < "$repo_root/script/validation_20261018/derive_bien_roster.R")"
expected_output="$(cat "$receipt")"
if [[ "$remote_output" != "$expected_output" ]]; then
  echo "BIEN_TOTORO_ROSTER_MISMATCH" >&2
  diff -u <(printf '%s\n' "$expected_output") <(printf '%s\n' "$remote_output") || true
  exit 1
fi
printf '%s\n' "$remote_output"
echo "BIEN_TOTORO_ROSTER_VERIFIED"
