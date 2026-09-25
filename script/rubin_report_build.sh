#!/usr/bin/env bash
# Build the campaign report page from the aggregator's JSON.
#   bash script/rubin_report_build.sh [AGG_DIR] [OUT_HTML]
set -euo pipefail
AGG="${1:-$HOME/pigauto_rubin_pool/agg}"
OUT="${2:-docs/dev-log/arc/2026-09-25-rubin-campaign-report.html}"
python3 - "$AGG/report.json" script/rubin_report_template.html "$OUT" <<'PY'
import sys
data, tpl, out = sys.argv[1:4]
t = open(tpl).read(); d = open(data).read()
assert "__DATA__" in t
open(out, "w").write(t.replace("__DATA__", d))
print(out, len(t) + len(d))
PY
