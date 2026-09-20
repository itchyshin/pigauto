#!/usr/bin/env bash
# Build the results page from the template and the aggregated data.
#   bash script/campaign_sim_page_build.sh <agg-prefix> <out.html> "<status line>"
set -euo pipefail
AGG="${1:?agg prefix}"; OUT="${2:?out html}"; STATUS="${3:-Campaign complete.}"
HERE="$(cd "$(dirname "$0")" && pwd)"
TMP="$(mktemp -d)"
Rscript "$HERE/campaign_sim_page_data.R" --agg "$AGG" --out "$TMP/data.json" >&2
python3 - "$HERE/campaign_sim_page.template.html" "$TMP/data.json" "$OUT" "$STATUS" <<'PY'
import json, sys
tpl, data, out, status = sys.argv[1:5]
html = open(tpl).read()
html = html.replace("__DATA__", open(data).read().strip())
html = html.replace("__STATUS__", json.dumps(status))
open(out, "w").write(html)
print("wrote", out, len(html), "bytes")
PY
