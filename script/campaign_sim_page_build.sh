#!/usr/bin/env bash
# Build the results page from the template and the aggregated data.
#   bash script/campaign_sim_page_build.sh <agg-prefix> <out.html> "<status line>" [<avonet-agg-prefix>]
# The fourth argument is optional: give it the AVONET case-study aggregate prefix and the page gains
# its "AVONET 300" tab. Omit it and that tab reports that the case study has not landed yet.
set -euo pipefail
AGG="${1:?agg prefix}"; OUT="${2:?out html}"; STATUS="${3:-Campaign complete.}"; AVO="${4:-}"
HERE="$(cd "$(dirname "$0")" && pwd)"
TMP="$(mktemp -d)"
Rscript "$HERE/campaign_sim_page_data.R" --agg "$AGG" --out "$TMP/data.json" >&2
if [ -n "$AVO" ] && [ -f "${AVO}_summary.csv" ]; then
  Rscript "$HERE/campaign_sim_page_data.R" --agg "$AVO" --out "$TMP/avonet.json" >&2
else
  echo 'null' > "$TMP/avonet.json"
fi
python3 - "$HERE/campaign_sim_page.template.html" "$TMP/data.json" "$OUT" "$STATUS" "$TMP/avonet.json" <<'PY'
import json, sys
tpl, data, out, status, avo = sys.argv[1:6]
html = open(tpl).read()
html = html.replace("__DATA__", open(data).read().strip())
html = html.replace("__AVONET__", open(avo).read().strip())
html = html.replace("__STATUS__", json.dumps(status))
open(out, "w").write(html)
print("wrote", out, len(html), "bytes")
PY
