#!/usr/bin/env bash
# Build the campaign pages (report, summary, accuracy) from the continuous aggregator's JSON (__DATA__) and the
# discrete results' JSON (__DISC__, written here by script/rubin_disc_pages_data.py from script/rubin_study/data/).
#   bash script/rubin_report_build.sh [AGG_DIR] [OUT_HTML]
set -euo pipefail
AGG="${1:-$HOME/pigauto_rubin_pool/agg}"
OUT="${2:-docs/dev-log/arc/2026-09-25-rubin-campaign-report.html}"
DISC="script/rubin_study/data/discrete/pages.json"
python3 script/rubin_disc_pages_data.py script/rubin_study/data "$DISC"
python3 - "$AGG/report.json" script/rubin_report_template.html "$OUT" "$DISC" <<'PY'
import sys
data, tpl, out, disc = sys.argv[1:5]
t = open(tpl).read(); d = open(data).read()
assert "__DATA__" in t and "__DISC__" in t
t = t.replace("__DATA__", d).replace("__DISC__", open(disc).read())
assert "__DATA__" not in t and "__DISC__" not in t
open(out, "w").write(t)
print(out, len(t))
PY
# one-page summary from the same JSON
python3 - "$AGG/report.json" script/rubin_summary_template.html docs/dev-log/arc/2026-09-25-rubin-summary.html "$DISC" <<'PY'
import sys
data, tpl, out, disc = sys.argv[1:5]
t = open(tpl).read(); assert "__DATA__" in t and "__DISC__" in t
t = t.replace("__DATA__", open(data).read()).replace("__DISC__", open(disc).read())
assert "__DATA__" not in t and "__DISC__" not in t
open(out, "w").write(t); print(out)
PY
# accuracy companion from the same JSON (per-value and downstream accuracy tables, rounded to 5 places)
python3 - "$AGG/report.json" script/rubin_accuracy_template.html docs/dev-log/arc/2026-09-25-rubin-accuracy.html "$DISC" <<'PY'
import sys, json
data, tpl, out, disc = sys.argv[1:5]
r = json.load(open(data)); t = open(tpl).read(); assert "__DATA__" in t and "__DISC__" in t
rd = lambda rows: [{k: (round(v, 5) if isinstance(v, float) else v) for k, v in x.items()} for x in rows]
D = {"cl": rd(r["cells_l"]), "cn": rd(r["cells_n"]), "dl": rd(r["down_l"]), "dn": rd(r["down_n"])}
dd = json.load(open(disc)); DD = {"cl": dd["cl"], "pr": dd["pr"]}   # the accuracy page needs only these
t = t.replace("__DATA__", json.dumps(D, separators=(",", ":"))).replace("__DISC__", json.dumps(DD, separators=(",", ":")))
assert "__DATA__" not in t and "__DISC__" not in t
open(out, "w").write(t); print(out)
PY
