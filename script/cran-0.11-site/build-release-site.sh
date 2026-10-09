#!/usr/bin/env bash
set -Eeuo pipefail

repo_root=$(git rev-parse --show-toplevel)
source_commit=$(git -C "$repo_root" rev-parse HEAD)
site_root=$(mktemp -d "/private/tmp/pigauto-cran011-site-${source_commit:0:7}.XXXXXX")
cache_root=$(mktemp -d /private/tmp/pigauto-cran011-cache.XXXXXX)
cache_source=${PIGAUTO_SITE_CACHE_SOURCE:-/private/tmp/pigauto-cran011-site-final-01a11136/_site}
build_log="${site_root}.build.log"
crawl_log="${site_root}.crawl.log"

if [[ ! -d "$cache_source" ]]; then
  printf 'Verified pkgdown asset cache source is missing: %s\n' "$cache_source" >&2
  exit 2
fi

git -C "$repo_root" archive "$source_commit" | tar -x -C "$site_root"
printf '%s\n' "$source_commit" > "$site_root/SOURCE_COMMIT.txt"
python3 "$site_root/script/cran-0.11-site/seed_pkgdown_cache.py" \
  --source-site "$cache_source" --cache-root "$cache_root"

if ! (cd "$site_root" && \
  R_USER_CACHE_DIR="$cache_root" \
  Rscript --vanilla script/cran-0.11-site/build-offline.R > "$build_log" 2>&1); then
  tail -n 100 "$build_log" >&2
  exit 1
fi

PIGAUTO_SITE_DIR="$site_root/_site" \
  Rscript --vanilla "$site_root/pkgdown/clean-internal-pages.R"
if ! (cd "$site_root" && \
  python3 script/cran-0.11-site/crawl.py _site > "$crawl_log"); then
  cat "$crawl_log" >&2
  exit 1
fi
grep -F SITE_CRAWL_OK "$crawl_log"

(cd "$site_root" && Rscript --vanilla -e 'pkgdown::check_pkgdown()')
python3 - "$site_root/_site" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
retired = (
    "articles/simulation-study.html",
    "articles/simulation-study.md",
    "articles/articles/simulation-study.html",
    "articles/articles/simulation-study.md",
)
assert all(not (root / route).exists() for route in retired)
for name in ("search.json", "sitemap.xml"):
    text = (root / name).read_text()
    assert not any(route in text for route in retired), f"retired route found in {name}"
search_entries = json.loads((root / "search.json").read_text())
print(f"SITE_BUILD_AND_RETIREMENT_OK: source={Path(sys.argv[1]).parent.joinpath('SOURCE_COMMIT.txt').read_text().strip()} html={len(list(root.rglob('*.html')))} search={len(search_entries)}")
PY

grep -F 'Finished building pkgdown site for package pigauto' "$build_log"
tail -n 12 "$build_log"
printf 'SITE_ROOT=%s\nBUILD_LOG=%s\nCRAWL_LOG=%s\n' "$site_root" "$build_log" "$crawl_log"
