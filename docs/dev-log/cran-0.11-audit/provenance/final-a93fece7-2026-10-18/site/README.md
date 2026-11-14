# Rose post-merge Pages and live-site verification

**Verdict: PASS** for the merged source commit `b2e6054685e6bdedb073ba46e058af2c0a619a18` (PR #238), limited to the GitHub Pages deployment and live documentation site. This does not establish CRAN readiness or acceptance.

GitHub Pages workflow run `39090594794` completed successfully at the exact merged SHA. Its pkgdown job `117566315006` uploaded Pages artifact `12122696328`; the deployment log recorded `pages_build_version` and deployment ID `b2e6054685e6bdedb073ba46e058af2c0a619a18` and reported success.

The existing `verify_live_site.py` verifier exited 0 against `https://itchyshin.github.io/pigauto/`. It checked 44 retired routes and 62 retained routes with no failures; sitemap returned 200 with 62 entries and no retired targets; search returned 200 with 614 entries and no retired targets. The expected utility 404 page returned 200. One search result (`news/index.html`) contains historical `/dev/` text, but no retired target is indexed as a route.

Reviewed live pages in Chrome: `reference/impute.html`, `reference/pigauto_report.html`, and adjacent `articles/getting-started.html`. The references state baseline-only behavior by default (`gnn=FALSE`) and explicit GNN opt-in (`gnn=TRUE`); the getting-started guidance agrees. Rendered text was legible without overlap; long reference headings wrap beside the logo. Screenshots were inspected but not exported. No local HTTP server was started.

Evidence files:

- `verdict.json` — machine-readable verdict, exact run/deployment binding, scope and limits.
- `live-check/discovery-check.json` — route/search/sitemap summary.
- `live-check/retired-route-statuses.tsv` and `live-check/retained-route-statuses.tsv` — per-route receipts.

The generic main CI run `39090594763` was outside this review's ownership. No source or site files were edited, no package build or CRAN action was performed. Highest proven rung here is merged-main Pages deployment plus live-site content and retirement checks; CRAN artifact/submission/acceptance remain unproven.
