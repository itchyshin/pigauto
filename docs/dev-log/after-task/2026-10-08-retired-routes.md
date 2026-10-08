# After-task report: live retired-route check

## 1. Goal

Refresh the deployed-site evidence for the retired pigauto pages and record the result in the CRAN 0.11 evidence branch.

## 2. Implemented

Added a route-by-route receipt for the 34 retired `/dev/` pages and four legacy coordination pages. The Codex In-app Browser displayed pigauto's `Page not found (404)` page for all 38 routes.

## 3a. Decisions and Rejected Alternatives

The receipt records the visible page title. This browser pass did not measure HTTP response codes. I left search and sitemap claims open because this pass did not verify their current state.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/provenance/live-retired-routes-2026-10-08.tsv`
- `docs/dev-log/cran-0.11-audit/retirement-manifest.md`
- `docs/dev-log/after-task/2026-10-08-retired-routes.md`
- `docs/dev-log/after-task/2026-10-08-retired-routes-assessment.json`

## 5. Checks Run

- Opened every route in the Codex In-app Browser with `?audit=20261008-live-routes`: 38 of 38 showed the visible title `Page not found (404) • pigauto`.
- Receipt structure check: 38 unique routes, 34 under `/dev/`, four legacy routes, six tab-separated fields per row, and all rows marked `PASS_RENDERED_NOT_FOUND`.
- `git diff --check` passed for the receipt.
- Search probe: the deployed search box accepted text, but exposed no verifiable result list. Search remains open.
- Sitemap: not rechecked in this slice.

## 6. Tests of the Tests

The receipt validator checks the route count and required fields. That catches omitted or malformed rows. It cannot establish the live result; the browser observations provide that evidence. No planted route control was run.

## 7a. Issue Ledger

- Closed for this slice: live rendered-page check of the 34 manifested retired routes and four named legacy routes.
- Open: current search-index and sitemap checks; direct HTTP status for the routes; local visual review of the candidate site.

## 8. Consistency Audit

Compared the 34 `/dev/` paths with `retirement-manifest.md`, then opened those exact paths plus `/VALIDATION_LEDGER.html`, `/AGENTS.html`, `/CLAUDE.html`, and `/goodagents.html`. All displayed the same not-found page. This receipt says nothing about current sitemap or search-index contents.

## 9. What Did Not Go Smoothly

The search input did not expose a result list in its accessibility state after text entry. I kept that check unresolved. A file-patch encoding mistake was caught by the row validator and corrected before staging.

## 10. Known Residuals

This is deployed route evidence only. It does not close the local site visual review, source and documentation fixes, BirdTree redistribution documentation, frozen-tarball checks, platform results, or independent review. The release remains `NOT_READY`.

## 11. Team Learning

Record the observed page title separately from the HTTP response status. A custom not-found page is useful route evidence, but it does not by itself prove the server's status code.

Memory receipt: loaded pigauto's `route.py` manifest and the project LOAD-FIRST conditions. Tree provenance and shared-lane ownership shaped the scope. The Golden Set was not checked. This slice covered live route behavior.

## 12. Cross-Product Coverage

Covered: live direct-page rendering for 34 retired benchmark paths and four legacy paths.

Does NOT cover: HTTP status codes, sitemap entries, search-index entries, retained-page content, local visual layout, package source or exact-artifact checks, Windows results, rights clearance, or release readiness.
