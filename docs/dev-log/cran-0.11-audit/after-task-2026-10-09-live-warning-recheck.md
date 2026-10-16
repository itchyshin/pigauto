# After-task: live homepage warning recheck

## 1. Goal

Confirm the current Pages deployment and determine whether the homepage warning-rendering defect remains in the live reader surface.

## 2. Implemented

Added a current deployment and browser observation to the CRAN audit ledger. No package source, website source, deployment, or exact release artifact was changed.

## 3a. Decisions and Rejected Alternatives

The homepage display defect remains separate from G7’s bounded route-verification pass. The latest deployment matches merged main, but the warning still renders as a literal marker. Assumption: the package’s current Pages deployment is represented by GitHub’s latest successful `github-pages` deployment; if another deployment completes, repeat the browser check against it.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-live-warning-recheck.md`

## 5. Checks Run

- Chrome GitHub Deployments page: latest Pages success is pkgdown run #674 for merge commit `d76804e768bf77f430f43dcf591633f9cd900dab`.
- Chrome cache-busted homepage: version 0.11.0 is served and the page displays literal `[!WARNING]` text.
- `git show origin/main:README.md`: merged README still uses `[!WARNING]` callout syntax.
- Existing local fix evidence remains the exact-head pkgdown build at `45df106ba58f5bc97d7a274bc0a22a19d202fa17`, with the corrected `<blockquote>` render and build-log SHA-256 `ba7e2f361fa0c6980bbbbcca9fb342388602968edeee60963450e540f3c49db6`.
- This check did not rerun the retired-route, sitemap, or search-index audit; the existing G7 receipt is bounded to the same deployment commit.
- `slop_check.py`: 0 findings. `check-after-task.R`: required report structure passed; the overall Unlazy audit remains open on exact-artifact and independent-review gates, and the repo-wide checker also reports separate imputation-sim leaves.

## 6. Tests of the Tests

No tests changed. This is a live reader-surface observation cross-checked against merged source and a separately tested local rendering fix.

## 7a. Issue Ledger

- Confirmed current defect: the live homepage exposes the Markdown warning marker.
- Open: source fix is not merged or deployed.
- Open release gates remain: exact frozen-artifact Windows results and independent approval for that exact artifact.

## 8. Consistency Audit

The deployed content identifies version 0.11.0, and the latest successful Pages deployment corresponds to the merged commit. The source still contains the unsupported callout, matching the observed rendering. The local candidate fix has a separate successful render receipt.

## 9. What Did Not Go Smoothly

Nothing blocked this read-only recheck. GitHub’s deployment page was needed to bind the live page to a source commit.

## 10. Known Residuals

This turn did not re-audit all live routes, search results, or sitemap entries. It does not establish a new deployment, exact-tarball validation, Windows support, or CRAN readiness. The README source PR still awaits selection of its proposed title and body.

## 11. Team Learning

A route audit can pass while a separate visible rendering defect remains. Keep structural route/discovery checks and reader-facing rendering checks as distinct evidence.

## 12. Cross-Product Coverage

This covers the current homepage, merged README markup, and latest Pages deployment identity. It does NOT cover the full route set, retired URL behavior, sitemap/search completeness, Windows, the exact artifact, or CRAN acceptance.
