# After-task: candidate-site retirement and warning verifier

## 1. Goal

Harden the local post-cleanup site check for the unmerged homepage warning correction after independent review found false-pass cases. Keep local candidate evidence separate from deployed-site and artifact gates.

## 2. Implemented

Hardened the standard-library Python verifier and expanded its suite from five to sixteen tests. After Pat's independent review found two remaining false-pass cases, the verifier now catches decimal zero opacity values and expands comma-separated alternatives inside `:is()` and `:where()` selectors. It also pins the retired-route receipt hash; checks exact sitemap, HTML, and search counts; requires discovery URLs to use the configured origin and site prefix; verifies the expected warning body in a blockquote with a bold label; scans supported HTML/CSS hiding mechanisms, attribute selectors, and same-site CSS imports; ignores script/style text; repeatedly decodes paths; and rejects traversal. Search path count and sitemap membership are checked, while route multiplicities are not pinned. The cleaned candidate still passes. The durable release ledger records that this is local candidate evidence only, with no change to the deployed-site verdict or release readiness.

## 3a. Decisions and Rejected Alternatives

The check examines generated files, sitemap targets, search-index targets, article-index links, and the homepage warning callout. It runs after `pkgdown/clean-internal-pages.R`, because raw pkgdown output can contain routes that the production cleanup intentionally removes. It does not treat a local build as evidence of deployment. The 44-route count is pinned to the current reviewed receipt so a shortened manifest cannot silently weaken this release check.

The static parser verifies blockquote/strong-label structure, required body text, hidden attributes/classes, and supported hiding declarations in matching inline or same-site CSS rules and imports. It does not implement a browser's complete CSS cascade or evaluate accessibility and rendered layout. These remain outside this local structural check.

## 4. Files Touched

- `script/cran-0.11-site/verify_candidate_site.py`
- `script/cran-0.11-site/test_candidate_site.py`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/candidate-site-warning-callout-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/provenance/candidate-site-verifier-tests-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-candidate-site-verifier.md`

## 5. Checks Run

- `python3 script/cran-0.11-site/test_candidate_site.py`: 16 tests passed in 0.929 seconds, with 35 negative-control scenarios.
- `python3 script/cran-0.11-site/verify_candidate_site.py --site-dir /private/tmp/pigauto-warning-callout-site --manifest docs/dev-log/cran-0.11-audit/provenance/live-retired-routes-2026-10-09-v2.tsv`: exit 0; 44 unique retired routes absent, 62 sitemap entries matching all 62 HTML files, 613 search entries, 568 in-sitemap paths, `CANDIDATE_SITE_OK`.
- `git diff --check`: passed after the verifier changes.
- After commit `1872bd7` was pushed to `release/cran-0.11-gate`, Chrome verified PR #228 remained Draft at 50 commits with no deployment. Workflow run `38001692992` was skipped under the configured pull-request pkgdown guard. It produced no site build or deployment result.
- Follow-up commit `74d5a9a` was pushed to the same branch. Chrome verified it as the current PR head; PR #228 remains Draft and unmerged, with all checks passing and one skipped check. The PR page showed no deployments.
- The evidence was later refreshed through commit `8fa149a`; Chrome showed PR #228 at 52 commits, Draft and unmerged. Workflow run `38003919177` was skipped under the pull-request pkgdown guard, and the PR page listed no deployment. This confirms an evidence push only; it does not establish a new site build or deployment.
- Candidate build log: 46,943 bytes, SHA-256 `b776a7a5a88dd6c61b07295dc5ffb5d2cca6a5960758f2639c1ab268b4a3188d`. Candidate source is `45df106ba58f5bc97d7a274bc0a22a19d202fa17`, parent `d76804e768bf77f430f43dcf591633f9cd900dab`.
- Captured verifier output: SHA-256 `444b18fca2d969928bea89643c803009d63379226d3de97e3240689f6dbf75f4`.

## 6. Tests of the Tests

Thirty-five negative-control scenarios passed across retired output and discovery targets, old warning syntax, unrelated or incomplete warning text, hidden attributes/classes/stylesheets, decimal zero opacity, `:is()` and `:where()` selector alternatives, script-only warning text, empty/truncated inventories, external discovery URLs, truncated/duplicate/substituted route receipts, encoded retired URLs, and traversal routes. The candidate output is a positive control for the pinned 44-route receipt, page/search counts, and expected warning text.

## 7a. Issue Ledger

- Fixed for the candidate build: static output now gives a repeatable assertion for the homepage callout and for route cleanup across generated files and discovery indexes.
- Still open: the corrected README source is unmerged; the live homepage is unchanged.
- Still open: the exact frozen tarball's Windows checks fail and are not bound to its SHA-256. G8 and G9 remain unmet.

## 8. Consistency Audit

Compared the local candidate build after the same cleanup used by the production site workflow. The verifier checks all 44 unique receipt entries against the generated tree, sitemap, search index, and article index; exact page/search inventory counts; discovery origin/prefix; and homepage warning structure/body against common hiding mechanisms. The website candidate is separate from the frozen tarball and does not alter package code, defaults, or multiple-imputation behavior.

## 9. What Did Not Go Smoothly

The full GitHub Pages wrapper remained in mixed-types vignette rendering beyond its four-minute estimate and was interrupted. The direct `pkgdown::build_site()` call and subsequent cleanup had completed. Chrome rejected the local `file:` preview URL, and the browser policy prohibited substituting a local server, so visual review of this candidate was not performed. Unlazy's runner passed both candidate checks but could not persist its updated ignored ledger because a lock-file write returned `EPERM`; its measured outcomes were recorded manually in that ledger.

## 10. Known Residuals

The current live homepage still displays `[!WARNING]`. The candidate output is local and unmerged; it has no deployment receipt. Pat's final review passes the documented bounded static-check claim and confirms the 16 tests and candidate audit. A general visibility checker could still miss CSS expressions that compute opacity to zero, such as `opacity: calc(0)` or a custom property resolving to zero. The verifier does not calculate the browser's full CSS cascade. The wrapper build, computed browser appearance, accessibility review, and local visual review remain unverified. This work does not close the exact-artifact or independent-release-review gates, and nothing was submitted to CRAN.

## 11. Team Learning

Run route-retirement checks against final output after production cleanup because a raw pkgdown build can include routes that the deployment process removes.

Memory receipt: loaded the pigauto LOAD-FIRST manifest with `route.py pigauto` and searched the brain for current pigauto release decisions. No brain files were changed.

Golden Set: not run; no package behavior or statistical result changed.

## 12. Cross-Product Coverage

This check covers local static HTML output, sitemap targets, search targets, article-index links, and the homepage callout after cleanup. It does NOT cover visual appearance or accessibility, an actual Pages deployment, CRAN tarball contents, Windows checks, or CRAN acceptance.
