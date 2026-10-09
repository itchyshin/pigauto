# After-task: candidate-site retirement and warning verifier

## 1. Goal

Make the local post-cleanup site check repeatable for the unmerged homepage warning correction and confirm it does not restore retired pages or discovery entries.

## 2. Implemented

Added a standard-library Python verifier and five offline tests. Ran them against the candidate pkgdown output after the production cleanup script. Recorded the result and updated the CRAN audit ledger without changing the deployed-site verdict or release readiness.

## 3a. Decisions and Rejected Alternatives

The check examines generated files, sitemap targets, search-index targets, article-index links, and the homepage warning label. It runs after `pkgdown/clean-internal-pages.R`, because the raw pkgdown output can contain routes that the production cleanup intentionally removes. It does not treat a local build as evidence of deployment.

Assumption: the rendered visible label `Warning:` is an adequate static assertion for this correction. A visual and accessibility review after deployment could require a stronger criterion.

## 4. Files Touched

- `script/cran-0.11-site/verify_candidate_site.py`
- `script/cran-0.11-site/test_candidate_site.py`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/candidate-site-warning-callout-2026-10-09.log`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-candidate-site-verifier.md`

## 5. Checks Run

- `python3 script/cran-0.11-site/test_candidate_site.py`: 5 tests passed in 0.045 seconds.
- `python3 script/cran-0.11-site/verify_candidate_site.py --site-dir /private/tmp/pigauto-warning-callout-site --manifest docs/dev-log/cran-0.11-audit/provenance/live-retired-routes-2026-10-09-v2.tsv`: exit 0; 44 retired routes absent, 62 sitemap entries, 613 search entries, 568 entries with paths, 62 HTML pages, `CANDIDATE_SITE_OK`.
- `git diff --check`: passed before adding this report.
- Candidate build log: 46,943 bytes, SHA-256 `b776a7a5a88dd6c61b07295dc5ffb5d2cca6a5960758f2639c1ab268b4a3188d`. Candidate source is `45df106ba58f5bc97d7a274bc0a22a19d202fa17`, parent `d76804e768bf77f430f43dcf591633f9cd900dab`.
- Captured verifier output: SHA-256 `444b18fca2d969928bea89643c803009d63379226d3de97e3240689f6dbf75f4`.

## 6. Tests of the Tests

Five negative controls passed: the tests reject a retired output file, sitemap target, search target, article-index link, and the literal `[!WARNING]` marker. These controls show each principal retirement surface can fail independently.

## 7a. Issue Ledger

- Fixed for the candidate build: static output now gives a repeatable assertion for the homepage callout and for route cleanup across generated files and discovery indexes.
- Still open: the corrected README source is unmerged; the live homepage is unchanged.
- Still open: the exact frozen tarball's Windows checks fail and are not bound to its SHA-256. G8 and G9 remain unmet.

## 8. Consistency Audit

Compared the local candidate build after the same cleanup used by the production site workflow. The verifier reads the deployed retired-route receipt and checks all 44 entries against the generated tree, sitemap, search index, and article index. It also checks the homepage output. The website candidate is separate from the frozen tarball and does not alter package code, defaults, or multiple-imputation behavior.

## 9. What Did Not Go Smoothly

The full GitHub Pages wrapper remained in mixed-types vignette rendering beyond its four-minute estimate and was interrupted. The direct `pkgdown::build_site()` call and subsequent cleanup had completed. Chrome rejected the local `file:` preview URL, and the browser policy prohibited substituting a local server, so visual review of this candidate was not performed. Unlazy's runner passed both candidate checks but could not persist its updated ignored ledger because a lock-file write returned `EPERM`; its measured outcomes were recorded manually in that ledger.

## 10. Known Residuals

The current live homepage still displays `[!WARNING]`. The candidate output is local and unmerged; it has no deployment receipt. The wrapper build and local visual review remain unverified. This work does not close the exact-artifact or independent-release-review gates, and nothing was submitted to CRAN.

## 11. Team Learning

Run route-retirement checks against final output after production cleanup because a raw pkgdown build can include routes that the deployment process removes.

Memory receipt: loaded the pigauto LOAD-FIRST manifest with `route.py pigauto` and searched the brain for current pigauto release decisions. No brain files were changed.

Golden Set: not run; no package behavior or statistical result changed.

## 12. Cross-Product Coverage

This check covers local static HTML output, sitemap targets, search targets, article-index links, and the homepage callout after cleanup. It does NOT cover visual appearance or accessibility, an actual Pages deployment, CRAN tarball contents, Windows checks, or CRAN acceptance.
