# Pigauto 0.11.0 post-merge deployment and artifact evidence

Source: `bb5835d1b214d783da7b0b414df99aa6ba926bc7`. PR #226 was merged into
`main`; Pages run #570 and R-CMD-check run #37526171411 succeeded on that exact
merge commit. The audit tarball built from that source was checked on macOS.

The public site now serves the merged content. The home page has the four current
navigation labels and no per-trait benchmark links. The getting-started page shows
`phylo_signal_gate = FALSE`; the multiple-imputation article and
`multi_impute()` reference show the corrected current text. All 34 retired
`/dev/` HTML pages and four retired walkthrough URLs return HTTP 404. None of
the 34 retired benchmark routes appears in the 64-location sitemap or 574-path
search index. Exact live responses are in
`provenance/public-deployment-live-check-2026-10-06.tsv`.

The exact audit artifact is `pigauto_0.11.0.tar.gz`, SHA-256
`d34d469981386546b4ad03b9c7aab81674556708eeaf524c1c7365b6c68dec26`,
5,122,013 bytes and 247 archive members. Its DESCRIPTION reports version 0.11.0;
the package-path scan found no development, website, BACE, Git, macOS metadata,
or RStudio project paths. The artifact is retained locally at
`/private/tmp/pigauto-cran-011-gate-artifacts/pigauto-0.11.0-bb5835d-d34d4699/pigauto_0.11.0.tar.gz`;
the exact identity and full scan are in
`provenance/post-merge-artifact-identity.json`.

`R CMD check --as-cran --no-manual` on that exact tarball completed on macOS
Tahoe with R 4.6.0: `Status: OK`. Testthat reported 3,066 passes, zero failures,
161 warnings and 86 skips; the package check itself reported no ERROR or WARNING.
The retained full check log is
`provenance/post-merge-R-CMD-check-macos.log` (SHA-256
`657a19dc6584f1475f78f35d96824304be45b1279459f72bfa5f1c72113d41a0`);
the test output is retained alongside it. The merged source commit also passed
Ubuntu R release, Ubuntu R-devel, and macOS R release jobs in run #37526171411.
Those jobs are source-commit checks; the full check of the frozen tarball is the
macOS run above.

G6 is complete. G7 is partial: the exact frozen tarball has its full macOS
check, and the exact merged source commit passed Ubuntu release, Ubuntu devel,
and macOS release checks. The exact tarball is awaiting Windows R-release and
R-devel results from win-builder. G8's bounded audit
criteria are now complete: all three fresh reviewers returned READY after the
cache-busted public-page check, the CRAN release-gate selftest passed its planted
negative controls, and artifact/check hashes were re-verified. The review
outcomes, control output and live-page recheck are retained in
`provenance/review-panel-2026-10-06.md`,
`provenance/release-gate-selftest-2026-10-06.log`, and
`provenance/post-review-live-page-recheck-2026-10-06.md`.

This closes the audit evidence gate only. BirdTree data redistribution rights
remain unresolved, so the release ledger and CRAN authorization stay
`NOT_READY`; no upstream contact or CRAN submission has occurred.
