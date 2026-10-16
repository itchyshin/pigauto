## 1. Goal

Resolve the candidate-site verifier's unreproducible search-path expectation and validate the full fresh website built from the exact source correction under the deployed pkgdown version.

## 2. Implemented

- Changed the candidate-site verifier and fixture from 568 to 567 non-empty search-path references.
- Built the full PR #234 source snapshot with pkgdown 2.2.1 in a temporary R library, matching the version reported by the deployed Pages workflow.
- Retained the site build, crawl, unit-test, and strict-verifier outputs with checksums; added a provenance receipt and current release-ledger status.
- Retained a separate `pkgdown::check_pkgdown()` transcript, including R 4.6.0 and pkgdown 2.2.1, after the independent reviewer noted that its output was absent from the original build log.
- Confirmed the corrected Getting Started page has not deployed; Chrome still shows the old 300 × 4 statement for seven traits and 14 latent columns.

## 3a. Decisions and Rejected Alternatives

I did not change the expected count until I checked pkgdown version drift. Both the fresh pkgdown 2.2.0 and 2.2.1 builds produced 567 non-empty references; the retained historical receipt's 568 count was not reproducible. I kept the package behavior and user documentation unchanged because this issue was limited to the audit verifier and evidence.

## 4. Files Touched

- `script/cran-0.11-site/verify_candidate_site.py`
- `script/cran-0.11-site/test_candidate_site.py`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/site-search-index-candidate-2026-10-16.md`
- `docs/dev-log/cran-0.11-audit/provenance/site-build-2026-10-16-6ebc676-pkgdown-2.2.1.log`
- `docs/dev-log/cran-0.11-audit/provenance/site-crawl-2026-10-16-6ebc676-pkgdown-2.2.1.log`
- `docs/dev-log/cran-0.11-audit/provenance/candidate-site-verifier-tests-2026-10-16.log`
- `docs/dev-log/cran-0.11-audit/provenance/candidate-site-verifier-2026-10-16-pkgdown-2.2.1.json.log`
- `docs/dev-log/cran-0.11-audit/provenance/pkgdown-check-2026-10-16.log`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-16-site-search-verifier.md`

## 5. Checks Run

- Installed CRAN's pkgdown 2.2.1 source tarball (SHA-256 `9253a08bfbdddfa9f7f1404c5d1e5c9bc349845b078d06b68ef9288ece2db11c`) into a temporary library; no change to the regular R library.
- Full release-style site build from source commit `6ebc67636f31fa8ff324ddd4a327d5cd808fba99` completed. A separate retained `pkgdown::check_pkgdown()` run reported no problems; its transcript is SHA-256-bound in the receipt.
- Local crawl passed: 3,584 references checked, no errors.
- Strict candidate-site verifier passed: 44 retired routes absent, 62 sitemap entries matching 62 HTML pages, 613 search records, 567 non-empty paths, and no retired search target.
- Opened the deployed homepage, Getting Started, and MI article in Chrome. The MI page describes both model packages as optional; Getting Started remains out of date until PR #234 merges and deploys.
- At the last successful live check (2026-10-16, 02:08 UTC), PR #234 was open with Ubuntu release, Ubuntu devel, and macOS R-release checks queued. A fresh GitHub API refresh at 02:26 UTC failed to connect, so its current status remains unverified. This is a separate release gate and did not pass in this slice.

## 6. Tests of the Tests

All 16 candidate-verifier tests pass, including the positive fixture and the existing negative controls for retired routes, invalid discovery indexes, hidden warning callouts, external URLs, and unsafe paths. The test output is retained and SHA-256-bound in the provenance receipt.

## 7a. Issue Ledger

- Resolved for the current builder: hardcoded search-path expectation contradicted fresh output. The observed value is 567 under both pkgdown 2.2.0 and 2.2.1.
- Historical receipt: still records 568; the new dated receipt supersedes that count without altering the historical evidence file.
- Open: corrected documentation remains undeployed; final artifact, Windows/platform, network-capable `--as-cran`, and exact-hash reviewer gates remain unmet.

## 8. Consistency Audit

The verifier expectation, synthetic fixture, result log, provenance receipt, and ledger now agree on 567. Sitemap and HTML inventories both contain 62 routes. All non-empty search paths are in the sitemap, and retired paths are absent. The prior frozen tarball predates the documentation correction and remains ineligible as the final artifact.

## 9. What Did Not Go Smoothly

The first CRAN Archive URL returned 404; the current CRAN source path supplied pkgdown 2.2.1. Chrome blocked direct navigation to the live search JSON, so its raw production contents were not independently counted. The automated closeout template generator could not obtain `git status` because Git invokes Xcode and this workstation has not accepted the Xcode license. I did not accept that legal agreement. This report was written manually; current branch status remains unverified.

## 10. Known Residuals

The deployed homepage/site has not been rebuilt from the correction. There is no local full-site visual review in this slice. The candidate verifier checks that search paths are valid and retired routes are absent; it does not cover whether every distinct content page is independently searchable. No release readiness claim is made.

## 11. Team Learning

Pat's earlier independent source/render review passed the bounded Getting Started correction. Rose independently reran the 16 verifier tests and strict candidate check, recounted the candidate inventories, verified all retained hashes, and found no correctness blocker in the 568-to-567 correction. Rose identified one missing durable proof: the standalone `pkgdown::check_pkgdown()` output. I retained that transcript with its hash; the independent reviewer has not yet rechecked this added transcript. Memory receipt: used the project AGENTS instructions supplied in this task and its candidate-site verification requirements; the current LOAD-FIRST manifest was refreshed with `route.py`. Golden Set: not checked; this tooling correction was tested with its own negative-control suite. The project instructions prohibit writing to the brain vault without approval, so no brain files were changed.

## 12. Cross-Product Coverage

Covers: local candidate-site structure, retirement manifest enforcement, search-path validity, and verifier negative controls for the PR #234 source snapshot.

Does NOT cover: a deployed corrected site, full visual inspection, final source tarball, R-release/R-devel or Windows results bound to that tarball, CRAN incoming checks, or CRAN submission.
