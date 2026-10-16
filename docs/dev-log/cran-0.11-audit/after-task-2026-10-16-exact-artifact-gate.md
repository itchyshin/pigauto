# After-task: exact artifact and site evidence refresh, 2026-10-16

## 1. Goal

Refresh the exact-artifact/site evidence for the pigauto 0.11 release audit and prepare an unmerged evidence PR. Preserve unresolved platform status and do not submit to CRAN.

## 2. Implemented

- Updated the release gate record and JSON ledger for source commit `780c08470b263cc26b4139c711a5c144fd29389c`, exact tarball SHA-256 `0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d`, and the completed source-CI result.
- Retained fresh site build/crawl logs and passed the updated candidate-site verifier against the current build.
- Corrected the candidate search-path fixture expectation from 568 to 567 after the actual fresh site and prior clean build independently showed 567 path-bearing records among 613 search records.
- Kept `G8` and `G9` open: Win-builder results are still absent and the service acknowledged filename/size without a server-side digest.

## 3a. Decisions and Rejected Alternatives

- Kept the exact tarball frozen at the clean merged source commit; evidence-only changes do not modify the package artifact.
- Did not claim platform-clean from successful source CI or Win-builder upload acknowledgement.
- Did not infer a separate BirdTree redistribution licence; the ledger retains the maintainer warranty with citation/credit and its limitation.
- Assumptions made without asking: none.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/provenance/site-build-780c084.log`
- `docs/dev-log/cran-0.11-audit/provenance/site-crawl-780c084.log`
- `docs/dev-log/cran-0.11-audit/provenance/candidate-site-verifier-tests-780c084.log`
- `docs/dev-log/cran-0.11-audit/provenance/candidate-site-verification-780c084.log`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-16-exact-artifact-gate.md`
- `script/cran-0.11-site/verify_candidate_site.py`
- `script/cran-0.11-site/test_candidate_site.py`

## 5. Checks Run

- `python3 script/cran-0.11-site/test_candidate_site.py`: 16 tests passed.
- `python3 script/cran-0.11-site/verify_candidate_site.py ...`: passed; 44 retired routes, 62 sitemap entries, 613 search entries, 567 path-bearing search entries, and 62 HTML pages.
- Unlazy `--approve --reverify` on `GATES.md`: G5 and G5b passed; overall 9/11 met, G8/G9 open.
- `python3 .../cran_release_gate.py release-ledger.json`: `READY FOR CLAIMED RUNG` (tarball-clean only).
- Main source-CI run `38970604972`: Ubuntu R-release, Ubuntu R-devel, and macOS R-release passed. These checks are not exact-tarball platform results.
- Gmail search for current Win-builder result messages: none returned.
- PR #236 is open as a draft; its pkgdown job is intentionally skipped on pull requests by the repository workflow. The PR is unmerged.
- Staged source/document edits pass `git -c core.whitespace=-blank-at-eol diff --cached --check`. Raw R/test and pkgdown logs retain their original trailing whitespace and were kept byte-for-byte.

## 6. Tests of the Tests

The candidate-site verifier tests had previously failed against the real 613-entry index because the test expected 568 path-bearing records; inspection confirmed 567. The fixture and expectation were corrected, then all 16 verifier tests passed. The verifier was also run against the actual fresh site, not only fixtures.

## 7a. Issue Ledger

- Fixed: candidate-site verifier off-by-one expectation for path-bearing search records (568 expected, 567 observed).
- Open: Win-builder result arrival and checksum binding for the exact tarball.
- Open: final exact-artifact independent verdict and current platform-clean rung.
- Open: draft evidence PR #236 is published and intentionally remains unmerged for review.

## 8. Consistency Audit

Checked the live-site retirement/search/sitemap evidence against the fresh local site counts, the exact tarball hash and source commit, the updated ledger claims, and all three source-CI jobs. The user-facing site evidence remains bounded: Pat reviewed the NEWS exception and route receipts, not all 106 routes or all retained page contents. No package-source or installed-byte change was made in this evidence refresh.

## 9. What Did Not Go Smoothly

The first terminal invocation hit the Xcode Git shim's unaccepted-license message; the direct Command Line Tools Git binary worked without changing machine configuration. The repository closeout helper tried to resolve its output under the brain vault, outside this task's writable roots, so this report was written directly in the isolated evidence checkout. The route helper had no manifest for either temporary checkout; the supplied repository LOAD-FIRST manifest was used as the project routing source.

## 10. Known Residuals

The exact tarball is locally clean and its macOS `R CMD check --as-cran` result is `Status: OK`, but Windows results are not present or digest-bound. Exact-artifact platform coverage and the final independent panel remain open. Draft evidence PR #236 is open and unmerged. No CRAN submission or merge is claimed by this report.

## 11. Team Learning

Site discovery counts need two explicit denominators: total search records and records containing paths. A fixture should model both, including pathless entries, so a total-record count is not mistaken for the route-bearing count. Source-CI success and a platform upload acknowledgement remain separate from exact-tarball platform evidence.

Memory receipt: used the supplied pigauto LOAD-FIRST manifest (effective defaults, preserve `r_cal = 0`, use recovery-to-truth, and diff main before building) and the repo's release ledger; the applicable exact-artifact and platform-binding lessons shaped the rung claims. The brain route helper returned no manifest for the checkout paths. Golden Set: not run; this evidence-only refresh did not modify package behavior or exercise a known code-defect class.

## 12. Cross-Product Coverage

This follow-up changed the exact-artifact matrix runner and its documentation, plus retained evidence; it did not alter the package archive, package behavior, trait families, backend adapters, or defaults. It covers bounded exact-artifact optional-backend MI interoperability, including all four backend-availability states and termwise fixed-effect pooling, as well as the previously recorded local-site generation, retirement controls, and search/sitemap checks. It does NOT cover clean-room dependency installation, broad recovery or interval coverage, joint covariance/loading/variance-component pooling, checksum-bound Windows results, every retained public page's content, or CRAN acceptance.


## Follow-up: complete optional MI matrix receipt

Repaired the matrix runner receipt after independent review found that it omitted hashes for indirect R sources and did not record its working directory. Focused Python runner tests pass (2/2). Reran the exact frozen tarball SHA-256 `0fccf4bb0a179434fe15fa3a7ff514112bc5220e60c4b3e025fa00536f5ca29d`; all four cells and 12 steps passed, with unchanged archive hash throughout. Manifest and logs are retained under `provenance/mi-matrix-0fccf4bb-provenance-complete/`. The bounded claim and exclusions are recorded in `GATES.md` and the release ledger. Independent reviewer recheck passed with no material discrepancy after rehashing the artifact, runner, all five executed sources, and all 12 logs. Shared dependency symlinks mean this is not a clean-room installation test. This adds bounded adapter evidence only; platform-clean and overall release readiness remain unmet.
