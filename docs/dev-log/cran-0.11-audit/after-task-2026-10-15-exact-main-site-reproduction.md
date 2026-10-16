# After-task: exact-main local site reproduction and evidence repair

## 1. Goal

Reproduce the site checks from the exact merged `main` source and reconcile the CRAN 0.11 site evidence without freezing or submitting an artifact.

## 2. Implemented

- Rebuilt the full offline site from a `git archive` of `origin/main` at `8082ef8273177ff664d49ac4bd643f7def671f20`.
- Preserved the source marker, build/crawl/pkgdown-check/verifier logs, rendered-content assertions, and checksum-bound receipt.
- Reconciled the release ledger and gate narrative. The earlier receipt's ledger hash did not match its bytes; it is now explicitly superseded. The current build's actual pkgdown version is recorded as 2.2.0.
- Confirmed the three pages in the prior Chrome review are byte-identical in the exact-main reproduction; the review is bounded to those three pages.
- Incorporated independent Rose and site-auditor reviews: both verified the source marker, receipt/log hashes, current pkgdown identity, and bounded page-hash transfer. Rose flagged two stale narrative fields; both are corrected. The independent site audit found no other current-field contradiction.
- Committed the reviewed receipt and checkpoint locally on `codex/cran-011-gate-refresh`. Push attempts via SSH (twice) and HTTPS (once) failed with GitHub `remote: fatal: early EOF` / `index-pack failed`; a read-only remote-ref check found no branch, so the checkpoint remains local. No final evidence PR was opened because deployment and exact-artifact gates remain unmet.

## 3a. Decisions and Rejected Alternatives

- Kept the previous receipt as historical evidence and added a new receipt with an external ledger hash rather than preserving a self-hash field.
- Did not claim a new Chrome review: Chrome blocked both `file://` and localhost for the local build. The earlier three-page review transfers only because the three rendered HTML files are byte-identical.
- Did not treat the queued Pages or R-CMD-check runs as passing, and did not freeze a tarball while deployment is stale.
- Assumptions made without asking: none.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `.unlazy/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/exact-main-site-reproduction-2026-10-15.json`
- `docs/dev-log/cran-0.11-audit/provenance/site-build-exact-main-8082ef8-2026-10-15.log.gz`
- `docs/dev-log/cran-0.11-audit/provenance/site-crawl-exact-main-8082ef8-2026-10-15.log`
- `docs/dev-log/cran-0.11-audit/provenance/pkgdown-check-exact-main-8082ef8-2026-10-15.log`
- `docs/dev-log/cran-0.11-audit/provenance/candidate-site-tests-exact-main-8082ef8-2026-10-15.log`
- `docs/dev-log/cran-0.11-audit/provenance/candidate-site-verifier-exact-main-8082ef8-2026-10-15.log`
- `docs/dev-log/cran-0.11-audit/provenance/rendered-content-check-exact-main-8082ef8-2026-10-15.log`
- `docs/dev-log/cran-0.11-audit/provenance/site-source-commit-exact-main-8082ef8-2026-10-15.txt`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-15-exact-main-site-reproduction.md`

## 5. Checks Run

- Fresh offline build from archived `origin/main` `8082ef8`: 62 HTML pages and 62 sitemap entries.
- Crawler: `SITE_CRAWL_OK`, 3,584 local references, 34 retired pages, zero errors; its six crawler tests passed.
- Candidate verifier: 44 retired routes absent, 613 search records, 567 non-empty paths, `CANDIDATE_SITE_OK`.
- Candidate verifier tests: 16 passed.
- Reverification of the existing unlazy commands: retired-route verifier tests 3 passed; candidate verifier tests 16 passed; release-ledger JSON parsed; exact-main candidate verifier returned `CANDIDATE_SITE_OK`.
- `pkgdown::check_pkgdown()`: “No problems found.” Environment: R 4.6.0, pkgdown 2.2.0, Pandoc 3.12.
- Four rendered-content assertions passed: 300 × 14 dimensions, partial-overlap behavior, and absence of the stale 300 × 4 / dropped-species claims.
- The three reviewed article HTML files were byte-identical to the prior Chrome-reviewed build.
- The pkgdown build transcript is stored losslessly as deterministic gzip because its raw output contains two trailing-space lines; decompression reproduces the reviewed 32,940-byte log exactly.
- Independent review confirmed the checksum-bound receipt and referenced log hashes. Reviewers explicitly limited visual coverage to the three byte-identical pages and kept deployed-site verification open.
- Fresh Chrome checks at 2026-10-15 22:44 MDT: R-CMD-check run `38843648553` remained at 0/3 jobs and Pages run `38843648852` remained queued after more than 90 minutes; GitHub Status reported Actions operational. The live Getting Started page still rendered the stale `X_scaled (300 × 4)` output and unmatched-species-dropped wording.
- Unlazy read-only status: 4 gates met, 3 unmet (G3 exact-artifact R-release/R-devel, G4 exact artifact and deployed-site independent approval, G5 unmerged final evidence PR). G7 now points to the exact-main output and its checked command passes.

## 6. Tests of the Tests

The candidate-site suite retains negative controls for invalid discovery targets and retirement states. The earlier verifier's 568-path assumption was shown to fail against the actual 567-path output; the current 567-path verifier passes the exact-main site.

## 7a. Issue Ledger

- Fixed: receipt hash mismatch and incorrect current pkgdown-version claim; both are superseded by the checksum-bound exact-main receipt.
- Confirmed: local structural/site checks pass for exact `main` source `8082ef8`.
- Open: deployed-site G7, queued main checks/deployment, frozen exact-artifact G8, and final independent release closure G9.
- Deferred: no additional site or package changes were made.

## 8. Consistency Audit

Checked the exact source marker, generated pkgdown config, build/crawl/check logs, verifier output, retired manifest, rendered Getting Started text, live page state, and the historical Chrome-review page hashes. Rose identified and we corrected the stale “below” checkpoint reference and pre-merge PR #234 wording in the panel field. The site auditor independently found no other current-field contradiction. The current evidence consistently separates local build success from the stale deployed site. The fresh live check confirms the deployment queue is still the concrete G7 blocker.

## 9. What Did Not Go Smoothly

The first build lacked Pandoc on `PATH`; adding the installed `/opt/homebrew/bin/pandoc` enabled the successful run. Local Chrome access to file URLs and localhost was blocked by browser policy. The brain closeout generator resolves paths from its own repository root, so its initial template landed in the adjacent directory; that empty file was removed and no brain-vault content was retained or changed. The system Git shim was blocked by the unaccepted Xcode licence; the bundled Git runtime worked without changing that setting. The raw build log contains two trailing-space lines, so the committed evidence stores a deterministic lossless gzip copy and checks its decompressed hash. GitHub rejected three push attempts with `early EOF` while unpacking the pack; no remote ref was created.

## 10. Known Residuals

The local build used pkgdown 2.2.0. No new full-site browser review was possible; only the three byte-identical pages inherit the earlier visual review. Live deployment, queued CI, checksum-bound platform checks, a final tarball, and CRAN submission remain outstanding. The queue cause is unknown. The working tree also retains unrelated pre-existing untracked `graft/` and `script/cran-0.11-site/__pycache__/` paths, which were not staged.

## 11. Team Learning

Bind site receipts to an exact source marker and record the actual generated pkgdown version. Store the receipt hash outside the receipt itself. A rendered-page review should name the pages inspected and should never imply whole-site visual coverage.

Memory receipt: Loaded the routed pigauto manifest and the prior CRAN 0.11 release-audit memory. Exact source binding, G7 separation, and no-tarball-before-deployment shaped this work. Golden Set: not in scope; no package behavior or scientific default changed.

## 12. Cross-Product Coverage

Covers: exact-main local pkgdown build, sitemap/search structure, retired-route checks, rendered Getting Started text, and a bounded three-page visual review by byte identity.

Does NOT cover: successful live deployment, current live sitemap/search/retired-route checks, Windows/macOS release checks on a final tarball, all-page visual/accessibility review, statistical recovery, final independent exact-artifact review, or CRAN submission.
