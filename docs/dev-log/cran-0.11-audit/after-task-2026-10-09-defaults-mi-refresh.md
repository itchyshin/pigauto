## 1. Goal

Refresh PR #228's G1/G2 receipts with current-source effective-default and optional-model multiple-imputation evidence. Keep the final artifact and deployment gates open.

## 2. Implemented

Updated G1 to distinguish inventory-wide declared-formal/default matching from selected runtime-route checks, recorded the latest 198-pass defaults run and current 234-pass focused MI suite, and marked older receipts historical. Updated G2 with fresh source-mode real-object integration for drmTMB and gllvmTMB. Preserved all three raw logs with SHA-256 hashes.

## 3a. Decisions and Rejected Alternatives

Followed Shinichi's direction that the checkouts belong to one pigauto audit lane. The branch/worktree is the evidence checkout for this slice, not a competing lane. Kept drmTMB and gllvmTMB optional, retained automatic extraction adapters and termwise fixed-effect pooling, and recorded the earlier installed-backend matrix as supporting evidence. Did not claim joint covariance pooling, broad inferential validity, scientific optimality of defaults, or exact-tarball verification. Did not merge, deploy, or submit.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/defaults-inventory.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-defaults-mi-refresh.md`
- `docs/dev-log/cran-0.11-audit/provenance/defaults-current-2026-10-09-cf88d78.log`
- `docs/dev-log/cran-0.11-audit/provenance/mi-real-current-2026-10-09-cf88d78.log`
- `docs/dev-log/cran-0.11-audit/provenance/mi-focused-current-2026-10-09-cf88d78.log`

## 5. Checks Run

- `python3 ~/shinichi-brain/tools/route.py pigauto` loaded the repository LOAD-FIRST manifest.
- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla script/cran-0.11-defaults/run.R`: 198 passed, 0 failed, 0 warnings, 0 skipped; `CRAN_DEFAULTS_AUDIT_OK` on source commit `cf88d78e37d08130bd33823e1fd4594f3afa2d4f`.
- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla script/cran-0.11-integration/check-adapters.R real`: 50 expectations passed across six actual fits, with convergence 0 and `pdHess = TRUE` for all fits; `REAL_ADAPTERS_OK`.
- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla -e 'devtools::test(filter = "multi-impute-posterior|multi-impute-trees", stop_on_failure = TRUE)'`: 234 passed, 0 failed, 9 expected small-sample warnings, 0 skipped.
- Retained log hashes: defaults `94596524c0dfb08da789080ff8609eee282a789b4aad5be3a83709d6bcb445de`; real adapters `6f324f0429683b03af062fdced124c907e8d7d34ce74987635bcb9f93585ee94`; focused MI `71006a3a975e48d0f90a7a3cfc99908bfff534530305a090c01728e0f1b19b1b`.
- G1's current evidence distinguishes inventory-wide formal/default comparison from runtime verification of named critical routes. Older 860/905 results remain unreconciled historical receipts and are not combined with the current focused-suite count. The after-task structure checker passed; its full closeout then listed five unmet `.unlazy/imputation-sim` gates outside this CRAN audit.

## 6. Tests of the Tests

The defaults runner's assertion count is higher than the prior receipt and passed all current assertions. The real-backend harness fits both model classes and compares against independent Gaussian and Rubin-rule calculations; it is not a mock-only adapter test. The focused MI suite exercises posterior and tree routes with small fixtures. No deliberate fault injection was performed in this evidence-refresh slice.

## 7a. Issue Ledger

Resolved for this slice: G1 latest defaults receipt and G2 current source-mode object integration are now attached to PR #228's evidence ledger.

Still open: reconcile the historical 860-versus-905 count difference if it affects a future claim; exact final tarball validation; deployed-site verification; checksum-bound platform results; independent final-artifact verdict.

## 8. Consistency Audit

Checked `GATES.md` against the current source commit and raw logs. Confirmed the current runner tests selected effective routes while the inventory checks declared formals/defaults; the ledger no longer implies all argument interactions were executed. Preserved the distinction between `multi_impute()` automatic posterior/conformal routing and `multi_impute_analysis()` analysis-aware `lm`/`bayes_norm` defaults. Preserved the source-check versus exact-tarball boundary and termwise pooling scope. The three copied log hashes match their originating files.

## 9. What Did Not Go Smoothly

The shared preflight census counts each active checkout as a lane, while Shinichi clarified that these checkouts all belong to the same pigauto driver. I followed that ownership decision and claimed a narrow lease on the evidence files. The brain closeout helper resolves its report root to the parent Shinichi repository, so its `new` command created a blank template there instead of in pigauto. Auto-review blocked removing that out-of-scope file; it remains untouched. The correct pigauto report was authored directly in the release-evidence checkout. The report structure check passes, but full closeout is withheld by five unmet `.unlazy/imputation-sim` gates outside this audit.

## 10. Known Residuals

This slice does NOT establish exact-tarball installation, deployed website behavior, Windows/macOS results bound to the final archive, or CRAN acceptance. The historical 860/905 count difference is still unexplained. The earlier installed backend matrix is not rerun in this slice. No new recovery campaign was run.

## 11. Team Learning

When the maintainer identifies multiple checkouts as one project lane, record that ownership and use exact-path leases to protect the evidence files. In release ledgers, distinguish complete inventory comparison from runtime route coverage, and keep old counts visibly historical until their source snapshots are reconciled.

Memory receipt: loaded the pigauto LOAD-FIRST manifest with `route.py`; its release-boundary, evidence, and lane guidance shaped this update.

Golden Set: no package source behavior changed; no source known-mistake class was in scope.

## 12. Cross-Product Coverage

Covers current-source default inventory comparison, named effective routes, optional drmTMB/gllvmTMB fixed-effect MI extraction and pooling, and focused posterior/tree MI tests.

Does NOT cover every interaction among public formals, scientific optimality of all defaults, joint covariance pooling, broad inferential validity, final-artifact installation, deployed-site behavior, platform release checks, or CRAN submission.
