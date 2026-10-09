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

## Follow-up: exact-head CI ledger reconciliation, 2026-10-09

Two independent reviewers found that the PR descriptions had advanced beyond the durable ledger: GATES.md still stopped at run #37913439649 on `017e0a6`, and its G0 scope statement stopped at the earlier site commit. I updated GATES.md to record the passing #37916716258 receipt on `cf88d78` and the latest #37924160198 receipt on `dbe1248`, with the three platform durations and the candidate-source/force-Suggests limits. I also bound the G0 carry-forward to the exact six-file `cf88d78..dbe1248` audit-record delta. That delta contains no package, user-facing documentation, data, generated-help, or website-input changes.

Validation: `git diff --check` passed. The after-task structure check passed; its overall exit remains 1 because five unmet gates under `.unlazy/imputation-sim/gates/` belong to another acceptance ledger and were not changed here. The naturalness checker passed with zero findings. No package/site tests were rerun because this was an evidence-only correction. The PR description text previously verified in Chrome already records #37924160198 and retains the 8-of-11 tally. G7 (deployed-site verification), G8 (final post-merge tarball and checks), and G9 (independent final artifact review) remain open. Merge, deployment, and submission remain Shinichi's decisions.


## Follow-up: exact-head CI completion, 2026-10-09

Run #37927503719 completed successfully for PR #231 source head `7c2a3c361de64709864cfd2f6f63264cd91d2f11`, through synthetic merge `b9c86d8` into base `0b0f71f`. Ubuntu R-release, Ubuntu R-devel, and macOS arm64 R-release all reported `R CMD check Status: OK`; the full suite reported 3,196 passes, 0 failures, 175 warnings, and 83 skips on every platform. The macOS focused MPS test reported 201 passes, 0 failures, 50 warnings, and 0 skips in 273 seconds. The preserved consolidated log is 5,971,990 bytes with SHA-256 `1353b224fc8c7a4ba7ecf1936833e2d9b547d1fdc328eb3eeb8b5735b4d361f5`. This workflow used `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`; it is source-candidate evidence only and does not satisfy G8.

Files touched in this follow-up: `GATES.md`, this after-task report, and the retained consolidated workflow log `provenance/source-ci-37927503719.log`. The source delta after the previously checked `dbe1248` contains only audit-evidence records; no package behavior, bundled data, user-facing pages, generated help, or website inputs changed. G0–G6 remain met; G7–G9 remain open. No merge, deployment, or submission occurred. `git diff --check` passed. No package or site checks were rerun because this follow-up records completed CI evidence only.


## Follow-up: current-head CI completion, 2026-10-09

Run #37930221718 passed on PR #231 source head `2b56c64592e4920403fa5cdd78656cd805f2e0e5`, through synthetic merge `5b4f979` into base `0b0f71f`. Ubuntu R-release, Ubuntu R-devel, and macOS arm64 R-release all reported `R CMD check Status: OK`; each full suite reported 3,196 passes, 0 failures, 175 warnings, and 83 skips. The focused macOS MPS test reported 201 passes, 0 failures, 50 warnings, and 0 skips in 350.7 seconds. The retained consolidated log is 5,971,540 bytes with SHA-256 `d3d94372cbcc503f2d0dd4ec84ba577d075320bb7727121b83ef0726ac27b1e6`. This workflow set `NOT_CRAN=true` and `_R_CHECK_FORCE_SUGGESTS_=false`; it does not satisfy G8.

Files touched in this follow-up: `GATES.md`, this report, and `provenance/source-ci-37930221718.log`. These are evidence-only updates; no package code, data, user-facing pages, generated help, or website inputs changed. G0–G6 remain met and G7–G9 remain open. `git diff --check` passed for the two Markdown records; the verbatim raw log retains upstream trailing whitespace. No package or site checks were rerun because this follow-up records completed CI evidence only.
