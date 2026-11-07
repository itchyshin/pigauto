## 1. Goal

Address the ten open pigauto issue threads, verify actionable repairs and
reply to each reporter for confirmation. Leave issue closure to the authors.

## 2. Implemented

Shared observed-row rank handling fixes factor-covariate lambda estimation
and prediction. Unavailable estimates retain a finite replay value with status
metadata. Prediction rejects unused arguments. CSV ordinal input accepts an
explicit level order, and data-only species warnings explain retained rows.

The README and articles distinguish released evaluation from development,
explain repeated observations and uncertainty, and add worked missingness
plots. Baseline reports show offline metrics and conformal scores, including
bare fits. Daniel Noble's accessor and quick-start commits are preserved.

Draft PR 262 contains the implementation. All ten replies were posted and
read back, and all ten issues were still open at 22:56:45 UTC on 7 November.

## 3a. Decisions and Rejected Alternatives

Keep full one-hot covariates for the GNN and remove aliases only in the
baseline design. Warn on unidentified missing-row effects rather than assert
that every extrapolation is coding-invariant. Preserve the legal zero gate.

Do not replace the imputation model, invent a cause for private-data ties or
interval narrowing, submit to CRAN, or rerun a broad benchmark campaign.
CI expansion in PR 255 remains separate. PRs 252 and 253 remain open.

Assumption made without asking: interpret the closure instruction as asking
reporters to close issues if satisfied. This was stated before work began.

## 4. Files Touched

- `DESCRIPTION`
- `NAMESPACE`
- `NEWS.md`
- `R/bm_internal.R`
- `R/check_pigauto.R`
- `R/fit_baseline.R`
- `R/impute.R`
- `R/multi_impute.R`
- `R/pagel_lambda.R`
- `R/plot_missingness.R`
- `R/predict_pigauto.R`
- `R/preprocess_traits.R`
- `R/read_traits.R`
- `R/report.R`
- `R/result_accessors.R`
- `README.md`
- `_pkgdown.yml`
- `docs/dev-log/issue-review-20261107/AFTER-TASK.md`
- `docs/dev-log/issue-review-20261107/AFTER-TASK.md.assessment.json`
- `docs/dev-log/issue-review-20261107/PLAN.md`
- `docs/dev-log/issue-review-20261107/PLAN.md.assessment.json`
- `docs/dev-log/issue-review-20261107/REPLIES.md`
- `docs/dev-log/issue-review-20261107/REVIEW.md`
- `docs/dev-log/issue-review-20261107/REVIEW.md.assessment.json`
- `docs/dev-log/issue-review-20261107/ci-status.json`
- `docs/dev-log/issue-review-20261107/final-check.log`
- `docs/dev-log/issue-review-20261107/posted-comments.json`
- `docs/dev-log/issue-review-20261107/replies-assessment.json`
- `docs/dev-log/issue-review-20261107/smoke.R`
- `docs/dev-log/issue-review-20261107/smoke.log`
- `docs/dev-log/issue-review-20261107/source-test-summary.json`
- `man/check_pigauto.Rd`
- `man/completed_data.Rd`
- `man/fit_baseline.Rd`
- `man/get_imputed_datasets.Rd`
- `man/impute.Rd`
- `man/multi_impute.Rd`
- `man/pigauto_report.Rd`
- `man/plot_missingness.Rd`
- `man/predict.pigauto_fit.Rd`
- `man/read_traits.Rd`
- `tests/testthat/fixtures/cran-audit-exported-formals.tsv`
- `tests/testthat/test-community-surface.R`
- `tests/testthat/test-cran-audit-defaults.R`
- `tests/testthat/test-docs-phytools-warning.R`
- `tests/testthat/test-issue-input-clarity.R`
- `tests/testthat/test-issue-visuals.R`
- `tests/testthat/test-lambda-covariates.R`
- `tests/testthat/test-lambda-default.R`
- `tests/testthat/test-result-contract.R`
- `vignettes/common-pitfalls.Rmd`
- `vignettes/getting-started.Rmd`
- `vignettes/mixed-types.Rmd`
- `vignettes/multiple-imputation.Rmd`
- `vignettes/visualisation.Rmd`

## 5. Checks Run

Final focused tests passed, including the bare-fit metric comparison,
full-dummy/treatment-design equivalence, sparse fallback replay, explicit
ordinal levels and prediction argument rejection. The public n=200, seed=1
smoke gave lambda 0.5015374 with both encodings and completed-value difference
2.47e-10; five data-only rows stayed unchanged. See smoke.R and smoke.log.

The full pkgdown site built and check_pkgdown reported no problems. The new
vignette rendered and its plot labels were inspected; the offline report SVG
was rendered separately and inspected. Pandoc printed a MathJax deprecation
warning. Documentation regeneration and git diff --check passed.

Final source suite at implementation commit 7e542057: 854 tests, 3,761
passes, zero failures/errors, 184 warnings and eight skips. The warnings and
skips are retained; they are not passes. The full pre-final-report package
check completed in 10m 5s with zero errors, warnings or notes. The final rebuilt
check (tests disabled because the final suite ran separately) completed in
2m 55.9s with zero errors, warnings or notes. Both built vignettes. Exact check
logs are retained with this report. GitHub CI at 23:04:41 UTC remains in progress for Linux release, Linux
devel and macOS release on implementation commit 7e542057. The pkgdown
PR job is deliberately skipped. No CI success is claimed. The named unblock
event is completion of run 42422826371; the next eligible status check is
23:19 UTC unless a completion or failure notification arrives sooner.

## 6. Tests of the Tests

Before repairs, the factor design regression failed on the singular objective,
its shape and SE equivalence. The input tests failed for the absent ordinal
argument and silent prediction typo. Sparse lambda replay failed before the
finite fallback repair. Report comparison tests catch absent held-out metrics;
visual tests cover one-trait inputs and restored graphics settings.

## 7a. Issue Ledger

- 259: rank repair and unestimated-status handling proposed; reporter asked to test.
- 261: CRAN evaluation workaround and version distinction proposed; CRAN release outstanding.
- 258: warning corrected; data-only row preservation checked.
- 257: four missingness views and worked guide proposed.
- 254: README, report, ordinal input and draw-role clarifications proposed.
- 256: point, draw-distribution and downstream-inference questions answered.
- 250: posterior fields documented; private-data ties await a shareable reproducer.
- 251: single-fill distinction clarified; narrower correlation intervals remain unexplained.
- 260: multi-observation and seed/interval limitations clarified.
- 135: existing scoped evidence linked; broad comparator and graph-specific studies remain open.

Posted comment URLs are retained in posted-comments.json. No issue was closed.

## 8. Consistency Audit

Read all issue bodies/comments, existing overlapping PRs and current main before
editing. Checked prediction replay, full-rank and aliased covariate designs,
export/manual/formal inventories, navigation, novice examples, report paths
and uncertainty wording. Independent reviewers found and helped correct the
missing-only-level, sparse replay, plot-input, clipping and bare-fit gaps.
The final statistical reviewer found no remaining P1 or P2 source finding.

## 9. What Did Not Go Smoothly

The original checkout was dirty and on an old handover branch, so work used an
isolated managed checkout. The initial preflight invocation used a repository
name where a directory was needed; the corrected invocation succeeded.
Connector writes returned explicit 403 errors; authenticated CLI writes
succeeded and connector reads independently verified them. One smoke invocation
used the wrong report argument name and was corrected before its retained run.

The initial full suite retained two assertions for removed CRAN 0.10 guidance.
Those assertions were updated to test the current version mismatch and passed
focused verification. A fresh final-source full suite was launched.

## 10. Known Residuals

The PR stays draft while CI is running. The final evidence-only commit skips
CI because no package source, tests, manuals or public site inputs changed;
the implementation run remains identified by its exact source commit.

The patch is proposed in a draft PR, not merged, deployed or released to CRAN.
Private mammal ties, interval narrowing and broad external-comparator evidence
remain unresolved. GNN report charts still depend on their existing CDN.
Unidentified missing-only factor effects are warned about, not identified by
this repair. Optional route-status reconciliation has no demonstrated failure.

## 11. Team Learning

Memory receipt: loaded the canonical contract, codex-efficiency procedure,
repository manifest from route.py, R engineering, validation, debugging,
planning, review and verification guidance. Prediction correctness, scoped
scientific claims, ownership and the zero-gate fallback shaped this work.
Golden Set: not run; relevant repository evaluation and replay regressions
were used. No brain-vault files were written. Daily efficiency audit found
compaction breaches in other recorded sessions; it did not establish such a
breach for this chat. Model receipts are in REVIEW.md; no savings claim is made.

## 12. Cross-Product Coverage

Covers: baseline factor/numeric design rank, sparse baseline replay, prediction
argument handling, CSV ordinal order, input-row warning, missingness views,
baseline report metrics, documentation and all ten issue replies.

Does NOT cover: BACE source, new posterior/GNN campaigns, scientific interval
calibration, general downstream inferential validity, private-data diagnosis,
CRAN submission, a production deployment, or closure of any issue.
