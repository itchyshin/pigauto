## 1. Goal

Resolve issue #241's misleading held-out scoring workflow and draft an unposted maintainer reply.

## 2. Implemented

Added `predict(fit, evaluation = TRUE)`. It selects the fit's held-out baseline and masks validation/test values from prediction inputs. Default production prediction is preserved. New prediction objects record which latent indices are safe to score; `evaluate_imputation()` rejects scored indices outside that set. Corrected both vignette workflows, the extracted R scripts, and the monomorphic-trait test that used the same unsafe workflow.

## 3a. Decisions and Rejected Alternatives

Kept production as the default because ordinary imputation should use observed information. Used explicit provenance rather than guessing leakage from perfect scores, which can occur honestly for trivial traits. Retained the private masking argument for existing internal callers. Raw matrices and older objects remain caller-verified because their fitting history is unavailable.

Assumption: the request authorizes a local fix and a draft, not posting, merging, releasing, or deploying. The fix is based on current main `78b1306d`, after checking the intervening merge of #240.

## 4. Files Touched

- R/predict_pigauto.R and R/evaluate_imputation.R
- man/predict.pigauto_fit.Rd and man/evaluate_imputation.Rd
- vignettes/mixed-types.Rmd and vignettes/mixed-types.R
- vignettes/getting-started.Rmd and vignettes/getting-started.R
- tests/testthat/test-evaluation-mode.R and tests/testthat/test-evaluation-gnn.R
- tests/testthat/test-monomorphic-discrete.R
- tests/testthat/fixtures/cran-audit-exported-formals.tsv
- This directory: GOAL.md, AFTER-TASK.md, REPLY-DRAFT.md, assessments, and evidence logs.

## 5. Checks Run

The reporter's full Parts A, B and C reproduced on macOS R 4.6.0 against the reported source. The public evaluation argument then reproduced the masked scores in all three examples. The mixed-types article rendered with mass test RMSE 0.9284875, binary accuracy 0.5 and categorical accuracy 0.3333333. `pkgdown::check_pkgdown()` reported no problems. Documentation regenerated successfully.

Final full suite: 3,672 passed, zero failed, 184 warnings and eight skips. The separate package check, with vignette building and repeated tests disabled, had zero errors, two warnings about unbuilt vignettes and one online-clock note. These are development checks, not release certification. Commands and receipts are in CHECKS.tsv.

## 6. Tests of the Tests

On unchanged source, the new regression tests failed because production scoring was accepted, the public argument did not mask predictions, and invalid argument values were accepted. The GNN test forces the calibrated GNN weight to one and verifies that perturbing held-out truth by 100 leaves evaluation predictions unchanged. Formal metric comparisons verify agreement with `evaluate(fit)` for five trait types.

## 7a. Issue Ledger

Issue #241: fixed locally; GitHub remains unchanged. The published site and CRAN package still require a separate update. No new release was prepared.

## 8. Consistency Audit

Inspected internal evaluation, prediction baseline selection, input restoration, predictive mean matching donors, and conformal locality masks. Searched vignette and test scorer calls. The monomorphic-trait test was a sibling of the reported workflow and now requests evaluation mode. The API default inventory now includes the new argument. The closeout wrapper passed report structure, then failed on unrelated, inherited acceptance ledgers. Those ledgers were preserved.

## 9. What Did Not Go Smoothly

The initial test command used the original checkout by mistake and failed to locate the new test. It was rerun in the isolated copy before implementation evidence was accepted. The first broad run caught the API inventory mismatch and the sibling test workflow. Some progress messages were briefly misread as GNN errors; they were small-validation warnings, and the correction was stated in chat. Scratch copies were used for execution because the managed worktree was outside the writable sandbox; final files were copied selectively into the isolated worktree.

## 10. Known Residuals

Raw matrices and older prediction objects cannot supply trustworthy provenance. This work does not establish general calibration or inference validity. No website deployment, public reply, package release, or change to another lane was made.

CARRIED-OVER: local branch `codex/issue-241-evaluation`, uncommitted for review. Resume at the attached worktree and inspect its diff before integration. The original dirty checkout was preserved.

## 11. Team Learning

Memory receipt: loaded the routed pigauto manifest, R-package-engineer, validation-harness, systematic-debugging, test-driven-development and verification-before-completion instructions. The production/evaluation separation and gate-zero fallback shaped the fix. Brain search found no issue-specific prior decision. No brain vault writes were made.

Golden Set: not run; no hub guard or retrieval mechanism changed. Regression evidence directly targets the reported scoring failure.

## 12. Cross-Product Coverage

Covers: deterministic GNN-off mixed-type scoring for continuous, count, binary, categorical and ordinal traits; active GNN prediction input masking; production-object rejection; partial masks; old-object compatibility; missing splits and incompatible newdata; corrected vignette workflows.

This does NOT cover new statistical validation, multi-observation recovery, all eight trait families individually under the new argument, stochastic multiple-imputation evaluation, legacy saved-model reconstruction, release approval, or website deployment.
