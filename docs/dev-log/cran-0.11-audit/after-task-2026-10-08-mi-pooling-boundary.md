## 1. Goal

Add regression coverage for captured multiple-imputation fit failures that leave fewer than two poolable fits.

## 2. Implemented

Added a focused regression test through `with_imputations()` and `pool_mi()`. With two imputations and one synthetic fit failure, it checks the failure denominator, retained `n_fits` and `n_failed` metadata, captured `pigauto_mi_error`, the drop warning, and refusal to pool one remaining fit. No production code changed because existing behavior already enforced the minimum-two-fits boundary.

## 3a. Decisions and Rejected Alternatives

Kept this as a regression test rather than changing pooling logic. The audit found a missing test for a one-valid-fit outcome, not a defect in the established rule. The test uses `lm()` to exercise package-independent pooling and makes no claim about optional backend integration.

## 4. Files Touched

- `tests/testthat/test-mi-provenance.R`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-mi-pooling-boundary.md` (force-staged because the broad `docs/` ignore rule excludes new files)
- `.unlazy/mi-pooling-boundary/GATES.md` (ignored local acceptance ledger)

## 5. Checks Run

- `Rscript --vanilla -e 'devtools::test(filter = "mi-provenance")'`: PASS, 55 expectations, 0 failures, 0 warnings, 0 skips.
- Unlazy re-verification: PASS, 1/1 gate. The exact approved command ran `devtools::test(filter = "mi-provenance", stop_on_failure = TRUE)` and emitted `MI_PROVENANCE_PASS` with exit status 0.
- `git diff --check`: PASS.
- Independent read-only review confirmed the test covers the real `with_imputations()` to `pool_mi()` path and checks both failure-count attributes. The reviewer did not run tests; the commands above are the execution evidence.
- Chrome inspection showed PR #231 remains a draft with Ubuntu release and devel checks passing, macOS release still running, and pkgdown skipped. PR #228 remains a draft with conflicts and no checks. These are release-state observations, not validation of this uncommitted test.

## 6. Tests of the Tests

The synthetic second-fit failure exercises captured-error handling. Assertions require the one-of-two warning, failure metadata, dropped-error warning, and exact fewer-than-two stop. If the failure were hidden, counts changed, warning omitted, or pooling allowed one valid fit, at least one assertion would fail. No mutation test was run.

## 7a. Issue Ledger

- Fixed: added direct coverage for the boundary where two requested fits yield only one poolable fit.
- Open: real drmTMB and gllvmTMB object integrations remain separate checks in the release audit.
- Open: current Windows results are not verified against a final exact tarball.
- Open: the release evidence PR needs corrected merged and deployed source, a new frozen artifact, and conflict resolution.

## 8. Consistency Audit

The test is adjacent to existing MI provenance tests and exercises the `with_imputations()` and `pool_mi()` workflow rather than internal helpers. It checks that captured errors remain visible in result metadata and are not silently pooled. It does not change coefficient extraction, Rubin arithmetic, or optional model adapters.

## 9. What Did Not Go Smoothly

The closeout generator resolved a relative report path against the brain repository and created a template there. I verified that location and removed only the generated file, confirmed its absence, then recreated the report with an absolute pigauto worktree path. No brain note remains from that attempt.

## 10. Known Residuals

This test proves one failure boundary only. It does not prove external backend extraction, convergence or Hessian checks, saved-object portability, or recovery properties. At this report’s verification point, the test is staged but not committed, and branch CI has not run against a commit containing it. No tarball or CRAN gate is validated by this slice.

## 11. Team Learning

Memory receipt: loaded the repository `AGENTS.md`, pigauto LOAD-FIRST context, and unlazy acceptance instructions. The relevant invariant was that captured fit failures stay visible and pooling requires at least two valid fits. An independent reviewer checked the test scope.

Golden Set: not run because this narrow status-boundary regression does not match a recorded Golden Set class.

## 12. Cross-Product Coverage

Covers: package-independent `with_imputations()` failure capture, fit-count metadata, `pool_mi()` failure dropping, and the minimum-two-fits guard.

This does NOT cover: drmTMB or gllvmTMB installed-package integrations, saved and reloaded model objects, other extraction methods, broad Rubin arithmetic correctness, downstream inferential validity, or exact release-artifact checks.
