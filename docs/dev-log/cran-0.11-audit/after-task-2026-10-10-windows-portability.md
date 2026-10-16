## 1. Goal

Make the Windows test suite tolerate an installed torch R package without the libtorch runtime, while preserving tests on systems where the runtime works.

## 2. Implemented

Extracted the existing guarded torch installation check into a private helper so accelerator-selection tests can mock that prerequisite. Changed three tests that call torch tensor operations to use the existing libtorch-aware skip helper.

## 3a. Decisions and Rejected Alternatives

Kept torch as an imported runtime package and retained these tests when libtorch works. The default `gnn = FALSE` paths bypass torch, so their tests remain unguarded. This patch changes test portability only; runtime availability messages stay as they are. A different Windows failure could change that scope.

## 4. Files Touched

- `R/check_pigauto.R`
- `tests/testthat/test-check-pigauto.R`
- `tests/testthat/test-covariate-alignment.R`
- `tests/testthat/test-monomorphic-discrete.R`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-10-windows-portability.md`

## 5. Checks Run

- Focused red check before implementation: `devtools::test(filter = "check-pigauto")` failed in the two accelerator-selection tests because the runtime prerequisite had no mockable binding.
- Focused post-change tests: `devtools::test(filter = "check-pigauto|covariate-alignment|monomorphic-discrete")` passed 140 expectations, with 11 fixture warnings and no failures.
- Full suite: `devtools::test()` passed 3,556 expectations, with 182 warnings and 8 skips, in 244.8 seconds.
- `devtools::document()` completed without generated help changes.
- `git diff --check` passed.

## 6. Tests of the Tests

The pre-change red run failed at the intended point: testthat could not mock `.check_pigauto_torch_is_installed` because the seam did not exist. After extracting that helper, the focused tests passed. The no-libtorch behaviour still needs confirmation on Windows; this Mac has a working libtorch runtime.

## 7a. Issue Ledger

- Fixed in this candidate: mock accelerator tests depended on the actual libtorch-installation result, and three tensor-operation tests skipped only when the torch R package was absent.
- Deferred: validate R-release and R-devel Win-builder runs against a new exact artifact.

## 8. Consistency Audit

Reviewed the existing `skip_if_no_libtorch()` helper and adjacent tests. The three changed tests execute normally when a real tensor probe succeeds. The later monomorphic imputation tests remain unguarded because the default `gnn = FALSE` path bypasses torch. An independent code review found no required correction.

## 9. What Did Not Go Smoothly

The first closeout-generator invocation resolved paths relative to the brain repository rather than this linked worktree and failed before writing a file. The report was then created directly in the pigauto worktree. The structural after-task check passed. The full closeout check also listed five open gates under `.unlazy/imputation-sim/gates/`; they name the separate `pigauto-imputation-sim` worktree, so I left them untouched.

## 10. Known Residuals

Release evidence still requires a Windows result on a newly frozen artifact. This patch leaves the existing tarball unchanged. The global closeout gate remains open because the separate imputation-simulation ledger has five unmet gates.

## 11. Team Learning

Memory receipt: loaded the pigauto `LOAD-FIRST` manifest, R package guidance, validation harness, test-driven development, and after-task instructions. The prediction-correctness and optional-runtime guidance shaped the decision to preserve runtime behaviour and guard only backend-dependent tests. Golden Set: not run; no known-mistake class beyond the reported test-environment mismatch was in scope.

## 12. Cross-Product Coverage

Covers the pigauto R test suite and the candidate source changes. Does NOT cover Windows execution, a rebuilt release tarball, deployment, publication, or CRAN submission.
