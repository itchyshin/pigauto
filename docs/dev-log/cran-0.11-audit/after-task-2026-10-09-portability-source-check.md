# After-task: full package check of the portability patch source

## 1. Goal

Broaden source-level verification of the Windows portability patch without editing its active branch or implying Windows validation.

## 2. Implemented

Built a temporary source archive from commit `45bf230c7a10f56d0357b7fd5f5403a263eb72d4` and ran the full runnable package test suite on macOS. Added raw logs, a hash-bound receipt, and the bounded result to the CRAN audit ledger.

## 3a. Decisions and Rejected Alternatives

The check used `--no-build-vignettes` to avoid a docs build in this source-only patch verification. Its two vignette warnings remain visible; the result is not called a clean CRAN check. The `--as-cran` attempt was retained as a network-blocked attempt, not treated as a test result. Assumption: the source patch commit remains the reviewed candidate; if its tests or implementation change, this receipt applies only to commit `45bf230c7a10f56d0357b7fd5f5403a263eb72d4`.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-portability-source-check.md`
- `docs/dev-log/cran-0.11-audit/provenance/winbuilder-portability-check-45bf230-2026-10-09.md`
- Four build, failed incoming-check, package-check, and test-output logs under `docs/dev-log/cran-0.11-audit/provenance/`
- No package source, tests, live site, or frozen release artifact was edited.

## 5. Checks Run

- `git archive 45bf230c7a10f56d0357b7fd5f5403a263eb72d4` created the temporary source snapshot.
- `R CMD build --no-build-vignettes .` succeeded. The resulting archive SHA-256 is `58e4b5628f0f58339fd49e50f5aa8dabed9ed8e635b7b0564d310a5fd185d405` and size is 5,052,462 bytes.
- `env -u NOT_CRAN _R_CHECK_FORCE_SUGGESTS_=true R CMD check --as-cran --no-manual --run-donttest pigauto_0.11.0.tar.gz` stopped at incoming feasibility because CRAN and Bioconductor indexes could not resolve. It did not run package tests.
- `_R_CHECK_FORCE_SUGGESTS_=true R CMD check --no-manual --run-donttest pigauto_0.11.0.tar.gz` completed with exit 0 and `Status: 2 WARNINGs`. The warnings are the expected absent vignette outputs from `--no-build-vignettes`. Testthat reported 3,190 passes, 0 failures, 161 warnings, and 86 skips.
- Environment: macOS arm64, R 4.6.0. Raw logs and checksums are listed in the provenance receipt.
- `slop_check.py`: 0 findings. `check-after-task.R`: required report structure passed; overall Unlazy verification remains red because release and other audit gates remain open. `git diff --check`: passed after preserving raw logs in compressed form.

## 6. Tests of the Tests

No tests were changed. The full runnable test suite in the temporary package archive passed. The changed portability test files were included in that suite; this host result does not prove the Windows-specific missing-runtime case.

## 7a. Issue Ledger

- Source-patch evidence strengthened: all runnable tests passed in the full package check.
- Open: Win-builder R-release and R-devel results for a final artifact, with hash-bound upload evidence.
- Open: libtorch-enabled Windows coverage for the tensor/model tests and real `check_pigauto(gnn = TRUE)`.
- Open: exact release artifact Windows failures remain unbound to the frozen archive SHA-256. G8 remains unmet.

## 8. Consistency Audit

The source archive was built from the exact reviewed patch commit and kept separate from the frozen release tarball. Raw logs were compressed without changing their bytes because their terminal output contains meaningful trailing whitespace; the receipt records both compressed and uncompressed hashes. The log distinguishes the network-blocked `--as-cran` attempt from the offline-capable package check and records the two vignette warnings caused by the build option.

## 9. What Did Not Go Smoothly

The first raw-directory check invocation stopped on derived DESCRIPTION fields. The corrected build-then-check flow reached package checks. The `--as-cran` command then stopped before installation because this environment could not resolve CRAN and Bioconductor indexes.

## 10. Known Residuals

This is macOS source-snapshot evidence with two warnings. It does not cover Windows runtime behavior, a libtorch-enabled Windows run, the frozen tarball, final site deployment, or CRAN acceptance. The separate README source PR still awaits the user’s choice on its exact title and body.

## 11. Team Learning

For package-source validation, build the archive before `R CMD check`; checking the raw directory can fail on fields synthesized during build. When network-dependent CRAN incoming feasibility is unavailable, preserve that failure and distinguish it from a standard package check.

## 12. Cross-Product Coverage

This covers the built package source snapshot and its runnable test suite on macOS. It does NOT cover Windows, libtorch-enabled Windows execution, the frozen release artifact, live website deployment, or CRAN acceptance.
