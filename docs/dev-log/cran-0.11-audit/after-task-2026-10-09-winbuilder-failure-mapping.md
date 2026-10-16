# After-task report: map Win-builder failures to the candidate patch, 2026-10-09

## 1. Goal

Check whether the five reported Win-builder test failures correspond to the changes in the separate Windows-portability patch, without modifying that patch or the frozen release archive.

## 2. Implemented

Compared the R-release and R-devel failure records with source commit `45bf230c7a10f56d0357b7fd5f5403a263eb72d4` in its isolated review worktree. All five reported failures have a direct source-level counterpart in that patch. This establishes a plausible, bounded mapping, not a Windows fix result.

## 3a. Decisions and Rejected Alternatives

The two accelerator failures arise in mocked `check_pigauto()` tests. The patch adds a mock for the new torch-runtime probe so those tests exercise their intended CUDA and MPS branches independently of the host's Lantern installation. The three Lantern errors occur in tensor tests. The patch changes their guard from checking only whether the `torch` R package exists to checking whether a tensor operation works.

Keep the patch in its separate lane. Its tests and package check were run on macOS, while the reported Windows logs have no source-archive checksum. Neither source correspondence nor macOS success binds the failure to the frozen archive or proves a Windows pass.

## 4. Files Touched

- Added this audit report.
- Updated `docs/dev-log/cran-0.11-audit/GATES.md` with the bounded mapping and remaining platform limitation.
- No package source, tests, frozen archive, or separate portability branch was changed.

## 5. Checks Run

- Read both retained compressed Win-builder logs. Each reports the same five failures: two expectations in `test-check-pigauto.R` (CUDA and MPS device selection), one Lantern error in `test-covariate-alignment.R`, and two Lantern errors in `test-monomorphic-discrete.R`.
- Compared those five sites with the four-file patch from `d76804e768bf77f430f43dcf591633f9cd900dab` to `45bf230c7a10f56d0357b7fd5f5403a263eb72d4`.
- The existing candidate-source receipt records a focused macOS test run and a full runnable package check. The latter reported 3,190 passes, 0 failures, 161 warnings, and 86 skips, then ended with 2 warnings because vignettes were omitted from its build. These are source-snapshot results, not Windows or frozen-archive checks.
- No new R tests or platform checks were run in this slice.

## 6. Tests of the Tests

The patch's CUDA and MPS tests keep their branch assertions while mocking the torch-installed probe. The three tensor tests retain their assertions when Lantern is usable and skip when only the R package is present. This is a source review of test routing; no absent-Lantern Windows runtime was available for direct confirmation.

## 7a. Issue Ledger

- Mapped: all five reported failure locations correspond to the separate patch.
- Open: the old logs are not checksum-bound to the 0.11.0 archive.
- Open: no Windows result exists for the patch source or a newly frozen archive.
- Open: the patch has not been merged into `main` and is not part of the frozen archive.

## 8. Consistency Audit

The frozen archive remains `pigauto_0.11.0.tar.gz`, SHA-256 `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4`, source commit `d76804e768bf77f430f43dcf591633f9cd900dab`. The candidate patch is a distinct commit and its separate source archive has a different hash. G8 and G9 remain open; the release verdict remains NOT READY.

## 9. What Did Not Go Smoothly

The first targeted search used a guessed provenance path that did not exist. The retained Windows result receipt identified the correct logs, and the follow-up search was limited to the five failure classes.

## 10. Known Residuals

This mapping does not show that all five failures are repaired on Windows, that either old Win-builder run used the frozen archive, or that the patch is safe to include without the owning lane's review. Exact-artifact platform checks and independent exact-artifact closure remain outstanding.

## 11. Team Learning

Matching a reported failure to a candidate test guard is useful triage, but it does not substitute for rerunning the corrected source on the affected platform with archive provenance.

## 12. Cross-Product Coverage

Covers the five reported test failures and their source-level mapping to a separate candidate patch.

Does not cover a Windows rerun, a patched exact artifact, deployment, merge, or CRAN acceptance.
