## 1. Goal

Bind current release evidence to one frozen pigauto 0.11.0 archive, refresh defaults, optional-model MI, website, and Windows checks, and keep the draft release PR unmerged.

## 2. Implemented

Updated the `simulate_benchmark()` help example to use `gnn = FALSE` so it works when optional torch is absent; the public function default remains `TRUE`. Added a regression assertion requiring nonempty benchmark output. Refreshed the release ledger and gate record for the frozen archive, current defaults audit, exact-hash MI matrix, R-hub Windows checks, and fresh local site review. The candidate remains `NOT_READY` because the deployed site and final release-validator closeout are open.

## 3a. Decisions and Rejected Alternatives

Kept the public `gnn = TRUE` default unchanged; the failure was in a documentation example's runtime assumptions, and the baseline-only example demonstrates that path without optional torch. Rebuilt and rechecked a new frozen archive after the example/test change rather than reusing older hash-bound results. Kept the live deployment gate open because source commit `9a457eb` has not been merged or deployed. Did not merge PR #236 or submit to CRAN.

Assumptions made without asking: treated the user's public-upload authorization as covering the R-hub public package upload; it does not authorize a merge, deployment, or CRAN submission. Kept the current package version at 0.11.0 based on the prior publication-history check recorded in the ledger.

## 4. Files Touched

- `R/simulate_benchmark.R`
- `man/simulate_benchmark.Rd`
- `tests/testthat/test-new-features.R`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/provenance/defaults-current-2026-10-17-9a457eb.log`
- `docs/dev-log/cran-0.11-audit/provenance/site-9a457eb-review.md`
- `docs/dev-log/cran-0.11-audit/provenance/windows-rhub-c2e768ad-2026-10-17.md`
- `docs/dev-log/cran-0.11-audit/provenance/mi-matrix-c2e768ad/manifest.json` and its 12 step logs
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-17-current-artifact-closeout.md`
- Frozen artifact receipts and logs: `exact-artifact-c2e768ad-source-binding.md`, `exact-artifact-c2e768ad-local-check.md`, `exact-artifact-c2e768ad-as-cran-check.log`, and `exact-artifact-c2e768ad-as-cran-testthat.Rout`

## 5. Checks Run

- `devtools::test(filter = "new-features")`: 112 expectations passed, 0 failures.
- One-seed `simulate_benchmark(..., gnn = FALSE)` smoke: completed and returned nonempty results.
- Current-source `devtools::test(filter = "cran-audit-defaults")`: 202 passed, 0 failed, 0 warnings, 0 skips; retained log hash `4c4d20917e712f09e93fab8bb62707edfce28c71ddae70af5c71df5a85d1cc99`.
- Exact frozen archive `R CMD check --as-cran --run-donttest --timings`, `NOT_CRAN` unset, force-Suggests enabled: `Status: OK`; testthat 3,195 pass, 0 fail, 161 warnings, 86 skips.
- Exact-hash optional-backend MI matrix: 12/12 steps passed across neither, drmTMB-only, gllvmTMB-only, and both-installed cells; hash unchanged before and after each step.
- R-hub Windows R-release 4.6.1 and R-devel r90686: both `Status: OK`, each 2,442 pass and 0 failures after R-hub rebuilt the source extracted from the byte-identical uploaded object.
- Fresh local pkgdown build and crawler: 62 HTML pages, 614 search entries, 3,588 references, 34 retired routes, zero crawl errors; `pkgdown::check_pkgdown()` clean.
- Browser content/accessibility review of the rendered `simulate_benchmark()` page completed. The independent site review confirmed retired route output and discovery checks.
- Unlazy `--reverify` approved and reran four read-only checks: current defaults, exact-hash MI manifest, current site build, and current crawl all passed. The repository ledger has 7 met gates and 4 open gates: local pixel-level review, deployed-site verification, post-merge artifact freeze, and final independent closure.
- `cran_release_gate.py release-ledger.json`: `READY FOR CLAIMED RUNG` at `platform-clean`; `--selftest` passed all 14 planted controls. The whole release verdict remains `NOT_READY`.
- `check-after-task.R` on this report: structure passed. `closeout.py check` still exits on the brain repository evidence gate because it sees nine unrelated unfinished brain ledgers; the repo-scoped structure validator passed.
- `slop_check.py` on the after-task, site, and Windows receipts: 0 findings. The checker emitted path warnings for ignored files under the brain source root, but its findings count was zero.

## 6. Tests of the Tests

The `simulate_benchmark()` regression now asserts `nrow(result$results) > 0`; it would fail if the example returned an empty result table as the previous Windows run did. The exact-artifact MI runner checks all four dependency configurations, all three steps per cell, and verifies the archive SHA before and after each step. The R-hub object readback was byte-compared with the frozen archive before the platform checks.

## 7a. Issue Ledger

- Fixed: the `\donttest` benchmark example assumed an optional torch runtime; it now selects the supported GNN-off route.
- Fixed: the old Windows failure caused by the empty benchmark output table; the new assertion guards against recurrence.
- Deferred: deployed-site parity; it requires the authorized merge and deployment path.
- Deferred: final release-validator closeout and any maintainer release decision.

## 8. Consistency Audit

Checked the function default, roxygen example, generated Rd, regression test, exact tarball example run, and rendered candidate reference page together. The package behavior default remains `gnn = TRUE`; the example explicitly uses `FALSE`. The current archive identity, local check logs, current-source defaults receipt, MI matrix manifest, R-hub object identity, Windows logs, site build/crawl receipts, and ledger all refer to the stated commit/hash. Independent reviewers confirmed the platform and local-site claims within their bounded scopes. The current live website remains at an earlier source state.

## 9. What Did Not Go Smoothly

The first R-hub Windows run failed because the example implicitly required torch. The example and test were corrected, then rebuilt into a new artifact. The lane-lease helper could not write its registry outside the task's writable roots. The closeout generator's `new` command targets the brain repository root, outside this task's writable roots. Its combined `check` also runs the brain-wide evidence gate and stops on nine unrelated unfinished brain ledgers; the package report passed the repo-scoped structure validator. The system Git shim also requires the Command Line Tools binary because the Xcode license is not accepted.

## 10. Known Residuals

The fresh site is a local candidate. Its deployed routes and search/sitemap remain unverified until merge and deployment. Pixel-level appearance was not assessed; browser inspection covered rendered text and accessibility-tree content. The MI matrix reuses shared dependency symlinks and is not clean-room dependency installation. It supports termwise fixed-effect pooling, not joint covariance, loading, or variance-component pooling and does not establish broad statistical recovery or interval coverage. R-hub rebuilt the extracted archive source before checking, so its results follow the exact uploaded bytes through extraction and source checkout rather than checking the original compressed file in place. The release remains `NOT_READY` and nothing was submitted to CRAN.

## 11. Team Learning

An exported example should state optional-runtime requirements explicitly. A small benchmark example can use a deterministic runtime-free path while preserving the public function default. Keep package-object readback, extraction, platform rebuild, and platform check as separate links in the artifact provenance chain. Memory receipt: used the pigauto LOAD-FIRST constraints and the exact-artifact release lesson from the memory registry; these kept the local site, deployed site, and frozen archive claims separate. Golden Set: not run because this slice covered docs/example portability and release evidence; no known-mistake-class audit was requested.

## 12. Cross-Product Coverage

Covers ✓: `simulate_benchmark()` documentation example, its generated help, and the tested baseline-only example route; current source defaults audit; optional-backend installed MI workflows on the exact archive; local candidate website routes and discovery; macOS exact-archive check; Windows checks of rebuilt source extracted from the exact uploaded object.

Does NOT cover ✗: current deployed-site behavior; visual pixel appearance; clean-room optional dependency installation; statistical recovery beyond the retained bounded DGPs; joint pooling of covariance, loadings, or variance components; every possible default interaction; CRAN submission, acceptance, or publication; merge or deployment authorization.
