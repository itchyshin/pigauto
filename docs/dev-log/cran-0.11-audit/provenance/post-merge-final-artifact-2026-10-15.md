# Exact-artifact refresh, 2026-10-15

## Frozen source and archive

The final candidate was built from clean merged commit `5e46c46b6970bf153a3b4bd00d55c9ace69225b0` in a detached worktree. The built package is `pigauto` 0.11.0. Archive SHA-256 is `4b32f669bb92474f6ec9acb45eb7c7806467ed104bb4338371e4b14369cf230f`; size is 5,142,383 bytes; inventory has 251 entries with SHA-256 `3934bd7d60734627060bd952980db374541890428d42499cae73e62e48a1a2f6`. No forbidden or excluded development paths were found. The freeze receipt and ordered inventory are retained beside this report.

## Exact archive checks

The full local check used the exact archive with `NOT_CRAN` unset and force-Suggests enabled. The offline-capable `R CMD check --run-donttest --timings` completed with `Status: OK`; testthat reported 3,190 passes, 0 failures, 161 warnings, and 86 skips. Its log is `post-merge-final-offline-check-2026-10-15.log`.

The separate `R CMD check --as-cran --run-donttest` attempt stopped during CRAN and Bioconductor incoming-feasibility checks because the host could not resolve their package-index hosts. It did not reach installation or tests. This is an environment-blocked attempt, not an as-cran pass. Its log is `post-merge-final-as-cran-network-block-2026-10-15.log`.

## Optional model-package workflows

The exact archive installed successfully in four isolated-library cells: neither optional model package installed, drmTMB only, gllvmTMB only, and both installed. Fixture and real-object checks passed in every cell. Real adapter runs checked fixed-effect extraction, component labels, saved/reloaded fit lists, independent Gaussian-regression arithmetic for drmTMB, and independent Rubin pooling arithmetic for both backends. The package itself remains independent of both packages. This establishes bounded termwise fixed-effect pooling interoperability; it does not establish joint covariance pooling or general inferential validity. The harness manifest and every cell log are retained with this report.

## User authorization and data credit

Shinichi directly confirmed in this task that the bundled trees may be used with citation and credit. This is recorded as the maintainer's authorization and warranty for this package release. The source notice and user documentation provide attribution. BirdTree's reviewed downloads page requires citation but does not publish separate redistribution terms. This record does not claim an independent legal opinion or a BirdTree-issued licence.

## Remaining release gates

The previous Windows logs are bound to an older candidate and cannot be reused for this hash. The current R-release and R-devel uploads have not completed. Chrome's native file picker could not open while the office Mac was locked. The exact artifact therefore has no checksum-bound Windows result. The offline local check does not close the CRAN-shaped incoming check, and no independent panel has reviewed this exact archive and updated receipt. Do not claim platform-clean or submission-ready. No CRAN submission occurred.
