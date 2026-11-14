# Exact-archive installed-help and novice workflow check

**Verdict: PASS** for the installed help and the bundled novice example, using the exact frozen archive `a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb` (5,137,837 bytes), source commit `b2e6054685e6bdedb073ba46e058af2c0a619a18`, installed in `successor-b2e6054/mi-matrix/neither`.

The final attempt used `Rscript --vanilla` through the owned process-group supervisor, with one-thread caps, `NOT_CRAN` unset, `R_LIBS_USER` set to the exact `neither` library, `R_LIBS_SITE=/nonexistent`, and a 300-second hard cap. It completed in 3.065 seconds with exit 0. The archive SHA and size were unchanged.

Both installed help topics passed: `impute` shows `gnn = FALSE` in usage and describes `gnn = TRUE` as the explicit GNN opt-in; `pigauto_report` distinguishes baseline-only reports from GNN-enabled fits. The novice example filled its two missing cells, kept all observed values and column classes unchanged, and produced a fit with `model_config$gnn == FALSE`. The browser callback was captured only through this R process's `options(browser=...)`; no browser settings were changed.

The example emitted a small-validation warning because the tiny data set has one held-out cell for affected traits. This smoke check does not support a coverage or general-performance claim. The package imports torch and its namespace loaded; the evidence establishes GNN training is disabled by the fit configuration, not that torch is absent or never initialized during package loading.

The first two assertion failures are preserved with their logs and receipts. Attempt 1 used an incorrect overstrike-cleanup expression; attempt 2 did not normalize wrapped help lines. Attempt 3 fixed those check-only issues and passed. No source or archive was edited and no external action was taken.

See `verdict.json` for paths, command, script hashes, exact artifact binding, output receipts, and limitations. The complete PASS log is `attempt3/run.log`; the full command/run metadata is `attempt3/run-receipt.json`.
