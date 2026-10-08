# Current source package check

- Source commit: `21a8b348fc532b26f78e5a259b7d77bc81442f03`
- Command: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla -e 'devtools::check()'`
- Platform: R 4.6.0, macOS Tahoe 26.7, aarch64-apple-darwin23
- Options: `--no-manual --as-cran`; `NOT_CRAN=true`; `_R_CHECK_FORCE_SUGGESTS_=FALSE`; remote incoming checks disabled by `devtools::check()`.
- Duration: 9m02s in the unlazy-bound rerun on 2026-10-07. The source-only commit `d2c2736` had an earlier 6m37s run with one remote-clock NOTE; this later check is bound to `21a8b34`.
- Result: 0 errors, 0 warnings, 0 notes. The embedded test phase and vignette rebuild completed successfully.
- Limitation: this is current-source validation. Developer-check settings differ from the final frozen-tarball gate, which must run with `NOT_CRAN` unset and force-Suggests enabled.
