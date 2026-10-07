# Current source package check

- Source commit: `d2c2736d54b3a83dcc626a023ba3f44ae2c71e01`
- Command: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla -e 'devtools::check()'`
- Platform: R 4.6.0, macOS Tahoe 26.7, aarch64-apple-darwin23
- Options: `--no-manual --as-cran`; `NOT_CRAN=true`; `_R_CHECK_FORCE_SUGGESTS_=FALSE`; remote incoming checks disabled by `devtools::check()`.
- Duration: 6m37.4s; check completed on 2026-10-07 at 22:15:53 UTC.
- Result: 0 errors, 0 warnings, 1 NOTE. The NOTE reports that remote system time could not be verified. The package tests and vignette rebuild completed successfully; the embedded test phase took about 234 seconds.
- Limitation: this is current-source validation. Developer-check settings differ from the final frozen-tarball gate, which must run with `NOT_CRAN` unset and force-Suggests enabled.
