# Bounded source-site review: baseline-default help correction

Verdict: **PASS for this committed source snapshot and local site build only.** This does not validate a new release archive, merged main, or production deployment.

## Source and command

- Source branch: `docs/cran-011-baseline-help`; PR #238.
- Full source commit: `816987fd160c1d2c5c58b020b07357ff17c077e5`.
- Worktree was clean before and after verification; the source files were not edited.
- Estimate before running: 5–10 minutes; hard cap 15 minutes. The successful builder took roughly 6 minutes wall time and was not killed.
- Executed `script/cran-0.11-site/build-release-site.sh` with `PIGAUTO_SITE_CACHE_SOURCE=/private/tmp/pigauto-final-closure-20261017/site-source/_site`, `/usr/bin/git` first on PATH and `/usr/local/bin` included for Rscript, and `OPENBLAS_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`, `OMP_NUM_THREADS=1`, `JULIA_NUM_THREADS=1`.
- First attempt stopped before building because the restricted PATH omitted `/usr/local/bin/Rscript`. Its complete build log is retained as `site-build-first-attempt.log`. The second, corrected attempt succeeded (exit 0); its complete build and crawl logs are retained as `site-build.log` and `site-crawl.log`.

## Build and page review

The built site is `/private/tmp/pigauto-cran011-site-816987f.DEUFDP/_site`, with `SOURCE_COMMIT.txt` matching the full commit above. pkgdown completed and `pkgdown::check_pkgdown()` reported no problems. The builder's retirement check passed. The crawl reported 62 HTML pages, 3,588 local references checked, 34 retired pages, and zero errors. The build reported 62 HTML pages and 614 search entries; the sitemap has 62 locations. The four retired simulation-study paths were absent from the build and discovery targets.

I reviewed `reference/impute.html`, `reference/pigauto_report.html`, and adjacent `articles/getting-started.html` in a Chrome tab at a 1200-pixel-wide viewport. The impute page states that the default `gnn = FALSE` uses the baseline alone and `gnn = TRUE` trains and calibrates the GNN; the usage signature displays `gnn = FALSE`. The report reference now calls the output an imputation report and describes baseline-only reports for default fits, with GNN comparison/training history when `gnn = TRUE`. The neighboring guide agrees: its one-call example uses the default, its model-training section labels `gnn = TRUE` as the optional GNN example, and its installation/runtime guidance says torch is only needed when opting into GNN. The screenshots showed readable body text and code without overlap; reference titles wrap to multiple lines beside the logo but remain legible. No exported screenshot files are claimed.

The temporary loopback server was started on `127.0.0.1:8765` for these pages and stopped with Ctrl-C; the owned server process exited 0. Its access log is preserved in the stop output, not as a separate file.

## Scope limits

This is source-branch and local rendered-site evidence only. No R CMD check, new tarball, source-CI completion, merge, live-site check, or production deployment was performed or inferred. The earlier 14ee archive is withdrawn and does not receive validation from this build. No scientific behavior changed in this verification.
