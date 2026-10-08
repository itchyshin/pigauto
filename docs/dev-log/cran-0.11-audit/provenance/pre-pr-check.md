# Pre-PR source check, 2026-10-06

The audit-branch build is a pre-PR rehearsal. The exact release candidate is to be frozen from the merged source after the maintainer's merge gate.

The worktree was copied to `/private/tmp/pigauto-cran-011-src/pigauto` with `.git`, `graft`, `_site`, and `.unlazy` excluded. The copy had no symlinks. `R CMD build pigauto` produced `pigauto_0.11.0.tar.gz` (SHA-256 `a9c7ce7b65aec751edd0c8246697ecb51eda9ffe19c874566ad246b4b67c4852`). Its DESCRIPTION says version 0.11.0 and omits drmTMB and gllvmTMB dependencies. The tarball contains no `graft/`, `.unlazy/`, `_site/`, `script/`, `dev/`, `docs/`, `pkgdown/`, or `BACE/` subtree.

On macOS Tahoe with R 4.6.0, `R CMD check --as-cran --no-manual pigauto_0.11.0.tar.gz` finished with **Status: OK**. The retained `pre-pr-check.log` is the check's `00check.log` (SHA-256 `1b3d54d54d7c4cc7341ead4feb3b2ef320b4bc1ea0c1d5cc423bdc6968610346`). It includes successful incoming-feasibility, dependency, installation, documentation, example, test, and vignette checks. The source defaults audit also passed 126 assertions. The four isolated installed-library cells passed their fixture and real adapter checks; those cells installed from the worktree rather than this tarball.

The local pkgdown build passed a crawl of 65 rendered pages and 3,621 local references, with all 34 retired direct pages absent. Desktop and mobile views were inspected. Public deployment will be checked after merge.

The remaining gates are the nine gllvmTMB recovery seeds, which require approval for the measured run above three hours; the maintainer's audit-PR merge; a frozen post-merge tarball and platform checks; the open BirdTree redistribution-rights point; and the separate fresh release panel. Nothing has been submitted to CRAN.

**Superseded note (2026-10-07):** later bounded recovery evidence records all ten gllvmTMB seeds in `../recovery/RESULTS.md`. The old nine-seed count above is retained as a historical pre-PR rehearsal statement, not a current gate. The remaining merge, exact-artifact, rights, and release-panel gates must still be verified against current state.

## Fresh source and site checks after reader fixes, 2026-10-07

`devtools::check()` ran from a clean archive of implementation commit `ce438ff` with the current `vignettes/multiple-imputation.Rmd`. The matching raw output is retained in `devtools-source-check.log` (SHA-256 `87e421580f428b56877df6215e8ed6a23654a8ff5edc30ffad406713572fa7e4`). It completed in 6m03.3s with 0 errors, 0 warnings, and one environment NOTE because the remote system clock could not be verified. The test suite, package vignettes, and vignette rebuild all passed. A separate four-NOTE log found during review came from a different temporary copy containing generated site and log files; it is not the log for this clean-archive run. This source check used the local development configuration; it does not replace the exact post-merge tarball check.

A fresh standard pkgdown build used the current vignette and `_pkgdown.yml` in an isolated source copy. The normal Pages cleanup and crawler then passed with 64 HTML pages, 3,590 local references, all 34 archived benchmark pages absent, and zero structural errors. A planted `VALIDATION_LEDGER.html` and `AGENTS.md` made the updated crawler fail as expected. The generated multiple-imputation article contains the separate optional-package installation instructions. The homepage HTML retains its full experimental warning in the main text and omits the duplicate sidebar warning.

The live deployment still shows the old sidebar and `VALIDATION_LEDGER.html` route until merge and deployment. The live sitemap and local candidate visual review remain open. The CRAN index currently lists 0.10.0, published 2026-07-30, so 0.11.0 remains the unpublished candidate version at this check. No tarball was frozen from these source edits, and no release submission occurred.
