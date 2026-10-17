# Frozen 0.11.0 platform checks

This bounded package-check workflow checks the original `c2e768ad` compressed source archive on Windows R-release, Windows R-devel, and Linux R-release. It does not rebuild the archive. The prior R-hub checks remain evidence for rebuilt source and retain that qualification.

The upstream public R-hub object uses HTTP. Its content is untrusted until the pinned SHA-256, size, inventory, package, and version checks pass. The helper fails before dependency installation or checking if any identity differs. It verifies the archive again immediately before and after checking, including a failed check. Dependency discovery reads a separate DESCRIPTION copy.

The workflow installs every declared dependency, Pandoc, TeX and libtorch, exercises a real torch tensor, and runs `R CMD check --as-cran --run-donttest --timings` on the original absolute archive path. `NOT_CRAN` is unset, force-Suggests is true, and per-stage timing and vignette output are retained. A nonzero exit, unfinished log, or final status other than OK fails the job. A NOTE or WARNING requires explicit review; none is silently waived.

Estimated cold-run time is 45–75 minutes per parallel job, with a 90-minute job timeout. This is package validation on GitHub runners, with no simulation campaign. Only this named audit branch triggers automatic PR runs. Evidence commits outside the workflow paths do not trigger this workflow.

Before trusting results, download all three `frozen-c2e768ad-*` artifacts from the same completed workflow run, run `python script/cran-0.11-platform/verify_results.py RESULT_DIRECTORY`, and inspect raw check/test logs, dependency/session manifests, skip reasons, manuals and stage timings. The verifier proves declared mechanical checks; independent review still adjudicates the logs and scope. Retain the run/job URLs and archive identities beside the release ledger.

`artifact.py selftest` requires a verified candidate at `frozen-artifact/pigauto_0.11.0.tar.gz`. It runs one positive and five negative identity controls without changing that candidate. Missing or altered archives fail closed.

These checks do not close the deployed website, final post-merge freeze, rights/consent, or submission gates. No submission or merge is performed by this workflow.
