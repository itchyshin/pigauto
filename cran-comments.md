## Update scope

pigauto 0.11.0 updates the published CRAN 0.10.0 release. The default
imputation route now uses the phylogenetic baseline with the GNN off;
discrete traits estimate their own phylogenetic signal, and the optional
safety floor and signal gate default off. For eligible continuous,
single-observation data without covariates, automatic multiple imputation
selects posterior draws with a separation prior for residual covariance.
Other automatic completions remain prediction diagnostics and are refused
by the downstream fitting and pooling interfaces.

The optional drmTMB and gllvmTMB adapters remain available for termwise
fixed-effect pooling. Neither external model package is a pigauto installation
dependency; users install the package needed by their chosen model. Their
real-package integration checks live in a build-excluded harness, with
package-independent extractor and arithmetic fixtures in the shipped tests.

No change to the software licence is intended. Third-party data attribution
and phylogeny provenance remain in inst/NOTICE.

## Release evidence

The audit remains in progress. The original candidate archive, SHA-256
c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2,
passed direct R CMD check --as-cran --run-donttest --timings checks on
macOS arm64/R 4.6.0, Windows R-release 4.6.1, Windows R-devel r90686,
and Linux R-release 4.6.1. NOT_CRAN was unset and force-Suggests enabled.
All four final check logs report Status: OK.

Final Windows release/devel check totals were 536.636 and 629.689 seconds;
the tests took 333 and 407 seconds, and vignette rebuilds took 19 and 21
seconds. Earlier runs and their failures remain in the evidence record.
Each final check reports 3,195 passing test expectations, zero failures,
161 test warnings and 86 skips. These are not warning-free or complete
feature-validation claims.

A local replay with the same CRAN-mode reporter retained all 161 warning
objects. Independent review classified them as sample-size/calibration,
constant-trait/correlation, covariate-path, ordinal-imbalance, few-tree-draw,
signal-gate-fallback and species-mismatch diagnostics in the bounded test
fixtures. The Windows warning objects have not been captured; equal counts
do not establish equal messages.

Installation guidance has since been clarified to distinguish the CRAN and
GitHub versions. The checked archive is predecessor evidence for that
updated source. The final post-merge archive, deployed-site verification,
warning review, timing decision and independent release panel remain open.
These records will be replaced with final-artifact results before submission.

Primary CRAN package-index metadata refreshed on 2026-10-17 still lists
pigauto 0.10.0. No direct reverse dependencies were found in Depends,
Imports, LinkingTo, Suggests or Enhances. The package contains no compiled
source. No change to the licence or maintainer is proposed.

## Submission status

Nothing has been submitted for version 0.11.0. Exact-candidate check results,
platform logs, timing, spelling, URL checks and inventory will replace the
pending release-evidence paragraph before any maintainer submission.
