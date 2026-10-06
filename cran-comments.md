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

The pre-CRAN audit is in progress. No final candidate tarball has been frozen,
and no platform, incoming, or submission-ready claim is made here. The source,
site deployment, and final installed-artifact identities will be recorded
separately after the audit PR is merged by the maintainer.

Current CRAN package-index metadata (2026-10-06) lists no direct reverse
dependencies in Depends, Imports, LinkingTo, Suggests or Enhances.

## Submission status

Nothing has been submitted for version 0.11.0. Exact-candidate check results,
platform logs, timing, spelling, URL checks and inventory will replace the
pending release-evidence paragraph before any maintainer submission.
