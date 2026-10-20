## Update scope

pigauto 0.11.0 updates CRAN 0.10.0. The default imputation route uses the
phylogenetic baseline with the GNN off. Discrete traits estimate their own
phylogenetic signal; the optional safety floor and signal gate default off.
For eligible continuous, single-observation data without covariates,
automatic multiple imputation selects posterior draws with a separation
prior for residual covariance. Other automatic completions remain
prediction diagnostics and are refused by downstream fitting and pooling.

The optional drmTMB and gllvmTMB adapters support termwise fixed-effect
pooling. Neither model package is an installation dependency of pigauto.
Users install their chosen model package. Extraction and arithmetic tests
ship with pigauto; real-package integration checks use a build-excluded
harness. No joint covariance, loading or variance-component pooling is claimed.

No licence or maintainer change is proposed. Third-party data credits and
phylogeny provenance remain in inst/NOTICE. The primary CRAN index checked
on 18 October 2026 contained no reverse dependencies across Depends,
Imports, LinkingTo, Suggests or Enhances. The package has no compiled source.

## Exact archive checks

The candidate was built from clean merged source
b2e6054685e6bdedb073ba46e058af2c0a619a18. The archive is 5,137,837 bytes,
with SHA-256
a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb.

The original compressed archive was checked directly, without rebuilding,
on macOS arm64/R 4.6.0, Windows R-release 4.6.1, Windows R-devel and
Linux R-release 4.6.1. All four reported Status: OK, with zero errors,
zero warnings and zero NOTEs. Checks used --as-cran --run-donttest
--timings, NOT_CRAN unset, force-Suggests enabled and all 29 declared
dependencies available. PDF and HTML manuals were included.

Native check times were 452.674 seconds on macOS, 627.141 seconds on
Windows release, 516.487 seconds on Windows devel and 534.689 seconds on
Linux. The respective vignette rebuild stages reported 12/15 seconds
CPU/elapsed, 20 seconds, 17 seconds and 19/19 seconds CPU/elapsed.
Windows single-value timings are reported as logged. These observations
are not timing guarantees for CRAN machines.

Each native test suite recorded 3,195 passing expectations, zero failures,
161 fixture warnings and 86 skips. Separate native-layout replays capture
individual diagnostic objects for consequence review; their times are not
part of native check timing. The fixtures exercise small validation sets,
constant traits, route-specific ignored covariates, ordinal imbalance,
few tree draws, documented fallbacks and unmatched species. Passing these
checks does not establish general inferential validity or nominal coverage.

Installed help and the bundled beginner workflow passed with the GNN off;
observed values and classes were preserved and both missing values filled.
The optional-model installation matrix passed all 12 steps across neither
backend, each backend separately, and both. All six vignette render
diagnostics passed with authored disabled teaching chunks still disabled.
Installed size was INFO, not WARNING, at 7.0 MB locally and 5.9 MB remotely.

## URLs, spelling and external services

The URL check found one distinct HTTP 403 at
https://doi.org/10.1126/science.1157704, repeated in Getting Started source
and rendered documentation. Fresh primary Crossref metadata returned HTTP
200 and matched the cited Hackett et al. 2008 Science article. The DOI is
retained; automated resolver access remained denied during this check.

All 404 flagged spelling tokens were reviewed in 1,886 source contexts.
No ordinary-word typo was identified among those flags. DESCRIPTION flags
were technical terms, names or DOI fragments.

An isolated installed-package check without rgbif or terra passed 60
expectations with zero failures or warnings and eight explicit skips.
This verifies the tested missing-package and offline paths; it does not
establish live external-service availability.
