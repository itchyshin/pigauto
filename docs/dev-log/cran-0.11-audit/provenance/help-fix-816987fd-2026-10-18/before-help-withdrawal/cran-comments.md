## Update scope

pigauto 0.11.0 updates CRAN 0.10.0. The default imputation route uses the
phylogenetic baseline with the GNN off. Discrete traits estimate their own
phylogenetic signal; the optional safety floor and signal gate default off.
For eligible continuous, single-observation data without covariates,
automatic multiple imputation selects posterior draws with a separation
prior for residual covariance. Other automatic completions remain
prediction diagnostics and are refused by the downstream fitting and
pooling interfaces.

The optional drmTMB and gllvmTMB adapters support termwise fixed-effect
pooling. Neither model package is an installation dependency of pigauto.
Users install the package needed for their chosen model. Package-independent
extraction and arithmetic tests ship with pigauto; real-package integration
checks are retained in a build-excluded harness.

No licence or maintainer change is proposed. Third-party data credits and
phylogeny provenance remain in inst/NOTICE. No direct reverse dependencies
were found in the primary CRAN index across Depends, Imports, LinkingTo,
Suggests or Enhances on 17 October 2026. The package has no compiled source.

## Exact archive checks

The candidate was built from clean merged source
8e275ca738129c8d608f06ffc407b8141c1a2636. The archive is 5,137,648 bytes,
with SHA-256
14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb.
Native checks used this original gzip archive directly, without rebuilding:
R CMD check --as-cran --run-donttest --timings, NOT_CRAN unset, and
_R_CHECK_FORCE_SUGGESTS_=true. All 29 declared dependencies were available.

Windows R-release 4.6.1, Windows R-devel r90686 and Linux R-release 4.6.1
reported Status: OK. The independent local macOS arm64/R 4.6.0 check
completed with zero errors, zero warnings and one timing NOTE:

    compare_methods: user 2.532 s, system 0.992 s, elapsed 18.406 s

The complete native Windows release/devel checks took 645.813 and 622.432
seconds. Tests took 415 and 404 seconds, and vignette rebuilds took 22
seconds in each. Linux took 336.330 seconds overall; its logged CPU/elapsed
stage pairs were 206/206 seconds for tests and 11/11 seconds for vignettes.
The shared Mac check took 3,564.464 seconds overall, with tests logged at
5 minutes CPU/31 minutes elapsed and vignette rebuilds at 17 seconds
CPU/120 seconds elapsed. These are retained observations, not timing
guarantees for CRAN's machines. No run was repeated to replace its NOTE
or select a faster result.

Each native check reported 3,195 passing expectations, zero failures,
161 test warnings and 86 skips. Separate replays on the same checked
installations matched these totals. Their retained warning objects, messages,
source locations, callers and skip identities passed the four-platform audit.
The diagnostics concern small calibration fixtures, constant traits, ignored
covariates on the selected route, ordinal imbalance, few tree draws,
documented fallback and deliberately unmatched species. They retain their
stated limitations; this does not certify general coverage or inferential
validity. The objects come from the separate replays; the native reporter
retains totals.

All six individual vignette render diagnostics also completed; disabled
teaching chunks remained disabled. Installed size was reported as 5.9 MB
on the remote platforms and 7.0 MB locally, as INFO rather than WARNING.

## URLs, spelling and external services

The URL check found one distinct HTTP 403 at
https://doi.org/10.1126/science.1157704, repeated in the Getting Started
source and rendered guide. Primary Crossref metadata resolved with HTTP
200 and matched the cited Hackett et al. 2008 Science article. The DOI is
retained; automated access to that resolver remained denied during this
check. This is not a claim that the access failure has been fixed.

All 404 flagged spelling tokens were reviewed in their source contexts.
No ordinary-word typo was identified among them. DESCRIPTION flags were
technical terms, names or DOI fragments: al, compositional, de, doi, et,
sysbio, syy, Villemereuil and Vovk.

A separate installed-package run without rgbif or terra passed the
applicable failure-path fixtures: 60 expectations, zero failures or
warnings, and eight explicit skips. This verifies the recorded offline
and missing-package behaviour; it does not establish live external-service
availability. The optional-model installation matrix passed all 12 steps
across neither backend, each backend separately, and both backends, using
fresh pigauto libraries and recorded dependency symlinks.
