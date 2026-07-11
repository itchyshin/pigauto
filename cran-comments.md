## Release scope

This is the first CRAN submission of pigauto. Version 0.10.0 provides
phylogenetically informed imputation for mixed species-trait data. The
analysis-aware multiple-imputation backend is experimental and deliberately
limited to one incomplete continuous covariate under MAR for Gaussian `lm`,
binomial-logit `glm`, and Gaussian one-random-intercept `lmer`; fixed effects
only. Conformal-width, Brownian/MC-dropout, PMM, and posterior-tree draws are
prediction diagnostics and are unsupported for downstream inference.

A warning-free 6,000-task package-level campaign at clean source SHA
`2e3809d` passed all 24 fixed-effect cells, with 93.9%-96.3% coverage,
pooled-SE/empirical-SD ratios of 0.942-1.030, and 100% finite results. This
validates only the documented narrow scope. The final cross-platform checks are
recorded below; the corrected Windows rerun is the last outstanding release
gate.

## Test environments

Release-candidate checks completed so far:

- local macOS, R release: 0 errors under `_R_CHECK_FORCE_SUGGESTS_=false`; the
  only note is the normal `New submission` note and unavailable-Suggests info;
- GitHub Actions, Ubuntu, R release: 0 errors, 0 warnings, 0 notes;
- GitHub Actions, Ubuntu, R-devel: 0 errors, 0 warnings, 0 notes;
- GitHub Actions, macOS, R-release: 0 errors, 0 warnings, 0 notes;
- win-builder, R-devel: corrected rerun submitted from merged main `ddd16c5`;
  result pending.

## R CMD check results

The exact 4.85 MB tarball was checked locally and completed with no errors. The
authoritative GitHub R-CMD-check matrix is green on Ubuntu release, Ubuntu
R-devel, and macOS R-release. The local check used
`_R_CHECK_FORCE_SUGGESTS_=false` because several optional packages are not
installed on this machine; unavailable-Suggests messages are recorded as
information, not package errors.

## Documentation and website audit

The release audit covers the authored source and rendered output for:

- README and pkgdown home page;
- DESCRIPTION, LICENSE, NEWS, and these CRAN comments;
- all exported-function roxygen and Rd topics, including examples;
- every vignette under vignettes/;
- pkgdown navigation, reference indexing, articles, news, and methodology
  links;
- generated docs/ and inst/doc/ artifacts from a clean source export;
- version numbers, support claims, retired BACE and TabPFN material,
  acknowledgement wording, and broken internal or external links.

Audit inventory at the release-candidate documentation freeze:

- the public namespace and generated Rd topics were refreshed after the
  `multi_impute_analysis()` implementation landed;
- five authored vignettes: getting started, common pitfalls, mixed types,
  tree uncertainty, and GNN architecture;
- no hand-maintained files under `inst/doc/`; vignette artifacts are generated
  from the five source files during the package build;
- all five vignette destinations are indexed in the pkgdown navigation, and
  every namespace export is assigned to a pkgdown reference section;
- a clean pkgdown build produced 91 public HTML files after removing three
  internal agent-coordination pages; a local href/src inventory
  found no missing internal targets;
- current user-facing pages contain no live BACE integration or TabPFN support
  claim. Historical NEWS entries remain as an accurate release record, and the
  README acknowledgement retains TabPFN only to describe contributor work;
- five historical benchmark pages had their missing image assets restored;
  they remain outside the release navbar until rerun against release code;
- BACE-named static pages remain only as historical head-to-head benchmark
  evidence. They are not indexed in the navbar or reference, and do not
  represent a supported BACE integration or package dependency;
- `urlchecker::url_check()` examined 19 package URLs with no errors or
  redirections after canonicalizing the pkgdown URL.

The final submission will add the win-builder result and the exact CRAN
submission timestamp to this record.

## Bundled-data licensing gate

`inst/NOTICE` separates the package's MIT-licensed code from bundled
third-party data. AVONET is attributed under its Creative Commons Attribution
licence. The bundled example phylogenies are pruned derivatives of the
MIT-licensed `megatrees` 1.0.0 release asset; its exact SHA-256 digest and the
underlying Jetz et al. (2012) citation are recorded in the notice. The
previously bundled Delhey-derived object was removed because an explicit data
redistribution licence could not be verified.

## Downstream dependencies

There are currently no known reverse dependencies on CRAN.
