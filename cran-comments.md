## Release scope

This is the first CRAN submission of pigauto. Version 0.10.0 provides
phylogenetically informed imputation for mixed species-trait data. The planned
analysis-aware multiple-imputation backend is experimental and deliberately
limited to one incomplete continuous covariate under MAR for Gaussian `lm`,
binomial-logit `glm`, and Gaussian one-random-intercept `lmer`; fixed effects
only. Conformal-width, Brownian/MC-dropout, PMM, and posterior-tree draws are
prediction diagnostics and are unsupported for downstream inference.

A warning-free 6,000-task controls-only campaign at clean source SHA
`430b2c9` demonstrated that the narrow
analysis-aware target is attainable. This does not validate the package
implementation. Submission remains blocked until `multi_impute_analysis()`
passes its package-level validation campaign on a frozen clean SHA.

## Test environments

Final results are pending the release-candidate checks:

- local macOS, R release: pending
- GitHub Actions, Ubuntu, R release: pending
- GitHub Actions, macOS, R release: pending
- win-builder, R-devel: pending

## R CMD check results

Pending. Release requires 0 errors, 0 warnings, and 0 notes on R release and
R-devel before this section is finalized.

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

- namespace-export and Rd-topic counts will be refreshed after the
  `multi_impute_analysis()` implementation and generated documentation land;
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

The final submission will replace all pending check entries above with the
frozen release-candidate evidence.

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
