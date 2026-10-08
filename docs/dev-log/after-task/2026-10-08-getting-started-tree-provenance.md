## 1. Goal
Correct the false MCC description of `tree300` in the Getting Started vignette and verify the rendered reader-facing text.

## 2. Implemented
Replaced the MCC claim with the accurate description that `tree300` is posterior member 69 of the BirdTree Hackett-backbone sample, pruned to the 300 AVONET species. Added the Hackett et al. (2008) reference and DOI.

## 3a. Decisions and Rejected Alternatives
Kept the change limited to the demonstrably inaccurate vignette description and its source citation. Did not alter or remove bundled phylogenies because redistribution rights remain unresolved and require a separate release decision. Did not change historical NEWS claims.

## 4. Files Touched
- `vignettes/getting-started.Rmd`
- `docs/dev-log/after-task/2026-10-08-getting-started-tree-provenance.md`

## 5. Checks Run
- `git diff --check`: passed.
- Rendered `vignettes/getting-started.Rmd` to `/private/tmp/getting-started-birdtree.html`: passed; rendered text contains the posterior-member description and the Hackett citation. Pandoc emitted the existing `--mathjax` deprecation warning.

## 6. Tests of the Tests
This is a documentation-only correction, so no software test was added. The rendered HTML was inspected for both the corrected provenance sentence and the cited paper. Rendering confirms the source reaches the reader surface; it does not independently establish the tree's provenance or redistribution rights.

## 7a. Issue Ledger
- Fixed: Getting Started called posterior sample member 69 an MCC phylogeny.
- Open: rights to redistribute BirdTree-derived objects bundled in the package are not established by the citation instructions alone.
- Open: the deployed Getting Started page still needs a fresh build and deployment after the source change is merged.

## 8. Consistency Audit
Checked the current `origin/main` data generator, `R/data.R`, `inst/NOTICE`, and `vignettes/tree-uncertainty.Rmd`. The generator selects `megatrees::get_tree_bird_n100()[[69L]]`; the data help and notice describe a posterior sample rather than an MCC summary; the tree-uncertainty article already distinguishes the posterior trees. The false MCC wording was isolated to Getting Started. The exact composition and redistribution licence for `trees300` remain separate open items in the release audit.

## 9. What Did Not Go Smoothly
The first render tried to write knitr's intermediate Markdown beside the protected worktree source and stopped. Redirecting the intermediate directory and output to temporary storage made the render pass. The closeout generator defaults to the brain repository, so its first draft landed there accidentally; I removed that single generated file and regenerated the report at this pigauto worktree's absolute path. The pigauto after-task checker passed the structural report checks but stopped on five unmet gates in the tracked `.unlazy/imputation-sim` campaign ledger; that separate campaign was not changed.

## 10. Known Residuals
This local branch is not merged or deployed. No CRAN tarball has been frozen or checked by this slice. Bundled-tree redistribution rights, all other release gates, and live-page verification remain unproven.

## 11. Team Learning
The prior site audit approved this wording, with the Hackett citation. I verified the publication identity, journal details, pages, and DOI against PubMed. Memory receipt: the routed LOAD-FIRST manifest and CRAN release-gate protocol shaped the work. I did not write to the brain vault.

Golden Set: not checked; this prose-only provenance correction did not match a software-regression case.

## 12. Cross-Product Coverage
Covers the Getting Started R Markdown source and this local rendered HTML. It does NOT cover the bundled data files' redistribution rights, their generated help or NOTICE wording, the deployed site, the frozen CRAN artifact, defaults, optional-model multiple-imputation workflows, or other trait families and product surfaces.
