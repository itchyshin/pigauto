# BirdTree provenance correction: after-task report

## 1. Goal

Make the provenance and attribution for pigauto's bundled BirdTree objects accurate and consistent across source help, generated manuals, tree-generation scripts, and the installed notice.

## 2. Implemented

Corrected the tree descriptions to identify the relevant posterior sample, its backbones, the deterministic 50-tree composition, and the example-only use. Added citations and clarified that the `megatrees` software licence does not establish redistribution rights for the underlying BirdTree data.

## 3a. Decisions and Rejected Alternatives

- Retained the bundled examples for this source change while making their provenance and rights status explicit. Removal or replacement requires a separate implementation decision.
- Did not claim that public download availability or citation guidance grants redistribution permission. The BirdTree pages inspected require citation but do not state a data redistribution licence.
- Recommended user-directed download and explicit tree input as the release path while rights remain unresolved. An automatic fetch helper would add network and source-maintenance behaviour without resolving that rights question.

## 4. Files Touched

- `R/data.R`
- `data-raw/make_avonet300.R`
- `data-raw/make_avonet_full.R`
- `data-raw/make_trees300.R`
- `inst/NOTICE`
- `man/tree300.Rd`
- `man/tree_full.Rd`
- `man/trees300.Rd`
- `docs/dev-log/cran-0.11-audit/tree-provenance-2026-10-08.md`

## 5. Checks Run

- A deterministic R check confirmed the seed-42 sample contains 28 Ericson and 22 Hackett trees.
- Fresh roxygen regeneration matched all three generated tree help pages.
- Source, generation comments, and `inst/NOTICE` passed provenance and rights-caveat assertions.
- The Getting Started vignette rendered with the corrected posterior-member description.
- A fresh pkgdown reference build passed text assertions for `tree300`, `tree_full`, and `trees300`.
- `R CMD build --no-build-vignettes --no-manual --no-resave-data` succeeded. The fresh candidate source archive included the three tree manuals and `inst/NOTICE`. Its SHA-256 was `04a916fb317a2ff25150eb4c381a6ee7df34a11676c59c49b8e37e52b67d62f8`, size 5,048,555 bytes; it was a working-tree smoke build, not the frozen release artifact.
- The approved unlazy ledger passed all eight gates, then passed `--reverify` with 8/8 runnable checks rerun successfully. The checks covered sample composition, roxygen output, source/NOTICE agreement, Getting Started rendering, patch scope, freshly rendered reference pages, archive contents, and this report's structure.
- The project-wide `check-after-task.R` passed report structure validation, then returned exit 1 because five unrelated `.unlazy/imputation-sim/gates/leaf-*.md` ledgers remain unmet. Those belong to another active campaign lane and were left unchanged.
- The saved report passed `slop_check.py` with 0 findings. Its pattern scan reported 2 X-not-Y matches per 1,000 words, below the configured threshold.
- `git diff --check` passed before this report was added.
- Independent source review approved the documentation and provenance changes. It did not independently compare saved trees with upstream files or decide redistribution rights.

## 6. Tests of the Tests

The first fresh roxygen comparison failed because `tree_full.Rd` had not yet received the updated rights caveat. After updating it, the same exact source-to-generated-page comparison passed. This demonstrates that the check detected a stale generated manual.

## 7a. Issue Ledger

- Fixed: tree descriptions and installed attribution did not consistently identify the selected BirdTree posterior sample and mixed-backbone composition.
- Fixed: `tree_full` lacked the explicit underlying-data rights caveat present in the revised source wording.
- Open: a written permission or licence basis for redistributing the underlying BirdTree data has not been identified.
- Open: local visual review, deployment verification, and the exact frozen-artifact release gate remain outside this slice.
- Open: the repository-wide after-task closeout cannot pass while five unrelated imputation-simulation ledgers are unmet; this lane's scoped eight-gate ledger passed.

## 8. Consistency Audit

Checked the three tree help topics, their generated Rd files, all three tree-generation scripts, and the installed notice. Checked the rendered help pages and Getting Started output. The account now distinguishes the MIT licence for `megatrees` software from rights in BirdTree data. The official BirdTree download page provides archives and a species-subsetting tool and requires citation of Jetz et al. (2012) and BirdTree.org; the page inspected does not state redistribution terms ([BirdTree downloads](https://birdtree.org/downloads/), checked 2026-10-07 in America/Edmonton).

## 9. What Did Not Go Smoothly

The first regeneration check exposed stale generated help for `tree_full`; it was corrected and the check passed on rerun. The pkgdown build reported deprecated MathML options and could not write a Sass cache in the protected user cache, but it completed and produced the requested reference pages. The brain-specific closeout helper resolves paths against the brain repository, so it was not used to create this repository's report.

## 10. Known Residuals

This slice does not establish redistribution rights for BirdTree data. The candidate package build is not a frozen CRAN artifact. The reference build received structural text review but not browser visual review. The repository-wide after-task validator remains red because of five unrelated active imputation-simulation ledgers. No site deployment was performed or verified. No merge or CRAN submission occurred.

## 11. Team Learning

An attribution requirement is not itself a redistribution licence. Keep the software licence and the underlying dataset's rights separate in package notices and release evidence. Source and generated manuals need a fresh regeneration comparison because a correct roxygen source can coexist with stale checked-in help.

Memory receipt: loaded the pigauto LOAD-FIRST routing manifest and release-audit state. The instructions to diff against main, preserve `r_cal = 0`, and separate source evidence from exact-artifact evidence shaped this bounded documentation and package-build check. Golden Set: not run; this slice concerns data provenance wording and generated help rather than a known package-code regression class.

## 12. Cross-Product Coverage

- BirdTree provenance: covers source help, generated manuals, data-generation comments, and the installed notice. Does NOT cover upstream tree-byte verification, written redistribution permission, every public site route, deployment, or the frozen release artifact.
- User tree acquisition: covers the recommendation to direct users to BirdTree's download and subset tools and then accept an explicit tree object. Does NOT cover an automated download helper or the legal terms for republishing BirdTree data.

### Editorial assessment

Style: 2/10, medium confidence; this is a short evidence report with concrete checks and bounded claims.

Genre and coverage: after-task report, full text reviewed on 2026-10-08.

Evidence and repair: “working-tree smoke build, not the frozen release artifact” preserves the evidence boundary; no wording change needed. The rights statement links to the official download page.

Gates: science = not applicable to this documentation slice; facts = unknown because the saved tree bytes were not independently compared with the upstream asset; references = pass for the official BirdTree download and citation page cited here.

Provenance: self-review of this report; report version 2026-10-07; no style controls or held-out examples were reviewed.
