# Live documentation recheck, 2026-10-08

## Scope and method

This is a cache-busted read-only check of the deployed reader pages against the committed source at local `origin/main` (`0b0f71fee838c6ed51ef832ed819270eeafaf29b`, merged PR #230). A second reviewer checked the live rendered pages in a browser. Current lane worktrees and their uncommitted changes were not treated as committed or deployed evidence.

## Findings

| Surface | Deployed observation | Committed-source comparison | Status |
|---|---|---|---|
| [Home page](https://itchyshin.github.io/pigauto/?audit=20261008r11) | The cache-busted page badge is 0.11.0. | `origin/main:DESCRIPTION` is 0.11.0; `origin/main:NEWS.md` calls this a local candidate and says no CRAN submission or public release is implied. | Add a clear candidate-status note near the visible version/install guidance so the site badge is not mistaken for the CRAN version. |
| [Getting started](https://itchyshin.github.io/pigauto/articles/getting-started.html?audit=20261008r8) | The article calls `tree300` a “matching pruned Hackett MCC phylogeny from BirdTree.org”. | `origin/main:vignettes/getting-started.Rmd:161-164` makes that claim. `origin/main:R/data.R:32-40` and `origin/main:inst/NOTICE:34-37` identify `tree300` as posterior member 69 from the Hackett backbone, selected as an example tree, not an MCC tree. | Correct the reader-facing description and any neighbouring generated documentation. |
| [Tree sensitivity article](https://itchyshin.github.io/pigauto/articles/tree-uncertainty.html?audit=20261008r7) | The article describes supplied posterior trees and bundled `trees300`; it gives no explicit retrieval steps and does not say pigauto downloads trees. | `origin/main:vignettes/tree-uncertainty.Rmd:27-36` has the same user-supplied-tree scope and identifies the bundled sample's source collection, but gives no user-initiated retrieval example. | Add clear instructions for users who choose to obtain BirdTree samples themselves, with source and citation details. Do not imply pigauto fetches them. |
| [Multiple-imputation article](https://itchyshin.github.io/pigauto/articles/multiple-imputation.html?audit=20261008r9) | The page says `draws_method = "auto"` chooses posterior draws for eligible continuous-trait inputs and otherwise uses conformal draws. | This matches `origin/main:vignettes/multiple-imputation.Rmd:26-36,56-62` and the resolver in `origin/main:R/multi_impute.R:415-428,438-444,908-950`. | No mismatch identified in this slice. |
| [CRAN package page](https://cran.r-project.org/web/packages/pigauto/index.html) | The package listing reports 0.10.0, published 2026-07-30, with an older vignette set. | The website documents the local 0.11.0 candidate, while `NEWS.md` explicitly distinguishes it from a public release. | Keep candidate status visible and check that retained article links do not silently serve stale versions. |

The deployed tree article is linked from the current home page, so its stale content remains discoverable alongside current candidate material. An initial cache-backed page fetch returned an older 0.11.0.9002 homepage and 0.10.0.9000 article; the browser reviewer used cache-busted URLs and saw 0.11.0 on both current pages. The current source candidate also has separate uncommitted tree-retrieval prose; it is not committed or deployed evidence.

## Tree provenance and rights boundary

The team-reviewed source audit found that the BirdTree citation and mixed-backbone descriptions require consistent updates across `R/data.R`, generated help, the data generator, `inst/NOTICE`, and the tree article. At committed `release/cran-0.11-audit` HEAD `6482992`, the generator and NOTICE describe the mixed Ericson/Hackett sample, while the help and generated `man/trees300.Rd` still describe all 50 trees as Hackett-only. The same branch is based on `aa2ade4`, before current `origin/main` `0b0f71f`; its uncommitted edits do not count as resolved. These corrections need reconciliation in the source lane and a fresh generated-help check.

Rights remain an independent open item. BirdTree's download page requires the Jetz et al. citation and a BirdTree citation when the web tool is used. The upstream `megatrees` package reports MIT + file LICENSE, but that package metadata does not itself establish redistribution terms for the underlying BirdTree data in pigauto. No upstream permission request is recorded here. Keep the release rights gate open until the basis for redistributing the bundled derivatives is clear.

## Evidence limits

This check does not validate package behavior, built HTML from the current source, the deployment commit, tree-data redistribution rights, or any CRAN artifact. The exact candidate source and affected pages must be rebuilt and reviewed after the source lane commits its corrections. The website's candidate badge is not evidence that version 0.11.0 is on CRAN.

## References

- BirdTree downloads and citation instructions: <https://birdtree.org/downloads/>
- CRAN listing: <https://cran.r-project.org/web/packages/pigauto/index.html>
- CRAN Repository Policy: <https://stat.ethz.ch/CRAN/web/packages/policies.html>
- `megatrees` source repository and package description: <https://github.com/daijiang/megatrees>
