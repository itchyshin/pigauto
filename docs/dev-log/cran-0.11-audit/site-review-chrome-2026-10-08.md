# Chrome recheck of deployed tree and import pages, 2026-10-08

## Scope

Read the current public pages in Chrome with cache-busting query strings, then compared the most relevant source and generated manual on GitHub `main`. No local source, deployed site, or other lane was changed.

## Deployment binding

GitHub Actions [pkgdown run #583](https://github.com/itchyshin/pigauto/actions/runs/37716522073) completed successfully after a push to `main` at `0b0f71fee838c6ed51ef832ed819270eeafaf29b` on 2026-10-07 20:09 MDT. The run links to the public site. PR #231 remains draft and its page says the branch has not been deployed. The live pages show version 0.11.0. This binds the latest successful pkgdown workflow observed here to merged `main`; it does not make pending PR #231 changes live.

## Live-page findings

| Page | Chrome observation | Release implication |
|---|---|---|
| [tree300 reference](https://itchyshin.github.io/pigauto/reference/tree300.html?audit=20261008) | Describes a posterior Hackett-backbone sample distributed by the MIT-licensed `megatrees` package. It does not identify member 69. | Stale against the reviewed acquisition provenance and still overstates what the `megatrees` MIT licence establishes about underlying BirdTree data. |
| [trees300 reference](https://itchyshin.github.io/pigauto/reference/trees300.html?audit=20261008c) | Calls all 50 trees a random sample from the BirdTree Hackett-backbone posterior and attributes the MIT licence to `megatrees`. | Stale against the candidate's mixed Ericson/Hackett composition and rights distinction. |
| [tree_full reference](https://itchyshin.github.io/pigauto/reference/tree_full.html?audit=20261008c) | Describes the data as the same posterior Hackett sample used for `tree300` and repeats the MIT-licence attribution. | Same provenance and rights mismatch as `tree300`. |
| [read_tree reference](https://itchyshin.github.io/pigauto/reference/read_tree.html?audit=20261008b) | Documents Newick/NEXUS parsing through `ape::read.tree()` and `ape::read.nexus()` only. | Does not tell BirdTree users where to obtain trees, how to import one or many, or which access-tool citations to include. |
| [Getting Started](https://itchyshin.github.io/pigauto/articles/getting-started.html?audit=20261008b) | Calls `tree300` a matching pruned Hackett MCC phylogeny from BirdTree.org; the generic local-file example has no BirdTree retrieval guidance. | Reader-facing provenance remains stale and acquisition guidance is missing from the main novice route. |

GitHub `main` [R/data.R](https://github.com/itchyshin/pigauto/blob/main/R/data.R) and [man/tree300.Rd](https://github.com/itchyshin/pigauto/blob/main/man/tree300.Rd) still contain the older `tree300` wording. The live page agrees with that committed source. `tree300` was not corrected in this check; an earlier interim statement in chat that it was corrected was a misread and is withdrawn.

The live pages' displayed package version is 0.11.0, while the official [CRAN package record](https://cran.r-project.org/web/packages/pigauto/index.html?audit=20261008) lists 0.10.0, published 2026-07-30. No version change is inferred from the site badge.

## Limits

This check covered five reader pages and the linked `tree300` source/manual, plus the successful pkgdown workflow record. It did not verify every retained or retired URL, search results, sitemap contents, a fresh local candidate build, local visual rendering, data redistribution rights, or the final release tarball. The sitemap remains unverified because Chrome blocked the prior direct sitemap navigation. Keep G4, G6, and G7 open.
