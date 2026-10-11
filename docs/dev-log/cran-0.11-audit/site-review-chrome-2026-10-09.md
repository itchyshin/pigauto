# Live site check in Chrome, 2026-10-09

This is a pre-merge baseline for the deployed pigauto site. Each URL was opened with a cache-busting query in Chrome. These checks establish what the browser served at the time; they do not bind the deployment to a Git commit.

| Page | Current live observation | Candidate source comparison |
|---|---|---|
| [Getting Started](https://itchyshin.github.io/pigauto/articles/getting-started.html?cran-audit=20261009) | Still calls `tree300` a matching pruned Hackett MCC phylogeny. | Current `vignettes/getting-started.Rmd` identifies posterior member 69 and says it is not an MCC or consensus tree. |
| [tree300 reference](https://itchyshin.github.io/pigauto/reference/tree300.html?cran-audit=20261009) | Calls it a posterior Hackett-backbone sample distributed by `megatrees`; its source line cites Li (2026), `megatrees` 1.0.0 and Jetz et al. (2012). | Current `man/tree300.Rd` identifies member 69 and explains that the `megatrees` software licence does not establish rights for BirdTree data. The deployed page is partly corrected but not yet aligned in detail. |
| [Multiple-imputation article](https://itchyshin.github.io/pigauto/articles/multiple-imputation.html?cran-audit=20261009) | Describes the inverse-Wishart prior as remaining available through `posterior_control`. | Current `vignettes/multiple-imputation.Rmd` limits that option to reproducing legacy results. |
| [Article index](https://itchyshin.github.io/pigauto/articles/index.html?cran-audit=20261009) | Lists “Four ways to impute a phylogenetic trait matrix” and links to `articles/simulation-study.html`. | Current candidate build removes the obsolete public route while preserving historical source and evidence in Git. |
| [Retired simulation-study route](https://itchyshin.github.io/pigauto/articles/simulation-study.html?cran-audit=20261009) | Returns the full historical article. The page itself says it is excluded from the 0.11.0 package build, but the deployed route still serves it. | Candidate fresh-site build and retirement manifest require that route to be absent. |

The current `tree300` reference no longer uses the old MCC label; this updates the 2026-10-08 observation for that page. Getting Started remains stale. The article index and directly served retired route prove that the historical page remains publicly discoverable.

Chrome's direct navigation to `https://itchyshin.github.io/pigauto/sitemap.xml?cran-audit=20261009` returned `ERR_BLOCKED_BY_CLIENT`. The web lookup tool also could not access the sitemap. Search-index membership was not checked during this pass. Those controls remain unverified.

**Gate result:** G7 remains unmet. The source PR is still unmerged, deployed reader pages do not fully match candidate source, the retired route remains linked and served, and sitemap/search closure has not been proven. Recheck after Shinichi authorizes merge and deployment.
