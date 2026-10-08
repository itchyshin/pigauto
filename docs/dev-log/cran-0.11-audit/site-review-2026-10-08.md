# Cache-busted public-site review, 2026-10-08

I opened the public pigauto site in the Codex In-app Browser, using cache-busted URLs. The pages displayed version 0.11.0.

- The homepage shows the baseline-first default workflow, optional GNN installation, and the distinction between completed data, diagnostic predictions, and analysis-aware imputation.
- The deployed tree-uncertainty article describes `trees300` as a sample drawn from the 100-tree BirdTree collection and cites Ericson, Hackett, and Jetz. The deployed page lacks user-run `megatrees` retrieval steps, BirdTree subset-tool instructions, and the BirdTree web-tool citation instruction. The candidate worktree has the retrieval guidance, but it has not been deployed.
- The deployed `tree300` and `tree_full` help pages still say the trees come from the MIT-licensed `megatrees` package. They do not distinguish that software licence from rights to the underlying tree data. The `tree_full` page describes the tree as a sample shared with `tree300` rather than identifying posterior member 69.
- The deployed getting-started article and current vignette source call `tree300` a Hackett MCC phylogeny. `R/data.R` in the candidate worktree identifies the bundled object as posterior member 69, a posterior sample rather than a consensus tree.
- The deployed multiple-imputation article includes optional installation instructions for both `drmTMB` and `gllvmTMB`. Its inverse-Wishart paragraph says the option remains available, but does not clearly limit it to reproducing the earlier analysis as the candidate source now does.
- Direct sitemap navigation returned `net::ERR_BLOCKED_BY_CLIENT`. The request did not reveal the sitemap contents.

These observations confirm G4 and G7 remain open. They describe the current deployment. Candidate-source updates need a later deployment check. Direct URLs opened with the query `audit=20261008`:

- `https://itchyshin.github.io/pigauto/`
- `https://itchyshin.github.io/pigauto/articles/tree-uncertainty.html`
- `https://itchyshin.github.io/pigauto/reference/tree300.html`
- `https://itchyshin.github.io/pigauto/reference/tree_full.html`
- `https://itchyshin.github.io/pigauto/articles/getting-started.html`
- `https://itchyshin.github.io/pigauto/articles/multiple-imputation.html`
- `https://itchyshin.github.io/pigauto/sitemap.xml`

## Follow-up cache-busted browser check

I reopened the tree-uncertainty, `tree300`, `tree_full`, getting-started, and multiple-imputation pages in the Codex In-app Browser with `audit=20261008r2`.

- The deployed tree-uncertainty article now gives the 50-tree mixed-backbone provenance and cites Ericson, Hackett, and Jetz. It still lacks the user-run `megatrees` retrieval steps, the BirdTree subset-tool alternative, and instructions to cite BirdTree when its web tool is used. The candidate worktree has the retrieval guidance, but it is not live.
- The deployed `tree300` and `tree_full` reference pages still connect the tree sample with the `megatrees` MIT licence. The `tree300` page calls the sample Hackett-backbone; the `tree_full` page says it uses the same sample, without identifying posterior member 69.
- The deployed getting-started page still calls `tree300` a Hackett MCC phylogeny. Its remaining workflow details were not part of this targeted tree-page check.
- The deployed multiple-imputation page still says the inverse-Wishart option remains available without the candidate source's historical-reproduction boundary. It presents `drmTMB` and `gllvmTMB` as optional model packages installed separately, matching the dependency decision.
- The live site search returned “No results found” for `bench_continuous`. The former `/dev/bench_continuous.html` route and `/VALIDATION_LEDGER.html` now return pkgdown 404 pages. These checks confirm those three discovery and retirement outcomes only.

These are live-page observations only. They do not verify the local rendered candidate, sitemap contents, or any later deployment.

## Chrome recheck, 2026-10-08

I reopened key public routes in Chrome with a fresh query string. The package badge remains 0.11.0, but this check did not identify the Pages deployment commit.

- `/reference/tree300.html` still says the tree is a posterior Hackett sample distributed by the MIT-licensed `megatrees` package. This conflicts with the candidate help, which identifies member 69 and states that the package licence does not establish rights to the underlying data.
- `/articles/tree-uncertainty.html` describes the 50-tree sample as mixed Ericson/Hackett posterior trees and correctly limits the workflow to descriptive point-prediction sensitivity. It still lacks the user-controlled retrieval steps and BirdTree web-tool citation guidance present in the candidate.
- `/articles/getting-started.html` still calls `tree300` a Hackett MCC phylogeny. The rest of the visible novice workflow separates completed data, prediction diagnostics, and narrow analysis-aware MI.
- `/articles/multiple-imputation.html` still says `residual_prior = "iw"` remains available without the candidate's historical-reproduction-only boundary. It presents `drmTMB` and `gllvmTMB` as optional packages installed separately, and describes the fixed-effect scope of their pooled outputs.

Direct cache-busted Chrome URLs: `https://itchyshin.github.io/pigauto/reference/tree300.html?audit=20261008lead`, `https://itchyshin.github.io/pigauto/articles/tree-uncertainty.html?audit=20261008lead`, `https://itchyshin.github.io/pigauto/articles/getting-started.html?audit=20261008lead`, and `https://itchyshin.github.io/pigauto/articles/multiple-imputation.html?audit=20261008lead`. This Chrome pass records several mismatches on the current deployment; it did not reveal the Pages commit or cover the full site.

<!-- Earlier assessment superseded by the complete assessment below. -->

## Correction to the saved site snapshot, 2026-10-08

A later independent review reran the crawler against `/private/tmp/pigauto-cran011-site-20261008` and found `articles/simulation-study.md` still present. That saved snapshot therefore fails the retired-route check; its earlier 3,548-reference, zero-error description is withdrawn as evidence of complete retirement. The current candidate was rebuilt and cleaned after the NEWS and reader edits. The approved G5/G5b reverify passed; the current crawler reports 62 HTML pages, 3,441 local references, 34 retired routes, and zero errors; search has 613 entries and 57 unique non-empty paths. Rendered-content assertions found no route file or sitemap/search entry for any of the four retired variants. An isolated planted `.md` route now gives the expected crawler failure, retained at `provenance/site-crawler-simulation-retirement-control.txt`.

## Historical simulation article leak, 2026-10-08

Chrome opened `https://itchyshin.github.io/pigauto/articles/simulation-study.html` and displayed the historical benchmark article. The page itself says it is excluded from the 0.11.0 package build and public site. Its source sits under the `.Rbuildignore` rule `^vignettes/articles$`, but pkgdown still renders vignettes from that source tree. The generated article index linked to it, the sitemap listed it, and the search index included it. Pkgdown also generated a nested redirect at `articles/articles/simulation-study.html`.

The cleanup implementation removes the article's HTML and Markdown outputs, pkgdown's nested HTML redirect, and the article-index entry, then filters the four routes from sitemap and search. The crawler checks those direct routes and both discovery indexes. The earlier saved build reported 62 HTML pages, 3,548 local references, 34 retired routes, and 606 search entries, but the independent recheck found `articles/simulation-study.md` still present. That run is not clean evidence for retirement and its zero-error claim is withdrawn. The citation `provenance/site-crawler-simulation-retirement-control.txt` was absent when the reviewer checked; the isolated negative-control output has since been recreated there. Use the current corrected metrics and results in the correction section above. The earlier `pkgdown::check_pkgdown()` result remains a local package-site check only.

The direct Chrome inspection of the public deployment confirms the current leak; it does not confirm a later deployment. Chrome's URL policy rejected the local `file://` preview, and the policy prohibits using a local server or another browser surface to reach the same files. Local candidate visual review therefore remains open. G4/G5 remain partial and G6 remains unmet until the source change is merged, deployed through the authorized process, and checked on the live routes and discovery surfaces.

Naturalness assessment: 2/10, moderate confidence; concise technical audit note. Coverage: the complete report, including the historical article leak and its bounded candidate verification. Science: pass for the described interface and build behaviour; no claim about the benchmark's scientific validity. Facts: pass against the live Chrome page, `.Rbuildignore`, fresh pkgdown output, cleanup, crawler, and pkgdown check; no deployment after the source correction is claimed. References: pass for the direct public route and repository paths described above. Reviewer: Codex self-review, 2026-10-08.

## Cache-busted live-route recheck in Chrome, 2026-10-08

Chrome opened cache-busted versions of the homepage, article index, and historical simulation-study route using `?audit=20261008codexlead`.

- The homepage returned the 0.11.0 badge and current reader guidance for the GNN-off default, explicit `gnn = TRUE`, continuous posterior MI, and analysis-aware MI limits. This checks visible text on the homepage; other retained pages still require review.
- The live article index still lists “Four ways to impute a phylogenetic trait matrix” and links to `/articles/simulation-study.html`.
- The direct route returned the full historical article. Its own text says it is excluded from the 0.11.0 package build and public site. This confirms that the route remains both discoverable and directly served on the live deployment.
- The browser returned `net::ERR_BLOCKED_BY_CLIENT` for `/sitemap.xml`; no conclusion about the live sitemap is drawn from that attempt.

These cache-busted responses confirm that the source branch's local route cleanup has not reached the deployed site. G7 remains unmet pending merge, deployment, and a fresh check of retained and retired routes and discovery files. The PR page showed source/docs PR #231 still Draft, with no reviews and no deployment. This observation does not identify the current Pages deployment commit.

## Candidate-source performance-claim correction, 2026-10-08

The live Getting Started page still described 10,000-tip support as tested on a standard laptop. The candidate source repeated that claim, called the large-tree path dense in all cases, and described sparse Lanczos as future work. Inspection of `R/build_phylo_graph.R` showed that trees above 7,500 tips use optional `RSpectra` Lanczos when available, with dense eigendecomposition as a fallback. It also showed that distance, adjacency, Laplacian, squared-distance, and phylogenetic-correlation matrices remain dense, while `k_eigen = "auto"` only caps the returned feature count. The vignette now describes these behaviors and the resulting memory burden without claiming laptop performance. A fresh local build rendered the corrected wording; the crawler and `pkgdown::check_pkgdown()` passed. The live page still needs the separate post-deployment check.

Assessment for the cache-busted Chrome addition: 2/10, medium confidence; technical audit record, limited to the added section. “Still lists” and “returned the full historical article” give specific browser outcomes; no repair needed. Science = not applicable; facts = pass for the observed pages and direct route; references = not applicable, no external citation added. Codex self-review, 2026-10-08; the earlier route-control evidence was already known, and these cache-busted responses are new.
