# BirdTree acquisition guidance: after-task report

## 1. Goal

Give users a cited, user-directed path for obtaining a BirdTree phylogeny and reading a selected tree into pigauto.

## 2. Implemented

Expanded `read_tree()` help to explain BirdTree's subset and full-tree downloads, the 2,500-species subset limit, local Newick/NEXUS import, and requested citations. A single-tree file returns `phylo` for `impute()`; a multi-tree file returns `multiPhylo` for `multi_impute_trees()`, whose current scope is experimental prediction-sensitivity checking only. Added a regression test for the multiple-Newick-file return shape. The help states that pigauto does not fetch BirdTree data automatically and that the multi-tree path is not validated for downstream inference.

## 3a. Decisions and Rejected Alternatives

- Kept tree retrieval user-directed. BirdTree's documented subset workflow requests a download and later email/job-status retrieval; the official pages reviewed do not document a stable package-facing API.
- Linked both the BirdTree download and subset pages. Did not add an automatic fetch helper or runtime network access.
- Shinichi instructed that BirdTree trees may be used with citation and credit. This edit records the source's citation request; it does not assign a formal licence to BirdTree data or replace the separate release rights review.

## 4. Files Touched

- `R/read_tree.R`
- `man/read_tree.Rd`
- `tests/testthat/test-read-tree-multiphlo.R`
- `.unlazy/birdtree-reader/GATES.md` (ignored local ledger)
- `.unlazy/birdtree-reader/check_guidance.py` (ignored local checker)
- `.unlazy/birdtree-reader/check-roxygen-copy.sh` (ignored local checker)
- `docs/dev-log/cran-0.11-audit/tree-acquisition-guidance-2026-10-08.md`

## 5. Checks Run

- Unlazy re-verification: **6/6 runnable gates passed; 1/1 manual gate reviewed; 7/7 total met**. `gate-check.mjs --reverify --approve --timeout 600` ran against a temporary copy of the local ledger because the managed worktree blocks writes to its `.unlazy` directory. Its `CWD` remained this worktree, so each oracle checked the current source. The 600-second limit was needed for the complete pkgdown build; the default 120 seconds timed out once. The reverified runnable checks covered the initial 24 help assertions, roxygen regeneration, Rd parsing, existing Newick/NEXUS smoke tests, the fresh site build and cleanup, and the new multiple-tree test. Following independent review, the strengthened source/Rd checker passed 28 assertions and a fresh complete site build passed the route and inference-limit assertions.
- Official [BirdTree downloads](https://birdtree.org/downloads/), [subset instructions](https://birdtree.org/subsets/), and [FAQ](https://birdtree.org/faq/) reviewed. The site requests citation of Jetz et al. (2012) for full or partial tree data, asks users of its web tool to cite BirdTree.org, documents the subset tool's 2,500-species limit and user-requested download workflow, says the tree distributions are Newick, and recommends a reasonable number of draws (more than 100) for full-tree analyses.
- Independent help review approved the wording after the Newick/FAQ addition and flagged the risk of implying a single tree is enough. A later rendered-page review found that the site checker tested keywords separately, not the route pairings or inference limit. I strengthened it to verify the single-tree-to-`impute()` route, the multi-tree-to-`multi_impute_trees()` route, and the downstream-inference warning; the refreshed full build passed.
- Opened the deployed [`read_tree` help](https://itchyshin.github.io/pigauto/reference/read_tree.html) in Chrome on 2026-10-08. It still shows version 0.11.0 and only the generic Newick/NEXUS reader description; the BirdTree guidance is absent, as expected before this PR is merged and deployed.
- Built the complete pkgdown site with pkgdown 2.2.0 in a clean source copy under `/private/tmp`. The offline build disabled CRAN date lookups in that temporary copy and supplied the package's GitHub link to pkgdown, avoiding blocked network requests without changing project configuration. It completed URL, favicon, OpenGraph, article/reference metadata, sitemap, redirects, and search-index checks. The current rendered `read_tree.html` contains the single/multiple-tree guidance, and the deployment cleanup removed four internal pages from files, sitemap, and search while retaining the reader routes.
- The raw build rendered `AGENTS.md`, `CLAUDE.md`, `VALIDATION_LEDGER.md`, and `goodagents.md` as public HTML and included them in `sitemap.xml` and `search.json`. The Pages workflow already runs `pkgdown/clean-internal-pages.R` after building. I ran that same script against the clean site, then verified all four files and routes were absent from the output, sitemap, and search index while `read_tree.html` remained. The cleanup script in the source copy matches the worktree version.
- A second independent review of the current PR diff verified that the official GitHub release metadata names `tree_bird_n100.rda` and lists the SHA-256 recorded in `inst/NOTICE`. The reviewer did not download the archive or compare its tree bytes with pigauto's saved objects. Reachable history records the selection of posterior member 69 and the safety-floor rationale in commit `9566630` (2026-07-10); the current CPU-pinned AVONET test exercises that safety-floor criterion. No reproducible Robinson-Foulds calculation or comparison record was found, so the phrase “tied for the smallest Robinson-Foulds distance” remains unsupported. The same NOTICE file has work on `release/cran-0.11-audit`, so this slice did not edit it.
- Local visual review remains unmet: Chrome refused the generated `file://` page because its browser policy permits only HTTP/HTTPS and blocks serving the file through a workaround. Generated HTML and sitemap/search entries were inspected as text.
- `git diff --check`: passed.
- A new `read_tree()` test writes a Newick file containing two trees and verifies that both are returned as `multiPhylo` and each entry inherits from `phylo`: **3 assertions passed**.
- The guidance checker now covers the single-tree and multi-tree return classes and their separate handoffs: **28 assertions passed**. The generated Rd passes `tools::checkRd()`.
- The managed worktree blocks R from writing generated files in place. A clean copy under `/private/tmp` was used to run roxygen with the current source overlaid; `cmp` confirmed exact agreement with `man/read_tree.Rd`: **`ROXYGEN_RD_MATCH_PASS`**.
- The after-task structure check passed. Six unrelated unmet gates remain in `.unlazy/imputation-sim/gates/leaf-campaign.md`, `leaf-env.md`, `leaf-prerun.md`, `leaf-results.md`, `leaf-runner.md`, and `.unlazy/tree-provenance/GATES.md`; they belong to other slices and remain untouched.
- The BirdTree-guidance source diff also includes one focused multi-tree regression test. No package tarball was built or frozen in this slice.

## 6. Tests of the Tests

The expanded checker first failed on the newly required multi-tree return and handoff statements, as intended. After those statements were added, all 24 initial assertions passed. Independent review then found the rendered-site checker did not bind each tree class to its correct function or require the inference warning. The strengthened checker first failed because it expected a hyphen in the rendered phrase “prediction sensitivity”; after aligning the assertion with the actual help text, the source/Rd checker passed 28 assertions and a fresh complete site build passed all route and limitation checks. The focused Newick regression test passed all three assertions. Direct roxygen output into the managed worktree was blocked by filesystem permissions; running roxygen in a clean temporary copy and comparing the generated Rd with the checked-in page passed exactly. The first site-gate rerun reached the verifier's 120-second timeout while rendering articles. A later run with a 600-second limit completed the whole build and passed cleanup, sitemap, search, and rendered-page assertions. The Newick/NEXUS smoke test ran with `stop_on_failure = TRUE` and emitted a success marker only after completion.

## 7a. Issue Ledger

- Fixed: `read_tree()` now directs bird-tree users to BirdTree's download workflow and states the requested citations.
- Fixed: the help distinguishes single-tree `phylo` from multi-tree `multiPhylo` results and documents the separate API routes.
- Open: the current README and Getting Started article already show generic local-file tree input. This help slice leaves those pages for the separate public-surface inventory and reconciliation.
- Open: remove or qualify the “tied for the smallest Robinson-Foulds distance” claim in `inst/NOTICE`, or add the missing reproducible comparison evidence. The upstream asset checksum is listed in official release metadata, but no downloaded-byte comparison was made. Ownership of `inst/NOTICE` is shared with the release-audit lane.
- Met for the local build: the Pages workflow's cleanup step removes all four internal coordination pages from rendered files, sitemap, and search index. Deployed-route verification remains open.
- Open: rendered HTML was inspected as text but not visually in a browser because Chrome blocks local `file://` pages and the permitted workflow offers no HTTP-serving workaround.
- Open: package-level release rights and exact-artifact gates remain separate and unmet.

## 8. Consistency Audit

Checked the README's six-line local-file workflow, the Getting Started vignette's custom-data example, `read_tree()` source help, generated `read_tree` Rd, `multi_impute_trees()` input contract, and existing tree-uncertainty guidance. A multi-tree file returns a `multiPhylo` list accepted by `multi_impute_trees()`; `impute()` requires one `phylo`. The updated help distinguishes those paths. The full site inventory will determine whether the README and vignette also need the download steps.

## 9. What Did Not Go Smoothly

The first gate run checked generated Rd before running roxygen and correctly caught the stale page. The test gate initially lacked a success-only output token; it was changed to use `stop_on_failure = TRUE` and rerun. The brain closeout generator targets the brain repository rather than this worktree, so its `new` command was not usable for this repo; no file was written by that failed attempt. The local site-output check became the sixth runnable gate after the initial source/help, Rd, and Newick/NEXUS checks and official-source review. A fresh reviewer check then caught the missing multiPhylo return case; the documentation and regression test were extended before the seven-gate re-verification.

## 10. Known Residuals

This slice updates help source and adds a multi-tree return-shape regression test. The deployed website remains unchanged until deployment. The local build and existing post-processing confirmed the four internal coordination pages are removed from deploy output and discovery files. The independent review verified the upstream release metadata for the asset name and SHA-256; it did not verify a downloaded archive against the bundled tree bytes. The Robinson-Foulds superlative still lacks a reproducible calculation. Visual review, post-deployment verification, formal redistribution terms, the full retired-URL audit, optional model-object MI workflows, and the final 0.11 tarball remain open. The live reader page is still the pre-change version until merge and deployment; the release ledger remains `NOT_READY` with no frozen artifact identity. The repo-wide after-task gate check is also still unmet because of the six separate ledgers listed above.

## 11. Team Learning

Memory receipt: loaded the pigauto `LOAD-FIRST` route, `AGENTS.md`, brain index and model-routing map; searched the brain for BirdTree acquisition and the current release decision. The route's emphasis on comparing source and generated help shaped the stale-Rd gate. Golden Set: not run because this was a help-only change with no known code regression class.

## 12. Cross-Product Coverage

Covers `read_tree()` source help, regenerated Rd, BirdTree download guidance, single- and multi-tree Newick return behavior, existing NEXUS smoke coverage, and a fresh pkgdown build with deployment cleanup and sitemap/search inspection. This does NOT cover visual page review, the README, complete Getting Started reader surface, the unsupported Robinson-Foulds superlative, deployed updates, retired-route verification, formal data-rights classification, optional drmTMB/gllvmTMB integration, or frozen-tarball checks.
