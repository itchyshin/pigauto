# BirdTree acquisition guidance: after-task report

## 1. Goal

Give users a cited, user-directed path for obtaining a BirdTree phylogeny and reading a selected tree into pigauto.

## 2. Implemented

Expanded `read_tree()` help to explain BirdTree's subset and full-tree downloads, the 2,500-species subset limit, local Newick/NEXUS import, and passing the resulting object to `impute()`. The help states that pigauto does not fetch BirdTree data automatically, records the requested citations, and cites BirdTree's FAQ for the Newick format and full-tree draw guidance. It clarifies that `read_tree()` reads one file and bounds pigauto's experimental multi-tree path to descriptive prediction-sensitivity checks; those stochastic datasets are not validated for downstream inference.

## 3a. Decisions and Rejected Alternatives

- Kept tree retrieval user-directed. BirdTree's documented subset workflow requests a download and later email/job-status retrieval; the official pages reviewed do not document a stable package-facing API.
- Linked both the BirdTree download and subset pages. Did not add an automatic fetch helper or runtime network access.
- Shinichi instructed that BirdTree trees may be used with citation and credit. This edit records the source's citation request; it does not assign a formal licence to BirdTree data or replace the separate release rights review.

## 4. Files Touched

- `R/read_tree.R`
- `man/read_tree.Rd`
- `.unlazy/birdtree-reader/GATES.md` (ignored local ledger)
- `.unlazy/birdtree-reader/check_guidance.py` (ignored local checker)
- `docs/dev-log/cran-0.11-audit/tree-acquisition-guidance-2026-10-08.md`

## 5. Checks Run

- `node ~/.codex/skills/unlazy/scripts/gate-check.mjs --reverify --approve --root <worktree> .unlazy/birdtree-reader/GATES.md`: **5/5 runnable gates passed; 1/1 manual gate reviewed; 6/6 total met**. The runnable checks covered 20 source/help assertions, roxygen regeneration, Rd parsing, existing Newick/NEXUS `read_tree()` smoke tests, and removal of internal pages from fresh site output.
- Official [BirdTree downloads](https://birdtree.org/downloads/), [subset instructions](https://birdtree.org/subsets/), and [FAQ](https://birdtree.org/faq/) reviewed. The site requests citation of Jetz et al. (2012) for full or partial tree data, asks users of its web tool to cite BirdTree.org, documents the subset tool's 2,500-species limit and user-requested download workflow, says the tree distributions are Newick, and recommends a reasonable number of draws (more than 100) for full-tree analyses.
- Independent help review approved the wording after the Newick/FAQ addition and flagged the risk of implying a single tree is enough. The help now distinguishes the one-file reader from BirdTree's tree distributions and states the current scope limit of pigauto's experimental multi-tree workflow.
- Opened the deployed [`read_tree` help](https://itchyshin.github.io/pigauto/reference/read_tree.html) in Chrome on 2026-10-08. It still shows version 0.11.0 and only the generic Newick/NEXUS reader description; the BirdTree guidance is absent, as expected before this PR is merged and deployed.
- Built the complete pkgdown site with pkgdown 2.2.0 in a clean source copy at `/private/tmp/pigauto-pr231-site.AAt8TD`. The build completed, including URL, favicon, OpenGraph, article/reference metadata, sitemap, redirects, and search-index checks. The copy's `R/read_tree.R`, `vignettes/getting-started.Rmd`, `_pkgdown.yml`, and `DESCRIPTION` match this worktree.
- The raw build rendered `AGENTS.md`, `CLAUDE.md`, `VALIDATION_LEDGER.md`, and `goodagents.md` as public HTML and included them in `sitemap.xml` and `search.json`. The Pages workflow already runs `pkgdown/clean-internal-pages.R` after building. I ran that same script against the clean site, then verified all four files and routes were absent from the output, sitemap, and search index while `read_tree.html` remained. The cleanup script in the source copy matches the worktree version.
- A second independent review of the current PR diff verified that the official GitHub release metadata names `tree_bird_n100.rda` and lists the SHA-256 recorded in `inst/NOTICE`. The reviewer did not download the archive or compare its tree bytes with pigauto's saved objects. Reachable history records the selection of posterior member 69 and the safety-floor rationale in commit `9566630` (2026-07-10); the current CPU-pinned AVONET test exercises that safety-floor criterion. No reproducible Robinson-Foulds calculation or comparison record was found, so the phrase “tied for the smallest Robinson-Foulds distance” remains unsupported. The same NOTICE file has work on `release/cran-0.11-audit`, so this slice did not edit it.
- Local visual review remains unmet: Chrome refused the generated `file://` page because its browser policy permits only HTTP/HTTPS and blocks serving the file through a workaround. Generated HTML and sitemap/search entries were inspected as text.
- `git diff --check`: passed.
- The after-task structure check passed. The BirdTree unlazy ledger reverified all six gates: five runnable checks passed and the official-source review remains the completed manual gate. Six other unmet gates exist in `.unlazy/imputation-sim/gates/leaf-campaign.md`, `leaf-env.md`, `leaf-prerun.md`, `leaf-results.md`, `leaf-runner.md`, and `.unlazy/tree-provenance/GATES.md`; they are outside this slice and remain untouched.
- This slice's source diff is limited to the `read_tree()` source help and generated help page. The ignored unlazy ledger and checker also record the official evidence. No package tarball was built or frozen in this slice.

## 6. Tests of the Tests

The updated checker first falsely failed because Rd wrapped the Newick sentence across lines. Normalizing whitespace fixed that. A later assertion caught the source's line-broken “prediction-sensitivity” phrase; the help was reworded plainly as “descriptive checks of prediction sensitivity.” After regeneration, all 20 assertions passed. The Newick/NEXUS smoke test ran with `stop_on_failure = TRUE` and emitted a success marker only after completion. The required report validator then exposed unrelated unmet ledgers elsewhere in the worktree, so its full command did not exit successfully even though this slice's own six gates passed.

## 7a. Issue Ledger

- Fixed: `read_tree()` did not direct bird-tree users to BirdTree's download workflow or state the citation request.
- Open: the current README and Getting Started article already show generic local-file tree input. This help slice leaves those pages for the separate public-surface inventory and reconciliation.
- Open: remove or qualify the “tied for the smallest Robinson-Foulds distance” claim in `inst/NOTICE`, or add the missing reproducible comparison evidence. The upstream asset checksum is listed in official release metadata, but no downloaded-byte comparison was made. Ownership of `inst/NOTICE` is shared with the release-audit lane.
- Met for the local build: the Pages workflow's cleanup step removes all four internal coordination pages from rendered files, sitemap, and search index. Deployed-route verification remains open.
- Open: rendered HTML was inspected as text but not visually in a browser because Chrome blocks local `file://` pages and the permitted workflow offers no HTTP-serving workaround.
- Open: package-level release rights and exact-artifact gates remain separate and unmet.

## 8. Consistency Audit

Checked the README's six-line local-file workflow, the Getting Started vignette's custom-data example, `read_tree()` source help, generated `read_tree` Rd, and existing tree-uncertainty guidance. The new help adds BirdTree-specific acquisition, format, credit, and bounded multi-tree context beside the file-import function. The full site inventory will determine whether the README and vignette also need these instructions.

## 9. What Did Not Go Smoothly

The first gate run checked generated Rd before running roxygen and correctly caught the stale page. The test gate initially lacked a success-only output token; it was changed to use `stop_on_failure = TRUE` and rerun. The brain closeout generator targets the brain repository rather than this worktree, so its `new` command was not usable for this repo; no file was written by that failed attempt. The local site-output check became the fifth runnable gate after the initial four source/help checks and the official-source review were recorded. The slice's six gates now consist of those five runnable checks and the reviewed-source gate.

## 10. Known Residuals

This slice updates help source. The deployed website remains unchanged until deployment. The local build and existing post-processing confirmed the four internal coordination pages are removed from deploy output and discovery files. The independent review verified the upstream release metadata for the asset name and SHA-256; it did not verify a downloaded archive against the bundled tree bytes. The Robinson-Foulds superlative still lacks a reproducible calculation. Visual review, post-deployment verification, formal redistribution terms, the full retired-URL audit, optional model-object MI workflows, and the final 0.11 tarball remain open. The live reader page is still the pre-change version until merge and deployment; the release ledger remains `NOT_READY` with no frozen artifact identity. The repo-wide after-task gate check is also still unmet because of the six separate ledgers listed above.

## 11. Team Learning

Memory receipt: loaded the pigauto `LOAD-FIRST` route, `AGENTS.md`, brain index and model-routing map; searched the brain for BirdTree acquisition and the current release decision. The route's emphasis on comparing source and generated help shaped the stale-Rd gate. Golden Set: not run because this was a help-only change with no known code regression class.

## 12. Cross-Product Coverage

Covers `read_tree()` source help, regenerated Rd, BirdTree download guidance, existing local Newick/NEXUS smoke coverage, and a fresh pkgdown build with the deployment cleanup and sitemap/search inspection. This does NOT cover visual page review, the README, complete Getting Started reader surface, the unsupported Robinson-Foulds superlative, deployed updates, retired-route verification, formal data-rights classification, optional drmTMB/gllvmTMB integration, or frozen-tarball checks.
