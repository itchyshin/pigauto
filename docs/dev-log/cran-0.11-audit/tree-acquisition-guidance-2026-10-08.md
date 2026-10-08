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

- `node ~/.codex/skills/unlazy/scripts/gate-check.mjs --reverify --approve --root <worktree> .unlazy/birdtree-reader/GATES.md`: **4/4 runnable gates passed; 1/1 manual gate reviewed; 5/5 total met**. The runnable checks covered 20 source/help assertions, roxygen regeneration, Rd parsing, and existing Newick/NEXUS `read_tree()` smoke tests.
- Official [BirdTree downloads](https://birdtree.org/downloads/), [subset instructions](https://birdtree.org/subsets/), and [FAQ](https://birdtree.org/faq/) reviewed. The site requests citation of Jetz et al. (2012) for full or partial tree data, asks users of its web tool to cite BirdTree.org, documents the subset tool's 2,500-species limit and user-requested download workflow, says the tree distributions are Newick, and recommends a reasonable number of draws (more than 100) for full-tree analyses.
- Independent help review approved the wording after the Newick/FAQ addition and flagged the risk of implying a single tree is enough. The help now distinguishes the one-file reader from BirdTree's tree distributions and states the current scope limit of pigauto's experimental multi-tree workflow.
- Opened the deployed [`read_tree` help](https://itchyshin.github.io/pigauto/reference/read_tree.html) in Chrome on 2026-10-08. It still shows version 0.11.0 and only the generic Newick/NEXUS reader description; the BirdTree guidance is absent, as expected before this PR is merged and deployed.
- `git diff --check`: passed.
- `Rscript ~/shinichi-brain/tools/check-after-task.R <report>`: structure check passed, then the repository-wide ledger recheck exited 1 because six unmet gates exist outside this slice's owned paths: `.unlazy/imputation-sim/gates/leaf-campaign.md`, `leaf-env.md`, `leaf-prerun.md`, `leaf-results.md`, `leaf-runner.md`, and `.unlazy/tree-provenance/GATES.md`. The BirdTree ledger itself reverified as 5/5 met. Those other ledgers were left untouched.
- The local source diff is limited to the `read_tree()` source help and generated help page. The ignored unlazy ledger and checker also record the official evidence. No package tarball was built or frozen in this slice.

## 6. Tests of the Tests

The updated checker first falsely failed because Rd wrapped the Newick sentence across lines. Normalizing whitespace fixed that. A later assertion caught the source's line-broken “prediction-sensitivity” phrase; the help was reworded plainly as “descriptive checks of prediction sensitivity.” After regeneration, all 20 assertions passed. The Newick/NEXUS smoke test ran with `stop_on_failure = TRUE` and emitted a success marker only after completion. The required report validator then exposed unrelated unmet ledgers elsewhere in the worktree, so its full command did not exit successfully even though this slice's own five gates passed.

## 7a. Issue Ledger

- Fixed: `read_tree()` did not direct bird-tree users to BirdTree's download workflow or state the citation request.
- Open: the current README and Getting Started article already show generic local-file tree input. This help slice leaves those pages for the separate public-surface inventory and reconciliation.
- Open: the source edit has not been built into a fresh pkgdown site. The deployed page was inspected as a baseline only; post-deployment verification remains open.
- Open: package-level release rights and exact-artifact gates remain separate and unmet.

## 8. Consistency Audit

Checked the README's six-line local-file workflow, the Getting Started vignette's custom-data example, `read_tree()` source help, generated `read_tree` Rd, and existing tree-uncertainty guidance. The new help adds BirdTree-specific acquisition, format, credit, and bounded multi-tree context beside the file-import function. The full site inventory will determine whether the README and vignette also need these instructions.

## 9. What Did Not Go Smoothly

The first gate run checked generated Rd before running roxygen and correctly caught the stale page. The test gate initially lacked a success-only output token; it was changed to use `stop_on_failure = TRUE` and rerun. The brain closeout generator targets the brain repository rather than this worktree, so its `new` command was not usable for this repo; no file was written by that failed attempt.

## 10. Known Residuals

This slice changes help source, not the deployed website. It does not decide a formal redistribution licence, audit every public page or retired URL, validate optional model-object MI workflows, or bind evidence to the final 0.11 tarball. The live reader page is still the pre-change version until merge and deployment; the release ledger remains `NOT_READY` with no frozen artifact identity. The repo-wide after-task gate check is also still unmet because of the six separate ledgers listed above.

## 11. Team Learning

Memory receipt: loaded the pigauto `LOAD-FIRST` route, `AGENTS.md`, brain index and model-routing map; searched the brain for BirdTree acquisition and the current release decision. The route's emphasis on comparing source and generated help shaped the stale-Rd gate. Golden Set: not run because this was a help-only change with no known code regression class.

## 12. Cross-Product Coverage

Covers `read_tree()` source help, regenerated Rd, BirdTree download guidance, and existing local Newick/NEXUS smoke coverage. This does NOT cover the README, complete Getting Started reader surface, pkgdown build, deployed pages, retired routes, formal data-rights classification, optional drmTMB/gllvmTMB integration, or frozen-tarball checks.
