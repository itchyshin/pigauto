## 1. Goal

Correct two reader-facing mismatches found in the deployed CRAN 0.11 site: the homepage warning rendered as a literal marker, and GPU advice did not distinguish the opt-in GNN from the default baseline-only route.

## 2. Implemented

- Commit `45df106ba58f5bc97d7a274bc0a22a19d202fa17` changes the README warning to a standard Markdown warning label.
- Commit `a4cd21e197922e4ccd8c6f2e7c38a74858ce6558` says CUDA, MPS, or CPU selection applies when `gnn = TRUE`, and that default `gnn = FALSE` fits the baseline without calling torch.
- The source correction is a candidate on branch `codex/cran-011-reader-surface`; it has not been merged or deployed.

## 3a. Decisions and Rejected Alternatives

- Keep the documented default `gnn = FALSE`; do not change package behavior to match ambiguous wording.
- Keep the GPU availability examples, but scope them to users opting into GNN training.
- Assumption made without asking: this warning and GPU copy is the full newly verified reader-surface defect. If the post-merge browser review finds another mismatch, reopen the documentation gate and repair it before freezing a new artifact.

## 4. Files Touched

- `README.md`
- `vignettes/getting-started.Rmd`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-10-reader-surface-fix.md`

## 5. Checks Run

- Source inspection in `R/impute.R` confirmed the `gnn = FALSE` formal default and baseline-only branch; `R/utils_torch.R` confirms device selection is CUDA, then MPS, then CPU. The graft graph refresh was unavailable because its cache directory returned `EPERM`.
- `git diff --check` passed.
- The report's required 12-section shape passed `check_after_task()`. The full `closeout.py check` correctly remains red because its repo-wide recheck finds five unmet gates in the separate, pre-existing `.unlazy/imputation-sim/` ledger; that arc's `GOAL.md` says it is complete while its last recorded state still leaves G13c open. This reader-surface slice does not own that ledger.
- `python3 ~/shinichi-brain/tools/slop_check.py <absolute getting-started.Rmd path>` reported 0 findings.
- Fresh `script/cran-0.11-site/build-release-site.sh` at exact commit `a4cd21e197922e4ccd8c6f2e7c38a74858ce6558` passed: 6 crawler tests, `SITE_CRAWL_OK`, `pkgdown::check_pkgdown()` with no problems, 62 HTML pages, 613 search entries, and retired-route assertions.
- Chrome reviewed the fresh local homepage and Getting Started page. The homepage shows a rendered bold “Warning:” label; Getting Started says the GPU choice is conditional on `gnn = TRUE` and the default does not call torch.
- Pat and Rose independently reviewed the correction. Both passed local source/render consistency; Rose confirmed the old archive cannot attest to the corrected docs.

## 6. Tests of the Tests

The release-site builder runs six crawler fixture tests and explicit assertions that retired routes are absent from generated output, sitemap, and search. The exact build passed those controls. This slice did not add a new test because the changed content is prose and the site build checks the rendered page and retirement output.

## 7a. Issue Ledger

- RESOLVED IN CANDIDATE SOURCE: homepage `[!WARNING]` marker displayed literally on the live site; standard Markdown warning rendering verified locally.
- RESOLVED IN CANDIDATE SOURCE: Getting Started GPU advice omitted that `gnn = FALSE` is the default; exact-source output and runtime defaults agree.
- OPEN: these commits are not merged or deployed. The live pages still require a post-deployment review.
- OPEN: the frozen tarball SHA `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4` predates both documentation commits and must not be described as containing them.

## 8. Consistency Audit

Checked the live homepage, live Getting Started page, `impute()` formal default, baseline-only execution branch, and `get_device()` selection order. Checked the fresh rendered homepage and Getting Started output, site crawler, sitemap/search retirement checks, and pkgdown diagnostics. The article's separate tree-uncertainty inference boundary remains an explicit validation limitation and was not changed by this documentation fix.

Memory receipt: loaded the pigauto `route.py` manifest. Its prediction-path and `r_cal = 0` guidance shaped the narrow documentation correction; no statistical behavior changed. The cross-repo Golden Set was not run because this prose-only correction did not touch a known-mistake class.

Golden Set: not run; no known-mistake class was in scope for this documentation-only slice.

## 9. What Did Not Go Smoothly

The first site build archived `HEAD` and therefore omitted an uncommitted vignette edit. I committed the edit and reran the builder at the exact candidate commit. The closeout helper resolved its relative report path to the separate brain vault; that single accidental file was removed, and this report is written in the pigauto worktree.

## 10. Known Residuals

- The merged site still shows the previous homepage warning and GPU wording until the candidate source change is merged and Pages deploys it.
- A new tarball must be built from the final merged source, hashed, inventoried, and checked. Existing platform evidence bound to the older hash does not transfer.
- The earlier R-release and R-devel Win-builder logs show five failures each and do not identify the frozen archive hash. Windows remains unproven for the final artifact.
- This documentation review does not close data rights, platform, or independent exact-artifact gates and does not authorize CRAN submission.

## 11. Team Learning

Pkgdown site builders that archive `HEAD` verify only committed source. Commit the candidate prose before the exact-head build, then inspect the generated page in a browser.

## 12. Cross-Product Coverage

- README homepage warning: covered locally in pkgdown output; does NOT cover the live site until a successful deployment is checked.
- Getting Started GPU/device advice: covered in source and exact-head rendered output; does NOT cover future changes to other GPU instructions or validate model performance.
- Package behavior: inspected only for `gnn = FALSE` and device fallback; this prose slice does NOT cover all defaults, all trait types, MI inference validity, or platform behavior.
- Release artifact: no replacement tarball was built in this slice; this evidence does NOT cover the older frozen tarball or any final CRAN check.
