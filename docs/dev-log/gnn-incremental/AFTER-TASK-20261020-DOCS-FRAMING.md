## 1. Goal

Reframe the GNN validation section in general terms while keeping its measured scope and limits clear.

## 2. Implemented

Renamed the article “Pigauto validation results and limits” and the GNN section “When does the neural component add predictive value?” The opening now explains the general question before naming the tested repeated-observation covariate setting. Added a practical interpretation: treat this as a candidate setting, compare pigauto with a covariate-aware non-neural model on task-matched held-out data, and report missingness strata separately.

The scenario table, measured results, inferential status, and scope limits remain explicit. An independent Astra review corrected the ablation description and added continuous-response, explicit `gnn = TRUE`, and graph-specificity limits. Corrected the linked report’s polynomial OLS pooled RMSE from 0.876 to 0.878, matching the defined pooled metric, sqrt(mean dataset MSE). No package method, API, default, or simulation changed.

## 3a. Decisions and Rejected Alternatives

- Interpreted the request as broader framing and practical guidance while keeping the measured limits explicit.
- Kept the tested design and comparator results in the section so readers can assess the limits of the guidance.
- No alternatives involving stronger claims or new simulations were used.

## 4. Files Touched

- `vignettes/validation-results.Rmd`
- `_pkgdown.yml`
- `docs/dev-log/gnn-incremental/RESULTS-20261020.md`
- `vignettes/validation-results.Rmd.assessment.json`
- `docs/dev-log/gnn-incremental/AFTER-TASK-20261020-DOCS-FRAMING.md`
- `docs/dev-log/gnn-incremental/AFTER-TASK-20261020-DOCS-FRAMING.md.assessment.json`

## 5. Checks Run

- `pkgdown::build_article("validation-results")`: rendered the article.
- `pkgdown::check_pkgdown()`: `No problems found`.
- Independent Astra review: verified the candidate-use interpretation against the retained report and paired summary; corrected the weight transfer description and stated that graph-specific value was not isolated.
- `python3 ~/shinichi-brain/tools/slop_check.py vignettes/validation-results.Rmd`: zero findings, 1,425 words, zero em dashes.
- `git diff --check`: passed.
- Independent arithmetic check: sqrt(0.771113038657619) rounds to 0.878, matching the pooled polynomial OLS RMSE in the retained summary.
- Exact-draft assessment completion checks: passed for the article and this report.
- `check-after-task.R` passed the report structure check. The full acceptance reverify remains non-green: the GNN ledger has seven met checkboxes without executed approval records, and unrelated imputation-sim and October-validation ledgers also have unresolved checks. I left those campaign records untouched during this documentation follow-up.
- Scoped unlazy status for `gnn-incremental` shows all eight gates marked met from retained evidence; seven runnable checks have no approval record in this clone, so they were not rerun.
- Rendered article inspection: the new heading and practical interpretation appear in the generated HTML.

## 6. Tests of the Tests

No package code or executable tests changed. The rendered article was inspected for the new heading and qualified use guidance; the paired summary table and full report were checked against the result claims.

## 7a. Issue Ledger

Fixed: the section heading described only the study design and did not foreground the user’s question about when the neural component may help. Added a candidate-use interpretation and a task-matched comparison recommendation. No simulation estimates changed.

## 8. Consistency Audit

Checked the article’s scenario values and limits against `RESULTS-20261020.md` and `paired_summary.csv`; recalculated polynomial OLS pooled RMSE as sqrt(0.771113) = 0.878. The section distinguishes within-fit neural-branch contribution from separately fitted non-neural comparisons, and states that graph message passing itself was not isolated. The generated table of contents follows the new heading.

## 9. What Did Not Go Smoothly

The lane-lease helper could not write to its central registry under this filesystem permission profile. Work stayed in an isolated temporary clone, and no leased campaign closeout file was edited. The repo route tool returned no LOAD-FIRST manifest for that temporary clone; the supplied pigauto AGENTS.md instructions remained the project guidance.

## 10. Known Residuals

The evidence still comes from one repeated-observation design with 500 species, about five observations per species, ten covariates, and the registered missingness pattern. The broader follow-up across sample sizes, repeat counts, signal strengths, and missingness patterns has not been run. The article recommends comparison on the intended prediction task; it does not recommend pigauto as generally superior.

## 11. Team Learning

Memory receipt: the supplied pigauto AGENTS.md rules on regime-specific claims and comparator scope shaped the wording. Brain search found no additional relevant note, and route.py returned no manifest for this temporary clone. No brain files were changed.

Golden Set: not in scope for this documentation-only edit.

## 12. Cross-Product Coverage

Covers: the pkgdown validation article section, its generated table of contents, and the article’s local render and structural check.

Does NOT cover: new simulation regimes, other trait types, single-observation prediction, package methods or defaults, uncertainty calibration, downstream inference, or a broad performance recommendation.
