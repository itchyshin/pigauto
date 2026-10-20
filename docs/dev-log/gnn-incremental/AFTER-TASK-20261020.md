## 1. Goal

Test whether the current gated neural correction improves held-out prediction in repeated-observation covariate scenarios, with a same-fit neural-off comparison, matched non-neural references, reproducible results, and a qualified pkgdown summary.

## 2. Implemented

Added a frozen five-scenario ADEMP contract, manifest, resumable campaign runner, same-fit ablation, matched GLS and regression references, fail-closed accounting, and negative-control checks. Fixed the optional phylogenetic signal gate to exclude held-out validation and test cells, with a direct `phytools::phylosig()` regression oracle. Ran 500 datasets in each scenario and both requested gate policies. Retained all replicate results in a checksummed archive on Totoro and the compact summary bundle in the repository. Wrote the full scientific report and updated the validation-results article after the article path was released by PR #239.

## 3a. Decisions and Rejected Alternatives

The confirmatory policy is the public `impute()` default, `phylo_signal_gate = FALSE`; `TRUE` is a sensitivity arm. Both use `lambda_mode = "estimate"`. The primary comparison is paired within the exact same fit and changes only the neural contribution. GLS and regression models are reference comparisons; they do not define a universal ceiling. We made no public API, default, or method change. We did not run on DRAC or Kohaku: the measured Totoro CPU lane completed the approved campaign within its cap, and this R/torch workload did not require a GPU. We did not deploy or publish the site.

## 4. Files Touched

- `R/fit_pigauto.R`, `R/phylo_signal.R`, `man/fit_pigauto.Rd`
- `tests/testthat/test-phylo-signal-training-mask.R`
- `.unlazy/gnn-incremental/PLAN.md`, `.unlazy/gnn-incremental/GATES.md`
- `docs/dev-log/gnn-incremental/ADEMP.md`, `GOAL.md`, `PILOT-20261019.md`, `RESULTS-20261020.md`, `PLAN-VS-ACTUAL-20261020.md`, `AFTER-TASK-20261020.md`
- `script/gnn_incremental/` runner, manifest, accounting, test, pilot and production outputs, including `results/production-20261019-v1-summary/` and its compact archive
- `vignettes/validation-results.Rmd`

The 278 MB full replicate archive remains on Totoro persistent home storage. Its path and SHA-256 are recorded in the results report. The old untracked article in the other active-validation worktree was not touched.

## 5. Checks Run

- Deterministic checks, same-fit ablation, accounting and production-launcher guards returned their expected PASS markers before production.
- Totoro strict production audit: `RUN_ACCOUNTING_PASS`; 2,500/2,500 datasets, 5,000/5,000 gate fits and 30,000/30,000 method attempts succeeded, with no failures, invalid results, nonconvergence, skips, missing outcomes, duplicates or unpaired records.
- Summary regeneration on Totoro and from the complete downloaded replicate archive: `SUMMARY_REGENERATION_PASS`.
- `devtools::test()` on the rebased branch: 3,655 passed, 0 failed, 8 skipped, 184 warnings. Skips are missing optional BIEN cache or dependency branches and are not counted as passes. Warnings are small-fixture calibration and other expected test diagnostics.
- `pkgdown::build_site(); pkgdown::check_pkgdown()`: rendered all articles including the updated validation article; `No problems found`. Nonfatal notices include deprecated Pandoc math flags and two noninteractive examples.
- `devtools::check(error_on = "warning")` on the rebased branch: `Status: OK`, 0 errors, 0 warnings, 0 notes; 9m 28s on R 4.6.0 macOS arm64.
- Exact-draft assessments: recorded beside each final draft, with `ASSESSMENT_RECORDED` for all four. The scope-specific closeout checker passed. The after-task structure check passed; its full repository acceptance scan also found unchecked leaves in the separate `.unlazy/imputation-sim/` pipeline. Those belong to other validation work and remain untouched here.

## 6. Tests of the Tests

The generator and GLS passed a hand-calculable covariance oracle; a deliberately broken comparator was detected. The ablation test checks that the mean-baseline weight and fitted linear covariate contribution remain fixed and verifies aligned rows. Accounting tests detect incomplete or duplicate manifests, stale fingerprints, missing and unknown method statuses, and invalid pilot classification. The signal-gate regression test checks a direct `phylosig()` oracle and uses extreme held-out values to ensure split-masked cells cannot influence the estimate. All negative controls were exercised and caught their intended faults.

## 7a. Issue Ledger

Fixed: validation/test leakage in optional phylogenetic-signal estimation; GLS comparator row selection; quadratic-basis handling; permissive/missing terminal method statuses and audit fingerprint/task validation. Regression tests now cover these fixes. Preserved: early invalid pilot archives and logs as historical failures. No statistical method or default was silently changed.

## 8. Consistency Audit

The report and pkgdown table use pooled RMSE defined as the square root of mean dataset-level MSE. Paired intervals and resampling use datasets as the independent units. Confirmatory tests are limited to the two prespecified default-policy scenarios; gate-on results are sensitivity evidence. The report identifies that the interactive polynomial basis exactly represents its generated covariate function and that the nonlinear aggregate GLS parity masks a weaker partly-observed-species comparison. The extracted production bundle had no Git metadata, so its SHA-256 and 66 per-file hashes are the source identity. Local strict audit appropriately rejected the differing laptop runtime; authoritative accounting passed on Totoro.

## 9. What Did Not Go Smoothly

Two early pilot issues were retained and diagnosed: the comparator used a quadratic basis incorrectly and a row selector could include the wrong observations. The campaign was corrected and the expanded pilot rerun before production. Early production-throughput estimates exceeded the approved 630-minute envelope; observed later throughput improved, and the final run completed in about 520 minutes. We corrected intermediate scenario-count reporting and waited for PR #239 to release the article path before editing it.

## 10. Known Residuals

The positive result applies to the two stated repeated-observation covariate scenarios, 500 species, ten standard-normal covariates, Gaussian observation errors, this tree process, missingness design, and current public `impute()` default. It does not establish general superiority over non-neural methods, other trait types, single-observation workflows, empirical datasets, uncertainty calibration, downstream inference, or distribution-free coverage. The no-signal control had a small within-fit loss. Eight package tests were skipped for unprovisioned optional data/dependency conditions. The local package check result is recorded above. The article source is local and has not been deployed.

## 11. Team Learning

Memory receipt: applied pigauto operating rules for isolated worktree use, current-source comparison, Totoro caps, split-leakage checks, dataset-level uncertainty, visible failures, and qualified scientific claims. No brain files were edited.

Golden Set: the external Golden Set was not in scope; adversarial controls were included in this campaign's deterministic and audit tests.

## 12. Cross-Product Coverage

Covers: R source and regression tests, reproducible simulation runner and manifest, Totoro production evidence, retained summaries and provenance, scientific report, and local pkgdown article source.

Does NOT cover: deployment or publication, public API/default/method changes, other trait families, single-observation prediction, uncertainty calibration, downstream analyses, external-data validation, or cross-machine independent reproduction of the neural fits.
