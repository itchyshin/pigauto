# CTmax seed-sensitivity validation

## Goal and ownership
Owner: this Codex lane. Base source 033649b93115c6d8191e791f06b1ddf52eb82d19. Scope is a current-source real-data validation harness, not a new imputation model. Package runtime, the frozen CRAN archive and the release owner's files are excluded.

Acceptance: pinned public inputs and source, fixed outer folds, a finite-output smoke, measured pilot time, complete attempted-fit accounting, paired point-prediction summaries and an independent review. Missing mammal inputs and pending Win-builder results remain separate dependencies.

## Design fixed before fitting
Use the experimental CTmax data and tree from BACE commit b4b3cde02c4b171947042100e26b13b9f7ac425d, as linked in issue #260: 2,661 observations and 524 species. Keep the six reported covariates; use the same body_water replacement for missing medium. Save input hashes, fold membership and software identity. Raw inputs and per-row predictions remain outside the package and public Git history.

Five outer folds are assigned to rows within species using split seed 2601109. These are new row holdouts on the same dataset, not an external replication, whole-species validation or a rerun of the original standard-condition folds. Hide CTmax before preprocessing. Covariates and the full tree remain available. The model never receives outer-test truth. Repeated rows of a species can remain in training by design.

For every outer fold, run one baseline and six GNN fits using seeds 51001 through 51006. All fits use the same outer hidden cells. Default 2,000 training epochs are requested with CPU execution and two torch threads per process. Internal validation splits change with fit seed, so this measures end-to-end seed sensitivity rather than neural initialization alone. Preserve actual training history and warnings.

Let e_is = prediction_is - truth_i. Report cell RMSE sqrt(mean(e_is^2)), species-balanced RMSE sqrt(mean_species(mean_within_species(e_is^2))), per-cell seed range max_s(prediction_is)-min_s(prediction_is), and the empirical distribution of those ranges. Repeat point metrics excluding Notophthalmus viridescens. Preserve all failed or missing attempts in the denominator.

Compare two predeclared three-seed averages, 51001:51003 and 51004:51006, against their constituent single fits on the same cells. No best-seed selection is allowed. These are exploratory paired point-prediction comparisons. Do not average interval endpoints or claim validated ensemble coverage. Report individual-fit interval coverage and widths descriptively; correlated rows do not support binomial independence assumptions. No default change or general calibration claim follows from this dataset.

## Compute and stop rules
Totoro is reachable and supplies torch and pigauto. Use an exported exact source rather than its preinstalled pigauto version. First run one epoch on one fold, then a 2,000-epoch pilot on fold 1 / seed 51001. Historical reference: the reporter measured about 12 minutes per fit with three torch threads on an eight-core laptop; Totoro runtime is unmeasured. Pilot estimate: 40 minutes, with a 40-minute process limit. Stop and report if it overruns or yields invalid/empty output.

Only after the pilot, estimate the full 30-GNN-fit run plus five baselines. Use at most six concurrent fits and two torch threads per fit, below Totoro's 150-core cap. A total estimate above three hours needs user approval after pilot results. Check source identity, output validity, memory and remaining shared-machine capacity before launch. No campaign runs on GitHub Actions or the laptop. Resume completed jobs rather than refitting them.

Pilot receipt: fold 1 / seed 51001 completed 2,000 epochs in 250.3042 seconds with finite output. Pilot metadata predates the strict per-result prepared-input and runner digests, so it remains timing evidence and is not admitted to the final result set. The final metadata-enabled smoke passed in 2.6272 seconds. The production estimate is 25–40 minutes at six concurrent processes, with a 40-minute launcher limit. Totoro launcher PID 1747891 started the 35-fit roster. Raw outputs remain under `~/sim-archive/pigauto/seed260-20261108/` on Totoro. No model settings or scoring rules were changed after inspecting the pilot.

Harness controls: complete synthetic 35-fit output summarized correctly. Wrong prepared-input identity, mixed runner identity, corrupt results and missing results were refused with an attempt ledger. A plan-only launcher test reused 35 valid results, then 34 with one missing job, and rejected a wrong source before fitting. The final remote plan-only check showed zero reused jobs and 35 queued. An initial cross-process synthetic test used an R temporary directory that disappeared at process exit; the reusable self-test now executes its subprocesses before R exits.

## Other two tracks
Mammal issues #250 and #251 have no new reporter material. Existing requests remain sufficient; no repeat reminders are needed. No causal diagnosis can be supplied without a reproducer or stored model/draw evidence.

The separate Pigauto chat owns the CRAN candidate. Its latest receipt reports local, Linux R-hub and GitHub Windows release/devel checks passed for the exact archive, while independent Win-builder results and the final release gate remain pending. Do not duplicate that archive or treat this validation as release authority. Two targeted Gmail searches returned no matching new results; that alone is not proof none exist elsewhere.
