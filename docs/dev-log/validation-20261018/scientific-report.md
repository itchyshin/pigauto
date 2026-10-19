# Comprehensive pigauto validation: implementation and production evidence

Status: Totoro production completed and independently reconciled. It establishes the measured outcomes below for the frozen source, data, settings and environment groups. It does not establish broad guarantees outside these regimes, nor authorize method, default, interface, release or publication changes.

## What is verified

The isolated lane starts at main b2e6054685e6bdedb073ba46e058af2c0a619a18. Package R source and public interfaces are unchanged. The original dirty checkout is preserved.

The manifest represents 151 configurations and 32,240 independent dataset attempts: prediction 16,400; ZI 1,400; composition 1,600; repeated observations 1,500; covariates 1,200; trees 60; posterior inference 4,000; analysis-aware inference 6,000; real data 80. Seeds, settings and estimands are frozen in the executable manifest and data-generating contract.

The new public contract tests pass four tests and 86 expectations, with zero failures, errors or skips. They cover all eight types separately and mixed, observed values, type constraints, simplex constraints, row alignment, deterministic seeds, default lambda handling, prediction and saved-object round trips. Two small-calibration warnings are retained.

The main-source Totoro development suite reports PASS 3554, FAIL 0, WARN 182 and SKIP 9. Elapsed time was 4m04s; maximum resident memory was 1.17 GiB. This supports those exercised tests only. Missing rgbif, unavailable BIEN cache in that frozen baseline, and installed-dependency branches explain skips; no skip is promoted to a pass.

## Production execution and accounting

The user-approved Totoro run used 48 single-thread workers, a 54,000-second total budget, and a 1,800-second per-job limit. The jobs completed in 6:12:02 wall time; GNU time reported average CPU use of 4,753% and maximum resident memory of 1,070,932 KiB. All 32,240 manifest jobs have exactly one retained result: 32,199 ok, one failed, and 40 unavailable. The independent review matched every result job ID, seed and scenario to the frozen job table. No job was silently reused from an earlier campaign.

The 1,122,060 result rows comprise 876,017 ok, 246,040 unavailable and three failed metric rows. These are metric-level counts, distinct from the 32,240 job-envelope statuses. Most unavailable metric rows are expected, including unsupported ordinal Brier scores and point-only floor intervals.

The one failed attempt is `prediction_factorial_018__0095` (n=100, OU, lambda 0.3, trait correlation 0.5, observed-variable MAR, 30% missing). Its GNN-off c2 conformal-coverage, width and interval-score rows were non-finite or absent. The fit itself returned without an adapter error; the failure remains in the denominator and is not interpreted as a package crash. The 40 unavailable attempts are the 20 PanTHERIA and 20 FishBase masks. Their frozen records state that the registered local source data were absent. AVONET and the configured BIEN subset each completed 20 masks.

The original production summariser stopped after the workers completed because it required one installed-package inventory hash across the whole run. The corrected post-hoc evaluator re-read all raw records, checked each record against its nested runtime hash and full job specification, verified a single exact captured `sessionInfo` within every scenario/adapter context, and regenerated the tables without altering raw records. The final evaluator regression fixture covers the matrix form returned by `installed.packages()`.

There are six installed-package inventory hashes across the run, reflecting changes in the host library inventory. Every one of the 151 scenario/adapter contexts has one captured session fingerprint. Only `posterior_001` spans multiple inventory hashes; its six groups have one loaded-session fingerprint. Runtime groups and context tables preserve this variation. The six inventory row counts are 343, 346, 350, 351, 352 and 335 packages. Thus the records establish which inventory and loaded session accompanied each result; this is not evidence that the host library remained globally unchanged.

All 32,384 retained worker RDS files have entries in `receipts/production-raw-rds.sha256`. Raw replicate envelopes remain on Totoro at `/home/snakagaw/hsq_work/pigauto_validation_20261018/source-5e25a197/script/validation_20261018/results/production/`; the compressed full result tables, approval receipts, runtime tables and run log are included locally. The source identity is `67ea5121b8e9e6db4009f562ebf3758384deeff7f16a01685e61f5e20761dcfa`; the job-list hash is `0adafdb698ecce426f75ca09e9c6ba6028829ae98495f5453ba919e5b61d8769`. The production approval and resume receipts preserve the user's direct approval and exact budgets.

## Package artifact

The package builds with its vignettes. The strict R CMD check --no-manual stops at unavailable suggested packages brms and rgbif. A separately labelled check of the same artifact with _R_CHECK_FORCE_SUGGESTS_=false returns Status: OK in 4m12.89s, with maximum resident memory 1,073,452 KiB. Its installed tests report PASS 3272, FAIL 0, WARN 163 and SKIP 87. These counts differ from the development suite because installed check conditions skip more tests.

Checked artifact: pigauto_0.11.0.tar.gz, SHA-256 348829845b658ee77c080723948e92a884943b6a4b2d959ebaea3e79ad89b7b2. This is not an all-dependency clean check, a manual build, a CRAN submission or cross-platform evidence.

## Original bounded pilot

All 18 requested attempts were retained. Execution finished in 9m11.71s, with maximum resident memory 899,972 KiB. Thirteen attempts were ok, two failed and three were unavailable. The frozen summariser then rejected explicit unavailable metrics because of an evaluator schema bug. Corrected reanalysis reports PILOT_ACCOUNTED; raw records are unchanged and their source identity remains intact.

The failed routes exposed test-adapter issues: tree data loading and a BIEN dataframe diagnostic. AVONET was also incorrectly marked unavailable because data(c(...)) loaded no named objects. All three root causes were fixed and covered by actual bundled-data and mask fixtures. A bounded follow-up pilot verifies the fits. PanTHERIA and FishBase remain unavailable because registered raw inputs are absent.

The original posterior n=300 and n=1000 cases converged and completed in 92.5s and 223.7s. Analysis lm, glm and lmer n=100 cases completed in 3.4s, 7.4s and 4.1s. These are route timing observations, not scientific calibration results.

## Follow-up pilot

Seventeen full-setting cases cover all twelve analysis model/size/mechanism combinations, corrected tree/AVONET/BIEN routes, a posterior denominator replay and a MAR prediction fit with probability calibration metrics. All seventeen envelopes are ok, with zero failed result rows. Execution finished in 3m04.37s at maximum resident memory 588,092 KiB. All twelve analysis cells used m=50 and retained 50/50 successful downstream fits. The posterior replay used four chains, 5,000 iterations, 1,000 burn-in, thin 20 and 1,000 retained draws; it converged and retained 50/50 downstream fits. Its x intervals account for all 83 requested masked cells; y had zero requested cells and is explicitly unavailable. Sixty reliability-bin rows were retained. No case has replicated calibration evidence.

## Scientific interpretation

The production summary contains 5,697 scenario/method/estimand rows: 5,630 recovery rows and 23 descriptive rows, plus 44 inferential rows. Of the inferential rows, none is supported; 12 are classified failed and 32 indeterminate under the prespecified bands and Monte Carlo uncertainty. This includes complete-data controls, so some failures arise even without imputation. These are findings for the frozen simulation regimes, not universal claims about all uses of pigauto.

Five of eight `posterior_full` slope regimes are failed and three indeterminate. Under the same independent truth, its reported-SE ratios range from 0.750 to 1.039 and empirical 95% coverage from 86.2% to 95.2%. The clearest failures are under-dispersed uncertainty in regimes `posterior_002` and `posterior_004`, where coverage is 86.2% and 88.6%, respectively, with SE ratios 0.750 and 0.813. Regimes `posterior_005` and `posterior_007` also fail the joint criteria despite coverage point estimates near 95%, because standardized bias or Monte Carlo bounds miss a prespecified band. Four of eight complete-data posterior controls fail, showing that the miss is not solely attributable to imputation. In analysis-aware inference, 23 of 24 rows are indeterminate and one `analysis_008` GLM row is failed; none is supported. Estimates, intervals and classification reasons are in the full summary table.

For the 82 prediction scenarios, each with 200 replicates, the equal-weighted mean of scenario-level GNN-off conformal coverage is 91.2% for continuous trait c1, 91.1% for c2, 92.3% for counts and 92.1% for proportions. Scenario means range from 69.5% to 97.7% across those four trait families. These descriptive averages fall below nominal 95% and are consistent with the stated limitation that phylogenetic dependence and structured missingness can violate exchangeability. They do not justify a distribution-free coverage claim. Discrete prediction, zero-inflated counts, compositions, repeated observations, covariates, trees and available real-data cases are summarized in retained per-scenario tables; their descriptive results do not establish inferential calibration.


Probability-score scope is explicit in the result rows: all 164 ordinal Brier-score keys are unavailable because Brier score is unsupported for ordinal outcomes; 82 of 164 ordinal ECE keys are unavailable. The other unavailable rows are control-arm interval metrics for the point-only floor. These expected unavailable values are retained rather than interpreted as failures. The production paired-contrast table contains 3,395 scenario/metric contrasts, including 1,655 with 200 independent pairs, 564 with 100 pairs and 28 with 500 pairs. It supports regime-specific paired estimates and Monte Carlo uncertainty, not one global superiority claim.

The runner-contract roster independently expected 5,697 scenario/method/estimand groups, representing 1,122,060 job-level keys. All were present, with zero absent, partially absent, unexpected or unrostered result rows. Three rows were failed and 246,040 were explicitly unavailable, matching their planned status accounting. The payload roster covered 504 reliability-bin and posterior-interval groups: none was missing or failed, and 16,400 GNN-off ordinal calibration-bin records were explicitly unavailable. Strict payload checks reject malformed calibration ranges, counts, statuses, coverage and interval widths. A fresh read-only derivation on Totoro matches all five BIEN trait names and both registered cache hashes. An independent replay of 4,000 posterior DGP cells matches all 8,000 saved x/y interval denominators; a deliberately corrupted denominator is detected. Only trait names, matched-species count, subset seed and cache hashes were retained locally from the BIEN inputs.

The inferential assessment uses the prespecified coverage, bias, SE-calibration bands, Monte Carlo uncertainty, independently defined truth and complete-data controls. It retains all attempts and convergence outcomes. For regimes without an independently justified inferential target, the report limits interpretation to paired differences and descriptive coverage.

Prediction error and discrete Brier scores are accompanied by confidence reliability bins and expected calibration error where probability evidence is supported. Brier scores alone are not described as proof of calibration. Probability evidence must cover every requested masked cell. Posterior cell interval summaries require their requested-cell denominators; legacy pilot interval metrics without them remain unavailable.

Conformal coverage is empirical under the tested regime. Phylogenetic dependence and structured missingness prevent a blanket distribution-free guarantee. Tree results concern draw-level recovery from the first ten of fifty bundled trees through the default shared route. BIEN concerns a deterministic 300-species subset, not full-dataset scaling.

## Provenance and retained evidence

Original pilot Git tree: a440c5ff42ca54a2680894a6aad22f1ee2e859eb. Source archive SHA-256: 8426b826a5904736fd7c5ca41c6d67fe2fc9398db5e602aba094f9fdd1fb1ace. Raw source identity with configured BIEN input bytes: 56dc8ba44821253a73f701c4fdffb329e70c526c8b906835e63fe06cf205b45d.

Follow-up Git tree: 5b858b7e8892473e0fcb2ea339aa1be94f818cad. Source archive SHA-256: fffa9a0aec4d8d78da954b909ae2fb41d34fb2f0f7903afd9a7ec318a1d18702. Runtime fingerprints now include RNGkind. Source, external input, job and runtime hashes must match before resume or production approval. Later evaluator corrections carry separate evaluator hashes in evaluation_provenance.rds; old raw identities are never rewritten.

Raw envelopes remain on Totoro under /home/snakagaw/hsq_work/pigauto_validation_20261018/. The frozen jobs, summaries, complete result rows, attempt tables, calibration bins, posterior interval tables, paired contrasts and runtime tables are retained as compressed receipts here. The raw RDS inventory and compressed receipt files have separate SHA-256 manifests. Working output directories are excluded from the package and Git history. Small pilot result bundles and their checksums are also retained here.

## Production decision

The user approved the exact Totoro production proposal after the pilot and refreshed host check. `receipts/production-approval.rds` and `receipts/production-approval-resume.rds` bind the approval to the source, initial runtime identity, job-list hash, worker count, time budgets, context and approval text. The original proposal JSON remains the pre-approval proposal artifact; its `approved: false` field is not the production authorization record.

The production estimate was 5.5 to 14.8 elapsed hours. The measured wall time was 6:12:02. Totoro remained under the 150-core ownership cap during the run; the run used 48 single-thread workers. All fits completed across the full 32,240-job manifest. The production shell exited with status one because the original summary step rejected mixed runtime inventory hashes after the fits had finished. The corrected evaluator regenerated and independently reconciled the results. No method, default, public interface or package source was changed.

## Independent review and writing assessment

A separate reviewer checked the retained job and attempt tables, metric counts, analysis m, sampler settings, reported interval denominators, calibration bins, source identities and inferential classifications. A later independent DGP replay verifies the posterior interval denominators against regenerated masks. The paired-contrast table contains 3,395 scenario/metric contrasts, including 1,655 with 200 pairs and 28 with 500 pairs. No inferential scenario met the prespecified support criteria; paired results remain regime-specific and do not support a global superiority claim.

Writing assessment: perceived AI-like style 2/10 (medium confidence), based on the report's specific counts and bounded claims. Evidence: “no inferential scenario met the prespecified support criteria” states the measured result directly; no change needed. “It does not establish broad guarantees outside these regimes” keeps the scope clear; no change needed. Scientific gate: pass, classifications and calibration summaries checked against retained production tables and prespecified criteria. Factual gate: pass, run accounting, source identities, BIEN roster and denominator replay checked against receipts and independent review. Reference gate: not applicable, no external references are cited. Assessment provenance: self-review by Codex, whole scientific report, independent mode, version v1, 2026-10-18; relevant human-control examples were seen before assessment.
