# Comprehensive pigauto validation: implementation and pilot evidence

Status: production scientific assessment pending approval. No production campaign has started. One pilot replicate per case cannot establish reliability, coverage or downstream inferential validity.

## What is verified

The isolated lane starts at main b2e6054685e6bdedb073ba46e058af2c0a619a18. Package R source and public interfaces are unchanged. The original dirty checkout is preserved.

The manifest represents 151 configurations and 32,240 independent dataset attempts: prediction 16,400; ZI 1,400; composition 1,600; repeated observations 1,500; covariates 1,200; trees 60; posterior inference 4,000; analysis-aware inference 6,000; real data 80. Seeds, settings and estimands are frozen in the executable manifest and data-generating contract.

The new public contract tests pass four tests and 86 expectations, with zero failures, errors or skips. They cover all eight types separately and mixed, observed values, type constraints, simplex constraints, row alignment, deterministic seeds, default lambda handling, prediction and saved-object round trips. Two small-calibration warnings are retained.

The main-source Totoro development suite reports PASS 3554, FAIL 0, WARN 182 and SKIP 9. Elapsed time was 4m04s; maximum resident memory was 1.17 GiB. This supports those exercised tests only. Missing rgbif, unavailable BIEN cache in that frozen baseline, and installed-dependency branches explain skips; no skip is promoted to a pass.

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

All production calibration assessments remain indeterminate because the replicated campaign has not run. Supported inference requires the prespecified coverage, bias and SE-calibration criteria with Monte Carlo uncertainty, independently defined truth and complete attempt accounting. Complete-data analyses are assessed too.

Prediction error and discrete Brier scores are accompanied by confidence reliability bins and expected calibration error where probability evidence is supported. Brier scores alone are not described as proof of calibration. Probability evidence must cover every requested masked cell. Posterior cell interval summaries require their requested-cell denominators; legacy pilot interval metrics without them remain unavailable.

Conformal coverage is empirical under the tested regime. Phylogenetic dependence and structured missingness prevent a blanket distribution-free guarantee. Tree results concern draw-level recovery from the first ten of fifty bundled trees through the default shared route. BIEN concerns a deterministic 300-species subset, not full-dataset scaling.

## Provenance and retained evidence

Original pilot Git tree: a440c5ff42ca54a2680894a6aad22f1ee2e859eb. Source archive SHA-256: 8426b826a5904736fd7c5ca41c6d67fe2fc9398db5e602aba094f9fdd1fb1ace. Raw source identity with configured BIEN input bytes: 56dc8ba44821253a73f701c4fdffb329e70c526c8b906835e63fe06cf205b45d.

Follow-up Git tree: 5b858b7e8892473e0fcb2ea339aa1be94f818cad. Source archive SHA-256: fffa9a0aec4d8d78da954b909ae2fb41d34fb2f0f7903afd9a7ec318a1d18702. Runtime fingerprints now include RNGkind. Source, external input, job and runtime hashes must match before resume or production approval. Later evaluator corrections carry separate evaluator hashes in evaluation_provenance.rds; old raw identities are never rewritten.

Raw envelopes, jobs, summaries, complete result rows, attempt tables, paired method contrasts and timing tables are under script/validation_20261018/results/. Persistent Totoro source and result directories are under /home/snakagaw/hsq_work/pigauto_validation_20261018/. Compressed exact logs and their raw-byte SHA-256 receipts are versioned in this report's receipts directory. Working output directories are excluded from the package and Git history. Small pilot result bundles and their checksums are retained with the report; raw outputs also remain on both hosts.

## Production decision

The revised measured extrapolation is 265.9 CPU-hours, or 5.5 to 14.8 elapsed hours with 48 single-thread workers, before queue time. There is no sample-size extrapolation. Twenty-five configurations still use timing from another missingness mechanism, so this remains a planning estimate. Seven configurations retain historical failed or unavailable timing attempts; the two raw-data configurations without successful measurements are PanTHERIA and FishBase. All timing attempts and source identities remain visible. The upper estimate uses a twofold timing margin and 75% worker utilisation.

The concrete proposal is all 32,240 attempts on Totoro, up to 48 single-thread workers, a fifteen-hour global stop and a thirty-minute per-job stop. Fresh aggregate CPU use must remain below 150 cores. The current reviewed source identity with the retained BIEN input bytes is 67ea5121b8e9e6db4009f562ebf3758384deeff7f16a01685e61f5e20761dcfa. The unapproved proposal is production-proposal.json; approved=false cannot pass the launch gate.

The approved plan and operating contract require explicit approval for a production run estimated above three hours. No production approval receipt has been created. Fresh Totoro load and aggregate user CPU/memory checks remain required before launch; DRAC remains the fallback. Scientific failures will be reported rather than changing methods silently.

## Independent review and writing assessment

A separate reviewer verified every raw job ID, status count, analysis m, sampler settings, complete interval denominator, calibration bins and separate source identities. Scalar scientific summaries are descriptive or indeterminate; paired comparisons have at most one independent pair. No supported replicated scientific conclusion is claimed.

Post-draft style assessment: perceived AI-like style 2/10, based on direct outcome statements, concrete counts and explicit limits. Scientific gate: pilot scope and unrun production distinguished. Factual gate: checked against retained logs and raw records. Reference gate: no new external factual claim depends on an unverified reference. The repository prose checker passes with zero findings and no em dashes.
