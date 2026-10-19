## 1. Goal

Implement and execute the approved comprehensive pigauto validation programme, preserve every outcome, assess the prespecified scientific targets, and deliver reproducible reports and receipts.

## 2. Implemented

An isolated source lane, 151-scenario manifest with 32,240 jobs, all-eight-type public regression fixtures, simulation adapters, retained attempt envelopes, uncertainty-aware summaries, paired comparisons and source/runtime/job approval controls. Totoro production completed the full manifest. No package R source, methods, defaults or public interfaces changed.

## 3a. Decisions and Rejected Alternatives

Use Totoro because measured load and memory show safe headroom. User approved the measured 48-worker production proposal. No campaign ran on GitHub Actions. Historical results remain historical until equivalence is demonstrated. Comparator scope is current GNN-off versus mean/mode, with targeted GNN covariate contrasts. No broad competitor campaign.

Assumptions made without asking: real-data BIEN uses a fixed 300-species subset to keep the pilot bounded; unavailable PanTHERIA and FishBase sources produce explicit unavailable attempts. Tree simulation uses the first ten bundled trees and five imputations per tree through the default shared route. These restrict claims and are documented deviations from broader possible readings of the plan.

## 4. Files Touched

- `.gitignore`
- `LOOP/lanes/validation-20261018/GOAL.md`
- `LOOP/lanes/validation-20261018/checkpoint.md`
- `docs/dev-log/validation-20261018/PLAN.md`
- `docs/dev-log/validation-20261018/after-task.md`
- `docs/dev-log/validation-20261018/data-generating-contract.md`
- `docs/dev-log/validation-20261018/plan-versus-actual.md`
- `docs/dev-log/validation-20261018/production-estimate.csv`
- `docs/dev-log/validation-20261018/production-estimate.md`
- `docs/dev-log/validation-20261018/production-proposal.json`
- `docs/dev-log/validation-20261018/receipts/aftertask-check.log.gz`
- `docs/dev-log/validation-20261018/receipts/archive-sha256.txt`
- `docs/dev-log/validation-20261018/receipts/base-suite-totoro.log.gz`
- `docs/dev-log/validation-20261018/receipts/contract-local.log.gz`
- `docs/dev-log/validation-20261018/receipts/jobs.csv.gz`
- `docs/dev-log/validation-20261018/receipts/mar-sentinel.log.gz`
- `docs/dev-log/validation-20261018/receipts/package-check-available.log.gz`
- `docs/dev-log/validation-20261018/receipts/package-check-strict.log.gz`
- `docs/dev-log/validation-20261018/receipts/package-tests.Rout.gz`
- `docs/dev-log/validation-20261018/receipts/pilot-a440c5f-results.tar.gz`
- `docs/dev-log/validation-20261018/receipts/pilot-a440c5f.log.gz`
- `docs/dev-log/validation-20261018/receipts/pilot-attempts.csv.gz`
- `docs/dev-log/validation-20261018/receipts/pilot-calibration_bins_by_replicate.csv.gz`
- `docs/dev-log/validation-20261018/receipts/pilot-paired_method_contrasts.csv.gz`
- `docs/dev-log/validation-20261018/receipts/pilot-posterior_intervals_summary.csv.gz`
- `docs/dev-log/validation-20261018/receipts/pilot-result_rows.csv.gz`
- `docs/dev-log/validation-20261018/receipts/pilot-summary.csv.gz`
- `docs/dev-log/validation-20261018/receipts/proposed-source.rds`
- `docs/dev-log/validation-20261018/receipts/prose-final.log.gz`
- `docs/dev-log/validation-20261018/receipts/raw-log-sha256.txt`
- `docs/dev-log/validation-20261018/receipts/supplement-5b858b7-results.tar.gz`
- `docs/dev-log/validation-20261018/receipts/supplement-5b858b7.log.gz`
- `docs/dev-log/validation-20261018/receipts/supplement-attempts.csv.gz`
- `docs/dev-log/validation-20261018/receipts/supplement-calibration_bins_by_replicate.csv.gz`
- `docs/dev-log/validation-20261018/receipts/supplement-paired_method_contrasts.csv.gz`
- `docs/dev-log/validation-20261018/receipts/supplement-posterior_intervals_summary.csv.gz`
- `docs/dev-log/validation-20261018/receipts/supplement-result_rows.csv.gz`
- `docs/dev-log/validation-20261018/receipts/supplement-summary.csv.gz`
- `docs/dev-log/validation-20261018/receipts/unlazy-final.log.gz`
- `docs/dev-log/validation-20261018/receipts/unlazy-prepilot.log.gz`
- `docs/dev-log/validation-20261018/reuse-audit.md`
- `docs/dev-log/validation-20261018/review-final.md`
- `docs/dev-log/validation-20261018/review-production.md`
- `docs/dev-log/validation-20261018/review-orchestration.md`
- `docs/dev-log/validation-20261018/review-scientific.md`
- `docs/dev-log/validation-20261018/posthoc_summarise.R`
- `docs/dev-log/validation-20261018/test_posthoc_summarise.R`
- `docs/dev-log/validation-20261018/receipts/production-approval.rds`
- `docs/dev-log/validation-20261018/receipts/production-approval-resume.rds`
- `docs/dev-log/validation-20261018/receipts/production-attempts.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-calibration_bins_by_replicate.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-evaluation_provenance.rds`
- `docs/dev-log/validation-20261018/receipts/production-paired_method_contrasts.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-posterior_intervals_by_replicate.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-posterior_intervals_summary.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-production-resume-20261018.log.gz`
- `docs/dev-log/validation-20261018/receipts/production-raw-rds.sha256`
- `docs/dev-log/validation-20261018/receipts/production-result_rows.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-runtime_contexts.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-runtime_groups.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-summary.csv.gz`
- `docs/dev-log/validation-20261018/validation-ledger.md`
- `docs/dev-log/validation-20261018/plan-versus-actual.md`
- `LOOP/lanes/validation-20261018/checkpoint.md`
- `.unlazy/validation-20261018/GATES.md`
- `docs/dev-log/validation-20261018/scenarios.csv`
- `.unlazy/validation-20261018/GATES.md`
- `docs/dev-log/validation-20261018/assessments/scientific-report.assessment.json`
- `docs/dev-log/validation-20261018/assessments/after-task.assessment.json`
- `docs/dev-log/validation-20261018/assessments/plan-versus-actual.assessment.json`
- `docs/dev-log/validation-20261018/assessments/review-production.assessment.json`
- `docs/dev-log/validation-20261018/scientific-report.md`
- `docs/dev-log/validation-20261018/validation-ledger.md`
- `script/validation_20261018/README.md`
- `script/validation_20261018/check_contract.R`
- `script/validation_20261018/drac.sbatch`
- `script/validation_20261018/estimate.R`
- `script/validation_20261018/evaluate.R`
- `script/validation_20261018/freeze.R`
- `script/validation_20261018/manifest.R`
- `script/validation_20261018/orchestration.R`
- `script/validation_20261018/package_check.sh`
- `script/validation_20261018/pilot.sh`
- `script/validation_20261018/production.sh`
- `script/validation_20261018/run.R`
- `script/validation_20261018/runners.R`
- `script/validation_20261018/summarise.R`
- `script/validation_20261018/test_approval.R`
- `script/validation_20261018/test_harness.R`
- `script/validation_20261018/test_orchestration.R`
- `script/validation_20261018/test_runners.R`
- `script/validation_20261018/test_summary.R`
- `script/validation_20261018/test_timeouts.sh`
- `script/validation_20261018/roster.R`
- `script/validation_20261018/test_roster.R`
- `script/validation_20261018/derive_bien_roster.R`
- `script/validation_20261018/audit_posterior_denominators.R`
- `script/validation_20261018/verify_bien_roster_totoro.sh`
- `script/validation_20261018/audit_roster_receipts.R`
- `docs/dev-log/validation-20261018/receipts/production-bien-source-roster.csv`
- `docs/dev-log/validation-20261018/receipts/production-bien-roster-verification.txt`
- `docs/dev-log/validation-20261018/receipts/production-posterior-denominator-audit.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-receipt-sha256.txt`
- `docs/dev-log/validation-20261018/receipts/production-roster-audit.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-metric-roster-audit.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-roster-summary.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-unexpected-rows.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-unrostered-rows.csv.gz`
- `docs/dev-log/validation-20261018/receipts/production-roster-issues.csv.gz`
- `tests/testthat/test-validation-contract.R`

## 5. Checks Run

The main-source Totoro baseline suite passed 3,554 expectations, with zero failures, 182 warnings and nine skips, in 4m04s at peak 1.17 GiB. New public contract fixtures passed four tests and 86 expectations with zero failures/errors/skips and two small-calibration warnings. Harness, adapter, approval, source, summary and expected-roster controls pass their bounded checks. The local roster test includes malformed probability-bin, interval, status and denominator controls, and passes. All 82 prediction scenarios construct without fitting, including MAR. Linux process-group timeout regression passes on Totoro. Production completed 32,240 attempts in 6:12:02 and retained 1,122,060 metric rows: 876,017 ok, 246,040 unavailable and three failed. The expected roster contains 5,697 groups and 1,122,060 planned job-level output keys, all present; zero were absent, partial, unexpected or unrostered. A fresh Totoro BIEN comparison matches the retained roster. Independent replay of 4,000 posterior DGP cells matches all 8,000 interval denominators, and the deliberate corruption is detected.

Original pilot: 18/18 accounted, thirteen ok, two failed and three unavailable in 9m12s. Corrected follow-up: 17/17 ok in 3m04s, full analysis and sampler settings independently verified. The production estimate was 5.5 to 14.8 hours at 48 workers; the run took 6:12:02. The original summariser failed after compute because it rejected mixed runtime inventory hashes. The corrected evaluator regenerated the summaries; its matrix row-count defect was fixed and retested. The complete roster, Totoro BIEN and posterior-denominator checks passed through the exact approved unlazy G9 gate. The repository-wide after-task acceptance recheck also reports unmet inherited imputation-sim ledgers owned by another lane; those files were preserved. A strict package check must not be called clean when suggested packages are unavailable.

## 6. Tests of the Tests

Negative controls reject duplicate, stale, malformed and missing records; false 100% coverage; altered source, data, job settings and approval budgets; absent DRAC allocations; orphan worker processes; lost prediction cells; mismatched RNG kind; and unsupported metrics incorrectly represented as scientific failures. Known-value fixtures verify Monte Carlo calculations, interval widths, complete-data paired differences and method pairing.

## 7a. Issue Ledger

Fixed harness issues: metric applicability, envelope schema consistency, MAR missingness driver, non-finite score denominators, posterior convergence and downstream attrition, runtime RNG identity, BIEN input binding, source/resume fingerprints, process-group shutdown, mixed runtime inventory evaluation and matrix-backed installed-package counts. The expected-key roster detects globally absent and partially absent outputs, rejects unrecognized result statuses, and validates probability-bin and posterior-interval payloads. Negative controls cover malformed payloads, missing bins, incomplete denominators and a corrupted expected denominator. One prediction replicate retains three failed c2 interval metrics. Package strict-check dependency absence remains explicit. No new package implementation bug has been established by these fixes.

## 8. Consistency Audit

Reviewed adjacent generators, scoring, manifests, timeout launchers, approval checks, summaries, expected-key roster and evidence labels. BIEN's five expected traits were derived on Totoro from the exact source caches; their names and hashes are compared byte-for-byte with the retained verification CSV and approved production input identity. Posterior interval denominators are compared with independently regenerated frozen masks. Current source and raw pilot source identities are separate where evaluator corrections followed the freeze. Existing active-observation brute-force checks and shared/per-tree tests are part of the baseline suite. No unvalidated optional integration is promoted by a skip.

## 9. What Did Not Go Smoothly

The first pilot was stopped because unsupported metrics were mislabelled failed. Its source and files remain superseded harness evidence. The shutdown audit found a child-process-group escape, fixed and tested on Linux. Torch needed an existing private native runtime; vignette building needed an existing private Pandoc path. No shared machine configuration or library was changed. The report scaffold helper resolved a relative path to the vault; the newly generated file was immediately moved to this repository and removed from that unintended location.

## 10. Known Residuals

Production is complete. No inferential row met the support criteria: 12 failed and 32 indeterminate. Four complete-data posterior controls also fail. The 40 unavailable real-data jobs are PanTHERIA and FishBase because inputs were absent. Six installed-package inventory identities occurred; each of the 151 scenario/adapter contexts has one captured session fingerprint, with `posterior_001` spanning six inventories. The independent runner-contract roster covers all 1,122,060 planned job-level method/estimand keys, with zero absent, partial, unexpected or unrostered rows; strict metric-payload checks found zero missing or failed groups. The fresh BIEN roster check and all 8,000 posterior mask-denominator comparisons pass. Independent review confirmed invalid statuses are rejected and Totoro proof content matches the retained roster. Optional brms and rgbif are absent on this host. DRAC was not used.

## 11. Team Learning

Memory receipt: loaded the pigauto route and LOAD-FIRST guidance, prior prediction-first audit, compute routing, uncertainty and fallback invariants, and exact-evidence guidance. These shaped the source-isolated lane, narrow comparator roster, retained failures and Totoro-first pilot. No vault or memory update is shipped without approval.

Golden Set: no package source bug was changed; test-harness negative controls and existing baseline regression tests were exercised. No claim that the full cross-project Golden Set passed. Independent review confirmed exact job accounting, observed-key result/summary parity, scientific classifications, malformed-payload rejection, status validation, the Totoro account and proof-to-roster binding. It identified the installed-package matrix count defect, which was fixed and retested. The global checker also lists inherited imputation-sim ledgers, which this lane does not own. Those ledgers were preserved.

Independent reviewers challenged scientific estimands, failure denominators, scenario constructors, source identity, RNG identity and process ownership. Real run accounting caught issues that syntax checks did not.

## 12. Cross-Product Coverage

Covers: all eight public trait contracts; available baseline package tests; all 151 planned configurations and 32,240 production attempts; explicit failure and unavailable accounting; Totoro production and isolated package-check evidence.

This does NOT cover inferential support for the tested regimes, an all-dependency clean package check, Windows/macOS parity, full BIEN scaling, all 50 posterior trees, a scientific shared-versus-per-tree comparison, refreshed historical sensitivity campaigns, CRAN readiness or publication.
