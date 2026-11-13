# Seed validation, mammal diagnostics and release follow-up

## 1. Goal
Do the feasible follow-up work on seed sensitivity, missing mammal evidence and the separate CRAN release, while preserving ownership and scientific claim boundaries.

## 2. Implemented
Completed 35 CTmax validation fits on Totoro and a separate ten-fit factor-versus-numeric coding check on the reporter's original five folds. Added reviewed, tested no-refit diagnostic helpers for mammal issues #250 and #251. Verified existing exact-archive release logs and prepared updated draft comments without changing the release owner's files. RESULTS.md gives the bounded scientific conclusions.

## 3a. Decisions and Rejected Alternatives
Fixed outer folds and seed trios before fitting. Preserved default package behaviour. Rejected a new default, ensemble intervals, downstream pooling and a post-hoc subgroup fallback presented as validated. Kept the first pilot separate because its metadata predates strict identity checks. The planned production repeat of its fold/seed reproduced point predictions exactly.

Assumptions made without asking: use a new stratified row-holdout design on the same public dataset, six fixed seeds and two fixed trios; keep internal calibration randomness as part of end-to-end seed sensitivity; use Totoro CPU at six fits and two torch threads each. These choices do not recreate the original 375-target experiment or prove behaviour under uniform whole-species holdouts.

## 4. Files Touched
New files under script/seed_validation_260/: prepare.R, run.R, summarize.R, launch.py, selftest.R, check_factor_coding.R and stratify.R. New files under script/issue_250_diagnostics/: diagnose_ties.R, selftest.R and README.md. New files under script/issue_251_diagnostics/: subset_intervals.R, selftest.R and README.md. New files under docs/dev-log/seed-validation-260/: PLAN.md, RESULTS.md, RELEASE-STATUS.md, cran-comments-followup-DRAFT.txt, source-manifest.json and this report.

Remote raw inputs, fits, logs and summaries remain on Totoro under ~/sim-archive/pigauto/seed260-20261108/. The original dirty checkout, package runtime, shared instruction files and frozen release archive were not edited. Local branch: codex/seed-validation-260. The initial study was kept local. On 9 November the user authorized pushing it, opening a PR and updating #259/#260. The dedicated RELEASE-STATUS.md and cran-comments-followup-DRAFT.txt remain local and are excluded from the PR's final diff; their earlier audit is recorded here as context only. Release operations remain separate.

## 5. Checks Run
The final metadata-enabled one-epoch smoke passed in 2.6272 seconds. The pilot completed 2,000 epochs in 250.3042 seconds. Production completed all 35 expected fits with finite point predictions and intervals, within its 40-minute limit. Actual GNN times were 225–260 seconds. Five fits used ordinary early stopping, so 2,000 was requested rather than universally completed. Five captured warnings were the expected GNN-off covariate-route advisory.

The 61-file runtime source/data/metadata manifest matched exact model source 033649b93115c6d8191e791f06b1ddf52eb82d19. The harness self-test, both mammal self-tests and an eight-tip real pigauto_result tie-helper check passed. Factor/numeric coding agreed on all five original CTmax folds under the declared 1e-5 tolerance; actual completed-value differences were at most 2.57e-10 degrees C.

Independent reviewers recomputed all primary and post-hoc metrics, row/truth and source/input/runner identities, and attempt counts. Agreement with summaries was within 1e-12. git diff --check passed. Prose checks passed with no em dashes. The after-task structural function was run separately; no claim is made that an unrelated historical acceptance ledger was reverified.

No exports, roxygen or installed-package paths changed. script/ and docs/ are explicitly excluded by .Rbuildignore. The unchanged package-source check evidence from the earlier source remains separate from these runner and scientific checks; no duplicate full package check or new release check was launched.

## 6. Tests of the Tests
The harness refuses wrong prepared-input identity, mixed runner identity, corrupt records and missing attempts while retaining a complete attempt ledger. Plan-only tests reused 35 valid records, then 34 with one queued job, and rejected wrong source before launching workers. Launcher limits and interruption process-group cleanup were independently reviewed.

The tie helper distinguishes exact from rounding-only ties and rejects invalid completions, a mismatched saved tree and changed observed values. The interval helper rejects invalid correlations, zero-width reference intervals and unavailable/nonfinite/missing m levels. A reviewer caught silent omission of unavailable levels; the corrected helper now refuses the grid.

## 7a. Issue Ledger
#260: fresh-seed point validation complete for this one row-holdout regime; averaging improved aggregate error but not the singleton subgroup, so the general remedy remains open. #259: coding invariance now checked on the original public application folds, GNN off. #250 and #251: helpers ready, actual manuscript reproducer/draws still unavailable. #261: reporter confirms the GitHub fix; CRAN availability remains separate. #135: broader mixed-type and multi-dataset comparative evidence remains open. No issue was closed or newly messaged.

## 8. Consistency Audit
Compared observation-weighted and species-balanced summaries, retained and entirely withheld target observations, newt inclusion/exclusion, individual and fixed-trio predictions, and individual interval summaries. No ensemble intervals were inferred. The post-hoc subgroup result is explicitly labelled and was not used to tune a model.

Read exact-archive local, Windows release/devel and R-hub Ubuntu release/clang logs and matching hash files. Found stale PENDING entries in the release owner's draft comments and prepared an update separately. Flagged the AGENTS.md versus current code/help covariate description for the owner's product-contract reconciliation.

## 9. What Did Not Go Smoothly
Review found missing source/input binding, mixed-runner acceptance and incomplete attempt accounting in early harness drafts; these were fixed before production. A cross-process fixture initially disappeared with R's temporary directory, so the reusable self-test now runs subprocesses before R exits. The old pilot was not admitted to final results. The local graph was stale for the active model source; the exported source and manifest governed the run. No new Win-builder result was found in targeted Gmail searches, which does not establish absence elsewhere.

## 10. Known Residuals
The true causes of mammal ties and interval narrowing remain unknown without inputs. The seed study is one dataset, one tree and one row-holdout design; singleton status and lack of target training values coincide here. No general calibration or universal ensemble improvement is established. No new default, package API, public scientific post, release tag or CRAN submission occurred. Independent Win-builder review and the final release gates remain with the existing Pigauto chat. Formal exact-draft prose-assessment compiler: NOTASSESSED.

## 11. Team Learning
Memory receipt: the canonical contract was already loaded; compute-routing, CRAN release gate, quantitative-analysis, validation-harness and current ownership checks shaped the work. Prior issue state and existing release work were reused. No brain memory was written. Golden Set: not run for a runtime repair because no package implementation was changed; planted experiment/diagnostic failure controls were run instead.

The independent numerical review passed, including the post-hoc partition. Reconciliation: the intended prospective seed comparison completed, with an explicitly post-hoc breakdown that limits the interpretation. Data-dependent investigations and release work remained separate rather than being treated as completed by a successful campaign.

## 12. Cross-Product Coverage
This work covers one continuous CTmax trait, repeated observations, six GNN seeds on fixed row folds, two fixed point averages, baseline factor/numeric coding on the original application folds, and synthetic/real-object tests of no-refit helpers. It does NOT cover all trait types, a uniform species/clade-holdout campaign, new ensemble uncertainty or inferential validity, the actual private mammal examples, general calibration, public release or CRAN submission. Those are distinct unfinished targets.

Review note: self-review of this whole factual report on 8 November 2026, perceived AI-like style 3/10 with moderate confidence. Located distinctions “2,000 was requested” and “explicitly post-hoc” prevent overclaiming. Science and factual gates pass within independently verified artifacts; reference identity is given by exact source/archive hashes and the linked sibling reports. Formal exact-draft compiler remains NOTASSESSED.
