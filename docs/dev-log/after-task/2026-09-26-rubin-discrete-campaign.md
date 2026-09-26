# After-task: discrete traits added to the freq vs BACE Rubin campaign (n = 100 and 300)

Lane `claude:pigauto-rubin-freq-bace` · branch `arc/rubin-freq-bace` · 2026-09-26 · Claude Code (Opus 5.5)

## 1. Goal

Add the discrete traits (`bin`, `ord`, `cat3`) to the finished continuous campaign as one simulation for one
publication: the same cells, seeds and datasets, M = 20, a proper and an improper frequentist arm beside the same BACE
arms, per-value scores and a Rubin-pooled downstream estimate. Shinichi approved the plan and then the launch ("let's
go ahead and do the simulation"); n = 100 and 300 first, n = 1000 decided later.

## 2. Implemented

- `script/rubin_discrete.R`: an exact joint draw of the missing tips of an Mk trait (forward filtering, backward
  sampling, transition probabilities by the matrix exponential); castor arms mirroring freq A (parametric bootstrap
  of the rates, a refit rejected and redrawn when the observed tips are impossible under it) and freq B (fitted
  rates); flexible rate models (ARD for `bin` and `cat3`, SRD for `ord`) by default and v1's equal-rates models for
  comparison arms; per-value scores (modal-class accuracy, Brier, ECE, 95% class-set coverage and size, ordinal class
  distance, all ties fractional); the Rubin-pooled PGLS slope of c1 on `bin` with the complete-data interval and a
  per-dataset target.
- `script/rubin_cell.R`: default-off `--save_imp` (keeps every arm's M datasets) and `--discrete` (castor fills the
  freq datasets after their continuous draws under seeds 505/606, comparison arms under 707/808; every arm scored on
  discrete traits; per-value scores and the estimand recorded separately).
- `script/rubin_discrete_truth.R` (Monte Carlo population target, 5,000 datasets per cell),
  `script/rubin_discrete_rescore.R` (scores saved imputations afterwards), `script/rubin_disc_pool.sh`,
  `script/rubin_discrete_aggregate.R` (built by a workflow agent, checked by a second, revised by me),
  the discrete section of `script/rubin_study/build_data.R`, and the discrete section and findings of
  `script/rubin_study/study.qmd` (rendered `study.html`).
- `script/rubin_campaign_nibi.sh`: `RUBIN_OUT`, `CELL_FLAGS`, `TAG` in the job name; the task list is built in a
  temporary file so a dry run never rewrites the list queued tasks read.

## 3a. Decisions and Rejected Alternatives

- Frequentist arm: castor two ways, as approved. After the review, flexible rates became the primary model
  (measured: `bin` accuracy 0.60 against 0.50 for equal rates at λ = 0.3, n = 100) and v1's equal rates a comparison
  arm; both run in every dataset, so the labelling can be reversed without compute.
- Downstream coverage read against the per-dataset target, because complete-data intervals cover the population
  target in only 0.303 to 0.557 of datasets at λ = 1, ρ = 0.5 (0.905 to 0.985 for the per-dataset target).
- BACE only on the DRAC clusters (the stored fits' BACE build and R); Totoro ran the frequentist arms.
- The full arrays were submitted at once instead of a 20-fit pilot (Shinichi: time matters), with the first
  completions checked. nibi's 440 queued n = 100 tasks were moved to rorqual (400) and fir (480) at 10:12.
- Rejected: a `glmmTMB` probit arm for `bin` (one trait only, breaks the mirror); retrying BACE failures that the
  stored campaign also had (they are deterministic).

## 4. Files Touched

`git diff --stat f6f2406..HEAD`: 27 files, 12,482 lines added (most in the rendered `study.html` and the result
tables under `script/rubin_study/data/discrete/`), all under `script/rubin_*`, `script/tests-rubin/` and
`docs/dev-log/`. Nothing in `R/` or `BACE/` (`git diff --stat f6f2406..HEAD -- R/ BACE/` is empty). 15 commits.

## 5. Checks Run

- Completeness against the design: BACE 2,400 of 2,400 and frequentist 2,400 of 2,400 result files; none missing,
  duplicated or unexpected.
- Aggregator sanity: all 2,400 datasets identical between the BACE and frequentist files (mask and discrete traits;
  c1 within 3.7e-13); mode floor and complete-data estimate identical across a dataset's fits.
- `build_data.R` rebuild: the continuous tables are byte-identical to the committed ones (git shows no change).
- Lane suite `testthat::test_dir("script/tests-rubin")`: FAIL 0, PASS 466 (the 12 warnings are BACE's own).
- `study.qmd` rendered; every quoted number computed from the final aggregates (a script, not by hand);
  `slop_check.py` 0 findings on `study.qmd`, the plan and this report.
- `check-after-task.R`: structure check passed. Its acceptance-ledger step reports UNMET for this lane's six gate
  files (the continuous slices) and for another lane's five (`imputation-sim`, protected). `gate-check.mjs --status`
  reports all six of this lane's files ALL MET (16 gates). The ledger step re-runs every CHECK with `--reverify`,
  and the stored approvals are keyed to an earlier session's PATH and scratch paths, so they do not run here. Not
  re-approved: re-approval is a deliberate act, and one check needs a 2026-09-24 scratch folder.

## 6. Tests of the Tests

- The joint-draw test first compared with castor's own marginals and failed for SUEDE; brute-force marginals showed
  castor was wrong, not the sampler, and the test now compares with brute force.
- Boundary rate matrix (a rate exactly 0): analytic 2-state transition probabilities, semigroup property, and joint
  draws against brute force.
- Mocked castor failures: one trait failing leaves the others filled; bootstrap refits under which the data are
  impossible are rejected and counted, never returned as NA.
- A review workflow (four reviewers, a skeptic per finding) returned 15 findings; 13 survived and are fixed.

## 7a. Issue Ledger

- FIXED: defective rate matrices crashed the castor arm in about 30% of λ = 1, n = 100 freq A datasets (review).
- FIXED: bootstrap refits with zero-likelihood rates returned NA draws silently (found in my own smoke).
- FIXED: equal-rates castor imputes uniform classes at low λ (review; flexible rates now primary).
- FIXED: the c1 ~ bin estimand stopped when `bin` had one observed class and took the per-value rows with it.
- FIXED: a phantom arm "_er" in the aggregator (`paste0(character(0), "_er")` returns "_er").
- FOUND: castor's hsp marginals are wrong for non-reversible models (up to 0.47 for SUEDE); v1's ordinal
  probabilities were approximate.
- FOUND: BACE returns NA for 1.6% of its `cat3` draws (non-finite class probabilities in `.predict_bace`).
- FOUND: BACE fit success on fragile λ = 1 datasets depends on the machine (6 datasets differ between the runs).
- Scheduler: no failed, timed-out or out-of-memory tasks on nibi, fir or rorqual; rorqual refused one 880-task
  submission (1,000-task per-user cap) and took it split.

## 8. Consistency Audit

- BACE re-runs reproduce the stored continuous fits to the last digit on the same cluster (417 of 423); across
  clusters they are independent replicates (228 of 1,977 identical; median difference 0.004 to 0.016).
- Frequentist continuous rows: 2,236 of 2,400 files identical to the results of record within 1e-5; the rest differ by
  cross-host numerical paths (Totoro R 4.5.3 against fir R 4.5.0), up to 0.17 for one dataset.
- BACE failures: 62 of the 63 datasets without a BACE fit also failed in the continuous campaign.

## 9. What Did Not Go Smoothly

- My first BACE smoke was lost at the final save: I edited `rubin_cell.R` while Rscript was reading it.
- I read a stored failure record as a successful fit and resubmitted 8 known failures on fir; cancelled after
  4 minutes, outputs moved out of the pool.
- The monitor needed four fixes (macOS bash 3.2 has no associative arrays; a filter captured the whole status line;
  cancellations counted as failures; moved jobs' logs not matched).
- The session's permission check blocked the first job submission; Shinichi's explicit go-ahead cleared it.
- The destructive-command guard blocked an `rm -rf` of a scratch folder; used a new folder instead.

## 10. Known Residuals

- BACE has no fit on 63 datasets (all λ = 1: 44 at n = 100, 19 at n = 300) and no chained arm on 11 more; castor has
  none missing. Reported, not scored.
- Totoro's result files record `git_hash` "HEAD" (its copy sits in a git repository without commits); they ran
  55a0998. Fixed in the code for future runs (95dce1a), not synced mid-run.
- n = 1000 not run (about 3,400 core-hours of BACE); on hold (Shinichi, 2026-09-26). ECE is computed but not shown
  (biased at 30 to 90 cells).
- Pages updated with the discrete results (Shinichi, 2026-09-26; commit 648aa21): report
  https://claude.ai/artifact/PFkoRFtTjox4tndEbuBtPQ (v6), summary https://claude.ai/artifact/9vNrmz2PUfjW9vZwuKmvFc
  (v5), accuracy https://claude.ai/artifact/CGGmi4Km6unGPhxEiv94qk (v3); still private. Built from the templates via
  `script/rubin_report_build.sh` (discrete data: `script/rubin_disc_pages_data.py` to `data/discrete/pages.json`).
  The page build exposed a double rounding in `study.qmd` (largest lead 0.1446 written as 0.15), now 0.14.
- Brain notes not written (need Shinichi's approval): castor is pigauto's frequentist discrete arm; the finding of
  record; the lessons below.
- Dan not contacted (Shinichi's hold).

## 11. Team Learning

Memory receipt: pigauto AGENTS.md LOAD-FIRST (no `R/` change, recovery to truth); D-287 (estimate, pre-run, approval
before the full run); D-143 (Totoro 40 of 150 cores, checked by the launcher); D-64 (sockets only); no DRAC
login-node compute (only version queries and log reads there).

Golden Set: not in scope (no package code changed).

Durable lessons (candidates for LESSONS, pending Shinichi's approval): (1) never edit a script while Rscript or bash
is running it: both read it incrementally; run from a frozen copy; (2) before retrying a "stored success", read its
errors: a result file can be a failure record; (3) exact reproduction holds within a cluster, not across; place
re-runs on the stored fit's cluster when exactness matters; (4) macOS ships bash 3.2 (no associative arrays);
(5) check per-user submit caps (`sacctmgr ... MaxSubmitJobs`) before moving work between clusters.

## 12. Cross-Product Coverage

Covers: the continuous campaign's datasets at n = 100 and 300 (λ 0.3/0.7/1, ρ 0/0.5, 200 seeds), M = 20, three BACE
arms, four castor arms (flexible and equal rates, proper and improper), per-value discrete scores and a Rubin-pooled
c1 ~ bin slope with two targets.

Does NOT cover: n = 1000; downstream estimands for `ord` or `cat3`; a frequentist discrete imputer that uses the
other traits; non-threshold discrete DGPs; MAR or clade missingness; real trees.
