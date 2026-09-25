# After-task: freq vs BACE campaign under Rubin's rules (option B)

Lane `claude:pigauto:rubin` · branch `arc/rubin-freq-bace` · 2026-09-24 to 2026-09-25 · Claude Code (Opus 5.5)

## 1. Goal

Run the approved campaign (option B) to completion: frequentist arms (freq A proper, freq B improper) and BACE arms
(as shipped, chained, residual control) on the 18 core cells, 200 datasets per cell (BACE 100 at n = 1000), and deliver
a complete frequentist-vs-BACE report. Shinichi approved decisions 1 to 3 and asked to hold Dan until the campaign is in.

## 2. Implemented

- Chained BACE arm (`mi_bace_chain`, BACE's own `bace_final_imp` run M times from the previous draw; gate G-S3b).
- BACE input cleaning (`fit_bace_mi`: drop unused factor levels, leave out a discrete trait with fewer than 2 observed
  classes; recorded per fit; gate G-S3c). Without it BACE as shipped stops on most lambda = 1 datasets.
- Frequentist robustness (`rubin_freq.R`): Cholesky solves, bounds on every conditional distribution, a plausibility
  check on every Rphylopars fit (implied tip variance), bootstrap redraws counted (gate G-S2d).
- Fast exact PGLS (`est_pgls_slope_fast`, eigenbasis REML; equal to `nlme::gls` to 1e-5; gate G-S1d).
- prp per-cell scoring on the logit scale (Meng N10).
- Multi-cluster campaign tooling: `rubin_campaign_nibi.sh` (grid or RETRY_FILE, task cap), `rubin_totoro.sh`,
  `rubin_status.sh`, `rubin_retry_prep.R`, `rubin_pool.sh`, `rubin_missing.R`, `rubin_campaign_aggregate.R`,
  `rubin_report_template.html` + `rubin_report_build.sh`, `rubin_convcheck.R`.

## 3a. Decisions and Rejected Alternatives

- BACE settings runs 5, nitt 50,000 (Shinichi): BACE's own convergence check was shown not to measure mixing.
- Rerun only BACE failures that ran on the old code: failures on the fixed code are deterministic given the seed.
- Frequentist plausibility thresholds from measurement: refit vs original 50-fold (normal 0.5 to 1.7); original vs
  observed 1e4-fold (prp legitimately inflated up to 204-fold at lambda = 1 under the shared-lambda model).
- Compute moved between clusters by queue behaviour; nibi's queue and fir's low fair share were worked around.
- Rejected: installing packages into the shared v1 environment on nibi; rebuilding the library on narval.

## 4. Files Touched

`git diff --stat 3f28e26..HEAD`: 37 files, about 3,450 lines added, all under `script/rubin_*`, `script/tests-rubin/`,
`docs/dev-log/` and `.unlazy/rubin-freq-bace/`. Nothing in `R/` or `BACE/` (`git diff --stat 3f28e26..HEAD -- R/ BACE/` is
empty). The close commit adds the n = 1000 headline and data-driven freq A wording to both templates, the `bc` fix in
the report template (its recommendations section had never rendered), and the test fixes below.

## 5. Checks Run

- Completeness (`rubin_missing.R`): BACE 1200 / 1200 / 600 and frequentist 1200 / 1200 / 1200 files at n = 100 / 300 / 1000;
  none missing.
- Sanity line: complete-data estimates identical between BACE and frequentist files on all 6,000 matched datasets
  (max |diff| 0).
- Lane suite `testthat::test_dir("script/tests-rubin")`: FAIL 0, PASS 234 after the fix below. The first close run
  gave FAIL 3 in G-S1d (fast PGLS vs `nlme::gls`). Cause: two test files set `RNGkind("L'Ecuyer-CMRG")` and never
  reset it, so G-S1d built different datasets; on one (lambda 0.3, n 60, seed 502) the REML optimum is at lambda = 0
  and `gls` with free starts stops at 0.108 with a lower restricted log-likelihood (-66.474 against -66.450; a grid
  over [0, 1] agrees with the fast estimator). The fast estimator used for all campaign scoring is the correct one.
  Fixed: the two files restore the RNG kind, G-S1d pins it, and a new test asserts the fast estimator reaches at least
  `gls`'s log-likelihood at that boundary case.
- Both pages rendered headless in Chrome with no console errors and no NaN; `slop_check.py` 0 findings on both.

## 6. Tests of the Tests

- G-S3b negative control: c1's design matrix identical across shipped final runs, varying across chained runs.
- G-S3c positive control: the test dataset really has an empty level; as shipped BACE fails on it in 1 s.
- G-S2d: a non-PD parameter set and a 1e6-inflated Sigma_p are rejected; a mocked failing draw is redrawn and counted.
- G-S1d: the fast estimator reproduces `nlme::gls` on 36 datasets including lambda at the boundary.

## 7a. Issue Ledger

- FIXED: BACE as shipped fails on empty discrete levels (input cleaning).
- FIXED: frequentist blow-ups (about 1 in 300 fits; degenerate bootstrap refits and a silent ginv fallback).
- FIXED: PGLS scoring too slow at n = 1000 (58 s per fit with gls; milliseconds with the eigenbasis estimator).
- FOUND: BACE's convergence check tests the deterministic convergence phase and is scale-dependent.
- FOUND: the frequentist model's shared lambda inflates prp's variance at lambda near 1 (c1, c2 unaffected).
- Crashes and scheduler outcomes (sacct, rubin_ jobs since 2026-09-24 13:30): nibi 699 completed, 1 timeout; fir
  1162 completed, 1 failed, 1 cancelled (a duplicate), 2 redundant duplicates cancelled at close; rorqual 524
  completed, 1 failed, 4 timeouts, 1 cancelled. Every lost fit was rerun; the completeness check is the tally of record.
- FIXED: G-S1d test leaked RNG state and treated `gls` as ground truth at a boundary (see Checks).
- FIXED: report template referenced an undefined `bc`, so the recommendations, design note and n = 1000 limitation
  never rendered in versions 1 to 4 of the published report.

## 8. Consistency Audit

- Complete-data estimates identical between the frequentist and BACE result files for every matched dataset (same
  simulated data across arm sets and clusters).
- BACE fitting functions byte-identical on the Mac, nibi, fir, rorqual and Totoro despite version labels 0.0.0.9000
  and 0.1.0.
- Superseded results kept, not used: frequentist v1 (fir pigauto_rubin/results/freq), v2 (pigauto_rubin_f2), the nibi
  comparison run (pigauto_rubin_freq).

## 9. What Did Not Go Smoothly

- Frequentist blow-ups needed three passes (moment bounds, then plausibility, then a lambda-aware yardstick).
- `closeout.py`-style path handling aside, three of my own slips were caught: vectorised `format()` padding in the
  completeness check, a retry-prep rename that overwrote 6 pass-1 records (now per-pass folders), and a duplicate
  resubmission of 2 running fits (cancelled before they started).
- Queues: nibi held BACE for 2+ hours; fir started 1 of 354 n = 1000 tasks at 40 GB until resubmitted at 24 GB.
- narval's CPUs cannot run the fir-built library (illegal instruction).

## 10. Known Residuals

- BACE fits that error inside a file even after input cleaning are recorded, not scored: 77 datasets (47 at n = 100,
  20 at n = 300, 10 at n = 1000), so BACE rows rest on 1153, 1180 and 590 datasets rather than 1200, 1200 and 600.
  Frequentist: 3 Rphylopars failures of 3,600.
- At lambda = 1 and n = 1000, freq A and chained BACE both cover slightly above 0.95 (0.959, 0.958); downstream slope
  coverage at n = 1000 is 0.967 (freq A) and 0.975 (chained). Not investigated; the pages state it.
- The accuracy companion (CGGmi4Km6unGPhxEiv94qk) is rebuilt on the final data from `script/rubin_accuracy_template.html`
  via `rubin_report_build.sh`; its summary paragraph is now computed from the data. With the final data, chained BACE
  is as accurate as freq A at n = 300 and 1000 (0.565, 0.518 against 0.564, 0.525) and less accurate only at n = 100.
- Dan not contacted (Shinichi's hold). Pages are private until Shinichi shares them.

## 11. Team Learning

Memory receipt: pigauto AGENTS.md LOAD-FIRST (recovery to truth, prediction path first, r_cal untouched: no R/ change),
compute guards D-139 (estimate first, stop over 30%), D-143 (150 cores per user, checked on Totoro), D-64 (sockets only).

Golden Set: not in scope (no package code changed).

Durable lessons (candidates for LESSONS, pending Shinichi's approval): (1) a silent pseudo-inverse fallback hides
numerical failure; validate against hard bounds instead; (2) seeded BACE failures are deterministic, so retries only
help when the code changed; (3) request memory from measured peak plus margin: 40 GB blocked fir scheduling that 24 GB
unblocked.

## 12. Cross-Product Coverage

Covers: types_mixed DGP, BM with Pagel's lambda, MCAR 30%, n 100/300/1000, lambda 0.3/0.7/1, rho 0/0.5, M = 20,
five MI arms plus complete data, per-cell and two downstream estimands.

Does NOT cover: MAR or clade missingness, OU evolution, real trees, discrete-trait Rubin scoring, a pigauto arm,
BACE with other chain settings in the campaign.
