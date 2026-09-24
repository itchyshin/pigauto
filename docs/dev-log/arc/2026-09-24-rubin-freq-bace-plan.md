# Plan: frequentist stack vs BACE, re-run with proper multiple imputation and Rubin's rules

Lane `arc/rubin-freq-bace` (worktree `../pigauto-rubin-freq-bace`), branched from `arc/imputation-sim` at
`14b919d` so the v1 runner, design table and aggregator come with it. Scope: the frequentist stack and BACE
only, for Dan Noble's BACE paper. No pigauto arm (a pigauto posterior is being built in another lane and can
plug into the same estimands later). Owns `script/rubin_*`, `docs/dev-log/arc/2026-09-24-rubin*`,
`.unlazy/rubin-freq-bace/`. Never edits `R/` or `BACE/`.

## Why (from v1, 2026-09-24)

1. v1 scored per-cell prediction intervals, not multiple imputation. BACE's interval was the 2.5 and 97.5
   percentile of 20 imputed datasets, whose coverage ceiling is 0.872 under a correct model; the Rubin
   interval from the same 20 draws, mean plus or minus `qt(0.975, M - 1) * sqrt((1 + 1/M) * B)`, covers 0.950
   (200,000 simulated sets).
2. **BACE's continuous imputations are posterior means, not predictive draws.** Verified in the source:
   `.predict_bace()` (BACE/R/model_functions.R) takes `pred_prob[,1]`, "Extract posterior mean", from
   `.pred_cont()`, which returns the posterior mean, SD and credible interval of the *fitted values*
   `X beta + Z u`. No residual variance is drawn. So the spread across BACE's `n_final` datasets comes from
   refitting and from the chained imputations of other traits, and excludes the residual term entirely. As MI
   in Rubin's sense this is improper for continuous and count traits: Rubin intervals built on them should
   under-cover, and more so where the residual share is large (low lambda). Discrete traits are sampled from
   predicted probabilities (`sample = TRUE` path), which is a draw. How large the continuous shortfall is has
   not been measured (AGENT-INFERRED direction, not magnitude).
3. v1 measured no downstream estimand, which is the thing Rubin's rules exist for.

Consequence for the paper: "BACE vs frequentist" under Rubin's rules must say what each arm's imputations
are. Changing BACE is out of scope for this repo (and is Dan's code); options for the BACE arm are in Q1.

## Arms

| arm | imputations (M = 20) | what it tests |
|---|---|---|
| BACE as shipped | its `n_final` datasets (posterior means for continuous) | what a BACE user gets today |
| BACE + residual draw (if Q1 = yes) | each posterior-mean value plus a draw from N(0, sigma_units) using that fit's residual variance, taken from the returned model objects | whether proper predictive draws fix BACE's coverage |
| freq A: parametric bootstrap then draw | refit Rphylopars `model = "lambda"` (and castor Mk, phyloglm) on a parametric bootstrap sample, then draw missing cells from that fit's conditional normal | proper MI: carries Sigma and lambda uncertainty |
| freq B: draws from one fit | M draws from N(prediction, anc_var) of a single fit | improper contrast: parameters fixed |
| floor | mean / mode, no MI | reference |

`freq` at the BM default is dropped from the headline: v1 already showed the specification gap, and the
paper reports both specifications in its v1 tables.

## Estimands and scoring (all Rubin-pooled across M)

1. Per-cell: missing value; interval mean plus or minus `t_{M-1} * sqrt((1 + 1/M) B)`; coverage, width,
   interval score. (Within-imputation variance is zero for the value itself.)
2. Downstream slope: phylogenetic GLS of continuous trait c2 on c1 (`nlme::gls`, `corPagel` on the true
   tree) fitted on every completed dataset; pooled with `W + (1 + 1/M) B`, Barnard-Rubin df. Score bias,
   RMSE and 95% coverage against the same model fitted on the complete true data.
3. Downstream correlation: Pearson correlation of c1 and c2 on Fisher's z, pooled the same way.
4. Diagnostics: fraction of missing information per estimand; BACE convergence verdict and ESS as in v1.

## Design

Core cells of v1 only (BM, MCAR 30%, lambda {0.3, 0.7, 1}, rho {0, 0.5}, n {100, 300, 1000}); rho = 0.5 is
where a downstream slope exists, rho = 0 checks that the pooled slope stays near zero. 200 replicates for the
frequentist arms, 100 for BACE, as v1. Factorial deferred until the core answers.

BACE settings are not fixed yet: step 2 below chooses them.

## Steps

1. **Smoke (Mac, minutes):** one cell, n = 60, both frequentist arms and BACE with M = 20; confirm every
   estimand is produced and the Rubin pooling reproduces a hand computation.
2. **BACE settings pre-run (Totoro, D-139):** core cells at lambda 0.3 and 0.7, n = 100 and 300, 5 seeds,
   grid `runs` {5, 10, 15} x `nitt` {50k, 100k} (burnin 20%, thin keeping about 1,600 samples); record BACE's
   convergence verdict, drift, median ESS, wall time and the per-cell Rubin coverage. Estimate: 24 cells x 6
   settings x 5 seeds = 720 fits at a mean of about 25 min = 300 core-hours, about 2 h at 150 cores.
   **Stop and show Shinichi the numbers and the re-derived budget.**
3. Campaign on the approved settings. Budget re-derived from step 2; BACE at n = 1000 dominates it.
4. Aggregate, report for Dan: a results section in the methods note, board update.

## Open questions for Shinichi (drafted)

- Q1: *"Add a 'BACE + residual draw' arm built in our runner from BACE's returned model objects, so the paper
  can separate BACE-as-shipped from BACE-with-proper-draws? yes / no"* (Tell Dan either way: this is a
  property of BACE's imputations that his paper should state.)
- Q2: *"Frequentist MI: parametric bootstrap (A) as the headline, draws-from-one-fit (B) as the contrast?
  yes / no"*
- Q3: *"Downstream estimands: PGLS slope of c2 on c1 plus their correlation, on the core cells? yes / no"*

## Compute and constraints

Totoro is shared: read per-user RSS and `free -g` before every launch; the lambda-default lane's benchmark
was using about 140 cores on 2026-09-23. nibi for anything queued. No DRAC login-node compute; no Duo.
