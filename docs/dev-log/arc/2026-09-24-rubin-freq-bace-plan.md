# Plan: frequentist stack vs BACE, re-run with proper multiple imputation and Rubin's rules

Lane `arc/rubin-freq-bace` (worktree `../pigauto-rubin-freq-bace`), branched from `arc/imputation-sim` at
`14b919d` so the v1 runner, design table and aggregator come with it. Scope: the frequentist stack and BACE
only, for Dan Noble's BACE paper. No pigauto arm (a pigauto posterior is being built in another lane and can
plug into the same estimands later). Owns `script/rubin_*`, `docs/dev-log/arc/2026-09-24-rubin*`,
`.unlazy/rubin-freq-bace/`. Never edits `R/` or `BACE/`.

## Why (from v1, 2026-09-24)

> **Correction (2026-09-24, later the same day; supersedes point 2 below).** The installed BACE, the
> code `BACE::bace()` actually runs (built 2026-08-09, `~/Library/R/arm64/4.6/library/BACE`), does draw a
> residual. `bace_final_imp()` calls `.predict_bace(..., sample = TRUE)`, and for a gaussian trait that
> branch takes one posterior iteration's fitted value and adds `rnorm(0, sqrt(sigma2_units))` from the same
> iteration before back-transforming. Each final dataset is a posterior predictive draw, but only conditional
> on a shared anchor: every final run starts from the same converged dataset (see the review note below). Point 2
> was verified against the in-tree `BACE/` clone, which is stale (commit `de87d8c`, 2026-04-01) and still
> takes the posterior mean. Checked by deparsing the installed namespace, after the S3 builder flagged it.
> Consequences: the "BACE + residual draw" arm adds a second residual on top of BACE's own, so it is kept
> in the smoke only as a measured contrast and its place in the campaign goes back to Shinichi; nothing
> should reach Dan saying BACE imputes posterior means; and every compute host must run the same BACE build
> (the pre-run checks for `sample = TRUE` in `bace_final_imp`). Point 1 (the percentile interval's 0.872
> ceiling) is unaffected.

> **Review note (Meng, 2026-09-24; `docs/dev-log/arc/2026-09-24-rubin-review.md`).** In the installed build
> every final run starts from the same `last_data`, the last convergence iterate, whose fills came from
> `sample = FALSE` runs. Predictors later in the formula order are therefore identical across the M datasets
> (c1's design matrix was identical in all 20 final runs of the smoke fit). In a tree-free toy, this scheme
> covered the slope 0.900 and cells 0.88 to 0.90, against 0.95 for independent chains; chaining final run i
> from dataset i-1 restored 0.950. The size in the phylogenetic DGP is not measured. A `bace_chain` arm built
> from BACE's own functions is the proposed replacement for `bace_resid`; that is Shinichi's call.

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
| BACE as shipped | its `n_final` datasets (installed build: one posterior draw plus residual per final run, all runs starting from one converged dataset) | what a BACE user gets today |
| BACE + residual draw (negative control after the correction) | BACE's datasets plus a second draw from N(0, sigma_units) | how coverage responds to a doubled residual; a chained BACE arm is the proposed replacement (review B2) |
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
3. Downstream correlation: GLS-whitened phylogenetic correlation of c1 and c2 on Fisher's z (normal
   reference), pooled the same way. A tip-level Pearson correlation was replaced in planning because its
   1/(n-3) variance is wrong under phylogeny.
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

## Decisions (Shinichi, 2026-09-24: "yes yes yes" to Q1, Q2, Q3 below)

- Q1: *"Add a 'BACE + residual draw' arm built in our runner from BACE's returned model objects, so the paper
  can separate BACE-as-shipped from BACE-with-proper-draws? yes / no"* (Tell Dan either way: this is a
  property of BACE's imputations that his paper should state.)
- Q2: *"Frequentist MI: parametric bootstrap (A) as the headline, draws-from-one-fit (B) as the contrast?
  yes / no"*
- Q3: *"Downstream estimands: PGLS slope of c2 on c1 plus their correlation, on the core cells? yes / no"*

## Compute plan: Totoro + DRAC (Shinichi 2026-09-24: "plan simulations using DRAC + totoro")

Split by what each machine is good at, measured in v1:

| work | where | why | v1 measurement it rests on |
|---|---|---|---|
| step 1 smoke (n = 60) | Mac | minutes | |
| frequentist arms, all cells (A needs M = 20 refits per replicate) | **Totoro**, `script/campaign_sim_totoro.sh` pattern, <= 150 cores (D-143; 250 was a temporary allowance) | no queue; v1's freq_lambda wave did 11,200 replicates in 12 min, so A at 20 refits is roughly 4 h for the core | freq_lambda 2.0 / 3.9 / 10.6 s per replicate at n = 100 / 300 / 1000 |
| step 2 BACE settings pre-run (n = 100, 300) | **Totoro** | fits are 20 to 50 min, RAM modest, results within hours | BACE 1,161 / 3,050 s at n = 100 / 300 |
| BACE campaign at n = 100 and 300 | **nibi** job arrays (`script/campaign_sim_nibi_array.sh` pattern), 1 core, 32 GB per task | many independent long single-threaded fits; nibi drained v1's BACE fastest | median 1,161 / 3,050 s at v1 settings; scale by the chosen `nitt` x (`runs` + `n_final`) |
| BACE campaign at n = 1000 | **nibi** arrays, 1 core, 32 to 48 GB, `--time` 12:00:00 | 3 h per fit at v1 settings, 5 to 6 h for clade-masked; fir as backup only | 11,019 s median at n = 1000; clade cells timed out at 5 h |

Rules carried from v1 (see `docs/dev-log/arc/2026-09-23-simulation-v1-summary.md`, operational lessons):
read `free -g` and per-user RSS on Totoro before every launch (another user held 930 GB on 2026-09-22; the
lambda-default lane's benchmark used about 140 cores on 2026-09-23); nibi caps a user at 1,000 submitted array
tasks, so split arrays by n; resume-skip is per host, so seed a host's results directory from the pool before
any recovery array; per-host pool subdirectories and (filename, arm set) dedupe; never aggregate while an rsync
runs; a segfault or a timeout leaves no rds, so the runner cannot floor it (record by hand and label it). No
DRAC login-node compute; attach through the `~/.ssh/cm-*` ControlMaster sockets, never trigger Duo.

The full-campaign budget is re-derived from step 2 and approved by Shinichi before anything beyond the
pre-run is submitted (D-139).
