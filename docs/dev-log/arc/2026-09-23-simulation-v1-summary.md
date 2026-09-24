# Simulation v1: results and what the redo must change

Four-arm phylogenetic imputation study, run 2026-09-20 to 2026-09-23 on branch `arc/imputation-sim`
(draft PR #184). Stopped by Shinichi on 2026-09-23 to redo with longer MCMCglmm settings and, after the
lambda-default lane (D-278), a new pigauto baseline. This note is the hand-off to v2; the full record is
the methods note (`2026-09-20-simulation-methods-bace.md`), the article
(`vignettes/articles/simulation-study.Rmd`), the committed aggregates (`script/campaign_sim_results/`)
and the private board (https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24, v10).

## What was run

Seven traits of mixed type plus an always-observed MAR driver, simulated on `ape::rcoal` trees with
Pagel's lambda on the tree correlation and cross-trait correlation rho. Arms: the frequentist stack at
Rphylopars' `model = "BM"` default (`freq`) and at `model = "lambda"` (`freq_lambda`); BACE; pigauto
with the GNN off (in-house and Rphylopars joint solvers); pigauto with the GNN on; a mean or mode floor.

| stage | cells | replicates | complete |
|---|---|---|---|
| core: BM, MCAR 30%, lambda {0.3, 0.7, 1}, rho {0, 0.5}, n {100, 300, 1000} | 18 | 200 (BACE 100) | yes, G10 |
| factorial: BM/OU x 4 mechanisms x lambda {0.3, 1} x rho x n {100, 1000} | 56 | 200 (BACE 100 or 30 to 100) | yes, G11 |
| AVONET300 case study | 1 | 20 | yes |
| covariate sensitivity, 2 covariates on the core cells | 18 | 200 | yes |

Three of 600 core BACE replicates are hand-written failure records (one MCMCglmm segfault, two
four-hour timeouts), labelled in the files.

## Results

**1. The ranking inverts with phylogenetic signal.** Core slice, z-RMSE on the continuous family
(lower is better):

| n, lambda | BACE | freq (BM) | freq_lambda | pigauto GNN off | pigauto GNN on | floor |
|---|---|---|---|---|---|---|
| 100, 0.3 | 1.004 | 1.158 | 0.913 | 1.031 | 1.022 | 1.014 |
| 1000, 0.3 | 0.891 | 1.106 | 0.878 | 1.001 | 1.003 | 1.006 |
| 100, 1.0 | 0.794 | 0.584 | 0.566 | 0.487 | 0.541 | 1.018 |
| 1000, 1.0 | 0.647 | 0.475 | 0.464 | 0.379 | 0.407 | 1.002 |

At low signal only the two arms that estimate the phylogenetic share (BACE, `freq_lambda`) beat the
floor. At high signal pigauto with the GNN off leads and BACE is last.

**2. The frequentist stack's model choice matters more than the method.** Rphylopars at
`model = "lambda"` beats its BM default by 0.23 to 0.25 z-RMSE at lambda 0.3 and by 0.103 (MCSE 0.0009)
across the factorial. Against `freq_lambda`, BACE's continuous-trait advantage disappears. At the BM
default the stack sits above the mean floor at lambda 0.3.

**3. BACE's real advantage is discrete traits at low signal.** At lambda 0.3 it leads by 11 to 17
accuracy points (0.587 against 0.422 at n = 1000), and every other arm is at or below the mode floor.
At lambda 1 the frequentist Mk path and pigauto lead (0.96 to 0.97 against BACE's 0.85).

**4. Only conformal intervals reach nominal coverage.** pigauto: 0.95 to 0.96 at n >= 300. BACE 0.80 to
0.89, `freq_lambda` 0.88 to 0.91. Under clade-biased missingness pigauto's coverage drops to 0.86,
because masked cells are no longer exchangeable with the calibration split; the model-based Rphylopars
interval is the one that improves there.

**5. BACE is costly and fragile.** Errors on 16 to 20% of core replicates and up to 63% in the Brownian
clade-biased lambda = 1 cells; by its own `assess_convergence()` verdict only 15 to 20% of fits converge
at lambda 0.3, though median ESS is 572 to 1,600; about three hours per fit at n = 1000.

**6. Under clade-biased missingness BACE is inaccurate at low signal and fragile at high signal.** At
lambda 0.3 it sits above the floor with no failures (1.181 at n = 100); at lambda 1 most of its
replicates fail.

**7. Smaller results.** OU behaves like BM throughout. Two covariates improve the frequentist stack by
0.02 to 0.05 at lambda <= 0.7 and change no ranking. On AVONET the estimated lambda is about 1, so
`freq_lambda` equals `freq`, and pigauto's Rphylopars solver is the most accurate arm (0.498). That
solver returned absurd values on 84 simulated replicates, and the in-house solver on none.

## Why v2 is needed

1. **BACE settings were too short.** The current settings are `nitt = 50000`, `burnin = 10000`,
   `thin = 25`, `runs = 5`, `n_final = 20`. The low convergence-pass rate with healthy ESS points at
   drift across the five sequential imputation iterations, so `runs` needs to rise as well as `nitt`.
   BACE's vignette uses `runs = 15`, `nitt = 100000` for a real analysis, which is about 2.9x the cost.
   Pre-run a grid (`runs` 5/10/15 x `nitt` 50k/100k) on low-signal cells, measure the pass rate and wall
   time, then set the campaign.
2. **pigauto assumed lambda = 1** in its joint baseline (no residual component). D-278 orders a
   lambda-estimating default (Shinichi's lane, `../pigauto-lambda-default`). Redo the pigauto arms on
   that baseline.
3. **The monomorphic-discrete failures at lambda 1** (fixed thresholds) hit castor and pigauto GNN on.
   The GNN-on half is fixed on main (#185, keep the K axis for a single-level categorical trait); castor
   still fails by design and is scored at the floor.

## Correction found after close (2026-09-24): BACE's coverage is mostly a construction artefact

BACE's 95% interval in v1 is the 2.5 and 97.5 percentile over its `n_final = 20` imputed datasets
(`run_bace()` in `script/campaign_gnn_off_lib.R`). With 20 draws those quantiles cover a new exchangeable
draw only **0.872** of the time even if BACE's model were exactly right (measured by simulation, 200,000
replicates). The observed BACE coverage of 0.80 to 0.89 therefore sits at or near that ceiling in most cells,
so v1's statement that BACE's intervals fall well short of 0.95 is mostly our artefact, not BACE's. The plan
had flagged this (a reviewer's HIGH finding: build the interval from all retained MCMC samples, gate G6 asked
for at least 500) and the fix was never implemented; G6's resolution in the ledger waived that clause
wrongly. Two further unknowns: whether each imputed dataset is a posterior predictive draw (with residual
noise) or a draw of the conditional mean, which would push coverage lower still.

What still stands: pigauto's conformal coverage (0.95 to 0.96 at n >= 300) and the frequentist stack's
model-based coverage (0.80 to 0.91) are unaffected, and the accuracy results are unaffected. What does not
stand: any comparison of BACE's coverage with the other arms.

None of the v1 intervals use Rubin's rules, and none should: Rubin pools a downstream estimate (a regression
coefficient fitted on each completed dataset) as total variance `W + (1 + 1/M) B`, whereas v1 scores the
per-cell prediction interval for the missing value itself. v1 never measured a downstream estimand, so the
multiple-imputation-plus-Rubin workflow (`multi_impute()` -> `with_imputations()` -> `pool_mi()`) was not
evaluated for any arm.

**v2 must:** (1) build BACE's interval from enough draws that the construction ceiling is at 0.95 (many more
final imputations, or quantiles over the retained posterior predictive samples if BACE exposes them), and
confirm what one imputed dataset is; (2) add a downstream estimand (for example the slope of one continuous
trait on another, fitted with a phylogenetic GLS on every completed dataset) and score its bias and 95%
coverage after Rubin pooling for pigauto (`multi_impute`, M = 20 or more) and BACE (its `n_final` datasets),
with the complete-data fit as the reference.

## Operational lessons (keep these)

- Divergence rule in the aggregator: z-RMSE above 3x the floor or interval score above 1e3 is a failure
  that did not throw; floor it and count it (`n_divergent`).
- Resume-skip is per host: seed a host's results directory with the pool before submitting recovery
  arrays, or the recovery recomputes the campaign.
- Pools use per-host subdirectories; dedupe on (filename, arm set); never aggregate while an rsync runs.
- Clade-masked BACE at n = 1000 needs a 12 h Slurm limit, not 5 h; nibi caps submissions at 1,000 array
  tasks per user.
- Totoro is shared: another user's process held 930 GB on 2026-09-22; read per-user RSS before blaming
  your own wave, and guard RAM.
- Run `freq_lambda` from the start; it became the fifth arm mid-campaign.

## Decisions recorded (G13c, Shinichi 2026-09-22 and 2026-09-23)

The BACE paper reports both frequentist specifications with the gap as a finding and leads with the
discrete result, with the failure and convergence costs beside it. The article stays unlisted until
Szymek signs off on the corrected Pagel-lambda form. PR #184 stays a draft until Shinichi has read the
board. The joint baseline estimates lambda by default (D-278).
