# Simulation plan: phylogenetic trait imputation, frequentist stack vs BACE vs pigauto

Status: **plan only, nothing launched** (Shinichi, 2026-09-19: "we plan for it"). Framework: ADEMP (Morris, White and Crowther 2019, *Stat. Med.* 38: 2074) with the 11 reporting items of Williams et al. (2024, *MEE* 15: 1926). Feasibility evidence: today's with/without-GNN campaign (200 cells), the solver diagnostic (200 cells), and the per-type small test (15 cells; appendix).

Two papers from one set of runs. The BACE paper reports arms 1 and 2. The pigauto paper, after BACE is published, reports arms 1 to 4. Every replicate stores every arm's output, so the second paper needs no new compute.

## A. Aims

- **Primary (BACE paper).** Does joint Bayesian imputation (BACE) improve accuracy, probability calibration and interval coverage over the standard frequentist stack (Rphylopars for continuous-family traits, castor's Mk model for discrete traits), across trait types, phylogenetic signal, cross-trait correlation, missingness mechanism and sample size?
- **Primary (pigauto paper).** Does pigauto with the GNN off match or exceed the better of arms 1 and 2 on the same three targets at a fraction of the cost, and what, if anything, does the GNN add?
- **Secondary.** Robustness to a mis-specified evolutionary model (OU traits, low signal), to MAR missingness, and to tree size; the price of joint modelling when traits are independent.

## D. Data-generating mechanism

Levels: one tree per replicate; one row per species (single observation); K traits per species.

1. Tree: `ape::rcoal(n)` scaled to unit height (ultrametric; required by BACE). Sensitivity: `rlineage`-style birth-death trees via `BACE::sim_bace`'s generator.
2. Latent traits: `L ~ MVN(0, Sigma_phylo (x) V_lambda + Sigma_resid (x) I)`, with `V_lambda = lambda * V_tree + (1 - lambda) * I` (Pagel's lambda on the tree covariance), `Sigma_phylo` a K x K matrix with unit diagonal and off-diagonal `rho`, and residual variance `(1 - lambda)`. This is the mixed-signal, cross-correlated liability model both BACE and pigauto assume; OU traits replace `V_tree` by the OU covariance with `alpha = 2`.
3. Trait construction from the latents (one latent per trait): continuous = latent; count = `Poisson(exp(1.5 + 0.8 L))`; proportion = `plogis(L + N(0, 0.3^2))`; binary = `L > median`; ordinal = quartile cuts of `L` into 4 levels; categorical = tercile cuts into 3 levels. Type mix per replicate: 2 continuous, 1 count, 1 proportion, 1 binary, 1 ordinal, 1 categorical (K = 7), the `types_mixed` DGP already in `script/campaign_gnn_off_lib.R`.
4. Missingness applied to the complete data, the same mask for every arm within a replicate: MCAR at a fraction `m`; MAR where the probability a cell is missing depends on an always-observed continuous trait (logistic, odds ratio 3 per SD), same marginal fraction.
5. Factors and levels:

| factor | levels | why |
|---|---|---|
| n (species) | 100, 300, 1000 | the campaign's ladder; 1000 is where BACE costs 20 min per fit |
| lambda (signal) | 0.3, 0.7, 1.0 | low-signal regime is where the mean floor and lambda shrinkage decide (solver diagnostic) |
| rho (cross-trait) | 0, 0.5 | joint methods should gain only when rho > 0 |
| missingness | MCAR 0.1, MCAR 0.3, MAR 0.3 | 0.3 matches the campaign; MAR is the realistic case |
| evolutionary model | BM, OU (alpha = 2) | sensitivity |
| covariates | none; 2 environmental covariates (effect 0.5 SD each on every trait's latent, phylogenetic signal 0.3 in the covariates themselves, fully observed) | BACE and pigauto can use them; the frequentist stack needs phylolm / phyloglm to use them at all |

Full factorial without covariates: 3 x 3 x 2 x 3 x 2 = 108 cells. The BM / MCAR-0.3 / no-covariate slice (18 cells) is the core; the rest are sensitivity. The covariate factor doubles the design; the recommendation is to run it on the core slice only (18 more cells) unless the BACE paper wants covariates as a headline.

Covariate handling per arm, to be settled with Szymek:

- BACE: covariates enter every trait's fixed formula (`y ~ x1 + x2 + other traits`), its native use.
- pigauto: `covariates =` enter through the GNN (`obs_refine`, `cov_linear`), so arm 4 uses them and arm 3 (`gnn = FALSE`) ignores them with a warning. A covariate-aware baseline (phylogenetic GLS with fixed effects) exists in the codebase (`bench_covariate_sim` era) but is not the default `gnn = FALSE` path; if the covariate factor is in, arm 3 needs that route or is reported as "no covariates" honestly.
- Frequentist stack: Rphylopars and castor take no covariates. The fair comparator with covariates is `phylolm::phylolm(y ~ x1 + x2, model = "lambda")` per continuous trait, `phylolm::phyloglm` per binary trait, and castor for multi-state traits (no covariates), i.e. a second version of arm 1. Without that, a covariate DGP would measure "has covariates" rather than "is Bayesian".

6. Replicates per cell: **200** for arms 1, 3, 4 (coverage MCSE `sqrt(0.95 * 0.05 / 200)` = 1.5 points, so a 5-point coverage difference is 3 MCSE); **100** for BACE (MCSE 2.2 points) unless the pre-run shows it affordable at 200. Reported per cell with MCSE.

## E. Estimands and targets

For each masked cell (species i, trait k) the truth `y_ik` is stored with the replicate. Targets, per trait type:

- Point accuracy: z-RMSE on continuous-family traits (`(y - yhat) / sd_train`), accuracy and macro-F1 on discrete traits. Both are computed on the masked cells only.
- Probability calibration (discrete traits): Brier score and expected calibration error (10 bins) of the predicted class probability against the outcome. Arm 1 supplies castor's normalised state likelihoods; arm 2 the posterior class frequencies over draws; arms 3 and 4 pigauto's `probabilities`.
- Interval coverage and width (continuous-family traits): 95% interval coverage `mean(lo <= y <= hi)` and mean standardised width. Arm 1: `yhat +/- 1.96 * sqrt(anc_var)` from Rphylopars; arm 2: 2.5 and 97.5 percentiles over BACE draws (n_final = 20); arms 3 and 4: pigauto's conformal interval (production mode) and, as a secondary, the BM-SE interval.
- Cost: wall time per fit; failure rate per arm (a failed fit counts as a failure, never dropped).

## M. Methods (arms)

1. **Frequentist stack.** `Rphylopars::phylopars(model = "BM", phylo_correlated = TRUE, pheno_correlated = TRUE, REML = TRUE)` on the continuous-family columns jointly (count on log1p, proportion on logit, back-transformed); `castor::hsp_mk_model` per discrete trait (ER rates for binary and categorical, SUEDE for ordinal, 3 trials). Standard practice; no cross-type information.
2. **BACE** (Bayesian, joint, MCMCglmm): `nitt = 50,000`, `burnin = 10,000`, `thin = 25`, 2 chains, `n_final = 20`, OVR on. The chain length is fixed by a convergence pre-check (Gelman-Rubin on the pre-run cells), not by cost.
3. **pigauto, GNN off** (`gnn = FALSE`, default safety machinery, `joint_solver` as decided in arc C; both solvers if the default is still open).
4. **pigauto, GNN on** (defaults, 2000 epochs) and, stored for free, GNN on predicting from the full baseline.
Reference: mean/mode floor.

## P. Performance measures and MCSE

| measure | formula | MCSE |
|---|---|---|
| z-RMSE (per trait, per cell) | `sqrt(mean(z^2))` over masked cells, then mean over replicates | `sd / sqrt(n_sim)` over replicate means |
| accuracy | `mean(yhat == y)` | `sd / sqrt(n_sim)` |
| Brier | `mean(sum_k (p_k - 1[y = k])^2)` | `sd / sqrt(n_sim)` |
| ECE | `sum_b (n_b / N) * abs(acc_b - conf_b)`, 10 bins | bootstrap over replicates |
| coverage | `mean(lo <= y <= hi)` | `sqrt(p (1 - p) / n_sim)` |
| width | `mean((hi - lo) / sd_train)` | `sd / sqrt(n_sim)` |
| wall, failures | seconds; `mean(failed)` | `sd / sqrt(n_sim)`; binomial |

Reporting: one table row per cell x arm x trait type with metric (MCSE); figures as in `script/campaign_gnn_off_figures.R` (dot and MCSE bars, facetted by factor); a worked case study on AVONET300 (real data, item 9), where today's numbers already exist for all four arms.

## Compute estimate (D-139)

Measured per-fit walls today (4 threads per cell): BACE 124 / 477 / 1,131 s at n = 100 / 300 / 1000; pigauto GNN on 118 / 165 / 526 s; pigauto GNN off 0.5 / 4.5 / 117 s (the n = 1000 figure is the gate calibration; the pure arm takes 2 s); frequentist stack under 3 s at every n.

Per replicate, arms 1 + 3 + 4: about 122 / 173 / 646 s; BACE alone: 124 / 477 / 1,131 s.

Full factorial without covariates, 36 cells per n, arms 1, 3, 4 at 200 replicates and BACE at 100:

| n | arms 1, 3, 4 | BACE | total |
|---|---:|---:|---:|
| 100 | 36 x 200 x 122 s = 244 core-h | 36 x 100 x 124 s = 124 core-h | 368 core-h |
| 300 | 36 x 200 x 173 s = 346 core-h | 36 x 100 x 477 s = 477 core-h | 823 core-h |
| 1000 | 36 x 200 x 646 s = 1,292 core-h | 36 x 100 x 1,131 s = 1,131 core-h | 2,423 core-h |
| all | | | **about 3,600 core-hours** |

Totoro at 144 threads (D-143 cap): about 25 h wall. A DRAC job array at 500 cores: about 8 h. The core slice (18 cells, no covariates) is a third of that; the covariate factor on the core slice adds about the same again, plus phylolm / phyloglm fits (seconds).

Pre-run before any submit (D-139): one replicate per core cell at n = 100 and n = 1000, all arms, with Gelman-Rubin diagnostics on BACE's two chains, output inspected, estimate re-stated, then Shinichi approves. Target: Totoro for the core slice, a DRAC array for the factorial.

## Williams et al. (2024) reporting items, self-audit

| item | covered by |
|---|---|
| 1 aims | section A, two primary aims and secondary |
| 2 DGP | section D, math and factor table |
| 3 estimands | section E, truth stored per replicate |
| 4 methods | section M, four arms plus floor, settings stated |
| 5 performance measures | section P, formulas |
| 6 software and versions | `sessionInfo()` saved beside results; pigauto commit, BACE, castor, Rphylopars versions in each rds (runner already records pigauto version and host) |
| 7 code availability | `script/` on pigauto main; BACE paper repo links to it |
| 8 seeds and reproducibility | master seed per cell; masks derived from `seed + 1000` (already the campaign convention) |
| 9 real-data case study | AVONET300, all arms, in the main text |
| 10 full results incl. failures | per-cell table with failure rate; no fit dropped |
| 11 MCSE on every number | section P table; aggregator computes MCSE for every summary |

## Open before launch (decide with Shinichi and Szymek)

- Phylogenetic signal: lambda in {0.3, 0.7, 1.0} on the tree covariance is one parameterisation; BACE's own simulator uses a "phylo_signal" fraction of variance. Agree one definition and state it in both papers.
- Covariates: in or out for the BACE paper; if in, arm 1 becomes the phylolm / phyloglm / castor stack and arm 3 needs a covariate-aware non-GNN route or is reported as covariate-free.

- Arc C default: `joint_solver` for arm 3 (Rphylopars solver inside pigauto, 30% better on real data at 50x the fit time and a Suggests dependency, or the in-house solver repaired first).
- BACE chain length: fix by convergence, then re-estimate compute.
- Whether to add a single-trait-type sensitivity (each type alone) for the BACE paper's per-type claims; the `types_mixed` DGP already yields per-type metrics from one run.
- Compute target: Totoro core slice first, DRAC array for the factorial.
