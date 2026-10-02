# Defaults audit for the user-facing arguments of pigauto

2026-10-02. Worktree `pigauto-defaults-audit`, branch `docs/defaults-audit`, based on `origin/main` b565cad (PR #189 merged). No default is changed. The branch carries the documentation fixes, the brms guard and the tests described under Decision 4 and inconsistencies 5 and 6; every verdict below is a proposal and Shinichi decides. **Line numbers are origin/main b565cad** unless marked "(branch)"; on the branch `NEWS.md` shifts by about 25 lines and `R/pool_mi.R` by up to 13.

**Aim.** A user should be able to use pigauto's defaults well without knowing the alternatives. For each argument that offers alternatives, the tables give the current default (read from the function signature in `R/`, with file and line), the evidence, and a proposed verdict:

- **R**: recommended default.
- **A**: documented alternative, with caveat.
- **D**: candidate for deprecation (Shinichi decides).

Evidence paths are relative to the repository. "No benchmark found" means a search of `NEWS.md`, `docs/dev-log/` and `script/*.md` found none. Evidence on unmerged branches is cited as "branch X, not on main" and was read with `git show`, not checked out. All numbers state their regime; if a regime is missing from a quote it was missing from the source.

## Decisions and proposals

Decisions below were made by the maintainer (Shinichi, 2026-10-02). Items marked "proposal" are follow-up work he asked to have written down; none is implemented here.

1. **`draws_method` stays `"conformal"` for now (decision).**
   - Key caveat: conformal multiple-imputation (MI) *draws* bias downstream slopes and lose coverage under Rubin's rules. Branch `arc/mi-gls-attenuation`, not on main (`docs/dev-log/mi-gls/results.md`): across 16 phylogenetic generalised least squares (PGLS) regimes (two traits, bivariate Brownian motion with correlation 0.7, true slope 0.70, lambda 1 or 0.5, n 300 or 1000, MCAR or clade-structured missingness, 120 replicates, m = 20, 500 GNN epochs), conformal draws biased the pooled slope by -0.20 to -0.46 and the pooled 95% interval covered the truth in 0 to 17% of replicates; MC-dropout draws biased it by -0.03 to -0.38.
   - Per-cell conformal *intervals* are a different object and stay calibrated: coverage 0.954 to 0.971 in the gated MCAR rows of the PR #189 sweep (`docs/dev-log/mi-posterior/results.md`, Headline 1; `impute(gnn = FALSE)`, 20 regime x trait rows). Use them for per-cell uncertainty.
   - Already enforced in code: `with_imputations()` refuses `multi_impute()` output that is not posterior (`R/with_imputations.R:83-109`), and `pool_mi()` refuses legacy and unknown-provenance objects (`R/pool_mi.R:170-192`). A list of fits built by hand is still pooled, with only a provenance warning (`R/pool_mi.R:190`, `:403-407`), so a user who bypasses `with_imputations()` can still pool conformal draws.
   - `draws_method = "posterior"` (PR #189; continuous traits only, no GNN, no `species_col`, no `covariates`) fixes the inference problem. Paired bias against complete-data slopes is -0.014 to +0.009 in the in-model regimes (+0.017 is reached only in the stress regimes) (`docs/dev-log/mi-posterior/results.md`, Headline 2).
   - The cost of the refusal is a usability trap: a user who calls `multi_impute()` with its defaults and then follows the documented `with_imputations()` / `pool_mi()` path is stopped at the last step. The default route does not lead to a supported analysis. This is the main argument for changing the default.
   - Proposal: a follow-up PR adds `draws_method = "auto"`: posterior when every trait is continuous and there is no `species_col` or `covariates`; otherwise conformal plus a one-line message saying why the draws cannot be pooled.
2. **`gnn`: proposal to make `gnn = FALSE` the default.**
   - A one-line message when `covariates` are supplied: "covariates are used only by the GNN; set `gnn = TRUE`". `gnn = "auto"` was considered and rejected: it would hide a large runtime difference (2026-09-19 campaign, default arms, 4 threads: GNN on 118 s against GNN off 0.5 s at n = 100, 165 against 4.5 s at n = 300, 526 against 117 s at n = 1000) and a torch requirement, and there is no measured evidence that the GNN uses covariates well.
   - Rerun at current main (`docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`, branch; 200 cells, 20 seeds, same masks as 2026-09-19): GNN on as shipped is worse than GNN off in all ten DGP x n cells, by 12 to 23% in continuous z-RMSE on the BM, OU and AVONET300 data (AVONET300 +0.101, MCSE 0.011) and 1 to 2% on the low-signal DGP. GNN on predicting from the full baseline is within two MCSEs of GNN off everywhere (AVONET300 -0.018, MCSE 0.011). Fit time: GNN off 1.4 to 24 s, GNN on 122 to 496 s (n = 100 to 1000, 4 threads). The evidence supports the proposal; the flip remains Shinichi's decision.
   - Follow-ups: (a) the flip itself, with a NEWS entry, a minor version bump, and the line that restores the old behaviour (`gnn = TRUE`); (b) `gnn = TRUE` should predict from `baseline_full` and report both answers. Today only `gnn = FALSE` predicts from the refit without held-out cells (`R/fit_pigauto.R:551-567`; `R/predict_pigauto.R:349-357`). The 2026-09-19 campaign attributes the 15 to 25% worse z-RMSE of shipped GNN-on to this held-out-cell cost, not to the network (`docs/dev-log/arc/2026-09-19-campaign-gnn-off-results.md`, item 1). (c) Move torch from Imports to Suggests (`DESCRIPTION:37-38` lists it under Imports).
   - Side effect to decide: `conformal_method = "mondrian"` requires `gnn = TRUE` (`docs/dev-log/handover/2026-09-22-claude-handover-mondrian-conformal.md`, Key Decisions), so it would become opt-in on top of an opt-in.
3. **`joint_solver`: flagged.** The auditor's proposal, not a maintainer decision, is to treat this question before the GNN default.
   - On 2026-09-19 (`docs/dev-log/arc/2026-09-19-avonet-gap-results.md`, main ebbf63e) the in-house solver lost about 30% of continuous z-RMSE on AVONET300 against `joint_solver = "rphylopars"` (0.790 against 0.433).
   - Rerun at current main (`docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`, branch; GNN off, default safety machinery, 20 seeds, same masks): the in-house default is now 0.539 on AVONET300, level with raw Rphylopars (0.547), so most of that gap closed with PRs #187, #191 and #192. The Rphylopars solver is still lower: paired difference -0.098 (MCSE 0.025, -14%) on AVONET300 and -7% on the BM and OU simulations at n = 1000, about -1 to -2% at n = 300, and +6.5% worse on `ou_mixed` at n = 100. It costs 7 to 22 times the fit time (104 s against 4.8 s on AVONET300). Discrete accuracy and interval coverage do not differ.
   - Proposal: keep `"inhouse"` as the default for speed and self-containment; document `"rphylopars"` as the more accurate option for n of a few hundred or more when fit time is acceptable.
4. **`pool_mi`: brms fits are refused, with a pointer to `brms::brm_multiple()`** (posterior concatenation, not Rubin's rules). Added on this branch (`R/pool_mi.R:240-250`, branch); on main a `brmsfit` fell through to `coef()` / `vcov()`, and `coef()` on a brmsfit does not return its fixed effects. `brms`, `drmTMB` and `gllvmTMB` are added to Suggests on this branch; `pool_mi()` already had adapters for the last two (`R/pool_mi.R:412-474`). End-to-end tests through `multi_impute(draws_method = "posterior")`, `with_imputations()` and `pool_mi()` now cover glmmTMB, lme4, drmTMB and gllvmTMB (`tests/testthat/test-pool-mi-backends.R`).

5. **New candidate question (auditor's, not a decision): `safety_floor` and `phylo_signal_gate`.** The rerun found the pure baseline lower in error than the default safety machinery at small n and on low-signal simulated data, and 3 to 6 times faster, because lambda estimation now does the shrinkage these switches were added for. Simulations only; worth a targeted check on real data before any change.

## Small inconsistencies found

Verified by reading the code. Items 5 and 6 are fixed in the roxygen on this branch (NEWS entries are left as history); the rest are not changed here.

1. `epochs`: `impute()`, `multi_impute()` and `multi_impute_trees()` default to 2000 (`R/impute.R:386`, `R/multi_impute.R:380`, `R/multi_impute_trees.R:258`); `fit_pigauto()` defaults to 3000 (`R/fit_pigauto.R:372`).
2. `impute()` passes `phylo_signal_method = "lambda"` as a bare string (`R/impute.R:403`) with no `match.arg()`; only `fit_pigauto()` validates it (`R/fit_pigauto.R:393`, `:411`), after preprocessing and graph building. A typo errors late.
3. `multi_impute_trees()` defaults `draws_method` to `"mc_dropout"` and offers only `"mc_dropout"` and `"conformal"` (`R/multi_impute_trees.R:262`); `multi_impute()` defaults to `"conformal"` and also offers `"posterior"` (`R/multi_impute.R:372-373`). The two multiple-imputation entry points disagree on the default.
4. Roxygen of the internal `fit_joint_mvn_baseline()` and `fit_joint_threshold_baseline()` says `lambda_mode` is `"fixed_1"` "(default)" (`R/joint_mvn_baseline.R:35`, `R/joint_threshold_baseline.R:284`). Those functions' own signature defaults are indeed `"fixed_1"` (`:58`, `:309`), but `fit_baseline()` always passes the user's value, whose default is `"estimate"` (`R/fit_baseline.R:298`). The comment is misleading, not wrong about the signature.
5. `gate_method`: the signature default is `"cv_folds"` (`R/fit_pigauto.R:387`; introduced in commit 9d4782e, PR #102), but the roxygen at `R/fit_pigauto.R:170-171` and NEWS (`NEWS.md:1573`) still call `"single_split"` the default.
6. `min_val_cells`: signature default 20 (`R/fit_pigauto.R:394`); roxygen says "Default 10" (`R/fit_pigauto.R:217-220`).
7. `n_imputations` roxygen on `impute()` says "MC-dropout imputation sets" (`R/impute.R:37`); under `gnn = FALSE` the draws are baseline-posterior draws (`NEWS.md`, 0.11.0 section).
8. NEWS entries for `pool_method` (Phase H, `NEWS.md:1310`) and an older "auto" pooling note (`NEWS.md:927`) describe different arguments (`predict.pigauto_fit()` against a benchmark helper); the current enum is `"median"`, `"mean"`, `"mode"` with `"median"` the default (`R/predict_pigauto.R:214`).
9. `NEWS.md:2448-2457` still presents `share_gnn = TRUE` as a tree-uncertainty feature; the 0.11 compatibility section supersedes this (completions are descriptive sensitivity only). Readers of the older entry can be misled.
10. `docs/dev-log/handover/2026-09-22-claude-handover-mondrian-conformal.md` cites `docs/dev-log/2026-08-16-mechanism-coverage-results.md`, which is not on main (not found in `docs/dev-log/`).

## impute()

Signature: `R/impute.R:380-406`. Primary entry point; chains preprocessing, baseline, GNN, calibration and prediction.

| arg | default | alternatives | evidence | verdict | note |
|---|---|---|---|---|---|
| `gnn` (L405) | `TRUE` | `FALSE` | Rerun at current main (`docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`, branch, 200 cells, 20 seeds): GNN on as shipped is 12 to 23% worse than GNN off on BM, OU and AVONET300, 1 to 2% on the low-signal DGP; GNN on from the full baseline ties GNN off. Earlier: 2026-09-19 campaign, 200 cells, 20 seeds, 30% MCAR, n 100/300/1000: on `bm_mixed` and `ou_mixed` GNN off, raw Rphylopars and GNN on from the full baseline are within one MCSE; shipped GNN-on is 15 to 25% worse in z-RMSE. Wall time at n = 1000: 2.2 s (pure) to 117 s (default safety machinery) for GNN off, 526 s GNN on (`docs/dev-log/arc/2026-09-19-campaign-gnn-off-results.md`). On AVONET300 GNN on from the full baseline is 0.723 against 0.790 GNN off (same note, table). Rerun at current main: results pending, `docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`. | A now; proposal in Decision 2 is to make FALSE the default | The 2026-09-19 AVONET300 gain (0.723 against 0.790) did not survive the new baseline defaults: at current main the full-baseline GNN arm is 0.520 against 0.539, 1.6 MCSEs. Covariates act only through the GNN (`R/fit_pigauto.R:555-558`). |
| `lambda_mode` (L388) | `"estimate"` | `"fixed_1"`, `"cv"`, `"bayes"` | 13 real-data cases (continuous traits, 5 seeds, GNN off, Totoro): estimate helps or ties in 12, never worse by more than 1%; PanTHERIA -14.8% (n 300) (`docs/dev-log/lambda-default/real-data.md`). 18 core simulation cells, 200 seeds: z-RMSE falls in every cell, 0.021 to 0.086 (`docs/dev-log/lambda-default/benchmark.md`). `cv` and `bayes` cost 10 to 50% more for little extra gain (`real-data.md`). Full-REML fix: bias -0.061 to -0.022 at true lambda 0.3 (n 300, 30% MCAR, 200 seeds) (`docs/dev-log/2026-09-24-lambda-reml-report.md`). | R for `"estimate"`; A for `"cv"`, `"bayes"`; A for `"fixed_1"` | Pre-registered gate G12 was missed at lambda 0.3, n 100 (17% of the gap to a freq-lambda reference closed against 50% asked); Shinichi kept the default on 2026-09-23. A known case against it: LepTraits flight duration, n 2000, 1.000 against 0.968. `"cv"`/`"bayes"` are per-column only and reroute the pipeline. |
| `joint_solver` (L389) | `"inhouse"` | `"rphylopars"` | Current main, GNN off, 20 seeds: AVONET300 in-house 0.539, Rphylopars 0.440 (paired -14%); simulations at n = 1000 -7%, at n = 100 tie or worse; 7 to 22 times slower (`docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`, branch). Earlier gap 0.790 against 0.433 (`docs/dev-log/arc/2026-09-19-avonet-gap-results.md`). Under `model = "lambda"` Rphylopars is 26 times slower, with 1 to 8 explosive seeds per cell before the plausibility guard (`docs/dev-log/lambda-default/benchmark.md`). | R for `"inhouse"`; A for `"rphylopars"` (Decision 3) | The 30% gap of 2026-09-19 is no longer current. Rphylopars is a Suggests dependency; the fallback to in-house exists. |
| `predict_method` (L390) | `"auto"` | `"exact"`, `"per_column"` | 18 core cells, 200 seeds, GNN off, two cross-trait correlations: z-RMSE of continuous-family traits falls 3.1 to 6.6% at lambda 0.3, 6.1 to 9.4% at 0.7, 1.2% (n 100) and 0.3% (n 300) at lambda 1, and rises 0.4% at n 1000. Discrete accuracy +0.009 to +0.032 at 0.3 and 0.7. 13 real cases, 5 seeds, no MC error: -17 to -51% on AVONET and PanTHERIA, -0.9 to +0.8% on seven others (`NEWS.md`, 0.11.0.9000; `docs/dev-log/exact-default/after-task.md`). | R | Coverage at lambda 1, n = 100 falls a further 0.006 to 0.018. Not re-benchmarked with the GNN on. `"exact"` alone was rejected (+2.0 to +2.6% z-RMSE at lambda 1). |
| `multi_obs_aggregation` (L387) | `"hard"` | `"soft"` | No benchmark at current main. A small script run (150 species, 2 reps, 100 epochs, 2026-04-18) gives soft minus hard accuracy of +0.016, +0.014, +0.024 for the binary trait in three scenarios and -0.016 to +0.003 for the three-level trait (`script/bench_multi_obs_mixed.md`). NEWS prose: binary +1.4 to 2.4 pp, categorical mixed (`NEWS.md:2495-2515`). | A | Default kept `"hard"` for back-compatibility only. Predates the lambda, exact and GNN-off work. Evidence is two replicates. |
| `phylo_signal_gate` (L401), `phylo_signal_threshold` (L402), `phylo_signal_method` (L403) | `TRUE`, `0.2`, `"lambda"` | `FALSE`; other thresholds; `"blomberg_k"` | In the 2026-09-19 campaign the default GNN-off arm (gate and safety floor on) has the lowest z-RMSE of the pigauto arms on low-signal `bace_dgp` (0.956 against 1.034 pure, n 1000, 20 seeds, mean-floor 1.0). It costs about 0.02 at n = 100 on `bm_mixed` (0.270 against 0.251) and `ou_mixed` (0.431 against 0.409), and is identical at n >= 300. The gate fires on raw trait values, so it blocks covariate lift for weak-signal traits (`NEWS.md:1924-1957`). No benchmark found for the threshold 0.2 itself or for `"blomberg_k"`. | R for `TRUE`/`"lambda"`; A for `FALSE` (covariates, or a pure-BM comparison); `"blomberg_k"`: D | `"blomberg_k"` is "not dimensionally comparable" to the threshold (`R/fit_pigauto.R:212-216`) and is not validated by `impute()` (inconsistency 2). |
| `safety_floor` (L400) | `TRUE` | `FALSE` | Floor guarantees validation RMSE no worse than the grand mean by construction; plants, weak-signal traits fixed (`NEWS.md:2050-2085`). In the campaign the safety machinery is what separates default from pure GNN-off on low-signal data (same low-signal numbers as above). Cost: `calibrate_gates()` is 117 s at n = 1000 against 2.2 s without it (campaign, wall-time table). | R, to be re-examined | At current main the pure arm is lower than the default GNN-off arm on the low-signal DGP (-0.8 to -2.7%) and at n = 100 (-2.7 to -4.5%), identical at n = 1000, and 3 to 6 times faster (`docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`, item 3): lambda estimation now does the shrinkage the safety machinery used to provide. Simulations only. The pure traditional-statistics arm is `gnn = FALSE, safety_floor = FALSE, phylo_signal_gate = FALSE`. The cost at n >= 1000 is a usability item, not an accuracy item. |
| `conformal_split_val` (L404) | `FALSE` | `TRUE` | Opt-in; forcing the split regressed AVONET300, OVR categorical and BIEN smoke benches by 2 to 26% on small validation sets (`NEWS.md:1802-1818`). The `"auto"` route-choice split already separates route choice from calibration for traits with >= 38 held-out rows (`NEWS.md`, 0.11.0.9000). No benchmark found at current main. | A | Use `TRUE` when accurate 95% intervals matter more than point error, per NEWS. |
| `epochs` (L386) | `2000L` | any | No benchmark found for the value. Early stopping (`patience = 10`, `eval_every = 100`, `R/fit_pigauto.R:381-382`) usually ends training earlier. Irrelevant under `gnn = FALSE`. | R, provisional | Differs from `fit_pigauto()` (inconsistency 1). |
| `missing_frac` (L384) | `0.25` | `0` to skip splitting | Held-out fraction for validation and testing, not the missingness of the data. At n = 100 a trait has about 21 validation cells and production-interval coverage is 0.89 to 0.93; at n >= 300 it is 0.96 to 0.98 (campaign, item 4). | R | `0` disables evaluation and gate calibration. `min_val_cells` (default 20) warns when too few cells remain. |
| `pool_method` (L395) | `"median"` | `"mean"`, `"mode"` | Ordinal, AVONET Migration, 3 seeds, n = 1500, 30% held out, N_IMP = 20: median 0.713, mode 0.779, class-mode baseline 0.800 (`NEWS.md:1310-1360`). Applies only when `n_imputations > 1`. | R for `n_imputations = 1`; A (`"mode"`) for ordinal with `n_imputations > 1` | NEWS itself recommends `"mode"` for ordinal with many draws; the default was kept for back-compatibility. Not re-measured at current main. |
| `n_imputations` (L384) | `1L` | integers > 1 | These are stochastic prediction draws (baseline posterior plus dropout), not analysis-aware MI (`R/predict_pigauto.R:82-85`). | R | For downstream inference use `multi_impute(draws_method = "posterior")` or `multi_impute_analysis()`. |
| `clamp_outliers` (L396), `clamp_factor` (L397) | `FALSE`, `5` | `TRUE`; other factors | AVONET Mass, seed 2030: 24,330 to 6,273 (-74%); seed 2031: 312 to 318 (`NEWS.md:1239-1275`, N_IMP = 20). Two seeds only. | A | Affects log-transformed continuous, count and zero-inflated magnitude traits only. |
| `match_observed` (L398), `pmm_K` (L399) | `"none"`, `5L` | `"pmm"` | PMM "is not a tail-safety tool for pigauto" and failed its redesign pilot; it is an experimental stochastic mechanism (`NEWS.md:1112-1140`, `NEWS.md:586`). | A for `"none"`; `"pmm"`: D | Candidate for deprecation: no supported use case remains in the docs. |
| `log_transform` (L383) | `TRUE` | `FALSE` | Auto-logs strictly positive continuous traits. No benchmark found comparing the two. On AVONET300 Mass the MCSE is 0.24 to 0.35 and the log transform is one of three unseparated suspects for the real-data gap (campaign, "What this does not say"). The Rubin campaign on branch `arc/rubin-freq-bace` (not on main) ran posterior MI with `log_transform = FALSE`. | R, provisional | `draws_method = "posterior"` draws live on the log scale; analysing a logged trait on the raw scale is outside the supported regime (`NEWS.md`, 0.11.0.9000). |
| `joint_refine_iter` (L391), `em_iterations` (L392), `em_offdiag` (L394), `em_tol` (L393) | `0L`, `0L`, `FALSE`, `1e-3` | positive integers; `TRUE` | EM refinement was disabled on 2026-05-17 after divergence (`docs/dev-log/arc/2026-09-19-avonet-gap-results.md`, "What this means"). `joint_refine_iter` has no effect on traits routed through `"exact"` (`R/impute.R:143-144`). No benchmark found. | A (off); candidates for deprecation if the solver question is settled in favour of Rphylopars | Opt-in controls with a rollback guard; left off. |
| `species_col`, `trait_types`, `multi_proportion_groups`, `covariates`, `seed`, `verbose` (L380-386) | `NULL` or `TRUE` | n/a | Structural arguments, not tuning choices. | R | `covariates` are ignored with a warning under `gnn = FALSE`. |

## multi_impute()

Signature: `R/multi_impute.R:371-384`. Produces `m` completed datasets.

| arg | default | alternatives | evidence | verdict | note |
|---|---|---|---|---|---|
| `draws_method` (L372-373) | `"conformal"` | `"mc_dropout"`, `"posterior"` | Conformal draws: pooled PGLS slope bias -0.20 to -0.46, coverage 0 to 17% in 16 regimes (branch `arc/mi-gls-attenuation`, not on main, `docs/dev-log/mi-gls/results.md`, Table 1; 120 replicates, m = 20). MC-dropout: -0.03 to -0.38 (same table). Posterior: paired bias -0.014 to +0.009 (gls, in-model regimes) and -0.026 to +0.017 across gls and phylolm in the stress regimes; pooled 95% interval coverage 0.845 to 0.975 against 0.76 to 0.97 for complete data (the 48 in-model rows, regimes 17 to 40, 200 reps, n 300/1000, 30% missing; `docs/dev-log/mi-posterior/results.md`). Branch `arc/rubin-freq-bace`, not on main (`docs/dev-log/arc/2026-10-02-rubin-pigauto-campaign.md`, 3,600 fits, 200 datasets per cell, `m = 20`): posterior MI matches proper frequentist MI except n = 1000, lambda = 1, where slope coverage is 0.860 against 0.963 (paired difference -0.102, SE 0.016); 12 sampler errors (Cholesky) and 56 + 51 unconverged fits at n = 100. | Decision 1: stay `"conformal"` for now; proposal `"auto"` | Conformal: R for per-cell uncertainty, not for downstream inference. `"mc_dropout"`: A (draws from the baseline posterior when `gnn = FALSE`). `"posterior"`: A with caveats (continuous traits only; small negative bias at lambda = 1; failures at n = 100). | 
| `m` (L371) | `100L` | integers >= 2 | Posterior campaign used `m = 20` (`docs/dev-log/mi-posterior/design.md:218`). No benchmark found for the default of 100. | R, provisional | 100 is large for the posterior route (n = 1000 took about 1,800 s per fit, single core, 36-fit pre-run on the arc/rubin branch; the cost is mostly sampler length, not `m`). |
| `gnn` (L381), `epochs` (L380), `missing_frac` (L378), `lambda_mode` (L382) | `TRUE`, `2000L`, `0.25`, `"estimate"` | as for `impute()` | See `impute()`. Under `"posterior"` all four are ignored and a message lists any that were supplied (`R/multi_impute.R:388-391`). | as for `impute()` | `gnn` flip would change which draws are produced under the default `"conformal"` route. |
| `posterior_control` (L383) | `list()` | chains, burn-in, `n_iter`, `thin`, `keep_draws`, `param_uncertainty`, `seed`, `auto_extend`, `max_extend` | Convergence rule: split R-hat < 1.05 and bulk ESS > 400, up to 3 extensions; 7,999 of 8,000 simulation fits converged (`docs/dev-log/mi-posterior/results.md`). `param_uncertainty = "none"` is for validation only and is refused downstream. | R for defaults | Do not change `param_uncertainty` for analysis. |
| `log_transform` (L377) | `TRUE` | `FALSE` | See `impute()`. | R, provisional | Posterior MI is linear-Gaussian on the logged scale. |

## multi_impute_trees()

Signature: `R/multi_impute_trees.R:251-264`. Posterior-tree sensitivity completions. Not a supported downstream-inference route (`NEWS.md`, 0.11.0 compatibility section; `R/with_imputations.R:92-98` refuses its output).

| arg | default | alternatives | evidence | verdict | note |
|---|---|---|---|---|---|
| `share_gnn` (L260) | `TRUE` | `FALSE` | Documented as a 10 to 15 times speed-up at n = 10,000, not measured in this audit (`R/multi_impute_trees.R:83-88`). No benchmark found comparing the two on accuracy or on tree sensitivity. Under `gnn = FALSE` the argument only changes whether the baseline is refit per tree without held-out cells. | R, provisional | Only meaningful with `gnn = TRUE`. Per-tree replay carries `joint_solver`, `predict_method` and `joint_refine_iter` since 0.11 (`docs/dev-log/after-task/2026-09-18-gnn-off.md`, section 8). |
| `draws_method` (L262) | `"mc_dropout"` | `"conformal"` | No `"posterior"` option (inconsistency 3). Same attenuation evidence as `multi_impute()`. For tree-sensitivity use the draw method hardly matters; the tree is the quantity of interest. | A | Default disagrees with `multi_impute()`. Candidate to align once Decision 1 is settled. |
| `gnn` (L263) | `TRUE` | `FALSE` | See `impute()`. | as `impute()` | Under `gnn = FALSE` `"mc_dropout"` means baseline-posterior draws (one message). |
| `m_per_tree` (L251) | `1L` | integers >= 1 | No benchmark found. | R | Does not create an inference-supported workflow (`NEWS.md`, older entries marked superseded). |
| `reference_tree` (L261) | `NULL` (maximum clade credibility tree) | any `phylo` | No benchmark found. | R | Requires `phangorn`; otherwise falls back to `trees[[1]]` with a warning. |
| `epochs` (L258), `missing_frac` (L256), `log_transform` (L255) | `2000L`, `0.25`, `TRUE` | as for `impute()` | See `impute()`. | as `impute()` | |

## fit_pigauto()

Signature: `R/fit_pigauto.R:350-401`. Pipeline-level entry. Only arguments with alternatives that a user might change are tabled; hyper-parameters (`hidden_dim`, `lr`, `dropout`, `corruption_*`, `refine_steps`, `lambda_shrink`, `lambda_gate`, `warmup_epochs`, `edge_dropout`, `eval_every`, `patience`, `clip_norm`) have no recorded benchmark in the dev-log and are left at their defaults (verdict R, "no benchmark found").

| arg | default | alternatives | evidence | verdict | note |
|---|---|---|---|---|---|
| `gnn` (L356) | `TRUE` | `FALSE` | See `impute()`. `baseline_full` (L357) is built only when `gnn = FALSE` (`R/fit_pigauto.R:559-570`). | A now; see Decision 2 | Follow-up (b) in Decision 2 changes this. |
| `use_transformer_blocks` (L363), `n_heads` (L364), `ffn_mult` (L365) | `TRUE`, `4L`, `4L` | `FALSE` (legacy attention) | Small script run, 200 species, 2 reps, 100 epochs, 2026-04-16: identical to legacy in high-signal scenarios; RMSE 0.8349 against 0.8539 (low signal) and 0.7174 against 0.7509 (moderate signal); about 25 to 30% slower (`script/bench_discriminative_phase9.md`; `NEWS.md:2635-2655`). No benchmark found at current main or at n >= 300. | A (transformer is a documented default with two replicates behind it) | Irrelevant when the gate is closed and under `gnn = FALSE`. |
| `use_attention` (L362) | `TRUE` | `FALSE` | No benchmark found. Applies only to the legacy path (`use_transformer_blocks = FALSE`). | D | Candidate for deprecation: it matters only for a legacy architecture that is itself an opt-out. |
| `use_trait_attention` (L366), `n_trait_heads` (L367), `trait_embed_dim` (L368) | `FALSE`, `2L`, `32L` | `TRUE` | Negative result: WorldClim plus `use_trait_attention = TRUE` was +0.11 worse in the W-series (`NEWS.md:780-797`). | A | One data set (BIEN plants); candidate for deprecation if never re-benchmarked. |
| `conformal_method` (L384) | `"split"` | `"bootstrap"`, `"mondrian"` | Mondrian addresses undercoverage under clade-structured missingness (-3.4 pp for split at n 300/1000); verified on simulated mechanisms only, never re-run on the real data where undercoverage appeared (fishbase 0.89 to 0.91; pantheria 0.87 to 0.94) (`docs/dev-log/handover/2026-09-22-claude-handover-mondrian-conformal.md`; `NEWS.md:391-410`). Needs about 38 validation cells per trait, single-observation data and `gnn = TRUE`. No benchmark found for `"bootstrap"`. | R for `"split"`; A for `"mondrian"`; `"bootstrap"`: D | Roxygen: these "support nominal held-out diagnostics, not package-certified coverage" (`R/fit_pigauto.R:150-159`). |
| `conformal_split_val` (L386) | `FALSE` | `TRUE` | See `impute()`. | A | |
| `gate_method` (L387), `gate_cv_folds` (L389), `gate_splits_B` (L388) | `"cv_folds"`, `5L`, `31L` | `"median_splits"`, `"single_split"` | Smoke on AVONET300 (60 epochs, 60 Mass + 60 Migration cells masked, 1 replicate): 82.6 s single, 77.2 s median, 75.9 s cv_folds. Synthetic mixed-type (n = 300, 5 reps): continuous RMSE 1.0778 (single), 1.0720 (median), 1.0611 (cv_folds) (`NEWS.md:1527-1620`). | R for `"cv_folds"` | Docs still say `"single_split"` is the default (inconsistency 5). Five replicates; not re-run at current main. |
| `safety_floor` (L390), `phylo_signal_gate` (L391), `phylo_signal_threshold` (L392), `phylo_signal_method` (L393) | `TRUE`, `TRUE`, `0.2`, `"lambda"` | see `impute()` | See `impute()`. | see `impute()` | |
| `min_val_cells` (L394) | `20L` | any | 19 cells is the smallest set whose conformal ceiling n/(n+1) reaches 0.95 (`NEWS.md:391-410`; `NEWS.md`, silent-fallback section). | R | Docs say 10 (inconsistency 6). |
| `lambda_mode` (L395), `joint_solver` (L396), `predict_method` (L397), `joint_refine_iter` (L398) | `"estimate"`, `"inhouse"`, `"auto"`, `0L` | see `impute()` | See `impute()`. | see `impute()` | |
| `epochs` (L372) | `3000L` | any | No benchmark found. | R, provisional | Differs from `impute()` (inconsistency 1). |
| `gate_cap` (L361) | `0.8` | 0 to 1 | "Safety comes from regularisation, not the cap" (`R/fit_pigauto.R:91-92`). No benchmark found. | R | |

## fit_baseline()

Signature: `R/fit_baseline.R:295-307`. The phylogenetic baseline alone.

| arg | default | alternatives | evidence | verdict | note |
|---|---|---|---|---|---|
| `lambda_mode` (L298) | `"estimate"` | `"fixed_1"`, `"cv"`, `"bayes"` | See `impute()`. Per-type dispatch: continuous-family columns estimate their own lambda; binary, categorical, zero-inflated gate and ordinal columns stay at lambda = 1 (`NEWS.md`, "Default flip"). | R | Covariate-aware baselines estimate lambda; `"cv"` and `"bayes"` fall back to lambda = 1 there, with a warning. |
| `lambda_fixed` (L299) | `NULL` | named numeric vector | Rebuilds a baseline at stored lambdas; reproduces mu and se within 1e-10 over 12 seeds (`docs/dev-log/exact-default/after-task.md`, section 6). | R | Mechanism, not a tuning choice. |
| `joint_solver` (L303) | `"inhouse"` | `"rphylopars"` | See `impute()` and Decision 3. | R / A, see `impute()` | |
| `predict_method` (L304), `predict_route` (L306) | `"auto"`, `NULL` | `"exact"`, `"per_column"`; named route | See `impute()`. `predict_route` takes precedence over `predict_method` and replays a stored choice. | R | A trait with fewer than 5 validation cells, or no split, uses `"exact"`. `"exact"` falls back to per-column above about 20,000 unknown cells. |
| `multi_obs_aggregation` (L297) | `"hard"` | `"soft"` | See `impute()`. | A | |
| `model` (L295) | `"BM"` | other strings | Only `"BM"` is the documented path; no benchmark found for any alternative. | R | Candidate for a docs check: confirm what other values do or remove the argument. |
| `em_iterations` (L300), `em_tol` (L301), `em_offdiag` (L302), `joint_refine_iter` (L305) | `0L`, `1e-3`, `FALSE`, `0L` | positive integers; `TRUE` | See `impute()`. Fixed in 0.11.0.9000: `em_iterations >= 1` ignored `predict_method` inside the EM loops. | A (off) | |

## predict.pigauto_fit()

Signature: `R/predict_pigauto.R:212-219`.

| arg | default | alternatives | evidence | verdict | note |
|---|---|---|---|---|---|
| `n_imputations` (L213) | `1L` | integers > 1 | Draws are stochastic prediction variation, "not validated analysis-aware multiple imputations" (`R/predict_pigauto.R:82-85`). | R | |
| `pool_method` (L214) | `"median"` | `"mean"`, `"mode"` | See `impute()`. | R for 1 draw; A (`"mode"`) for ordinal with many draws | |
| `clamp_outliers` (L215), `clamp_factor` (L216) | `FALSE`, `5` | `TRUE` | See `impute()`. | A | |
| `match_observed` (L217), `pmm_K` (L218) | `"none"`, `5L` | `"pmm"` | See `impute()`. | `"pmm"`: D | |
| `return_se` (L212) | `TRUE` | `FALSE` | Baseline SE is exact under Brownian motion and model-dependent otherwise; discrete-trait values are uncertainty scores, not standard errors (`CLAUDE.md`, uncertainty section). The per-cell conformal interval is the primary 95% interval. | R | Do not feed discrete `se` into Rubin arithmetic. |
| `baseline_override` (L213) | `NULL` | `list(mu, se)` | Internal, used by `multi_impute_trees()`. | R | |
| `newdata` (L212) | `NULL` | data.frame | Not benchmarked here. | R | |

## with_imputations()

Signature: `R/with_imputations.R:72-74`. Fits one analysis per completed dataset.

| arg | default | alternatives | evidence | verdict | note |
|---|---|---|---|---|---|
| `.on_error` (L74) | `"continue"` | `"stop"` | Behavioural, not a statistical choice: with `"continue"` failures are captured, a warning summarises them, and `pool_mi()` drops error elements (`R/with_imputations.R:51-57`). No benchmark found. | R | With `"continue"`, check the warning: pooling over fewer than `m` fits silently changes the Rubin variance denominator. |
| `.progress` (L73) | `interactive()` | `TRUE`, `FALSE` | None needed. | R | |
| `mi` (input) | none | `pigauto_analysis_mi`, or posterior `pigauto_mi` | Refuses conformal and MC-dropout `multi_impute()` output, tree-sensitivity output, plug-in posterior draws and bare dataset lists (`R/with_imputations.R:83-109`). | R | The refusal is the practical form of Decision 1's caveat. |

## pool_mi()

Signature: `R/pool_mi.R:131-136`. Rubin's rules with optional Barnard-Rubin degrees of freedom.

| arg | default | alternatives | evidence | verdict | note |
|---|---|---|---|---|---|
| `conf.level` (L132) | `0.95` | (0, 1) | Posterior-MI pooled 95% coverage 0.845 to 0.975 against 0.76 to 0.97 complete-data, within 0.05 in 48 of 48 rows (`docs/dev-log/mi-posterior/results.md`). | R | |
| `df_fun` (L135) | `NULL` | function returning complete-data df | With a finite complete-data df the Barnard-Rubin (1999) correction applies; otherwise the large-sample df (`R/pool_mi.R:73-80`). No benchmark found for the effect of supplying it. | R | Supply it for small n. |
| `coef_fun` (L133), `vcov_fun` (L134), `tidy_fun` (L136) | `NULL` | functions | Automatic adapters cover `lm`, `glm`, `gls`, `lme`, `merMod`, `glmmTMB` (conditional), `drmTMB` and `gllvmTMB_multi` fixed effects (`R/pool_mi.R:18-19`, `:90-96`). A 6,000-task campaign passed 24 of 24 predeclared cells for the analysis-aware backend, within its narrow scope (`NEWS.md:575-583`). | R | Extractability does not mean the model has passed a validation gate (`R/with_imputations.R:38-41`). |
| `fits` (input) | none | list from `with_imputations()` | `brmsfit` objects are refused with a pointer to `brms::brm_multiple()` and `brms::combine_models()` (new on this branch, `R/pool_mi.R:240-250`, branch); MCMCglmm was already refused (`R/pool_mi.R:218-229`). Worked examples for each backend: `vignettes/multiple-imputation.Rmd`. | R | Decision 4. |

## suggest_next_observation()

Signature: `R/active_impute.R:417-423`.

| arg | default | alternatives | evidence | verdict | note |
|---|---|---|---|---|---|
| `by` (L418) | `"cell"` | `"species"` | No benchmark found. `"species"` sums available reductions across the species' missing traits. | R | Mixing variance and entropy in the species sum has no stated common scale. |
| `types` (L419-423) | all eight types | subsets | No benchmark found. | R | |
| `top_n` (L417) | `10L` | any | None needed. | R | |
| (method) | model-based proxy | none | Ranks cells by modelled variance or entropy reduction; the 0.11.0 NEWS describes it as a model-based proxy ranking under assumptions, with no sampling-design guarantee and no demonstrated field gain. It scores with the raw tree correlation (lambda = 1) and does not read the fit's `lambda_per_trait` (`NEWS.md`, "Known limitation"). | A | Under the new default `lambda_mode = "estimate"` the helper and the baseline can disagree on the correlation structure. Not benchmarked. |

## What this audit does not cover

- Whether any proposed default change improves outcomes for users. Decisions 2 and 3 rest on a rerun that is on this branch, not on main; it covers single-observation MCAR data only.
- The GNN-on arms under the current baseline defaults: the `predict_method = "auto"` benchmark used GNN off, and the lambda benchmark ran GNN on only for 100 seeds at n <= 300 and 50 at n = 1000.
- Multi-observation, `multi_proportion`, `zi_count`, covariate and missing-not-at-random regimes. Most benchmarks cited are single-observation MCAR.
- Hyper-parameters of the GNN (no benchmark found) and the arguments of `multi_impute_analysis()`, `cross_validate()`, `simulate_benchmark()` and `compare_methods()`.
- Branch-only evidence (`arc/mi-gls-attenuation`, `arc/rubin-freq-bace`) was read but not re-run; its numbers are as reported there.
