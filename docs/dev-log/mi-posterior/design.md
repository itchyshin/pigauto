# Posterior multiple imputation: design and symbolic alignment

2026-09-24. Branch `arc/mi-posterior` (stacked on PR #187). Plan:
`~/.claude/plans/read-agents-md-and-docs-dev-log-handover-curious-eagle.md` (approved 2026-09-24).
Locked by Shinichi: full Sigma_E; GNN ignored; weak inverse-Wishart priors with parameter expansion
for Sigma_P; option name `draws_method = "posterior"`.

## 1. Model

Continuous traits only, on pigauto's latent scale (the z-scored, possibly log-transformed columns of
`preprocess_traits()$X_scaled`). Single observation per species. For tip i (n tips) and trait k
(K traits):

    y_ik = mu_k + a_ik + e_ik

- a (all N nodes of the tree, tips and internal, by K traits): vec(a) ~ N(0, Sigma_P %x% A), where A is
  the phylogenetic covariance over nodes built so that its tip block equals pigauto's correlation
  matrix R = cov2cor(vcv(tree)). Use `build_henderson_S_inv()` (R/henderson_s_inv.R): its sparse
  precision Q = A^{-1} over the extended tree with the `tip_sqrt_d` correlation scaling, as
  `draw_conditional_bm()` already does. Never form A or R densely in the sampler.
- e_i ~ N(0, Sigma_E) independently across tips, with a full K x K Sigma_E.
- mu (K) has a flat prior.
- Implied marginal: vec(Y) ~ N(mu, Sigma_P %x% R + Sigma_E %x% I_n).
- Implied per-trait Pagel lambda: lambda_k = Sigma_P[k,k] / (Sigma_P[k,k] + Sigma_E[k,k]). Implied
  phylogenetic and residual correlations: cov2cor(Sigma_P), cov2cor(Sigma_E).
- Reduces to #187's single-lambda Kronecker model when Sigma_P = lambda * Sigma and
  Sigma_E = (1 - lambda) * Sigma, and to the prototype when Sigma_E = 0.

Missing cells y_mis are treated as unknowns (data augmentation).

## 2. Gibbs sampler (one sweep)

Revised 2026-09-24 after the design review (`review-design.md`: B1, B2, B3 and R1 to R6). An
orchestrator finding, recorded here as B4, is also folded in.

### 2.0 The one precision matrix (review B1)

`build_henderson_S_inv(tree)` returns Q, the precision of the raw node covariance: tips first, then
the non-root internal nodes, with the root excluded. The root state is therefore fixed at 0, and mu is
the root state, so mu is identified under a flat prior. The tip-block Schur complement of Q is
`vcv(tree)^{-1}`, not `R^{-1}`.

Define `D = diag(c(tip_sqrt_d, rep(1, n_internal_nonroot)))` and `Qc = D Q D`. D is diagonal, so Qc
has Q's sparsity pattern exactly. Its tip-block Schur complement is `D_T vcv^{-1} D_T = R^{-1}`, so
a ~ N(0, Sigma_P %x% Qc^{-1}) has tip marginal Sigma_P %x% R. Internal nodes sit on the rescaled
(irrelevant) scale, and tip values of a are in the units of y.

Qc is built once per tree. It is the only node precision used anywhere in the sampler: the (a, mu)
draw, the sufficient statistic, and the tests. Bare Q never appears.

### 2.1 State and priors (review B2)

The state holds the working effects xi (N x K), the working covariance Sigma_W (K x K), the
expansion parameters alpha (K), Sigma_E (K x K), mu (K), and y_mis.

Real parameters: a = xi diag(alpha) (row i is diag(alpha) xi_i) and
Sigma_P = diag(alpha) Sigma_W diag(alpha).

Priors, following MCMCglmm's parameter-expanded G prior:
- Sigma_W ~ IW(nu_W = K + 1, S_W = I_K).
- alpha ~ N(0, V_alpha I_K) with V_alpha = 1000.

Together these imply a scaled-F prior on each phylogenetic SD, flat near 0, and a uniform marginal
prior on the phylogenetic correlations. Other priors:
- Sigma_E ~ IW(nu_E = K + 1, S_E = 0.01 x diag of the observed-cell variances on the latent scale).
  Sigma_E is not expanded, as in MCMCglmm's R structure.
- mu is flat.

### 2.2 One sweep

1. **(a, mu, y_mis) | y_obs, Sigma_P, Sigma_E, drawn as one block (B4).** Alternating a and y_mis
   mixes badly when Sigma_E is small (lambda near 1), because y_mis nearly equals a_tip. So y_mis is
   integrated out of the first draw.
   - (a) Draw (a, mu) | y_obs. The joint precision over (vec(a) node-major, mu) is the sum of two
     parts:
     - the prior, `Qc %x% Sigma_P^{-1}` in node-major order (K x K block per node pair);
     - the likelihood, summed over tips i with observed set O_i, where
       `P_i = E_{O_i} (Sigma_E[O_i, O_i])^{-1} E_{O_i}'`. P_i adds to the (a_i, a_i), (a_i, mu),
       (mu, a_i) and (mu, mu) blocks.

     The linear term is the sum over tips of `P_i y_i`, with zeros in the unobserved slots. A tip with
     no observed trait adds nothing. The sparsity pattern is fixed across sweeps, so the symbolic
     Cholesky factorisation is done once and each sweep refactors numerically
     (`Matrix::Cholesky` + `update`). Draw by sample-from-precision: solve `L' z = N(0, I)` and add
     the mean.
   - (b) Draw y_mis | a, mu, y_obs, Sigma_E per tip. The residual e_i = y_i - mu - a_i is
     N(0, Sigma_E). Condition its missing components on its observed ones; a tip with nothing
     observed draws from N(0, Sigma_E).

   With the parameters fixed, (a) then (b) is an exact independent draw of (a, mu, y_mis) from its
   full conditional. That is what makes G2 an exactness test.
2. **Working effects.** xi = a diag(1 / alpha).
3. **alpha | xi, mu, y (completed), Sigma_E.** At tips, y_i - mu = diag(xi_i) alpha + e_i with
   e_i ~ N(0, Sigma_E); internal nodes carry no likelihood. This is Gaussian and conjugate:
   - precision `Pa = (Sigma_E^{-1} * (sum_i xi_i xi_i')) + I / V_alpha`, where * is the elementwise
     product and the sum runs over tips;
   - mean `Pa^{-1} sum_i diag(xi_i) Sigma_E^{-1} (y_i - mu)`.

   Then a = xi diag(alpha).
4. **Sigma_W | xi.** IW(nu_W + N, S_W + xi' Qc xi), using the WORKING xi over all N nodes, never the
   real a. Then Sigma_P = diag(alpha) Sigma_W diag(alpha).
5. **Sigma_E | e.** IW(nu_E + n, S_E + e'e), with e = Y_completed - 1 mu' - a_tip.

Reported per kept sweep: Sigma_P, Sigma_E, lambda_k = Sigma_P[k,k] / (Sigma_P[k,k] + Sigma_E[k,k]),
mu, and y_mis. The correlations cov2cor(Sigma_P) and cov2cor(Sigma_E) are not stored; compute them
per draw from the stored `params$Sigma_P` and `params$Sigma_E` arrays (K x K x draws), for example
`array(apply(Sigma_P, 3, cov2cor), dim(Sigma_P))`. alpha and Sigma_W are not identified and are
neither reported nor diagnosed.

### 2.3 Chains and diagnostics

- 4 chains by default, with dispersed starts. One chain starts at #187's REML lambda split; the
  others start from that point perturbed.
- Burn-in and thinning are set from the smoke.
- Split R-hat and bulk ESS (rank-normalised; Vehtari et al. 2021) are computed internally, with no new
  dependency, for every element of Sigma_P (upper triangle), Sigma_E (upper triangle) and lambda_k.
- `converged` = TRUE when every max R-hat < 1.05 and every min ESS > 400.

### 2.4 Lambda cross-check (review R4, verified)

#187's per-column model (`R/joint_mvn_solver.R`, `.mvn_gls_mean_at_lambda`) uses
`R_oo(lambda) = lambda R + (1 - lambda) I` on the correlation-scale R. Here the marginal for trait k
is Sigma_P[k,k] R + Sigma_E[k,k] I = (Sigma_P[k,k] + Sigma_E[k,k]) (lambda_k R + (1 - lambda_k) I),
so the two lambdas are the same quantity.

### 2.5 Changes made during the build (S1, 2026-09-24)

1. **Collapsed Metropolis moves before step 1.** With the section 2.2 sweep alone, mixing near the
   boundaries was too slow to use: at n = 1000, bulk ESS for lambda was 20 to 41 from 4 x 10,000
   sweeps at truth lambda = 1 and 0.95. So Metropolis moves run on the expanded state with
   (a, mu, y_mis) integrated out. Per trait k, three multiplicative moves:
   - alpha_k scaled by c;
   - row and column k of Sigma_E scaled by d;
   - a "ridge" move with log c = -log d, which moves lambda_k.

   Each proposal is symmetric on the log scale, so the acceptance ratio carries the Jacobians c and
   d^(K+1). These three moves hold Sigma_W fixed. Two further kinds of move follow in the same block.
   - Per trait k, a prior-only rebalance: row and column k of Sigma_W are scaled by d and alpha_k by
     1/d, so Sigma_P and the likelihood do not change. The step is fixed (log d ~ N(0, 0.5)) and the
     Jacobian is d^K. This keeps the expanded chain irreducible in diag(Sigma_W) when the Gibbs
     steps are switched off (a test hook).
   - Per trait pair (k, l), three off-diagonal moves: Fisher-z random walks on the Sigma_E and on the
     Sigma_W correlation, with the diagonals fixed and Jacobian (1 - r_new^2)/(1 - r_old^2); and a
     covariance ridge, Sigma_P[k,l] += delta and Sigma_E[k,l] -= delta (implemented as
     Sigma_W[k,l] += delta/(alpha_k alpha_l)), with delta = eps sqrt(v_k v_l),
     v = diag(Sigma_P + Sigma_E), and unit Jacobian. A proposal that is not positive definite is
     rejected. The pair moves were added because the total cross-trait covariance is well identified
     but its split between Sigma_P and Sigma_E is not, and Gibbs steps 4-5 traverse that split
     slowly: bulk ESS for Sigma_P[1,2] was 58 to about 140 at n = 1000 (the builder's
     measurement, not committed; the code comment gives 140 and an earlier version of this section
     gave 144).

   The three per-trait moves target p(alpha, Sigma_E | Sigma_W, y_obs). The rebalance, the Sigma_W
   correlation walk and the covariance ridge also change Sigma_W, so the whole Metropolis block
   leaves p(alpha, Sigma_W, Sigma_E | y_obs), the (a, mu, y_mis)-marginal, invariant. Step 1 then
   redraws (a, mu, y_mis) from its full conditional, so the sweep is a valid partially collapsed
   Gibbs sampler (van Dyk and Park 2008). Step sizes of the per-trait and pair moves adapt during
   burn-in only, so the kept phase is a fixed kernel.

   Evidence that the posterior is still the right one:
   - The committed test "the Metropolis-only and Gibbs-only kernels target the same posterior"
     (`tests/testthat/test-mi-posterior.R`) runs both kernels on one fixture. It requires agreement
     within 3.5 combined MCSE for lambda, the log variances and both correlations; the clean code
     reaches max |z| about 1.8. The S5 checker confirmed it catches mutants that break a prior, a
     degrees-of-freedom term or a Jacobian.
   - G2 (exactness) passes.
   - The G3 calibration run: 95% intervals for lambda and rho_P cover 92 to 98% over 150 fits
     (`evidence/README.md`).

   During the build, S1 also reported that the moves alone, the Gibbs steps alone and both together
   agree within |z| <= 1.5, and that the sampler matches MCMCglmm under the same priors. Those scripts
   were not committed, so the claims are recorded only as the builder's report. The committed test
   above replaces the first one.

   Cost: 5 to 7 times more per sweep. Measured on the Mac Studio at K = 2: 4.9 ms per sweep at
   n = 300 and 10.2 to 10.7 ms at n = 1000.
2. **Defaults.** 4 chains, each 1,000 burn-in plus 5,000 sweeps, thinned so that 1,000 draws are kept
   in total (thin 20). With 3,000 sweeps per chain the worst G3 setting already had min ESS 473.
3. **Chains run one after another within a fit.** `parallel` is not in DESCRIPTION, so campaigns
   parallelise across cells instead.
4. **`cell_interval$row`** is the integer row of the user's `traits`. A test covers rows shuffled
   relative to the tree.
5. **Start values.** Chain 1 starts at #187's REML lambda split when n <= 2,000 tips. Above that, the
   REML start would need a dense n x n matrix, so chain 1 starts at lambda = 0.9.
6. **Improper mode** fixes Sigma_P and Sigma_E at the posterior means of the same run. #187's REML
   gives only the per-trait lambdas, not the two full K x K matrices. Every improper draw is an exact
   independent draw. An internal `param_uncertainty = "both"` returns the proper and improper draws
   from one run, and the simulation uses it.
7. **Zero-length branches** are floored at 1e-6 of the tree height before Qc is built.
8. **`with_imputations()` and `pool_mi()`** accept the new provenance marker
   `pigauto_posterior_mi_v1` (orchestrator decision, option B). Conformal and mc_dropout objects are
   still refused. The docs state the congeniality scope (narrowed after the S5 review):
   - covered: analyses that are linear in the imputed traits on the imputation scale (the log
     scale for traits that `log_transform` logged), with every analysis variable among the imputed
     traits;
   - not covered: external covariates, nonlinear terms or interactions among the imputed traits, or
     a log-transformed trait analysed on its raw scale.

   Draws from `param_uncertainty = "none"` carry the marker `pigauto_posterior_plugin_diagnostic`,
   which `with_imputations()` and `pool_mi()` refuse.

## 3. Outputs

- **m imputations** (default 20) for Rubin pooling: completed latent matrices from m kept sweeps spaced
  across chains, decoded to the original scale with pigauto's existing decode path.
- **Per-cell predictive draws** (default at least 1,000 kept sweeps across chains) of each missing
  cell, from step 1(b); per-cell 95% interval = 2.5% and 97.5% quantiles. Never taken from the m
  imputations.
- **Parameter draws**: Sigma_P, Sigma_E, lambda_k, mu per kept sweep.
- **Improper mode** (`param_uncertainty = "none"`): the full sampler runs first. Sigma_P and Sigma_E
  are then fixed at their posterior means from that run (not #187's REML; see 2.5 item 6), and only
  step 1 is redrawn, as exact independent draws, so mu is still drawn. For comparison in validation
  only.

## 4. Frozen API (S3 and S4 code against this; S1 and S2 implement it)

```r
mi <- multi_impute(traits, tree, m = 20L, draws_method = "posterior",
                   posterior_control = list(
                     n_chains = 4L, n_iter = NULL, burnin = NULL, thin = NULL,
                     keep_draws = 1000L, param_uncertainty = c("full", "none"),
                     seed = NULL))
mi$datasets                       # list of m completed data.frames (original scale)
mi$posterior$cell_interval        # data.frame: row, trait, lower, upper, median (95%)
mi$posterior$diagnostics          # data.frame: parameter, rhat, ess_bulk; attr "converged"
mi$posterior$params               # list: Sigma_P, Sigma_E (arrays K x K x draws), lambda (draws x K), mu
mi$draws_method == "posterior"
```

Errors (clear messages naming the problem): any non-continuous trait (binary, categorical, ordinal,
count, proportion, zi_count, multi_proportion); multi-observation data; covariates (not modelled;
say so). `gnn` is ignored with a message. `pool_mi()` and `with_imputations()` work unchanged.

## 5. Estimands and acceptance (see `.unlazy/mi-posterior/GATES.md`)

- Downstream: slope of y on x from `nlme::gls(y ~ x, corBrownian)` and from a lambda-estimating model
  (`phylolm(model = "lambda")`, with `nlme::gls(corPagel)` as a check). Primary metric: paired bias,
  mean over replicates of (MI pooled slope - complete-data slope) in the same replicate. Also SE ratio
  (mean pooled SE / empirical SD of pooled slope) and 95% coverage, compared with complete data under
  the same analysis model.
- Per cell: coverage of the masked truth by the 95% predictive interval, in exchangeable and in
  clade-biased masks, beside conformal. Clade-biased coverage is descriptive only (review R5). It is
  expected to fall below the exchangeable G7 band when the tree model is misspecified, and results.md
  reports it as such.
- Convergence (review B3): every campaign fit records max R-hat, min ESS and `converged`. G6 and G7
  score converged fits only. A regime where more than 2% of fits are non-converged fails its gate and
  is named in the output.

### 5a. What G1 must cover (review R1, R6, S3)

1. The Qc tip-block Schur complement equals `cov2cor(vcv(tree))` to 1e-8 on an ultrametric and on a
   non-ultrametric tree.
2. `xi' Qc xi` equals the dense `t(xi) solve(dense Qc^{-1}) xi` on a 10-tip tree.
3. The (a, mu) precision is positive definite with mu included (flat prior; the root is excluded),
   including when some tips have no observed traits.
4. The PX mapping `Sigma_P = diag(alpha) Sigma_W diag(alpha)` and `a = xi diag(alpha)`, checked on a
   hand-computed K = 2 case.
5. Reduction to #187: with Sigma_P = lambda Sigma and Sigma_E = (1 - lambda) Sigma, the implied dense
   marginal covariance equals Sigma %x% (lambda R + (1 - lambda) I).
6. Decode indexing: only tip rows of a enter the residuals and the completed data; the internal rows
   never do.
7. Observed cells are never altered in any imputation.
8. Seed reproducibility.
9. API errors: non-continuous traits, multi-observation data, covariates; the `gnn` message.
10. Output shapes match section 4, and `pool_mi()` / `with_imputations()` run on the result.
11. The improper mode (`param_uncertainty = "none"`) holds the parameters fixed across sweeps.

### 5b. G2 and G3 as specified (review R2, R3)

- G2a: Sigma_P and Sigma_E positive definite and fixed, Sigma_E not small; n = 40, K = 2; 20,000
  draws of y_mis. Compare against the dense conditional of y_mis | y_obs under
  N(1 mu', Sigma_P %x% R + Sigma_E %x% I_n), with mu also held fixed (an internal test hook that fixes mu along with the covariances).
  - mean within 4 MCSE per cell;
  - covariance entries within 0.03.
- G2b: Sigma_E = 1e-6 I. The conditional mean matches `draw_conditional_bm()` within 1e-3 (on the
  scale of the data).
- G3: n = 1000, 3 seeds per truth setting. Truth lambda pairs (1, 1), (0.5, 0.5), (0.3, 0.9) and
  (0.05, 0.95); the phylogenetic correlation truth is 0.7 or 0.
  - posterior-mean lambda within 0.1 of truth for every trait;
  - every off-diagonal element of cov2cor(Sigma_P) within 0.1 where lambda_k >= 0.3 for both traits;
  - posterior lambda within 2 posterior SD of #187's REML lambda.

  Boundary bias (posterior mean minus truth at lambda 0.05 and at 1) is written to the gate output
  and reported, even when it passes.

### 5c. Harness decisions (2026-09-24)

Orchestrator decisions on the simulation harness (`script/mi_gls/`), made after the verified review
of the v2 harness (findings dgp#0, estimands#0-5, gates#0-4, cluster#0-6). They change how G6 and
G7 are computed, not the thresholds of the approved plan.

- **D1. Downstream-coverage truth.** Regimes 1-16 keep 0.7. It is exact for any analysis model there,
  because the DGP covariance is proportional, so E[y | x] = 0.7 x. In regimes 17-24,
  `Sigma_P[1,2] + Sigma_E[1,2]` is the OLS estimand only; neither corBrownian GLS nor
  phylolm(lambda) targets it, and complete-data coverage of it is near 0. There the truth is a
  pseudo-truth: the mean complete-data slope over the expected reps of that (regime, analysis
  model), computed in `03_summarise_v2.R`. `covered` is recomputed in 03 for complete and every
  method from the saved estimate, se and df (t quantile with the saved df; normal if df is
  missing). The `true_beta_pop` coverage is kept only as a descriptive column (`coverage_pop`).
- **D2. Complete-data reference rows.** 03 writes `method == "complete"` rows per (regime, analysis
  model) with SE ratio, coverage and R. G6 rule 3 (SE ratio in [0.90, 1.15] under phylolm) is
  selectable by env `MI_SE_RULE`: `relative` (default since CP1, section 5d.3: MI SE ratio divided
  by the complete-data SE ratio) or `absolute` (the original approved plan). Both numbers are
  printed per regime, plus an `ANALYSIS_MODEL_SE_RATIO` line whenever the complete-data ratio is
  itself outside the band. Shinichi chose `relative` at CP1 (5d.3).
- **D3. Completeness, fail-closed.** Expected reps come from env `MI_N_REPS` (default 200). The
  expected grid is built from `regimes.R`, never from the files. A missing rep file counts as a fit
  failure and as non-converged. A missing (regime, analysis model, method) row, R or n <= 0, or any
  non-finite gated quantity fails; nothing is skipped. G7 does the same over the expected MCAR set
  (trait x in x_only regimes; x and y in both-missing regimes), and also fails a row where a
  converged rep has no cell intervals for that trait (`n_reps_scored < n_converged`). The 2% limits
  (fit failures, non-convergence) fail only above 2%, compared on counts: exactly 4/200 passes and
  5/200 fails. A run restricted with `MI_REGIMES` or fewer than 200 reps can pass its rules but
  never prints the G6 or G7 token.
- **D4. G6 rule 4 (proper vs improper SE ratio).** Superseded by 5d.4: reported, not gated. The
  mean posterior_full vs posterior_none SE ratio is printed per analysis model within each
  both-missing block (9-16 stress test, 17-24 Kronecker, 33-40 twins), and per-regime pairs are
  printed too. Neither decides the gate, and a non-finite pair is counted, not failed. A missing
  posterior_none row still fails, through the completeness check (D3).
- **D5. Mean shortfall** is the positive part, `mean(pmax(complete_coverage - coverage, 0))` over the
  posterior_full rows, so over-coverage in one row cannot cancel under-coverage in another. The
  per-row rule, coverage >= complete - 0.05, is unchanged.
- **D6. posterior_none runs only where both traits are missing** (it exists for rule 4), and the
  expected grid reflects that. It comes from the same chain run as posterior_full, via
  `param_uncertainty = "both"` (section 2.5), so it adds no second 4-chain run.
- **D7. Conformal arm.** Each cell also runs `impute(<same masked df>, tree, gnn = FALSE, seed =
  <cell seed>)` with pigauto defaults otherwise, and scores `prediction$conformal_lower/upper` on
  exactly the masked cells, against the same truth on the original scale as the posterior per-cell
  coverage. Coverage and mean width are recorded per trait, next to the posterior widths.
  `05_cell_coverage.R` prints conformal beside posterior_full per regime x trait x mechanism,
  descriptive and never gated. A conformal failure is recorded and does not stop the cell.
- **D8. Frozen code and routing.** Each cell output records `code_sha` (env `MI_POST_SHA`), the
  pigauto version, R version and host. The cell script keeps `devtools::load_all()` from its working
  directory, which in the campaign is a `git archive` of one SHA. Totoro is the primary route
  (`12_totoro_campaign.sh`: at most 150 parallel cells, single-threaded BLAS, n = 1000 regimes
  first, completed cells skipped, `CELL_FAILED` lines without stopping the rest, one code SHA per
  output directory). Both runners refuse a git checkout as the code directory, since `load_all()`
  would also source untracked files. 03 writes the G6/G7 CSVs into `docs/dev-log/mi-posterior/`,
  which are copied back to the worktree where the ledger CHECK lines read them.
  `11_fir_array.sbatch` is the DRAC fallback, fixed per the cluster findings.

### 5d. CP1 decisions (Shinichi, 2026-09-24, after the smoke run)

1. The full campaign was launched on Totoro from commit 69670d4:
   - simulation: 24 regimes x 200 reps on 130 cores;
   - real data: 10 cells on 10 cores;
   - G3 calibration rerun: 200 fits on 10 cores.
2. **G3 becomes a calibration check** (`script/mi_gls/gate_calibration.R`). The per-fit 0.1 tolerance
   is about 1.7 posterior SD at lambda = 0.5, so honest fits fail it. 150 fits at commit 7a0f470
   showed the posterior is calibrated (`evidence/README.md`). Rules:
   - 50 or more seeds per setting, all from one code SHA;
   - converged in at least 98% of fits per setting;
   - 95% interval coverage in [0.88, 1.00] for lambda_k where the truth is inside (0, 1) (settings
     B-D; the lambda = 1 truths of A sit on the boundary), and for rho_P where both lambda_k >= 0.3
     (settings A-C);
   - |mean bias| <= 0.03 for every lambda_k, and for rho_P where its coverage is gated;
   - REML agreement in at least 90% of fits.
3. **The SE ratio is judged relative to complete data** under the same analysis model:
   (MI SE ratio) / (complete SE ratio) in [0.90, 1.15]. The complete-data ratio itself can sit near
   0.80 where the analysis model is misspecified (regimes 21 and 23 under phylolm, per the harness
   review).
4. **Proper vs improper SE ratio is reported, not gated.** On the 30-tip test fixture (S1
   measurement; `tests/testthat/test-multi-impute-posterior.R`), fixing Sigma at its posterior mean
   gave a per-cell predictive variance 1.4-1.6% larger than the proper one (intervals about 0.7-0.8%
   wider), a Jensen effect. No measurement at the campaign sizes is recorded. So "proper > improper"
   need not hold when everything is right.
5. **The real-data 5% slope criterion is reported, not gated** (decision R3, unchanged).

### 5e. CP2 follow-up decisions (Shinichi, 2026-09-24, after the diagnosis)

Evidence: `diagnosis.md`.

1. **In-model twin regimes 25 to 40.** Regimes 1 to 16 simulate from the raw covariance of
   non-ultrametric trees, which is outside the model the sampler fits. So each of them gets an
   in-model twin, regimes 25 to 40:
   - same seeds, trees and masks as the source regime;
   - each tip's row divided by sqrt(diag(V_sim)), so vec(Y) ~ N(0, Sig %x% (lambda R + (1 - lambda) I))
     with R = cov2cor(vcv(tree)).

   The G6 and G7 gates apply to the in-model regimes (17 to 40). Regimes 1 to 16 are kept, unchanged,
   as a reported misspecification stress test.
2. **Automatic chain extension.** When the convergence rule fails after the default run, every chain
   continues from its saved state and RNG state for another n_iter kept-phase sweeps. This repeats
   up to `max_extend` times (default 3, so at most 4x the default length) before the draws are
   returned.
   - Adaptation stays confined to the burn-in, so an extended chain equals a longer run.
   - A fit that converges first time is byte-identical to the previous code.
   - Motivation: all 44 non-converged campaign fits in regimes 1 to 24 failed on ESS only (min bulk
     ESS 190 to 395, max split R-hat at most 1.028; `evidence/diagnosis/campaign_cells.csv`). The 32
     in regimes 5, 21 and 23, the regimes over the 2% limit, were re-run at 4x length (burn-in 4,000
     and 20,000 sweeps). That re-run also lengthened burn-in, which the extension does not do. All 32
     converged, and the pooled slopes did not move (`diagnosis.md`, Failure 2). The other 12 were not
     re-run at 2x or 4x.

   Users pay the extra time only when a fit needs it (usability, D-139).
3. **The re-run uses the new code for:**
   - regimes 25 to 40 (3,200 fits);
   - the 44 previously non-converged cells of regimes 1 to 24.

   The 4,756 converged cells of regimes 1 to 24 keep their 69670d4 results. The new code must
   reproduce them byte for byte; this is checked on a sample before the campaign.
4. **Not changed:** the priors, including the Sigma_E prior (its smaller-scale variant failed
   numerically in 2 of 24 fits), and the `cov2cor(vcv(tree))` convention. Both are listed as open
   items.

## 6. Out of scope here

Discrete, mixed and multi-observation data; GNN blending; covariates in the imputation model; changing
the default draws method; editing `R/joint_mvn_solver.R` (its plug-in covariance shrinkage goes to the
lambda lane as a note: `script/mi_gls/` six-tree check, 0.67 to 0.46).
