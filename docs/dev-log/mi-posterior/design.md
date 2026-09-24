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
mu, and y_mis. Also stored: cov2cor(Sigma_P) and cov2cor(Sigma_E). alpha and Sigma_W are not
identified and are neither reported nor diagnosed.

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

## 3. Outputs

- **m imputations** (default 20) for Rubin pooling: completed latent matrices from m kept sweeps spaced
  across chains, decoded to the original scale with pigauto's existing decode path.
- **Per-cell predictive draws** (default at least 1,000 kept sweeps across chains) of each missing
  cell, from step 2; per-cell 95% interval = 2.5% and 97.5% quantiles. Never taken from the m
  imputations.
- **Parameter draws**: Sigma_P, Sigma_E, lambda_k, mu per kept sweep.
- **Improper mode** (`param_uncertainty = "none"`): parameters fixed at their posterior mean (or #187's
  REML estimate); only steps 1-2 run. For comparison in validation only.

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

## 6. Out of scope here

Discrete, mixed and multi-observation data; GNN blending; covariates in the imputation model; changing
the default draws method; editing `R/joint_mvn_solver.R` (its plug-in covariance shrinkage goes to the
lambda lane as a note: `script/mi_gls/` six-tree check, 0.67 to 0.46).
