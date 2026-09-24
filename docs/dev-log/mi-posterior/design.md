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

1. **(a, mu) | y, Sigma_P, Sigma_E.** Joint Gaussian. With Z mapping nodes to tips, the precision of
   vec(a) is Sigma_P^{-1} %x% Q + Sigma_E^{-1} %x% (Z'Z) (sparse), and the linear term is
   (Sigma_E^{-1} %x% Z') vec(y - 1 mu'). Draw mu | a, y then a | mu, y, or block them; draw a with one
   sparse Cholesky of the (K N) system (Matheron / sample-from-precision), reusing the symbolic
   factorisation across sweeps.
2. **y_mis | a, mu, y_obs, Sigma_E.** Per tip: e_i = y_i - mu - a_i is K-variate normal N(0, Sigma_E);
   condition its missing components on its observed components. Tips with no observed trait draw
   e_i from N(0, Sigma_E).
3. **Sigma_P | a.** Inverse-Wishart: IW(nu_P + N, S_P + a' Q a) (a as an N x K matrix), with
   parameter expansion (below).
4. **Sigma_E | e.** IW(nu_E + n, S_E + e'e).

Priors: nu = K + 1; S = 0.01 x diag of the observed-cell variances on the latent scale (weak).
Parameter expansion for Sigma_P (Gelman 2006; MCMCglmm's `alpha.mu`/`alpha.V`): a = diag(alpha) a*,
with working parameters alpha drawn each sweep; this gives scaled-F priors on the phylogenetic
standard deviations and improves mixing when a variance is near zero (lambda near 0) or when
Sigma_E is near zero (lambda near 1). Record the implementation choice in the code header.

Chains: 4 by default, dispersed starts (one at #187's REML lambda split, others perturbed); burn-in and
thinning set from the smoke. Diagnostics: split R-hat and bulk ESS for each element of Sigma_P,
Sigma_E and lambda_k, computed internally (no new dependency).

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
  clade-biased masks, beside conformal.

## 6. Out of scope here

Discrete, mixed and multi-observation data; GNN blending; covariates in the imputation model; changing
the default draws method; editing `R/joint_mvn_solver.R` (its plug-in covariance shrinkage goes to the
lambda lane as a note: `script/mi_gls/` six-tree check, 0.67 to 0.46).
