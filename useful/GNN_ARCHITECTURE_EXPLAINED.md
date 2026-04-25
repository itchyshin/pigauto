# pigauto GNN architecture — explained, with diagnosed bottlenecks

> Generated 2026-04-25 to make the architecture transparent and identify
> why the GNN keeps losing to `phylolm-lambda` BLUP on simulated linear
> covariate-effect data even after Fixes A through D were applied.

## End-to-end data flow

```
USER INPUT:
  traits      data.frame  n_species (or n_obs in multi-obs) × p_traits
  tree        phylo       n_species tips
  covariates  data.frame  n_species (or n_obs) × n_user_cov   (OPTIONAL)

STEP 1: preprocess_traits()
  X_scaled      n_obs × p_latent    (z-scored / log-transformed / one-hot
                                     encoded depending on trait type)
  trait_map     list of descriptors per trait (type, levels, mean, sd)
  obs_to_species n_obs              (only present in multi-obs mode)
  covariates    n_obs × n_user_cov  (z-scored)

STEP 2: build_phylo_graph()
  coords  n_species × k_eigen        Laplacian eigenvectors (spectral
                                     phylogeny features)
  adj     n_species × n_species      Gaussian kernel on cophenetic dist.
  D_sq    n_species × n_species      squared cophenetic distance (for
                                     rate-aware attention)

STEP 3: fit_baseline()                ← FIXED at this point, NOT trained
  mu  n_species × p_latent           BM kriging predictions per trait,
                                     IGNORING covariates entirely.
  se  n_species × p_latent           Conditional-MVN standard errors.

  For continuous: each column independently fits BM(tree) via
    cov2cor(vcv(tree)) and conditional E[y_missing | y_observed].
  This is essentially `phylolm(y ~ 1, model="lambda")` — INTERCEPT-ONLY,
  no cov terms.

STEP 4: ResidualPhyloDAE training
  See "Model architecture" below.
  The GNN is the only part of pigauto that sees user covariates.

STEP 5: r_cal calibration on val set
  Per-trait grid search of r_cal that minimises pigauto val RMSE.
  Optionally bounded below by safety floor: r_cal must satisfy
    pigauto_val_RMSE ≤ mean_val_RMSE.

STEP 6: predict.pigauto_fit()
  pred = (1 − r_cal) * mu + r_cal * delta + cov_linear_fixed_effects
```

## Model architecture (`ResidualPhyloDAE` in `R/model_residual_dae.R`)

```
Input:
  x       (n_obs × p_latent)              corrupted trait input
  coords  (n_species × k_eigen)           spectral phylo features
  covs    (n_obs × cov_dim)               cov_dim = p_latent + 1 + n_user_cov
                                          [baseline_mu | NA-mask | user_covs]

If multi-obs: coords_obs = coords[obs_to_species]
Else:        coords_obs = coords

combined = concat(x, coords_obs, covs)    # 2D tensor, big

ENCODER (2 linear layers, single MLP):
  h = enc1(combined) → ReLU → dropout
  h = enc2(h)        → ReLU
  # h is now n_obs × hidden_dim (e.g. 32)

# Aggregate to species level if multi-obs:
if multi_obs: h_species = scatter_mean(h, obs_to_species, n_species)
else:         h_species = h                 # n_species × hidden_dim

GNN MESSAGE PASSING (n_gnn_layers transformer blocks):
  for layer in 1..n_gnn_layers:
    h_species = GraphTransformerBlock(h_species, adj, D_sq)
    # Each block:
    #   - multi-head attention with phylo prior (log(adj)+gauss bandwidth)
    #   - FFN with hidden_dim → ffn_mult*hidden_dim → hidden_dim
    #   - LayerNorm pre-each-block, residual skip
    #   - FFN output init to ZERO so block ≈ identity at step 0

# Broadcast back to obs level if multi-obs:
if multi_obs: h = h_species[obs_to_species]
else:         h = h_species

COVARIATE REFINEMENT (post-Fix B, fires when n_user_cov > 0):
  user_covs = covs[:, last n_user_cov columns]
  cov_h     = cov_encoder(user_covs)          # 4-layer MLP, hidden_dim out
  h         = h + obs_refine(cat(h, cov_h))   # residual injection

DECODER (2 linear layers):
  h     = dec1(h) → ReLU
  delta = dec2(h)                             # n_obs × p_latent

LINEAR COV HEAD (Fix C+D, post-2026-04-25, fires when n_user_cov > 0):
  fixed_effects = cov_linear(user_covs)       # direct linear regression
  # Initialised at 0.01x scale to start small.

GATE:
  rs = sigmoid(res_raw) * gate_cap            # per-trait, in (0, gate_cap)
  res_raw initialises at -1 (continuous) or near 0 (discrete)

OUTPUT:
  pred = (1 − rs) * baseline_mu + rs * delta + fixed_effects
```

## Loss & training

```
Per minibatch (full batch in pigauto, no minibatching):
  pred = forward(x_corrupted, ..., baseline_mu = baseline_mu)
  loss_main   = compute_mixed_loss(pred, y_true)        # MSE/BCE/CE per trait
  loss_shrink = MSE(delta − baseline_mu)                # push delta → baseline
  loss_gate   = MSE(rs)                                 # push rs → 0
  total       = loss_main + 0.03*loss_shrink + 0.01*loss_gate
  total.backward(); optimiser.step()
```

**Three regularisations push delta toward baseline:**

1. `lambda_shrink = 0.03` weight on MSE(delta − baseline_mu) — actively
   prevents delta from deviating from the BM baseline.
2. `lambda_gate = 0.01` weight on MSE(rs) — actively pushes the gate
   toward zero (i.e. baseline-only blend).
3. `gate_cap = 0.8` upper-bounds rs to 0.8 even when training fully
   opens the gate.

## Diagnosed bottlenecks

### 1. The BM baseline IGNORES user covariates

`fit_baseline()` fits `phylolm(y ~ 1, model="lambda")` per trait — pure
intercept + BM, no cov terms. Compare to `phylolm-lambda BLUP`:

```
phylolm BLUP:
  pred = X * beta_hat        +    V_lambda BLUP correction
         (linear cov part)        (phylo residual part)

pigauto baseline:
  mu = 0                     +    V_lambda BLUP correction
       (no cov part)              (phylo residual part)
```

The entire linear-cov contribution `X * beta_hat` has to be re-derived
inside the GNN's `delta`. The GNN starts from scratch, learning by
gradient descent through 4-6 nonlinear layers what `phylolm` recovers
analytically with one solve.

**This alone explains most of the linear-effect underperformance.**

### 2. Three regularisations push delta toward the BM baseline

Even if the GNN figured out the right cov coefficients, three penalties
pull it back:

- `lambda_shrink * MSE(delta − baseline)` actively penalises any cov
  contribution to delta.
- `lambda_gate * MSE(rs)` actively closes the gate so delta isn't used.
- `gate_cap = 0.8` caps the gate even when training pushes it up.

These exist to prevent the GNN from over-fitting tiny datasets — they
ARE useful safety properties. But they're set at fixed values that may
be too aggressive for cov-rich, well-identified problems.

### 3. The validation calibration plus safety floor

After training, `r_cal` is grid-searched per trait on the val set to
minimise val RMSE, then optionally floored by safety_floor=TRUE so
`pigauto_val_RMSE ≤ mean_val_RMSE`. If the GNN is undertrained or noisy
on val, the calibration gives `r_cal ≈ 0` and predictions collapse to
the baseline (which has no cov info).

### 4. The corrupted-input task is denoising, not regression

The GNN trains on a denoising autoencoder objective: input has random
masks + corrupted cells, output predicts uncorrupted y. This is good
for trait imputation but doesn't directly optimise the cov→trait
regression. phylolm directly fits regression coefficients with
maximum-likelihood; the GNN finds them via gradient descent through
many indirect layers.

### 5. Single hidden_dim bottleneck

Default `hidden_dim` is small (typically 32). The encoder concatenates
~50–80 input dimensions and projects to 32, losing information. The
cov_encoder added in Fix B helps, but the main GNN stack still operates
at hidden_dim=32.

### 6. Information dilution across the GNN stack

Even with Fix A/B, the user covariates only enter at:
- The encoder input (mixed with everything else, projected to hidden_dim)
- The obs_refine residual (after the GNN, mixes covs back in)

The GNN message-passing layers themselves never see covariates. So the
GNN propagates phylogenetic similarity but cannot use covariates for
the propagation. A "phylogenetically-similar species lives in a
covariate-different climate" signal is impossible to encode in the
current message passing.

## What fixes E and beyond might look like

Now that the architecture is mapped, here are the higher-impact
interventions in increasing order of intrusiveness:

### Fix E: lower lambda_shrink and lambda_gate when covariates are present

Make those constants adaptive — when n_user_cov > 0, halve the shrink
penalty. The user paid for a covariate signal; we shouldn't penalise
deviation from the cov-less baseline.

### Fix F: train cov_linear on a regression-only loss head

Add an auxiliary loss `MSE(cov_linear(covs), y − baseline_mu)` so the
linear path is fitted directly to the residual-from-baseline regression
target. This converges to phylolm's `beta_hat` quickly without the
indirect blend gradient.

### Fix G: feed user_covs to baseline_mu — fit phylolm internally

Replace `fit_baseline` with `fit_baseline_with_covs` that runs `phylolm(
y ~ user_covs, model="lambda")` per trait. The baseline becomes
`X*beta_hat + BM_BLUP_residual` — the same as phylolm's full prediction.
The GNN's role becomes adding NONLINEAR corrections only. This is the
proper decomposition: linear cov + BM phylo + GNN nonlinear.

### Fix H: cov-aware GNN layers

Inject covariates into every transformer block via cross-attention or
concat-and-project. Each layer attends to (h_species, covariates)
jointly. This lets the GNN propagate covariate-aware signals through
the phylogenetic graph.

## Recommendation

I recommend trying **Fix G first** because it directly addresses the
root cause: `pigauto`'s baseline is pure-phylo, while `phylolm-lambda`'s
baseline is linear-cov + phylo. By upgrading the baseline to include
covariates linearly, the GNN's job becomes learning *only* the
nonlinear / interactive corrections — exactly the regime where it
should win.

Concretely Fix G means: replace `fit_baseline_continuous` with a call
to `phylolm::phylolm(y ~ ., data = X_with_covs, phy = tree,
model = "lambda")` and use BLUP for prediction. Then the GNN delta is
honestly a residual model — and the lambda_shrink penalty toward
baseline becomes a DESIRABLE inductive bias rather than a handicap.

## Status of fixes

- ✅ Fix A: `n_user_cov = n_cov_cols` always (single-obs + multi-obs both create obs_refine)
- ✅ Fix B: dedicated cov_encoder MLP for richer cov capacity
- ✅ Fix C: cov_linear direct linear cov head added to delta
- ✅ Fix D: cov_linear contribution moved OUTSIDE the (1−r)/r blend gate
- 🔜 Fix E (proposed): adaptive lambda_shrink/lambda_gate when covs present
- 🔜 Fix F (proposed): auxiliary regression loss on cov_linear
- 🔜 Fix G (proposed, recommended next): use `phylolm-lambda` as baseline
  when covs supplied — pigauto's GNN becomes pure nonlinear residual
- 🔜 Fix H (proposed): cov-aware GNN layers via cross-attention
