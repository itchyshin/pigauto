# GNN architecture findings — covariates are starved in single-obs mode

> **Status:** confirmed bug / design limitation. Found 2026-04-25 while
> diagnosing why v1 GNN earnings sim showed pigauto LOSING to phylolm-BLUP
> on linear and nonlinear effects. The user explicitly flagged this as
> the right thing to investigate.
>
> **Action:** address AFTER current Day 1 v2 sims and Day 2 multitrait sim
> finish. The user said: "of course we can work on GNN and non-linear stuff
> once we finish what we are doing now so I want you to remember this is
> an important point". Hence this file.

## The smoking gun

In `R/fit_pigauto.R` line 398:

```r
n_user_cov <- if (multi_obs) n_cov_cols else 0L
```

In **single-observation mode** (the typical comparative-biology setup
where each species has one row of trait data), `n_user_cov` is forced to
**zero regardless of how many covariate columns the user passes**.

This silently disables the only architectural path that re-injects user
covariates after the GNN message-passing stack:

```r
# In ResidualPhyloDAE$forward(), R/model_residual_dae.R line 268:
if (self$n_user_cov > 0L && !is.null(self$obs_refine)) {
  ...                                  # FALSE in single-obs mode
}
```

## What pigauto's GNN actually sees with single-obs covariates

With `n_user_cov = 0`, user covariates enter the model **exactly once**, in
the input encoder (line 200 of `model_residual_dae.R`):

```r
combined <- torch::torch_cat(list(x, coords_obs, covs), dim = 2L)
h <- self$enc1(combined)              # one linear projection
h <- self$act(h); h <- self$drop(h)
h <- self$enc2(h); h <- self$act(h)   # one more linear projection
```

After that, the signal goes through:

1. ReLU + dropout (kills small linear contributions)
2. ~~A few~~ several GNN / transformer layers operating on the
   **species-level phylogenetic adjacency only** — no covariate info
   propagates through these
3. A 2-layer decoder

By the time the decoder produces `delta`, the covariate signal has been
washed through several layers that didn't see it again. The GNN message
passing is doing pure phylogenetic graph propagation — the same thing
phylolm-BLUP does analytically.

**This is exactly what the v1 sim shows:** pigauto behaves like
"phylolm-BLUP-but-noisier" because architecturally it almost is one in
single-obs mode. The covariates barely participate in the nonlinear
GNN computation.

## CLAUDE.md confirms this is documented (and stale)

`CLAUDE.md` actually says, in the Architecture section:

> **Input**: ... `covs` (`n_obs × cov_dim` — currently baseline mean +
> NA-mask indicator; **no user covariates yet**) ...

That `no user covariates yet` was true historically. The `obs_refine`
multi-obs path was added in v0.6.0+. But the **single-obs path was never
upgraded** to feed user covariates back in.

## What the fix probably looks like

There are progressively-richer fixes:

### Fix A (minimal): always set `n_user_cov = n_cov_cols`

```r
# R/fit_pigauto.R line 398
n_user_cov <- n_cov_cols   # was: if (multi_obs) n_cov_cols else 0L
```

And let `obs_refine` fire in single-obs mode too. In single-obs the
`index_select(1L, obs_to_species)` step is a no-op (obs == species), so
the residual injection at line 273 just does:

```r
h <- h + obs_refine(cat(h, user_covs))
```

i.e. one MLP that takes (post-GNN species hidden state, user covs) and
adds a covariate-conditioned correction. This gives the model **a
nonlinear path from covariates to delta that doesn't go through the GNN
phylogenetic propagation**.

This is the smallest change with the biggest architectural payoff. It
should be the first thing we try.

### Fix B: covariate-aware GNN layers

Instead of one residual injection, make every GNN layer covariate-aware:

```r
# Each GNN layer gets cross-attention with covs as keys/values
m <- attention(h_species, K = cov_proj(covs), V = cov_proj(covs))
```

Richer but bigger change to the model.

### Fix C: dedicated cov MLP that bypasses GNN entirely

```r
delta = decoder(GNN(h_phylo) + cov_MLP(covs))
```

Adds a parallel non-phylogenetic path the gate can blend in.

### Recommendation

Try **Fix A first** because:
1. It's a one-line config change in `fit_pigauto.R` plus a guard removal in `model_residual_dae.R`
2. The `obs_refine` MLP already exists and is well-tested in multi-obs
3. We can immediately re-run the GNN earnings sim to see whether the
   nonlinear / interactive losses go away

If Fix A doesn't move the needle, Fix B or C in turn.

## Why this matters for the paper

The v1 sim says **"pigauto's GNN does not extract nonlinear cov signal
beyond phylolm"** — but that's because the GNN is *architecturally
prevented from doing so* in single-obs mode, not because GNNs can't do
that in principle. The earlier `bench_multi_obs.R` results (10–19 % lift
on simulated CTmax with `acclim_temp` covariate) showed exactly this:
when `obs_refine` IS active (multi-obs mode), the GNN extracts covariate
signal. Single-obs is the broken case.

**Re-running the GNN earnings sim after Fix A will tell us whether the
real architectural claim** ("the GNN extracts nonlinear cov signal that
phylolm cannot") holds when the GNN is actually given access to the
covariates at every layer.

If post-fix pigauto beats phylolm-lambda BLUP on nonlinear / interactive
cells, the paper has a real story:
> "We identified and fixed a path-blocking issue in pigauto's single-obs
> mode that prevented user covariates from re-entering the GNN after
> phylogenetic message passing. With the fix, pigauto extracts
> nonlinear covariate-effect signal that phylolm-lambda BLUP cannot
> capture, by [X] % RMSE on [regimes]."

If post-fix pigauto STILL doesn't beat phylolm, then we know the
limitation is not just plumbing — it's the architecture itself.

## Plan

1. ⏸ Wait for current Day 1 v2 sims to finish (~3:30 PM today).
2. ⏸ Run Day 2 multitrait sim (~30 min).
3. ⏸ Harvest both, write decisive writeup.
4. ▶ **THEN** apply Fix A. Single-line change. Re-run v2 tree300 sim
   only (fastest, n=300, ~4 hours).
5. Compare post-fix v2 results to current v2 results. If post-fix
   pigauto wins on nonlinear / interactive: write up the architectural
   bug as a paper finding. If not: try Fix B or C. If still nothing:
   the architecture itself genuinely doesn't capture nonlinear cov
   signal beyond phylolm and the paper claim must be scoped to
   safety-floor / multi-obs / mixed-types / UQ.

This file is the persistent record of the finding. Don't lose it.
