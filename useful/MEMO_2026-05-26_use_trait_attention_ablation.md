# `use_trait_attention` ablation — 60-rep, N=2000

Run: 2026-05-26.  Lead: @b1805 (cluster run).  Code:
[`script/sim_bench/ablation_use_trait_attention.R`](../script/sim_bench/ablation_use_trait_attention.R)
(landed via PR #116, after rebase onto the post–PR #117 main).
Raw data archived as the `results_full_60.xlsx` attachment on
[issue #106](https://github.com/itchyshin/pigauto/issues/106).

## Inferential target

Is `use_trait_attention = TRUE` (added in v0.9.3, default `FALSE`)
genuinely redundant against the joint-MVN baseline at BIEN scale, or is
the encoder taking a "lazy" shortcut by routing `baseline_mu` through the
residual path and ignoring the attention tokens?

## Data-generating process

This iteration **fixes** the single-trait DGP flaw that Shinichi flagged
on the earlier ablation attempt. The corrected DGP injects two
correlated continuous traits `y1`, `y2` (latent cross-trait Σ), so the
joint-MVN baseline (`fit_baseline()` → joint MVN path) genuinely fires
and `baseline_mu` carries real multi-trait structure. Without that,
`use_trait_attention` has nothing to attend across.

Other DGP knobs (per `ablation_use_trait_attention.R`):

* `N_SPECIES = 2000`
* 4 process scenarios: `bm_strong`, `ou_strong`, `nonlinear_cov`,
  `weak_signal`
* 3 missingness mechanisms: `phylo_MAR`, `trait_MAR`, `trait_MNAR`
* `MISS_RATE = 0.30`, `DEP_STRENGTH = 1.5`
* 5 seeds per cell × 12 cells = 60 replicates total

## Three arms

| arm             | `use_trait_attention` | extras                                        |
|-----------------|-----------------------|-----------------------------------------------|
| `pigauto_OFF`   | `FALSE`               | v0.9.3 default                                |
| `pigauto_ON`    | `TRUE`                | full B3 path                                  |
| `pigauto_ON_L0` | `TRUE`                | `baseline_mu` masked + `lambda_shrink = 0`    |

The `pigauto_ON_L0` arm **disarms the lazy-optimizer trap**: masking
`baseline_mu` and dropping the L2 shrink penalty forces the encoder to
rely entirely on the attention tokens. If `pigauto_ON` were quietly
copying `mu` and ignoring attention, `pigauto_ON_L0` should reveal new
useful signal under the attention path.

## Result

z-RMSE, conformal coverage and interval width, averaged across all 12
scenario × mechanism cells (60 replicates each method):

| method            | z_rmse | coverage | width |
|-------------------|--------|----------|-------|
| `column_mean`     | 1.291  | —        | —     |
| `bm_kriging`      | 1.163  | —        | —     |
| **`pigauto_OFF`** | **1.038** | **0.887** | **2.65** |
| `pigauto_ON`      | 1.056  | 0.884    | 2.67  |
| `pigauto_ON_L0`   | 1.057  | 0.887    | 2.68  |

## Reading

1. **pigauto comfortably beats `column_mean` and `bm_kriging`** —
   z-RMSE drops from 1.16 (BM kriging) to 1.04 (pigauto_OFF). The GNN
   correction + calibrated gate is doing useful work at this scale.
2. **`use_trait_attention` regresses pigauto** — `pigauto_ON` lifts
   z-RMSE from 1.038 to 1.056 (+1.7%) and slightly widens conformal
   intervals (2.65 → 2.67). Coverage is statistically indistinguishable.
3. **Even disarmed, the attention path doesn't help** — `pigauto_ON_L0`
   is no better than `pigauto_ON`. The "lazy optimizer" hypothesis is
   refuted: the network is not hiding useful structure behind `mu`. At
   N=2000 with a working joint-MVN baseline, there is simply no
   additional cross-trait signal for the attention mechanism to recover.
4. **Conclusion.** The joint-MVN / threshold-joint baseline already
   encodes the dominant cross-trait Σ at BIEN scale. The within-row
   attention is mathematically redundant against it and adds noise.

## Implication for the codebase

`use_trait_attention` stays default-`FALSE`. It is kept as opt-in
because b1805's local smoke test at N=300 (smaller than the baseline's
convergence regime) did show attention dropping validation loss before
the joint-MVN catches up — there may be a small-N regime where it
helps. At BIEN scale it does not, and that's documented here and in the
GNN-architecture vignette.

The autoresearch sweep finding ("attention regresses on BIEN") is now
independently corroborated by a 60-rep ablation on a corrected
multi-trait DGP. No further architecture knob is on the table for BIEN
at this scale; the next direction is environmental covariates
(WorldClim).
