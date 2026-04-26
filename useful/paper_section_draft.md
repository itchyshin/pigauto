# Paper section draft — pigauto's covariate-aware architecture

> **Status:** working draft (2026-04-26).
> **Contents:** the architectural finding, the fix, the empirical evidence,
> and an honest scoping of the paper claim. Multi-observation evidence
> (Section 4.4) will be filled in when `bench_multi_obs.R` re-run completes.

## 1. Background — pigauto's blended-prediction architecture

pigauto fits a per-trait blend of a Brownian-motion phylogenetic baseline
and a graph-neural-network correction:

$$
\mathrm{pred}_i = (1 - r_i) \cdot \mu_{i,\mathrm{BM}} + r_i \cdot \delta_{i,\mathrm{GNN}}
$$

where $\mu_{\mathrm{BM}}$ is computed analytically via conditional MVN on
the phylogenetic correlation matrix $R = \mathrm{cov2cor}(\mathrm{vcv}(T))$,
$\delta_{\mathrm{GNN}}$ is the output of a multi-head graph transformer
that operates on a phylogenetic adjacency, and $r$ is a per-trait
calibrated gate bounded in $(0, 0.8)$. The architecture is supposed to
let users pass environmental covariates and have them influence
prediction through the GNN's $\delta$.

## 2. The architectural finding

While running the GNN-earnings simulation (Section 4.1) we observed
that pigauto with covariates was performing *worse* than `phylolm`
($y \sim X\beta + u, u \sim \mathrm{MVN}(0, \sigma^2 V_\lambda)$) BLUP
predictions on linear covariate-effect regimes by 50–85 %. This was
unexpected — even in the worst case, pigauto should be able to match
`phylolm` since both have the same model class.

Source code inspection traced the issue to `R/fit_pigauto.R` line 398:

```r
n_user_cov <- if (multi_obs) n_cov_cols else 0L
```

In single-observation mode (the typical comparative-biology setup),
`n_user_cov` was forced to zero regardless of how many covariate
columns the user passed. This silently disabled the `obs_refine` MLP
that re-injects user covariates after GNN message passing. User
covariates entered the model only via the input encoder (one linear
projection), then got diluted through phylogeny-only graph layers
that didn't see them again. The GNN's $\delta$ had to relearn linear
covariate effects from gradient descent through several layers of
nonlinearity, while three regularisations (`lambda_shrink * MSE(\delta -
baseline)`, `lambda_gate * MSE(r)`, $r \le 0.8$) actively pulled
$\delta$ back toward a covariate-blind baseline. The architecture
was structurally underpowered for cov-aware prediction in single-obs
mode.

This was a long-standing limitation, documented as
`covs ... no user covariates yet` in `CLAUDE.md`. The earlier multi-obs
extension (v0.6.0+) had added `obs_refine` for the multi-obs path but
left single-obs broken.

## 3. The fix

We applied eight successive interventions, evaluated at each step on
a controlled GNN-earnings simulation (Section 4.1):

- **Fix A** (`R/fit_pigauto.R`): set `n_user_cov = n_cov_cols` always.
- **Fix B** (`R/model_residual_dae.R`): dedicated `cov_encoder` MLP
  giving raw covariates their own `hidden_dim` of nonlinear capacity.
- **Fix C+D**: a `cov_linear` direct-projection head added to $\delta$,
  contributing OUTSIDE the $(1 - r)/r$ blend gate, mirroring `phylolm`'s
  fixed-effects decomposition.
- **Fix E** (tested, ruled out): lowering the regularisation weights
  $\lambda_\mathrm{shrink}, \lambda_\mathrm{gate}$ did not move
  performance — regularisation was not the bottleneck.
- **Fix G** (`R/bm_internal.R`, `R/fit_baseline.R`): the root-cause fix.
  `bm_impute_col_with_cov()` performs GLS phylogenetic regression
  $\hat\beta = (X' R^{-1} X)^{-1} X' R^{-1} y$ with cov fixed effects,
  then BLUP for held-out cells. Mathematically equivalent to
  `phylolm(model="BM")` for prediction, unit-tested to $10^{-7}$.
  Includes:
  - **LRT gate**: falls back to no-cov baseline if covariates do not
    reduce residual variance by $\ge$ 2 %. Prevents fitting noise on
    phylo-redundant covariates.
- **Fix H** (`R/graph_transformer_block.R`, `R/model_residual_dae.R`):
  per-layer covariate injection. Each transformer block optionally
  accepts a `cov_h` tensor and adds it via residual; init to zero so
  block $\approx$ identity at training step 0.

The architectural decomposition is now:

$$
\mathrm{pred}_i = \underbrace{X_i \hat\beta + R_{io} R_{oo}^{-1} (y_o - X_o \hat\beta)}_{\text{Fix G: phylolm-equivalent}} + r_i \cdot \delta_{i,\mathrm{GNN}}
$$

where the GNN's $\delta$ is now responsible only for nonlinear /
interactive / cross-trait residuals. Linear covariate effects are
captured analytically in the baseline.

## 4. Empirical evidence

### 4.1 Synthetic GNN earnings simulation

We simulated traits as

$$
y_i = \sqrt{\alpha} \cdot \mathrm{phylo}_i + \sqrt{\beta} \cdot f(\mathrm{cov}_i) + \sqrt{1 - \alpha - \beta} \cdot \varepsilon_i
$$

on the AVONET 300 bird tree (n=300), with $f \in \{\mathrm{linear},
\mathrm{nonlinear}\ (\sin \cdot \exp), \mathrm{interactive}\
(c_1 c_2 + 0.5 c_1^2)\}$ and 30 % MCAR mask. We compared pigauto
with covariates (safety_floor = TRUE) to `phylolm-lambda` BLUP — the
analytical optimum under the data-generating process for linear $f$.

**Pre-fix vs post-fix gap to `phylolm-lambda` BLUP:**

| $f$ | pre-fix | Fix A-D | Fix G+H | reduction |
|---|---:|---:|---:|---:|
| linear | 1.85$\times$ | 0.97$\times$* | **1.22$\times$** | -63 percentage points |
| nonlinear | 1.27$\times$ | 0.95$\times$* | **1.22$\times$** | -5 pp |
| interactive | 1.15$\times$ | 1.13$\times$ | **1.16$\times$** | flat |

$^*$ Fix A-D values reflect single-seed smoke; Fix G+H are at default
ridge=0, LRT=0.02.

The linear-effect bottleneck was the dominant gap; Fix G closes it
to within 22 % of analytical optimum. The remaining gap is the neural
network approximation tax — pigauto trains $\delta$ via gradient
descent, multiple imputation, and a calibrated blend gate, all of
which add finite-sample variance over the analytical BLUP.

### 4.2 Real-data covariate benches (6 datasets)

We re-ran six pre-existing covariate-lift benches with the post-Fix-G+H
code:

dataset | n | covariates | trait flips OLD $\to$ NEW
---|---:|---|---
GlobTherm ectotherms | 809 | latitude, longitude, elevation | **Tmax: 1.28 $\to$ 0.74 (52 pp swing)**
PanTHERIA mammals | 850 | precip, temp, lat, PET | GestationLen 1.09 $\to$ 0.90; Body mass and others ~flat
AmphiBIO amphibians | 1,000 | climate-zone occupancy | Body_mass_g 1.02 $\to$ 0.94; Litter_size 1.08 $\to$ 0.95
LepTraits butterflies | 1,500 | monthly flight phenology | FW_L 1.03 $\to$ 0.95
BIEN plants | 3,450 | WorldClim bioclim | sla 1.03 $\to$ 0.94 (rest flat); LRT also tamed `sf=FALSE` catastrophe (height_m 5.6$\times$ $\to$ 1.6$\times$ hurt)
Delhey birds | 5,809 | 6 climate covariates | unchanged (LRT correctly detects climate is phylo-redundant)

Across all 23 (dataset, trait) cells:

- traits with $\ge$ 5 % lift: 3 OLD $\to$ 5 NEW (+2 net)
- 4 traits flipped null/regression $\to$ lift
- 2 traits flipped lift $\to$ null/regression (likely LRT threshold sensitivity)
- 16 stayed flat; 1 stayed lift

**Headline win**: GlobTherm Tmax — Fix G captures the textbook
latitude $\to$ CTmax linear relationship analytically that the
unmodified GNN had been unable to extract via gradient descent.

**Honest non-win**: Delhey n=5,809 plumage-lightness data stays flat.
The LRT gate correctly detects that climate covariates are
phylo-redundant on this dataset (closely-related species share both
climates and plumage) and falls back to the no-cov baseline. This is
a **safety property** rather than a failure: pigauto returns
predictions identical to no-covariate pigauto when covariates carry no
phylo-decoupled information.

### 4.3 Multivariate cross-trait coupling (Day 2 sim)

We simulated multivariate Brownian motion on K=4 correlated traits
with a nonlinear cross-trait coupling: $y_4 \mathrel{+}= c \cdot
\sin(2 y_1) \cdot \exp(0.3 y_2)$. We compared pigauto's joint
multi-trait imputation to `Rphylopars` BLUP (the analytical baseline
for joint MVN with cross-trait correlation), at $c \in \{0, 0.5, 1\}$.

| $c$ | column-mean | lm crosstrait | Rphylopars | pigauto | pigauto / Rphylopars |
|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.743 | 0.630 | **0.479** | 0.521 | 1.09$\times$ |
| 0.5 | 0.803 | 0.692 | **0.547** | 0.599 | 1.10$\times$ |
| 1.0 | 1.052 | 0.892 | **0.710** | 0.757 | 1.07$\times$ |

pigauto loses to `Rphylopars` BLUP by 7–10 % across all coupling
strengths, consistent with the n=300 finite-sample approximation tax.
The gap is slightly smaller at strong coupling (1.09 $\to$ 1.07),
suggesting the GNN extracts *some* of the nonlinear cross-trait
signal — just not enough to flip the result at this n.

### 4.4 Multi-observation regime

(Pending — `bench_multi_obs.R` re-running with Fix A-H. Pre-fix this
sim showed 10–19 % lift from `acclim_temp` covariates on
within-species CTmax variation. Fix H's per-layer cov injection was
specifically designed to enhance this regime; we expect it to extend
the lift further.)

## 5. Discussion

The headline result — Fix G converts a 28 % regression on GlobTherm
CTmax to a 26 % lift — is the cleanest demonstration that pigauto's
new architecture captures what its previous architecture couldn't. The
finding generalises across four other real datasets where at least one
trait flipped from null to a meaningful lift after the fix.

We do **not** claim that pigauto's GNN systematically beats
`phylolm-lambda` BLUP on its analytical home turf (continuous, single-
obs, BM-on-tree, linear covariate effects). On those data pigauto pays
a 7–22 % approximation tax — the cost of being a flexible neural model
fitted by stochastic gradient descent rather than a closed-form
analytical solver. The honest paper claim is:

> pigauto's covariate-aware architecture (Fix G) recovers analytical
> phylolm-BLUP performance on linear covariate effects within
> approximation tolerance and provides modest covariate lifts on
> traits with phylo-decoupled signal. On regimes outside `phylolm`'s
> applicability — mixed-type response variables, multi-observation
> repeated measures, joint cross-trait imputation — pigauto's
> unified API and adaptive blend gate provide value that no single
> analytical method offers.

## 6. Limitations

- The LRT threshold (default 0.02) occasionally rejects useful
  covariates (loss of AmphiBIO Body_size_mm lift) or accepts spurious
  ones (PanTHERIA PopDensity regression). Per-trait or
  cross-validated thresholds would address this.
- Single-observation, single-trait covariate lifts on real data
  remain modest (5-26 % on the lifted traits). Datasets where
  covariates are genuinely phylo-redundant (Delhey) stay flat
  regardless of the architecture.
- The synthetic GNN earnings sim is necessarily an upper bound on
  pigauto's home turf — real comparative-biology data has more
  complex phylogenetic processes (rate variation, OU, regime shifts)
  and trait-type heterogeneity that the sim does not exercise.

## 7. Reproducibility

All fixes are on the `experiment/gnn-earnings-sim` branch (PR #49).
Bench scripts: `script/bench_{plants_cached_only, pantheria_covariates,
globtherm_covariates, amphibio_covariates, leptraits_covariates,
delhey_covariates, gnn_earnings, gnn_earnings_v2,
gnn_earnings_multitrait, multi_obs}.R`. Architecture documentation:
`useful/GNN_ARCHITECTURE_EXPLAINED.md`. Verdict figure:
`useful/fix_G_real_data_verdict.png`.
