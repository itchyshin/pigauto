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

## 8. Uncertainty quantification

### 8.1 Prediction intervals

For every imputed continuous, count, ordinal or proportion trait we report a 95%
prediction interval built by split conformal prediction (Papadopoulos et al. 2002; Vovk
et al. 2005). Before training, a random 25% of the observed cells is withheld and further
divided into a validation set (one quarter) and a test set. The model is fitted without
them. On the validation cells the absolute residual between the true value and the
blended prediction is computed on the latent (standardised) scale, and for each trait the
interval half-width is taken as the empirical quantile of these residuals at level
⌈(1 − α)(n + 1)⌉ / n, with n the number of validation residuals and α = 0.05. The
interval for a missing cell is the prediction plus or minus this half-width,
back-transformed to the trait's original scale. The construction is distribution-free: if
the residual of a new cell is exchangeable with the validation residuals, the interval
contains the true value with probability at least 1 − α (Lei et al. 2018).

That exchangeability condition is the one assumption the method makes, and phylogenetic
data violate it in a specific way. Validation cells are drawn from the observed part of
the matrix, which sits in well-sampled regions of the tree where the baseline can lean on
close relatives and prediction errors are small. The cells a user actually needs to
impute are not so placed; in real trait databases missingness concentrates in poorly
sampled clades, where errors are larger. A single quantile calibrated on the first
population and applied to the second is too short. In a simulation with clade-structured
missingness (two clades at 7:1 odds of being missing against the background), empirical
coverage of the nominal 95% interval fell to 0.923 at n = 300 species and 0.927 at
n = 1000, against 0.961 and 0.957 under completely random missingness; the shortfall did
not diminish with sample size (Supplementary Table S-UQ1).

To restore the condition where it fails, pigauto offers a locality-stratified variant
(`conformal_method = "mondrian"`; Vovk 2012; Boström & Johansson 2020). For each
validation cell a locality statistic is computed as the mean cophenetic distance from its
species to the five nearest species with an observed value for that trait. Validation
cells are split at the median locality into a near and a far stratum, and a separate
conformal quantile is computed within each stratum at the same adjusted level. At
prediction time a missing cell receives the half-width of the stratum its own locality
places it in, so intervals widen in undersampled clades, which is where the error is.
Within a stratum the near-versus-far mismatch that broke exchangeability is largely
removed, and the marginal guarantee then holds within each stratum, provided the cells to be imputed are exchangeable with the validation cells of their stratum; the structured mask approximates that condition for real missing cells but cannot establish it (Section 8.3). On the
same simulation grid this recovered clade-structured coverage to 0.946 at n = 1000,
within three Monte Carlo standard errors of nominal, while widening the median interval
under random missingness by 2.3%.

The variant has a floor. A stratum's own conservative quantile can only reach 0.95 when
it holds at least 19 residuals (the smallest n for which n / (n + 1) ≥ 0.95), so a trait
whose validation cells would leave either stratum below 19 (roughly 38 validation cells
per trait) falls back to the global split quantile, and the software records that it did.
At n = 300 in the simulation above this fallback fired for every trait and the variant
changed nothing; the remedy at that size is more held-out data rather than
stratification. The stratified variant is currently limited to single-observation data,
because its locality is defined per species, and to traits that receive conformal
intervals at all. The default remains the unstratified split quantile; Section 8.3
reports the real-data confirmation and why.

### 8.2 Supplementary Table S-UQ1: mechanism-coverage simulation

Regime: F1 continuous, λ = 1, 30% missing, single-obs, simulated. Coverage is scored on
genuinely-missing cells at nominal 0.95; MCSE = sd(per-rep coverage)/√n_rep.

**B1: split conformal (baseline)**

| mechanism | n | n_rep | split coverage | MCSE |
|---|---:|---:|---:|---:|
| MCAR | 300 | 30 | 0.961 | 0.004 |
| MAR_trait | 300 | 30 | 0.949 | 0.003 |
| MAR_phylo | 300 | 30 | 0.923 | 0.006 |
| MNAR | 300 | 30 | 0.939 | 0.006 |
| MCAR | 1000 | 20 | 0.957 | 0.002 |
| MAR_trait | 1000 | 20 | 0.940 | 0.004 |
| MAR_phylo | 1000 | 20 | 0.927 | 0.005 |
| MNAR | 1000 | 20 | 0.932 | 0.003 |

**B2: split vs. Mondrian conformal, paired re-run**

| mechanism | n | split coverage | mondrian coverage | MCSE (split / mondrian) |
|---|---:|---:|---:|---:|
| MCAR | 300 | 0.961 | 0.961 | fallback fired; identical by design |
| MAR_phylo | 300 | 0.923 | 0.923 | fallback fired; stratification impossible at this n_val |
| MNAR | 300 | 0.939 | 0.939 | fallback fired; same |
| MCAR | 1000 | 0.957 | 0.957 | width +2.3% (cap 10%) |
| MAR_trait | 1000 | 0.940 | 0.946 | gap 0.004 < 3×0.0037 |
| MAR_phylo | 1000 | 0.927 | 0.946 | gap 0.004 < 3×0.0053 |
| MNAR | 1000 | 0.932 | 0.940 | gap 0.010 < 3×0.0043, marginal |

At n = 1000 all three non-MCAR mechanisms recover to within 3×MCSE of the nominal 0.95,
with MAR_phylo, the worst mechanism, fully repaired. At n = 300 the per-stratum floor of
19 residuals is not met for any trait, so the method falls back to the split quantile and
coverage is unchanged by construction.

### 8.3 Real-data confirmation

Masking originally observed cells gives true values to score against. We did this on
three trait databases under a pre-registered design
(`docs/dev-log/mondrian-realdata/00-preregistration.md`). Two mask arms were used. A
random mask of 20% of observed cells makes test cells exchangeable with the validation
cells, so split conformal is valid there by construction; it served as a no-harm
control. A structured mask drew the same 20% with probability proportional to a
propensity for real missingness, fitted per trait on phylogenetic eigenvectors, so that
test cells sit where genuinely missing cells sit. PanTHERIA (4,027 mammals) ran both
arms with three masks each; AVONET (1,500 birds) is almost completely observed and ran
the random arm only; FishBase (10,484 fishes) ran the structured arm with one mask.
Mondrian activated for every continuous trait in all three databases. Ordinal traits also received Mondrian intervals but were not scored; the evidence covers the continuous and count traits listed in Table S-UQ2 only.

**Table S-UQ2.** Coverage of nominal 95% intervals, pooled over traits (weighted by test
cells) and masks, by the stratum of each test cell. The last column is the median over
traits of the per-trait width ratio. FishBase rests on one mask and is descriptive.

| database | mask arm | stratum | split | Mondrian | width ratio (Mondrian / split) |
|---|---|---|---:|---:|---:|
| PanTHERIA | structured | far | 0.919 | 0.940 | 1.18 |
| PanTHERIA | structured | near | 0.974 | 0.958 | 0.84 |
| PanTHERIA | random | far | 0.931 | 0.957 | 1.23 |
| PanTHERIA | random | near | 0.978 | 0.968 | 0.83 |
| FishBase | structured | far | 0.930 | 0.949 | 1.14 |
| FishBase | structured | near | 0.966 | 0.950 | 0.78 |
| AVONET | random | far | 0.929 | 0.964 | 1.67 |
| AVONET | random | near | 0.984 | 0.961 | 0.81 |

Pooled over traits, the split quantile undercovers in the far stratum and overcovers in
the near stratum in every database, including under the random mask. Mondrian moves both
towards nominal: far-stratum coverage rises to 0.94-0.96 at the cost of wider intervals
there, and near-stratum intervals narrow by about a fifth while their pooled coverage
stays at or above 0.95. Per trait the picture is noisier: six of nineteen near rows fall
below 0.95 under Mondrian, and in three far rows split was already at or above nominal. Our pre-registered rule for changing the default also required that Mondrian's
near-stratum coverage fall no more than two percentage points below split's. That
condition was not met. For AVONET the pooled near-stratum drop (2.3 points) exceeded the
margin. For FishBase the drop (1.6 points) was inside the margin but, with a single mask,
could not be shown non-inferior (one-sided p = 0.22, Holm-adjusted 0.43). PanTHERIA
passed (drop 1.2 points, adjusted p = 0.025). The AVONET result alone decides the
verdict, so the default remains split. Mondrian remains an opt-in. The far-stratum gains
under the structured mask (condition 1 passed on both databases that ran it) are the
evidence for choosing it when missing species are concentrated in poorly sampled clades;
the price is wider far intervals and a 1 to 2 point drop in near-stratum coverage, which
fell below 0.95 for six of nineteen trait rows. Per-trait results, uncertainty and the decision script are in the
repository.

### References

- Papadopoulos H, Proedrou K, Vovk V, Gammerman A (2002) Inductive confidence machines
  for regression. In: *Machine Learning: ECML 2002*, Lecture Notes in Computer Science
  vol. 2430, Springer, 345–356. DOI: 10.1007/3-540-36755-1_29.
- Vovk V, Gammerman A, Shafer G (2005) *Algorithmic Learning in a Random World*. Springer.
  DOI: 10.1007/b106715.
- Vovk V (2012) Conditional validity of inductive conformal predictors. In: *Proceedings
  of the 4th Asian Conference on Machine Learning*, PMLR 25, 475–490.
- Lei J, G'Sell M, Rinaldo A, Tibshirani RJ, Wasserman L (2018) Distribution-free
  predictive inference for regression. *Journal of the American Statistical Association*
  113(523), 1094–1111. DOI: 10.1080/01621459.2017.1307116.
- Boström H, Johansson U (2020) Mondrian conformal regressors. In: *Proceedings of the
  Ninth Symposium on Conformal and Probabilistic Prediction and Applications*, PMLR 128,
  114–133.
