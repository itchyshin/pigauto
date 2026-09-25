# Diagnosis of the two campaign failures (posterior MI, campaign 69670d4)

2026-09-24. Branch `arc/mi-posterior`. Campaign: 24 regimes x 200 reps at code commit
69670d44f96e5934d194abac8c6653ee82e11811 (`results_tables.md`, `sim_summary.csv`). This note
diagnoses (1) the negative paired slope bias in the lambda = 1 regimes and (2) the excess
non-convergence in regimes 5, 21 and 23. No file under `R/`, `tests/`, `NEWS.md` or `GATES.md` was
changed.

Provenance. Every model fit below ran on Totoro from the frozen campaign copy
(`/home/snakagaw/pigauto_mi_posterior/69670d44f9/code`, package loaded with
`devtools::load_all()` from there), driven by the scripts in `script/mi_gls/diag/` at commit
3f0e0fbd3c (archived to `/home/snakagaw/pigauto_mi_posterior/diag_3f0e0fbd3c/code`, SHA in
`code/SHA`). The one exception is the prior-sensitivity arm, which ran from a scratch copy of the
frozen tree with one line changed (`evidence/diagnosis/PRIOR_CHANGE.diff`). BLAS and OpenMP threads
were pinned to 1 in the launching shell; at most 60 cores were used (briefly 64 for about 3 minutes
while two queues overlapped). Wall time about 55 minutes. Every table here is printed by
`script/mi_gls/diag/analyse_diag.R` into `evidence/diagnosis/diag_tables.md` from the CSVs in the
same folder. "Measured" marks a number read from those files; "Inference" marks a reading of them.

## Summary

**Failure 1 (lambda = 1 slope bias).** Not mainly a property of the estimand, and not a sampler bug.
The dominant cause, for phylolm, is a mismatch between the data-generating process (DGP) of regimes 1
to 16 and the sampler's model: those regimes simulate on a non-ultrametric `ape::rtree()` with the
raw covariance `vcv(tree)`, so tip variances differ up to about 7-fold, while the sampler (by
design) uses the correlation matrix `cov2cor(vcv(tree))`, which forces equal tip variances. A
second, smaller component comes from the posterior's residual covariance Sigma_E at lambda near 1,
where the data say little about Sigma_E and its prior pulls the residual correlation towards 0. At
n = 300 with x only missing, a third component of about -0.003 remains even with exact imputations.

| hypothesis | verdict | decisive numbers (phylolm; regime 1, n = 300 / regime 3, n = 1000) |
|---|---|---|
| H2: exact imputations from the true parameters still give the paired bias at finite n | Refuted as a main cause; a small n = 300 component is measured | oracle_true -0.0029 (MCSE 0.0012) / -0.0001 (0.0014), against campaign -0.0259 (0.0016) / -0.0175 (0.0019) |
| New, H4: the regime 1 to 16 DGP lies outside the sampler's model family (non-ultrametric tree; covariance vs correlation) | Supported; the dominant phylolm component | package draws at the model's best attainable (KL) parameters minus oracle_true: -0.0163 (0.0005) / -0.0137 (0.0007); on in-model twin data the full sampler's regime 3 bias is -0.0025 (0.0029) against -0.0173 (0.0029) |
| H1: the posterior point estimate of the covariances is off | Partly supported, for Sigma_E only | lambda and phylogenetic correlation match the KL pseudo-truth (0.994 vs 0.993; 0.708 vs 0.700); residual correlation 0.17 / 0.27 against 0.70; campaign minus the KL oracle -0.0068 (0.0010) / -0.0037 (0.0015); a 100-fold smaller Sigma_E prior scale moves the slope by +0.0033 (0.0019) / +0.0034 (0.0012) but made 2 of 24 fits fail numerically |
| H3: pipeline (decode, z-scoring, latent scale, flat mu) | Refuted for decode, z-scoring and the latent scale; flat mu not isolated | observed cells unchanged in all 138 re-runs (max change 0), no log transform, no NA; the same pipeline on in-model twin data removes 86% of the phylolm bias |

Recommendation: no change to the sampler's algorithm or defaults on this evidence. The main lever
is a harness decision about regimes 1 to 16 (below). A smaller Sigma_E prior scale would remove
about 0.003 to 0.005 of the bias, but it is a default change (campaign re-run) and failed
numerically in 2 of 24 test fits, so it is not recommended as is.

**Failure 2 (non-convergence in regimes 5, 21, 23).** All 32 non-converged fits fail on bulk ESS
(190 to 395), none on R-hat (max split R-hat 1.005 to 1.028, threshold 1.05). The low-ESS set is
the whole weakly identified block: in regime 5 (lambda = 1) both lambdas and all of Sigma_E (with
Sigma_P[1,2] just under 400), in regimes 21 and 23 (lambda 0.7, residual correlation 0) five to all
eight of the eight monitored quantities; the single lowest is
Sigma_E[1,2] in 6 of the 8 fits re-run at 1x. At 2x chain length 30 of 32 converge (the two that do
not are regime 21 reps 20 and 68, min ESS 345 and 375 on lambda[1]); at 4x all 32 converge (min ESS
472 to 2,283, max R-hat at most 1.011). The pooled slopes do not move systematically (mean signed
change 4x minus 1x: phylolm +0.0029 (0.0026), gls -0.0007 (0.0042)); the mean absolute change
(0.012, 0.018) equals the change between the 2x and 4x runs (0.011, 0.017), that is, Monte Carlo
noise of m = 20 draws, against a pooled SE of 0.080 and 0.105. This is an efficiency problem, not a
bias. Meeting the 2% rule at default settings needs a longer default chain (a default change,
campaign re-run) or a harness change.

## Failure 1: negative paired slope bias at lambda = 1

### What the campaign measured

Paired bias = mean over reps of (pooled MI slope minus complete-data slope in the same rep).
posterior_full, all 200 reps (`results_tables.md`): phylolm -0.026 (MCSE 0.002) in regime 1,
-0.019 in regime 3, -0.008 to -0.021 in the other lambda = 1 regimes; gls -0.002 to -0.013. The
plug-in arm (covariances fixed at their posterior means) has the same bias as the full posterior in
all four both-missing lambda = 1 regimes (converged reps, `campaign_cells.csv`: phylolm full vs
plug-in -0.0205 vs -0.0218 in regime 9, -0.0162 vs -0.0162 in 11, -0.0130 vs -0.0143 in 13, -0.0082
vs -0.0076 in 15), so parameter uncertainty is not the cause.

### A fact about the DGP that the hypotheses did not name

Regimes 1 to 16 reuse the v1 DGP (`dgp_v2.R`, verbatim from `01_cell.R`): `tree <- ape::rtree(n)`,
`V_sim <- vcv(sim_tree) / max(vcv(sim_tree))`, `vec(Y) ~ N(0, Sig %x% V_sim)` with
`Sig = [[1, 0.7], [0.7, 1]]`. `ape::rtree()` returns a non-ultrametric tree, so diag(V_sim), each
tip's root-to-tip depth over the maximum depth, is not constant. Measured (Table 1b): mean
diag(V_sim) 0.54 to 0.56 and minimum 0.13 to 0.15, against a maximum of 1.

The sampler's model (`design.md` section 1; the header of `R/mi_posterior.R`) is
`vec(Y) ~ N(1 mu', Sigma_P %x% R + Sigma_E %x% I_n)` with `R = cov2cor(vcv(tree))`, which gives
every tip the same marginal variance. For a non-ultrametric tree no (Sigma_P, Sigma_E, mu)
reproduces `Sig %x% V_sim`: regimes 1 to 16 lie outside the sampler's model family. Regimes 17 to 24
simulate from `R` itself (the Kronecker branch of `dgp_v2.R`), inside the family; their bias is
within 0.006 of 0 except the two n = 300 regimes with residual correlation 0 (phylolm -0.015 in 21
and -0.010 in 23, MCSE 0.005). The two analysis models differ in the same way: `gls(corBrownian)` uses the correlation
matrix (`ape:::corMatrix.corBrownian()` calls `vcv.phylo(tree, corr = TRUE)`), the same `R` as the
sampler, while `phylolm(model = "lambda")` uses the covariance `vcv(tree)`, which matches the DGP at
lambda = 1.

So "the true (Sigma_P, Sigma_E, mu) of the DGP" does not exist in the sampler's parametrisation for
regimes 1 to 16, and the requested oracle was split into three arms.

### Experiment

Same regimes, reps and seeds as the campaign (the DGP is seed-deterministic; the complete-data
slopes recomputed by the oracle script equal the campaign's to < 1e-8 in every rep used), m = 20,
downstream fits and pooling copied verbatim from `01_cell_v2.R` (`script/mi_gls/diag/oracle_cell.R`):

- `oracle_true`: exact draws of the missing cells from their conditional under the DGP's own
  covariance `Sig %x% V_sim` (dense conditional, mu = 0 known). Tests H2.
- `oracle_kl`: the package's fixed-parameter draws (`.mip_fixed_draws()`, the function behind the
  G2 exactness gate) under the sampler's model at its best attainable parameters,
  Sigma_P = a* Sig and Sigma_E = b* Sig + 1e-6 I, where (a*, b*) minimise
  KL(N(0, V_sim) || N(0, a R + b I)) for that rep's tree; mu = 0 fixed. The DGP is separable and the
  model family is closed under linear maps of the traits, so this is the KL projection of the DGP
  onto the model, the parameters a large-n posterior would settle on. It needs no estimation and
  no prior.
- `oracle_nominal`: the same package draws at the literal DGP values in the sampler's
  parametrisation (lambda = 1: Sigma_P = Sig, Sigma_E = 1e-6 I, as G2b does; lambda = 0.5:
  Sigma_P = Sigma_E = 0.5 Sig), mu = 0. Reported because it was requested; it is not a fair oracle
  here (unit tip variance against a mean true tip variance of 0.55).

Three sampler arms re-run the campaign's `multi_impute(draws_method = "posterior")` call on the
frozen code (`script/mi_gls/diag/rerun_cell.R`):

- `h1`: regimes 1 and 3, reps 1 to 12, re-run to record the posterior summaries the campaign output
  did not keep. They reproduce the campaign exactly (pooled slopes and min ESS identical).
- `twin`: regime 3, reps 1 to 20, with the truth replaced by its in-model twin: each tip's row
  divided by sqrt(diag(V_sim)[i]), so the data are exactly `N(0, Sig %x% R)` (same seed, tree,
  noise and masks). Same sampler, priors and pipeline, on data inside the model family.
- `prior`: regimes 1 and 3, reps 1 to 12, one prior change in a scratch copy of the package:
  `S_E = 1e-4 diag(obs var)` instead of `0.01 diag(obs var)` (chain settings unchanged). Run because
  Table 2 showed the residual correlation shrunk towards 0.

Rep counts differ by arm because of the one-hour budget (regime 1 oracle: all 200 reps; other
regimes: the first 24 to 60 reps; regimes 7, 15 and 4, and the lambda = 0.5 regimes other than 2,
were not run).

### Table 1. Paired slope bias by arm (mean over reps, MCSE in brackets)

"Campaign" is the campaign's posterior_full paired bias on the same reps as the oracle arms. The last
two columns are paired differences within rep.

| regime | analysis | reps | campaign posterior | oracle_true (DGP covariance) | oracle_kl (sampler model, KL pseudo-true) | oracle_nominal (sampler model, nominal) | campaign minus oracle_kl | oracle_kl minus oracle_true |
|---|---|---|---|---|---|---|---|---|
| 1: lambda 1, n 300, MCAR, x | gls | 200 | -0.0134 (0.0015) | -0.0030 (0.0012) | -0.0056 (0.0012) | -0.0780 (0.0016) | -0.0078 (0.0010) | -0.0026 (0.0004) |
| 1: lambda 1, n 300, MCAR, x | phylolm | 200 | -0.0259 (0.0016) | -0.0029 (0.0012) | -0.0192 (0.0012) | -0.0967 (0.0018) | -0.0068 (0.0010) | -0.0163 (0.0005) |
| 3: lambda 1, n 1000, MCAR, x | gls | 40 | -0.0054 (0.0020) | 0.0003 (0.0015) | -0.0013 (0.0016) | -0.0791 (0.0026) | -0.0041 (0.0014) | -0.0016 (0.0004) |
| 3: lambda 1, n 1000, MCAR, x | phylolm | 40 | -0.0175 (0.0019) | -0.0001 (0.0014) | -0.0138 (0.0015) | -0.0964 (0.0029) | -0.0037 (0.0015) | -0.0137 (0.0007) |
| 5: lambda 1, n 300, MAR, x | gls | 60 | -0.0017 (0.0031) | -0.0011 (0.0026) | 0.0014 (0.0028) | -0.0753 (0.0030) | -0.0031 (0.0015) | 0.0025 (0.0013) |
| 5: lambda 1, n 300, MAR, x | phylolm | 60 | -0.0132 (0.0028) | -0.0013 (0.0024) | -0.0111 (0.0025) | -0.0892 (0.0030) | -0.0021 (0.0015) | -0.0098 (0.0013) |
| 9: lambda 1, n 300, MCAR, both | gls | 60 | -0.0236 (0.0047) | -0.0025 (0.0034) | -0.0027 (0.0037) | -0.0469 (0.0039) | -0.0209 (0.0038) | -0.0002 (0.0012) |
| 9: lambda 1, n 300, MCAR, both | phylolm | 60 | -0.0305 (0.0047) | -0.0019 (0.0031) | -0.0120 (0.0033) | -0.0597 (0.0037) | -0.0185 (0.0039) | -0.0101 (0.0014) |
| 11: lambda 1, n 1000, MCAR, both | gls | 24 | -0.0099 (0.0043) | -0.0007 (0.0029) | -0.0019 (0.0028) | -0.0526 (0.0030) | -0.0079 (0.0033) | -0.0013 (0.0012) |
| 11: lambda 1, n 1000, MCAR, both | phylolm | 24 | -0.0158 (0.0043) | 0.0003 (0.0029) | -0.0092 (0.0029) | -0.0625 (0.0031) | -0.0066 (0.0034) | -0.0095 (0.0012) |
| 13: lambda 1, n 300, MAR, both | gls | 60 | -0.0102 (0.0052) | 0.0012 (0.0029) | -0.0008 (0.0027) | -0.0434 (0.0031) | -0.0094 (0.0040) | -0.0020 (0.0015) |
| 13: lambda 1, n 300, MAR, both | phylolm | 60 | -0.0150 (0.0051) | -0.0003 (0.0025) | -0.0080 (0.0025) | -0.0517 (0.0033) | -0.0070 (0.0040) | -0.0077 (0.0014) |
| 2: lambda 0.5, n 300, MCAR, x | gls | 60 | 0.0000 (0.0042) | -0.0076 (0.0036) | -0.0047 (0.0037) | -0.0852 (0.0044) | 0.0047 (0.0026) | 0.0029 (0.0012) |
| 2: lambda 0.5, n 300, MCAR, x | phylolm | 60 | -0.0142 (0.0026) | -0.0008 (0.0024) | -0.0151 (0.0023) | -0.0989 (0.0034) | 0.0009 (0.0018) | -0.0143 (0.0010) |

Measured. At lambda = 1 `oracle_true` is within 2.5 MCSE of 0 in every regime (largest -0.0030
(0.0012), regime 1 gls). The package's own draws at the model's best attainable parameters
(`oracle_kl`) are biased for phylolm in every lambda = 1 regime (-0.0080 to -0.0192). The gap
`oracle_kl minus oracle_true`, the cost of the model family alone, is -0.0077 to -0.0163 for phylolm
and -0.0026 to +0.0025 for gls. `oracle_nominal` is far more biased (-0.043 to -0.099), as expected
from draws over-dispersed by a factor of about 1 / 0.55.

### Table 1b. KL pseudo-true parameters and tip variances of the DGP (raw scale, mean over reps)

| regime | reps | a* (Sigma_P scale) | b* (Sigma_E scale) | lambda* = a*/(a*+b*) | mean diag(V_sim) | min diag(V_sim) | optim failures |
|---|---|---|---|---|---|---|---|
| 1 | 200 | 0.537 | 0.0035 | 0.993 | 0.557 | 0.149 | 0 |
| 2 | 60 | 0.247 | 0.2808 | 0.467 | 0.544 | 0.144 | 0 |
| 3 | 40 | 0.521 | 0.0024 | 0.995 | 0.540 | 0.128 | 0 |
| 5 | 60 | 0.531 | 0.0034 | 0.993 | 0.552 | 0.143 | 0 |
| 9 | 60 | 0.540 | 0.0034 | 0.994 | 0.561 | 0.139 | 0 |
| 11 | 24 | 0.520 | 0.0023 | 0.996 | 0.537 | 0.137 | 0 |
| 13 | 60 | 0.534 | 0.0034 | 0.994 | 0.555 | 0.136 | 0 |

At lambda = 1 the model's best fit puts almost all variance in Sigma_P (lambda* 0.993 to 0.996) at
a scale close to the mean tip variance (a* 0.52 to 0.54 against 0.54 to 0.56). It cannot represent
the spread of tip variances (0.13 to 1).

### Table 2. Posterior summaries (H1)

Posterior means per fit on the sampler's latent (z-scored) scale, averaged over reps, SD over reps in
brackets. lambda and correlations are scale free, so they compare directly with Table 1b and with
the DGP's 0.7.

| arm | regime | reps | lambda_x | lambda_y | corr_P | corr_E | beta_P = SP12/SP22 | beta_E = SE12/SE22 | converged |
|---|---|---|---|---|---|---|---|---|---|
| h1 | 1 | 12 | 0.994 (0.001) | 0.992 (0.005) | 0.708 (0.040) | 0.168 (0.113) | 0.738 (0.121) | 0.178 (0.145) | 12/12 |
| h1 | 3 | 12 | 0.996 (0.002) | 0.997 (0.001) | 0.704 (0.010) | 0.266 (0.107) | 0.705 (0.132) | 0.306 (0.141) | 12/12 |
| prior | 1 | 11 | 1.000 (0.000) | 0.999 (0.002) | 0.694 (0.044) | -0.015 (0.037) | 0.722 (0.126) | -0.015 (0.055) | 10/11 |
| prior | 3 | 11 | 1.000 (0.000) | 1.000 (0.000) | 0.698 (0.013) | 0.002 (0.057) | 0.693 (0.133) | -0.002 (0.068) | 11/11 |
| twin | 3 | 20 | 0.997 (0.001) | 0.997 (0.001) | 0.703 (0.011) | 0.231 (0.120) | 0.696 (0.106) | 0.263 (0.145) | 20/20 |

Measured. Posterior lambda (0.992 to 0.997) matches lambda* (0.993, 0.995) and the phylogenetic
correlation (0.708, 0.704) matches 0.7. The residual correlation does not: 0.168 and 0.266 against
0.7 in the KL pseudo-truth (Sigma_E there is b* Sig). Inference: Sigma_E is about 0.5% to 0.8% of
each trait's latent variance here, so the data carry little information about it, and its
inverse-Wishart prior IW(3, 0.01 diag(obs var)), whose marginal prior on the correlation is uniform
on (-1, 1), pulls the correlation towards 0. This is the H1 mechanism the brief named, located in
Sigma_E, not in Sigma_P or lambda. The in-model twin shows the same shrinkage (0.231), as expected,
since its true Sigma_E is 0.

### Prior-sensitivity arm

Pooled slope with `S_E = 1e-4 diag(obs var)` minus the default-prior slope of the same rep:

| regime | reps | gls | phylolm | corr_E default | corr_E prior arm | Sigma_E[1,1] default | Sigma_E[1,1] prior arm |
|---|---|---|---|---|---|---|---|
| 1 | 11 | 0.0046 (0.0013) | 0.0033 (0.0019) | 0.169 | -0.015 | 7.80e-03 | 3.81e-04 |
| 3 | 11 | 0.0044 (0.0012) | 0.0034 (0.0012) | 0.270 | 0.002 | 4.78e-03 | 2.56e-04 |

Measured. The smaller prior scale shrinks Sigma_E about 20-fold (below the KL value), drives
lambda to 1.000 and the residual correlation to about 0, and raises the pooled slope by +0.0033 to
+0.0046. That is 50% to 110% of the "campaign minus oracle_kl" gap in these two regimes (Table 1:
regime 1 phylolm -0.0068, gls -0.0078; regime 3 -0.0037, -0.0041), but only 13% to 22% of the phylolm
bias. Two of 24 fits failed (regime 1 rep 3, regime 3 rep 10). A serial re-run of regime 1 rep 3
stops with `.updateCHMfactor(...): leading principal minor of order 1198 is not positive`
(`evidence/diagnosis/logs/logs/prior_serialcheck_r1_3.log`): the sparse Cholesky update fails once
Sigma_E may become nearly singular. One further fit of the 22 did not converge.

### Table 3. The sampler on in-model twin data (regime 3), and the prior arm

| arm | regime | analysis | reps | paired bias (MCSE) | campaign bias, same reps |
|---|---|---|---|---|---|
| h1 | 1 | gls | 12 | -0.0135 (0.0060) | -0.0135 (0.0060) |
| h1 | 1 | phylolm | 12 | -0.0248 (0.0066) | -0.0248 (0.0066) |
| h1 | 3 | gls | 12 | -0.0048 (0.0039) | -0.0048 (0.0039) |
| h1 | 3 | phylolm | 12 | -0.0168 (0.0038) | -0.0168 (0.0038) |
| twin | 3 | gls | 20 | -0.0045 (0.0029) | -0.0057 (0.0031) |
| twin | 3 | phylolm | 20 | -0.0025 (0.0029) | -0.0173 (0.0029) |
| prior | 1 | gls | 11 | -0.0073 (0.0058) | -0.0119 (0.0063) |
| prior | 1 | phylolm | 11 | -0.0207 (0.0068) | -0.0241 (0.0072) |
| prior | 3 | gls | 11 | 0.0003 (0.0038) | -0.0041 (0.0042) |
| prior | 3 | phylolm | 11 | -0.0124 (0.0039) | -0.0158 (0.0040) |

Measured. On in-model data the full sampler's phylolm bias is -0.0025 (0.0029), against -0.0173
(0.0029) for the campaign on the same 20 reps (the twin removes 86% of it). The gls bias does not
change (-0.0045 against -0.0057); at n = 1000 the gls bias is small and of the size the Sigma_E
prior arm moves (+0.0044). (The twin's complete-data slope differs from the campaign's, so the
comparison is between paired biases, not slopes.)

### Reading the components (inference)

Campaign bias = oracle_true + (oracle_kl minus oracle_true) + (campaign minus oracle_kl):

- Finite n with exact imputations (H2): about -0.003 at n = 300 with only x missing; within MCSE of
  0 in every other lambda = 1 regime examined.
- Model family, no estimation, no prior (H4): the dominant phylolm component (-0.008 to -0.016),
  small for gls (-0.003 to +0.003). This fits the two analysis models: phylolm weights tips by the
  DGP's own covariance, in which low-variance tips carry the most weight, and the sampler
  over-disperses exactly those tips (it gives every tip the same variance); gls(corBrownian) uses the
  sampler's own correlation matrix, so it is nearly congenial with the imputation model.
- The sampler's estimates versus the KL values (campaign minus oracle_kl): -0.002 to -0.009 in
  regimes 1, 3, 5, 11 and 13, of which roughly half to all is the Sigma_E prior in regimes 1 and 3.
  It is larger in regime 9 (n = 300, both traits missing: gls -0.021 (0.004), phylolm -0.019
  (0.004), 60 reps). The plug-in arm's equality with the full posterior places this component in the
  posterior MEAN of the covariances, not their uncertainty. Whether the Sigma_E prior explains all of
  the regime 9 gap is unresolved: the prior arm ran only on regimes 1 and 3, and the flat-mu versus
  fixed-mu difference between the sampler and oracle_kl was not isolated.

lambda = 0.5 (regime 2 only, 60 reps): the phylolm bias (-0.0142) is all model-family component
(oracle_kl minus oracle_true -0.0143 (0.0010)), and oracle_true is 0. The gls numbers there
(campaign 0.0000, oracle_true -0.0076 (0.0036)) are not decisive, and the positive gls bias of the
other lambda = 0.5 regimes (+0.006 to +0.017) was not examined. The lambda = 0.5 DGP is further
outside the model: Pagel's transform of a non-ultrametric tree adds a tip-specific nugget
(1 - lambda) diag(vcv(tree)), not (1 - lambda) I.

### Recommendation for failure 1

1. No change to the sampler's algorithm or defaults is indicated for the dominant component. The
   sampler does what its model (`design.md` section 1) says; regimes 1 to 16 simulate from a
   different model. For ultrametric (dated) trees `vcv(tree)` is proportional to
   `cov2cor(vcv(tree))`, so the issue does not arise; it matters only for non-ultrametric trees.
2. A harness decision for Shinichi, one of:
   - (a) Keep regimes 1 to 16 as they are and report them as a stress test of a misspecified
     imputation model (tip-variance heterogeneity from a non-ultrametric tree), and gate G6 rule 1
     (bias) only on in-model regimes. No re-run; a gate change, which is his call.
   - (b) Make regimes 1 to 16 in-model, either as the twin does (divide each tip's row by
     sqrt(diag(V_sim)), equivalently simulate from `cov2cor(V_sim)` as regimes 17 to 24 already do)
     or with ultrametric trees (for example `ape::rcoal()`). Expected effect, measured on the twin:
     regime 3 phylolm bias from -0.017 to about -0.003 (not measured for regime 1). What remains is
     the smaller Sigma_E and n = 300 components, about -0.003 to -0.008 in regimes 1 and 3 and
     possibly up to -0.02 in regime 9 (unresolved). Cost: re-run regimes 1 to 16 (3,200 fits at
     about 250 to 500 s each, roughly 220 to 450 core-hours). A harness change, not a package
     change.
3. Not recommended on this evidence: switching the sampler to the covariance `vcv(tree)` instead of
   the correlation matrix (dropping the tip scaling D in Qc). It would remove the model-family
   component for non-ultrametric trees, but it changes the algorithm, departs from pigauto's
   `R = cov2cor(vcv(tree))` convention used by every other imputation path, and needs a full
   campaign re-run.
4. The Sigma_E prior: a 100-fold smaller `S_E` reduced the bias by 0.003 to 0.005 in the two regimes
   tested, but it is a default change (campaign re-run), it was not tested where Sigma_E is truly
   non-zero (regimes 17 to 24), and it made 2 of 24 fits fail in the sparse Cholesky update. Not
   recommended as is. If the residual-correlation shrinkage is to be addressed, an untested
   alternative is an `S_E` whose off-diagonal follows the observed trait covariance rather than a
   diagonal matrix, so that the prior does not centre the residual correlation on 0. It would need
   its own numerical-stability check and a campaign re-run.

## Failure 2: non-convergence above 2% in regimes 5, 21 and 23

Measured from the campaign output (`evidence/diagnosis/campaign_cells.csv`): the non-converged fits
are regime 5 reps 59, 65, 120, 182, 183; regime 21 reps 15, 20, 28, 68, 73, 146, 157, 175, 182, 188,
190; regime 23 reps 5, 6, 37, 59, 76, 79, 87, 88, 111, 144, 153, 170, 174, 186, 196, 198 (5, 11 and
16 of 200). Every one fails on ESS only: min bulk ESS 190 to 395 (threshold > 400), max split R-hat
1.005 to 1.028 (threshold < 1.05). The campaign output kept only the max R-hat and the min ESS, so
8 of the 32 (the lowest ESS in each regime) were re-run at the campaign length to get per-parameter
diagnostics; the re-runs reproduce the campaign exactly (slopes and diagnostics identical). All 32
were re-run at 2x (burn-in 2,000, 10,000 sweeps) and 4x (burn-in 4,000, 20,000 sweeps), same seed,
keep_draws 1,000, with the 4 chains of each fit run in parallel (the chains are identical to a
serial run; see the note in `rerun_cell.R`; the plug-in draws are not, and are not used here).

Per-parameter bulk ESS at 1x (parameters below 400; from `rerun_diag_long.csv`):

| regime | rep | parameters with ESS <= 400 (of 8 monitored) |
|---|---|---|
| 5 | 182 | 6: lambda[1] 327, lambda[2] 320, Sigma_E[1,1] 326, Sigma_E[1,2] 293, Sigma_E[2,2] 319, Sigma_P[1,2] 382 |
| 5 | 183 | 6: lambda[1] 236, lambda[2] 235, Sigma_E[1,1] 240, Sigma_E[1,2] 227, Sigma_E[2,2] 235, Sigma_P[1,2] 398 |
| 21 | 20 | 5: lambda[1] 227, Sigma_E[1,1] 231, Sigma_E[1,2] 318, Sigma_P[1,1] 247, Sigma_P[1,2] 284 |
| 21 | 68 | 7: all but Sigma_E[2,2]; lowest Sigma_P[1,2] 227 |
| 21 | 175 | 8: all; lowest Sigma_E[1,2] 270 |
| 23 | 174 | 8: all; lowest Sigma_E[1,2] 226 |
| 23 | 186 | 8: all; lowest Sigma_E[1,2] 209 |
| 23 | 196 | 8: all; lowest Sigma_E[1,2] 190 |

Inference: in regime 5 (lambda = 1) the slow direction is the split of variance into a Sigma_E near
zero (lambda and all of Sigma_E move together); in regimes 21 and 23 (lambda 0.7 for both traits,
residual correlation 0, n = 300) it is the whole phylogenetic versus residual split.

### Table 4. The 32 fits at 1x (campaign), 2x and 4x

| regime | rep | campaign max R-hat | campaign min ESS | 2x max R-hat | 2x min ESS (param) | 4x max R-hat | 4x min ESS (param) | phylolm slope 1x / 2x / 4x | gls slope 1x / 2x / 4x |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 59 | 1.024 | 377 | 1.005 | 1055 (Sigma_E[1,2]) | 1.003 | 2053 (Sigma_E[1,2]) | 0.740 / 0.740 / 0.743 | 0.790 / 0.791 / 0.793 |
| 5 | 65 | 1.007 | 371 | 1.006 | 1034 (Sigma_E[1,2]) | 1.002 | 1778 (Sigma_E[1,2]) | 0.709 / 0.705 / 0.701 | 0.726 / 0.725 / 0.717 |
| 5 | 120 | 1.013 | 296 | 1.003 | 769 (Sigma_E[1,2]) | 1.001 | 1566 (Sigma_E[1,2]) | 0.699 / 0.698 / 0.694 | 0.725 / 0.724 / 0.721 |
| 5 | 182 | 1.018 | 293 | 1.007 | 663 (Sigma_E[1,2]) | 1.004 | 1041 (Sigma_E[1,2]) | 0.690 / 0.689 / 0.699 | 0.706 / 0.708 / 0.718 |
| 5 | 183 | 1.028 | 227 | 1.015 | 738 (Sigma_E[1,2]) | 1.002 | 1551 (Sigma_E[1,2]) | 0.782 / 0.785 / 0.781 | 0.798 / 0.799 / 0.796 |
| 21 | 15 | 1.012 | 339 | 1.003 | 718 (Sigma_P[1,2]) | 1.003 | 1342 (Sigma_P[1,2]) | 0.097 / 0.137 / 0.110 | -0.043 / -0.015 / -0.040 |
| 21 | 20 | 1.016 | 227 | 1.009 | 345 (lambda[1]) | 1.005 | 682 (lambda[1]) | 0.342 / 0.343 / 0.360 | 0.155 / 0.134 / 0.162 |
| 21 | 28 | 1.012 | 336 | 1.004 | 671 (lambda[2]) | 1.002 | 1353 (lambda[2]) | 0.217 / 0.205 / 0.183 | 0.017 / -0.007 / -0.023 |
| 21 | 68 | 1.020 | 227 | 1.012 | 375 (lambda[1]) | 1.007 | 472 (lambda[1]) | 0.315 / 0.315 / 0.311 | 0.212 / 0.225 / 0.162 |
| 21 | 73 | 1.015 | 359 | 1.009 | 630 (Sigma_E[1,2]) | 1.004 | 1307 (Sigma_P[1,2]) | 0.350 / 0.364 / 0.365 | 0.202 / 0.222 / 0.245 |
| 21 | 146 | 1.020 | 318 | 1.007 | 445 (Sigma_E[2,2]) | 1.003 | 1313 (lambda[2]) | 0.201 / 0.225 / 0.216 | 0.136 / 0.151 / 0.156 |
| 21 | 157 | 1.012 | 296 | 1.012 | 509 (Sigma_P[1,2]) | 1.011 | 499 (Sigma_P[1,2]) | 0.229 / 0.213 / 0.224 | 0.103 / 0.083 / 0.089 |
| 21 | 175 | 1.012 | 270 | 1.002 | 1002 (Sigma_P[1,2]) | 1.002 | 1824 (Sigma_P[1,2]) | 0.396 / 0.411 / 0.395 | 0.167 / 0.212 / 0.183 |
| 21 | 182 | 1.016 | 356 | 1.008 | 657 (Sigma_P[1,2]) | 1.001 | 836 (Sigma_P[1,2]) | 0.412 / 0.400 / 0.412 | 0.265 / 0.266 / 0.258 |
| 21 | 188 | 1.013 | 372 | 1.008 | 757 (Sigma_P[1,2]) | 1.004 | 1157 (Sigma_P[1,2]) | 0.278 / 0.284 / 0.293 | 0.042 / 0.022 / 0.016 |
| 21 | 190 | 1.005 | 305 | 1.004 | 778 (lambda[2]) | 1.003 | 1495 (lambda[2]) | 0.227 / 0.227 / 0.238 | 0.105 / 0.113 / 0.128 |
| 23 | 5 | 1.007 | 372 | 1.005 | 636 (Sigma_P[1,2]) | 1.010 | 1257 (Sigma_P[1,2]) | 0.312 / 0.309 / 0.290 | 0.087 / 0.113 / 0.093 |
| 23 | 6 | 1.012 | 229 | 1.010 | 482 (Sigma_P[2,2]) | 1.006 | 1099 (Sigma_P[1,2]) | 0.289 / 0.285 / 0.290 | 0.164 / 0.146 / 0.147 |
| 23 | 37 | 1.011 | 289 | 1.002 | 550 (lambda[1]) | 1.001 | 1056 (lambda[1]) | 0.273 / 0.278 / 0.264 | -0.015 / -0.029 / -0.054 |
| 23 | 59 | 1.017 | 320 | 1.008 | 647 (Sigma_P[1,2]) | 1.004 | 1346 (lambda[1]) | 0.311 / 0.300 / 0.321 | 0.226 / 0.188 / 0.219 |
| 23 | 76 | 1.015 | 261 | 1.002 | 487 (Sigma_P[1,2]) | 1.003 | 886 (Sigma_E[1,2]) | 0.323 / 0.354 / 0.360 | 0.080 / 0.130 / 0.115 |
| 23 | 79 | 1.009 | 387 | 1.004 | 955 (Sigma_P[1,2]) | 1.000 | 2283 (Sigma_P[1,2]) | 0.123 / 0.117 / 0.117 | -0.060 / -0.088 / -0.075 |
| 23 | 87 | 1.021 | 301 | 1.003 | 475 (Sigma_P[1,2]) | 1.004 | 1285 (Sigma_P[1,2]) | 0.280 / 0.284 / 0.294 | 0.195 / 0.181 / 0.188 |
| 23 | 88 | 1.006 | 378 | 1.003 | 519 (Sigma_E[1,2]) | 1.010 | 1058 (lambda[1]) | 0.273 / 0.291 / 0.281 | 0.061 / 0.060 / 0.071 |
| 23 | 111 | 1.006 | 395 | 1.008 | 523 (lambda[2]) | 1.003 | 1513 (lambda[2]) | 0.172 / 0.186 / 0.184 | 0.108 / 0.129 / 0.120 |
| 23 | 144 | 1.016 | 287 | 1.004 | 653 (Sigma_P[1,2]) | 1.005 | 1338 (Sigma_P[1,2]) | 0.129 / 0.109 / 0.103 | 0.030 / 0.009 / 0.001 |
| 23 | 153 | 1.017 | 340 | 1.004 | 865 (lambda[2]) | 1.003 | 1509 (lambda[2]) | 0.308 / 0.322 / 0.317 | 0.262 / 0.291 / 0.275 |
| 23 | 170 | 1.011 | 232 | 1.006 | 518 (Sigma_P[1,2]) | 1.003 | 1267 (lambda[1]) | 0.343 / 0.335 / 0.323 | 0.205 / 0.211 / 0.203 |
| 23 | 174 | 1.010 | 226 | 1.006 | 669 (Sigma_P[1,2]) | 1.003 | 1337 (Sigma_P[1,2]) | 0.335 / 0.334 / 0.342 | 0.174 / 0.161 / 0.170 |
| 23 | 186 | 1.024 | 209 | 1.007 | 555 (Sigma_E[1,2]) | 1.002 | 1225 (Sigma_E[1,2]) | 0.384 / 0.371 / 0.397 | 0.193 / 0.168 / 0.188 |
| 23 | 196 | 1.005 | 190 | 1.003 | 530 (Sigma_P[1,2]) | 1.004 | 915 (Sigma_P[1,2]) | 0.297 / 0.291 / 0.325 | 0.176 / 0.177 / 0.244 |
| 23 | 198 | 1.022 | 284 | 1.013 | 684 (Sigma_P[1,2]) | 1.004 | 1310 (Sigma_E[1,2]) | 0.445 / 0.451 / 0.440 | 0.231 / 0.263 / 0.219 |

Summary (measured):

- 2x: 30 of 32 converged; max R-hat 1.002 to 1.015; min ESS 345 to 1,055. Not converged: regime 21
  reps 20 and 68 (lambda[1], 345 and 375).
- 4x: 32 of 32 converged; max R-hat 1.000 to 1.011; min ESS 472 to 2,283.
- Pooled slopes, 32 reps: signed change 2x minus 1x phylolm +0.0024 (0.0024), gls +0.0012 (0.0039);
  4x minus 1x phylolm +0.0029 (0.0026), gls -0.0007 (0.0042). Mean absolute change: 2x minus 1x
  0.0097 and 0.0174; 4x minus 1x 0.0119 and 0.0175; 4x minus 2x 0.0109 and 0.0172. Pooled SE in
  these reps (4x): phylolm median 0.080 (0.043 to 0.126), gls median 0.105.
- Wall time per fit, 4 chains in parallel on Totoro: median 210 s at 2x, 432 s at 4x (the 1x serial
  re-runs: median 276 s).

Inference: the mean absolute slope change between the 1x and the longer runs is no larger than the
change between the 2x and 4x runs, so it is the Monte Carlo noise of m = 20 draws, not a
convergence effect; there is no systematic shift (signed changes within 1.1 MCSE of 0). The failures
are a shortfall in effective sample size in weakly identified directions, not a failure to find the
posterior (no R-hat above 1.03). The bias in these regimes is unaffected: in the campaign, the
phylolm paired bias of converged and non-converged fits is -0.0121 vs -0.0013 in regime 5, -0.0146 vs
-0.0027 in regime 21 and -0.0104 vs -0.0127 in regime 23 (`campaign_cells.csv`).

### Recommendation for failure 2

- If G6 rule 6 (non-convergence <= 2% at default settings) is to hold: double the default chain
  (burn-in 2,000, 10,000 kept-phase sweeps). Expected effect, inferred from the 32 re-runs and
  assuming fits that converged at 1x still converge at 2x: regime 21 at 2 of 200 (1.0%), regimes 5
  and 23 at 0. Cost: about twice the sampler time for every user (measured n = 300 serial fit about
  250 to 280 s at 1x on Totoro). This changes a default, so the campaign must be re-run.
- Alternative without a default change: an automatic extension when the rule fails (continue the
  same chains until ESS > 400 or a cap). This is an algorithm change and also needs a re-run.
- Alternative without a package change: keep the defaults and the existing warning, and make rule 6
  a harness decision (for example, pre-register "re-run a non-converged fit at 2x and report both").
  A gate change, which is Shinichi's call.

The slopes show that none of these changes the scientific result; the choice is about guarantees
and user wall time.

## What this does not cover

- The lambda = 0.5 gls sign (positive in regimes 4, 6, 8, 10, 12, 14, 16) was not examined; only
  regime 2 was run, where gls is 0.
- Regimes 7, 15 and 4 were not run in the oracle arms (budget); regimes 11 (24 reps) and 3 (40 reps)
  are partial; the twin ran only on regime 3; the prior arm only on regimes 1 and 3 (22 of 24 fits).
- The remaining "campaign minus oracle_kl" gap in regime 9 (about -0.02) is unresolved.
- The flat-mu versus fixed-mu difference between the sampler and oracle_kl was not isolated.
- Regime 21's own phylolm bias (-0.014, MCSE 0.005; in-model DGP) was not investigated; it is not
  driven by the non-converged fits (converged fits -0.0146).
- The prior arm's numerical failure was confirmed on one of the two failed fits only.

## Files

- `script/mi_gls/diag/00_collect_campaign.R`: one row per campaign (regime, rep) from the frozen
  outputs -> `evidence/diagnosis/campaign_cells.csv`.
- `script/mi_gls/diag/oracle_cell.R`, `rerun_cell.R`: the arms above. `make_jobs.sh`,
  `run_queue.sh`: Totoro job lists and queue runner. The queues actually run are in
  `evidence/diagnosis/logs/` (`jobs_*.txt`; the three large lists are kept as their first lines,
  `*_head.txt`, and are the `grep` filters of `jobs_serial.txt` by regime described there).
- `script/mi_gls/diag/collect_diag.R` -> `evidence/diagnosis/oracle_long.csv`, `rerun_summary.csv`,
  `rerun_diag_long.csv`.
- `script/mi_gls/diag/analyse_diag.R` -> `evidence/diagnosis/diag_tables.md` (every table here).
- `evidence/diagnosis/PRIOR_CHANGE.diff`: the one-line change of the prior arm.
- `evidence/diagnosis/logs/`: queue logs and sample job logs, including the two failed prior-arm fits
  and the serial re-run that shows the Cholesky failure.
