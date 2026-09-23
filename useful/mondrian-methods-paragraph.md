# Methods paragraph: locality-stratified (Mondrian) conformal intervals

Drafted 2026-09-22 for the pigauto paper. Every number below is read from `R/fit_helpers.R`
(`compute_conformal_scores`, `mondrian_locality`), `NEWS.md`, and
`docs/dev-log/2026-08-16-mechanism-coverage-results.md` on `origin/main`. Prose in Shinichi's
register: plain, precise first, no em dashes.

## Suggested placement

A new subsection under the uncertainty-quantification part of the Methods, after the paragraph
describing the split conformal interval and before the multiple-imputation draws.

---

### Prediction intervals

For every imputed continuous, count, ordinal or proportion trait we report a 95% prediction
interval built by split conformal prediction (Papadopoulos et al. 2002; Vovk et al. 2005). Before
training, a random 25% of the observed cells is withheld and further divided into a validation set
(one quarter) and a test set. The model is fitted without them. On the validation cells the
absolute residual between the true value and the blended prediction is computed on the latent
(standardised) scale, and for each trait the interval half-width is taken as the empirical quantile
of these residuals at level ⌈(1 − α)(n + 1)⌉ / n, with n the number of validation residuals and
α = 0.05. The interval for a missing cell is the prediction plus or minus this half-width,
back-transformed to the trait's original scale. The construction is distribution-free: if the
residual of a new cell is exchangeable with the validation residuals, the interval contains the
true value with probability at least 1 − α (Lei et al. 2018).

That exchangeability condition is the one assumption the method makes, and phylogenetic data
violate it in a specific way. Validation cells are drawn from the observed part of the matrix,
which sits in well-sampled regions of the tree where the baseline can lean on close relatives and
prediction errors are small. The cells a user actually needs to impute are not so placed; in real
trait databases missingness concentrates in poorly sampled clades, where errors are larger. A
single quantile calibrated on the first population and applied to the second is too short. In a
simulation with clade-structured missingness (two clades at 7:1 odds of being missing against the
background), empirical coverage of the nominal 95% interval fell to 0.923 at n = 300 species and
0.927 at n = 1000, against 0.961 and 0.957 under completely random missingness; the shortfall did
not diminish with sample size (`docs/dev-log/2026-08-16-mechanism-coverage-results.md`).

To restore the condition where it fails, pigauto offers a locality-stratified variant
(`conformal_method = "mondrian"`; Vovk 2012; Boström et al. 2021). For each validation cell a
locality statistic is computed as the mean cophenetic distance from its species to the five nearest
species with an observed value for that trait. Validation cells are split at the median locality
into a near and a far stratum, and a separate conformal quantile is computed within each stratum at
the same adjusted level. At prediction time a missing cell receives the half-width of the stratum
its own locality places it in, so intervals widen in undersampled clades, which is where the error
is. Within a stratum the near-versus-far mismatch that broke exchangeability is largely removed, and
the coverage guarantee holds per stratum rather than only on average. On the same simulation grid
this recovered clade-structured coverage to 0.946 at n = 1000, within three Monte Carlo standard
errors of nominal, while widening the median interval under random missingness by 2.3%.

The variant has a floor. A stratum's own conservative quantile can only reach 0.95 when it holds at
least 19 residuals (the smallest n for which n / (n + 1) ≥ 0.95), so a trait whose validation cells
would leave either stratum below 19 (roughly 38 validation cells per trait) falls back to the global
split quantile, and the software records that it did. At n = 300 in the simulation above this
fallback fired for every trait and the variant changed nothing; the remedy at that size is more
held-out data rather than stratification. The stratified variant is currently limited to
single-observation data, because its locality is defined per species, and to traits that receive
conformal intervals at all. It has been evaluated on simulated missingness mechanisms; the default
remains the unstratified split quantile until the real-data benchmarks in which the undercoverage
was first observed have been re-run with it.

---

## References to add if not already present

- Papadopoulos H, Proedrou K, Vovk V, Gammerman A (2002) Inductive confidence machines for
  regression. ECML 2002, LNCS 2430, 345–356.
- Vovk V, Gammerman A, Shafer G (2005) Algorithmic Learning in a Random World. Springer.
- Vovk V (2012) Conditional validity of inductive conformal predictors. ACML 2012, PMLR 25,
  475–490. (Mondrian / label- and taxonomy-conditional conformal.)
- Lei J, G'Sell M, Rinaldo A, Tibshirani RJ, Wasserman L (2018) Distribution-free predictive
  inference for regression. JASA 113(523), 1094–1111.
- Boström H, Linusson H, Löfström T, Johansson U (2021) Mondrian conformal regressors. COPA 2021,
  PMLR 152, 24–41.

Citation years and venues are from memory and must be checked against the DOIs before submission;
Garfield's rule applies.
