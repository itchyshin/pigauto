# Mechanism-coverage campaign (B1) — results

Campaign complete 2026-08-16 on Totoro. **200 jobs, 0 failures** (~25 min total across
two dispatch waves; throttled mid-run from 40 to 15 workers for D-143 — see the
compute-lane lesson doc). Design + pre-registered decision rule:
`2026-08-16-mechanism-coverage-design.md`. Runner `mech_cell.R`, summariser
`summarise_mech_cov.R`, artifacts `~/pigauto_regime_map/results/mech_cov/`.

Coverage is conformal-interval coverage at nominal 0.95, scored on **genuinely-missing
cells** — the real-usage surface, where calibration cells come from the observed
complement. MCSE = sd(per-rep coverage)/√n_rep.

## The verdict: exchangeability failure CONFIRMED

| mechanism | n | n_rep | coverage | MCSE | verdict vs MCAR |
|---|---:|---:|---:|---:|---|
| MCAR | 300 | 30 | 0.961 | 0.004 | control ok (≥0.94) |
| MAR_trait | 300 | 30 | 0.949 | 0.003 | not distinguishable |
| **MAR_phylo** | 300 | 30 | **0.923** | 0.006 | **BELOW MCAR (>3 MCSE)** |
| **MNAR** | 300 | 30 | **0.939** | 0.006 | **BELOW MCAR (>3 MCSE)** |
| MCAR | 1000 | 20 | 0.957 | 0.002 | control ok |
| **MAR_trait** | 1000 | 20 | **0.940** | 0.004 | **BELOW MCAR (>3 MCSE)** |
| **MAR_phylo** | 1000 | 20 | **0.927** | 0.005 | **BELOW MCAR (>3 MCSE)** |
| **MNAR** | 1000 | 20 | **0.932** | 0.003 | **BELOW MCAR (>3 MCSE)** |

The pre-registered rule (any MAR/MNAR arm below MCAR by >3× pooled MCSE, with MCAR ≥
~0.94) fires at both n. Three secondary observations:

1. **The pipeline itself is healthy**: under MCAR the intervals slightly over-cover
   (0.957–0.961), as split conformal should (the ⌈(1−α)(n+1)⌉ quantile is conservative).
   The defect is specifically distributional, not mechanical.
2. **Phylo-structured missingness hurts most** (−3.4 pp), value-MNAR next (−2.5 pp),
   observed-trait-MAR least (−1.2 pp, only separable at n=1000). Mechanistically:
   calibration residuals come from well-sampled regions of the tree; genuinely-missing
   cells concentrate in poorly-sampled clades where prediction errors are larger, so the
   calibration quantile is too small for the target population.
3. **It does not wash out with n** — coverage at n=1000 is no better than at n=300
   (MAR_phylo 0.927 vs 0.923). Consistent with the real-data record: the arithmetic
   ceiling explains nothing at large n; the mechanism does.

## Relation to the real-data undercoverage

Direction and mechanism match the committed benches (fishbase 0.89–0.91 at n=10,654;
pantheria 0.87–0.94), but the simulated magnitudes (2–3.5 pp) are smaller than the worst
real cases. Real trait databases plausibly combine the mechanisms (clade-structured AND
value-dependent sampling) and add model misfit; this campaign establishes the mechanism
exists at realistic strengths, not that it accounts for every observed gap. The 0.65
pantheria `gestation_d` case (n_val small) is likely ceiling + mechanism stacked.

## Decision: B2 proceeds

Per the pre-registered rule, the weighted/Mondrian conformal fix is warranted. The
MAR_phylo ranking constrains the design: the correction must condition on **phylogenetic
sampling locality** (e.g. distance-to-nearest-calibration-cell or local observed density),
not merely on trait values. Candidate estimators, to be specified before implementation:

- **Mondrian (stratified) conformal**: stratify calibration cells by a locality statistic
  (e.g. quantile bins of mean phylogenetic distance to observed conspecific cells);
  per-stratum quantiles. Guarantees per-stratum validity if strata are exchangeable
  within themselves; costs width where strata are small.
- **Weighted conformal (Tibshirani et al. 2019)**: likelihood-ratio weights from an
  estimated missingness propensity (logistic on phylo eigenvectors + observed traits).
  More elegant, needs a propensity model good enough to trust.

Regime: F1 continuous, λ=1, 30% missing, single-obs, simulated; mechanisms as
parameterised in the design doc. Discrete types not measured (no conformal intervals).

## Verification target for B2 (pre-registered)

On this same campaign grid, the fixed method must bring every MAR/MNAR arm to within
3×MCSE of 0.95 **without** dropping MCAR below 0.94 or inflating median width by more
than ~10% under MCAR. Re-run via `mech_cell.R` with the new option; same seeds.
