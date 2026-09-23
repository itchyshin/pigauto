# Mondrian-scaled multi_impute() draws: does the stratum score make Rubin SEs more honest?

2026-09-23.

## What pigauto does today

`multi_impute(draws_method = "conformal")` runs `impute()` once, then for
each of `m` originally-missing cells samples a Normal draw centred on the
point estimate with SD = conformal_score / 1.96, on the appropriate
transformed (latent) scale (`.sample_conformal_draw()`,
`R/multi_impute.R`). Discrete types draw from Bernoulli / Categorical
instead. Until this change, that conformal score was always the single
global split-conformal score for the trait, even when the fit itself was
calibrated with `conformal_method = "mondrian"` and therefore carries
per-trait near/far stratum scores (`fit$conformal_mondrian`).

`predict.pigauto_fit()` already used the per-row near/far stratum score
for its own conformal interval output (B2, 2026-08-16): a target cell's
half-width depends on its phylogenetic locality, via
`mondrian_locality()`, relative to the training-time observed set. This
change extracts that per-row logic into `mondrian_cell_scores()`
(`R/predict_pigauto.R`) and threads it into `.sample_conformal_draw()` so
a Mondrian-calibrated fit's conformal draws use the same per-cell score
that its own conformal intervals already report. Split / bootstrap
calibration is untouched: `mondrian_cell_scores()` returns the constant
global score for those, so `.sample_conformal_draw()`'s behaviour is
identical to before this change.

## The argument

Split-conformal calibrates on a validation set that is, by construction,
close (in cophenetic distance) to the rest of the observed data. Under
clade-structured missingness (`MAR_phylo`), the genuinely-missing cells
are systematically farther from the observed set than the validation
cells were from each other. A single global conformal score therefore
undercovers exactly where it matters most -- the isolated clade -- and
overcovers near well-sampled clades. `conformal_method = "mondrian"`
addresses this for the *point* interval by splitting calibration residuals
into a near stratum and a far stratum and reporting the per-row applicable
one.

The argument for threading this into `multi_impute()`'s draws is that the
same mismatch propagates into the SD used for stochastic completions: if
draw SD tracks the (undercovering) global score, downstream Rubin pooling
underestimates between-imputation variance for isolated cells and,
mechanically, the pooled total SE for any coefficient whose leverage
comes disproportionately from that clade. Using the per-cell Mondrian
score as the draw SD -- far cells get a wider Normal, near cells a
narrower one -- is a plausible fix in the same spirit as B2's fix for the
point interval.

## Limits

Four reasons this is not, on its own, a validated inferential procedure:

1. **A 95% quantile is not a standard deviation.** `score / 1.96` assumes
   the calibration residual distribution is close enough to Gaussian that
   its 95th percentile and its SD are related by the Normal `z_0.975`
   constant. Conformal calibration gives no such guarantee; it is
   distribution-free by construction, which is precisely what makes the
   `/1.96` conversion a heuristic rather than a derived quantity.

2. **Conformal validity is marginal within a stratum, and two strata are
   a coarse step function.** Mondrian conformal guarantees marginal
   coverage *within each stratum*, treated as its own exchangeable pool.
   That is weaker than pointwise validity at every phylogenetic distance,
   and collapsing a continuous quantity (cophenetic distance to the
   nearest `k` observed tips) into two bins (near / far) at a single
   threshold is a coarse approximation to whatever the true relationship
   between locality and residual dispersion is. A cell just inside the far
   stratum and a cell just inside the near stratum can be phylogenetically
   almost equidistant from the observed set yet receive different SDs.

3. **The draws are not proper multiple imputations.** A proper MI draw
   should sample from (an approximation to) the posterior predictive
   distribution implied by re-fitting the *entire* generative model per
   draw. Here the model (BM baseline + GNN, gate, conformal calibration)
   is fit once; only the missing-cell values are perturbed post hoc around
   a fixed point estimate. This is the same limitation `multi_impute()`'s
   own roxygen already flags for the non-Mondrian case; Mondrian scaling
   changes the SD used, not this structural limitation.

4. **The stratum score is fixed noise across all `m` draws.** Because the
   near/far score for a given cell is a deterministic function of the
   fitted model and the tree (not a random quantity redrawn per
   imputation), every one of the `m` draws for a given missing cell uses
   the identical SD. The only randomness across draws is the Normal noise
   itself; there is no uncertainty in the scale parameter propagated
   into the pooled variance, unlike, say, a Bayesian bootstrap over the
   calibration set.

## Pre-registered claim

Under `MAR_phylo` missingness (clade-structured, as in
`~/pigauto_regime_map/mech_cell.R`), for a downstream `gls(y ~ x,
correlation = corBrownian(tree))` slope with `x` subject to missingness
and `y` fully observed:

- Mondrian-scaled draws give a Rubin SE ratio -- mean pooled total SE over
  the empirical SD of the pooled point estimate across reps -- closer to
  1 than split-scaled draws.
- Far-stratum cells' draw coverage of the masked truth (fraction of reps
  where the truth falls within the empirical 2.5-97.5 percentile range of
  that cell's `m` draws) is within 2 percentage points of 0.95.

Margins are fixed now, before the campaign runs: "closer to 1" means the
Mondrian arm's `|SE ratio - 1|` is smaller than the split arm's by at
least 0.03 (a fixed, pre-declared margin, not a post hoc significance
test), and "within 2 pp of 0.95" means far-stratum coverage falls in
`[0.93, 0.97]`.

## Simulation design

DGP reused from `~/pigauto_regime_map/mech_cell.R` (fetched via `ssh -o
BatchMode=yes snakagaw@totoro.biology.ualberta.ca 'cat
~/pigauto_regime_map/mech_cell.R'`), simplified to two correlated BM
traits (`x`, `y`, correlation 0.7) instead of that script's four, and its
`MAR_phylo` mechanism (two randomly-sized clades, each 15-35% of tips,
concentrate missingness at a clade:background odds ratio of 7:1) applied
to `x` only; `y` stays fully observed.

- `n = 1000` tips, `MAR_phylo` missingness on `x`, target `m_miss ~ 0.3`
  overall (matching `mech_cell.R`'s default).
- Downstream model: `nlme::gls(y ~ x, correlation =
  ape::corBrownian(phy = tree), data = dat, method = "ML")`, fit once per
  imputed dataset.
- `m = 20` imputations per method per rep; `500` paired reps (same tree,
  mask, and truth per rep across the split and Mondrian arms, so the
  comparison is a paired one).
- Three arms per rep: complete-data reference (truth `x`, no
  missingness -- establishes the target sampling SD of the slope with no
  imputation uncertainty at all), split-scaled MI, Mondrian-scaled MI.
- Metrics per rep, saved to one `.rds`:
  - pooled slope estimate, total SE, `df`, `fmi` (from `pool_mi()`) for
    each of the split and Mondrian arms;
  - the complete-data arm's slope estimate and model SE;
  - bias = pooled estimate - true slope, for each arm;
  - per-missing-cell draw PIT (`mean(draws <= truth)`) and empirical
    2.5-97.5 percentile coverage indicator, split by near/far stratum
    (stratum assignment from `mondrian_cell_scores()`'s threshold logic
    applied to the Mondrian fit).
- Cross-rep summary (`06b_summarise_mi_se.R`): SE ratio (mean pooled SE /
  empirical SD of the pooled estimate across the 500 reps) with a
  delta-method Monte Carlo SE on that ratio; coverage of the pooled
  slope's own CI against the true slope; bias; mean FMI; per-cell PIT and
  coverage aggregated by stratum.

Estimate before running the full campaign: `mech_cell.R`'s own single-cell
single-method fits at `n = 1000`, `epochs = 500` run in the range of
several minutes each on Totoro (per-node wall time noted in that
script's companion logs); this design fits `impute()` (or `multi_impute()`
internally calling it once) **twice** per rep (split arm, Mondrian arm) at
that same `n` and `epochs`, so a single paired rep is of that same order,
times two. 500 reps at even 2-3 minutes per rep would be several
compute-days if run serially; this is a DRAC-job-array candidate ("Totoro
or DRAC?" — see `~/shinichi-brain/projects/COMPUTE-PLAYBOOK.md`), not a
laptop or single-Totoro-process campaign. No campaign has been run yet;
see script/mondrian_confirmation/06_mi_se_sim.R's own header for the
one-rep-per-invocation contract and script/mondrian_confirmation/
06b_summarise_mi_se.R for the aggregator. A smoke run at `n = 200`,
`epochs = 20` (far below the pre-registered `n = 1000`, `epochs = 500`,
and therefore not informative about the claim above) confirmed the
pipeline runs end to end and produces both methods' output in one `.rds`;
wall time is reported in the after-task summary for this task, not
repeated here since it does not extrapolate linearly to the full design.

## Result

Pending: pre-run then campaign.
