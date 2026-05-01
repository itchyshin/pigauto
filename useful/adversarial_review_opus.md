# Adversarial review of pigauto — Opus

> 2026-04-28. Opus 4.7 1M reviewing the prior Sonnet 4.5 session's
> overnight Phase-1 work and the wider codebase. Source files cited
> are absolute paths within the working tree at `experiment/gnn-earnings-sim`.

## Bottom line

The prior session's high-level conclusion ("pigauto's regime of advantage
is real, but is the safety-gated pipeline rather than the transformer
architecture") is broadly defensible — but the **statistical confidence
attached to specific numbers is overstated** (3 reps per cell with CIs
that cross 1.0 at the headline crossover thresholds), and the prior
session's bug-hunting **missed two more bugs of similar pattern** to the
two it fixed. Most importantly, I confirmed by direct simulation a
silent **HIGH severity bug in `multi_impute(draws_method = "conformal")`**:
the default conformal sampler produces **identical datasets across all M
draws** (between-imputation variance literally zero) when the input
data.frame is not in tree-tip order. This is the same row-alignment bug
the prior session fixed in `build_completed`, but in a sister code path
that escapes the existing regression test because the test only checks
correlation with truth, not variance across draws. Multi-imputation
results that fed any Rubin's-rules pooling on shuffled inputs are
unreliable.

## A. Empirical claims

### Lambda threshold (0.15 / 0.30): PARTIALLY agree

The verdict's headline numbers reproduce within rounding from
`script/bench_lambda_sweep.rds`. The monotone trend (lower λ → higher
pigauto/lm_NL ratio) is real and clean. But the **claimed "first wins"
threshold is not statistically distinguishable from a tie**:

| cell | mean ratio | 95% CI (paired) |
|---|---|---|
| λ=0.15 nonlinear (claimed crossover) | 0.914 | **[0.802, 1.025]** straddles 1 |
| λ=0.30 interactive (claimed crossover) | 0.934 | **[0.852, 1.015]** straddles 1 |
| λ=0.20 nonlinear | 0.838 | [0.798, 0.879] clear win |
| λ=0.30 nonlinear | 0.825 | [0.800, 0.851] clear win |

The verdict's "Pigauto first wins at λ=0.15" claim should be re-stated
as "Pigauto first wins definitively at λ ≈ 0.20 (nonlinear), the
interaction interpretation is shakier (definitive only at λ ≈ 0.50)."
This isn't a fatal flaw in the story, but the precise crossover number
is being treated as known when it isn't.

The DGP-fairness concern is real. The interactive DGP is *exactly*
`Z[,1]*Z[,2] + 0.5*Z[,3]^2 + ...` — the saturated polynomial OLS has
zero misspecification. Calling pigauto a winner because it eventually
beats a perfectly-specified competitor is asking pigauto to do
distinctive work the competitor isn't even trying. The framing
"pigauto needs more phylo signal to overcome a saturated correct OLS
than a misspecified one" is honest, but it does highlight that the
crossover threshold is a property of *how mis-specified* lm_NL is, not
a property of pigauto's regime.

The prior session did NOT test stronger off-the-shelf nonlinear
methods (random forest, xgboost, BART). With ncov=10 polynomial features
exploding to 65, these would likely outperform `lm_nonlinear` at low λ
and possibly squeeze pigauto's advantage further. This is a real gap in
the comparison set.

### n-scaling: DISAGREE with the headline interpretation

The verdict claims pigauto's advantage grows with n because "more data
→ better learned function" (the AE story). The actual data tells a
different story:

```
nonlinear DGP, absolute RMSE:
            n=200   n=500   n=1000  n=2000
phylolm     1.19    1.22    1.22    1.20    (FLAT)
pigauto     1.09    1.09    1.10    1.08    (FLAT)
lm_NL       1.17    1.13    1.30    1.25    (DEGRADES with n)
```

Pigauto's RMSE is **flat** across n. Phylolm is **flat** across n.
**lm_nonlinear's RMSE actually gets worse** with n. The ratio "improves"
because the denominator gets worse, not because pigauto does better.
The "AE needs more data" interpretation is not what the data shows.
A more accurate interpretation: at large n, lm_nonlinear's polynomial
extrapolation pathology compounds (more rare-tail covariate values
appear, polynomial predictions blow up at extremes), while
phylogenetic methods are robust to this.

The interactive DGP shows a similar pattern: lm_NL goes from 1.04 to
1.23 as n grows; pigauto improves only marginally (1.22→1.18). The
"pigauto wins at n=2000 even on the interactive DGP" claim is
technically true but it's lm_NL's failure, not pigauto's success.

This is a meaningful re-interpretation. The paper figure should *not*
caption this as "AE methods extract more structure with more data" —
it should be honest that the relative advantage comes from
lm_nonlinear's polynomial extrapolation pathology at large n.

### Beta sweep: AGREE on direction, doubt on precision

The pattern (pigauto's advantage decreases with β) is real, monotone
in 4/5 cells, and consistent with an honest interpretation: at large β,
the saturated polynomial OLS captures most of the signal cleanly and
phylo regularization becomes a drag. The crossover at β=1.0 for
interactive (1.092 vs 1.0) has CI [0.99, 1.21] so the "lm_NL wins by
9%" is not statistically distinguishable from a tie.

The β=0.7 interactive cell has CI [0.62, 1.08] — extremely wide
because reps disagree. This and the lambda-15-nonlinear cell are the
biggest noise points; both are at the edges of the regime statement.

### Missingness: AGREE with caveats

This is the cleanest of the four sweeps. Pigauto/phylolm-λ ratio is
near-deterministic at ~0.91 across all 18 cells with very tight CIs
(SE on the order of 0.001-0.01). The "9/9 nonlinear wins, 9/9
interactive wins, vs phylolm" is robust.

Pigauto/lm_NL on nonlinear is also generally consistent (8/9 wins),
though the absolute values have wider CIs. The headline "pigauto's
advantage GROWS with missingness on nonlinear" is correct in
direction but the trend is partly driven by lm_NL's degradation as
training data shrinks (when total_miss=0.85 there are very few obs
left to fit poly+pairwise on).

Overall verdict: the missingness story is the strongest empirical
result, particularly the "always beats phylolm-λ" claim, which is the
cleanest characterisation in the entire campaign.

## B. Architecture analysis: DISAGREE with framing

The verdict says "transformer is incidental, drop the claim." But:
- Transformer wins **18/18 fits** (paired sign test p ≈ 7.6e-6),
  not "6/6 cells" (binomial p ≈ 0.03 the verdict suggests).
- The per-rep ratio CIs are extremely tight: [0.985, 0.999] for
  the worst cell; [0.989, 0.999] for the median cell. The
  transformer reliably outperforms at the 0.5-1.1% level.
- This **is not noise** in any statistical sense. The verdict
  conflates "small effect size" with "within MC noise." The signal
  is real and trivially statistically significant; it's just
  *practically* small relative to the 8-25% pipeline-level lift.

The defensible framing is: "transformer reliably beats GAT/GCN by a
small (0.5-1.1%) margin in cells where pigauto's pipeline already wins.
The architecture is not the source of pigauto's headline advantage —
the analytical pipeline contributes 95-99% of the gain — but the
direction of the architecture comparison is unambiguous."

Decision (keep transformer as default) is correct. But the rhetorical
framing "architecture is incidental" is too strong.

## C. Code correctness — bugs found

### C.1 [HIGH] `multi_impute(draws_method = "conformal")` between-imputation variance is zero on shuffled input

**File**: `R/multi_impute.R:302-371` (`.sample_conformal_draw`)
**Severity**: HIGH (default code path, silently corrupts MI output, undermines downstream Rubin's rules)
**Evidence**: I ran a direct verification script:
- 15-tip multi-obs tree, shuffled species column, 15 masked cells
- `multi_impute(m=5, draws_method="conformal")` (the default)
- Variance across the m=5 draws, **per masked cell: ALL ZEROS** (15/15)
- Single-obs shuffled: 4/10 cells have zero variance, others have outlier draws

**Root cause** (same pattern as the row-alignment bug fixed in `a3e6d39`):
- `pred$imputed` rows are in **internal (tree-tip) order**
  (`row_labels = object$obs_species`, set by predict.pigauto_fit:317)
- `imputed_mask` rows are in **user-input order**
  (`build_completed` writes `dimnames = list(rownames(original), trait_cols)`,
  `R/impute.R:330-331`)
- `.sample_conformal_draw` does
  `rows <- which(imputed_mask[, nm])` (user-input position) and then
  `imp[[nm]][rows] <- rnorm(...)` (writes to internal-order positions)
- `build_completed(traits, imp_df, species_col, input_row_order=...)` then
  reads `imputed[[nm]][imp_row]` where `imp_row = match(seq_len(n_row), input_row_order)`,
  re-aligning *as if `imp_df` were in internal order*. So the position
  the user sees as "user_2" reads from `imp[input_row_order==2]` — but
  the sampler wrote to `imp[2]` (a different species). The sampled
  values are written to wrong rows and the user never reads them back.

**What downstream breaks**: any user who calls
`with_imputations()` → `pool_mi()` on data with shuffled rows gets
between-imputation variance ≈ 0, FMI = 0 (artificial), total variance
underestimated, Rubin's-rules CIs too narrow. Single-obs partially
affected; multi-obs completely broken.

**Why the test missed it**: `tests/testthat/test-multi-impute.R:507`
checks `cor(pred, truth) > 0.85`. When the point estimate is already
high-quality (the test simulates strong species effects), every
"misrouted" sample lands on essentially the right value (the species
mean), so correlation stays high. The bug only manifests as zero
*sampling* variance, which the test does not check.

**Suggested fix**: in `.sample_conformal_draw`, take `imp_df` as
internal-order, build a parallel `imputed_mask_internal` reordered via
`input_row_order`, and use *that* to find the rows to sample. Or, more
robustly, perform the entire sampling on the internal-order matrix
*before* `build_completed` is called.

### C.2 [MEDIUM] `compute_conformal_scores` uses calibration validation set — selection bias

**File**: `R/fit_helpers.R:428-485`
**Severity**: MEDIUM (affects coverage guarantee, not point predictions)
**Evidence**: The same val_mask is used for `calibrate_gates`
(line 100, in `R/fit_helpers.R`) and `compute_conformal_scores`
(line 452). The gate is selected to minimise validation MSE, so
post-selection residuals on val are biased downward. The conformal
quantile inherits this bias and produces undercoverage.

This is the standard "post-hoc model selection on the same set"
issue. Theoretically: the conformal coverage guarantee assumes
exchangeability of val residuals with test residuals. Selecting `r_cal`
to minimise val loss makes val residuals smaller in expectation than
test residuals, breaking exchangeability.

The user's own coverage probe (memory note "coverage_investigation")
found 0.93 average coverage at n=150 with extreme variance. The
double-dipping is a candidate root cause. Standard fix: split val
into "calibration" half and "conformal" half. The bootstrap conformal
option (Phase 3a) doesn't fix this — it just smooths the quantile
estimator.

### C.3 [MEDIUM] `resolve_one_split` half-A still has NA propagation risk

**File**: `R/fit_helpers.R:244-264`
**Severity**: MEDIUM (similar pattern to the bug fixed in 573decd, but
in a sister branch)
**Evidence**: Line 244 sets `loss_a_pure_bm <- cal_mean_loss(ref_w, ha)`,
line 246 sets `best_la <- loss_a_pure_bm`, line 260 does `if (la < best_la)`.
If `loss_a_pure_bm` is NA (multi_proportion CLR with too few half-A
cells, exact analog of the bisected bug), then `la < best_la` evaluates
to NA and `if (NA)` errors.

The `573decd` fix only patched the half-B half. The half-A logic still
relies on `loss_a_pure_bm` being finite. This may be why the prior
session noted the categorical / binary regressions weren't fully
explained by the NA-fix — the NA path may be triggering on other
DGPs and corrupting the gate via the half-A path before the half-B
guard ever fires.

### C.4 [LOW] predict.pigauto_fit shape mismatch (still unfixed)

**File**: ambiguous; reproduced in `script/phase1_gnn_ablation.R`,
crashes at `R/model_residual_dae.R:249` (`enc1`)
**Severity**: LOW (only affects manual gate override use case, not
production paths)
**Evidence from log**: `linear(): input and weight.T shapes cannot be
multiplied (1473x22 and 23x64)`. The shape arithmetic at training
gives 1+15+7=23 (1 trait + k_eigen=15 + cov_dim=7); at predict the
input has 22 cols. So the predict path is short by 1 column when
modifying the calibrated gates pre-predict.

I could not reproduce on a fresh fit-then-modify-gates flow — the
shape mismatch only fires under a specific condition (likely involving
some interaction between manual override and the predict path's
covariate handling), which I did not have time to fully isolate. The
prior session's note that this is a real bug but only affects manual
gate-zeroing experiments is correct. It blocks the AE-attribution
ablation but not the main paper claims.

## D. Test quality — gaps

1. **The row-alignment bug stayed hidden for years.** The bundled
   `bench_multi_obs.R` simulator uses `rep(tree$tip.label, each = n_obs_per)`
   which puts data in tree-tip order *by accident*. There was no test
   exercising shuffled multi-obs input until 2026-04-26. **Same gap
   exists in tests/test-multi-impute.R:507** — it tests row alignment
   but not draw variance.

2. **The conformal MI variance bug** would be caught by a single-line
   test: assert that for `multi_impute(m=5)` on a missing cell,
   `var(d$trait[mask_idx[1]] across datasets) > 0`. This test does
   not exist.

3. **Coverage / calibration tests are minimal.** `tests/test-safety-floor.R`
   has 23 test_that blocks but they are mostly about
   simplex_grid mechanics and gate calibration math at the API level.
   None of them test the conformal coverage guarantee at the level it
   purports to make ("≥95% marginal coverage"). The user's separate
   `dev/coverage_sim/` is the real probe but is not in the test suite.

4. **No property-based tests for "predictions stable under input
   permutation."** The natural invariant is: for any permutation π
   of the input rows, `impute(traits[π,])` followed by inverse-π
   produces an identical `completed` to `impute(traits)`. Such a test
   would have caught the row-alignment bug (and the conformal
   sampling bug) immediately.

5. **`test-ovr-categorical.R:94` was failing on master before this
   work** and remained skipped/red. The prior session noted but did
   not fix it. I did not run the full test suite to verify whether
   this is genuinely broken or test fragility, but a pre-existing red
   test in a release-bound branch is poor hygiene. It deserves a
   triage at minimum.

6. **No regression bench or test catches discrete-trait performance
   drops < 10pp.** The binary signal_0.6 -7pp regression was found
   only by re-running the full per-type bench on a different machine.
   A "smoke" version of each per-type bench (1 cell, 1 rep, ~30s)
   could be a fast-CI gate. Not running because each per-type bench
   takes 20+ minutes is the proximate reason; there's no good
   alternative architecture for catching slow-drift performance
   regressions.

## E. Apr 17 → Apr 27 untraced regressions — top hypotheses

I cannot do a full bisect in this review (per the time budget), but
based on the file diffs in the window:

### binary signal_0.6 -7pp:
1. **Most likely: `5d83528` (GNN Fix A)** — this commit forces the
   `obs_refine` MLP to fire even in single-obs mode and re-injects
   user covariates after the GNN message passing. For binary traits
   with high phylo signal (signal_0.6 = strong species effect), the
   covariate path adds noise to a baseline that was already very
   close to optimal. The Fix A change is plausibly the proximate
   cause of binary regression.
2. **Second: `faf29e51` (safety_floor)** — the safety floor adds
   the mean-baseline corner to the simplex grid. For high-phylo-signal
   binary traits, the BM-LP corner should always win cleanly, but
   the new simplex search may pick a pure-mean fallback in some
   reps where a few half-B cells happen to favour the constant.
3. **Third: `b4d103d` (median_splits gate calibration)** — changes
   the gate selection from a single split to a median over 31 splits.
   This is meant to reduce bimodal flip-flop, but at small validation
   sets (typical for binary discrete) the median may pull toward 0 in
   cases where a single-split would have picked the right gate.

### categorical K_3 -8pp:
1. **Most likely: `4719ee2` (build_liability_matrix sd_prior_vec arg)** —
   threshold-joint baseline now uses a per-cell prior variance vector.
   For K=3 categorical OVR, each fit operates on a 1-vs-rest split
   where the "1" class has very few examples. A per-cell prior may
   destabilise the EM in low-data regimes that K=2 categorical
   doesn't see.
2. **Second: `62dabd7` (fit_ovr_categorical_fits_em + sd_prior_k
   plumbing)** — direct change to OVR categorical with EM iterations.
   If `em_iterations > 0` is now default for categorical, that's a
   regime shift.
3. **Third: `5d83528` (GNN Fix A)** — same hypothesis as binary;
   covariate injection for categorical OVR may be too aggressive.

### zi_count zf_0.2 +15%:
1. **Most likely: `faf29e51` (safety_floor)** — zi_count uses two
   latent columns (gate + magnitude). The mean baseline path may
   produce a degenerate signal for the gate column (mean of 0/1
   indicators is just the marginal frequency). At high zero-fraction
   (zf_0.2 means 80% zeros), the mean path may dominate via the
   simplex search and lock the gate column to constant prediction.
2. **Second: `b4d103d` (median_splits)** — the discrete-aware
   absolute-cell-floor logic interacts subtly with the gate column
   for zi_count. With MEDIAN of B=31 splits, fewer of those splits
   pass the absolute floor, pushing the gate to 0 more aggressively.

In all three cases, the fastest verification would be: revert the
single most-likely commit on a branch, re-run just the regressed
cell of the per-type bench, and check whether the regression
disappears. ~1-2 hours per regression as the prior session estimated.

## F. Safety floor critique

The safety floor mechanism is **statistically sound at the gate
level** (it is just constrained model selection with the convex hull
of `{BM, GNN, mean}` as the candidate space — there is nothing
unsound about restricting to a convex set), but it has **two
genuine concerns**:

1. **It is double-dipping with the conformal score** (covered in C.2
   above). The same val set is used for gate selection AND conformal
   quantile estimation. This breaks the conformal exchangeability
   guarantee. The "≥95% coverage" claim is no longer formally
   guaranteed once the gate is post-hoc selected on val.
2. **It papers over a training-objective mismatch.** Without the
   safety floor, the GNN delta wanders enough to harm single-obs
   accuracy. The fact that the safety floor "rescues" pigauto from
   this means the model is trained in a way that allows it to learn
   bad deltas, and we then post-hoc constrain the calibration to
   undo them. A better-posed model would not need post-hoc rescue.
   This is more an aesthetic concern than a fatal flaw — many
   regularised methods work this way (lasso, early stopping). But
   the "safety floor IS the contribution" framing in the paper
   should not over-claim that this is novel methodology — it's a
   standard convex-hull projection.

The safety floor is NOT a fatal flaw in the design. It is a sensible
backstop. But the paper should be careful not to position it as an
independent methodological contribution worth claiming credit for.

## G. AE story — NOT supported by direct evidence

The "the autoencoder is doing real work" claim rests on three legs:

1. **n-scaling.** AS DETAILED IN A.2: pigauto's RMSE is *flat* with n
   on both DGPs. The relative advantage grows because lm_nonlinear
   degrades. This is *not* AE-style learning — it's the AE staying
   robust where the polynomial competitor breaks.
2. **Missingness.** Pigauto's advantage grows with missingness, which
   IS consistent with AE literature (MIDAS/MIWAE/GAIN). However,
   pigauto's *absolute* RMSE also degrades with missingness — just
   less than lm_NL's. So the AE isn't "extracting more structure from
   sparse data" — it's just less brittle.
3. **Architecture ablation, gate-zeroing.** The prior session was
   blocked on the predict() shape bug, so could not directly measure
   how much the GNN delta contributes vs the analytical baseline.
   The architecture-ablation finding (transformer ≈ GAT ≈ GCN within
   1%) is itself evidence that the GNN delta is contributing very
   little — if it mattered, different architectures would show bigger
   differences.

**The claim "the AE is doing real work" is NOT empirically supported.**
What IS supported: "pigauto's analytical pipeline (BM/GLS baseline +
multi-obs aggregation + safety-gated calibration) is robust where
polynomial-feature OLS is brittle." The paper should make this exact
claim and not claim more. The "AE story" framing is over-reach.

This is the single most important re-interpretation in this review.

## H. Verdict-log audit

I spot-checked numbers from each of the 5 verdict files against the
underlying RDS:

1. **Lambda sweep**: claimed nonlinear λ=0.15 ratio 0.907; computed
   0.914. Close enough for verdict purposes (rounding from rep-level
   ratios vs cell-mean ratios). ✓
2. **n-scaling interactive n=2000**: claimed ratio 0.958; computed
   0.959. ✓
3. **Beta sweep nonlinear β=0.3**: claimed ratio 0.712; computed 0.715. ✓
4. **Missingness nonlinear sp_miss=0.7 within=0.5**: claimed ratio
   0.814; computed 0.831 (mean of 3 reps). The verdict appears to
   have used the single-rep value 0.814 (cell at total_miss=0.85 in
   the table), not the 3-rep mean. Minor inaccuracy but noticeable.
5. **Architecture ablation**: claimed max spread 1.15%; computed
   1.144%. ✓ Claimed 6/6 cells; computed 18/18 fits. The verdict
   collapsed reps into cells silently, which is conservative for
   the cell-level claim but understates the strength of the
   per-fit signal (binomial p ≈ 0.03 per the verdict's logic vs
   p ≈ 7.6e-6 per per-fit logic, 200,000× different).

Numbers are mostly accurate. Minor sleights:
- Verdict 4 used a single-rep snapshot in one cell instead of
  the 3-rep mean.
- Verdict 6 claims 6/6 (cell level) without acknowledging that
  18/18 (rep level) is stronger statistical evidence. The verdict
  *interprets* this as weaker than it is.
- All verdicts treat 3 reps as if it gave reliable cell-level CIs
  but at the headline crossover thresholds, the CIs cross 1.0
  (lambda nonlinear 0.15: [0.80, 1.03], etc).

## Disagreements with prior session

1. **"AE story is real" — DISAGREE.** The n-scaling pattern is driven
   by lm_NL's degradation, not pigauto's improvement. The architecture
   ablation suggests the GNN delta contributes little. The direct
   ablation that would prove the AE is doing work (gate-zeroing) was
   blocked by a bug and is unfinished. The "AE" framing is unsupported.

2. **"Architecture is incidental" — PARTIALLY DISAGREE.** The
   transformer's 18/18-fit win at 0.5-1.1% margin is a clean
   statistically significant signal, just a *small* one. The
   verdict's framing as "within MC noise" is misleading (the signal
   is reliably above noise; it's just modest). The right framing
   is "small but reliable architecture lift on top of a much larger
   pipeline lift."

3. **"Crossover thresholds at 0.15 / 0.30" — PARTIALLY DISAGREE.**
   With CIs straddling 1.0 at exactly those values, "first wins"
   is being claimed where there's only a tie + monotone trend.
   The defensible threshold (3 reps) is closer to 0.20 nonlinear /
   0.50 interactive.

4. **"NA-handling bug was the multi_proportion crash root cause" —
   AGREE for the crash, DISAGREE for the broader pattern.** The
   half-A path in `resolve_one_split` (line 244-264) has the same
   NA-propagation pattern that was patched in half-B, and may be a
   contributing factor to the binary/categorical/zi_count
   performance regressions that the NA-fix did not resolve.

5. **"Tests cover the row-alignment fix adequately" — DISAGREE.**
   The fix to `build_completed` was validated, but the same row-
   alignment pattern in `.sample_conformal_draw` was missed
   entirely. The test for multi_impute alignment only checks
   correlation, not draw variance, missing the conformal-MI bug.

## Severity ranking — top issues to fix

1. **[HIGH] Conformal MI between-imputation variance is zero on
   shuffled inputs.** `R/multi_impute.R::.sample_conformal_draw`.
   Default code path for `multi_impute()`. Silently corrupts
   downstream Rubin's-rules variance and FMI. Reproducer:
   any shuffled-row multi-obs input. Fix: do the sampling on
   internal-order matrices, then re-align via input_row_order
   AFTER sampling.

2. **[HIGH] `predict.pigauto_fit` shape mismatch** when calibrated
   gates are manually overridden. Blocks any AE-attribution
   ablation experiment, weakens the paper's empirical claim
   ("we never directly measured the GNN's contribution"). The
   prior session never root-caused this. Worth ~30 min to find
   the dimension-count divergence point in the predict path.

3. **[MEDIUM] Conformal scores use calibration validation set
   (double-dipping).** `R/fit_helpers.R::compute_conformal_scores`.
   Breaks the formal coverage guarantee. Standard fix: split val
   into calibration and conformal halves. The user's own coverage
   probe shows symptoms (extreme variance at small n).

4. **[MEDIUM] `resolve_one_split` half-A path has same NA-
   propagation risk** as the half-B path that was patched in
   `573decd`. May be silently contributing to the unresolved
   discrete-trait regressions. ~10 min to add the same `is.finite`
   guard.

5. **[MEDIUM] Paper claims "AE is doing work" without evidence.**
   The n-scaling pattern is lm_NL's degradation, not pigauto's
   improvement. The architecture ablation is consistent with the
   GNN contributing very little. Either reframe the paper to claim
   "robust analytical pipeline" instead of "denoising AE", or
   complete the gate-zeroing ablation (blocked on issue #2 above)
   to provide the missing direct evidence.
