# Rose review: `discrete_lambda` and the new gate / floor defaults (2026-10-06)

Scope: `git diff origin/main...feat/discrete-lambda` (2051126, e51c32c, 4b1e37c). I re-ran `test-discrete-lambda.R` and `test-lambda-dispatch.R` (all pass, 1 CRAN skip) and three small probes (below). I did not re-run the full suite.

**Verdict: PROCEED AFTER FIXES** (two blocking, both shipped text; the code is sound).

## What checks out

- Threading: `impute()` and `fit_pigauto()` (both baseline calls, both `model_config` sites) pass it to `fit_baseline()`, then `.fit_baseline_dispatch` / `_route` / `_auto` / `_core`, then threshold-joint (+ `_em`) and OVR (+ `_em`). `multi_impute()` passes it through `...`. The `multi_impute_trees()` shared-GNN replay falls back to `"fixed_1"` for old fits. I found no call site that drops it.
- `"fixed_1"` gives an empty `discrete_idx`, so the old code path runs unchanged. The ab02e31 fixture still matches at 1e-12 (per_column route).
- `options(pigauto.discrete_lambda)` is gone from `R/`, tests and vignettes. The parked "auto" and cumulative-ordinal commits are not on this branch.
- `predict()` never rebuilds the baseline; it uses the stored one. So the OVR lambdas that are not stored cannot make a prediction differ from the fit.
- Tests: no assertion was weakened. Every pinned test now passes `"fixed_1"` (or gate/floor TRUE) explicitly, and new tests cover the new defaults.

## Findings

1. **BLOCKING. Gain attributed to the wrong cause.** NEWS.md:25, R/impute.R:218, R/fit_pigauto.R:232 and R/fit_baseline.R:107 say "estimating [discrete lambda] raised discrete accuracy by 0.04 to 0.12". That range compares gate + floor off with lambda estimated against the old default. Discrete lambda on its own (gate and floor off in both arms) gives +0.02 to +0.04: screen 2 at λ = 0.3 gives 0.523 − 0.502 and 0.542 − 0.511; screen 6 at λ = 0 and 0.1 gives +0.032 to +0.044. *Fix:* "Together with gate and floor off, the new defaults raised discrete accuracy by 0.04 to 0.12 at λ ≤ 0.3 (simulated BM data, n = 100 to 1000); discrete lambda alone contributed +0.02 to +0.04."

2. **BLOCKING. Costs not fully disclosed.** NEWS.md:26-28 and the three roxygen blocks give only "a small cost in probability calibration". *Fix:* state the measured numbers. At λ = 1 with n ≥ 1000: accuracy −0.002 to −0.005 and Brier +0.010 to +0.016. On AVONET (real data, the only real discrete test): categorical Brier +0.014 to +0.020, accuracy unchanged or better. Ordinal still trails BACE at λ = 0.3 (0.356 / 0.379 vs 0.410 / 0.443). Also give the regime of "matched or beat BACE": BM, MCAR 30%, n = 100 and 300, accuracy.

3. NON-BLOCKING. **The zi_count gate is also affected.** Its liability is typed `"binary"` (joint_threshold_baseline.R:227), so it gets an estimated lambda too. In my probe `z_gate` was 0.005 under "estimate" and 0.01 under "fixed_1". This is undocumented and has no evidence. Either document it or exclude it.

4. NON-BLOCKING. **`"fixed_1"` does not mean λ = 1 on the exact route.** There, binary and ordinal columns share `lambda_block` (the old behaviour). In my probe the default "auto" route chose exact, and `bin` reported 0.01. Reword impute.R:216 and fit_pigauto.R:230 as "the previous behaviour (λ = 1, or `lambda_block` on the exact route)".

5. NON-BLOCKING. **Categorical lambdas are not reported.** Under "estimate", `$lambda_per_trait` shows 1 for every `cat3=*` column, yet the predictions do change (max |Δ| 0.79 at missing cells). NEWS.md:12 and fit_baseline.R:112 claim the fitted discrete lambdas appear there. A `lambda_fixed` replay reproduces the fit only because OVR re-estimates on the same data. *Fix:* report NA for categorical columns, or store the per-class values, and reword.

6. NON-BLOCKING. **Partial `lambda_fixed` overrides the default.** With a partial `lambda_fixed`, binary and ordinal columns missing from it silently go to 1 even under "estimate", while categorical is still estimated. Document this.

7. NON-BLOCKING. **No effect without a continuous-family trait.** `discrete_lambda` does nothing for binary or ordinal traits unless at least one continuous-family trait fires the threshold-joint baseline (fit_baseline.R:672). Without one they stay on label propagation. The docs say it applies unconditionally.

8. NON-BLOCKING. **test-exact-default.R:142 pin.** I reproduced the default-regime result at 10 seeds: exact 0.396 vs per_column 0.422. That gap is −0.025, which would fail the test's own 0.02 tolerance. At 50 seeds the gap is −0.014 (SE 0.010) under "estimate" versus 0.000 (SE 0.014) under "fixed_1". So the gap is not significant, but its sign is the same at both seed counts and it never appeared before. "Noise-level" overstates what is known. Keep the pin and update the comment with the 50-seed figures. Open a follow-up: under the new default, exact may lose discrete accuracy to per_column on no-signal data (and "auto" decides between them per trait).

9. NON-BLOCKING. **Misleading test title.** The first test in test-discrete-lambda.R says "byte-identical to the old path", but it only compares `"fixed_1"` with itself. The real old-path pin is the ab02e31 fixture, which covers the per_column route only. Rename the test, or add an exact-route reference.

10. NON-BLOCKING. **Stale text.** Code comments at fit_baseline.R:656-657 and 966-968 still say discrete/OVR columns are "always at lambda = 1". The dev-log README:4 still describes the opt-in option, and its "Does NOT cover" section (README:161) lists n = 1000, real data and the LP path, which later screens covered. The vignettes cite `docs/dev-log/...`, which is excluded from the built package by .Rbuildignore; use a GitHub URL.

11. NON-BLOCKING. **Things users will notice.** `fit$phylo_signal_per_trait` is now all NA by default; say so in NEWS. Multi-obs data is untested; say so. The floor still lowered Brier where discrete lambda costs it (alldisc λ = 1: 0.167 vs 0.173; screen 2, n = 300, λ = 1: 0.165 vs 0.170). "No longer improved accuracy" is true but incomplete; add "calibration slightly worse at λ near 1". The BIEN result is correctly a dev-log measurement (screen 8), not a test.
