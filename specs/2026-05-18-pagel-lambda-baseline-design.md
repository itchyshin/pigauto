# Pagel's λ for the BM baseline — design spec

**Status:** **APPROVED 2026-05-18.** All 4 open questions resolved (§8) — see decisions below.
**Branch (target):** `feature/pagel-lambda-baseline` (off `main` post-v0.9.4).
**Target version:** v0.10.0 (minor bump; structural baseline change).
**Motivated by:** `useful/sla_diagnosis_2026-05-18.md` (created in the v0.9.4 W-sweep): on the BIEN h2h bench, pigauto's joint MVN baseline beats grand-mean on `sla` by only ~7% (z-scored RMSE 0.91 vs 0.98) and is bit-identical to per-column BM. BACE's MCMCglmm achieves z-scored ≈ 0.69 on the same trait. Root cause: both the per-column BM kernel (`R/bm_internal.R`) and the Rphylopars joint MVN path are fixed at Pagel's λ = 1; neither can shrink toward the grand mean when phylogenetic signal is weak.

---

## 1. Goal

Enable Pagel's λ < 1 in the BM baseline so weak-phylo-signal traits get appropriate shrinkage toward the grand mean. Concretely: lift `sla` (and similar weak-signal continuous traits like the BIEN plant suite) toward BACE-level RMSE without touching the GNN.

Hard target: on the v0.9.4 BIEN h2h bench, `sla` z-scored val RMSE drops from 0.91 to ≤ 0.75 (BACE-level), keeping all other traits at their v0.9.4 numbers or better.

Soft target: ≥ 30% of the per-trait BACE-vs-pigauto gap on `sla` / `height_m` / `wood_density` closes; total BIEN_RMSE_TOTAL drops from 7.236 to ≤ 6.8.

## 2. Background

### 2.1 Where λ=1 lives

Two independent code paths both assume λ=1 silently:

1. **Per-column BM (`R/bm_internal.R::bm_impute_col`)** — the BIM (Brownian Imputation Model) used directly when joint MVN isn't available, and as the gradient kernel inside every downstream baseline. Builds `R = cov2cor(vcv(tree))`, then conditional MVN with that R. No λ parameter.
2. **Joint MVN (`R/joint_mvn_baseline.R::fit_joint_mvn_baseline`)** — calls `Rphylopars::phylopars(model = "BM")` per fit. phylopars *does* support `model = "lambda"`, but pigauto doesn't expose the choice.

Both downstream paths inherit the limitation:

- Threshold-joint (binary / ordinal liability): builds liability matrix → phylopars → λ=1
- OVR categorical: K independent threshold-joint fits → λ=1 each
- GNN's gated residual baseline (`baseline_mu` argument): uses the joint MVN output → λ=1

### 2.2 Why λ=1 is wrong for weak-signal traits

Under BM with λ=1, the conditional posterior at a missing tip is
$\hat{\mu}_i = \beta + R_{iO} R_{OO}^{-1}(y_O - \beta)$
which is the kriging predictor at full phylogenetic weight. When λ < 1,
$R(\lambda) = \lambda R + (1-\lambda) I$
which interpolates between full BM ($\lambda=1$) and uncorrelated diagonal ($\lambda=0$). For traits with weak phylo signal (Pagel's λ posterior near 0), the right move is large shrinkage of $y_O - \beta$ toward zero — i.e. predict the grand mean for missing cells, not the BM kriging value. That's exactly what BACE's MCMCglmm does via its variance-component posterior.

`R/build_phylo_graph.R` already computes per-trait Pagel's λ estimates (used by the `phylo_signal_gate`). So the package knows the right λ values per trait; it just doesn't feed them back into the baseline kernel.

### 2.3 What BACE does that pigauto doesn't

BACE fits MCMCglmm with phylogenetic random effect + residual variance + Bayesian priors. The posterior over the variance ratio (phylo / phylo + residual) is BACE's analogue of (1−λ). On weak-signal traits the posterior pulls toward 50/50 or smaller; the resulting kriging predictions look like Pagel's-λ-BLUP at the posterior λ.

## 3. Non-goals

- **OU, EB, kappa, delta models.** phylopars supports them; pigauto could too. Out of scope. λ is the single-parameter generalization that addresses the demonstrated weak-signal failure mode.
- **Joint estimation across traits.** Each trait gets its own λ (or a single λ for the joint MVN block). We do NOT estimate a multivariate-λ in the sense of separate λ per off-diagonal Σ entry.
- **Multi-obs aggregation changes.** The B1 soft path is unchanged.
- **Binary / categorical liability λ.** Threshold-joint baseline gets λ via the underlying phylopars call (free upgrade), but we don't separately tune the liability prior.
- **API removal.** Existing callers must continue to work. `bm_impute_col()` keeps λ=1 default.

## 4. Architecture

### 4.1 Per-column BM (`R/bm_internal.R`)

Add a `lambda` parameter to `bm_impute_col()`:

```r
bm_impute_col <- function(y, R_phylo, lambda = 1.0, ...)
```

- `lambda = 1.0` → behaviour unchanged (back-compat).
- `lambda` numeric in `[0, 1]` → use $R(\lambda) = \lambda R + (1-\lambda) I$ in the conditional MVN.
- `lambda = "estimate"` → fit ML λ on the observed cells of `y` using the standard Pagel transform. Falls back to the package's existing helper if there's one; otherwise inline ML over a 1-D grid + `optim`.

Per-column ML λ is cheap: closed-form for the profile likelihood at fixed λ + a 1-D optimisation. The full fit is matrix-inverse dominated. Expected runtime per column: ~10–50 ms at n=2000.

### 4.2 Joint MVN (`R/joint_mvn_baseline.R`)

Add `lambda` parameter to `fit_joint_mvn_baseline()` and pass through to phylopars:

```r
fit_joint_mvn_baseline <- function(..., lambda = 1.0)
```

- `lambda = 1.0` → `phylopars(model = "BM")`, unchanged.
- `lambda = "estimate"` → `phylopars(model = "lambda")` (phylopars estimates λ jointly with Σ).
- Numeric λ ∈ [0, 1] → use `phylopars(model = "BM")` after pre-transforming the tree's branch lengths by Pagel's λ. (phylopars does not accept a fixed λ directly; the tree-transform trick is standard.)

### 4.3 Dispatcher (`R/fit_baseline.R`)

Add `lambda_mode` to `fit_baseline()`:

```r
fit_baseline <- function(..., lambda_mode = c("fixed_1", "estimate"))
```

- `"fixed_1"` (default) → preserves all v0.9.x behaviour exactly.
- `"estimate"` → estimate per-trait λ via the joint MVN path (preferred) or per-column ML as fallback.

The dispatcher chooses:
- joint MVN available + ≥ 2 BM-eligible cols → `phylopars(model = "lambda")` once for the full block, gives one λ per trait (phylopars's diagonal Σ structure)
- joint MVN not available OR 1 BM-eligible col → per-column ML λ

### 4.4 Threading through `fit_pigauto()` and `impute()`

Add `lambda_mode = "fixed_1"` to both. Plumbs to `fit_baseline()`. Default keeps current behaviour; users opt in.

Future-friendly default question: should the default be `"estimate"` for new fits? Discussed in §7 (rollout).

### 4.5 Caching / reuse

Computed λ is stored on the `pigauto_fit` so `predict.pigauto_fit()` can rebuild the same baseline at prediction time. Add `model_config$lambda_per_trait` (named numeric vector keyed by trait latent column).

### 4.6 What the GNN sees

The GNN's `baseline_mu` input changes when λ ≠ 1, but it's still the same shape (n × p_latent). No GNN code changes. The gate calibration will adjust automatically to the new baseline quality — if the λ-baseline is better, some gates will close further; if it changes the residual structure, others may open.

## 5. Testing

### 5.1 Unit tests (`tests/testthat/test-pagel-lambda.R`, new file)

1. **λ=1 produces identical results to v0.9.4.** Run a fit with `lambda_mode = "fixed_1"` and verify the predictions match the current `fit_baseline()` exactly. Guards back-compat.
2. **λ=0 produces grand-mean predictions.** For continuous traits with `lambda = 0`, every missing cell's predicted z-score should be 0 (grand mean).
3. **λ=0.5 falls between λ=0 and λ=1.** Sanity check the interpolation.
4. **ML λ matches a known reference.** Simulate a known-λ trait under BM(λ), recover the ML λ to within 0.05.
5. **Joint MVN with λ-estimate matches per-column ML when Σ is diagonal.** Edge case where joint MVN shouldn't add anything.

### 5.2 Smoke test on AVONET 300

A fast trait-imputation roundtrip on the bundled AVONET 300 dataset with `lambda_mode = "estimate"`. RMSE on continuous traits should be within 5% of v0.9.4 (AVONET has strong phylo signal; λ should converge near 1). This is the safety regression check.

### 5.3 BIEN h2h CI benchmark

The acceptance criterion. Pre/post comparison of `script/gha/_bace_compat.R` BIEN dataset:

```
v0.9.4 baseline:
  height_m     1.285  (gate 0.68)
  leaf_area    1.943  (gate 1.00)
  sla          2.041  (gate 0)
  seed_mass    1.728  (gate 0)
  wood_density 0.238  (gate 0.40)
  TOTAL        7.236

Target with lambda_mode = "estimate":
  sla          ≤ 1.60  (closes ≥ 80% of the gap to BACE's 1.55)
  height_m     ≤ 1.28  (no regression; possibly small lift)
  others       ≤ baseline value + 0.05  (no regression)
  TOTAL        ≤ 6.85   (closes ≥ 0.4 of the 0.49 sla gap)
```

If `sla` regresses or the total moves the wrong way, the spec is wrong and we don't ship.

## 6. Risks

### 6.1 ML λ instability at small n

Per-column ML λ can be ill-conditioned when n_obs is small (< 30) for that trait. Mitigation: clip the optimisation to `[0.01, 0.99]` (avoid the singular endpoints), set a sensible prior toward 0.5 if data is sparse.

### 6.2 Existing callers break

`bm_impute_col()` is called from at least 6 internal sites. Adding an optional parameter with `lambda = 1.0` default is fully back-compat, but I need to audit each callsite to confirm none of them pass arguments positionally past the existing 4th arg. Action item: read every callsite before the PR.

### 6.3 Multi-tree workflow interaction

`multi_impute_trees()` re-fits per tree. λ should be re-estimated per tree (the tree changes). Confirm the multi-tree fit pipeline doesn't cache λ from the MCC tree.

### 6.4 Phylopars model = "lambda" might be slow

phylopars's λ-mode adds an ML estimation loop. Worst-case 2-5x slower than `model = "BM"`. Budget: at n=2000 the W0 fit was ~9 min; the λ-mode could go to 15-25 min. Acceptable for batch but expensive for the h2h CI runs across all 6 datasets. Need to time on a sample dataset before committing.

### 6.5 Default change in v0.10

If we change the default from `"fixed_1"` to `"estimate"`, existing user fits get implicitly upgraded. Per the user's standing rule in CLAUDE.md ("Options over forced changes — prefer opt-in flags with backward-compat defaults; never break existing fits silently"), we ship the feature as opt-in in v0.10.0 and revisit the default after a minor release cycle of feedback.

## 7. Rollout

### 7.1 v0.10.0 (this PR)

- `bm_impute_col()` gains `lambda` parameter, default 1.0
- `fit_joint_mvn_baseline()` gains `lambda` parameter, default 1.0
- `fit_baseline()` gains `lambda_mode` parameter, default `"fixed_1"`
- `fit_pigauto()`, `multi_impute()`, `impute()` gain `lambda_mode`, default `"fixed_1"`
- Unit tests + AVONET smoke + BIEN h2h verification
- NEWS entry describing the option

### 7.2 v0.10.1 or v0.11 (later, after feedback)

Consider flipping the default to `"estimate"` if user testing shows the BIEN-style lift generalizes and the AVONET-style strong-signal case stays a no-op (λ → 1).

### 7.3 Files changed (estimate)

- `R/bm_internal.R` — ~50 LOC additions for λ parameter
- `R/joint_mvn_baseline.R` — ~30 LOC for the phylopars dispatch
- `R/fit_baseline.R` — ~20 LOC for the dispatcher
- `R/fit_pigauto.R` — ~5 LOC for the model_config field
- `R/predict_pigauto.R` — ~10 LOC for re-using the stored λ
- `tests/testthat/test-pagel-lambda.R` — ~200 LOC new file
- `NEWS.md`, `DESCRIPTION` — version bump + entry
- `man/*.Rd` — regenerated via document()

Estimate: half-day to full-day of focused work, dominated by audit + testing of the callsites of `bm_impute_col()`.

## 8. Resolved decisions (locked 2026-05-18)

1. **Default value:** `lambda_mode = "fixed_1"` in v0.10.0. Back-compat-preserving. Revisit flipping to `"estimate"` in v0.11 after a release cycle of stability data. (Aligns with `feedback_options_over_breaking_changes.md`.)
2. **API name:** Two-tier naming:
   - Kernel layer (`bm_impute_col`, `fit_joint_mvn_baseline`): `lambda` — takes numeric `[0, 1]` OR `"estimate"`.
   - Dispatcher / user-facing layer (`fit_baseline`, `fit_pigauto`, `multi_impute`, `impute`): `lambda_mode` — enum `c("fixed_1", "estimate")`.
3. **Joint MVN λ-mode:** Trust phylopars's `model = "lambda"` (one joint fit returns per-trait λ̂_t along with Σ̂). Per-column ML λ is the fallback on phylopars convergence failure. Timing instrumented so we know when phylopars is slow.
4. **Threshold-joint and OVR paths:** `lambda_mode` flows all the way through. Discrete traits (binary / categorical / ordinal liability) inherit λ-estimate automatically as a free upgrade — these paths are deterministic delegates of the joint MVN.

## 9. Success criteria summary

- Hard: BIEN h2h `sla` z-scored RMSE ≤ 0.75 (BACE-level)
- Soft: BIEN_RMSE_TOTAL ≤ 6.85 (closes 80% of the sla-only gap)
- Back-compat: all v0.9.4 tests still pass with `lambda_mode = "fixed_1"` (default)
- AVONET regression check: < 5% RMSE change on continuous traits with `lambda_mode = "estimate"` (because AVONET has strong phylo signal; λ should converge near 1)
- R CMD check clean
- Tests added: 5+ new unit tests, 1 smoke test, 1 verification benchmark
