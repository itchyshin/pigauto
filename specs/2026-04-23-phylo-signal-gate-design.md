# Phylogenetic-signal gate — design spec

**Status:** approved (brainstorming 2026-04-23).
**Branch:** `feature/phylo-signal-gate` (stacked on `feature/safety-floor-mean-gate`).
**Plan:** `plans/2026-04-23-phylo-signal-gate.md` (to be written next).
**Stacks on:** `specs/2026-04-23-safety-floor-mean-gate-design.md` (PR #43).

---

## 1. Goal

Add a per-trait diagnostic that explicitly reports whether each trait has
enough phylogenetic signal for the Brownian-motion baseline to help, and
if not, forces the calibrated blend to the grand-mean corner
`(r_bm = 0, r_gnn = 0, r_mean = 1)` before any 3-way simplex search.
Makes the "pigauto knows when phylogeny helps" claim explicit and
scientifically interpretable.

## 2. Background

PR #43 (safety floor) ensures `pigauto_val_RMSE ≤ mean_val_RMSE` by
letting the calibrator pick `r_mean > 0` whenever the grand mean beats
the BM/GNN blend on validation. This works but is *implicit* — the
user sees `fit$r_cal_mean = 0.8` without any explanatory "why".

The phylogenetic-signal gate replaces the implicit mechanism with an
*explicit* upstream decision per trait: compute Pagel's λ on the
training cells; if `λ < threshold`, this trait is flagged weak-signal
and the blend goes straight to grand mean. No wasted BM fit, no
wasted GNN training effort on that trait. The fit object carries
`phylo_signal_per_trait` (named numeric of λ values) and
`phylo_gate_triggered` (named logical) so users / reports can
see which traits were routed this way.

**Scientific framing.** On plants with BIEN × V.PhyloMaker2 we already
know from the bench (commit `16ce036`) that four of five traits have
`r ≤ 0.21` — not because pigauto is broken, but because the
phylogenetic signal on pooled species means is weak. The phylo-signal
gate puts that fact on the label, as a Pagel's λ diagnostic reported
per trait.

## 3. Non-goals

- Changing the safety-floor calibration logic (PR #43 stays).
- Adding Blomberg's K as an alternative threshold (future; the spec
  accepts a `method = "lambda"` default with method-dispatch room, but
  only the `lambda` path ships in this PR).
- Environmental covariates (separate spec: GBIF centroids).

## 4. Architecture

### 4.1 Two-stage decision

```
for each trait t:
    λ_t ← phytools::phylosig(tree, trait_values_t, method="lambda")$lambda
    if λ_t < threshold (default 0.2):
        r_cal_bm_t = 0
        r_cal_gnn_t = 0
        r_cal_mean_t = 1   (forced, no calibration search)
        skip BM fit on trait t   (compute saving)
        skip GNN training contribution on trait t (via mask on loss)
    else:
        normal 3-way simplex calibration (PR #43 path)
```

### 4.2 Latent-column-to-trait mapping

A single trait may have multiple latent columns (categorical, zi_count,
multi_proportion). The gate decision applies to the **whole trait**,
not per-latent-column. For a categorical trait, if λ on the dominant
class-frequency signal is weak, ALL K latent columns get
`(0, 0, 1)`.

For multi-obs mode, λ is computed on the species-level aggregated
trait values (same aggregation rules as `aggregate_to_species()` in
`R/joint_threshold_baseline.R`).

### 4.3 When to compute λ

At the top of `fit_pigauto()` (before `fit_baseline()`). Needs: tree
+ the training-observed slice of `data$X_scaled` (or the raw values
for discrete traits). One-time cost per fit: `phylosig()` on an
n-tip tree is O(n²) for the REML likelihood; at n=10k this is
~2-5s per trait. Cached on the fit via `fit$phylo_signal_per_trait`.

### 4.4 Discrete traits

`phytools::phylosig()` supports `method = "lambda"` only for
continuous input. For binary / categorical / ordinal traits we
compute λ on the z-scored LATENT liability (which is continuous)
via:

- binary: λ on `estep_liability_binary(observed, mu_prior = 0,
  sd_prior = 1)$mu`.
- categorical: max of λ over K one-vs-rest liability columns.
- ordinal: λ on the z-scored integer class (already continuous in
  the latent matrix).

The zi_count gate column uses the binary path; the magnitude column
uses the continuous path. If their λs disagree (gate λ > threshold,
magnitude λ < threshold or vice versa), treat the WHOLE trait as
weak-signal iff BOTH columns are weak. This matches the Phase 5
coupled-gate convention in pigauto.

### 4.5 Storage on `pigauto_fit`

```r
fit$phylo_signal_per_trait        # named numeric, length = length(trait_map)
                                  # value is Pagel's λ (or NA if computation failed)
fit$phylo_gate_triggered          # named logical, same length
                                  # TRUE = this trait was routed to grand-mean
fit$phylo_signal_method           # character scalar, "lambda" or (future) "blomberg_k"
fit$phylo_signal_threshold        # numeric scalar, the threshold used
```

Traits with `phylo_gate_triggered[t] = TRUE` have
`fit$r_cal_bm[latent_cols(t)] = 0`,
`fit$r_cal_gnn[latent_cols(t)] = 0`,
`fit$r_cal_mean[latent_cols(t)] = 1` by construction — the safety
floor's stored weights are consistent with the gate decision.

### 4.6 Fit path optimisation

When a trait is gated, skip:

- its row in `fit_baseline()` (no BM conditional-MVN compute).
- its GNN loss contribution (`loss_mask[:, gated_latent_cols] = 0`
  in the per-batch loss).
- its entry in the calibration 3-way search (forced to `(0,0,1)`).

These are performance optimisations, not correctness; if any is
hard to implement they can be deferred to a follow-up without
breaking the spec.

### 4.7 API

```r
impute(..., phylo_signal_gate = TRUE,
            phylo_signal_threshold = 0.2,
            phylo_signal_method = "lambda")
fit_pigauto(..., phylo_signal_gate = TRUE,
                 phylo_signal_threshold = 0.2,
                 phylo_signal_method = "lambda")
```

Default: `phylo_signal_gate = TRUE` in v0.9.1.9003.

`phylo_signal_gate = FALSE` reproduces #43's behaviour (safety floor
only, no phylo-signal pre-filter).

### 4.8 Reporting

`print.pigauto_fit()` gains one line when any gate triggered:

```
Phylogenetic signal (Pagel's λ, threshold 0.2):
  gated (BM skipped → grand mean): sla (λ=0.04), leaf_area (λ=0.11),
  seed_mass (λ=0.09), height_m (λ=0.17)
  passed (BM + GNN active):        wood_density (λ=0.38)
```

If zero traits are gated, print a compact one-liner: `"All 5 traits
have phylogenetic signal ≥ 0.2; safety floor remains the fallback."`

## 5. Testing

### 5.1 Smoke canary (new file: `tests/testthat/test-phylo-signal-gate.R`)

1. **`phylosig_per_trait()` returns λ for each trait.** Synthetic: 3
   traits with known λ ∈ {0.05, 0.5, 0.95}; verify returned values
   are within ±0.1 of truth at n=200.

2. **Gate triggers iff λ < threshold.** Same synthetic data;
   threshold = 0.2 should gate trait 1 (λ=0.05), not trait 2 or 3.

3. **Gated trait uses pure-mean prediction.** Fit with gate ON;
   verify `fit$r_cal_mean[latent_cols(trait_1)] = 1` exactly, and
   predictions equal `mean_baseline_per_col` (on latent scale).

4. **AVONET 300 regression: no trait gated.** All AVONET continuous
   traits have λ > 0.9; gate does not trigger; predictions
   bit-identical to a `phylo_signal_gate = FALSE` run.

5. **Plants smoke: weak-signal traits gated.** 300-species BIEN
   subset (cached data); verify at least 2 of 5 traits flagged
   (height_m and one of sla/leaf_area/seed_mass).

6. **Binary / categorical path via latent liability.** Synthetic
   binary trait with strong-signal and weak-signal variants; verify
   threshold fires correctly on the latent liability.

7. **Backward compat: `phylo_signal_gate = FALSE`.** Reproduces #43
   exactly. Same `r_cal_*` slots.

8. **Print method shows gate info.** `capture.output(print(fit))`
   contains the "Phylogenetic signal" section when gate triggered.

### 5.2 Full canary update

Extend `script/regress.R` to include a per-trait λ readout in the
JSON output. No new hard gates (the invariant is that the fit
produces finite, safety-floored predictions; phylo-signal diagnostic
is reporting, not pass/fail).

## 6. Risks and mitigations

### 6.1 `phytools` is a heavy dependency

`phytools` imports a long chain (`maps`, `scatterplot3d`,
`numDeriv`, `optimParallel`, ...). To keep pigauto's Imports light,
add `phytools` to **Suggests**, and degrade gracefully when absent:

```r
if (phylo_signal_gate && !requireNamespace("phytools", quietly = TRUE)) {
  warning("phylo_signal_gate = TRUE requires the 'phytools' package; ",
          "falling back to phylo_signal_gate = FALSE.")
  phylo_signal_gate <- FALSE
}
```

### 6.2 λ is ill-defined for small trees

`phylosig()` is unstable below ~20 tips. Add a `min_tips = 20L`
guard; below that, skip computation and set `λ = NA`, do not
gate.

### 6.3 λ is ill-defined for constant traits

If a trait has zero variance on training cells, λ is NaN.
`phylosig()` already handles this gracefully (returns NA); set
`λ = NA` and don't gate.

### 6.4 Threshold of 0.2 is a judgment call

0.2 is Blomberg et al. 2003's "very weak" cutoff re-expressed in
λ. Configurable via `phylo_signal_threshold`. Document the
default choice in the roxygen with a citation.

### 6.5 Interaction with safety floor when gate triggers at edge

When a trait's λ is right at the threshold, adding noise could
flip the gate on successive fits. Mitigation: document that the
gate is deterministic given fixed training splits + seed; users
seeing borderline λ values can adjust the threshold.

## 7. Rollout

Single PR on `feature/phylo-signal-gate`. Stacks on #43's work.
After #43 merges:

1. Rebase `feature/phylo-signal-gate` onto `main`.
2. Ensure regression canary green.
3. Merge.

Version bump: `0.9.1.9002` → `0.9.1.9003`.

## 8. Success metric

After merge:

- On AVONET / PanTHERIA / FishBase / AmphiBIO: zero traits gated
  (λ > 0.2 everywhere). pigauto_on results bit-identical to #43.
- On plants (BIEN × V.PhyloMaker2): 3+ of 5 traits gated.
  `fit$phylo_gate_triggered` names them. Grand-mean RMSE matches
  what the gated path produces.
- Paper can claim: "pigauto explicitly reports per-trait
  phylogenetic signal (Pagel's λ). When λ < 0.2, BM is skipped and
  the prediction reverts to grand mean + covariates + GNN."
