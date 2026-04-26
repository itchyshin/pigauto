# Multi-obs row-alignment bug — diagnosis & fix

> Found and fixed 2026-04-26.
> Commit: `a3e6d39 fix(multi-obs): row-alignment bug in build_completed when input is shuffled`
> Follow-up commit: `5ee4b5e test(multi-impute): regression test`

## Symptom

Smoke tier of `script/bench_sim_bace_pigauto.R` (run 2026-04-26 morning,
post Fix A-H) showed pigauto **WORSE than the global column-mean baseline**
on all multi-obs cells:

| beta | pred_phylo | column_mean | pigauto_no_cov |
|------|-----------|-------------|----------------|
| 0.0  | 0.0       | 1.71        | 2.30 (34 % worse) |
| 0.0  | 0.6       | 1.73        | 2.38 (38 % worse) |
| 0.5  | 0.6       | 1.67        | 2.38 (43 % worse) |

A column-mean baseline is the trivial "predict the global mean for every
held-out cell" approach. Pigauto being worse than that tells you the
predictions are *systematically wrong*, not merely under-fit.

## Investigation

I tried four hypotheses before finding the real bug.

1. **Pagel's λ shrinkage** — pigauto's BM baseline always assumes 100 %
   phylo signal; thought it might be over-shrinking species means toward
   the phylogenetic prior when species RE dominates. Tested by
   monkey-patching `phylo_cor_matrix` to return identity (= λ=0
   extreme); pigauto RMSE went from 2.23 to 2.24. **Not the bug.**
2. **Predictor multicollinearity** — BACE generates predictors with
   their own phylo signal; thought double-counting through GLS-aware
   baseline. Swept `predictor_phylo_signal ∈ {0, 0.3, 0.6}`; deficit
   appeared at 0.0 too. **Not the bug.**
3. **Fix A-H regression** — re-ran `bench_multi_obs.R` (the bundled
   bench that previously showed 10–19 % cov-lift) with current code.
   **Lifts INCREASED**: 0.81–0.99 → 0.64–0.99. Fix A-H didn't break
   anything. **Not the bug.**
4. **DGP architecture mismatch** — BACE's `sample(species, n_cases,
   replace = TRUE)` produces randomly-ordered observation rows.
   `bench_multi_obs.R` uses `rep(tree$tip.label, each = n_obs_per)`
   which is already tree-tip-sorted. Hypothesis: the row order matters
   for pigauto. **THIS WAS THE BUG.**

The decisive test was inspecting per-row predictions on a single BACE
cell. Pigauto's predictions had correlation 0.0042 with the true species
mean — predictions were essentially random with respect to species
identity. Looking at consecutive rows showed an obvious shift pattern:
two observations of the same species got two completely different
predictions, while different species got similar ones.

## Root cause

`R/preprocess_traits.R` reorders the input data.frame's rows to
tree-tip order so the GNN's `scatter_mean(obs_to_species)` and
`index_select(obs_to_species)` operations work cleanly:

```r
# R/preprocess_traits.R, multi_obs branch
sp_order <- tree$tip.label[tree$tip.label %in% unique_species]
sp_rank  <- setNames(seq_along(sp_order), sp_order)
row_order <- order(sp_rank[obs_species_raw])
traits <- traits[row_order, , drop = FALSE]
obs_species_raw <- obs_species_raw[row_order]
rownames(traits) <- NULL                    # erases input-row provenance
```

The model produces predictions in this internal (tree-tip) order. Then
`R/impute.R::build_completed` merged predictions back into the user's
input data.frame:

```r
# OLD code
} else {
  # Multi-obs: imputed rows align 1:1 with original rows in input order
  imp_row <- seq_len(n_row)
}
```

The comment was wrong — internal order ≠ input order — but the bug was
hidden by `bench_multi_obs.R`, whose DGP happens to be already
tree-tip-sorted (the reorder is a no-op there). The bug only surfaced
when input rows were genuinely shuffled.

## Fix

Three small changes in `R/`, plus regression tests:

1. **`R/preprocess_traits.R`** stores `input_row_order` in the returned
   `pigauto_data` object. Length = `n_obs`. `input_row_order[k] = i`
   means "internal position k holds the i-th original-input row". NA
   entries flag synthetic rows that `preprocess_traits` adds for tree
   tips with no input data.

2. **`R/impute.R::build_completed`** accepts `input_row_order` and uses
   `match(seq_len(n_row), input_row_order)` to invert the reorder. The
   legacy `seq_len(n_row)` and `match(rownames(original),
   rownames(imputed))` paths are kept for any external callers that
   don't pass `input_row_order`, so this is non-breaking.

3. **`R/multi_impute.R`** + **`R/multi_impute_trees.R`** (4 call-sites
   total) pass `input_row_order = res$data$input_row_order` to
   `build_completed`.

The single-obs path of `preprocess_traits` also populates
`input_row_order` (= `match(tree$tip.label, rownames(traits))`). This
is a no-op when input rownames already match tree tips, and produces
equivalent output to the legacy rowname-match path otherwise.

## Verification

### Spot-check on a single BACE-sim multi-obs cell

|                          | before fix | after fix |
|--------------------------|-----------|-----------|
| pigauto RMSE             | 2.23      | **1.31**  |
| vs column_mean (1.89)    | 18 % worse | **31 % better** |
| vs species_mean (1.38)   | 62 % worse | **5 % better** |
| cor(pred, sp_mean)       | 0.004     | **0.986** |

### Full 12-cell BACE-sim sweep (`script/diag_bace_predictor_phylo.R`)

Pigauto / column_mean ratio (< 1.0 = pigauto wins):

| beta | pred_phylo | OLD ratio | NEW ratio |
|------|-----------|-----------|-----------|
| 0.0  | 0.0       | 1.344     | **0.648** |
| 0.0  | 0.3       | 1.360     | **0.736** |
| 0.0  | 0.6       | 1.378     | **0.729** |
| 0.5  | 0.0       | 1.276     | **0.764** |
| 0.5  | 0.3       | 1.351     | **0.772** |
| 0.5  | 0.6       | 1.427     | **0.766** |

Pre-fix: 27–43 % WORSE. Post-fix: 22–35 % BETTER. Median swing: ~62
percentage points. Every cell flipped from clear loss to clear win.

### Test suite

- 1016 tests pass (1 pre-existing failure in
  `test-ovr-categorical.R:94`, identical with or without this fix).
- New regression tests:
  - `test-preprocess.R`: `input_row_order` permutation-recovery semantics.
  - `test-fit-predict.R`: `impute()` with shuffled multi-obs input —
    correlation > 0.95 with species truth, max abs diff < 1.0.
  - `test-multi-impute.R`: `multi_impute()` with shuffled multi-obs
    input, both `conformal` and `mc_dropout` draws methods —
    correlation > 0.85 with species truth in every dataset.

### Bundled `bench_multi_obs.R` (no regression)

Re-running the previously-passing bench against the fix confirms no
regression and shows the cov-lift actually *improved* slightly because
of Fix A-H (commits earlier yesterday, unrelated to row-alignment):

| β | OLD lift | POST-FIX lift |
|---|---|---|
| 0.0 | 0 % | 0 % |
| 0.5 | 12–15 % | **22–25 %** |
| 1.0 | 15–19 % | **28–36 %** |

## What this fixes in production

Anyone who calls `pigauto::impute(traits, tree, species_col = ...)`
with a `traits` data.frame where the species column isn't already in
tree-tip order — i.e., the common case for any real-data ingestion.
The same bug propagates through `multi_impute()` and
`multi_impute_trees()` (Rubin's-rules + tree-uncertainty paths).

## What this does NOT change

- **Single-obs imputation** with rownames matching tree tips. The new
  `input_row_order` is `seq_len(n)` here, no-op behaviour change.
- **Internal latent-space predictions**. The model output was always
  correct; only the merge-back step was broken.
- **Existing on-disk `pigauto_fit` objects.** Model state is untouched.
- **Yesterday's real-data covariate-lift results** (PanTHERIA,
  GlobTherm, AmphiBIO, LepTraits, Delhey). All 5 of those benches use
  the single-obs `impute()` path — confirmed by inspection. Results
  stand as reported.

## Why we caught it now

We caught it because we ran a comprehensive simulation study on a
genuinely-shuffled DGP (`BACE::sim_bace`) for the first time. Every
prior benchmark either pre-sorted input rows in tree-tip order
(`bench_multi_obs.R`, `bench_multi_obs_real_tree.R`,
`bench_multi_obs_mixed.R`) or used the single-obs path (real-data
benches).

The takeaway is to keep simulating beyond familiar DGPs. The smoke
tier was originally meant to confirm pigauto was within 10 % of
phylolm-λ on linear gaussian; instead it surfaced a year-old (estimated)
production bug. Worth its weight.

## Follow-up work registered

- Re-run smoke + medium tiers of `bench_sim_bace_pigauto.R` post-fix to
  produce the comprehensive verdict the original sim study aimed for.
- Build pigauto's own multi-obs simulator (`useful/multi_obs_simulator_spec.md`).
  Independent of this bug, but the bug confirmed the value of having an
  in-package, well-controlled multi-obs simulator that we can audit
  end-to-end.
- Pagel's λ in BM baseline: deferred. The original motivation was
  wrong (bug, not shrinkage), but the feature might still be useful in
  regimes where species RE dominates phylo signal. Verify post-fix
  smoke-tier results before committing to building it.
