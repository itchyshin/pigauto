# share_gnn = TRUE as default in `multi_impute_trees()`

**Status:** Approved
**Author:** Claude Opus 4.7 (working with Shinichi Nakagawa)
**Date:** 2026-04-18

---

## Problem

`multi_impute_trees(traits, trees, m_per_tree, ...)` currently runs the
full pigauto pipeline — preprocess → graph → baseline → GNN training →
predict — once per posterior tree. At n=10,000 species with T=50 trees,
that is 17–33 hours, a cost that puts the Nakagawa & de Villemereuil
(2019) tree-uncertainty workflow out of reach for most users.

The expensive part is GNN training. Everything else is cheap relative
to it.

## Goal

Make tree-uncertainty-aware multiple imputation cheap enough that the
N&dV 2019 canonical workflow (T=50 posterior trees, one imputation per
tree, M=50 downstream model fits pooled with Rubin's rules) is the
default path.

## Non-goals

- No change to `multi_impute()` (single-tree MI). Users with one tree
  stay on that function unchanged.
- No new scientific method. This is an engineering change that reuses
  the existing gate-safety property: when phylogenetic signal is strong
  the gate closes, the GNN contributes nothing, and predictions fall
  back to the per-tree baseline — which is already tree-uncertainty
  aware by construction.

## Design decisions

### 1. `share_gnn = TRUE` is the new default

Under `share_gnn = TRUE`:

1. Pick a reference tree (see §2).
2. Run the full pigauto pipeline ONCE on the reference tree. Save the
   fitted GNN (weights), its spectral features, its calibrated gate
   `r_cal`, and its conformal scores.
3. For each posterior tree `tree_t`:
   - Recompute the BM / joint-MVN baseline using `tree_t`'s phylo
     covariance. This is the only per-tree work.
   - Predict by blending: `pred_t = (1 - r_cal) * baseline_t + r_cal *
     gnn_shared`. The GNN contribution is a tree-invariant constant
     computed once on the reference tree.
   - Draw `m_per_tree` imputations from `pred_t`.

Under `share_gnn = FALSE` (opt-out):
The full pipeline runs per tree, as today.

### 2. Reference tree selection: MCC via `phangorn`, fallback to `trees[[1]]`

New argument `reference_tree = NULL`:
- `NULL` (default): compute the **maximum clade credibility** tree
  from `trees` via `phangorn::maxCladeCred()`.
- If `phangorn` is not installed: warn (`"phangorn required for MCC
  tree selection; falling back to trees[[1]]. Install with
  install.packages('phangorn')"`) and use `trees[[1]]`.
- `reference_tree = <phylo>`: use the user-supplied tree directly.

`phangorn` is added to `Suggests` in `DESCRIPTION`.

### 3. `m_per_tree` default flips from 5L to 1L

Rationale: under `share_gnn = TRUE`, between-tree variance already
carries tree uncertainty, so per-tree over-sampling is redundant. T=50
trees × m_per_tree=1 = M=50, which matches the N&dV 2019 canonical
workflow exactly.

Users with T < 20 trees will be warned at runtime:
`"Low total imputations (T * m_per_tree = {M}). Consider m_per_tree = 5L
for stable Rubin's rules pooling when T < 20."`

### 4. Single-tree users are unaffected

`multi_impute()` remains the entry point for single-tree MI. The
existing error message in `multi_impute_trees()` already directs
single-tree users to the right place:

> "trees must contain at least 2 phylogenies for tree uncertainty
> propagation. For single-tree MI, use multi_impute() instead."

No behaviour change there.

### 5. Correctness claims (tree-uncertainty propagation)

| Gate state | Prediction | Tree uncertainty |
|---|---|---|
| Closed (r_cal ≈ 0) | `pred = baseline(tree_t)` | Fully preserved — GNN contribution is zero, baseline varies per tree |
| Partially open (r_cal ≈ 0.3) | `0.7·baseline_t + 0.3·gnn_shared` | 70% preserved — baseline portion varies, GNN portion is constant |
| Fully open (r_cal ≈ 1) | `pred ≈ gnn_shared` | Lost in GNN channel; still present in baseline channel which contributes (1-r_cal) |

On every real dataset benchmarked during the v0.9.0 campaign the gate
closed partially or fully (Delhey, binary, categorical, BM-continuous).
The bias direction is conservative: tree uncertainty is slightly
under-estimated, not over-estimated.

## API surface

```r
multi_impute_trees(
  traits, trees,
  m_per_tree = 1L,              # was 5L
  species_col = NULL,
  trait_types = NULL,
  multi_proportion_groups = NULL,
  log_transform = TRUE,
  missing_frac = 0.25,
  covariates = NULL,
  epochs = 2000L, verbose = TRUE,
  seed = 1L,
  share_gnn = TRUE,             # NEW
  reference_tree = NULL,        # NEW
  ...
)
```

Return object (`pigauto_mi_trees`) gains:
- `share_gnn`: logical (whether sharing was used)
- `reference_tree`: the tree actually used for training (resolved from
  MCC / user / fallback)
- `fit`: single `pigauto_fit` when `share_gnn = TRUE`; existing
  `fits` list of T fits when `share_gnn = FALSE`

## Implementation sketch

New internal helper in `R/multi_impute_trees.R`:

```r
predict_on_new_tree <- function(fit, new_tree, data, splits, graph_ref) {
  # Recompute baseline on new_tree only
  baseline_new <- fit_baseline(data, new_tree, splits = splits,
                               graph = graph_ref)
  # Reuse trained GNN delta from fit (tree-invariant under share_gnn)
  # Blend with calibrated gate
  predict(fit, baseline_override = baseline_new)
}
```

Requires modest surgery to `predict.pigauto_fit()` to accept an
override baseline while keeping all other state from the fit. That is
a single-argument addition, not a rewrite.

MCC tree helper in `R/multi_impute_trees.R`:

```r
resolve_reference_tree <- function(trees, reference_tree = NULL) {
  if (!is.null(reference_tree)) return(reference_tree)
  if (requireNamespace("phangorn", quietly = TRUE)) {
    return(phangorn::maxCladeCred(trees))
  }
  warning("phangorn required for MCC tree selection; falling back to ",
          "trees[[1]]. Install with install.packages('phangorn') for ",
          "MCC-based selection.", call. = FALSE)
  trees[[1]]
}
```

## Testing

1. **Unit test** (`tests/testthat/test-share-gnn.R`, new):
   - tiny synthetic dataset (n=20, T=3 trees)
   - `share_gnn = TRUE` produces a `pigauto_mi_trees` with correct shape
   - `share_gnn = FALSE` path still works (regression guard)
   - `reference_tree = <phylo>` is honoured
   - phangorn-missing fallback triggers warning
   - gate-closed case: predictions match per-tree BM baseline exactly

2. **Smoke benchmark** (`script/bench_share_gnn.R`, new):
   - Non-blocking, produces a small table: wall time and pooled FMI for
     `share_gnn = TRUE` vs `FALSE` on AVONET 300 with T=10
   - Expected: ~10× speedup, FMI within ±5% of the slow path
   - Output linked from the new pkgdown article

3. **Existing tree_uncertainty bench** (`script/bench_tree_uncertainty.R`):
   - Keeps `share_gnn = FALSE` explicitly to preserve the v0.9.0
     baseline comparison. We do NOT re-run it with new defaults.

## Documentation deliverables

1. `R/multi_impute_trees.R` roxygen:
   - New `@param share_gnn` and `@param reference_tree` entries.
   - New `@section When to use this:` block with the one-tree /
     many-tree decision guide.
   - Updated `@examples` shows T=50, m_per_tree=1 workflow and the
     Map() over `mi$tree_index` for the Rubin's-rules step.
2. `R/multi_impute.R` roxygen:
   - Same `@section When to use this:` block (for discoverability).
3. New pkgdown article `vignettes/tree-uncertainty.Rmd`, rendered at
   `articles/tree-uncertainty.html`. Navbar entry under Articles.
   Contents:
   - "Do I need this article?" — decision guide
   - What tree uncertainty is
   - Why the gate-closed regime makes `share_gnn = TRUE` safe
   - End-to-end worked example with BirdTree posteriors + Rubin's rules
   - Timing table from the smoke benchmark
4. `NEWS.md` entry for v0.9.1 (or v0.10.0 if we also land Phase 6 EM).

## Backward compatibility

- Existing user code running `multi_impute_trees(df, trees, m_per_tree=5)`
  continues to work but now runs ~10× faster with the new defaults.
- Users who need exact numerical equivalence to the pre-v0.9.1 slow
  path set `share_gnn = FALSE`.
- The `m_per_tree = 5 → 1` default change means users who call
  `multi_impute_trees(df, trees)` without specifying `m_per_tree` will
  get fewer datasets. For T ≥ 20 this is desirable (canonical workflow).
  For T < 20 the runtime warning points them to bump `m_per_tree`.

## Out of scope (deferred)

- `share_gnn = "weights_only"` (share weights, recompute per-tree
  graph): considered and deferred. Adds code complexity for a
  correctness benefit that is small when the gate closes. Revisit if
  empirical FMI loss exceeds 5% on any real dataset.
- Phase 6 EM for threshold baseline: separate design, not this spec.
- Multi-observation + tree uncertainty: the current
  `multi_impute_trees()` does not support multi-obs mode; this spec
  does not change that. Separate follow-up if requested.

## Success criteria

- [ ] `multi_impute_trees(df, trees)` with default arguments completes
      in ~10% of the pre-v0.9.1 wall time on AVONET 300 × T=50.
- [ ] Pooled FMI from `share_gnn = TRUE` is within ±5% of
      `share_gnn = FALSE` on the smoke benchmark.
- [ ] `phangorn`-missing environment triggers the fallback warning
      cleanly; rest of the pipeline proceeds.
- [ ] New pkgdown article renders without errors; navbar link works.
- [ ] Test suite stays green (558 → ~565 tests with the new unit tests).
- [ ] `R CMD check` stays clean.
