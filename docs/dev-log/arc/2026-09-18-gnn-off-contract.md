# gnn-off contract (frozen 2026-09-18; branch feat/gnn-off off origin/main 7f2aff3, v0.11.0)

Changes to this contract go through the orchestrator; builders do not drift it.

Amended after V2 (Rose): no-split fallback = pure baseline (r_bm = 1, no conformal); user covariates
are ignored under gnn = FALSE with a warning; `gnn` is validated as a single TRUE/FALSE; `latent_runs`
stores the raw blend at every cell (observed cells not overwritten), matching GNN-on.

## Semantics of `gnn = FALSE`

- `fit_pigauto(gnn = TRUE, baseline_full = NULL, ...)` and `impute(gnn = TRUE, ...)`. When `gnn = FALSE`
  the GNN is never constructed or trained and the function makes **zero `torch::` calls** (no
  `get_device()`, no `torch_manual_seed()`, no tensors). Everything below happens on plain R matrices.
- The gate calibration keeps GNN-on semantics: `calibrate_gates()` runs with `delta_cal = mu_cal`
  (the GNN corner degenerates to BM). Afterwards `r_cal_gnn` is folded into `r_cal_bm`
  (`r_cal_bm <- r_cal_bm + r_cal_gnn; r_cal_gnn[] <- 0`). `safety_floor` and `phylo_signal_gate`
  behave exactly as with the GNN. The pure traditional-stats arm is
  `gnn = FALSE, safety_floor = FALSE, phylo_signal_gate = FALSE` (existing args).
- Conformal scores: `compute_conformal_scores(trait_map, calibrated_gates = r_cal_gnn (zeros), mu_cal,
  delta_cal = mu_cal, r_cal_bm, r_cal_gnn, r_cal_mean, mean_baseline_per_col, fixed_cal = NULL,
  X_truth_r, val_mask_mat = val_mask_conf, method = conformal_method, bootstrap_B)`. Methods "split"
  and "bootstrap" supported; "mondrian" stops with a clear message under `gnn = FALSE`.
- Two baselines: `baseline` = the val+test-masked fit (unchanged; every scorer reads it);
  `baseline_full` = `fit_baseline(data, tree, splits = NULL, graph = graph, <same args>)`. `impute()`
  computes both BEFORE `graph$D <- NULL` and passes both in; `fit_pigauto()` computes `baseline_full`
  itself only when `gnn = FALSE` and `baseline_full` is NULL.
- `fit_baseline()` returns `path`: a named character vector (one entry per trait_map name) naming the
  dispatch that produced it: "joint_mvn", "threshold_joint", "ovr_categorical", "per_column_bm",
  "label_propagation", "multi_proportion_bm", "zi_gate_lp", or "zi_mag_constant" (as implemented in
  S1 and documented in the roxygen @return; zi_count reports the gate column's dispatch).

## `pigauto_fit` slots under `gnn = FALSE`

| slot | value |
|---|---|
| `model_state` | `list()` (length 0, never NULL) |
| `model_config` | all existing fields (`hidden_dim`, `k_eigen` (= ncol(graph$coords)), `n_gnn_layers`, `gate_cap`, `use_attention`, `use_transformer_blocks`, `n_heads`, `ffn_mult`, `use_trait_attention`, `n_trait_heads`, `trait_embed_dim`, `lambda_mode`, `dropout`, `refine_steps`, `cal_refine_steps`, `train_mask_heldout`, `cov_dim`, `input_dim = p`, `per_column_rs = TRUE`, `n_user_cov`, `seed`, plus whatever main v0.11.0 stores: `joint_solver`, `predict_method`, `joint_refine_iter`, …) **plus `gnn = FALSE`**. GNN-on fits get `gnn = TRUE`. |
| `graph` | as today (D stripped) |
| `baseline` | held-out `list(mu, se, path, ...)` |
| `baseline_full` | `list(mu, se, path, ...)` fit with `splits = NULL`; `NULL` on GNN-on fits |
| `calibrated_gates`, `r_cal_gnn` | named numeric `rep(0, p)` (names = latent names) |
| `r_cal_bm`, `r_cal_mean`, `mean_baseline_per_col` | from `calibrate_gates()` as above (length p; `mean_baseline_per_col` NULL when `safety_floor = FALSE`) |
| `val_rmse`, `test_rmse` | numeric scalar (plain-R RMSE of the blended latent prediction on val / test cells; `NA_real_` when no split); never NULL |
| `history` | 0-row data.frame with the same columns as a GNN-on history |
| `conformal_scores`, `conformal_method` | as computed; `conformal_mondrian = NULL` |
| everything else | identical to the GNN-on fit (`norm`, `species_names`, `obs_*`, `X_scaled`, `n_species`, `n_obs`, `multi_obs`, `trait_names`, `latent_names`, `trait_map`, `splits`, `safety_floor`, `phylo_signal_*`, `covariates`, `cov_*`) |

Shared builder: factor the existing `structure(list(...), class = "pigauto_fit")` block into an
internal `build_pigauto_fit(...)` used by both paths, so slot names cannot diverge.

## `predict.pigauto_fit` under `gnn = FALSE`

- `gnn_off <- isFALSE(cfg$gnn)`. When TRUE: no `ResidualPhyloDAE`, no `load_state_dict`, no torch.
- Effective baseline: `baseline_override` if supplied; else `object$baseline_full` iff
  `is.null(.mask_observed_idx)` (production mode); else `object$baseline` (evaluation mode).
- Latent prediction per column: `pred = r_bm * MU + r_mean * mean_baseline_per_col` (r_gnn = 0),
  matching the three-way blend; observed cells (per `X_scaled`, minus `.mask_observed_idx`) are
  overwritten with observed values, exactly as the GNN path's `X_seed` logic does.
- `n_imputations > 1`: draw `m` latent matrices `MU_draw = MU + rnorm * se` (se from the effective
  baseline, expanded to obs level in multi-obs; se = 0 at observed cells), blended the same way,
  observed cells restored; RNG seeded with `cfg$seed + 100000L` via `withr::with_seed` or a local
  `set.seed` that restores the previous RNG state.
- `rs_val <- rep(0, p)`. Everything downstream (decode, PMM, pooling, SE, conformal intervals,
  `imputed_datasets`) is the existing code, untouched.
- `plot_pigauto.R`'s "No training history available" message says `gnn = FALSE` when that is why.

## Wiring

- `impute(gnn = TRUE)`: forward to `fit_pigauto(gnn = gnn, baseline_full = baseline_full, ...)`; skip
  the `torch::cuda_is_available()` reclaim calls when `gnn` is FALSE.
- `multi_impute(draws_method = "mc_dropout")` on a GNN-off run: one `message()` "gnn = FALSE:
  mc_dropout draws are BM-posterior draws (no dropout)"; no error.
- `multi_impute_trees(gnn = FALSE)`: per-tree `baseline_t` fit with `splits = NULL` and passed as
  `baseline_override`; the baseline-arg replay from `model_config` includes `joint_solver`,
  `predict_method`, `joint_refine_iter` (pre-existing omission, now load-bearing).
- `summary.pigauto_fit` / `print.pigauto_fit` print one line `GNN: off (baseline only)` when
  `cfg$gnn` is FALSE. `evaluate()` unchanged.

## Style
Match surrounding code; roxygen on every changed export; no new dependencies; no edits outside the
leaf's OWNS globs (`.unlazy/gnn-off/gates/leaf-*.md`).
