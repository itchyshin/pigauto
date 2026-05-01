# pigauto — future work register

> Persistent notes on follow-up directions raised during the GNN-fix work
> (Apr 25-26, 2026). Each item explicitly out of scope for the current
> simulation study (see `useful/paper_section_draft.md`,
> `useful/fix_G_real_data_verdict.md`) but worth picking up after.

> **Update 2026-04-26 PM:** items below are reordered after the
> `multi_obs_row_alignment_diagnosis.md` (commit a3e6d39) was found and
> shipped. The smoke-tier "pigauto loses 6/8 cells" verdict was caused
> by that bug, not by the architectural concerns originally listed
> here. Re-running smoke tier in progress; medium tier queued for
> overnight. New top item: pigauto's own multi-obs simulator. The
> transformer-vs-GNN ablation and Pagel's λ items still stand on their
> own merits but are no longer "needed to explain results".

## 1. Transformer vs GNN architecture comparison

**User raised 2026-04-26 (during sim study running):** "once we do this
whether the use of transformer rather than GNN can improve things."

### Context
pigauto's current architecture (post Fix A-H) already uses a **graph
transformer** (`R/graph_transformer_block.R`) with multi-head self-attention
plus a learned phylogenetic Gaussian-bandwidth bias on attention scores.
So the question isn't strictly "transformer vs GNN" — both branches are
transformer-flavoured. The question is more nuanced:

### Concrete hypotheses to test

1. **Pure self-attention (no phylo bias)** vs current rate-aware bias.
   Does removing the phylo prior on attention scores hurt a lot, or
   does the model relearn the prior from data given enough $n$? This
   tests whether the inductive bias from the cophenetic-distance
   Gaussian is doing real work.

2. **Cross-attention with covariates as separate tokens.** Currently
   covariates enter via `cov_encoder + obs_refine` (additive residual)
   and via `cov_inject` per-layer (Fix H). A pure transformer
   alternative: treat each covariate as a separate token, let the
   trait token attend to it via cross-attention. This is the
   "perceiver"-style architecture and may scale better with many
   covariates.

3. **Pre-trained foundation transformer.** Train a single graph
   transformer on a large corpus of phylogenies (TimeTree, Open Tree,
   etc.) without trait labels (mask-language-model style on simulated
   BM traits), then fine-tune on the user's tree+traits. This is the
   "self-supervised pretraining" path mentioned in CLAUDE.md as
   "future work option B". Big lift potential but big infrastructure
   investment.

4. **Set transformer / DeepSets** as alternative to scatter_mean
   pooling in multi-obs mode. Currently obs-to-species pooling uses
   plain mean. A learned attention pool over within-species
   observations might capture heterogeneous within-species structure
   (e.g., when one observation is more informative than others due
   to acclimation context).

### Recommended testing approach

- Treat as a sweep dimension in the sim_bace_pigauto bench. Add a
  `model_arch` argument to `pigauto::impute()` (or a new exported
  helper) with options:
  - `"graph_transformer"` (current default)
  - `"plain_transformer"` (no phylo bias, otherwise identical)
  - `"perceiver_cov"` (covariates as cross-attention tokens)
- Run on the same DGP grid. Add as Tier-4 of the comprehensive sim.
- Decisive test: at strong covariate effect ($\beta = 0.6$) and weak
  phylo signal ($\alpha = 0.2$), does pure self-attention catch up
  to graph-transformer? If YES, the phylo bias is decorative. If
  NO, the bias is doing real work.

### Files that would change

- `R/graph_transformer_block.R` — add `bias_type` parameter
  (`"gauss" | "log_adj" | "none"`).
- `R/model_residual_dae.R` — pass through architecture choice.
- New `R/perceiver_cov_block.R` if pursuing option (2).

### Estimated effort

- Option 1 (no phylo bias): trivial — 1-line change to skip bias term.
- Option 2 (perceiver-style): moderate — new module, ~200 lines.
- Option 3 (pretrained foundation): major — separate pretraining
  infrastructure, ~2-4 weeks.
- Option 4 (set-transformer pooling): moderate — ~150 lines.

**Priority:** start with option 1 (1-line ablation). It will
disambiguate "graph-aware-attention helps vs neutral" cleanly.

---

## 1b. pigauto's own multi-obs simulator (`sim_pigauto_multi_obs()`)

**Approved 2026-04-26 by user.** Spec lives in
`useful/multi_obs_simulator_spec.md`. Build Phase 1 (continuous
response, single predictor, four-share variance decomposition) as a
deliberate package feature, exported with full roxygen + tests.

### Why now

- Independence: pigauto should not rely on `BACE::sim_bace` for its own
  sim studies (BACE is `Suggests:`-only).
- Transparency: the four-share decomposition (`phylo_signal`,
  `species_re_share`, `within_species_share`, residual) is an
  immediately-readable knob for ecologists.
- Coverage: BACE's DGP produces predictors with species + phylo + obs
  noise. There's no easy way to make a `bench_multi_obs.R`-style pure
  within-obs predictor (acclim_temp-like) inside BACE. Our simulator
  should make both regimes one-line each.
- The row-alignment bug confirmed the value of an in-package simulator
  that we can audit end-to-end, instead of having BACE be a black-box
  source of truth.

### Effort

Phase 1: ~half a day (continuous response only, single predictor,
fixed obs/species, gaussian only).
Phase 2-5: ~3 more days (multi-predictor, mixed types, interactions,
vignette).

---

## 1c. Pagel's λ in BM baseline (DEFERRED — wrong original motivation)

**Investigated 2026-04-26 as a potential fix for the multi-obs
"pigauto worse than column-mean on BACE-sim" symptom.** Turned out to
be a spurious motivation: a monkey-patched λ=0 (R = identity) gave
identical RMSE to λ=1 (R = full BM) on the failing cell. The actual
bug was row-alignment, not phylo shrinkage.

### Status

Building λ-aware `bm_impute_col` is still defensible work, just not
urgent. After post-fix smoke tier finishes, we'll know whether
pigauto's BM baseline still has any visible deficit. If yes, this is
back on the table.

### Concrete plan if revived

- Add `lambda ∈ [0, 1]` to `bm_impute_col` and `bm_impute_col_with_cov`,
  with `R_λ = λ R + (1-λ) I`. Default `lambda = 1` for backward compat.
- Add internal `estimate_pagel_lambda(y, R)` — REML profile over a
  grid `c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)`, optionally refined via
  `optimize()`.
- Wire into `fit_baseline.R`: estimate λ per BM column when ≥ 10
  observed species.
- Tests: at λ=0 baseline reduces to global mean; at λ=1 matches
  current behaviour byte-for-byte.

### Effort

~half a day TDD if needed.

---

## 2. LRT threshold tuning per dataset

**Raised in real-data bench analysis (2026-04-26).** The default LRT
threshold (0.02) caused PanTHERIA PopDensity to fit spurious GLS
coefficients (29 % regression), and may have rejected Body_size_mm's
useful coefficients (lost 17 % lift). A per-trait or
cross-validated threshold would address this.

### Approach options

- **F-statistic-based gate** instead of relative-variance reduction:
  if $F$ on the cov coefficient block is significant at $p < 0.10$,
  use cov-aware fit; else fall back.
- **Cross-validated threshold**: split observed data into $K$ folds,
  pick the threshold that minimises val RMSE per trait.
- **Per-trait threshold**: rather than one global default, tune
  per-trait based on data characteristics (n_obs, signal-to-noise).

### Files

- `R/bm_internal.R::bm_impute_col_with_cov` — add `lrt_method` parameter.

---

## 3. Mixed-types covariate test

Day-1 sim (continuous response) and the BACE-sim study cover gaussian
+ binary + threshold response types. But the unique pigauto advantage
that no analytical method offers is **multiple mixed-type traits
imputed jointly with shared covariates**.

### Concrete test

- Generate 4 traits: 2 gaussian, 1 binary, 1 categorical (K=3).
- Each trait has its own phylo signal and beta on the covariates.
- 30 % MCAR per trait.
- Compare pigauto+covariates against per-type analytical methods:
  phylolm-lambda for the gaussian traits, phylogenetic label
  propagation for binary/categorical.
- Question: does pigauto's joint imputation with cross-trait
  attention (the GNN's main role) improve over per-type fits?

This is where pigauto's architecture should genuinely shine, since
no other method handles all four trait types in one fit.

---

## 4. Multi-obs covariate-lift on real data

The bundled `bench_multi_obs.R` shows 10–19 % covariate lift on simulated
multi-obs CTmax. Pre-Fix-A-H. Re-run with current code (Fix H
specifically targets multi-obs cov injection).

Then test on a real multi-obs dataset:
- AnAge longevity (multiple cohort estimates per species)
- Marine fish thermal preference (multiple geographic measurements)

These are pigauto's true architectural home turf where neither
phylolm nor Rphylopars can natively apply.

---

## 5. Real-data multi-trait shared-covariate sim

Take a real comparative dataset (AVONET, PanTHERIA, AmphiBIO) with
multiple traits. Inject a synthetic covariate that affects 2 traits
the same way and 1 differently. See if pigauto correctly extracts
the trait-specific covariate sensitivities.
