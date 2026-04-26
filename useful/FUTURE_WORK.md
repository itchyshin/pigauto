# pigauto — future work register

> Persistent notes on follow-up directions raised during the GNN-fix work
> (Apr 25-26, 2026). Each item explicitly out of scope for the current
> simulation study (see `useful/paper_section_draft.md`,
> `useful/fix_G_real_data_verdict.md`) but worth picking up after.

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
