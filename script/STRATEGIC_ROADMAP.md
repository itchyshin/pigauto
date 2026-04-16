# pigauto strategic roadmap: benchmark datasets, simulation datasets, ML/DL improvements

Drafted 2026-04-15 after Phase 5 lands. Covers three question families the user raised:

1. What benchmark datasets should pigauto target beyond AVONET?
2. What simulation datasets should we build to discriminate methods?
3. How can we improve pigauto's machine-learning / deep-learning components?

The "priority" column reflects expected effort-to-impact ratio: **HIGH** = likely ≥2pp accuracy lift per 1 week of work; **MEDIUM** = ≥2pp per 2–3 weeks; **LOW** = exploratory, no guarantee.

---

## Part 1: Benchmark datasets

pigauto's only production benchmark is **AVONET** (9993 birds × 11 traits, Tobias et al. 2022). BACE also uses it. To validate pigauto as a general-purpose tool, we need a portfolio of biological datasets with different structure characteristics.

### 1.A Recommended additions (priority order)

| Dataset | Taxa | N species | Trait mix | Why it matters | Priority |
|---|---|---|---|---|---|
| **PanTHERIA** | Mammals | 5,416 | 5 continuous + 5 categorical + 3 binary + heavy missingness (~30–60%) | Gold standard for PCM method papers; mixed-type stress test | HIGH |
| **FishBase** core traits | Fishes | ~35,000 | 8 continuous + 4 categorical + 2 binary | Scale test (pigauto's largest validation to date has been ~10k) | HIGH |
| **TRY Plant Trait** | Plants | 300,000+ measurements, 2,000+ traits | Continuous-heavy, extreme sparsity | Stress-test on sparse-but-wide regime (unlike AVONET's dense-narrow) | MEDIUM |
| **AmphiBIO** | Amphibians | 6,776 | 4 continuous + 3 categorical + 5 binary | Smaller, mixed-type, IUCN-rich for evaluation | MEDIUM |
| **EltonTraits** | Birds + Mammals | 9,993 + 5,400 | Diet composition (multi_proportion!) + categorical | Tests `multi_proportion` trait type at scale | MEDIUM |
| **Tree of Sex** | Plants, animals | ~7,000 | Heavy categorical, ordinal | Exercises Phase 4 categorical once Phase 6 EM lands | LOW |

### 1.B Evaluation protocol per dataset

- **MCAR (random missingness)** at 10%, 25%, 50% — the default.
- **MAR (clade-correlated missingness)** — mask entire subclades at a time. Tests whether pigauto's phylogenetic prior actually helps in realistic sampling patterns (real databases are clade-biased, not random).
- **MNAR (trait-biased missingness)** — e.g. mask body-size values specifically for rare species. Tests robustness to non-ignorable missingness.

### 1.C Evaluation metrics (standardised)

- **Continuous**: RMSE, Pearson r, Spearman ρ, R², coverage of 95% conformal intervals
- **Binary / categorical**: accuracy, log-loss, Brier score, top-1 and top-3 accuracy (categorical)
- **Mixed-type**: combined score (e.g. macro-average of normalised metrics)
- **Compute**: wall-clock, peak RAM (already instrumented)

### 1.D Output artefacts

- `script/bench_{dataset}.R` — run pigauto + baselines (LP, BACE, Rphylopars) on same splits
- `script/bench_{dataset}.md` — per-scenario table  
- `pkgdown/assets/dev/bench_{dataset}.html` — linked from validation suite
- `script/compare_v080_bace.R` — aggregate comparison table

---

## Part 2: Simulation datasets

Real benchmarks answer "does this work on current data?" — but simulations answer "*when* does this work, and *when* does it fail?". Simulations are what make a paper rigorous.

### 2.A Recommended simulation studies

| Simulation | What it varies | What it discriminates | Priority |
|---|---|---|---|
| **Signal strength sweep** (Pagel's λ) | λ ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 1.0} | Low-λ: all methods converge to mean imputation. High-λ: phylo-aware wins. Sweet spot: moderate λ where methods discriminate | HIGH |
| **Evolutionary model** | BM / OU / ACDC / rate-heterogeneous | pigauto's BM baseline is optimal for BM but should gracefully degrade on OU; benchmark the gate-closure safety | HIGH |
| **Cross-trait correlation sweep** | ρ ∈ {0, 0.2, 0.4, 0.6, 0.8} across types | Level-C joint MVN should lift at high ρ; per-column BM at low ρ. Already tested for BM in `bench_joint_baseline.R`; extend to mixed types | HIGH |
| **Trait-type mixture** | % continuous vs binary vs categorical | Tests dispatcher logic end-to-end. Good for catching regressions | MEDIUM |
| **Clade-correlated missingness** | MAR pattern: entire subclades at various depths | Critical for real data — current benchmarks are MCAR | HIGH |
| **Rate shifts** | One clade evolves 3× faster | Tests whether attention GNN handles local rate heterogeneity | MEDIUM |
| **Phylogenetic noise** (tree error) | Branch-length perturbation, topology swaps | Already have `bench_tree_uncertainty.R`; extend to quantify breakdown at high tree error | MEDIUM |
| **Observation noise** | Add Gaussian / Laplacian / outlier noise to simulated values | Tests robustness — pigauto currently assumes clean observations | LOW |
| **Compositional data** | Multi_proportion with Dirichlet-BM evolution | Validates multi_proportion trait type against a proper compositional generator | MEDIUM |

### 2.B Simulation framework refactor

Currently `script/bench_*.R` has simulation code duplicated across files. Refactor to `R/sim_traits.R` (exported or `inst/`-only) with a clean API:

```r
sim_traits(
  tree, n_species,
  types        = c("continuous", "binary", "categorical", ...),
  evo_model    = "BM",          # or "OU", "ACDC"
  sigma        = identity_matrix,  # K x K cross-trait cov
  lambda       = 1.0,            # phylogenetic signal
  seed         = ...
)
```

One generator, all simulations. Makes it trivial for users to run their own power studies.

### 2.C Discriminative benchmark (Phase 8 planned)

`script/bench_discriminative.R` — low/mixed/moderate signal scenarios, mixed types, cross-trait correlation, realistic missingness. This is the benchmark BACE or pigauto will win or lose on. Master plan has this for Phase 8; the strategic framing: the benchmark MUST produce a non-flat ranking across methods, otherwise the benchmark doesn't discriminate and we're flying blind.

**Quality bar for the discriminative benchmark:**
- Top vs bottom method spread ≥ 5pp on ≥3 of 5 scenarios
- No method wins every scenario (tradeoffs visible)
- At least one scenario where LP (simple baseline) ties pigauto (shows pigauto doesn't regress)
- At least one scenario where BACE beats pigauto (honest about limitations)

---

## Part 3: ML / DL improvements

Phase A work from the old plans (curriculum masking, clade masking, SWA) showed flat results. Those were training-procedure tweaks. The higher-leverage changes are architectural and data-centric.

### 3.A Architecture improvements

**a) Graph Transformer + PhyloBias (HIGH priority)**

Current `ResidualPhyloDAE` uses scaled dot-product attention with a learnable log-adjacency bias. This is close to a graph transformer but uses only 2 layers. Scaling up:

- 4–6 transformer layers (was 2)
- Multi-head attention (was single-head)
- PhyloBias = `log(adjacency) + rate_feature` where `rate_feature` encodes local evolutionary rate per edge
- Careful initialisation to start near the attention-diagonal prior (keeps gate safety)

Expected lift: **+3–8pp on accuracy metrics**. Risk: model size grows, training slows ~2–4×.

**b) Rate-aware message passing (MEDIUM)**

Current graph Laplacian treats all edges equally. On real trees, some clades evolve faster (rate heterogeneity). Use per-edge rate estimates (from tree branch lengths + observed trait variance) to scale message magnitudes. Related: Pagel λ-aware adjacency.

**c) Hierarchical / multi-scale features (LOW)**

Coarse-grained clade-level features (mean trait per deep clade) as additional graph input. Low effort; modest gain.

### 3.B Training improvements

**a) Self-supervised pretraining on large phylogenies (HIGH)**

Pretext task: on a large unlabelled tree (TreeOfLife, Open Tree), mask random trait cells and predict. Pretrain the GNN backbone, then fine-tune on the target dataset. This is what modern NLP/CV does and it works.

Existing infrastructure: `fit_pigauto` already does MAE-style training (cell masking) — we'd just run it on 100× more data without target-dataset labels.

Cost: one-time pretraining run on a cloud machine. Deliverable: `pigauto::pretrained_backbone()` that users load like `torch::resnet50()`.

**b) Knowledge distillation from BACE (MEDIUM)**

BACE (MCMCglmm-based) is slow but well-calibrated on small datasets. Train pigauto to match BACE's posterior predictions on a held-out set → pigauto gets BACE's uncertainty calibration without the compute cost. Teacher–student setup standard in DL.

**c) Multi-task learning across trait types (MEDIUM)**

Currently the GNN predicts all traits from one shared encoder, but the loss is per-trait-type (MSE / BCE / CE) without explicit multi-task weighting. Uncertainty-weighted multi-task loss (Kendall & Gal 2018) should help when trait types have very different scales/difficulties.

### 3.C Loss / Objective improvements

**a) Evidential regression for uncertainty (MEDIUM)**

Replace the Gaussian-likelihood MSE loss with evidential output (predict μ, v, α, β for a Normal-Inverse-Gamma; Amini et al. 2020). Gives predictive + epistemic uncertainty in a single forward pass. Benefits conformal and MI downstream.

**b) Direct conformal-aware training (LOW)**

Current: conformal scores computed post-hoc on val-set residuals. Alternative: train the model with a conformal-aware loss (e.g. quantile regression + pinball loss) so the conformal guarantee is baked in. Modest expected lift.

### 3.D Data augmentation

**a) Tree bootstrapping during training (MEDIUM)**

For each epoch, sample a different posterior tree from `trees300` (or a user-supplied sample). The GNN sees slightly different graphs each epoch → learns to marginalise tree uncertainty into its predictions. Free uncertainty quantification at the imputation level (the "step 1" of the Nakagawa-de Villemereuil framework).

**b) Stochastic character mapping for categorical (LOW)**

Draw ancestral states via stochastic character mapping during training. Adds noise aligned with the evolutionary model. Exotic but potentially useful for deep categorical imputation.

### 3.E Architectural alternatives worth prototyping

- **Diffusion imputation**: treat missingness as noise to be denoised. Clean fit for deep generative models.
- **State-space / structural prior models**: instead of GNN, directly parameterise a linear-gaussian SSM on the tree with learnable per-edge dynamics.
- **LLM-augmented**: for species without tree placement, embed taxonomy as text and use an LLM to place them. Science-fiction today; plausible by 2028.

### 3.F ML/DL priority triage

Do first (quick wins):
1. **Graph transformer** (3.A.a) — biggest expected lift, relatively contained change
2. **Self-supervised pretraining** (3.B.a) — expensive once, free inference forever
3. **Rate-aware message passing** (3.A.b) — natural fit for phylogenetic data

Defer until Phase 6 and Phase 8 land:
4. Evidential regression
5. BACE distillation
6. Multi-task loss weighting

Experimental / low-priority:
7. Diffusion imputation
8. Tree bootstrapping during training
9. Stochastic character mapping

---

## Part 4: A concrete 6-month roadmap

### Months 1–2: Finish Level C (current)
- Phase 6: EM for Σ — replaces Rphylopars with stable custom solver. Unlocks categorical + multi_proportion joint fits.
- Phase 7: wire Level-C baseline into `fit_pigauto` (the GNN consumer changes, not the baseline producer).
- Phase 8: discriminative benchmark suite + AVONET head-to-head vs BACE.

### Months 3–4: Architecture upgrade
- v0.9.0: Graph transformer (3.A.a). Ships with strict back-compat (old `ResidualPhyloDAE` remains available behind a flag).
- Benchmark portfolio expansion: PanTHERIA, FishBase (1.A HIGH priority).

### Months 5–6: Pretraining + rigour
- v0.10.0: Self-supervised pretraining release. `pigauto::pretrained_backbone()` fine-tuning workflow.
- Simulation framework refactor (2.B): clean sim API.
- First manuscript draft using PanTHERIA + AVONET + discriminative benchmark.

### Months 7+: Polish
- Evidential uncertainty
- TRY / plant data
- CRAN-ready release

---

## Notes on BACE parity

The near-term competitive frame is "close the AVONET OVR gap against BACE". BACE with OVR achieved:
- Trophic.Level: 72% (pigauto 65.6%)
- Primary.Lifestyle: 72% (pigauto 61.6%)
- Migration: 81.6% (pigauto 68.8%, but Phase 3/5 threshold-joint will help Migration once AVONET re-runs)
- logMass: r=0.97 (pigauto 0.89)
- logWing: r=0.97 (pigauto 0.92)

Phase 6 + Phase 4 categorical activation should close most of the categorical gap. The continuous gap (logMass r=0.97 vs 0.89) requires architecture improvements (Phase 9+, Graph transformer).

**Pragmatic target**: by Phase 6 release (v0.8.x), match BACE on categorical. By v0.9.0 (architecture), match or beat BACE on continuous.

---

## Deferred items not in any phase

These emerged in prior plans but don't fit naturally into Level C → ML/DL:

- **Richer covariate support**: currently covariates must be fully observed. Allow NAs with an explicit model (e.g. mean imputation + observed indicator). Valuable for practical users.
- **Conditional imputation**: given some trait values, predict others. Currently all traits are optional; a "set trait X to this value, predict trait Y" API would help hypothesis testing.
- **Export to BEAST / MrBayes**: multi_impute output → Nexus / XML for downstream Bayesian PCM. Ecosystem play.
