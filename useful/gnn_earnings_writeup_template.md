# GNN earnings — does pigauto's neural architecture beat a smart linear baseline?

> **Status:** template. Numerical fields will be filled in by
> `script/analyze_gnn_earnings.R` once v2 phases finish.

## The question

The earlier real-data benches (PanTHERIA mammals, GlobTherm ectotherms, AmphiBIO amphibians, LepTraits butterflies) showed pigauto + climate / phenology covariates lifts trait imputation by 5–24 % on at least one trait per dataset. But a careful audit revealed that those lift numbers are largely consistent with what a **smart linear baseline** (phylogenetic linear regression with a fitted Pagel's λ — `phylolm(model="lambda")` — and BLUP correction) would already capture.

So the question this experiment answers is **narrower and harder**: does the GNN add value **beyond** what `phylolm-lambda` already extracts, on regimes where covariates carry nonlinear or interactive signal that a linear model cannot fit?

## Pre-registered decision rule

> pigauto + cov **earns its keep** if BOTH:
> 1. `mean(RMSE_pigauto_cov_sfT) < 0.90 · mean(RMSE_phylolm_lambda_blup)`, and
> 2. `|mean RMSE gap| > 1.96 · pooled simulation SE`.
>
> Evaluated on cells with `f ∈ {nonlinear, interactive}` and `β = 0.4`.

This rule was committed to git (`script/analyze_gnn_earnings.R`) before any v2 results were available.

## Experimental design (v2)

Semi-synthetic DGP on real phylogenies:

```
y_i = sqrt(α) · phylo_i + sqrt(β) · f(cov_i) + sqrt(1 − α − β) · ε_i

  phylo  ~ standardised BM(tree)
  cov    = 5 standard normals OR (cov_type=phylo_structured) cov1 ~ BM(tree)
  f(.)   ∈ {linear (=cov1), nonlinear (sin·exp), interactive (cov1·cov2 + 0.5·cov1²)}
  ε      ~ N(0, 1)
```

Sweep: `α ∈ {0.1, 0.4, 0.7}` (weak / mid / strong phylo signal) × `β ∈ {0, 0.2, 0.4}` (no / weak / mid covariate signal) × `f ∈ {linear, nonlinear, interactive}` × `cov_type ∈ {iid, phylo_structured}` × **5 reps**. Cells with `α + β > 0.9` (zero noise floor) are skipped.

Three trees with different structures and sizes:

- **`tree300`** — n=300 real bird phylogeny from AVONET 300
- **`mammal_tree`** — n=1000 random subsample of Bininda-Emonds 2007 (4,629-tip molecular)
- **`amphibio_tree`** — n=1500 taxonomic tree built from AmphiBIO Order/Family/Genus/Species with Grafen branch lengths

Six methods compared per cell:

| method | what it captures |
|---|---|
| `column_mean` | constant prediction (floor) |
| `lm(y ~ cov)` | covariates only, no phylogeny |
| `phylolm_lambda_blup` | linear cov + BM-with-noise (Pagel λ) BLUP |
| `pigauto_no_cov` | phylogeny only, no covariates |
| `pigauto_cov_sfT` | phylogeny + covariates, safety_floor = TRUE (default) |
| `pigauto_cov_sfF` | phylogeny + covariates, safety_floor = FALSE |

The smart linear baseline is `phylolm_lambda_blup`. **The decisive question is whether `pigauto_cov_sfT` beats it.**

## What v1 (the earlier flawed version) showed

The first round of phases used `phylolm(model="BM")` (no residual noise — misspecified for the iid-noise DGP) and i.i.d. covariates only (best case for pigauto). Those results are in `script/bench_gnn_earnings_{tree300,mammal_n1000,amphibio_n1500}.rds`.

**v1 is useful only as the upper bound for pigauto's apparent value-add** — what we got with the linear baseline handicapped. Whatever fraction of the v1 lift survives in v2 (with the misspecification fixed and phylo-redundant covariates included) is the **defensible** claim.

[v1 numbers TBD after v1 phases finish.]

## v2 results (decisive)

### Headline

[**[fraction]** of decisive cells (n=[total]) had pigauto + cov beat the smart linear baseline by the pre-registered rule.]

[Verdict: GNN earns its keep / mixed / scope claim down.]

### Per-cell verdict table

[Table from `useful/gnn_earnings_verdict.md`.]

### By covariate-structure regime

| cov_type | pigauto wins | median ratio | interpretation |
|---|---:|---:|---|
| `iid` (cov independent of phylogeny) | [n/N] | [ratio] | best case for pigauto; if it fails here, the GNN doesn't earn its keep at all |
| `phylo_structured` (cov has BM signal) | [n/N] | [ratio] | realistic comparative-biology case; this is where the paper claim must hold |

### By phylogenetic-signal regime (α)

| α | pigauto wins | interpretation |
|---|---:|---|
| 0.1 (weak phylo) | [n/N] | does pigauto over-extract noise when phylogeny is uninformative? |
| 0.4 (mid phylo) | [n/N] | the bulk of real-data regimes |
| 0.7 (strong phylo) | [n/N] | does the safety floor protect when phylogeny dominates? |

### By covariate-effect type (f)

| f | pigauto wins | interpretation |
|---|---:|---|
| `linear` | [n/N] | should be a tie — phylolm is the right model |
| `nonlinear` | [n/N] | the GNN should win here if it's doing real signal extraction |
| `interactive` | [n/N] | hardest case — interactions involve cross-feature structure |

### Heatmap

See `useful/gnn_earnings_figure.png`. Green tiles = pigauto better; red = phylolm better.

## Honest framing for the paper

Based on the v2 verdict, the paper's claim should be:

[Option A — if pigauto wins ≥50%]
> pigauto's neural architecture extracts nonlinear / interactive covariate-effect signal that linear phylogenetic methods (phylolm-lambda BLUP) cannot. On semi-synthetic data with realistic phylogenetic structure, pigauto + covariates beats the smart linear baseline by [X]–[Y] % RMSE on the regimes where covariates have nonlinear effects (`f ∈ {nonlinear, interactive}, β = 0.4`).

[Option B — if pigauto wins 25–50%]
> pigauto's neural architecture provides modest improvements over `phylolm-lambda` BLUP on specific regimes (α = [X], cov_type = [Y], f = [Z]). The architecture is most valuable when [scope conditions]. On regimes outside this scope, pigauto is comparable to or worse than the smart linear baseline; users should default to `phylolm-lambda` for those cases.

[Option C — if pigauto wins <25%]
> pigauto's primary value is **not** nonlinear signal extraction beyond `phylolm-lambda` BLUP. On semi-synthetic GNN-earnings sims, pigauto + covariates does not beat the smart linear baseline by a defensible margin. The architecture's documented benefits are: (i) the safety-floor calibration that adapts the BM/cov/mean blend per trait without manual tuning, (ii) multi-observation handling for repeated-measure ecological data, (iii) integrated uncertainty quantification (conformal intervals + multiple imputation), and (iv) a unified API across mixed trait types. These benefits are real and worth keeping; the "neural architecture extracts nonlinear signal" claim is **not supported by these simulations** and should not appear in the paper.

## Limitations

This experiment evaluates pigauto on a single-trait, additive (`phylo + f(cov)`), non-GxE DGP. Pigauto's architecture supports multi-trait sharing through the GNN, which a phylolm comparison cannot exercise (each trait is fit separately). A multi-trait extension would test:

- whether the GNN exploits cross-trait correlation that joint-MVN baselines (`Rphylopars`) miss
- whether nonlinear cross-trait interactions (e.g., `y_2 = sin(y_1) + ε`) are recovered

This would be a Phase 8 sim and is out of scope for this round.
