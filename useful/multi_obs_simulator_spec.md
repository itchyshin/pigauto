# pigauto multi-obs simulator — design spec

> Draft 2026-04-26.
> Status: design, not yet implemented. Pending user approval before
> writing code.

## Motivation

Two existing simulators each fail to give us what we need for testing
pigauto's multi-obs path:

| Simulator | Strength | Why it's not enough |
|---|---|---|
| `bench_multi_obs.R` (CTmax-style) | Within-observation predictor (acclim_temp) is i.i.d. across obs — pigauto's natural home turf. Showed 10–19 % cov-lift in April. | Hand-rolled inside one bench script. Not exported. Only one DGP family. No way to dial predictor structure between within-obs and between-obs. |
| `BACE::sim_bace()` | Maintained, supports gaussian/binary/poisson/multinomial/threshold response and predictor types, plus interactions. | Predictors carry phylo-signal + species-RE — between-species structure dominates. pigauto multi-obs collapses against column-mean (smoke tier 25–38 % regression at β=0). We can't tell if that's a pigauto bug or a DGP-architecture mismatch from this simulator alone. |

The user requested (2026-04-26) that we add a simulator to the pigauto
package itself — "people want to simulate reliably" — that exposes the
knobs that matter for evaluating multi-obs imputation methods.

## Goal

A function `sim_pigauto()` (exact name TBD; could be
`simulate_multi_obs()` or `pigauto_sim()`) that:

1. **Cleanly separates between-species from within-species variance** for
   every predictor and the response. This is the key knob the existing
   simulators conflate.
2. Supports the full pigauto trait-type vocabulary (continuous, count,
   ordinal, binary, categorical, proportion, zi_count, multi_proportion).
3. Returns the same standardised output shape as `script/sim_bace_dgp.R`
   so all bench scripts can switch between simulators with one line.
4. Is exported from the package (`@export`) with a thorough roxygen
   docstring and at least one runnable example.
5. Has unit tests in `tests/testthat/`.

## Design

### Public signature (proposed)

```r
sim_pigauto(
  n_species,
  tree              = NULL,           # phylo, NULL means simulate birth-death
  obs_per_species   = 1L,             # int, vector of length n_species,
                                      #   or function(n_species) returning vector
  predictors        = NULL,           # list of predictor specs (see below)
  response          = NULL,           # one response spec (see below)
  miss_frac         = 0.30,           # MCAR fraction on response (default)
  miss_predictors   = 0.0,            # MCAR fraction on each predictor
  birth             = 0.8,
  death             = 0.4,
  seed              = NULL
)
```

### Predictor / response specs

Each predictor / response is a list with fields that decompose its
variance into FOUR parts:

| Component | Knob | Meaning |
|---|---|---|
| Phylogenetic between-species | `phylo_signal` ∈ [0, 1) | Variance share that is BM-correlated across the tree |
| Non-phylo between-species | `species_re_share` ∈ [0, 1] | Variance share that is species-IID (no tree structure) |
| Within-species (between-obs) | `within_species_share` ∈ [0, 1] | Variance share at the observation level |
| Residual noise | implicit (`1 - sum`) | What's left |

Constraints: `phylo_signal + species_re_share + within_species_share ≤ 1`;
the leftover is residual noise scaled by `total_sd`.

```r
list(
  type             = "continuous",   # or "binary", "categorical", "count", "ordinal", "proportion"
  name             = "x1",            # column name in output
  total_sd         = 1.0,             # marginal SD on the latent / link scale
  phylo_signal     = 0.3,             # share of variance from BM
  species_re_share = 0.2,             # share from non-phylo species RE
  within_species_share = 0.4,         # share from observation-level signal
  # residual share: 1 - sum = 0.1
  link             = NULL,            # for binary: "logit" (default); ordinal: "probit" + thresholds
  thresholds       = NULL,            # for ordinal/multinomial
  K                = NULL             # categorical count
)
```

### Response spec

Same as predictor, plus:

```r
list(
  type          = "gaussian",       # or any pigauto type
  name          = "y",
  total_sd      = 1.0,              # marginal SD on the linear-predictor scale
  phylo_signal  = 0.4,              # response's own phylo signal (independent of predictors)
  species_re_share = 0.0,
  within_species_share = 0.0,
  beta          = c(x1 = 0.5, x2 = 0.3),   # regression coefficients on predictors
  interactions  = NULL,             # optional: list("x1:x2" = 0.2)
  link          = NULL,
  thresholds    = NULL,
  K             = NULL
)
```

### How variance gets generated

For a continuous predictor `x` at species `s`, observation `i`:

```
x[s, i] = total_sd * sqrt(phylo_signal)        * z_phylo[s]
        + total_sd * sqrt(species_re_share)    * z_species_re[s]
        + total_sd * sqrt(within_species_share)* z_within[s, i]
        + total_sd * sqrt(1 - phylo_signal - species_re_share - within_species_share)
                                              * eps_residual[s, i]
```

where `z_phylo` ~ MVN(0, R) under BM, the rest are IID N(0, 1).

Then add a constant intercept and possibly link-transform.

For the **response** the same decomposition applies, plus the systematic
`X β + interactions` term. The `total_sd` of the response represents the
*residual SD after the systematic part is removed*, so users can specify
"phylo signal explains 40 % of the residual variance" while the linear
predictor adds further variance on top.

### Semantics of the four shares for ecologists

These are easy-to-interpret defaults that correspond to common scenarios:

| Predictor scenario | Suggested config |
|---|---|
| Within-observation experimental treatment (e.g., acclimation temp) | `phylo_signal = 0`, `species_re_share = 0`, `within_species_share = 0.9` |
| Species-level climate niche (e.g., elevation midpoint) | `phylo_signal = 0.5`, `species_re_share = 0.3`, `within_species_share = 0` |
| Trait with environmental + phylo-conserved structure | `phylo_signal = 0.3`, `species_re_share = 0.2`, `within_species_share = 0.4` |
| Pure phylogenetic legacy | `phylo_signal = 0.9`, `species_re_share = 0`, `within_species_share = 0` |

This is much more transparent than BACE's mixed σ² components.

### Output

Returns a list with the same shape as `script/sim_bace_dgp.R`:

```r
list(
  tree            = phylo,
  df_complete     = data.frame,        # species + traits, no NAs
  df_observed     = data.frame,        # same with masking applied
  mask            = list(per-trait logical vectors),
  response_name   = "y",
  predictor_names = c("x1", "x2"),
  response_type   = "gaussian",
  n_species       = ...,
  n_cases         = ...,
  meta            = list(... full call args ...)
)
```

So existing bench scripts (`bench_sim_bace_pigauto.R`,
`bench_multi_obs.R`) can swap between simulators by changing one line.

## What this simulator answers that the others can't

- **Diagnose pigauto's multi-obs regression directly**: sweep
  `predictor.phylo_signal` from 0 to 0.6 holding all else fixed; see if
  pigauto's deficit appears only at high values (confirming the
  multicollinearity hypothesis).
- **Reproduce real-data scenarios**: papers often have datasets with
  mostly within-obs predictors (experiments) or mostly between-obs
  (climate niche). Users can simulate either by adjusting the four
  shares.
- **Test pigauto's ceiling** at unmixable settings (one within-obs
  predictor + one between-obs predictor + one mixed). This is where the
  pure GNN should genuinely add lift over phylolm.

## Open design questions

1. **Should `obs_per_species` allow correlation with phylo / species
   traits?** E.g., common species sampled more. Real data shows this.
   Default to no correlation; expose a `weights` argument later if
   needed.
2. **Should the simulator emit observation-level covariates** (e.g.,
   sampling site, latitude) as separate columns even when not used in
   the response? Yes — return them in `df_complete` for users to
   experiment.
3. **Mixed-type response handling**. Multi-trait response (multiple `y`
   columns) is needed for testing pigauto's joint imputation. Should we
   require users to call `sim_pigauto()` once per trait and merge, or
   allow `response = list(...)` taking a list of specs sharing one tree
   and one set of species? Vote: support `response = list(spec1, spec2,
   ...)` from day one — it's the same draw-from-MVN trick as for
   predictors, and it avoids users having to manually align observation
   indices.
4. **Auto-naming of predictors / response.** If `name` is omitted,
   default to `x1`, `x2`, ..., `y`. Match BACE's convention.
5. **Should the simulator return the random-effect draws** (phylo, species
   RE, within) for downstream model checking? Yes — under
   `result$random_effects`, parallel to BACE's API. Users will want to
   verify their imputer recovers them.

## Implementation plan

1. **Phase 1 — scalar continuous response, single predictor.** Build
   `sim_pigauto()` for the simplest case (gaussian response, gaussian
   predictor, one tree). Unit tests check variance shares match the
   knobs to within 5 % at n=500 species, 10 reps.
2. **Phase 2 — multi-predictor + multi-response.** Add list support;
   draw predictors jointly from a block-diagonal MVN respecting each
   predictor's phylo.
3. **Phase 3 — non-gaussian types** (binary, ordinal, count, categorical,
   proportion, zi_count, multi_proportion). Same liability-with-
   threshold approach BACE uses.
4. **Phase 4 — interactions** (`x1:x2` etc.). Optional.
5. **Phase 5 — exports + tests + vignette**. roxygen, testthat suite.
   New vignette section in `inst/doc/pigauto_workflow_mixed.html`
   showing a multi-obs simulation example.

Phase 1 alone is enough to repeat the diagnostic that resolves the
current question. Phases 2–5 can be done after we understand whether
pigauto multi-obs has a real bug.

## Estimated effort

- Phase 1 (continuous, single predictor): ~200 lines + tests, half a day
- Phase 2–3 (multi-predictor, all types): ~500 lines + tests, two days
- Phase 4–5 (interactions, vignette, polish): ~300 lines + docs, one day
- **Total: ~3–4 days of focused work for a fully exported, tested
  simulator.**

## Decision points for the user

1. **Approve Phase 1 only** (cheap, answers the immediate diagnostic
   question) and delay Phase 2+ until needed?
2. **Approve full Phase 1–5 build** as a deliberate package feature,
   committing to the simulator as a long-term pigauto deliverable?
3. **Go a different direction** — for example, rely on `BACE::sim_bace`
   indefinitely and patch pigauto's multi-obs path instead of
   simulating differently?

My recommendation: **Phase 1 first**, evaluated. If the diagnostic shows
pigauto multi-obs is fine on the new simulator (i.e., the BACE DGP was
the issue), commit to Phase 2–5 and ship as a v1.0 feature. If pigauto
multi-obs is broken on every simulator, fix the bug first and only then
build the full simulator.
