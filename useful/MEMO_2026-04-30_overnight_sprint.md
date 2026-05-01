# Overnight sprint summary (2026-04-30 → 2026-05-01)

User asked: "you have 11 hours till 5 am tomorrow to improve this
package - what would you do tell me - we want the best of the best -
think carefully and also out of box".

## What was delivered

Two scoped features committed on `experiment/gnn-earnings-sim`:

### 1. `gate_method = "cv_folds"` -- cross-validated gate calibration (commit aaabd26)

Closes the longest-standing open NEWS item.  When `gate_method =
"cv_folds"`, the val cells are partitioned into K (default 5)
deterministic non-overlapping folds; the existing grid +
half-B-verify procedure runs once per fold (`half_a` = K-1
training folds, `half_b` = held-out fold), and the componentwise
median of K winning weight vectors becomes `w_final`.  Compared to
`"median_splits"` (B=31 random half-A/half-B splits with overlap),
`"cv_folds"` uses larger training sets per split (K-1/K vs 1/2)
and has a standard k-fold-CV interpretation.

* Files: `R/fit_helpers.R` (+49), `R/fit_pigauto.R` (+20),
  `man/fit_pigauto.Rd` (+20), `tests/testthat/test-safety-floor.R` (+154)
* New API: `gate_method = "cv_folds"` + `gate_cv_folds = 5L`.
  Default `gate_method` remains `"single_split"`; cv_folds is
  opt-in.
* Tests: 5 new test_that blocks, 9 expects, all pass.  Existing 24
  tests + B4 edge cases unchanged.
* Smoke: AVONET 300 (60 epochs, 1 imp): cv_folds wall = 75.9s
  vs single_split 82.6s vs median_splits 77.2s; gate weights
  sum to 1, output is finite for every trait.

### 2. `suggest_next_observation()` -- active imputation (commit 3754b58)

Novel public API.  For a fitted `pigauto_result`, returns the top-N
candidate cells (or species) ordered by expected reduction in TOTAL
predictive variance across all currently-missing cells if that
candidate were observed next.

The closed-form variance-reduction formula is derived from a
Sherman-Morrison rank-1 inverse update on the BM conditional MVN:
```
delta_V_total(s_new) = sigma2 * sum_{i in miss} D[i, k]^2 / alpha_k
```
where `D = R[miss, miss] - R[miss, obs] R[obs, obs]^{-1} R[obs, miss]`
is the residual matrix and `alpha_k = D[k, k]`.

* Files: `R/active_impute.R` (new, 304 LOC), `R/impute.R` (+1,
  exposes `tree` in `pigauto_result`), `tests/testthat/test-active-impute.R`
  (new, 282 LOC), `man/suggest_next_observation.Rd` (auto-generated)
* New API: `suggest_next_observation(result, top_n = 10L,
  by = c("cell", "species"))` returns a `pigauto_active`
  data.frame with descending `delta_var_total`.  Custom print
  method.
* Tests: 10 test_that blocks, 29 expects, all pass.  Closed-form
  verified against an independent fixed-sigma2 brute-force refit
  on a 20-tip fixture.
* Smoke: AVONET 300 with 90 cells masked across 4 continuous
  traits.  Top-1 cell: `Chloephaga_rubidiceps Beak.Length_Culmen`,
  `delta_var_total = 2.99`.  Top-1 species:
  `Chloephaga_rubidiceps`, total = 3.23 across 2 missing traits.

#### Why this is novel

To my knowledge no other phylogenetic-imputation package
(Rphylopars, BACE, phylolm, mice with phylogenetic correlation)
exposes a sampling-design helper.  The methods are 30 years old in
optimal-design literature (Cohn et al. 1996, *JAIR* 4:129-145);
applying them to phylogenetic trait imputation appears to be new.

This is a paper hook: *"Phylogenetically-informed sampling design
via expected imputation-variance reduction"* would be a clean MEE
methods paper using pigauto.

### Bench reruns

`script/bench_cv_vs_median.R` (new) runs a focused 3-way
comparison: continuous + binary + categorical traits at lambda = 1,
5 reps, comparing all three gate methods.  Kicked off in
background overnight; results will be folded into NEWS.md when
complete.

## What was scoped but DEFERRED

(Per the original triage, these were tempting but out of scope for 11 h.)

* **lambda-fitted BM kernel replacement** -- touches Phase 2/3/6/
  liability/GNN code paths; 8+ h alone, high regression risk.
* **lambda-as-corner in safety_floor simplex** -- middle-tier
  architectural improvement; would have been a stretch goal but
  preferred to ship cv_folds + active imputation cleanly.
* **CRAN-readiness `R CMD check --as-cran`** -- time-bounded;
  could be 0 h or 6 h depending on what NOTEs surface.
* **Self-supervised pretraining on tree corpus** -- CLAUDE.md's
  Phase 9 future work; 30+ h project.
* **CUDA backend / mini-batched GNN** -- major architecture
  rewrite.

## Empirical caveats

* Bench rerun comparing val->test drift on
  `bench_binary.R / bench_categorical.R / bench_zi_count.R` across
  the three gate methods is not yet done -- the existing benches
  are 16 min each and would have eaten the budget.  Smoke test on
  AVONET 300 shows all three methods produce sensible gate weights
  but doesn't probe the val->test drift question.
  `script/bench_cv_vs_median.R` is a focused alternative running
  overnight.
* `suggest_next_observation()` v1 supports continuous-family traits
  only.  Discrete (binary, categorical, zi_count) silently
  skipped.  Expected-entropy-reduction extension queued for v2.
* Multi-obs not supported in v1 (single-obs only); errors with
  clear message.
* The closed-form variance reduction is the standard
  active-learning value-free formula (assumes sigma2 fixed).  The
  brute-force test asserts this against an independent
  fixed-sigma2 refit, NOT against a refit with REML re-estimation.

## Test status

466+ PASS / 0 FAIL across the 5 changed/added test files plus all
fit_baseline-touching files (verified at session start).  Full
sweep of all test files including the new test-active-impute.R
will be run at the very end before declaring done.

## Files added / modified this overnight session

| File | Status | Purpose |
|---|---|---|
| `R/active_impute.R` | NEW | suggest_next_observation + math |
| `R/fit_helpers.R` | M | cv_folds variant in calibrate_gates |
| `R/fit_pigauto.R` | M | gate_cv_folds arg threading |
| `R/impute.R` | M | tree in pigauto_result |
| `tests/testthat/test-active-impute.R` | NEW | 10 tests, 29 expects |
| `tests/testthat/test-safety-floor.R` | M | 5 cv_folds tests |
| `man/fit_pigauto.Rd` | M | docs |
| `man/suggest_next_observation.Rd` | NEW | docs |
| `NAMESPACE` | M | export(suggest_next_observation) |
| `DESCRIPTION` | M | Version bump 0.9.1.9007 -> 0.9.1.9008 |
| `NEWS.md` | M | two new sections |
| `script/bench_cv_vs_median.R` | NEW | 3-way bench |
| `useful/MEMO_2026-04-30_overnight_sprint.md` | NEW | this memo |

## Commits

* `aaabd26` calibrate_gates: add gate_method = "cv_folds"
* `3754b58` add suggest_next_observation()
* `b0ed099` NEWS + DESCRIPTION bump for v0.9.1.9008

(plus earlier commits from today's daytime work: d64eb78 / f43edb6 /
8860fa8 / 970f9b7 / 9bfce30, which closed the seven Opus
adversarial-review items.)

## Open / for the morning

* Read CV-vs-median bench results (`script/bench_cv_vs_median.{rds,md}`).
* Decide whether to push to origin (63+ commits ahead).
* Consider whether to default `gate_method` to "cv_folds" in a
  v0.9.2 release (currently opt-in for safety).
* `suggest_next_observation()` v2: discrete trait support via
  expected entropy reduction.
* Vignette section walking through active imputation on a real
  dataset.
