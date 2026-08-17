# Slice 1b — cross-trait refinement: what it fixes, and what it does not

2026-08-17 · Claude lane · G0 granted by Shinichi to cross the EM fence.
Branch `arc/inhouse-sigma-convergence`. **Mixed result, reported as measured.**

## The bug that was found and fixed

`.mvn_estep_refine()` combined the per-column BM posterior and the cross-trait
conditional posterior by **adding their precisions** (`1 / (prec_bm + prec_cross)`).
That is the variance of two **independent** estimates. Both are computed from the same
observed cells and are strongly dependent, so the rule double-counts the information and
understates the posterior variance. Measured consequence (recovery sim, 12 cells × 50
reps): SE coverage falls 0.925 → 0.856 (1 iteration) → 0.618 (3 iterations) while the
refined *mean* genuinely improves.

Two changes:

1. **`refine_variance = "conservative"` (new default)** — update the mean, hold the
   variance at the BM posterior value. The exact conditional variance under
   `vec(L) ~ MVN(0, Σ ⊗ R)` is what a full joint solve returns; neither rule computes it,
   and this one errs toward over-coverage, the safe direction for a package whose headline
   UQ claim is interval validity. The old rule remains as `"pooled"`, documented as
   overconfident.
2. **A divergence guard** — the Σ step must shrink each iteration; on growth or a
   non-finite step the loop rolls back to the last good iterate and reports `$diverged`.
   This is what makes `max_iter > 0` safe to expose after the 2026-05-17 disablement.

Exposed as `joint_refine_iter` on `fit_baseline()` / `fit_pigauto()` / `impute()`,
**default `0L` — no change to shipped behaviour.**

## What it fixes: moderate/low phylogenetic signal

Recovery sim, mean over cells (arms: C*n* = conservative at *n* iterations,
P*n* = pooled, E0 = current default, E2 = phylopars):

| arm | RMSE λ=0.2 | coverage λ=0.2 | RMSE λ=1.0 | coverage λ=1.0 |
|---|---:|---:|---:|---:|
| E0 (default) | 0.988 | 0.943 | 0.498 | 0.907 |
| **C1** | 0.920 | 0.958 | 0.497 | 0.907 |
| **C3** | 0.874 | 0.964 | 0.495 | 0.907 |
| **C5** | 0.863 | 0.966 | 0.502 | 0.902 |
| P1 (old rule) | 0.920 | **0.850** | 0.497 | 0.863 |
| P3 (old rule) | 0.869 | **0.670** | 0.500 | 0.838 |
| E2 phylopars | 0.812 | 0.936 | 0.413 | 0.935 |

At λ = 0.2 the fix improves **both** accuracy (−12.7% RMSE at 5 iterations) **and**
coverage (0.943 → 0.966, into the pre-registered [0.93, 0.97] band). The pooled rule is
confirmed as the coverage culprit. The guard fires on 82% of λ=1.0 reps and 0% of λ=0.2
reps — i.e. it abandons refinement exactly where refinement was known to be dangerous.

## What it does NOT fix: the AVONET gap

Real-data gate, 5 masks identical to `script/bench_external_comparators.md`, mean z-RMSE
(MCSE in parentheses):

| trait | `joint_refine_iter = 0` | `= 3` | joint_solver="rphylopars" | raw phylopars |
|---|---:|---:|---:|---:|
| Mass | 1.597 (0.510) | 1.568 (0.504) | 1.295 | 1.360 |
| Beak.Length_Culmen | 0.916 (0.165) | **0.961** (0.178) | 0.602 | 0.445 |
| Tarsus.Length | 1.199 (0.263) | 1.136 (0.273) | 0.873 | 0.639 |
| Wing.Length | 0.685 (0.085) | 0.649 (0.082) | 0.449 | 0.409 |

**Every change is well inside one MCSE, and Beak moves the wrong way.** The refinement
does not close the gap that `joint_solver = "rphylopars"` closes. **The pre-registered
real-data gate FAILS.**

### Why — measured, not assumed

AVONET's four continuous traits have Pagel's λ of **0.998, 0.993, 0.996, 0.998** —
essentially maximal phylogenetic signal. The sim's λ=1.0 block already said refinement
gives nothing there (RMSE 0.498 → 0.495). So the real-data null is the *predicted*
behaviour, not a new failure.

One prediction of mine was wrong and is corrected here: I expected the divergence guard to
be firing on AVONET. It does not (`diverged = FALSE` at 1, 3 and 5 iterations) — the loop
simply runs without converging and moves predictions negligibly. The mechanism is
"refinement has nothing to add at λ≈1", not "refinement is being aborted".

## What this means

- **Ship it as opt-in** — it is a genuine improvement in a genuine regime (moderate/low
  signal), it repairs a real statistical error in the variance, and it is guarded. Default
  stays `0L`, so nothing changes for existing users.
- **It is not the answer to the AVONET/BACE continuous gap.** That gap lives at λ≈1, where
  the difference between pigauto and phylopars is the **exact matrix-normal conditional**
  (Σ ⊗ R kriging over all cells jointly) versus per-column BM plus an approximate
  cross-trait correction. Closing it dependency-free means implementing that exact
  conditional — a real numerical project, not a parameter.
- **`joint_solver = "rphylopars"` therefore stays** as the opt-in yardstick. It remains
  a `Suggests:`-only, unreachable-by-default path (verified 2026-08-17: one call site,
  behind an explicit user choice), so pigauto still has no Rphylopars dependency.

## Regime

Sim: matrix-normal BM exactly, K=4, MCAR 30%, single-obs, n ∈ {100, 300}, λ ∈ {1.0, 0.2},
3 Σ designs, 50 reps, solver-level (no GNN). Real data: AVONET300, 5 masks, 30% MCAR,
full `impute()` path, single-obs. Coverage figures are BM analytic SE, not conformal.
