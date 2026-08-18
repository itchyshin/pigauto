# Exact Σ⊗R conditional — results

2026-08-17/18 · Claude lane · G0 granted by Shinichi. Design pre-registered at
`2026-08-17-exact-conditional-design.md` **before** the solver was written.
Branch `arc/exact-conditional`.

**Headline: this is the first change in the arc that improves the default
continuous path's accuracy** rather than improving diagnosis or honesty. It closes
50–60% of the AVONET gap to the `rphylopars` solver on three of four traits, with
no dependency.

## What was built

`predict_method = "exact"` (opt-in; default `"per_column"` unchanged) computes the
exact conditional mean/variance under `vec(L) ~ MVN(0, Σ ⊗ R)` instead of predicting
each trait independently. The trick is to condition in the **precision** form:

    precision(L_unk | L_obs) = P[unk, unk]      mean = −P[unk,unk]⁻¹ P[unk,obs] L_obs

with `P = Σ⁻¹ ⊗ S⁻¹`, so the `|O|×|O|` observed block is never inverted. `S⁻¹` is the
Hadfield & Nakagawa (2010) eq. 29 sparse precision over the extended tree — O(n)
nonzeros, already built by `build_henderson_S_inv()`. **This is the Kronecker lift of
`henderson_bm_predict()`, which has been doing exactly this for K = 1 in production.**

## Correctness — established three independent ways

1. **Dense oracle**: matches a conditional built directly from `Σ ⊗ R` to
   **7e-09 relative error** on mean *and* variance, on both the `A` and `R` scales.
   (Testing only one scale hides errors: means are scale-invariant, variances differ
   by tree height. Cost real debugging time; now a committed test.)
2. **Given the TRUE Σ it beats phylopars** (λ=0.2, n=300, 20 reps): exch07 0.639 vs
   0.645; hetero 0.870 vs 0.879. So the conditional is not merely correct, it is
   optimal — any residual gap is Σ *estimation*, not the conditional.
3. **Recovery sim**: at λ=1.0 the exact arm essentially equals phylopars
   (0.419 vs 0.413) against 0.498 for the shipped path. Coverage 0.963 / 0.934 —
   **G2 passes**, inside the pre-registered [0.93, 0.97] band.

G1 formally passes 8/12 cells; all four misses are at λ=0.2 and are attributable to
Σ estimation (see (2)), where phylopars' REML beats pigauto's closed form. **This
corrects an earlier claim of mine** that in-house Σ beat phylopars on Frobenius error
— that metric was ill-posed (different scale/parameterisation) and is retracted in
`2026-08-17-sigma-recovery-results.md`.

## Speed — why this is possible at all (K=4, 30% MCAR)

| n | trait-cells | H&N sparse | dense | speedup | nnz(S⁻¹) |
|---:|---:|---:|---:|---:|---:|
| 100 | 400 | 0.01 s | 0.02 s | 1.4× | 394 |
| 300 | 1,200 | 0.04 s | 0.35 s | 9× | 1,194 |
| 600 | 2,400 | 0.04 s | 3.26 s | **86×** | 2,394 |
| 1,200 | 4,800 | 0.10 s | infeasible | — | 4,794 |

`nnz(S⁻¹) = 4n − 6` exactly — linear, versus n² for a dense `A⁻¹`. Dense caps out near
n ≈ 600 at K = 4; fishbase (n = 10,654 × 5 traits ≈ 37,000 observed cells) is only
reachable sparsely. **Without Hadfield–Nakagawa this feature could not exist in-house.**

## G4 — the real-data gate: PASSES on 3 of 4 traits

AVONET300, 5 masks identical to `script/bench_external_comparators.md`, **paired**
per-mask deltas (the correct test — the marginal MCSEs in the raw table include
between-mask variance and overstate the noise):

| trait | per-column | exact | paired Δ | \|Δ\|/SE | gap closed |
|---|---:|---:|---:|---:|---:|
| Beak.Length_Culmen | 0.892 | **0.717** | −0.175 | **3.67** | 60% |
| Tarsus.Length | 1.170 | **1.023** | −0.148 | **2.95** | 50% |
| Wing.Length | 0.629 | **0.523** | −0.106 | **2.64** | 59% |
| Mass | 1.468 | 1.465 | −0.004 | 0.10 | 2% |

All four deltas negative; three at 2.6–3.7 SE. Reference rows on these masks:
`joint_solver = "rphylopars"` gives 0.602 / 0.873 / 0.449 / 1.295; raw phylopars
0.445 / 0.639 / 0.409 / 1.360.

**Mass is unmoved and that is not glossed**: it carries by far the largest
between-mask variance, and it is also the one trait where raw phylopars *loses* to the
rphylopars-solver configuration — so something other than the conditional governs it.
Unexplained, and flagged rather than buried.

## The void first run — recorded because it nearly became the answer

The **first** G4 run reported all-null effects and would have been reported as "the
exact conditional fails on real data." It was void: `predict_method` reached
`fit_joint_threshold_baseline()` but was never forwarded to `fit_joint_solver()`, so on
any dataset with binary/ordinal traits — i.e. all of AVONET — the option was a **silent
no-op**. Two hours of compute produced a believable false negative.

It was caught only by refusing to report a negative result without confirming the
treatment applied: tracing `exact_conditional_mvn` showed it was never called, and the
baseline was byte-identical to `per_column`. Log preserved as
`dev/bench_exact_avonet_VOID.log`.

Two lessons, both general:

- **This is the FOURTH silent-parameter-drop in `fit_baseline.R`'s call chain**
  (`lambda_mode` ×2, baseline config, covariates-into-joint, now `predict_method`).
  The Mission Control board carries a standing high-priority recruit item predicting
  exactly this. Four instances is a structural problem in how that file threads
  arguments, not a run of bad luck.
- **My end-to-end test was too simple to reach the broken branch.** It used two
  continuous traits, which take the *continuous-joint* path (which I had fixed); only
  mixed-type data exercises the threshold-joint path (which I had not). A plumbing test
  must exercise every dispatch branch the parameter can reach.

## Scope and what this does NOT claim

- Simulated evidence: matrix-normal BM exactly, K=4, MCAR, single-obs, n ≤ 300.
- Real data: ONE dataset, 5 masks, 30% MCAR, single-obs, continuous traits scored.
- **Does not reach parity** with the rphylopars solver — 50–60% of the gap, not 100%.
  The residual is Σ estimation at low signal, a separate slice.
- Default stays `per_column`. A default change should be proposed separately, on
  multi-dataset evidence, not inferred from this.
- `joint_solver = "rphylopars"` remains the opt-in yardstick and remains
  `Suggests:`-only, unreachable by default — **pigauto still has no Rphylopars
  dependency.**
