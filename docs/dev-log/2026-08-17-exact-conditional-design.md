# Exact matrix-normal conditional (Σ ⊗ R kriging) — design

2026-08-17 · Claude lane · **G0 granted by Shinichi** ("implement it in-house").
Pre-registered before implementation. Successor to
`2026-08-17-refinement-results.md`, which established that cell refinement cannot
close the gap at λ≈1 and that the remaining difference is the conditional itself.

## The problem, stated exactly

pigauto's joint baseline assumes

  vec(L) ~ MVN(0, Σ ⊗ R)

with L the n×K liability matrix, R = cov2cor(vcv(tree)), Σ the K×K trait covariance.
The correct baseline prediction is the conditional mean and variance of the missing
cells given the observed cells:

  E[L_M | L_O] = Σ_MO Σ_OO⁻¹ L_O          Var[L_M | L_O] = Σ_MM − Σ_MO Σ_OO⁻¹ Σ_OM

where the blocks are of the full nK×nK matrix Σ ⊗ R, and **M/O are arbitrary
per-cell index sets** — missingness is ragged, not block-structured, so the Kronecker
factorisation does not survive the conditioning.

What ships today instead: per-column BM (each trait conditioned independently, cross-trait
information discarded) plus an optional approximate refinement. Measured cost:
0.14–1.27 z-RMSE vs converged phylopars on AVONET
(`2026-08-16-continuous-gap-diagnosis.md`).

## Why this is tractable: condition in the PRECISION form

Naively, Σ_OO⁻¹ is |O|×|O| dense with |O| ≈ 0.7·nK — 37,000 for fishbase. Cubic. Dead.

But for a Gaussian with joint precision P, the conditional of the unobserved block is

  precision(L_M | L_O) = P_MM          mean(L_M | L_O) = −P_MM⁻¹ P_MO L_O

— **no inversion of the observed block at all**, just the M-block of the joint precision.
And the joint precision here is

  P = Σ⁻¹ ⊗ S⁻¹

where S⁻¹ is the Hadfield–Nakagawa sparse precision over the **extended** tree (tips +
non-root internal nodes), already built by `build_henderson_S_inv()` with O(n) nonzeros.
Σ⁻¹ is K×K dense but tiny. So P has ≈ K²·nnz(S⁻¹) nonzeros — sparse, and a sparse
Cholesky of P_MM is near-linear in n for a tree's fill-in structure.

**This is the same trick pigauto already uses for K = 1.** `henderson_bm_predict()`
(R/henderson_s_inv.R:260) does exactly this: it takes `no_q = c(missing tips, internals)`,
forms `Q[no_q, no_q]` and `Q[no_q, obs]`, solves `−Q_no_o y_o` through a sparse Cholesky
for the mean, and reads `diag(Q_no_no⁻¹)` for the variance. **The work here is the
Kronecker lift of a routine that already exists, is tested, and is in production for
single traits** — not a new numerical method.

## Algorithm

Cells are indexed (trait t, extended node i) → `(t−1)·N + i`, N = n_tips + n_internal.

1. `H <- build_henderson_S_inv(tree)` (cached; already done upstream).
2. `Sig_inv <- solve(Σ + eps·I)`; symmetrise.
3. `P <- kronecker(Sig_inv, H$Q)` (Matrix, sparse).
4. Index sets in the lifted space:
   - `obs`  = observed **tip** cells only,
   - `unk`  = missing tip cells **+ every internal-node cell, all traits**.
5. Mean: `chol <- Cholesky(forceSymmetric(P[unk, unk]))`;
   `v <- solve(chol, −P[unk, obs] %*% y_obs)`; missing-tip entries of `v` are the answer.
6. Variance: `diag(P_unk,unk⁻¹)` at the missing-tip positions.
7. **`cor_scale`**: R = D^(−1/2) A D^(−1/2), so observed values are pre-multiplied by
   `tip_sqrt_d` and predicted means divided by it; variances by its square. Identical to
   the K=1 path's convention — reuse `H$tip_sqrt_d`, do not re-derive.

**Root/mean convention.** S excludes the root, which is the improper-flat-root
parameterisation; the K=1 path documents that this assumes a zero root state and is
correct because `build_liability_matrix()` centres its inputs. The same holds here, and
the K=1 dense path stays as-is for back-compat (`bm_impute_col` estimates a GLS mean).

## The variance is the expensive part — and the honest risk

Step 6 needs the diagonal of a sparse inverse. The K=1 code solves `|M|` unit-vector
right-hand sides, which is O(|M|) solves — fine at |M| ≈ 100, **quadratic-ish and far too
slow at |M| ≈ 10,000**. Options, in the order they will be tried:

1. **Blocked RHS** — solve in chunks of ~256 columns. Simple, exact, ~|M|/256 sparse
   solves. Probably adequate to n ≈ 2,000.
2. **Takahashi selected inverse** — the standard sparse-inverse-subset recursion; exact,
   near-linear, but real work to implement correctly.
3. **Hutchinson stochastic estimator** — cheap, approximate; unacceptable for an
   interval-validity claim, so it is a last resort and would have to be labelled.

**Pre-registered scope decision:** ship (1) and gate the feature on problem size, with an
explicit, documented fallback to the current per-column path above the threshold, rather
than silently degrading. (2) is a follow-up if the size limit bites in practice.

## Verification — the gate is unusually strong here

The recovery sim (`dev/sigma_recovery_sim.R`) already carries **phylopars as arm E2**, and
phylopars computes this same conditional. So correctness is checkable directly, not merely
by "is it better than what we had":

- **G1 (correctness):** the new arm must match E2's RMSE to within ~3 × MCSE in every one
  of the 12 cells, and its Σ̂ must remain the existing estimate (this changes prediction,
  not estimation). A large residual gap to E2 means the implementation is wrong, not that
  the method is weak.
- **G2 (calibration):** SE coverage in [0.93, 0.97] — E2 achieves 0.935/0.936.
- **G3 (K=1 back-compat):** single-trait results unchanged (that path is untouched).
- **G4 (real data):** the AVONET 5-mask bench must move materially toward the
  `joint_solver = "rphylopars"` row (Mass 1.295 / Beak 0.602 / Tarsus 0.873 / Wing 0.449).
  **This is the claim that failed for refinement**, and it is the one that decides whether
  the dependency-free path is real.

Default stays opt-in until all four pass; only then is a default change proposed —
separately, with the evidence.

## Out of scope

Estimating Σ (unchanged — the sim showed the existing estimate is fine), multi-obs mode,
the GNN, and anything touching the conformal layer.
