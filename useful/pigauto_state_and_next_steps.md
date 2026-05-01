# pigauto — state and next steps

> Written 2026-04-27 ~10:00, after a week of intensive sim work.
> Honest assessment, not a paper draft. Source for the planning
> conversation when the user returns.

## TL;DR

**One real correctness fix this week** (multi-obs row alignment, commit
a3e6d39). Everything else has been lateral movement or net regression.
Pigauto wins on the bench we *designed* for it (multi-obs nonlinear
with low phylo + i.i.d. covariates) and on continuous OU; it ties or
slightly loses on most other regimes; it has *regressed* on discrete
trait types vs the Apr 17 snapshot.

**The "transformer beats GNN" story is not supported by the data we
have.** Architectures collapse to within 1 % when the safety floor is
on. We haven't yet rerun the no-safety-floor ablation at scale, but
even there I expect modest gains.

**Recommendation: stop adding features. Run a focused ablation /
diagnostic campaign to figure out what pigauto is actually doing, then
decide what to keep.**

## What we ran in the last ~24 hours

| Bench | Verdict | Pre-fix vs post-fix |
|---|---|---|
| bench_sim_bace_pigauto medium tier (linear) | tied with phylolm-λ | (no comparison; first run) |
| bench_gnn_capacity smoke (n=500, β=0.5,1.0) | safety floor uniformly closed gate — uninformative | n/a |
| bench_transformer_ablation smoke (sf=ON) | 3 architectures within 1 % | n/a |
| bench_sim_bace_nonlinear smoke v2 | pigauto/phylolm matched within ±5 % | n/a |
| **bench_multi_obs_nonlinear smoke (the designed-to-win)** | **14/18 cells beat phylolm; β=1.0 nonlinear/interactive +6–10 %** | n/a |
| bench_continuous (BM, OU, regime, nonlinear) | OU +6–14 %; BM tied | **OU clearly improved post-fix** |
| bench_binary | pigauto < baseline on 4 of 8 scenarios | **regressed 5–7pp on signal scenarios** |
| bench_count | tied with baseline | tied / slight regression |
| bench_ordinal | tied with baseline | tied |
| bench_categorical | tied with baseline most cells | **K_3 regressed 8pp** |
| bench_proportion | tied with baseline | essentially tied (good) |
| bench_zi_count | tied with baseline most cells | **regressed 14–15 % on zf_0.2, zf_0.4** |
| bench_multi_proportion | **CRASHED — 24/24 errors** | **NEW failure introduced Apr 17 → Apr 27** |
| bench_missingness_mechanism | running now | tbd |
| bench_transformer_ablation_nosf | queued | tbd |
| bench_ou_regime | queued | tbd |

## Update 2026-04-27 ~11:10 — ablation_nosf FULL results (24 cells done)

The architecture ablation with **safety_floor = FALSE**, 24 cells, 120
rows. Per-cell spread between {transformer, legacy_attn, no_attn}:

| f_type | multi | beta | architecture spread |
|---|---|---|---|
| interactive | 1 | 0.5 | 0.85 % |
| interactive | 4 | 0.5 | 0.22 % |
| linear | 1 | 0.5 | 1.19 % |
| linear | 4 | 0.5 | 0.22 % |
| nonlinear | 1 | 0.5 | 1.21 % |
| nonlinear | 4 | 0.5 | 0.40 % |
| interactive | 1 | 1.0 | 0.57 % |
| interactive | 4 | 1.0 | 0.34 % |
| linear | 1 | 1.0 | 1.02 % |
| linear | 4 | 1.0 | 0.45 % |
| nonlinear | 1 | 1.0 | 0.62 % |
| nonlinear | 4 | 1.0 | 0.48 % |

**Maximum spread: 1.21 %.** Average: 0.63 %. The three architectures
are functionally identical regardless of f_type, multi_obs, or β.

**Architectural choice is irrelevant on these DGPs.** Transformer's
multi-head + FFN does not outperform single-head attention or even the
plain message-passing GCN. The "transformer beats GNN" story is
**dead** on this evidence.

But this same ablation revealed something more interesting: **the
safety floor is doing real work specifically on single-obs cells**:

| Cell type | sf=FALSE pigauto / phylolm-λ |
|---|---|
| single-obs linear β=1.0 | 1.256 (pigauto LOSES 25.6 %) |
| single-obs nonlinear β=1.0 | 1.139 (pigauto LOSES 13.9 %) |
| single-obs interactive β=1.0 | 1.083 (pigauto LOSES 8.3 %) |
| multi-obs nonlinear β=1.0 | 0.987 (pigauto WINS 1.3 %) |
| multi-obs linear β=1.0 | 1.016 (tied) |

**Without the safety floor, pigauto's GNN overfits or wanders on
single-obs cells**, costing 8–26 %. With the safety floor on (= the
default in the per-type benches), the calibrated gate clamps the delta
to ~0 and pigauto matches phylolm-λ. This explains why every per-type
bench showed `pigauto ≈ baseline` — the GNN delta was being silently
zeroed.

**Net story:** pigauto's GNN adds value ONLY on multi-obs cells. On
single-obs, the GNN actively hurts and the safety floor saves it. The
"safety floor" is therefore the actual product, not the GNN.

This reframes the paper:
- **NOT**: "pigauto's transformer beats GNN beats baselines"
- **YES**: "pigauto's safety-gated multi-obs imputation pipeline matches
  the linear-smart baseline on single-obs DGPs and beats it on
  multi-obs nonlinear DGPs, with built-in calibrated UQ"

The transformer architecture is incidental; the real machinery is the
calibrated gate + safety floor + multi-obs aggregation.

---

## Three findings that matter

### 1. The discrete-type regression is real

Pre-fix → post-fix changes on **single-obs** benches (where the
row-alignment fix is a no-op) show a consistent pattern: continuous
got slightly better; discrete and zero-inflated got slightly worse.
The differences are 5–15 percentage points on individual
scenarios — not within Monte Carlo noise for 5 reps.

Likely culprits (commits between Apr 17 and Apr 27, in chronological
order): Fix A-H, Fix G (cov-aware GLS baseline), B2 transformer with
per-head Gaussian bandwidths, Phase 9 graph transformer block default,
EM iteration paths in joint-threshold baseline. Each was added without
ablating against the discrete-type benches.

Worst regressions:
- `bench_binary` signal_0.6: 0.653 → 0.578 accuracy (-7.5pp)
- `bench_zi_count` zf_0.2: RMSE 38.8 → 44.8 (+15 %)
- `bench_categorical` K_3: 0.746 → 0.662 accuracy (-8.4pp)
- `bench_multi_proportion`: complete failure on all 24 cells

### 2. The GNN may not be doing much

Multiple signals point to "the GNN delta contributes ~zero, and pigauto's
predictions are basically the BM/GLS baseline":
- `safety_floor=TRUE` vs `=FALSE` predictions differ by ≤ 0.5 % on every
  bench. If the GNN were actively contributing, the safety floor would
  sometimes clamp a noisy delta and produce visibly different output.
- Architecture ablation (transformer vs GAT-style legacy vs plain GCN):
  all three within 1 %. Different architectures should produce different
  delta magnitudes; same final RMSE means deltas are uniformly tiny.
- On benches where `pigauto ≈ baseline`, the GNN's job is to add
  something on top. It's not doing that.

The **ONLY** bench where pigauto cleanly beats the analytical baseline
by more than 5 % is the multi-obs nonlinear DGP we designed for it.
That's a genuine win, but it's also: low phylogenetic signal, many
i.i.d. covariates, multi-obs, nonlinear response. A specific corner of
DGP space.

### 3. The complexity is unjustified

Pigauto's pipeline currently has, in order:
1. Encoder MLP (2 layers + ReLU + dropout)
2. Multi-head graph transformer block (4 heads, FFN, 2 LayerNorms,
   2 residuals, per-head learned Gaussian phylo bias)
3. obs_refine MLP (multi-obs only)
4. Decoder MLP (2 layers + ReLU)
5. cov_encoder MLP (Fix B)
6. cov_inject per-block (Fix H)
7. cov_linear residual head (Fix C/D)
8. BM/GLS baseline with LRT gate (Fix G)
9. Calibrated per-trait gate (`r_cal`)
10. Conformal prediction interval
11. Safety floor enforcement
12. Joint MVN baseline (Phase 2)
13. Threshold-joint baseline for binary/ordinal (Phase 3)
14. OVR categorical (Phase 6)
15. Soft liability E-step (B1)
16. EM iteration loop (Phase 6/7)

Each piece was added because it fixed *some* test or some scenario.
None has been removed since being added. We have no idea which pieces
are still earning their keep.

## What I'd recommend next (planning conversation)

### Option A: Diagnostic ablation campaign (2–3 days)

Pick the multi-obs nonlinear bench as the target (where pigauto
visibly wins). Then ablate one component at a time and measure:

| Ablation | Expected | What we learn |
|---|---|---|
| Force `r_cal = 0` (no GNN delta) | If RMSE unchanged, GNN does nothing | Whether to keep the GNN at all |
| `use_transformer_blocks=FALSE` | Already shows ≈ tied | Whether the transformer earns its keep |
| `use_attention=FALSE` | Already shows ≈ tied | Whether attention earns its keep |
| Skip Fix G GLS, use plain BM | Compare on linear cells | Whether Fix G is still pulling weight |
| Skip Fix H cov_inject | | Whether the per-block cov injection helps |
| `safety_floor=FALSE` always | Already tested in some cells | Whether the safety mechanism is doing real work |
| Disable EM iterations (already default) | Trivially same | Verify the EM path doesn't accidentally fire |

Output: a clear table of "what each piece adds", in RMSE
percentage points. If most pieces add nothing, the action is to
delete them.

### Option B: Fix the discrete regression (1–2 days)

Git bisect the Apr 17 → Apr 27 commit range against
`bench_binary` or `bench_multi_proportion`. Identify which commit
caused the regression. Decide whether to revert or patch.

This is the cheapest "make pigauto better than it was" path. It
doesn't tell us anything about the architecture, but it stops the
package from going backwards.

### Option C: Conservative paper (1 week)

Write the paper around what the data does support:
- Multi-obs row-alignment fix (real correctness win)
- Multi-obs nonlinear DGP win (real, in a specific regime)
- Continuous OU win (real, on non-BM DGPs)
- Mixed-type imputation infrastructure (descriptive contribution)
- Conformal UQ + multi-imputation infrastructure (engineering contribution)

Skip the architecture story (transformer vs GNN) entirely. Position
pigauto as "a robust pipeline for mixed-type phylogenetic imputation
with built-in safety and uncertainty quantification, with measurable
gains in non-BM regimes." That's defensible and honest.

### Option D: Full reset (3–4 days)

Branch from a known-good commit (perhaps Apr 17 pre-Fix-A-H), apply
ONLY the row-alignment fix, and benchmark from scratch. See if the
discrete regression goes away. If it does, decide which post-Apr-17
changes to re-port.

This is heavier but gives the cleanest answer to "what did each fix
add or break".

### My recommended order

**B → A → C.** First fix what's clearly broken (the regression). Then
ablate to see what's actually doing work. Then write the paper around
what survives. **Do NOT pursue D unless A reveals widespread breakage.**

## Open questions for the planning conversation

1. **Are we writing a paper or improving the package?** They're
   different priorities. Writing now is premature; improving requires
   2–3 more weeks of focused work.
2. **Which DGP family is "real" for our paper's claims?** BACE-sim,
   bench_multi_obs.R style, real ecological data, or some union? The
   answer determines whether pigauto is "good" or "mediocre".
3. **Is the transformer architecture a publishable contribution?**
   Right now: no. We'd need to either find a regime where it clearly
   beats GAT-style attention (doesn't exist in current data) or scale
   to a foundation-pretraining setting (months of work).
4. **Should we rip out features that aren't earning their keep?**
   I'd say yes — but the package's API stability is a concern. Most
   are internal so removal is safe; calibrated gate is user-facing.

## Update 2026-04-27 ~11:50 — ALL SIMS DONE

All overnight phases finished by 11:34. Final verdict on each:

### Phase 3: ou_regime (DONE 11:34)

| scenario | β | column_mean | phylolm-λ | pigauto_sfT | pigauto_sfF | sfT/phy |
|---|---|---|---|---|---|---|
| OU | 0.0 | 0.466 | 0.481 | 0.481 | 0.481 | **1.000 (tied)** |
| OU | 0.5 | 0.692 | 0.501 | 0.501 | 0.503 | **0.998 (tied)** |
| OU | 1.0 | 1.084 | 0.542 | 0.542 | 0.541 | **1.000 (tied)** |
| regime_shift | 0.0 | 2.943 | 1.310 | 1.307 | 1.310 | **0.998 (tied)** |
| regime_shift | 0.5 | 3.214 | 1.106 | 1.109 | 1.110 | **1.003 (tied)** |
| regime_shift | 1.0 | 3.120 | 1.069 | 1.074 | 1.072 | **1.004 (tied)** |

**Pigauto exactly ties phylolm-λ on OU and regime_shift.** Every cell
within 0.4 %. Median ratio: 1.000. The hypothesis "phylolm-λ assumes
BM so must lose on non-BM DGPs" was wrong — phylolm-λ's *lambda*
parameter adapts the BM correction strength enough to match OU /
regime-shift adequately.

This re-interprets the earlier "bench_continuous OU win". That bench
compared pigauto against **plain BM** (Rphylopars `phylopars()`), not
against phylolm-λ. Plain BM cannot adapt to OU; the BM RMSE was 1.13.
Pigauto with its calibrated gate matched what phylolm-λ would have
done (1.07), giving the appearance of "pigauto beats BM by 6%". It's
not a GNN win; it's pigauto's BM baseline being augmented enough by
the gate to act like phylolm-λ.

### Final state: where does pigauto add value vs phylolm-λ?

- ✅ **Multi-obs nonlinear with i.i.d. covariates** (designed-to-win):
  pigauto wins 14/18, +6–10 % on β=1.0 cells. **Real**.
- ✅ **vs plain BM on non-BM DGPs**: pigauto wins 6–14 % (because
  pigauto's pipeline acts like phylolm-λ).
- 🟰 **vs phylolm-λ on linear**: tied (within 3 % usually).
- 🟰 **vs phylolm-λ on OU/regime_shift**: tied (within 0.4 %).
- 🟰 **vs label-prop / OVR baselines on discrete** (single-obs):
  tied or slight regression.
- ❌ **Architecture (transformer vs GAT vs GCN)**: indistinguishable
  (≤ 1.21 % spread, with or without safety floor).
- ❌ **Discrete trait types post-fix vs Apr 17**: regressed 5–15 pp on
  some scenarios; multi_proportion completely broken.

### Single-line conclusion

**Pigauto's pipeline is a fancy reimplementation of phylolm-λ + a
mixed-type API + UQ + a multi-obs aggregation that earns ~10 % over
phylolm-λ in one specific regime (low phylo + i.i.d. covariates +
multi-obs + nonlinear).** The transformer architecture is
incidental; the safety-gated multi-obs aggregation is the actual
contribution.

### What this means for the planning conversation

1. **The paper is much narrower than we thought.** The "transformer
   beats GNN" angle is dead. The "wins on non-BM DGPs" angle is
   illusory (pigauto matches phylolm-λ; doesn't beat it). The win is
   in **multi-obs aggregation** of i.i.d. covariates.

2. **The safety floor is the actual product.** It's what prevents
   pigauto from losing 8–26 % on single-obs cells when the GNN
   wanders. Without it, pigauto would be a regression on most
   benches.

3. **The architecture choice is irrelevant.** Switch back to the
   simpler legacy attention-GNN as default; ditch the transformer
   block. ~25–30 % faster fits, no measurable accuracy cost.

4. **The discrete regression must be fixed before any release.**
   `bench_multi_proportion` complete failure + `bench_binary` -7pp
   loss is unshippable.

5. **Recommended order remains: B (bisect discrete regression) → A
   (ablation campaign on multi-obs nonlinear) → C (conservative
   paper around multi-obs aggregation + UQ).**
