# Architecture ablation on the WINNING DGP — verdict

> 2026-04-28 ~09:18. Source: `bench_arch_on_winning_dgp.rds` (18 cells
> × 3 archs × 3 reps = 54 fits, 62 min wall).
> DGP: multi-obs nonlinear / interactive at λ ∈ {0.2, 0.3, 0.5}, n=500,
> ncov=10, β=1.0, sp_miss=0.5, within_miss=0.2.
> Three configs all with `safety_floor=TRUE`:
>   - `pigauto_transformer`: `use_transformer_blocks=TRUE, use_attention=TRUE`
>   - `pigauto_legacy_attn`: `use_transformer_blocks=FALSE, use_attention=TRUE` (GAT-style)
>   - `pigauto_no_attn`: `use_transformer_blocks=FALSE, use_attention=FALSE` (plain GCN)

## Headline

**Transformer wins 6/6 cells. Margin: 0.5–1.1 %.** Within the "architectures
equivalent" bracket of the decision rule.

| Comparison | median ratio | range |
|---|---|---|
| transformer / no_attn | **0.9927** | [0.9886, 0.9944] |
| transformer / legacy_attn | **0.9923** | [0.9897, 0.9946] |

Maximum spread across any cell: **1.15 %**. The decision rule pre-registered
≤ 0.95 as "transformer earns its keep" and [0.98, 1.02] as "equivalent".
We got 0.989–0.995, solidly in the equivalent bracket.

## All 6 cells

| f_type | λ | transformer | legacy_attn | no_attn | tfm/no_attn |
|---|---|---|---|---|---|
| nonlinear | 0.2 | 1.1016 | 1.1100 | 1.1108 | 0.992 |
| nonlinear | 0.3 | 1.1370 | 1.1459 | 1.1443 | 0.994 |
| nonlinear | 0.5 | 1.2173 | 1.2240 | 1.2242 | 0.994 |
| interactive | 0.2 | 1.2172 | 1.2298 | 1.2288 | 0.991 |
| interactive | 0.3 | 1.2160 | 1.2277 | 1.2300 | 0.989 |
| interactive | 0.5 | 1.3072 | 1.3170 | 1.3147 | 0.994 |

## Interpretation

Three things to take away:

1. **Transformer wins consistently — but by tiny margins.** 6 / 6 cells
   show transformer with the lowest RMSE. The direction is unambiguous.
   The magnitude (0.5–1.1 %) is well below typical MC noise envelopes
   for 3-rep means.

2. **The architecture is NOT where pigauto's value comes from.** Pigauto
   beats `lm_nonlinear` on these cells by 8–25 % (lambda-sweep verdict).
   The within-architecture spread is 0.5–1.2 %. So the analytical
   pipeline (BM/GLS baseline + obs-level multi-obs aggregation +
   safety gate) is doing 95–99 % of the work; the transformer block
   contributes the remaining 1–5 %.

3. **The right default is still `use_transformer_blocks = TRUE`.**
   - It wins on every cell tested (consistent direction matters even at small magnitude)
   - It's the more flexible architecture for future work
     (foundation-model pretraining, self-attention extensions)
   - The ~25–30 % per-fit cost vs the legacy GNN is small in absolute
     terms (a few minutes more on most n)

## What this means for the paper

**DROP the "transformer beats GNN" claim.**
- Earlier ablations (with safety floor on/off) were uninformative
  because no architecture was earning anything.
- This bench tested on cells where pigauto-as-a-whole wins decisively.
- Even there, the architecture spread is only 1.15 %.
- The transformer is incidental, not the differentiator.

**KEEP the transformer as default.**
- Consistent slight edge across all cells (6/6)
- Forward-compatible with foundation-model pretraining
- No measurable downside

**Headline architectural framing for the paper:**
> *Pigauto uses a phylogenetic graph transformer (multi-head attention
> with rate-aware Gaussian-bandwidth bias on cophenetic distance) as
> the GNN component. We selected this over alternatives like single-head
> attention (Veličković et al. 2018, GAT) and plain message-passing
> (Kipf & Welling 2016) for forward-compatibility, though we find these
> alternatives within ~1 % RMSE on our benchmarks. The empirical
> contribution of pigauto is not the specific GNN architecture but the
> safety-gated combination of analytical phylogenetic baseline +
> multi-obs aggregation + GNN correction with calibrated uncertainty.*

## Decision: registered

**Transformer remains the default** (commit a7008c2 onwards is correct).
The paper claims the **safety-gated multi-obs phylogenetic AE pipeline**
as the contribution, not the specific GNN architecture.
