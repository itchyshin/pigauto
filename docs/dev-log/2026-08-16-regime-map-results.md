# Regime map (#135) — results

Campaign complete 2026-08-16 ~07:45 MDT on Totoro. **5,400 jobs, 21,620 trait-rows,
0 failures.** Design and pre-registered interpretation:
`2026-08-15-regime-map-design.md` (ADEMP; Morris 2019 + Williams 2024).
Pipeline: `main` @ `c655d75` (post-#158/#159, leak-free training + two-layer #157 floor).

**Δ = per-rep paired (blend − baseline) test loss on the held-out surface. Negative = the
GNN helps.** MCSE = `sd(Δ)/√n_rep` on the paired difference. The design's detection rule is
`|Δ|/MCSE ≥ 3`.

## The headline

| family | λ | mean Δ | median \|Δ\|/MCSE | P(gate open) | P(floor fired) |
|---|---:|---:|---:|---:|---:|
| F1 linear | 0.2 | −0.104 | 8.94 | 0.54 | 0.046 |
| F2 nonlinear | 0.2 | −0.189 | 9.13 | 0.47 | 0.028 |
| F3 mixed | 0.2 | −0.273 | 9.86 | 0.24 | 0.033 |
| F1 | 0.5 | −0.029 | 4.36 | 0.87 | 0.033 |
| F2 | 0.5 | −0.086 | 5.40 | 0.89 | 0.035 |
| F3 | 0.5 | −0.152 | 5.90 | 0.67 | 0.035 |
| F1 | 0.8 | −0.006 | 2.80 | 0.87 | 0.002 |
| F2 | 0.8 | −0.028 | 4.17 | 0.95 | 0.007 |
| F3 | 0.8 | −0.060 | 2.82 | 0.62 | 0.043 |
| **F1** | **1.0** | **−0.001** | **1.41** | 0.86 | 0.002 |
| **F2** | **1.0** | **−0.017** | **3.43** | 0.96 | 0.000 |
| **F3** | **1.0** | **−0.036** | **1.16** | 0.58 | 0.013 |

## What it says — and the two readings must not be merged

### 1. The specificity control passes, and the claimed niche is real

λ = 1.0 is the one row-block where the **BM baseline is correctly specified** (the DGP is BM
on an untransformed tree, which is exactly what the baseline assumes). That is where the
design said the answer would actually discriminate. Detection rates there:

| family | mean \|Δ\|/MCSE | cells ≥ 3 MCSE |
|---|---:|---:|
| **F1 linear (specificity control)** | 1.51 | **4 / 36 (11%)** |
| **F2 nonlinear (claimed niche)** | 3.29 | **21 / 36 (58%)** |
| F3 mixed-type | 1.72 | 7 / 36 (19%) |

This is the pre-registered discrimination, and it came out the right way round:

- **F1 is null.** Where the baseline is correct and the DGP is linear-Gaussian, the GNN adds
  nothing detectable — several cells even go slightly positive. The design stated in advance
  that *a GNN win on F1 would be a red flag, not a success*. It did not fire.
- **F2 is not null.** The nonlinear cross-trait links (`sin(2·Z₁)`, `Z₂²·sign(Z₂)`) that a
  linear joint-MVN baseline structurally cannot represent are recovered by the GNN, in the
  majority of cells, with the baseline at its best.

That pairing is the evidence #135 was opened to get: not "the GNN is better" but **"the GNN
recovers structure the baseline cannot represent, and adds nothing when the baseline can."**

### 2. The large low-λ gains are baseline misspecification, NOT a second win

Δ grows monotonically as λ falls — F1 goes from −0.001 (λ=1) to −0.104 (λ=0.2), a hundredfold.
**This is not the GNN being cleverer at low signal.** pigauto's BM baseline assumes λ = 1. The
DGP at λ = 0.2 has far less phylogenetic signal than the baseline presumes, so the baseline is
misspecified, and the GNN is absorbing that misspecification. The effect appears in F1 —
the *linear* family — which is what identifies it as a baseline-model problem rather than a
nonlinearity the GNN uniquely captures.

**Reported honestly, that is a robustness claim, not an accuracy claim:** *the gated GNN
degrades gracefully when the phylogenetic model is wrong.* It is genuinely useful — real
trees rarely satisfy λ = 1 — but it must not be sold as the GNN outperforming a
correctly-specified baseline.

**Unrun comparison that bears directly on it:** pigauto already ships a Pagel-λ BM baseline
(NEWS, v0.10.0.9000, "New (opt-in): Pagel's lambda Brownian-motion baseline"). Some — possibly
most — of the low-λ gain may be obtainable by fitting λ in the *baseline* instead of routing
around it with a neural correction. **This campaign did not test that**, and no claim about
the low-λ regime should ship until it does. That is the single highest-value follow-up here.

### 3. Gate and floor behaviour (S1 / S2)

- **P(gate open) is high**: 0.47–0.96 across cells, and *rises* with λ (0.54 → 0.86 for F1).
  The gate opens most where the baseline is best — which is consistent with it opening on
  small, genuine refinements rather than on noise.
- **P(floor fired) is low**: 0.000–0.046, and lowest at high λ. The two-layer #157 floor
  shipped in PR #158 is therefore a rarely-exercised safety net, not a load-bearing crutch —
  the gate mostly earns its own opening. It is not decorative either: it fires on 3–5% of
  fits at λ ≤ 0.5.
- **F3's gate opens least** (0.24–0.67) despite F3 having the largest raw Δ. Discrete
  calibration is more conservative, as designed.

## Regime, and what this does not establish

- Simulated DGPs only (BM on λ-transformed trees, K = 4, Σ exchangeable r = 0.7); MCAR
  missingness only; single-obs; 60 reps (n = 100/300) and 30 (n = 1000); one machine.
- MAR/MNAR is a separate axis (`bench_missingness_mechanism.R`), untouched here.
- Conformal coverage was recorded opportunistically but is **not** designed for here — that is
  the coverage campaign (`2026-08-15-coverage-campaign-design.md`).
- No real-data anchor in this campaign. AVONET300 re-run post-#158 remains a Tier-3 item.
- The λ-misspecification finding above means **the absolute Δ values at λ < 1 are a property
  of the baseline's λ = 1 assumption**, not portable statements about the GNN.

## What this licenses for the paper (Shinichi's decision, not this lane's)

The roadmap (Tier 1.2) said the framing choice waits for this result and must not be made
retroactively by it. On the evidence:

- The **safety-gated architecture** framing is fully supported: the gate opens often, the
  floor almost never has to intervene, and the control family is null.
- A **regime-mapped lift** claim is supported *specifically for nonlinear cross-trait
  structure at correctly-specified λ* — F2 at λ = 1, 58% of cells.
- A general "GNN beats the phylogenetic baseline" claim is **not** supported, and the low-λ
  results must not be used to imply it until the Pagel-λ comparison is run.

## Artifacts

- `~/pigauto_regime_map/results/regime_map/` on Totoro — 5,400 per-cell RDS (+5 smoke).
- `~/pigauto_regime_map/regime_map_summary.{md,rds}` — full 108-cell × trait table with MCSEs.
- Runner `regime_cell.R`, dispatcher `run_campaign.sh`, summariser `summarise_regime_map.R`.
- Cost: ~2.5 h wall at ~100 workers (D-139 overrun logged in
  `2026-08-16-campaign-go-no-go.md`; the pre-run priced the fit, not the job).
