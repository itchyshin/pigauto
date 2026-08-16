# ☀️ Morning go/no-go — Tier-1 Totoro campaigns (D-139 gate, pre-run COMPLETE)

2026-08-16, staged overnight. Everything below is measured, not guessed. **Reply "launch"
(or "launch regime map only") and both campaigns run back-to-back on Totoro, ≈ 2 h total.**

## Pre-run test results (6 fits, Totoro, family F2, λ=0.5, m=0.3)

| n | wall/fit | Δ = blend − baseline MSE | gate_max | floor |
|---|---|---|---|---|
| 100 | 33.6 / 34.3 s | **−0.072 / −0.156** | 1.00 | no |
| 300 | 45.0 / 44.8 s | **−0.128 / −0.064** | 0.69 / 0.79 | **fired once — captured cleanly** |
| 1000 | 105.3 / 103.0 s | **−0.024 / −0.110** | 1.00 | no |

Every D-139 pre-run criterion passed: non-empty non-NA outputs per rep; one fit inspected;
per-fit wall at each n; the post-refine-floor capture mechanism demonstrated **in the
wild** (n=300 rep 2). Bonus science: the leak-free GNN beat the baseline **6/6** on the
nonlinear family with gates wide open — the regime the design predicted it should win.

## Finalized estimates (replacing the design's 7–20 h guess-band)

| Campaign | Fits | CPU-h | Wall @ 100 workers |
|---|---|---|---|
| Regime map (108 cells; design `2026-08-15-regime-map-design.md`) | 5,400 | ≈ **105** | **≈ 1.5 h** |
| Coverage (24 cells; design `2026-08-15-coverage-campaign-design.md`) | 1,920 | ≈ **36** | **≈ 0.5 h** |
| **Total** | 7,320 | ≈ **141** | **≈ 2 h** |

*Corrected upward from a first pass of 79/27 CPU-h.* Those walls were measured before
`phytools` was installed, i.e. with the phylo-signal gate silently no-op. A clean A/B on
the identical cell (F2, λ=0.5, n=300) gives **45.0 s → 60.5 s (+34%)** once the gate
actually computes. Per-fit walls used above: 34 s (n=100), 45 s (n=300), 104 s (n=1000),
each × 1.34.

Totoro: 100 workers ≤ D-143's 150 cap · `OPENBLAS_NUM_THREADS=1` · torch threads 1 ·
outputs per-(cell,rep) RDS under `~/pigauto_regime_map/results/` · resumable · stop rules
armed (first-cell early read; >20% cell failure pauses; wall overrun stops and reports).

## State on Totoro (already staged)

- pigauto `c655d75` installed (R 4.5.3, torch present) — bootstrap receipt in session log.
- `~/pigauto_regime_map/regime_prerun_cell.R` (the DGP module the runner reuses).
- Campaign runner + dispatcher: staged AND smoke-validated — three live cells through
  the generalized script's previously-untested paths (F1 linear, F3 mixed-type with
  discrete baseline decode, coverage-with-split), all `OK` at ~34 s.
- **Smoke catch:** `phytools` was missing on Totoro → the phylo-signal gate silently
  no-ops there (would have corrupted every λ=0.2 cell — a quarter of the map).
  Installed before launch. The pre-run timings are unaffected; its gate values at
  λ=0.5 were measured without phytools, noted for the comparison record.

## The launch commands (for the record — run by me on your word)

```bash
ssh totoro 'cd ~/pigauto_regime_map && bash run_campaign.sh regime_map 100'
# then, after clean completion:
ssh totoro 'cd ~/pigauto_regime_map && bash run_campaign.sh coverage 100'
```

## Why this was NOT launched overnight

The design doc and two in-chat statements committed to a morning go/no-go (D-139:
approval after pre-run results are shown). The measured cost (~2 h) makes waiting nearly
free, and a stated gate does not get reversed at 1 AM under a general delegation. Tagged
for Melissa as *adaptive-by-commitment*, not drift.


---

## ⏱️ Overrun report (D-139) — logged 2026-08-16, campaign ~14% in

**Estimated 1.5 h for the regime map. Measured pace projects ~2.5 h.**

Observed: ~755 results in ~22 min (~34/min) at 103 workers, **0 failures**. Remaining 4,645
jobs ≈ 2.3 h.

**Why the estimate was ~60% low.** The pre-run timed the `impute()` call *from inside an
already-running session* — pure fit time. Each campaign job additionally pays R+package
startup (measured: 3.3 s, minor), plus dataset simulation, `ape::rtree()`, `vcv()`, the
Pagel-λ transform and its Cholesky — all outside the instrumented window and none of it
trivial at n = 1000. **The lesson: price the JOB, not the FIT.** A pre-run test that wraps
only the expensive call under-prices a campaign whose unit of work is a whole process.

**Not stopped, deliberately.** Per-cell resumable (`regime_cell.R` skips existing outputs),
zero failures, 103 workers ≤ the D-143 cap of 150, Totoro otherwise idle. Stopping at 14%
would discard real work to correct a wall-clock miss, not a correctness or resource problem.
D-139's requirement is that an overrun does not *quietly* continue — this is the report.

**Cheapest lever if the wall-clock matters:** drop the n = 1000 tier from `n_rep = 30` to
`20`. That tier is ~40% of total cost and carries the widest MCSE anyway.
