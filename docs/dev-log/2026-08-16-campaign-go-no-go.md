# ☀️ Morning go/no-go — Tier-1 Totoro campaigns (D-139 gate, pre-run COMPLETE)

2026-08-16, staged overnight. Everything below is measured, not guessed. **Reply "launch"
(or "launch regime map only") and both campaigns run back-to-back on Totoro, ~2 h total.**

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

| Campaign | Fits | CPU-h (measured walls) | Wall @ 100 workers |
|---|---|---|---|
| Regime map (108 cells; design `2026-08-15-regime-map-design.md`) | 5,400 | ≈ 79 | **≈ 1–1.5 h** |
| Coverage (24 cells; design `2026-08-15-coverage-campaign-design.md`) | 1,920 | ≈ 27 | **≈ 0.5 h** |

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
