# Smoke run before CP1 (2026-09-24, Totoro)

- Code: commit 69670d4 (`git archive`), `/home/snakagaw/pigauto_mi_posterior/69670d44f9/code`.
- Thread caps exported at launch; 9 R processes, one thread each.
- Sampler defaults: 4 chains, each 1,000 burn-in plus 5,000 sweeps.

## Simulation (`sim/`, `sim_campaign.log`)

Five cells at rep 1, run 5 at a time through `script/mi_gls/12_totoro_campaign.sh`: 5 of 5 ok.

| Cell | n | Sampler wall (s) | Cell wall (s) | max R-hat | min ESS |
|---|---|---|---|---|---|
| regime 1 (lambda 1, MCAR, x only) | 300 | 89 | about 95 | 1.003 | 1,636 |
| regime 23 (two lambdas, MAR_phylo, both) | 300 | 95 | about 103 | 1.002 | 1,241 |
| regime 3 (lambda 1, MCAR, x only) | 1000 | 189 | about 225 | 1.006 | 2,221 |
| regime 12 (lambda 0.5, MCAR, both) | 1000 | 199 | about 237 | 1.008 | 1,699 |
| regime 22 (two lambdas, MCAR, both) | 1000 | 188 | about 223 | 1.001 | 3,371 |

Cell wall covers the sampler (proper and improper from one run), the complete, posterior and improper
downstream fits (gls and phylolm on 20 datasets each), and the conformal arm (`impute(gnn = FALSE)`,
2 to 16 s). One replicate per cell cannot measure bias or coverage; the per-cell coverage in the
`.rds` files only shows that the pipeline runs.

## G4 (`g4.log`)

`CONVERGENCE_OK`: regime 1 rep 1 max R-hat 1.003, min ESS 1,636; regime 23 rep 1 max R-hat 1.002, min
ESS 1,241. Wall 3:06.

## Real data (`real/`, `real_*.log`): SHORT chains for timing only

`MI_POST_NITER=100 MI_POST_BURNIN=50`, 600 sweeps per fit. None converged (expected at this length),
so the coverage numbers are not results.

| Cell | n tips | Traits kept | Sampler wall (s) | Cell wall | Peak RSS |
|---|---|---|---|---|---|
| avonet-mcar-m20260818 | 1,500 | 4 | 50 | 1:12 | 0.63 GB |
| pantheria-mcar-m20260818 | 4,027 | 4 | 60 | 1:09 | 0.71 GB |
| fishbase-structured-m20260818 | 10,484 | 5 | 688 | 12:00 | 1.5 GB |

These include set-up time. Scaling linearly to the default 24,000 sweeps gives about 30 to 35 min
(AVONET), about 40 min (PanTHERIA) and about 7.6 h (FishBase) per cell. That is an upper bound: fixed
set-up time is part of the 600-sweep figure.
