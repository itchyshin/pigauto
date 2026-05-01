# Active imputation demo: AVONET 300, 4 continuous traits

Run: 2026-04-30 19:55:13
n=300 species, mask=30%, 360 masked cells; observing 10 of them

## Test-set RMSE on REMAINING masked cells after observing 10 cells

Lower RMSE = better strategy.  Baseline = RMSE on ALL masked cells
(no observation) for reference.

| Trait | Baseline (no obs) | RANDOM | HIGH_SE | ACTIVE |
|---|---|---|---|---|
| Mass | 1.65e+03 | 2.39e+03 | 114 | 1.32e+03 |
| Beak.Length_Culmen | 20 | 20.2 | 16.5 | 18.3 |
| Tarsus.Length | 17 | 15.9 | 9 | 12.1 |
| Wing.Length | 41.4 | 50.7 | 43.7 | 39.6 |

## Interpretation

- **Mean RMSE across traits**: baseline=432, RANDOM=620, HIGH_SE=45.8, ACTIVE=347

If suggest_next_observation()'s closed-form variance-reduction
formula correlates with actual realised improvement, the ACTIVE
strategy should produce lower test RMSE than RANDOM and
(usually) than HIGH_SE.  The HIGH_SE heuristic is a sensible-
looking but suboptimal strategy: cells with high individual SE
may not have high TOTAL variance reduction (their SE is high
because they are far from observed species, but observing them
doesn't necessarily inform many other cells).

Top-3 ACTIVE recommendations (highest variance reduction):

```
Active-imputation suggestions (by = cell, n = 3)
               species              trait       type   metric    delta
 Chloephaga_rubidiceps Beak.Length_Culmen continuous variance 2.987119
 Somateria_spectabilis Beak.Length_Culmen continuous variance 2.987119
   Dendrocygna_guttata Beak.Length_Culmen continuous variance 2.473035
 delta_var_total delta_entropy_total
        2.987119                  NA
        2.987119                  NA
        2.473035                  NA

Observe the cell(s) with the largest delta_var_total to maximally reduce total imputation variance.
```
