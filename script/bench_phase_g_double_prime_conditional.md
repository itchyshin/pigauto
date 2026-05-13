# Phase G'' acceptance bench: conditional PMM (pmm_when = 'outside_observed')

Run: 2026-05-01 13:00:54
AVONET n=1500, miss=0.30, seeds=2030,2031,2032, N_IMP=20.  Synthetic heavy-tail: 3 reps.
Reference numbers for `none`, `clamp`, `pmm_always` come from PR #61's bench (script/bench_phase_g_prime_pmm.rds).

## Pre-registered acceptance criteria

1. pmm_outside seed-2030 Mass RMSE <= pmm_always seed-2030 RMSE
   (don't lose the Phase G' win).
2. pmm_outside doesn't regress > 5 % vs `none` on any cell
   (the criterion that pmm_always failed).

## Per-cell results across all configs

```
              dataset              trait seed      config         rmse
               AVONET Beak.Length_Culmen 2030       clamp    13.143423
               AVONET Beak.Length_Culmen 2030        none    13.458349
               AVONET Beak.Length_Culmen 2030  pmm_always    13.159920
               AVONET Beak.Length_Culmen 2030 pmm_outside    13.323408
               AVONET Beak.Length_Culmen 2031       clamp     9.128660
               AVONET Beak.Length_Culmen 2031        none     9.540088
               AVONET Beak.Length_Culmen 2031  pmm_always     9.507171
               AVONET Beak.Length_Culmen 2031 pmm_outside     9.190915
               AVONET Beak.Length_Culmen 2032       clamp    16.149460
               AVONET Beak.Length_Culmen 2032        none    15.681230
               AVONET Beak.Length_Culmen 2032  pmm_always    17.168591
               AVONET Beak.Length_Culmen 2032 pmm_outside    16.414227
               AVONET               Mass 2030       clamp  6276.758365
               AVONET               Mass 2030        none 29930.444672
               AVONET               Mass 2030  pmm_always  1274.832317
               AVONET               Mass 2030 pmm_outside  1272.467805
               AVONET               Mass 2031       clamp   303.280870
               AVONET               Mass 2031        none   308.802783
               AVONET               Mass 2031  pmm_always   259.235775
               AVONET               Mass 2031 pmm_outside   326.920616
               AVONET               Mass 2032       clamp   660.178865
               AVONET               Mass 2032        none   373.548031
               AVONET               Mass 2032  pmm_always  3386.361922
               AVONET               Mass 2032 pmm_outside  4965.576724
               AVONET      Tarsus.Length 2030       clamp     9.685146
               AVONET      Tarsus.Length 2030        none     8.776682
               AVONET      Tarsus.Length 2030  pmm_always    10.140568
               AVONET      Tarsus.Length 2030 pmm_outside     9.193638
               AVONET      Tarsus.Length 2031       clamp    21.538876
               AVONET      Tarsus.Length 2031        none    21.344997
               AVONET      Tarsus.Length 2031  pmm_always    22.764871
               AVONET      Tarsus.Length 2031 pmm_outside    21.823332
               AVONET      Tarsus.Length 2032       clamp     9.692533
               AVONET      Tarsus.Length 2032        none    10.173029
               AVONET      Tarsus.Length 2032  pmm_always    10.108823
               AVONET      Tarsus.Length 2032 pmm_outside    10.024632
               AVONET        Wing.Length 2030       clamp    30.101613
               AVONET        Wing.Length 2030        none    31.793467
               AVONET        Wing.Length 2030  pmm_always    30.827279
               AVONET        Wing.Length 2030 pmm_outside    29.739645
               AVONET        Wing.Length 2031       clamp    27.965588
               AVONET        Wing.Length 2031        none    28.069174
               AVONET        Wing.Length 2031  pmm_always    28.941890
               AVONET        Wing.Length 2031 pmm_outside    28.748704
               AVONET        Wing.Length 2032       clamp    28.616661
               AVONET        Wing.Length 2032        none    32.651438
               AVONET        Wing.Length 2032  pmm_always    31.957375
               AVONET        Wing.Length 2032 pmm_outside    29.791265
 synthetic_heavy_tail               mass    1       clamp  2972.871089
 synthetic_heavy_tail               mass    1        none  2972.871089
 synthetic_heavy_tail               mass    1  pmm_always  3095.271685
 synthetic_heavy_tail               mass    1 pmm_outside  2972.871089
 synthetic_heavy_tail               mass    2       clamp  1627.918522
 synthetic_heavy_tail               mass    2        none  1627.918522
 synthetic_heavy_tail               mass    2  pmm_always  1775.666319
 synthetic_heavy_tail               mass    2 pmm_outside  1627.918522
 synthetic_heavy_tail               mass    3       clamp  3159.778865
 synthetic_heavy_tail               mass    3        none  3159.778865
 synthetic_heavy_tail               mass    3  pmm_always  3219.737073
 synthetic_heavy_tail               mass    3 pmm_outside  3159.778865
 pct_in_observed_range    max_pred   obs_max
              99.77778    177.4366    414.00
              99.77778    188.5971    414.00
             100.00000    148.4000    414.00
             100.00000    172.2861    414.00
             100.00000    140.8917    319.50
             100.00000    138.5133    319.50
             100.00000    142.8000    319.50
             100.00000    158.2250    319.50
             100.00000    191.7014    355.50
             100.00000    240.2364    355.50
             100.00000    159.1000    355.50
             100.00000    214.1997    355.50
              99.55556 167846.5000  35000.00
              99.55556 669873.9424  35000.00
             100.00000   9320.6000  35000.00
             100.00000   9512.1000  35000.00
              99.77778   4567.1981  34093.30
              99.55556   4576.5051  34093.30
             100.00000   4874.0000  34093.30
             100.00000   5375.5055  34093.30
              99.55556  98371.7587 111000.00
              99.55556 107393.2979 111000.00
             100.00000  39500.0000 111000.00
             100.00000   7996.8319 111000.00
             100.00000    146.8418    269.90
             100.00000    153.7754    269.90
             100.00000    134.2000    269.90
             100.00000    140.9508    269.90
             100.00000    266.9236    382.80
             100.00000    249.6961    382.80
             100.00000    126.9000    382.80
             100.00000    263.1365    382.80
             100.00000    324.0102    475.00
             100.00000    327.9506    475.00
             100.00000    267.1000    475.00
             100.00000    299.3838    475.00
              99.55556    801.0297    759.80
              99.55556    802.1483    759.80
             100.00000    615.5000    759.80
              99.55556    759.8000    759.80
             100.00000    588.7479    760.60
             100.00000    597.7704    760.60
             100.00000    566.6000    760.60
             100.00000    569.6557    760.60
             100.00000    573.9100    636.50
              99.77778    663.6921    636.50
             100.00000    554.5000    636.50
             100.00000    587.5639    636.50
             100.00000   1215.6899  16204.05
             100.00000   1215.6899  16204.05
             100.00000   2015.9158  16204.05
             100.00000   1215.6899  16204.05
             100.00000   1318.7666  29290.88
             100.00000   1318.7666  29290.88
             100.00000   1094.6236  29290.88
             100.00000   1318.7666  29290.88
             100.00000   1101.7917  22934.01
             100.00000   1101.7917  22934.01
             100.00000   1009.2288  22934.01
             100.00000   1101.7917  22934.01
```

Crit 1 (pmm_outside seed-2030 Mass <= pmm_always): PASS
Crit 2 (pmm_outside no-regress vs none, threshold 5%): FAIL on 2 cells (worst 1229.3%)

## Verdict

**PHASE G'' PARTIAL** -- FAIL 2: pmm_outside regresses > 5 % on 2 cells (worst delta = 1229%)

