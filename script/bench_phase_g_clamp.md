# Phase G acceptance bench: clamp_outliers on AVONET Mass

Run: 2026-05-01 10:18:31
AVONET n=1500, miss_frac=0.30, seeds=2030,2031,2032, N_IMP=20, clamp_factor=5

## Pre-registered acceptance criteria

1. Seed-2030 Mass RMSE with clamp ON drops below 5,000 (currently 23,723).
2. Seeds 2031 and 2032 Mass RMSE change with clamp ON < 5 % of OFF RMSE.

## Per-seed results

```
 seed clamp       rmse  baseline  top1_pred
 2030   OFF 24330.0188 1819.6687 551050.271
 2030    ON  6273.0888 1819.6687 167846.500
 2031   OFF   312.2118  676.3354   3683.761
 2031    ON   317.6400  676.3354   3654.301
 2032   OFF   591.2584 5268.5304 100257.146
 2032    ON   432.5001 5268.5304 117486.928
```

Seed-2030 RMSE drop: 24330 -> 6273 (74.2% reduction)
Seed-2031 RMSE Δ: 312 -> 318 (1.74%)
Seed-2032 RMSE Δ: 591 -> 433 (26.85%)

## Verdict

**PHASE G FAILS criterion 1**: seed-2030 RMSE = 6273 >= 5000

