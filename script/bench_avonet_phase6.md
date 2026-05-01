# AVONET 300 post-Phase-6 benchmark

Run on: 2026-04-30 08:05:03
Level-C baseline (Phase 6 active) vs LP baseline.

RMSE: lift = lp - levelC (positive = levelC better).
Accuracy: lift = levelC - lp (positive = levelC better).

```
              trait        type   metric    levelC        lp      lift
               Mass  continuous     RMSE 0.2511250 0.6761676 0.4250426
 Beak.Length_Culmen  continuous     RMSE 0.3310674 0.6110457 0.2799783
      Tarsus.Length  continuous     RMSE 0.3376763 0.6727999 0.3351236
        Wing.Length  continuous     RMSE 0.2427396 0.5697892 0.3270496
          Migration     ordinal     RMSE 0.8898098 0.8898098 0.0000000
      Trophic.Level categorical accuracy 0.7704918 0.4918033 0.2786885
  Primary.Lifestyle categorical accuracy 0.8431373 0.5686275 0.2745098
```
