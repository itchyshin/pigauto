# Phase 6 A/B: OVR K-fits vs LP baseline on correlated categorical data

Run on: 2026-04-16 05:30:37
Species per scenario: 150, reps: 3, rho: 0.60
Data: continuous + binary + K-class categorical, all BM with rho-correlation.

Accuracy on held-out categorical test rows (argmax):

```
 K   acc_ovr    acc_lp        lift
 3 0.6083264 0.6278195 -0.01949318
 5 0.7905876 0.6605402  0.13004734
```

OVR lifts >= 2pp: 1 / 2
OVR lifts >= 5pp: 1 / 2
