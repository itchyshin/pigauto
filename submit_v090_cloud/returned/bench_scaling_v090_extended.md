# pigauto v0.9.0 scaling curve — extended (n = 5000-10000)

Run on: 2026-04-18 12:13:44
Machine: Linux 5.14.0-611.42.1.el9_7.x86_64 (x86_64), R 4.4.0
RSpectra available: FALSE

## Workload per n

- 2 continuous BM traits
- 1 binary trait (BM threshold)
- 1 categorical trait (K = 4, argmax of independent BM liabilities)
- 25% MCAR missingness
- Epochs: 100 for all n >= 5000

## Per-stage wall time

| N | stage | status | wall (s) | R heap max (MB) | error |
|---|---|---|---|---|---|
| 5000 | preprocess | ok | 0.0 | - |  |
| 5000 | graph | ok | 45.2 | - |  |
| 5000 | splits | ok | 0.0 | - |  |
| 5000 | baseline | ok | 2.5 | - |  |
| 5000 | train | ok | 393.1 | - |  |
| 5000 | predict | ok | 5.1 | - |  |
| 7500 | preprocess | ok | 0.0 | - |  |
| 7500 | graph | ok | 141.1 | - |  |
| 7500 | splits | ok | 0.1 | - |  |
| 7500 | baseline | ok | 6.8 | - |  |
| 7500 | train | ok | 835.1 | - |  |
| 7500 | predict | ok | 10.7 | - |  |
| 10000 | preprocess | ok | 0.0 | - |  |
| 10000 | graph | ok | 314.4 | - |  |
| 10000 | splits | ok | 0.1 | - |  |
| 10000 | baseline | ok | 14.3 | - |  |
| 10000 | train | ok | 1394.9 | - |  |
| 10000 | predict | ok | 15.8 | - |  |

## Accuracy on held-out MCAR cells

| N | trait | metric | value |
|---|---|---|---|
| 5000 | cont1 | rmse | 0.0329 |
| 5000 | cont2 | rmse | 0.0292 |
| 5000 | bin | accuracy | 0.9139 |
| 5000 | cat | accuracy | 0.9709 |
| 7500 | cont1 | rmse | 0.0265 |
| 7500 | cont2 | rmse | 0.0258 |
| 7500 | bin | accuracy | 0.7271 |
| 7500 | cat | accuracy | 0.8546 |
| 10000 | cont1 | rmse | 0.0225 |
| 10000 | cont2 | rmse | 0.0228 |
| 10000 | bin | accuracy | 0.6154 |
| 10000 | cat | accuracy | 0.6396 |
