# pigauto v0.9.0 scaling curve

Run on: 2026-04-17 13:43:11
Machine: Darwin 25.4.0 (arm64), R 4.5.2
Commit:  888d02052ff24513071283abb851d4ad64ebede7
RSpectra available: TRUE

## Workload per n

- 2 continuous BM traits
- 1 binary trait (BM threshold)
- 1 categorical trait (K = 4, argmax of independent BM liabilities)
- 25% MCAR missingness
- Epochs: 300 for n<=300, 200 for n<=3000, 100 for n<=10000

## Per-stage wall time

| N | stage | status | wall (s) | R heap max (MB) | error |
|---|---|---|---|---|---|
| 100 | preprocess | ok | 0.0 | 51 |  |
| 100 | graph | ok | 0.0 | 54 |  |
| 100 | splits | ok | 0.0 | 52 |  |
| 100 | baseline | ok | 1.9 | 84 |  |
| 100 | train | ok | 110.3 | 67 |  |
| 100 | predict | ok | 0.2 | 61 |  |
| 300 | preprocess | ok | 0.1 | 65 |  |
| 300 | graph | ok | 0.0 | 82 |  |
| 300 | splits | ok | 0.0 | 65 |  |
| 300 | baseline | ok | 3.3 | 94 |  |
| 300 | train | ok | 117.3 | 70 |  |
| 300 | predict | ok | 0.2 | 66 |  |
| 1000 | preprocess | ok | 0.0 | 60 |  |
| 1000 | graph | ok | 0.5 | 137 |  |
| 1000 | splits | ok | 0.0 | 93 |  |
| 1000 | baseline | ok | 12.1 | 156 |  |
| 1000 | train | ok | 89.7 | 93 |  |
| 1000 | predict | ok | 0.3 | 86 |  |
| 3000 | preprocess | ok | 0.0 | 62 |  |
| 3000 | graph | ok | 6.4 | 749 |  |
| 3000 | splits | ok | 0.0 | 356 |  |
| 3000 | baseline | ok | 38.1 | 774 |  |
| 3000 | train | ok | 127.5 | 283 |  |
| 3000 | predict | ok | 0.9 | 273 |  |
| 5000 | preprocess | ok | 0.0 | 63 |  |
| 5000 | graph | ok | 29.6 | 1971 |  |
| 5000 | splits | ok | 0.0 | 878 |  |
| 5000 | baseline | ok | 75.7 | 2019 |  |
| 5000 | train | ok | 719.5 | 656 |  |
| 5000 | predict | ok | 34.3 | 643 |  |

## Accuracy on held-out MCAR cells

| N | trait | metric | value |
|---|---|---|---|
| 100 | cont1 | rmse | 0.1527 |
| 100 | cont2 | rmse | 0.2722 |
| 100 | bin | accuracy | 0.9545 |
| 100 | cat | accuracy | 0.9286 |
| 300 | cont1 | rmse | 0.1795 |
| 300 | cont2 | rmse | 0.1290 |
| 300 | bin | accuracy | 0.9167 |
| 300 | cat | accuracy | 0.8506 |
| 1000 | cont1 | rmse | 0.0914 |
| 1000 | cont2 | rmse | 0.0778 |
| 1000 | bin | accuracy | 0.9662 |
| 1000 | cat | accuracy | 0.9835 |
| 3000 | cont1 | rmse | 0.0411 |
| 3000 | cont2 | rmse | 0.0411 |
| 3000 | bin | accuracy | 0.9673 |
| 3000 | cat | accuracy | 0.9851 |
| 5000 | cont1 | rmse | 0.0329 |
| 5000 | cont2 | rmse | 0.0291 |
| 5000 | bin | accuracy | 0.9951 |
| 5000 | cat | accuracy | 0.9921 |
