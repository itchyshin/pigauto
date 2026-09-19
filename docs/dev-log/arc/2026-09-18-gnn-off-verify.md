# gnn-off V1 smoke — AVONET300, 30% user-level MCAR, seed 2026, 1 rep (a smoke, not an interval)

Commit: 6745787

- GNN-off wall: 8.3 s
- pure GNN-off wall: 1.4 s
- GNN-on wall (epochs = 300): 114.1 s

Production-interval coverage (conformal_lower/upper on the user-masked cells; nominal 0.95):

| arm | trait | metric | value | coverage |
|---|---|---|---:|---:|
| GNN-off | Mass | zRMSE | 2.299 | 0.83 |
| GNN-off | Beak.Length_Culmen | zRMSE | 0.527 | 0.92 |
| GNN-off | Tarsus.Length | zRMSE | 0.569 | 0.99 |
| GNN-off | Wing.Length | zRMSE | 0.664 | 0.98 |
| GNN-off | Trophic.Level | accuracy | 0.800 |  |
| GNN-off | Primary.Lifestyle | accuracy | 0.711 |  |
| GNN-off | Migration | accuracy | 0.756 |  |
| pure GNN-off | Mass | zRMSE | 2.299 | 0.83 |
| pure GNN-off | Beak.Length_Culmen | zRMSE | 0.527 | 0.92 |
| pure GNN-off | Tarsus.Length | zRMSE | 0.569 | 0.99 |
| pure GNN-off | Wing.Length | zRMSE | 0.664 | 0.98 |
| pure GNN-off | Trophic.Level | accuracy | 0.800 |  |
| pure GNN-off | Primary.Lifestyle | accuracy | 0.711 |  |
| pure GNN-off | Migration | accuracy | 0.756 |  |
| GNN-on | Mass | zRMSE | 2.244 | 0.87 |
| GNN-on | Beak.Length_Culmen | zRMSE | 0.642 | 0.92 |
| GNN-on | Tarsus.Length | zRMSE | 0.543 | 0.99 |
| GNN-on | Wing.Length | zRMSE | 0.541 | 0.98 |
| GNN-on | Trophic.Level | accuracy | 0.800 |  |
| GNN-on | Primary.Lifestyle | accuracy | 0.678 |  |
| GNN-on | Migration | accuracy | 0.756 |  |

## Suite and check

**Test suite**: FAIL 1, ERR 0, SKIP 8, PASS 2446
- Failure: `test-community-surface.R::tree benchmark and interval help retain current boundaries`
  - Failure line: 130, Expected FALSE but got TRUE in grepl check

**R CMD check**: E 0, W 0, N 0 (clean)

