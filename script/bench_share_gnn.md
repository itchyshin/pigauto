# share_gnn smoke benchmark

n = 100, T = 5 trees, m_per_tree = 1, reps = 2.

| rep | method | wall (s) | M |
|---|---|---|---|
| 1 | share_gnn=TRUE | 32.2 | 5 |
| 1 | share_gnn=FALSE | 140.3 | 5 |
| 2 | share_gnn=TRUE | 26.3 | 5 |
| 2 | share_gnn=FALSE | 116.3 | 5 |

Mean wall time ratio (FALSE / TRUE): 4.39x
