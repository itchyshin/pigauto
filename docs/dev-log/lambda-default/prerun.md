# Pre-run: one core cell under the lambda default

Cell: types_mixed, BM, lambda 0.3, rho 0, MCAR 30%, n 300. Seeds 1 to 20. Arms: gnn_off, gnn_off_rphylopars,
gnn_on. Totoro, 20 parallel, pigauto 0.11.0.9000 from ~/R/lib-lambda. Launched 05:14 MDT 2026-09-23, done
05:18. 20 of 20 result files, 0 errors. Seconds per job (all three arms): median 221, max 286.

Paired against the committed run on the same 20 seeds (aggregated with the same script):

| arm | trait | committed zRMSE | lambda default | delta |
|---|---|---|---|---|
| gnn_off | c1 | 1.001 | 0.956 | -0.045 |
| gnn_off | c2 | 0.994 | 0.962 | -0.032 |
| gnn_off | cnt | 0.965 | 0.951 | -0.015 |
| gnn_off | prp | 1.003 | 0.942 | -0.061 |
| gnn_off_rphylopars | c1 | 0.968 | 0.919 | -0.049 |
| gnn_off_rphylopars | c2 | 0.996 | 0.916 | -0.081 |
| gnn_on | c1 | 1.009 | 0.953 | -0.056 |
| gnn_on | prp | 1.004 | 0.957 | -0.047 |

Binary and categorical accuracy: delta 0 in every arm (the discrete path is unchanged by design).
Ordinal accuracy moves by at most 0.004 in the in-house arms. Conformal coverage for c1 and c2 stays in
0.95 to 0.97.

Runtime decision: 221 s per job projects the full 3,600-job core run to about 9 h at 36 parallel, over the
5 h stop rule. The no-GNN arms use one thread each, so they run as a separate wave at 140 parallel (under
the 150-core cap); gnn_on runs as its own wave at 36 parallel x 4 threads.
