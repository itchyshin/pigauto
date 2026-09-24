# Round-2 smoke run (2026-09-24, Totoro, code 9597e18)

Before the round-2 campaign (design.md 5e), five cells ran at the package defaults.

| Cell | Purpose | Result |
|---|---|---|
| regime 1 rep 1 | byte-identity with 69670d4 | `IDENTITY_OK` (results, cell detail, cell coverage, core diagnostics identical); no extension |
| regime 12 rep 1 | byte-identity with 69670d4 (both-missing, plug-in arm) | `IDENTITY_OK`; no extension |
| regime 23 rep 196 | previously non-converged cell | extended once, converged; its numbers differ from 69670d4, as intended |
| regime 27 rep 1 | twin of regime 3 (n = 1000) | converged first time; per-cell coverage 0.962 (conformal 0.986) |
| regime 36 rep 1 | twin of regime 12 (n = 1000, both missing) | converged first time |

Wall time per cell under a Totoro load of about 180 (mostly other users): 4.5 min at n = 300,
about 11 to 12 min at n = 1000, and 8.3 min for the extended n = 300 cell. Peak RSS was 0.53 to
0.73 GB.

Files:
- `identity_check.R`, `identity_check.log`: the comparison;
- `cell_*.log`, `regime_*.rds`: the cells.
