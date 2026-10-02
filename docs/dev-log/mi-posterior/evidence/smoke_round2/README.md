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

## The 14 cells re-run by accident in the round-2 campaign

The campaign driver (`campaign_logs/driver.sh`) passed a single rep number for regimes 1 and 11
("118", "37"). `12_totoro_campaign.sh` reads a single number as a count, so it planned reps 1 to 118
and 1 to 37 (`campaign_logs/phase2_regime_1.log`, `phase2_regime_11.log`). Those two runs were
stopped after regime 1 reps 1 to 10 and regime 11 reps 1 to 4 had finished; reps 118 and 37 were then
run on their own. The 14 finished cells had all converged at 69670d4. Their 9597e18 copies are
byte-identical to the 69670d4 files on the five compared fields, so the merge (which keeps the later
copy) gives the same summary either way (`../summarise_round2.log`).

Files:
- `identity_check.R`, `identity_check.log`: the comparison for the smoke cells;
- `identity_extra.log`: the verdict line for each of the 14 accidental re-runs (regime 1 reps 1 to 10,
  regime 11 reps 1 to 4), all `IDENTITY_OK`; the loop that wrote it was not recorded;
- `identity_extra_recheck.log`: the same 14 comparisons re-run on Totoro with `identity_check.R`,
  with the command and the per-field output; all `IDENTITY_OK`;
- `campaign_logs/`: the round-2 driver and the two phase-2 logs above, copied from Totoro;
- `cell_*.log`, `regime_*.rds`: the cells.

Regime 1 rep 1 appears both here (smoke run) and among the 14, so the identity checks cover 16
comparisons of 15 distinct converged cells.
