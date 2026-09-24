# BACE settings pre-run: ready to launch, waiting for Shinichi (D-139)

Lane `arc/rubin-freq-bace`, step 2 of `docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md`. Nothing here
has been launched. The run takes about 2.5 to 3 hours on Totoro, so it needs approval first.

## Recommendation

Run the full grid on Totoro at 120 concurrent processes (if no other lane of this user is running there): 240 BACE fits, about 253 core-hours, wall-clock
**2.5 to 3 h** (the longest single fit sets the floor at about 2.4 h). A cheaper option drops n = 300:
70 core-hours, about 1 h, but convergence at n = 300 then goes unmeasured before the campaign.

## What it answers

Which BACE settings the campaign should use. v1 ran `nitt = 50000`, `runs = 5`, `n_final = 20` and passed
BACE's own convergence check too rarely. Per setting, the pre-run records BACE's convergence verdict, the
median ESS, the wall time, the failure count, and the per-cell Rubin coverage and width of both BACE arms.

Proposed selection rule: take the cheapest setting where (a) at least 80% of fits pass BACE's convergence
check, (b) median ESS is at least 100 (v1's G9b bar), and (c) per-cell Rubin coverage lies within Monte
Carlo error of the most expensive setting's. If no setting meets (a), report that rather than pick one.

## Grid

| factor | levels |
|---|---|
| `runs` | 5, 10, 15 |
| `nitt` | 50,000 and 100,000; burnin 20%, thin chosen to keep 1,600 samples |
| `n_final` = M | 20 |
| cells | lambda {0.3, 0.7} x rho {0, 0.5} x n {100, 300}; MCAR 30%; v1 core options (`--driver`, fixed thresholds) |
| seeds | 1 to 5 |

6 settings x 8 cells x 5 seeds = **240 fits**. The lane plan's "720 fits" was an arithmetic slip; its stated
grid gives 240. Arms per fit: `bace` and `bace_resid`, which share one BACE fit.

## Time estimate (derived, then checked)

BACE's cost scales with `nitt x (runs + n_final)`. v1's median wall times, 1,161 s at n = 100 and
3,050 s at n = 300, were measured at v1's settings, 50,000 x (5 + 20) = 1.25M units. That gives 0.93 ms per
unit at n = 100 and 2.44 ms per unit at n = 300. The Mac smoke (n = 60, 6,000 x 23 units, 88.5 s) measured
0.64 ms per unit, which fits between them.

| setting (runs, nitt) | units | n = 100 | n = 300 |
|---|---|---|---|
| 5, 50k | 1.25M | 19 min | 51 min |
| 10, 50k | 1.50M | 23 min | 61 min |
| 15, 50k | 1.75M | 27 min | 71 min |
| 5, 100k | 2.50M | 39 min | 102 min |
| 10, 100k | 3.00M | 46 min | 122 min |
| 15, 100k | 3.50M | 54 min | 142 min |

Total: 12,539 s per cell-seed at n = 100 and 32,940 s at n = 300, summed over settings; x 4 (lambda x rho)
x 5 seeds = **253 core-hours**, mean fit 63 min. At 120 processes the lower bound is 2.1 h, but the longest
fit is 2.4 h, so the wall-clock is about 2.5 to 3 h with longest-first ordering. Rubin scoring adds seconds
per fit. Uncertainty: v1's BACE timings come mostly from nibi, and Totoro's cores may be slower; if the first
n = 300 fits run more than 30% over their row above, the run stops and I re-report (D-139).

## Launch procedure

1. `rsync` the lane's `script/` to `totoro:~/pigauto_rubin/script/` through the `~/.ssh/cm-*` socket.
2. `bash script/rubin_prerun_totoro.sh 120` without `CONFIRM` first. It prints `free -g`, the heaviest
   users' memory, checks that the installed BACE draws posterior predictive values (`sample = TRUE` in
   `bace_final_imp`), checks the code is present, counts the cores this user already has busy, and stops.
   D-143's 150-core cap covers every lane run as this user: the lambda-default lane used about 140 cores
   on 2026-09-23, and if it is still running the script refuses 120 and names the largest PAR that fits.
   A smaller PAR stretches the wall-clock (for example about 5 h at 60), so wait for that lane or accept
   the longer run.
3. `CONFIRM=yes bash script/rubin_prerun_totoro.sh 120`. The run is detached with `setsid`; results land in
   `~/pigauto_rubin/prerun/runs<r>_nitt<nitt>/`, one directory per setting because the cell tag does not
   carry the BACE settings.
4. Read the first finished n = 100 rds before trusting anything else.

## Known risks

- **BACE can fail outright at small n.** In the Mac smoke, seed 1 at n = 60 failed deterministically with
  MCMCglmm's "Mixed model equations singular"; seed 2 ran. The runner records the error and still writes
  the rds, so failures are counted, not lost. The failure rate at n = 100 and 300 is one of the pre-run's
  outputs.
- **Same BACE build everywhere.** The in-tree `BACE/` clone (2026-04-01) imputes posterior means; the
  installed build (2026-08-09) draws posterior predictive values. The launch script refuses a host with the
  old build.
- **`bace_resid` is under review.** With the installed build it adds a second residual. It costs nothing to
  keep in the pre-run; whether it stays in the campaign is Shinichi's call.
