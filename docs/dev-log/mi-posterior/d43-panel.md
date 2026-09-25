# D-43 completion panel (2026-09-25, at commit 79aed52)

Three independent reviewers read the committed record and re-ran checks: two on Sonnet and one ceiling
reviewer on Fable. Rule: two or more NOT-DONE verdicts withhold a claim. **No claim was withheld.**

Claims judged:
- C1: implemented, tested and documented; G1 to G5c pass at the final code.
- C2: in-model per-cell coverage is about 95% (G7), with the stated under-coverage caveat.
- C3: the in-model downstream slope is close to unbiased with coverage near complete data; G6's relative SE-ratio rule fails in 3 rows, and not from noise.
- C4: every number in results.md traces to a committed file, and nothing is overstated.
- C5: the work is presented as not finished.

| Claim | sonnet-a | sonnet-b | ceiling (fable) |
|---|---|---|---|
| C1 | DONE | DONE | DONE |
| C2 | DONE | DONE | DONE |
| C3 | DONE | DONE | NOT-DONE |
| C4 | DONE | DONE | NOT-DONE |
| C5 | DONE | DONE | DONE |

## Dissents and how they were handled

- C3, ceiling-fable: Most of the claim is supported and reproduces from sim_summary.csv: 48 gated posterior_full rows, max |paired bias| 0.0142, all pass max(0.02, 2.5 MCSE); coverage 0.845-0.975 vs complete 0.760-0.970 with all 48 rows >= complete - 0.05 and mean positive shortfall 0.0016 (note the rule is one-sided: MI over-covers complete by up to +0.17 in the gls lambda = 0.5 rows, e.g. regime 36 gls 0.930 vs 0.760); all 16 lambda = 1 twin rows negative (-0.004 to -0.014, |z| 2.7-8.7) and stated as a caveat; 04_acceptance.R re-run fails exactly rows 35 (1.186), 36 (1.153), 38 (1.167) on the relative phylolm rule; 07_se_ratio_noise.R re-run is byte-identical to se_ratio_noise.txt (P(>= 3 rows out of band | noise) = 0.021, Cochran Q p = 0.012, bootstrap-checked MCSE), so 'not Monte Carlo noise' is supported at that strength. The overstatement is 'the pooled SEs are honest in absolute terms': over the 24 ga

- C4, ceiling-fable: Second half holds: results.md reports G6 'not met (3 in-model rows)' and G8 'pending', its status block says 'Not ready to merge', and its G5 row under-claims ('being re-run') relative to the ledger log in the same commit. Every simulation, calibration, real-data and diagnosis number I checked reproduces: ~60 numbers recomputed from sim_summary.csv and cell_coverage.csv (bias ranges, 16 lambda = 1 rows, 15/16 rows > 2 SE, 4-41% rescaled widths, 48/48 plug-in pairs, widths 0.002-0.18% narrower, 7,999/8,000, stress 72%/54%/16-75%, gls -0.0067 vs -0.0090); se_ratio_noise.txt byte-identical on re-run; G3 check reproduces A coverage 0.00 with bias -0.004/-0.003 and D rho bias -0.14; real_preview tables give PanTHERIA 0.87-0.96/0.90-0.98/0.91-0.99, means 0.931/0.948/0.953, 17 of 24, AVONET 0.957/0.959/0.962, 2.952 vs 2.858 at 4.05 ref SE, 0.148 vs 0.180 at -3.29, ESS 368/366; diagnosis.md line

Actions taken:
- C4: G1's counts are now committed in `evidence/ledger/g1_counts.log`.
- C4: the gate table in `results.md` cites the ledger re-verification log, and its stale "being re-run" and "not committed" notes are gone.
- C3: `results.md` already limits the absolute-calibration statement to the three failing rows. It reports the 0.80 to 1.06 range, below 0.90 in 5 rows, in the "other direction" bullet. The looser wording was in the panel's claim text, not in the documents.
