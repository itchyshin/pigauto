# Arcs: freq-vs-BACE Rubin CAMPAIGN (option B, approved by Shinichi 2026-09-24 evening)
Decisions: (1) drop BACE's convergence check from the selection rule, BACE runs 5, nitt 50,000, burnin 10,000, thin 25;
(2) arms freqA, freqB, bace, bace_chain, bace_resid (resid = negative control); (3) option B: BACE 200 reps at n 100/300,
100 reps at n 1000, after a 6-fit n = 1000 timing check; freq 200 reps everywhere; (4) hold Dan until the campaign is in.
Cells: types_mixed, BM, MCAR 30%, lambda {0.3, 0.7, 1} x rho {0, 0.5}, driver on, thresholds fixed, M = 20.
nibi root: ~/projects/def-snakagaw/snakagaw/pigauto_rubin (results/bace, results/freq). Launcher: script/rubin_campaign_nibi.sh.

- [x] A1 n = 1000 BACE timing check (gate passed on fit phase, 1-12% under): job 22634021 (6 tasks, seed 1, 48G, 14 h). Measure wall + MaxRSS; est 6.9 h/fit.
- [ ] A2 BACE n = 100 (job 22634022, 240 tasks x 5 seeds, 6G, 5:30) and n = 300 (job 22634023, 402 tasks x 3 seeds, 12G, 10:00).
- [x] A3 freq comparison on pre-run datasets (jobs 22633907/22633909, seeds 1-50, root pigauto_rubin_freq) -> interim page update.
- [x] A4 freq campaign (superseded twice; freq of record = fir pigauto_rubin_f3 v3, 3600/3600): copy A3 rds into pigauto_rubin/results/freq (same code, md5 verified), then submit remaining
       (n 100 seeds 51-200 BLOCK 50; n 300 seeds 51-200 BLOCK 25; n 1000 seeds 1-200 BLOCK 10, 8G, ~1 h tasks).
- [~] A5 BACE n = 1000 campaign LAUNCHED fir 61390620 + rorqual 21781308;, seeds 1-100 (seed 1 from A1), mem/time from A1. GATE: if A1 > 30% over 6.9 h, STOP + re-report.
- [ ] A6 aggregate all; complete freq-vs-BACE report (page + docs/dev-log/arc/2026-09-25-rubin-campaign-report.md);
       after-task; Melissa reconcile; commit + push; tell Shinichi. Nothing to Dan.
