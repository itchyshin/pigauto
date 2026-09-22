# Gates: pre-run (S4, Totoro, under 30 minutes)

OWNS: docs/dev-log/arc/*-simulation-prerun.md

Scope: 16 cells (12 core cells at n in {100, 1000} plus 4 sentinels at n = 1000: OU lambda 0.3, MAR, clade, MCAR 0.1), 1 replicate, every arm, BACE 2 chains, on Totoro; walls, Rhat and a re-derived budget written to the note; then STOP for Shinichi (G0).

- [x] G9a: 16 pre-run rds present on Totoro, none empty
  CHECK: ssh -o BatchMode=yes -o ConnectTimeout=12 snakagaw@totoro.biology.ualberta.ca 'n=$(ls ~/pigauto_sim/prerun_fast/*.rds 2>/dev/null | wc -l); [ "$n" -eq 16 ] && echo G9a PASS'
  EXPECT: G9a PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=8b01951b15dd/30 entries; output=G9a PASS

- [ ] G9b: every BACE cell converges by BACE's own assess_convergence() verdict, with the runs setting
  recorded (Gelman-Rubin across `runs` is NOT the diagnostic; see bace_diagnostics()).
  CHECK: ssh -o BatchMode=yes -o ConnectTimeout=12 snakagaw@totoro.biology.ualberta.ca 'cd ~/pigauto_sim && Rscript script/campaign_sim_checks.R --gate G9b --dir prerun_bace'
  EXPECT: G9b PASS
  EVIDENCE: probe measured runs 3 = 0/2 converged, runs 5 = 2/2, runs 10 = 2/2 at n = 100. The six
  stored BACE cells ran at runs = 2, below assess_convergence()'s min_iterations = 3, so this gate is
  expected to fail against prerun_bace and is re-run on the campaign at runs = 5.

- [x] G9c: arm 3b wall at n = 1000 measured and the budget re-derived in the note (slot-hours for core and trimmed factorial)
  EVIDENCE: arm 3b = 260 s at n = 1000 (10 cells, 209-337 s); budget re-derived at 11,564 slot-hours in
  docs/dev-log/arc/2026-09-20-simulation-prerun.md (commit ffd6d85).

- [x] G0: Shinichi approved the launch after reading the pre-run note (his words recorded here with date)
  EVIDENCE: Shinichi, 2026-09-20: "Go ahead" -- approving the pre-run note
  (docs/dev-log/arc/2026-09-20-simulation-prerun.md), runs = 5, n_final = 20, and the proposed
  allocation. Launched 13:42 Totoro core fast arms (3600 jobs, 62 slots); nibi core BACE n=100
  (job 22339703) and n=300 (22339823); rorqual core BACE n=1000 (21475732, one replicate per task).
