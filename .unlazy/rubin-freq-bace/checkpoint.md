GOAL: see GOAL.md (+ arcs.md).   STATE: campaign running on nibi + fir; rorqual/narval validating (2026-09-24 ~17:50 nibi).
ARCS DONE (verified): pre-run (240/240), freq timing (36/36), results page v1.
ARC IN PROGRESS: A1 (nibi 22634021), A2 (nibi 22634022/22634023), A3 (nibi 22633907/909, root pigauto_rubin_freq),
  A4 freq campaign on fir (61365818/19/20; n1000 seeds 1-200, n100/300 seeds 51-200; root ~/pigauto_rubin; env ~/pigauto_sim/env.sh).
  rorqual/narval: R lib copied from fir (BACE digest identical), env ~/projects/def-snakagaw/snakagaw/pigauto_rubin/env.sh,
  validation jobs rorqual 21771783, narval 3935076 -> if G-S4a PASS, they take part of A5.
STATUS: bash script/rubin_status.sh (Mac). Totoro: this user 155 cores busy in other lanes (>150, D-143) -> not used; recheck.
NEXT: A3 done -> pull nibi pigauto_rubin_freq/results/freq -> compare.R -> interim page; A1 done -> A5 split by seed across
  nibi / fir / rorqual / narval (per-host seed ranges, never overlapping; pool per host subdir, dedupe by filename).
OPEN GATES: A5 overrun gate (A1 > 30% over 6.9 h/fit). Never contact Dan; never edit R/ or BACE/.
TRUTH LIVES IN: branch arc/rubin-freq-bace; results on each host under pigauto_rubin/results/{bace,freq}; page https://claude.ai/artifact/WMHxQ8idafD5sLAVEJtg9Q
RESUME: Read .unlazy/rubin-freq-bace/GOAL.md, arcs.md, checkpoint.md; run script/rubin_status.sh; continue from NEXT.
HOST SPLIT (2026-09-24 19:5x nibi): BACE n300 seeds 1-102 on nibi (22634023, 204 tasks after scancel of the 198 PENDING
  seed 103-200 elements), seeds 103-200 on fir (61366948, 294 tasks, BLOCK 2). BACE n100 all seeds on nibi (22634022).
  freq: n100/300 seeds 1-50 on nibi root pigauto_rubin_freq (22633907/909), the rest on fir (61365818/19/20).
  A5 (n1000 BACE) will split seeds 2-100 across nibi/fir/rorqual/narval by disjoint ranges.
HOST SPLIT v2 (20:0x nibi; nibi BACE arrays still all PENDING after 2 h, fir/rorqual start in minutes):
  BACE n100: seeds 1-100 nibi (22634022, 120 tasks), 101-200 rorqual (21772063, 120 tasks).
  BACE n300: seeds 1-51 nibi (22634023, 102 tasks), 52-102 rorqual (21772064, 102 tasks), 103-200 fir (61366948, 294 tasks).
  rorqual validation 21771783 G-S4a PASS. narval validation FAILED: illegal instruction in ape (fir-built library not portable
  to narval CPUs) -> narval not used (rebuild would be needed). Totoro still over this user's 150 cap.
BACE LAMBDA = 1 (2026-09-24 ~20:50 nibi): BACE as shipped fails most lambda = 1 fits (fixed thresholds leave empty discrete
  levels; MCMCglmm "Mixed model equations singular", or "argument is of length zero"). v1 saw the same (114 bace failures n100).
  FIX (rubin_bace.R cc7daf75, gate G-S3c): fit_bace_mi drops unused factor levels and leaves out a discrete trait with < 2
  observed classes; recorded in diag$input_fix; no-op otherwise. Failed as-shipped rds moved to results/bace_asshipped_failed/
  (the as-shipped failure-rate record). Retry arrays (retry1): fir 61371997 (n300, 69), rorqual 21774491 (n100, 71),
  21774492 (n300, 27). nibi's queued lambda = 1 tasks run the fixed code. TODO after first pass: rerun
  rubin_retry_prep.R on every host (tasks that started before the deploy fail as-shipped) and resubmit (retry2).
FREQ BLOW-UPS: 5/600 freq-cmp fits (0.8%) have astronomically wide intervals on nibi compute; Mac reproduces the same
  cell cleanly; not the BLAS backend (fir BACE test). Diagnostic job nibi 22635770 pending. Fix still to decide.
FREQ FIX (2026-09-24 ~21:10 nibi): blow-ups were not nibi-only (fir: 6 rows up to 2e67, 8 fits with errors, of 3000).
  Cause: degenerate bootstrap refits (lambda* = 1 + near-singular Sigma_p) and a solve() + silent ginv() fallback. Fix in
  rubin_freq.R (0bca6ad1): Cholesky solves + hard bounds (cond var <= marginal var; |cond mean - mu| <= 20 marginal SD);
  freq A redraws the bootstrap sample up to 5 times (n_fail, n_degenerate); freq B records a failure. Gate G-S2d.
  rubin_cell.R 9dfb7926 records n_degenerate. FREQ RERUN of record: fir ~/pigauto_rubin_f2 jobs 61373104/05/06 (all 3600).
  Superseded (kept, not used): fir ~/pigauto_rubin/results/freq, nibi pigauto_rubin_freq/results/freq.
FREQ v3 (2026-09-24 ~21:50 nibi): v2 still had 2/2858 blow-ups: refits returning absurd but self-consistent Sigma_p
  (bounds scale with it). Added plausible_pars() on the implied tip variance: refit vs original fit factor 50 (normal
  0.5-1.7); original fit vs observed variance factor 1e4 (prp is legitimately inflated up to 204x at lambda 1: shared
  lambda + prp's non-phylogenetic noise, a frequentist-model finding for the report). rubin_freq.R 431cab36. Freq of
  record now fir ~/pigauto_rubin_f3 (61376909/10/11); v2 (f2) superseded. At final pooling: move the Mac pool's old
  freq/ aside before rubin_pool.sh.
FREQ v3 DONE (2026-09-24 ~22:25 nibi): 3599/3600 (1 running); max per-cell width 4.41 SD, max downstream CI 0.87 (no
  blow-ups); 12 degenerate draws + 59 refit failures redrawn; 3 fits failed in both arms (2 implausible original
  Rphylopars fits, 1 Rphylopars type error) = the frequentist failure rate to report.
  BACE lambda = 1 with the input fix: fir 47/47 new-code fits needed cleaning and completed.
REPORT PIPELINE READY (2026-09-24 ~23:00 nibi): bash script/rubin_pool.sh -> Rscript script/rubin_campaign_aggregate.R ->
  bash script/rubin_report_build.sh -> docs/dev-log/arc/2026-09-25-rubin-campaign-report.html (interim chip until
  bace >= 3000 and freq >= 3600). Publish as a NEW artifact (title "BACE and Frequentist MI"), update it at the end.
NEXT (in order): (1) n1000 timing check done (~02:40 nibi) -> gate (<= 9 h/fit incl chain, MaxRSS) -> launch A5 split:
  seeds 2-100 across fir / rorqual / nibi by disjoint ranges (seed 1 from timing check; rerun seed 1 for lambda = 1 only if
  it failed as-shipped). (2) first pass of n100/300 done -> retry-prep pass 2 on every host -> resubmit. (3) publish interim
  report once n100/300 BACE complete. (4) final: pool, aggregate, report, after-task, Melissa, commit, push, tell Shinichi.
