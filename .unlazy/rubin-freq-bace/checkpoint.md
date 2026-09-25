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
RETRY PASS 2 fir (2026-09-24 ~23:50 nibi): 8 lambda = 1 n300 failures moved, resubmitted as fir 61390210 (retry2). Pass-1
  as-shipped counts (from the retry-prep outputs, authoritative): fir 69 (68 lambda 1 + 1 lambda 0.7), rorqual 98 (all
  lambda 1). The pass-2 move on fir overwrote 6 pass-1 records = 6 fits that failed AGAIN after the fix (retry1). Fixed
  rubin_retry_prep.R: later passes use results/bace_failed_pass<k>/, never overwrite. Use dest bace_failed_pass3 next.
A5 LAUNCHED (2026-09-24 23:55 nibi): n1000 timing gate evaluated on the shipped-fit phase: 6/6 fits 12,136-13,540 s vs
  13,750 predicted (1-12% under), incl. both lambda = 1; chain phase = 20 of the same sweeps -> ~6.4-6.8 h/fit total vs
  6.9 h estimate; MaxRSS 16 GB so far. Launched without waiting for the chain phase (answers the time/memory question
  Shinichi's "timing check first" asked). n1000 BACE: fir seeds 2-60 (61390620, 354 tasks, 40G, 12h), rorqual seeds
  61-100 (21781308, 240 tasks). Seed 1 = timing check (nibi 22634021). Expected finish ~07:00-10:00 nibi if they start
  promptly. Watch the first finished fits' total wall and MaxRSS; stop and re-report if > 30% over (D-139).
RETRY RULE (2026-09-25 00:2x nibi): BACE fits are deterministic given seed + code (arm_seed 303), so retrying a failure
  on the FIXED code repeats it. Only failures from tasks that started on the OLD code (before the ~20:50 deploy) are
  retried; identify by the task log's "cell <tag>:" start time. rorqual: 17 old-code failures -> bace_failed_pass2/,
  resubmitted 21782231 (n100, 14) + 21782232 (n300, 3); 10 new-code failures stay in results/bace = after-fix failures.
  fir retry2 (61390210, 8) were mostly new-code failures (6 of 8 failed twice) and will repeat; harmless (fail in ~1 s).
  nibi: check old-code failures the same way at the end of its first pass. Pool now also pulls bace_failed_pass2.
N1000 HOST SPLIT v2 (2026-09-25 ~01:00 nibi): fir started 1/354 (40G/12h jobs queue slowly; low fairshare) -> cancelled
  180 PENDING fir tasks (seeds 31-60). Now: fir seeds 2-30 (61390620, 174), Totoro seeds 31-40 (60 fits, PAR 30,
  pgid 401996, ~/pigauto_rubin, driver script/rubin_totoro.sh, BACE digests identical), nibi seeds 41-60 (22645293, 120),
  rorqual seeds 61-100 (21781308, 240). Seed 1 = nibi timing check. Totoro counts in rubin_status.sh (totoro line).
  Pool: add Totoro results (rsync totoro:pigauto_rubin/results/bace -> pool/bace/totoro) at final pooling.
N1000 FIRST COMPLETE FIT (2026-09-25 01:36 nibi): timing task 5 total 5:55:49 (fit 12,136 s + chain 9,204 s), MaxRSS
  18.8 GB; estimate 6.9 h -> 14% under (D-139 pass). fir's 172 pending n1000 tasks (seeds 2-30 minus the 2 running)
  resubmitted at 24G/10h as 61453939 (two accidental duplicates of the running seeds 2,3 at lambda 0.3 rho 0 were
  cancelled before they started). Running n1000 at 01:40 nibi: fir 2, rorqual 78, nibi 70, Totoro 30.
INTERIM REPORT PUBLISHED (2026-09-25 ~03:35 nibi): https://claude.ai/artifact/PFkoRFtTjox4tndEbuBtPQ (file
  docs/dev-log/arc/2026-09-25-rubin-campaign-report.html; republish the same path to update). Rows with < 30 datasets are
  withheld as "running"; narrative is data-driven (accuracy claims need > 2 SE). Findings so far: BACE as shipped
  under-covers at lambda 0.3/0.7 (0.88-0.91 per cell), near nominal at lambda 1; chained BACE ~0.95 (over-covers
  slightly at lambda 1); freq A ~0.95 everywhere; freq B under-covers the correlation.
  Missing reruns: rorqual 21787578 (n100 seeds 51-55 lambda 0.7 rho 0.5, the nibi timeout block), 21788702 (n300:
  lambda 0.7 rho 0.5 seeds 55-57 = rorqual segfault block; lambda 0.7 rho 0 seeds 151-152 = fir segfault block).
N1000 TAIL CUT (2026-09-25 ~04:05 nibi): fir 24G resubmission 61453939 now 169 RUNNING (24G fixed fir scheduling).
  Totoro: stopped the xargs dispatcher (pid 402000) only; its 30 running fits (lambda 0.3 both rho + lambda 0.7 rho 0,
  seeds 31-40) continue; the unstarted second wave (lambda 0.7 rho 0.5, lambda 1 both rho, seeds 31-40) -> fir 61472904.
  rorqual: 59 PENDING n1000 tasks cancelled -> fir 61472913 (24G). rorqual keeps 180 running. nibi 118 running.
  All 594 n1000 fits now running or just submitted; expected done ~10:00-11:00 nibi (08:00-09:00 Edmonton).
BACE n100 + n300 COMPLETE (2026-09-25 ~06:00 nibi): 1200/1200 each (reruns of the timeout block and both segfault
  blocks all succeeded: segfaults were not deterministic). Report republished v3 (same URL) with n100/300 final and a
  per-n bias statement (BACE attenuates the slope at n = 100 vs freq A; within MC error at n = 300 for BACE as shipped).
N1000 (06:30 nibi): 75/600. Totoro batch 1 (30) done; Totoro batch 2 = 40 of fir's pending (61472913) moved,
  pgid 593856, file logs/n1000_from_fir.txt; fir keeps 22 pending (8 in 61472904, 14 in 61472913) + ~175 running.
