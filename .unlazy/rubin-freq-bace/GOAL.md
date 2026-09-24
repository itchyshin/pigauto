# GOAL (lane: freq vs BACE under proper MI + Rubin's rules; re-read every arc)
Worktree /Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace, branch arc/rubin-freq-bace from arc/imputation-sim 14b919d.
Plan: docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md. Lease claude:pigauto:rubin (script/rubin_, docs/dev-log/arc/2026-09-24-rubin, .unlazy/rubin-freq-bace/).
State 2026-09-24 close (steps 1-2 DONE, PAUSED for Shinichi): runner + four arms + Rubin pooling built; 12/12 runnable gates met;
Mac smoke green (n=60 seed 2; seed 1 BACE singular at n=60, reported). Pre-run plan READY, NOT LAUNCHED: 240 fits, ~253 core-hours,
2.5-3 h Totoro at 120 procs (docs/dev-log/arc/2026-09-24-rubin-prerun-plan.md; launcher refuses >150 cores for this user).
CORRECTION: installed BACE (2026-08-09) draws posterior predictive values; the in-tree BACE/ clone is stale. bace_resid = negative control.
OPEN for Shinichi: (1) approve the pre-run; (2) Meng B2: every BACE final run starts from one converged dataset (toy coverage 0.90 vs 0.95
chained); approve a bace_chain arm or not; campaign blocked on B2. Review: docs/dev-log/arc/2026-09-24-rubin-review.md.
After-task: docs/dev-log/after-task/2026-09-24-rubin-steps1-2.md. Never edit R/ or BACE/. Totoro shared: check free -g first.
