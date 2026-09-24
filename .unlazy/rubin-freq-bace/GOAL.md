# GOAL (lane: freq vs BACE under proper MI + Rubin's rules; re-read every arc)
Worktree /Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace, branch arc/rubin-freq-bace from arc/imputation-sim 14b919d.
Plan: docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md. Lease claude:pigauto:rubin (script/rubin_, docs/dev-log/arc/2026-09-24-rubin, .unlazy/rubin-freq-bace/).
State 2026-09-24 (ultra-plan approved; plan /Users/z3437171/.claude/plans/tidy-beaming-allen.md): Q1-Q3 YES. Locked in planning:
freq B = JOINT conditional draws at one fixed fit; coverage vs population rho (true slope = true phylo corr = rho in this DGP);
correlation = GLS-whitened phylogenetic correlation on Fisher z. Lease claimed (script/rubin_, script/tests-rubin/, docs/dev-log/arc/2026-09-24-rubin, .unlazy/rubin-freq-bace/).
Gate ledger: .unlazy/rubin-freq-bace/GATES.md + gates/leaf-*.md (untracked on purpose). S1 lib / S2 freq / S3 BACE dispatched in parallel;
next S4 rubin_cell.R + Mac smoke (n=60, M=20), S5 BACE pre-run plan, then STOP for Shinichi (D-139). Never edit R/ or BACE/. Totoro shared: check free -g first.
