# GOAL (lane: freq vs BACE under proper MI + Rubin's rules; re-read every arc)
Worktree /Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace, branch arc/rubin-freq-bace from arc/imputation-sim 14b919d.
Plan: docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md. Lease claude:pigauto:rubin (script/rubin_, docs/dev-log/arc/2026-09-24-rubin, .unlazy/rubin-freq-bace/).
State 2026-09-24 13:20 (APPROVED + RUNNING): Shinichi approved the pre-run and a chained BACE arm ("yes go ahead").
bace_chain built (gate G-S3b), five-arm smoke green (Mac + nibi job 22614229). Totoro refused (this user ~200 cores in other lanes,
D-143), so the pre-run runs on nibi: arrays 22614567 22614569 22614570 22614571 (240 BACE fits, 421 core-hours est.).
Results on nibi: ~/projects/def-snakagaw/snakagaw/pigauto_rubin/prerun/. Next: read the first finished rds early; stop and re-report
if the first fits run >30% over their row (D-139); when done, aggregate per setting (convergence, ESS, failures, bace vs bace_chain
coverage) and apply the selection rule. Campaign still needs Shinichi's approval of settings. Never edit R/ or BACE/.
