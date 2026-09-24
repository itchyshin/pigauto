# GOAL (lane: freq vs BACE under proper MI + Rubin's rules; re-read every arc)
Worktree /Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace, branch arc/rubin-freq-bace from arc/imputation-sim 14b919d.
Plan: docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md. Lease claude:pigauto:rubin (script/rubin_, docs/dev-log/arc/2026-09-24-rubin, .unlazy/rubin-freq-bace/).
State 2026-09-24 13:20 (APPROVED + RUNNING): Shinichi approved the pre-run and a chained BACE arm ("yes go ahead").
bace_chain built (gate G-S3b), five-arm smoke green (Mac + nibi job 22614229). Totoro refused (this user ~200 cores in other lanes,
D-143), so the pre-run runs on nibi: arrays 22614567 22614569 22614570 22614571 (240 BACE fits, 421 core-hours est.).
Results on nibi: ~/projects/def-snakagaw/snakagaw/pigauto_rubin/prerun/. Next: read the first finished rds early; stop and re-report
if the first fits run >30% over their row (D-139); when done, aggregate per setting (convergence, ESS, failures, bace vs bace_chain
coverage) and apply the selection rule. Campaign still needs Shinichi's approval of settings. Never edit R/ or BACE/.

OVERNIGHT 2026-09-24 -> 05:00 (Shinichi away; "work autonomously and keep working to finish what you need to do"):
 1. Pre-run completes -> aggregate (script/rubin_prerun_summary.R on nibi) -> docs/dev-log/arc/2026-09-25-rubin-prerun-results.md
    with the full table + settings recommendation + convergence-check recommendation.
 2. Read the installed BACE convergence check; document what it tests (why it fails 70-90% with ESS 750-1600).
 3. Freq-arm timing pre-run on nibi: 18 core cells x 2 seeds, arms freqA,freqB, <= 30 min wall (D-139 ok) -> measured rates incl. n=1000.
 4. Campaign plan + launch scripts (NOT launched): arms, reps (BACE 200 per Meng N12?), budget from measured rates (+25% at n=300).
 5. Meng N10: prp per-cell scoring on the logit scale (tests).
 6. Commit, push, handover note for 05:00. MUST NOT: launch the campaign, message Dan, edit R/ or BACE/, write to the brain vault.

STATE 2026-09-24 17:30 nibi (overnight items 1-5 DONE): pre-run 240/240, results page https://claude.ai/artifact/WMHxQ8idafD5sLAVEJtg9Q
(copy docs/dev-log/arc/2026-09-25-rubin-prerun-results.html). Freq timing 36/36 (freqA 31/69/268 s per rep at n 100/300/1000).
Fast exact PGLS in the runner; prp logit scoring; campaign launcher script/rubin_campaign_nibi.sh (dry run, NOT launched).
WAITING ON SHINICHI: (1) drop BACE convergence rule, runs 5 nitt 50k; (2) arms incl. bace_resid as control; (3) option B
~7,400 core-h after a 6-fit n = 1000 BACE timing check; (4) when to tell Dan. Nothing launched beyond the pre-run.
