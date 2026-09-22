# GOAL (re-read every arc; survives compaction)

Plan: /Users/z3437171/.claude/plans/distributed-snacking-turtle.md (approved by Shinichi 2026-09-20).
Worktree: /Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim, branch arc/imputation-sim from origin/main fd0b513.
Lease: claude:pigauto:84704 (release: LANE_ID='claude:pigauto:84704' ~/shinichi-brain/tools/lane_lease.sh --release pigauto).
Ledger: .unlazy/imputation-sim/gates/leaf-{env,runner,prerun,campaign,results}.md; verify with
  node ~/.claude/skills/unlazy/scripts/gate-check.mjs --reverify <leaf>.

Destination: the four-arm imputation simulation is computed on the corrected design, aggregated with paired MCSE,
and Shinichi has (a) a results Artifact, (b) a BACE-paper methods note, (c) an unlisted pkgdown article, plus
after-task, Melissa reconcile and handover.

Arcs (state on disk here; update the line when an arc closes):
- S1 ledger + worktree + lease ................ DONE 2026-09-20
- S2b Totoro env (agent totoro-env) ........... RUNNING
- S2c nibi env (agent nibi-env) ............... RUNNING
- S3 runner corrected design (agent runner-build) RUNNING
- S3-verify (fresh Opus reviewer, symbolic alignment) pending
- S5 drivers (design.R, totoro.sh, nibi_array.sh) DRAFTED 2026-09-20; smoke against the wrapper once S3 lands
- S4 pre-run Totoro 16 cells, note, STOP for G0 . pending
- G0 Shinichi .................................. pending
- S6a core Totoro / S6b factorial split nibi (half A) + rorqual (half B), fir hot backup / S6c AVONET / S6d covariates  pending
  (Shinichi 2026-09-20 "think backup plans, parallel"; measured: nibi 38k idle, rorqual 6k idle + fairshare 1.31,
   fir 8.6k idle, narval fairshare 0.09; Totoro 383/384 idle. Stall rule: a cluster with no task started 2 h after
   submit hands its remaining cells to the other or to fir; finished (cell, seed) rds are never re-run.)
- S2c' rorqual bootstrap: replay the nibi bootstrap script once nibi's smoke passes (own /project library)  pending
- MECHANICAL-VERIFY (Haiku) ..................... pending
- S7 aggregate; S7a Artifact (Shinichi decides); S7b BACE methods; S7c article  pending
- S8 after-task, Melissa, handover, PR .......... pending

Compute allowances: Totoro 250 cores for snakagaw (Shinichi, 2026-09-20, raising this lane above
D-143's standing 150) = 62 concurrent cells at 4 threads. nibi and rorqual: def-snakagaw_cpu arrays.
rorqual note: /project def-snakagaw is at its FILE-COUNT quota (about 499k/500k inodes), so its R
library and torch home live under /home there.

Pauses for Shinichi: G0 (after the pre-run note) and S7a (after the Artifact). Otherwise autonomous.
Must stop: edits to R/, BACE/, PR #175 files, the dirty main checkout; a Duo prompt; a merge; publishing the
article before Shinichi reads the Artifact.
Design facts locked: arms 1 freq (Rphylopars + castor + phyloglm Poisson), 2 BACE, 3a/3b pigauto GNN off, 4 GNN on;
fixed thresholds; trimmed factorial = 56 cells (clade included); primary contrast = BACE vs freq, z-RMSE + 95%
coverage (interval score), core slice, pooled over types.
