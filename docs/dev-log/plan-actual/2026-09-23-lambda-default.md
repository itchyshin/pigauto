# Plan vs actual: Pagel's lambda default (D-278/D-281), 2026-09-23

Plan: `~/.claude/plans/eager-weaving-pike.md`. Actual: worktree `pigauto-lambda-default`, branch
`feat/joint-lambda-default`, draft PR #187 (open, not merged). Reconciled at the point Rose's final
review fixes had landed (commit 39a97cf) and the GNN benchmark wave had finished, before Melissa's own
slice (S11) and close-out (S12) were declared complete.

Six commits ahead of origin/main: b8eb384 (feature), 18c83f2 (pre-run/plan), bcd8069 (real-data bench),
379a95d (Rphylopars guard), 6bd25b6 (gate decision), 39a97cf (final review fixes). 49 files changed,
6371 insertions, 312 deletions.

## Summary counts

adaptive: 4
drift: 4
unclear: 1

## Drift items (one line each)

1. After-task report section 5 claims "15 of 15 gates met" for leaves S2 to S6; the ledger's own leaf
   files list 14 gates for those leaves (G1, G2, G2b, G3, G4, G5, G5b, G6, G7, G7b, G7c, G8, G9, G10),
   all eventually MET. Owner: Ada (orchestrator, wrote the after-task report).
2. The plan's S9 (Haiku mechanical verify) was never dispatched; Node.js was missing so
   `gate-check.mjs` could not run, and Ada re-verified the ledger itself with a Python stand-in,
   covering only leaves S2 to S6. Leaves S8 and S12 still show `[ ]` unchecked with `EVIDENCE: pending`
   for G11, G13, G14 and G15 in `.unlazy/lambda-default/gates/leaf-S8.md` and `leaf-S12.md`; only G12
   carries an explicit `ABANDON` line. Owner: Ada (chose to self-verify rather than escalate the
   tooling gap or find an equivalent independent check).
3. PR #187's body is stale against the branch's current head. It still reads "Still running on
   Totoro: the GNN arm, 18 cells x 100 seeds" and cites the full-suite/check evidence from commit
   b8eb384 (PASS 2565), even though the GNN wave finished at 14:18 and the review-fix commit 39a97cf
   raised the final count to PASS 2594, 0 errors / 0 warnings / 1 note. Rose's final review flagged
   this exact staleness ("neither log is evidence for HEAD") and recommended waiting for the finished
   GNN wave before merge readiness; the body was not refreshed afterward. The body also still frames
   the G12 shortfall as one failed criterion (the lambda-0.3 half-gap) where `benchmark.md` documents
   two (the half-gap and the |delta| < 0.01 guard at lambda = 1, which moved 0.021 to 0.033).
   Owner: Ada (S12 owns the PR).
4. At the moment of this reconcile the worktree carries four uncommitted, unpushed changes: a fix in
   `R/fit_baseline.R` (a fully observed continuous column now reports its estimated lambda instead of a
   default 1), a matching new test in `tests/testthat/test-lambda-dispatch.R`, and the finalized prose
   in `docs/dev-log/lambda-default/benchmark.md` and `real-data.md` (GNN wave-2 results, the corrected
   AVONET 300 number). The after-task report's Issue Ledger (section 7a) already lists the
   fully-observed-column bug as "Fixed in this lane," but the fix exists only in the working tree, not
   on the branch or the PR. No AGENT_LOG entry for this lane's close and no lease-release evidence were
   found either. Owner: Ada (responsible for committing before declaring a fix landed, and for S12
   close-out).

## Adaptive items (plan explicitly allowed for these, or the process worked as designed)

1. G12 (ship only if lambda 0.3 closes half the committed gap) failed (17 to 28% closed against 50%).
   The plan's own rule for this case is escalate and pause, not silently ship or silently revert; that
   is what happened. Shinichi decided 2026-09-23 to keep `lambda_mode = "estimate"` as the default
   anyway, overriding the plan's stated default answer ("pause") to question 3. Decision owner:
   Shinichi.
2. The GNN benchmark wave was capped at 100 seeds per cell (50 at n = 1000), short of the committed
   yardstick's 200, because another Totoro user's job cut the lane's share of cores mid-run. Decision
   owner: Shinichi (the cap), trigger: external contention, not an agent choice.
3. A real-data benchmark (`script/bench_lambda_datasets.R`, seven trait datasets, 13 cases) was added
   beyond the plan's single Totoro simulation benchmark. This is additional evidence, not scope creep;
   nothing under the fenced paths changed.
4. Rose's D-43 final review returned NEEDS-CHANGES and a fresh Sonnet builder was dispatched to fix the
   findings (ordinal leak, partial `lambda_fixed`, redundant eigendecomposition, weak tests, stale
   roxygen, NEWS gaps). This slice was not named in the original table but is exactly the rework loop
   the plan's verification discipline implies after an adversarial review finds real bugs; the review
   worked as intended (it caught a real correctness bug, the ordinal leak).

## Unclear

1. Whether the plan's S12 close-out items (AGENT_LOG line, `lane_lease.sh --release`) have actually
   run. No matching AGENT_LOG entry for this lane's close was found (only the earlier D-278/D-281
   renumbering entries), and no lease file was found in a quick search. Given the uncommitted files
   above, S12 most likely has not finished, but this was not confirmed against a live lease registry.

## Axis-by-axis detail

| axis | plan | actual | verdict |
|---|---|---|---|
| scope | Joint solver + covariate path + dispatcher + default flip + Totoro benchmark; DEFER OU/EB/kappa/delta, Sigma_P/Sigma_E REML, GNN changes, script/campaign_*, BACE/, docs/dev-log/arc/2026-09-2*, arc/imputation-sim | Matches. `git diff --stat origin/main...HEAD` touches none of the fenced paths; Rose's independent scope check (section 5 of her review) confirms this. `R/predict_pigauto.R`, listed in the ledger's OWNS list, needed no edit (the fit already carries its own baseline). One added real-data benchmark beyond the plan (adaptive item 3). | matches, one adaptive addition |
| evidence/verification | S9 Haiku mechanical verify over all leaves before S10's adversarial pass; ledger CHECK/EXPECT/EVIDENCE lines filled before any "done" claim | S2 to S6 leaves reverified (13:30 to 14:45, some initially UNMET then fixed); S8's G11 substance is documented in `prerun.md` and `benchmark.md` but its formal EVIDENCE line is still "pending"; G12 is explicitly `ABANDON`ed (compliant with the plan's own escape hatch); S12's G13 to G15 (PR gate, after-task validator, tracked-docs count) were never run through their CHECK commands. See drift items 1 and 2. | drift |
| model routing | S0 Haiku, S1r Opus, S2/S3 Sonnet high (S4 to S6 reuse the same two agents), S9 Haiku, S10 Opus, S11 Sonnet | S0, S1r, S2, S3, S4 (reuse S2), S5/S6 (reuse S3), S10 and S11 ran as planned. S9 (Haiku) never dispatched; the orchestrator (Fable-tier) self-verified instead. A Sonnet fix builder for Rose's NEEDS-CHANGES findings was an unplanned but reasonable addition (adaptive item 4). | drift on S9 only |
| safety gates | G12 half-gap-or-pause; Totoro <= 150 cores (D-143); GNN/no-GNN arm seed counts implicitly 200 to match the committed yardstick; MUST STOP list (no merge, no edits to fenced paths, no lease violation) | G12 failed and was escalated and overridden by Shinichi (adaptive item 1); GNN wave capped at 100/50 seeds by Shinichi under Totoro contention (adaptive item 2); Totoro core cap respected throughout (140 single-thread wave 1, 120 threaded wave 2); no merge; no fenced-path edits. | adaptive, gates functioned as designed |
| public claims | PR body claims should match the evidence files; ledger EVIDENCE should be verbatim | Rose's final review caught several overclaims (discrete "identical" scoping, coverage range, ordinal leak) and most were corrected in `benchmark.md` and `real-data.md`. The PR body itself, however, is stale relative to the current branch head and still undercounts the G12 failure to one criterion. See drift item 3. | drift |
| handoff state | Draft PR open, after-task report on disk, Rose's review on disk, Melissa's reconcile on disk, lease released, AGENT_LOG line written | Draft PR open and targets main; after-task report and Rose's review are on disk; this reconcile is being written now. Four files are uncommitted and unpushed as of this reconcile, one of which the after-task report already describes as fixed; no AGENT_LOG close-out line or lease-release evidence found. See drift item 4 and the unclear item. | drift, plus one unclear item |
