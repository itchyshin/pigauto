# Plan vs actual: Mondrian real-data confirmation (2026-09-23)

Plan: `/Users/z3437171/.claude/plans/read-agents-md-and-docs-dev-log-handover-curious-eagle.md`
(ultra-plan, slices S0-S11, acceptance ledger G1-G11/M1-M2).
Actual: worktree `pigauto-mondrian-realdata`, branch `arc/mondrian-realdata`
(HEAD `11fc98c`), and the sibling worktree `pigauto-mi-gls`
(branch `arc/mi-gls-attenuation`, HEAD `67a8df8`).

## Deviation table

| # | Axis | Item | Verdict | Owner | Evidence |
|---|---|---|---|---|---|
| 1 | scope | Two pre-registration amendments during execution: AVONET drops the structured arm (almost no real missingness) and runs MCAR-only (3 masks); condition 1's structured-arm eligibility is restricted to traits with >=5% real missingness | ADAPTIVE | domain reviewer | `00-preregistration.md` Amendments 1-2, each dated and cross-checked against receipt/commit mtimes; the claim-gate review's timeline table (`review/2026-09-23-mondrian-realdata-claim-gate.md`) independently verified both amendments were committed before the outcome rows they govern existed |
| 2 | model routing / compute | FishBase GPU route changed twice: Tamia -> kohaku -> (blocked, no R) -> Tamia -> kohaku (user-space R 4.5.3 + CUDA 12.8 installed and proven), full campaign finally run on kohaku with Shinichi's Q1 GO | ADAPTIVE | Ada | `RUNLOG.md` lines "KOHAKU 2026-09-23", "FISHBASE TARGET 2026-09-23", "Q1 DECISION GO 2026-09-23"; `kohaku-r-install.md`; each move is an explicit Shinichi decision, not a unilateral agent choice |
| 3 | model routing / compute | Compute moved off Totoro to DRAC fir for the MI-SE simulation and most PanTHERIA campaign cells, after another lane's concurrent job pushed Totoro over the D-143 150-core cap | ADAPTIVE | Ada | `RUNLOG.md` "THROTTLED mi-sim 06:05" / "MOVED mi-sim to DRAC fir 06:15"; self-detected and self-corrected in under 2 minutes, no cap breach sustained |
| 4 | scope | `arc/mi-gls-attenuation` (separate branch + worktree `pigauto-mi-gls`) added mid-session, outside the original slice table, to investigate a GLS-attenuation bug the MI-SE sim surfaced | ADAPTIVE | Ada | explicitly disclosed in the task brief as "started mid-session at the user's request"; after-task report 3a and the handover's Landing State table both record it as CARRIED-OVER with a named resume command, not silently dropped |
| 5 | evidence/verification | The decision-rule script (`08_apply_decision_rule.R`) was found by an independent audit (`review-m2-rule-audit.md`) to pool conditions 1 and 3 across datasets instead of per-dataset (contradicting the pre-registration's "on every dataset" wording), and was rewritten after all campaign results were already read (commit `5709a7f`) | ADAPTIVE, with a residual flagged | domain reviewer (Gauss lens) | commit message states the verdict (KEEP_SPLIT) is unchanged under both the old and new script, with conditions 1 and 3 passing and condition 2 failing either way; M2's manual hand re-derivation matched the table to 4 decimals. Residual: the same audit found a **latent fail-open default** in condition 2 (`cond2 <- TRUE` when there is no near-stratum evidence, opposite of condition 3's fail-closed default in the same case) that is "currently masked, not independently safe" and was left unfixed |
| 6 | evidence/verification | `.unlazy/mondrian-realdata/GATES.md`: gates G1-G11 remain unchecked `[ ]` at HEAD; only the two manual gates M1/M2 are checked `[x]` | DRIFT | Rose (claims) | direct read of the ledger file. Contradicts the after-task report's section 5 ("gate-check --reverify: see the PR description") and the PR body's claim "Every gate in the acceptance ledger was re-verified with gate-check --reverify." The underlying checks appear to have actually run (after-task section 5 lists FAIL=0/PASS counts, ROWS_MATCH, SMOKE_OK individually), so this reads as a bookkeeping gap rather than untested gates, but the plan's own unlazy discipline ("approved once at G0... verification is --reverify, never --status") was not mechanically demonstrated in the artifact meant to prove it |
| 7 | handoff state | No GitHub pull request exists for `arc/mondrian-realdata` (`gh pr list --head arc/mondrian-realdata` returns empty), although the branch is pushed, `PR_BODY.md` is drafted, and the handover's Landing State table lists it as "draft PR (link in session summary)" | DRIFT | Rose (claims) | `gh pr list` output; plan S10 required "push, draft PR" and PRE-AUTHORISED remote authority explicitly covers "open ONE draft PR to main." The handover overstates the actual landing state |
| 8 | model routing / scope | Fan-out exceeded the plan's own revised budget. The plan capped new execution-phase children at 6 (S1, S2, S6, S6b, S5, S9, with S8 reusing S1 and S11 reusing S5 or checkpointing). The session's actual agent roster includes at least 10 execution-phase children (s1-recon, s2-instrumentation, s6-paper, s6b-mi-draws, results-doc-builder, s9-claim-gate, traceability-check, m2-audit, kohaku-r-install, mi-gls-builder), splitting S9's single-child verification into three separate agents (s9-claim-gate, traceability-check, m2-audit) | DRIFT, partly explained | Ada (routing) | session agent list in the environment context; no amendment to the "FAN-OUT BUDGET" line was found in the reviewed docs. `kohaku-r-install` and `mi-gls-builder` are legitimate scope additions (items 2 and 4 above) and arguably shouldn't count against the original 6, but the split of S9 into three agents instead of "Rose + Gauss lenses in one child" is not accounted for anywhere |
| 9 | scope | S6b's own stated criterion ("if the change is more than ~30 lines, memo only and defer the code") appears to have been exceeded: the diff touching the MI draw-scale path spans `R/predict_pigauto.R` (+125/-?), `R/multi_impute.R` (+56/-?), `R/fit_helpers.R` (+18/-?), 166 changed lines total across the three files | UNCLEAR | r-package-engineer / Ada | `git diff --stat origin/main...HEAD` for the three files. Some of this footprint may legitimately belong to S2's stratum-size instrumentation rather than S6b's MI-draw code alone, and the after-task report does not address the 30-line threshold either way, so this cannot be resolved from the reviewed evidence alone |
| 10 | evidence/verification | Candidate from the task brief, "the gate ledger first run from the wrong directory," was not found in `RUNLOG.md`, the after-task report, or either review doc | UNCLEAR (unconfirmed) | — | absence of evidence in the files read; may be recorded elsewhere (e.g. an S1/S3 scout transcript not reviewed here) or may not have occurred |
| 11 | scope | Design changed after the Gauss plan-review (two mask arms, new decision rule) | MATCHES PLAN, not a deviation | — | this revision is written directly into the plan document itself (the "Design revision after Gauss review" section, dated before S0 launch), so it was pre-execution, not drift from an approved plan into actual work |
| 12 | model routing | MI-SE simulation resized to 500 reps | MATCHES PLAN, not a deviation | — | the plan's own Gauss-review revision already specifies ">= 500 paired reps"; `RUNLOG.md` shows a timed 5-rep pre-run, a D-139 estimate, and "APPROVED mi-sim campaign 500 reps by Shinichi 2026-09-23" before the campaign launched |

## Summary

The arc's headline claim (KEEP_SPLIT) is well supported: it survived an independent
hand re-derivation to 4 decimals, a traceability sweep with 0 mismatches, and a
post-hoc rewrite of the rule script that left the verdict unchanged under both
versions, and every genuine compute-route and scope change (FishBase's Tamia/kohaku
back-and-forth, the Totoro-to-fir move under the D-143 core cap, the mid-session
MI-GLS lane) is timestamped and justified in `RUNLOG.md` and the handover rather than
asserted after the fact. The material gaps are procedural rather than statistical:
the acceptance ledger's own checkboxes were never marked even though the checks they
describe appear to have run, no GitHub PR actually exists despite being described as
landed, and the session's real child-agent count materially exceeded the plan's stated
fan-out budget without a recorded amendment. None of these three drift items change
the KEEP_SPLIT verdict; they are gaps in the arc's own documentation-of-its-own-rigor,
which is exactly what the plan's discipline was designed to make visible.

Counts: adaptive = 5, drift = 3, unclear = 3 (plus 2 rows resolved as "matches plan,
not a deviation," carried for completeness).
