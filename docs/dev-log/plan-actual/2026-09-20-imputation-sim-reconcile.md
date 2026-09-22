# Plan-vs-actual: the four-arm imputation simulation lane (2026-09-20)

The lane matches the approved plan on the axes that matter most. The DEFER fence held: `git diff
--stat fd0b513..HEAD` touches only `script/`, `docs/`, `.Rbuildignore`, `vignettes/articles/`, and
`.unlazy/`, with no changes under `R/` or `BACE/`, and the main checkout on
`handover/2026-08-09-cursor` is exactly as dirty as it was before this lane started. G0 is recorded in
Shinichi's own words ("Go ahead") before any compute committed, and no publish or merge has occurred,
consistent with the overnight grant of "everything except publishing and merging." The campaign is
reported as running, not claimed done, and a real mid-campaign bug (freq/gnn_on co-failure on 59 of
2,334 core cells, all in the lambda = 1 mixed cells) was recorded rather than patched, since fixing it
would mean editing `R/`, which is out of scope. Four material deviations are worth flagging, none of
which compromise the fence or the publication gate.

| axis | planned | actual | tag | owner |
|---|---|---|---|---|
| model routing / compute | rorqual named as the sole DRAC fallback if nibi's first array waited past 2 hours; fir appears in the sweep receipt as an environment, not a campaign host | rorqual retired mid-campaign (project inode quota, slower nodes); fir became the primary host for the n = 1000 BACE arm and factorial half A. The substitution is reasoned and logged (after-task section 3a, PROGRESS section 8a) | adaptive | Ada |
| safety gate / budget | G9c re-derived the whole-campaign budget at 11,564 slot-hours; G0 approved that number | measurement on fir shows the n = 1000 BACE arm alone would cost roughly 28,000 core-hours at the original 16 GB / 3:30 setting, about 2.4x the whole approved budget. No resubmission has been made at n = 1000 BACE; three costed options are written down and held for Shinichi | adaptive | Ada, pending Shinichi |
| evidence / verification | the plan names the `.unlazy` ledger as the reverify source of truth ("gate-check.mjs --reverify ... exits 0 before every 'done'") | `leaf-runner.md` (G2 to G6, Gold, Galign) and `leaf-results.md` (G12, G13a to d) still read unchecked with "EVIDENCE: pending", while PROGRESS.md and the after-task report assert specific PASS numbers for several of them. A rerun today of G13a (superlative grep) and G13b (`slop_check.py`) on the current files confirms both pass, so the claims check out where they can be checked independently, but a `--reverify` run today would report these leaves unmet | drift | Rose |
| handover state | overnight authority: "push arc/imputation-sim and open a DRAFT PR" | the branch is 25 commits ahead of `origin/main` locally; `git ls-remote origin arc/imputation-sim` returns nothing and no PR exists in any state. This is not logged anywhere as blocked or deferred | drift | Rose |

The ledger gap is worth separating from a false claim. Where I could check the underlying assertion
directly (G13a and G13b), it held up. The concern is procedural: the plan's own reverify mechanism
reads these gates as unmet, so a session that trusted only `gate-check.mjs` rather than the prose
would not see what PROGRESS.md and the after-task already report. The unpushed branch is a plainer
gap: overnight authority explicitly covered it, three close-out documents were written in its place,
and the push itself, which touches no other lane, was not done.

## Carried into the next session

- Push `arc/imputation-sim` and open the draft PR; authorized overnight, not yet done.
- Update `leaf-runner.md` and `leaf-results.md` with the gate evidence already reported in
  PROGRESS.md and the after-task, then run `gate-check.mjs --reverify` for real against both leaves.
- Shinichi's decision on the n = 1000 BACE arm: resubmit at 32 GB / 5 h, reduce replication to seeds
  1 to 30, or report the scaling failure as a measured result.
- Resubmit nibi's 224 TIMEOUT cell-seeds at doubled wall time and 32 GB once its two arrays drain.
- S6c (AVONET300) and S6d (covariate sensitivity), both correctly deferred until the core slice
  finishes.
- Shinichi's read of the S7a Artifact, and his call on whether and when S7c is published.

## Reconcile 2, 2026-09-22 (Melissa mode; six axes only)

The 09-20 carried-over items first, since they are what this reconcile owes: the draft PR (#184) turned
out to exist since the evening of 09-20 and the branch is now current on it; the ledger leaves now carry evidence for G14, G6c, G13a, G13b, G13d with
G10, G11, G6d, G12 recorded as open with their causes; the n = 1000 BACE decision was taken (30
replicates, Shinichi 2026-09-21); nibi's TIMEOUT cell-seeds were resubmitted and have drained; S6c is
done and S6d is at 3,584 of 3,600; the Artifact read and the S7c publication call remain with Shinichi.

| axis | planned | actual | tag | owner |
|---|---|---|---|---|
| scope | four arms plus floor (plan, Q1) | five: `freq_lambda` added as a fifth arm on Shinichi's word ("OK changed my mind. Yes, at the lambda stuff as the fifth arm.", 2026-09-21), run on core, factorial, AVONET and covsens | adaptive | Ada |
| scope | factorial "40 cells after the Q7 trim" (plan) | the design table has always emitted 56 (the plan's 40 omitted MCAR 0.30); G11 and the article say 56 | adaptive, and the plan text was wrong | Ada |
| evidence / verification | S6b factorial fast arms on nibi arrays, concurrent with the core slice | the fast-arm wave did not run for about 24 h: the Totoro launch was refused three times overnight and each refusal was logged without the completion summary saying the wave had NOT run; found 2026-09-21 morning when the factorial pool proved BACE-only; launched 09:19, landed 2026-09-22 09:02 on Totoro, not nibi | drift | Rose |
| evidence / verification | "every reported number carries MCSE" and "failures scored at the floor" (plan) | the aggregator had no rule for a fit that returns a finite absurd value; 15 + 34 + 2 such replicates dominated whole strata of the first factorial summary. Rule added and committed (1ec3e98) before any factorial number shipped; count reported as `n_divergent` | adaptive | Ada, verified by the per-replicate cross-check |
| evidence / verification | G14: 5 cells on two hosts agree on truth, mask and the frequentist arm | PASS on three host pairs for truth and mask; the frequentist clause is vacuous because nibi and fir hold BACE only. Stated as such in the ledger and the after-task | adaptive | Rose |
| public claims | methods note reports each arm at its documented default | note says "No arm reaches nominal 0.95" for a paragraph that covered only two arms; pigauto's conformal intervals do. Scoped the same day (bee0e87) | drift caught in-arc | Rose |
| public claims | the results Artifact "shown" to Shinichi | v1 and v2 published private; not visually inspected by the agent (browser not signed in; sign-in declined by rule). JS parse-checked, ids verified. Shinichi has not yet confirmed reading it | unclear | Shinichi |
| model routing | Fable orchestrates; Sonnet builders; Opus for load-bearing review | Shinichi corrected mid-arc ("do not use Fable for parallel work"; "you use Fable as appropriate and orchestrating jobs"); refuters pinned to Sonnet, critic to Opus | adaptive | Ada |
| safety gates | D-139 pre-run before any run over 30 min | held for the campaign; the `freq_lambda` waves (12 min, 7 min, 12 min) were under the line and launched on precedent | adaptive | Ada |
| handoff state | overnight authority: push and draft PR | draft PR #184 was opened 2026-09-20 23:52 UTC, so the 09-20 row above was written minutes before it happened and is stale, not a gap; the remote branch had then fallen 46 commits behind, and was brought current today with the PR body refreshed to the present state | adaptive (09-20 row corrected) | Rose |

```
DECISION RECEIPT
  Questions asked      — Q6 thresholds, Q7 factorial trim, Q8 primary contrast (Phase 0.4, 09-20);
                         G0 budget approval (after the pre-run, 09-20); n = 1000 BACE replication
                         (three costed options, 09-21); fifth arm yes/no (09-21); core-now vs
                         factorial-later for freq_lambda (09-21)
  Answers received     — "Go ahead" (G0); "option 2, reduce BACE to 30 seeds at n = 1000 in the
                         FACTORIAL"; "Keep model = BM only and disclose" then "OK changed my mind. Yes,
                         at the lambda stuff as the fifth arm."; "Core slice now, factorial after";
                         "totoro you can use up to 250 for snakagaw OK"
  Defaults accepted    — arms 3a and 3b both run (Q1 default); covariate sensitivity after the main
                         runs (Q2); phyloglm Poisson GEE for counts (Q3); fir substituted for the
                         retired rorqual without asking (reasoned, logged, reversible)
  Adaptive decisions   — divergence rule; per-host pool subdirectories; BACE convergence by BACE's own
                         verdict; freq_lambda run on AVONET and covsens for arm parity; the covsens
                         n = 1000 array cancelled and resubmitted when stuck
  Unresolved           — Shinichi's read of the Artifact and the six publication decisions on its last
                         tab (G13c); whether the pkgdown article goes public, unlisted, or stays
                         internal; Szymek's sign-off on the corrected Pagel-lambda form; the arm-3
                         default solver as a package change -> Rose
```
