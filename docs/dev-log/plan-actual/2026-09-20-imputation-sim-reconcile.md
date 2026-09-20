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
