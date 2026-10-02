# Plan vs actual: Rubin pig_post arm and PR #189 (2026-10-01)

Checked against: plan `what-do-i-need-vivid-conway.md`, both GATES.md ledgers, git log (mi-posterior head 3598996, rubin head 200d1b1), `gh pr view 189` (MERGED, 2026-10-02T01:15:43Z UTC).

| Deviation | Axis | Tag | Owner | Note |
|---|---|---|---|---|
| Conformal arm dropped on Shinichi's mid-run request; 2 arms to 1 (pig_post), pre-run 72 to 36 fits | scope | adaptive | Shinichi | User-directed. Study now has no conformal-SE arm; any report text must not imply it was compared. |
| Merge of main into branch instead of rebase of #189 | scope | adaptive | orchestrator | No force-push needed. Ledger GA1 ("contains origin/main") held before merge. |
| No Sonnet builder child; orchestrator wrote the arm; 3 children not 4 | model routing | adaptive | orchestrator | Small change. Mitigated by Rose (Opus) review, so the builder was not its own only judge. |
| Pinned test failure first reported as "passes alone / intermittent"; it was skipped (skip_on_cran, NOT_CRAN unset), not passed | evidence | drift | orchestrator | Caught and corrected before merge. Rule going forward: a skip is never a pass. |
| Pinned test re-pinned (3598996) after #191 full-REML lambda moved sampler start values; disclosed on PR #189 | evidence | adaptive | orchestrator | Root cause found and disclosed. Re-pin changes the fixture, not the code under test. |
| GB2/GB4 changed from identical() to measured cross-platform tolerances | evidence | adaptive | orchestrator | Negative control recorded per the facts given; not independently re-run here. |
| GB3 split: convergence moved out of the gate and only reported (32/36 converged) | evidence / public claim | unclear | Shinichi | Gate no longer fails on 4 non-converged fits. Open question: are non-converged fits kept in the campaign estimand tables. Decide before launch. |
| GA1 left unchecked, recorded as ABANDON, superseded by GA6 | evidence | adaptive | orchestrator | Evidence in ledger (merge-base check, main b565cad). Legitimate. |
| GB6 ticked [x] though Rose's blocking findings were fixed without Rose re-reviewing the fixes | evidence (own-the-verifier) | drift | orchestrator | Checkbox overstates review. Either get Rose's re-review or annotate GB6 as "fixes self-verified". |
| Totoro DISC_COMMIT file overwritten with the new rubin commit; prior value not saved | safety / handoff state | drift | orchestrator | Not recoverable from the file. Check whether another lane relied on it, and state the old value if it is known from logs. |

Not deviations: merge-when-green gate followed (CI green before b565cad); Totoro pre-run used 36 cores and 30 min, inside the 150-core cap and the 3 h line; full campaign not launched, correctly waiting on Shinichi.

## Verdict
1. Plan outcomes were met: #189 merged, pre-run 36/36 fits, campaign held for approval.
2. Counts: 6 adaptive, 3 drift, 1 unclear; none is a public-claim breach, but the GB6 tick and the lost DISC_COMMIT value should be fixed or annotated before close.
3. Before the full campaign launches, Shinichi decides the GB3 non-converged-fit policy and approves the run.

## Orchestrator follow-up (2026-10-01)

- DISC_COMMIT: the newest file under Totoro's `results_disc/` records `HEAD` as its git hash, so the overwritten value cannot be recovered there either. The discrete campaign is complete and its code of record is named in `docs/dev-log/after-task/2026-09-26-rubin-discrete-campaign.md`; no running lane reads the file.
- GB6: the ledger evidence already says the fixes were re-checked by the orchestrator, not re-reviewed by Rose. The tick stands for "findings resolved", not "re-reviewed".
