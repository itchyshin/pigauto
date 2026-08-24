# Pigauto evidence programme: plan versus actual (2026-08-18)

| Planned slice | Actual state | Evidence / next gate |
|---|---|---|
| Stage A ADEMP harness | Implemented | Scripted DGP, common 60/20/20 split, active/random/uncertainty/no-acquisition policies, receipts, summaries, and integrity tests are in this branch. |
| Stage A local smoke | Complete, non-evidential | 30-tip, 5-epoch continuous and binary runs produced finite protected-test predictions and exactly one restoration. |
| Stage A Totoro pilot | Complete, non-claim-bearing | 800 receipts / 3,200 policy rows at source SHA `beee5df`; independent receipt audit passed, zero failed fits, and the process group was reaped. The pilot is mixed: continuous lambda=1 is directionally favourable; binary lambda=1 is not yet reliably favourable. |
| Stage A full campaign | Complete, headline gate fails | Audited 8,000 receipts / 32,000 policy rows at pilot SHA `beee5df` plus extension SHA `ad3990c`; zero policy failures and exact treatment/provenance checks passed. Continuous lambda=1 active-minus-random intervals are below zero at n=100 and n=300, but both binary lambda=1 intervals include zero. The broad headline is unearned; only the continuous lambda=1 recovery result is supported. |
| Third-party prior art | Complete, audited | NotebookLM task `6d0404ae-9859-48ce-84c2-f190793421fe` was source discovery; four primary adjacent sources were checked. No novelty wording is authorised. |
| Stage B comparator | Blocked on merge | PR #173 is open, clean, and all CI checks are green. The comparator protocol is written, but the study cannot start until the PR is merged; no parity/default claim is authorised. |
| Stage C Mondrian | Protocol only | FishBase and PanTHERIA masked-observation confirmation remains independent future work. |

The programme is not scientifically complete: Stage A earns only its narrow
continuous lambda=1 recovery result, not the cross-family headline; Stage B
awaits PR #173's merge; and Stage C is protocol-only.

## Reconciliation update (2026-08-20)

PR #173 has since merged to `main`.  The post-merge Stage-B comparator
protocol is now recorded in
`docs/dev-log/2026-08-20-stage-b-exact-comparator-protocol.md`; it specifies
two labelled AVONET regimes, five shared masks per regime, retained failures,
and the no-parity/default-change fence.  No comparator fit has yet been run,
so Stage B remains a protocol, not a competitiveness result.

Stage C now has a completed observed-only, two-dataset 300-tip operational
pilot and a completed PanTHERIA full-data timing receipt, but neither is
calibration evidence.  The first full FishBase CPU receipt was stopped after
24 hours and 19 minutes with only its mask receipt retained.  The resulting
scalability correction and the separately bounded Tamia GPU ladder are
recorded in `docs/dev-log/2026-08-19-mondrian-scalability-correction.md`.
Tamia job `419940` allocated four H100s but failed before a task ran because
of an inconsistent nested GRES request. Jobs `419946` and `419948` completed
the four bounded inputs with CUDA available, but their four independent
single-rank Slurm steps all executed on physical GPU 0; they are CUDA/input
smokes, not four-GPU scaling evidence. The corrected one-step, four-rank job
`424950` then completed in 24 seconds, with all four split receipts successful
and four distinct physical GPU UUID/PCI pairs. Its 7.030--9.650 second
per-input elapsed times are valid bounded operational ladder evidence, but
are not a CPU/GPU comparison, calibration evidence, or authority for a full
campaign.

## Reconciliation update (2026-08-23)

The corrected Stage-B one-mask timing smoke is complete and recorded in
`docs/dev-log/2026-08-23-stage-b-timing-smoke.md`. Both pigauto default and
opt-in exact arms completed in the continuous and mixed regimes, together
with the installed external comparators; their runtime is now measured rather
than guessed. BACE was explicitly retained as unavailable in both receipts.
Consequently this is a valid operational/timing receipt but not a complete
five-mask comparator study or any competitiveness evidence. The next Stage-B
gate is a concrete BACE source and private-library verification, followed by
a revised full-run estimate and explicit approval.

The programme therefore remains incomplete: Stage A has the narrow
continuous-only result and its independent claim review; Stage B awaits its
approved compute decision; and Stage C awaits a completed feasibility ladder,
a revised full-campaign cost estimate, explicit campaign approval, and an
independent calibration-claim review.
