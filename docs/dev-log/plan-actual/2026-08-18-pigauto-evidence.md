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
