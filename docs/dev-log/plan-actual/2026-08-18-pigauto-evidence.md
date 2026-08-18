# Pigauto evidence programme: plan versus actual (2026-08-18)

| Planned slice | Actual state | Evidence / next gate |
|---|---|---|
| Stage A ADEMP harness | Implemented | Scripted DGP, common 60/20/20 split, active/random/uncertainty/no-acquisition policies, receipts, summaries, and integrity tests are in this branch. |
| Stage A local smoke | Complete, non-evidential | 30-tip, 5-epoch continuous and binary runs produced finite protected-test predictions and exactly one restoration. |
| Stage A Totoro pilot | Complete, non-claim-bearing | 800 receipts / 3,200 policy rows at source SHA `beee5df`; independent receipt audit passed, zero failed fits, and the process group was reaped. The pilot is mixed: continuous lambda=1 is directionally favourable; binary lambda=1 is not yet reliably favourable. |
| Stage A full campaign | Awaiting explicit approval | Every cell requires the registered 1,000-replicate minimum: 8,000 total / 7,200 additional replicates, projected at 1.16 / 1.04 hours respectively on 96 Totoro workers with 20% headroom. Retain all cells and errors. |
| Third-party prior art | Complete, audited | NotebookLM task `6d0404ae-9859-48ce-84c2-f190793421fe` was source discovery; four primary adjacent sources were checked. No novelty wording is authorised. |
| Stage B comparator | Blocked | PR #173 is open and unstable; its Ubuntu test commands reported 0 failures, but the workflow is not green. The comparator protocol is written; no parity/default claim is authorised. |
| Stage C Mondrian | Protocol only | FishBase and PanTHERIA masked-observation confirmation remains independent future work. |

The programme has cleared its pilot gate but is not scientifically complete:
the Stage-A headline remains unearned, Stage B awaits PR #173, and Stage C is
protocol-only.
