# Coverage qualification checkpoint

Frozen candidate `ae15cda0816c84d549b2aa623014c1fa4c26e5db` covers **12,897 of 13,000 executable locations (99.207692%)**, leaving **103 uncovered**. Independent Sol review verified the source-bound counters, witnesses and accounting. The result describes this combined candidate, separate from main and the coverage-only PR head.

Compared with 848a595b, five old misses gained named-test coverage: two floating-point fallback paths and three CUDA cleanup control-flow paths. The latter use explicit external and routing mocks while real scalar probes, fitting, predictions and losses run on CPU. They do not validate CUDA execution. The finite-variance repair #303 also adds one new executable location, covered by its regression. No old location was removed or newly held pending approval in this step.

| Original 1,904 misses | Locations |
|---|---:|
| Unchanged mapped locations covered by named test sets | 1,774 |
| Unchanged mapped locations uncovered | 103 |
| Changed or removed locations pending approval | 27 |
| Total | 1,904 |

The original unresolved remainder is **130**. For the prior 13,002-location catalog, 12,827 are named-covered, 103 uncovered and 72 changed. All 70 new current locations are covered: 13,002 minus 72 plus 70 equals 13,000. There are zero strict or mapped coverage regressions, and all 253 neural-method locations remain covered. Current coverage has 12,884 single-test witnesses and 13 collective witnesses.

The instrumented suite passes **5,591 assertions**, with 192 inventoried warnings and 33 classified skips. The installed suite passes **5,588 assertions**, with 192 warnings and 34 skips; NEWS accounts for the three-assertion difference. Package check reports zero errors, warnings and notes, with rebuilt vignettes and no PDF manual. The source-only gate passes 192 assertions without conditions. Collection took **18 minutes 7 seconds** on one CPU thread. The runtime receipt records torch 0.17.0, real CPU execution, unavailable CUDA/MPS and available zip, unzip and Pandoc.

Independent review recomputed all catalog identities, source hunks and unchanged text, counters, named witnesses and the new-location complement. It independently parsed both raw Rout logs and matched all condition rows. The remote raw RDS hash is receipt-bound; its bytes were not independently downloaded. Strict baseline-source identity remains unmet for modified files, and conservative reconciliation holds changed entries pending separate repair approval.

All 103 current misses retain source-bound triage: 71 failure/fallback, 14 dependency/platform, 11 suspected unreachable, four reachable and three legacy rows. Classification is a guide, not proof of unreachability. Numerical safeguards remain in place. No coverage exclusion was added.

Later final-mean and isolated-resource tests are outside this measured source. The resource prototype reached a real allocation-failure guard, but its R vector-heap limit is distinct from total process memory. Its owned-child admission, thread, portability and provenance checks require independent review before integration and full-profile credit. The final-mean fixture has a hand-calculated healthy control and refuses an unrepresentable double-valued conditional mean. Neither prototype establishes public-fit reliability or statistical calibration.

G7 and final G8 remain unmet. All PRs remain unmerged; no public issue update or CRAN release action has been taken.
