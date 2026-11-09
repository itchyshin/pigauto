# Coverage qualification checkpoint

Frozen candidate `848a595b5b6cf62df30c9d321ed566bec9e9aa15` covers **12,891 of 12,999 executable locations (99.169167%)**, leaving **108 uncovered**. Independent Sol review verified this combined-source result. The percentage is separate from main and the coverage-only PR head.

Compared with 3898a77b, **15 previously uncovered locations gained named-test coverage**. Two other misses moved to pending repair approval, and two previously covered catalog locations were removed in repair hunks. The net denominator decreases by four, the covered count increases by 13 and the missing count decreases by 17. The two removed misses are GBIF's empty-key fallback and the solver's stale rollback assignment. They receive no test-coverage credit.

| Original 1,904 misses | Locations |
|---|---:|
| Unchanged mapped locations covered by named test sets | 1,769 |
| Unchanged mapped locations uncovered | 108 |
| Changed or removed locations pending approval | 27 |
| Total | 1,904 |

The original unresolved remainder is **135**. For the prior 13,002-location catalog, 12,822 are named-covered, 108 uncovered and 72 changed. All 69 new current locations are covered: 13,002 minus 72 plus 69 equals 12,999. No unchanged previously covered mapped location regresses. All 253 neural-method locations remain covered. Current coverage has 12,878 single-test witnesses and 13 collective witnesses; the new collective witness is the ordinal CV summary branch.

Independent review checked every catalog identity, Git hunk and unchanged source text, counter value, named witness and new-location complement. It independently parsed both raw Rout logs and checked condition rows against the archived inventories. Raw RDS bytes remain remote, with their hash recorded in the source-bound receipt rather than independently downloaded.

The instrumented suite passes **5,518 assertions**, with 192 inventoried warnings and 33 classified skips. The installed suite passes **5,515 assertions**, with 192 warnings and 34 skips; NEWS accounts for the three-assertion difference. Package check reports zero errors, warnings and notes, with rebuilt vignettes and no PDF manual. The source-only gate passes 192 assertions without conditions. Collection took **21 minutes 48 seconds** on one CPU thread. The runtime receipt records torch 0.17.0, real CPU execution, unavailable CUDA/MPS, and available zip, unzip and Pandoc.

Every current missing location has source-bound triage in `residual-classification.csv`: 73 failure/fallback, 17 dependency/platform, 11 suspected unreachable, four reachable and three legacy rows. These labels guide the next work; they do not prove unreachability. Device-routing and CUDA cache mocks can test software cleanup control flow. Actual CUDA execution remains separate hardware evidence.

This snapshot includes the rollback repair #302 and cache-key cleanup #301. The later finite-variance repair #303 is outside it. That private kernel defect can return infinite variances under extreme tree scaling despite representable correlation-scale targets. Its separate regression fails before the repair, passes 25 assertions after it and preserves both healthy controls. New real floating-point fallback tests pass a 428-assertion focused prediction file; new mocked cleanup-boundary tests pass a 259-assertion GNN file. They carry no coverage credit in this snapshot. CPU boundary tests do not validate GPU execution or calibration.

G7 and final G8 remain unmet. No coverage exclusion was added. No PR was merged, public issue update posted or CRAN release action taken.
