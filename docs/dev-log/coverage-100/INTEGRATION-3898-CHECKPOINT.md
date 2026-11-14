# Coverage qualification checkpoint

Frozen candidate `3898a77b1a953341a414b377c16a3e0e4e95ba5c` covers **12,878 of 13,003 executable locations (99.038683%)**, leaving **125 uncovered**. The percentage describes this combined candidate, not main or the coverage-only PR head. Independent Sol review found no artifact blocker for this partial qualification checkpoint.

The preceding a89152a3 snapshot covered 12,892 of 13,028 locations with 136 uncovered. The net denominator decreases by 25, the covered count decreases by 14 and the uncovered count decreases by 11. The new percentage therefore reflects reviewed guard cleanup, not an increase in covered locations. No coverage exclusion was added. The separate repair approvals are still outstanding.

| Original 1,904 misses | Locations |
|---|---:|
| Unchanged mapped locations covered by named test sets | 1,754 |
| Unchanged mapped locations uncovered | 125 |
| Changed or removed locations pending approval | 25 |
| Total | 1,904 |

The original unresolved remainder is **150**, compared with 148 previously. Two previously covered original entries, R/draws_conditional.R lines 240 and 241 in the baseline source, now fall within changed hunks and are conservatively held pending approval. Lower current misses do not erase that approval requirement. For the prior 13,002-location catalog, 12,809 are named-covered, 125 uncovered and 68 changed. All 69 new current locations are covered: 13,002 minus 68 plus 69 equals 13,003. No unchanged previously covered mapped location regresses. All 253 neural-method locations are covered; current coverage has 12,866 single-test witnesses and 12 collective witnesses.

Independent review checked every catalog identity, mapped source text and hunk, hit counter, named witness and new-location complement. It independently compared the prior snapshot: eleven old misses became changed entries pending approval, two previously named-covered original entries became pending, and no old miss gained coverage in this snapshot. Raw RDS bytes were not independently downloaded.

The instrumented suite passes **5,426 assertions**, with 192 inventoried warnings and 33 classified skips. The installed suite passes **5,423 assertions**, with 192 warnings and 34 skips. The extra installed NEWS skip accounts for three assertions executed in the instrumented run. Package check reports zero errors, warnings and notes, with rebuilt vignettes and no PDF manual. The separate source gate passes 192 assertions without conditions. All condition rows reconcile. Missing optional fixtures and HTML Tidy remain disclosed.

Collection took **19 minutes 31 seconds**, using one CPU thread and torch 0.17.0. The runtime receipt records CPU tensor execution, unavailable CUDA/MPS and available zip, unzip and Pandoc. An SSH interruption delayed retrieval; the collector's terminal successful manifest and artifacts were recovered. Raw RDS remains remote with its byte hash in the named-evidence receipt. Strict verification continues to reject modified-source identity. Conservative reconciliation preserves changed original locations and always reports incomplete.

The source includes repairs through #298 and #300, plus the ambient-seed test. Later candidate `adad2c1d9e9542339fc4f5fd31f70c5b3e1f7c6f` adds the separately reviewed cache-key cleanup #301, joint-refinement rollback repair #302 and additional tests. It passes 1,716 focused assertions across 24 changed files; its full coverage is not yet measured. The rollback repair fixes a pre-existing opt-in path returning initialization after reporting one accepted refinement. A same-test exact-source comparison has six expected failures before the repair and 26 passes afterwards, including an independent dense oracle and public baseline dispatch.

The later metadata test exceeded the uncapped focused run's 120-second gate before reaching that new case. A smaller eight-tip fixture passes nine standalone assertions. The full focused prediction file then passes 406 assertions with explicit one-thread limits. Its executable gate includes those limits and passes reverification; no branch-hit claim follows from a backend-tolerant pass alone.

G7 remains unmet, and final G8 awaits the final candidate. Coverage does not establish calibration, GPU correctness, inferential validity or singleton resolution. Draft PRs remain unmerged, with no public issue update or CRAN release action.
