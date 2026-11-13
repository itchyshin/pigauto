# Prepared fixes, verified tests and exact remainder

The reviewed runtime combination `82ce4b0625a679026a6f9a9e0494aaa43da31ecf` passes full verification and records **12,905 of 13,002 executable locations covered (99.253961%)**. All 97 current misses are retained. Following Shinichi's approval, delivery finishes the prepared fixes and meaningful tests rather than pursuing 100%. The historical 100% gate remains unmet.

| Evidence | Covered / total | Scope |
|---|---:|---|
| Trusted starting Codecov artifact | 10,845 / 12,749 | Run42537328573, source033649b9 |
| Final CPU candidate | 12,905 / 13,002 | Frozen82ce4b06, expanded neural-method accounting |

The denominators reflect source changes and the additional neural-method inventory. The original 1,904 uncovered locations reconcile exactly: **1,780 covered with named tests, 95 still uncovered and 29 changed locations awaiting separate repair approval**. Thus 124 original entries remain unresolved. Two newly added defensive post-encoding guard lines are also uncovered, bringing current physical misses to 97. No exclusion was added to raise the percentage.

The complete 13,002-location catalog contains 12,827 named-covered entries, 95 misses and 80 changed entries. The current source adds 80 locations, of which 78 are covered and two are missed. This preserves the denominator accounting: 13,002 minus 80 plus 80 equals 13,002. Changed entries receive no new test-coverage credit. Compared with the preceding candidate, the positional covariate diagnostic gains a real witness; two formerly covered locations become approval-held through the count-validation repair.

Verification on the exact final source:

- Full instrumented suite: **5,740 passing assertions**, 192 inventoried warnings and 34 classified skips.
- Installed suite: **5,737 passing assertions**, 192 warnings and 35 skips. The source-only NEWS assertions explain the three-assertion difference.
- Package check: **zero errors, warnings or notes**, with rebuilt vignettes and the PDF manual excluded. This check-level result does not mean the test logs were warning-free.
- Source-only gate: **192 clean assertions**.
- All 31 changed test files: **1,960 clean assertions** locally, including 12 actual-MPS fixed-effect assertions on the same source.

Full CPU collection took 18 minutes 1 second, using one thread per job and a 30-minute cap. All suggested dependencies were available. The CPU host used torch0.17.0 without CUDA/MPS; its explicit MPS skip is classified. The separate actual-MPS test validates stored fixed-effect arithmetic and scalar bias. It is not aggregated into CPU coverage and does not establish general GPU support or statistical calibration.

Independent Sol and fresh bounded counter/oracle reviews confirm all 191 recorded source hashes, exact accounting, zero mapped coverage regressions and coverage of all 253 neural-method locations. There are 12,892 single-test witnesses and 13 collective witnesses. The raw RDS digest is retained through a remote hash receipt; the reviewers did not independently replay that full RDS.

The evidence bundle is `integration-82ce4b06/`, including raw logs, named evidence, complete old/new line ledgers, dependency versions, manifests and review receipts. Its residual CSV identifies every missed line; the triage CSV records owner, obstruction and needed evidence. Classification is not an impossibility proof. Reproduction commands are in `script/coverage_audit/README.md`; the raw RDS location and digest are in the bundle README.

`PR-INDEX.md` separates two coverage/evidence drafts from twenty repair drafts. These remain unmerged and require separate approval. The measured combined candidate is distinct from main and the coverage-only PR head. Scientific calibration, general GPU correctness, the singleton limitation and CRAN release remain separate. No issue update, tag, release or submission was made.

A later test-only CI portability follow-up supplies the shared-GNN fixture's reference tree explicitly and checks its identity. It reproduces and fixes a strict-warning failure when optional phangorn availability is controlled as absent, while preserving all numerical assertions and the warning whitelist. The focused MI suite passes 261 assertions. The original runner warning text was hidden, so the missing-dependency explanation is a supported reproduction rather than a confirmed runner diagnosis. Runtime source and coverage configuration remain unchanged; the counts above remain attributed to the exact82ce snapshot. Follow-up verification and current PR checks are recorded separately.
