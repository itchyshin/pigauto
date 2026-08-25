# pigauto 0.11 trust-and-usability product execution

## GOAL

```text
Solo platform: Codex
Deliverable: A source-identical local pigauto 0.11.0 trust-and-usability candidate.
HEADLINE: Make the scientifically safe path the easiest path; refuse known-invalid inferential use through pigauto's supported workflow.
IN PARALLEL: Read-only reconciliation, claim review, mechanical verification and community review only. Repository writing remains sequential.
DEFER: Exact-route default change; Stage C; BACE development; external competitiveness; broad interval certification; general active-sampling claims; push, PR, merge, release and CRAN submission.
DISCIPLINE: Prewritten acceptance ledger plus independent review; no new compute campaign; closure requires frozen source, tarball, tests, installed workflow, rendered documentation and plan-versus-actual receipt to agree.
```

## Approval and boundary

Shinichi approved the staged ultra-plan on 2026-08-25. Gate 4A source
provenance is already closed as a narrow pigauto-only subgate in the separate
`active-recovery` evidence worktree. The external-comparator gate remains open.
This product lane starts only after that evidence handover and does not import
its experimental results into package claims.

## Orientation receipt

- Worktree: `.worktrees/pigauto-0-11-trust-usability`
- Branch: `codex/pigauto-0-11-trust-usability`
- Base: refreshed `origin/main` at `7d775c06f1366d921f33add99f3bbe93dd0ca0b3`
- Baseline version: `0.10.0.9000`
- Lane: whole-repository Codex lease acquired; repository writers are sequential.
- Root checkout: unrelated and dirty; it is out of bounds.
- Prior work: current main already contains the narrow analysis-aware MI engine,
  result invariants, lower-level validators, read helpers and report machinery.
- Missing work: no reachable ref contains the requested diagnostic/tree MI
  refusals, shared `pigauto_check`, `completed_data()`, result summary contract,
  or canonical six-line novice journey.
- BACE history: `b615579` is a removal checklist only, not a cherry-pick; later
  GitHub development restored the bridge now being removed.
- Sister repositories: no package interface or evidence contract is reusable.
- Brain/repository reconciliation: the repository and the 2026-08-25 handover
  are the current technical authority; older project notes do not describe the
  requested 0.11 contracts.
- Search decision: no literature or novelty campaign is warranted because this
  lane authorizes no novelty, priority or competitiveness claim.
- Route decision: destination, interfaces, checks and output paths are known;
  no decision map or further user question is required.

### Clean baseline

Before product changes, the full suite ran with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` and completed in 600.2 seconds:
`FAIL 0 | WARN 336 | SKIP 8 | PASS 2180`. The warnings are the suite's
existing small-validation, zero-variance and bounded diagnostic warnings; no
baseline failure was carried into P1.

## Sequential slices

1. **P1 — inference boundary:** remove the installed BACE bridge; mark diagnostic
   and tree-sensitivity MI objects; admit only `pigauto_analysis_mi` through the
   supported workflow; propagate provenance into fit lists and enforce it when
   pooling.
2. **P2 — one diagnostics contract:** add fit-free `check_pigauto()`, persist it
   as `result$check`, add `completed_data()` and `summary.pigauto_result()`, and
   make the report consume the same contract without changing `$completed`.
3. **P3 — community surface:** teach the six-line workflow and bounded recipes;
   move solver/exact/calibration controls to an advanced section; sweep public
   claims to the evidence boundary.
4. **R1 — local candidate:** invoke the CRAN release gate, set the candidate
   version, generate documentation, run focused and full tests, build/install,
   execute the installed novice workflow, build pkgdown, run local
   `R CMD check --as-cran`, then freeze and inventory one tarball.
5. **V1/C1 — independent closure:** mechanical and adversarial review, rerun the
   acceptance ledger, reconcile plan versus actual, and write the after-task and
   local handover receipts. Any source repair invalidates the freeze and reruns
   candidate gates.

## Dispatch and ownership

One product builder owns the repository write surface across P1-P3 and is reused
sequentially. Read-only reviewers may inspect concurrently. Review agents never
write. This keeps every writing leaf disjoint in time and gives shared interface
files one active owner.

Model routing follows the approved plan: Terra/high for implementation, Luna for
mechanical reconciliation, and Sol/high for adversarial release and claim review.
No ultra effort is required.

## Verification contract

The prewritten local ledger is `.unlazy/pigauto-011/GATES.md`. Its commands are
not approved until every referenced script exists and has been read. The ledger
is ignored from Git; durable evidence is copied into the final candidate receipt.
No gate is silently removed. A source change after candidate freeze invalidates
the tarball and G8-G11.

## P1 receipt — accepted

P1 first went red in ten places against the permissive legacy behavior. The
initial implementation then failed adversarial review on mixed-class admission,
forged fit-list provenance, early warning placement and stale public guidance.
Those failures were reproduced before repair.

The accepted contract now requires both the `pigauto_mi_fits` class and
`mi_workflow = "pigauto_analysis_mi_v1"`. Diagnostic, tree-sensitivity,
legacy, unknown and inconsistent objects stop before callbacks or extractors.
Valid bare user fit lists remain an explicit escape hatch and emit exactly one
unverified-provenance warning only after a pooled result is successfully built.

The installed BACE bridge, dependency, export, help page and pkgdown reference
entry are removed. The in-tree BACE package and comparator helpers are
unchanged. Historical pooling directions in NEWS are explicitly superseded,
and current package and agent guidance names only the narrow
`multi_impute_analysis()` route for supported downstream pooling.

Final evidence:

- exact G5 command: `FAIL 0 | WARN 0 | SKIP 0 | PASS 46`;
- analysis-aware focused suite: 100 expectations, zero failures/warnings, one
  optional-package skip;
- diagnostic/tree generation suite: 200 expectations, zero failures and 32
  pre-existing small-fixture warnings;
- temporary installation: no BACE dependency, export, namespace binding or
  help entry;
- independent adversarial and mechanical re-reviews: pass;
- `git diff --check`: clean.

## P2 receipt — accepted

P2 adds a fit-free `check_pigauto()` contract, stores its compact summary on
`pigauto_result`, and gives users `completed_data()`, a result summary and a
report that share the same output-role language. Fully observed matched target
data now stop before graph construction or training and point to
`cross_validate()`.

Review exposed and repaired several neighboring correctness defects before the
slice was accepted: original row identities were lost around data-only species;
character predictions could be converted to factor codes; diagnostic draws used
an incomplete inverse row map; initial summary counts overlapped; and runtime
readiness reported accelerators without probing the selected device. The final
checker uses matched rows for type hints, emits parseable declarations even for
unusual names, preserves partial-species result invariants, and treats an
attempted but still-missing value as unresolved rather than filled.

Final evidence:

- exact G6 command: `FAIL 0 | WARN 0 | SKIP 0 | PASS 90`;
- exact G7 command: `FAIL 0 | WARN 0 | SKIP 0 | PASS 21`;
- neighboring `preprocess|shipping-coverage` suite: 106 expectations, zero
  failures and seven pre-existing small-fixture warnings;
- neighboring `new-features|multi-impute` suite: 445 expectations, zero
  failures, 148 pre-existing warnings and three expected skips;
- arbitrary backtick/backslash declaration actions independently round-trip;
- independent semantic and mechanical re-reviews: pass;
- documentation regeneration and `git diff --check`: clean.

## P3 receipt — accepted

P3 makes the supported path visible before advanced controls: the README and
getting-started vignette now teach one six-expression journey from files through
preflight, imputation, completed-data access and reporting. Installed novice
fixtures and an executable example preserve that journey for candidate testing.
Short recipes cover ordered states, integer measurements, proportions and
zero-inflated counts, compositions, repeated observations and species
reconciliation. Solver, EM, refinement, calibration and exact-prediction controls
remain available in a separately labelled advanced section; exact prediction
remains opt-in.

Claim review also repaired contradictory interval, tree-sensitivity, active-
sampling, BACE and historical PMM language. The retained tree benchmark path is
a descriptive-only tombstone, and current help distinguishes conditional
baseline standard deviations, diagnostic interval outputs and the one supported
analysis-aware MI route. Historical NEWS entries carry local withdrawal or
supersession boundaries so a deep link cannot revive an unsupported claim.

The claim gate normalizes Markdown whitespace and includes planted wrapped-text
negative controls. This closes a review-discovered bypass where ordinary regex
dot matching did not cross a line break.

Final evidence:

- exact `community-surface` suite: 79 expectations, zero failures, warnings or
  skips;
- exact six-expression installed example and all fixture paths parse and resolve;
- fit-free fixture preflight: `ready` for eight tree tips, three traits and two
  missing cells;
- README, vignette source and purl output are synchronized;
- seven changed roxygen source/help pairs regenerate byte-identically;
- pkgdown navigation, current tree article and descriptive-only tombstone agree;
- independent mechanical and adversarial re-reviews: pass, including planted
  wrapped BACE, PMM and tree-inference negatives;
- no `BACE/`, Stage C, default-route or exact-route implementation drift;
- `git diff --check`: clean.

## R1 pre-freeze receipt — accepted

The novice pre-run first exposed an invalid three-tip teaching fixture: the
default four spectral features require a larger tree. The fixture alone was
expanded to eight tips while retaining three traits and two missing cells. The
unchanged six-expression, default-argument CPU workflow then completed in 99.14
seconds and wrote its report. No package default was changed.

Reader-surface review found that unindexed `pkgdown/assets/dev` pages would
still be copied into a rendered site. Eleven stale comparator, coverage and
internal-test pages are therefore retained only as stable no-index tombstones;
their former contents remain recoverable from Git history. The current source
and rendered-site claim gate requires every exact tombstone and scans all HTML.

Candidate verification is one atomic, fail-closed run from a clean commit. It
remeasures HEAD and worktree status around every command; records exact
executables, arguments, directories, environments, timestamps and exit codes;
binds build, check, install, installed gates and rendered-site gates to one
tarball and staging root; and hashes the tarball plus hidden-aware installed-
library and site inventories. A partial freeze is never promoted to `CURRENT`.

Pre-freeze evidence:

- version set to local candidate `0.11.0`, with no submission or release claim;
- exact G5/G6/G7 focused gates remain 46/90/21 expectations;
- community-surface gate: 83 expectations, zero failures, warnings or skips;
- source G9 verifier and all six candidate scripts parse and pass light probes;
- component ledger enumerates every shipped data/extdata/example/logo/notice
  component and every asset-bearing excluded root;
- planted failures closed malformed dependency matching, wrapped prose, path
  spaces, stale logs, wrong source roots, different installed libraries,
  different site roots, clean HEAD switches, hidden-file drift and non-stopping
  source-test failures;
- independent mechanical and adversarial pre-freeze reviews: pass;
- `git diff --check`: clean.

## R1 candidate-attempt repair — accepted

The first atomic candidate attempt correctly stopped at the full source-test
gate and did not create `CURRENT` or a frozen candidate. Its 13 errors were all
one compatibility class introduced by the approved fully-observed-input stop:
older integration fixtures asked `impute()` to fit data with no missing target
cell. Two of those fixtures were invalid-argument tests whose expected
`clamp_factor` or `pmm_K` error had also been hidden behind the new preflight
stop.

The repair preserves the safety boundary. `impute()` now validates
`match_observed`, `clamp_factor` and `pmm_K` before preflight, while every active
integration fixture that genuinely tests fitting supplies at least one missing
matched target. A direct public test requires fully observed matched input to
stop and direct the user to `cross_validate()`.

Repair evidence:

- the failed predecessor remains only under ignored staging and was never
  promoted;
- all 13 mapped failures are repaired: two fail-fast precedence cases, eleven
  fitting fixtures, and the public stop contract;
- the expanded affected-neighbour suite covering preflight, clamp/PMM,
  fit/predict, joint solver, feature, safety-floor and zero-inflated-count paths
  exited zero; its 182 warnings are the pre-existing small-validation and
  zero-standard-deviation diagnostics exercised by those fixtures;
- documentation regeneration produced no generated-file drift;
- independent semantic and reconciliation reviews both accepted the repair,
  found no weakened assertion and found no active sibling fixture omitted;
- no default, exact-route, BACE, Stage-C or claim-scope drift;
- `git diff --check`: clean.
