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
