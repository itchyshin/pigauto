# Comprehensive validation

OWNS: script/validation_20261018/**, tests/testthat/test-validation-contract.R, docs/dev-log/validation-20261018/**, LOOP/lanes/validation-20261018/**

Scope: approved scientific validation implementation and bounded pilot; production above3h needs measured-pilot approval.

- [x] G0: Source mutation and exact approval controls reject stale evidence
  CWD: /Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto
  CHECK: Rscript script/validation_20261018/test_orchestration.R && Rscript script/validation_20261018/test_approval.R
  EXPECT: APPROVAL_PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto; path=d124951413a0/39 entries; output=ORCHESTRATION_PASS | APPROVAL_PASS

- [x] G1: Manifest and evaluator reject invalid evidence
  CWD: /Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto
  CHECK: Rscript script/validation_20261018/test_harness.R
  EXPECT: HARNESS_PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto; path=d124951413a0/39 entries; output=HARNESS_PASS

- [x] G2: Eight-type public contract
  CWD: /Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto
  CHECK: Rscript script/validation_20261018/check_contract.R
  EXPECT: CONTRACT_PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto; path=d124951413a0/39 entries; output=Auto-detected trait type "count" (log1p-encoded) for integer column(s): count. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(count = "continuous"). | Auto-detected trait typ

- [x] G3: Campaign adapters and dispatch
  CWD: /Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto
  CHECK: Rscript script/validation_20261018/test_runners.R
  EXPECT: RUNNERS_PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto; path=d124951413a0/39 entries; output=RUNNERS_PASS

- [x] G4: Pilot accounting
  CWD: /Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto
  CHECK: Rscript script/validation_20261018/summarise.R --pilot
  EXPECT: PILOT_ACCOUNTED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto; path=d124951413a0/39 entries; output=RESULT_FILES 18 EXPECTED 18 | PILOT_ACCOUNTED

- [ ] G5: Production scientific assessment
  EVIDENCE: pending

- [x] G6: Independent review and reconciliation
  EVIDENCE: Separate reviewer verified both raw result sets, full settings, denominators and identity separation; reports retain production as unrun.


- [x] G7: Replicate interval summary controls
  CWD: /Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto
  CHECK: Rscript script/validation_20261018/test_summary.R
  EXPECT: SUMMARY_PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto; path=d124951413a0/39 entries; output=SUMMARY_PASS

- [x] G8: Corrected full-setting follow-up accounting
  CWD: /Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto
  CHECK: Rscript script/validation_20261018/summarise.R --pilot --out script/validation_20261018/results/supplement
  EXPECT: PILOT_ACCOUNTED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto; path=d124951413a0/39 entries; output=RESULT_FILES 17 EXPECTED 17 | PILOT_ACCOUNTED
