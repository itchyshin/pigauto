# Comprehensive validation

OWNS: script/validation_20261018/**, tests/testthat/test-validation-contract.R, docs/dev-log/validation-20261018/**, LOOP/lanes/validation-20261018/**

Scope: approved scientific validation implementation, Totoro production, scientific assessment, and reproducible closeout.

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

- [x] G5: Production scientific assessment and attempt accounting
  CWD: /Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto
  CHECK: Rscript docs/dev-log/validation-20261018/test_posthoc_summarise.R; full post-hoc evaluator on the Totoro result tree; independent result-key and classification review
  EXPECT: 32,240 unique jobs accounted; source/job hashes and per-record runtime identities verified; every result key reconciled to the summary; inferential classifications reported by regime
  EVIDENCE: 32,199 ok, 1 failed, 40 unavailable; 1,122,060 result rows; 5,697 result keys; 12 failed, 32 indeterminate and 0 supported inferential rows. Runtime inventory variation is retained and reviewed. Scientific failures are findings, not passes.

- [x] G6: Independent review and reconciliation
  EVIDENCE: Independent reviewer matched all job IDs, seeds and scenarios; independently recomputed result key counts and assessment classifications; confirmed 12 failed, 32 indeterminate and zero supported inferential rows. It identified the matrix inventory-count defect, which was fixed and retested. Residual: the report retains all six environment hashes; only `posterior_001` spans multiple hashes, with one loaded-session fingerprint.


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

- [x] G9: Expected method/estimand roster and payload validation
  CWD: /Users/z3437171/.codex/worktrees/pigauto-comprehensive-validation/pigauto
  CHECK: Rscript script/validation_20261018/test_roster.R; bash script/validation_20261018/verify_bien_roster_totoro.sh; Rscript script/validation_20261018/audit_posterior_denominators.R; Rscript script/validation_20261018/audit_roster_receipts.R
  EXPECT: Every planned method/estimand key is present or explicitly unavailable; missing, partial and unexpected keys remain visible; probability bins and posterior interval denominators are checked.
  EVIDENCE: Local negative controls reject negative counts, out-of-range calibration values, invalid statuses, impossible coverage/width and tampered denominators. A fresh read-only Totoro derivation matches all five BIEN trait names, source-cache hashes, 6,425 matched species, subset size and seed. Independent replay of 4,000 posterior DGP replicates matches all 8,000 x/y interval denominators; the corrupted-denominator control fails as intended. Receipt audit: 5,697 method/estimand groups and 1,122,060 planned job-level keys all present; zero absent, partial, unexpected or unrostered rows. Of 504 metric groups, 16,400 ordinal GNN-off calibration-bin records are unavailable; zero groups are missing or failed. Three result rows failed and 246,040 were unavailable. Independent re-review passed the status validation, account target, and proof-to-roster binding.
