# Gates: results and deliverables (S7, S7a Artifact, S7b BACE methods note, S7c pkgdown article)

OWNS: script/campaign_sim_results/**, docs/dev-log/arc/*-simulation-results.md, docs/dev-log/arc/*-simulation-methods-bace.md, vignettes/articles/simulation-study.Rmd, .Rbuildignore, _pkgdown.yml

Scope: every reported number carries its regime and a finite MCSE; the primary contrast (BACE vs frequentist, z-RMSE and 95% coverage with interval score, core slice, pooled over types) is reported first; no superlatives.

- [ ] G12: results csv: every row has finite MCSE > 0 and filled regime columns (dgp, n, lambda, rho, mechanism, arm, reps); paired-difference rows vs the reference arm present
  CHECK: Rscript script/campaign_sim_checks.R --gate G12 --dir script/campaign_sim_results
  EXPECT: G12 PASS
  EVIDENCE: pending

- [ ] G13a: BACE methods note and pkgdown article contain no superlative
  CHECK: ! grep -iqE "\b(best|highest|definitive|state.of.the.art|unprecedented)\b" docs/dev-log/arc/*-simulation-methods-bace.md vignettes/articles/simulation-study.Rmd && echo G13a PASS
  EXPECT: G13a PASS
  EVIDENCE: pending

- [ ] G13b: shipped prose passes the slop check
  CHECK: python3 ~/shinichi-brain/tools/slop_check.py docs/dev-log/arc/*-simulation-methods-bace.md vignettes/articles/simulation-study.Rmd && echo G13b PASS
  EXPECT: G13b PASS
  EVIDENCE: pending

- [ ] G13c: Shinichi read the S7a Artifact and recorded what ships where (Dan's summary; article visibility)
  EVIDENCE: pending

- [ ] G13d: R CMD check unaffected: vignettes/articles is .Rbuildignore'd
  CHECK: grep -qE '^\^vignettes/articles\$$' .Rbuildignore && echo G13d PASS
  EXPECT: G13d PASS
  EVIDENCE: pending
