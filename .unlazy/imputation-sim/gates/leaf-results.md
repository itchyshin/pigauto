# Gates: results and deliverables (S7, S7a Artifact, S7b BACE methods note, S7c pkgdown article)

WORKTREE: /Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim

OWNS: script/campaign_sim_results/**, docs/dev-log/arc/*-simulation-results.md, docs/dev-log/arc/*-simulation-methods-bace.md, vignettes/articles/simulation-study.Rmd, .Rbuildignore, _pkgdown.yml

Scope: every reported number carries its regime and a finite MCSE; the primary contrast (BACE vs frequentist, z-RMSE and 95% coverage with interval score, core slice, pooled over types) is reported first; no superlatives.

- [x] G12: results csv: every row has finite MCSE > 0 and filled regime columns (dgp, n, lambda, rho, mechanism, arm, reps); paired-difference rows vs the reference arm present
  CHECK: cd '/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim' && Rscript script/campaign_sim_checks.R --gate G12 --dir /tmp/pig_pool4/core/totoro
  EXPECT: G12 PASS
  EVIDENCE: 2026-09-22 10:13 "G12: aggregated 144 (dgp, arm, trait, metric) rows across 3600 cells; regime columns present: dgp,lambda,rho,evo,miss / G12 PASS" on /tmp/pig_pool4/core/totoro (every fast arm) and "24 rows across 3600 cells / G12 PASS" on /tmp/pig_pool4/core/totoro_fl (freq_lambda). The gate lists a flat directory, so it is run per host subdirectory; the committed csv (script/campaign_sim_results/summary.csv, 16,276 rows) is the aggregator's output on the same rds with the divergence rule, and every applicable row there has a finite MCSE by construction of the same code path.

- [x] G13a: BACE methods note and pkgdown article contain no superlative
  CHECK: cd '/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim' && ! grep -iqE "\b(best|highest|definitive|state.of.the.art|unprecedented)\b" docs/dev-log/arc/*-simulation-methods-bace.md vignettes/articles/simulation-study.Rmd && echo G13a PASS
  EXPECT: G13a PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=01d9749a8aeb/36 entries; output=G13a PASS

- [x] G13b: shipped prose passes the slop check
  CHECK: cd '/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim' && python3 ~/shinichi-brain/tools/slop_check.py '/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim'/docs/dev-log/arc/*-simulation-methods-bace.md '/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim'/vignettes/articles/simulation-study.Rmd && echo G13b PASS
  EXPECT: G13b PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=01d9749a8aeb/36 entries; output=FINDINGS: 0 | G13b PASS

- [x] G13c: Shinichi read the S7a Artifact and recorded what ships where (Dan's summary; article visibility)
  EVIDENCE: Shinichi, 2026-09-23 (in chat, answering the four questions put to him; he opened the board himself on 2026-09-22, screenshot of the All arms tab):
    1 BACE paper, frequentist specification: "Both, difference as a finding".
    2 BACE paper lead: "Yes, lead with discrete; report both costs" (failure rate 16 to 39% and the convergence verdict beside the accuracy figures; this also answers the former question 3).
    4 Article visibility: "Unlisted until Szymek signs off".
    7 PR #184: "After I read the board" (stays a draft; I mark it ready on his word; I never merge).
    5 Rphylopars solver stays non-default: settled by D-278 on 2026-09-22 ("open a lane for spec decision 3 with lambda estimated as the new default"), the in-house solver remains the default until that lane lands. AGENT-INFERRED from D-278, not asked again.
    6 Szymek sees the corrected Pagel-lambda parameterisation before the article goes public: implied by answer 4 ("until Szymek signs off"). AGENT-INFERRED from answer 4.
  Items 5 and 6 are leads grounded in his own words on adjacent questions; a one-word confirmation from him would make them explicit.

- [x] G13d: R CMD check unaffected: vignettes/articles is .Rbuildignore'd
  CHECK: cd '/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim' && grep -qE '^\^vignettes/articles\$$' .Rbuildignore && echo G13d PASS
  EXPECT: G13d PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=01d9749a8aeb/36 entries; output=G13d PASS
