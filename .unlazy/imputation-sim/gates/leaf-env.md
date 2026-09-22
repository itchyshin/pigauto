# Gates: environments (Mac, Totoro, nibi)

OWNS: (no repo files; R libraries on three hosts)

Scope: pigauto 0.11 (gnn argument) plus castor, Rphylopars, BACE, MCMCglmm, missForest, phylolm, torch installed and smoke-proven on Mac, Totoro and nibi.

- [x] G1a: Mac has pigauto with the gnn formal
  CHECK: Rscript -e 'stopifnot("gnn" %in% names(formals(pigauto::impute))); cat("G1a PASS\n")'
  EXPECT: G1a PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=8b01951b15dd/30 entries; output=G1a PASS

- [x] G1b: Totoro has pigauto with the gnn formal plus castor and missForest
  CHECK: ssh -o BatchMode=yes -o ConnectTimeout=12 snakagaw@totoro.biology.ualberta.ca 'Rscript -e "stopifnot(\"gnn\" %in% names(formals(pigauto::impute)), requireNamespace(\"castor\"), requireNamespace(\"missForest\")); cat(\"G1b PASS\n\")"'
  EXPECT: G1b PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=8b01951b15dd/30 entries; output=Loading required namespace: castor | Loading required namespace: missForest

- [x] G1c: nibi has pigauto with the gnn formal in the /project library under the recorded module set
  CHECK: ssh -o BatchMode=yes -o ConnectTimeout=12 snakagaw@nibi.alliancecan.ca 'source ~/projects/def-snakagaw/snakagaw/pigauto_sim/env.sh && Rscript -e "stopifnot(\"gnn\" %in% names(formals(pigauto::impute)), requireNamespace(\"BACE\"), requireNamespace(\"castor\")); cat(\"G1c PASS\n\")"'
  EXPECT: G1c PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=8b01951b15dd/30 entries; output=Loading required namespace: BACE | Loading required namespace: castor

- [x] G7: Totoro smoke cell ran every arm
  CHECK: ssh -o BatchMode=yes -o ConnectTimeout=12 snakagaw@totoro.biology.ualberta.ca 'grep -c "done in" ~/pigauto_sim/smoke/smoke.log && ! grep -q ERROR ~/pigauto_sim/smoke/smoke.log && echo G7 PASS'
  EXPECT: G7 PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=8b01951b15dd/30 entries; output=6 | G7 PASS

- [x] G8: nibi smoke cell ran every arm on a COMPUTE node via salloc (hostname recorded, not l5/login)
  CHECK: ssh -o BatchMode=yes -o ConnectTimeout=12 snakagaw@nibi.alliancecan.ca 'f=~/projects/def-snakagaw/snakagaw/pigauto_sim/smoke/smoke.log; grep -c "done in" $f && ! grep -q ERROR $f && grep -q "^host=" $f && ! grep -qE "^host=(l[0-9]|login)" $f && echo G8 PASS'
  EXPECT: G8 PASS
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim/.unlazy/imputation-sim/gates; path=8b01951b15dd/30 entries; output=6 | G8 PASS
