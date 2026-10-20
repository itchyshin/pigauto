"""Verify a source-merge approval packet, without asserting artifact/release readiness."""
import csv, hashlib, json, subprocess, sys
from pathlib import Path
R=Path('/private/tmp/pigauto-pr236-validator-20261017')
T=Path('/private/tmp/pigauto-final-closure-20261017')
W=Path('/Users/z3437171/.codex/worktrees/cran-011-final-artifact/pigauto')
H=R/'docs/dev-log/cran-0.11-audit/provenance/help-fix-816987fd-2026-10-18'
F=R/'docs/dev-log/cran-0.11-audit/provenance/final-14ee7e5d-2026-10-17'
HEAD='816987fd160c1d2c5c58b020b07357ff17c077e5'
BASE='8e275ca738129c8d608f06ffc407b8141c1a2636'
SHA='14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb'
GIT='/Library/Developer/CommandLineTools/usr/bin/git'
def read(p):return json.loads(p.read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(ok,message):
 if not ok:raise ValueError(message)
def validate_record(x):
 require(x['headRefOid']==HEAD,'source head mismatch')
 require(x['baseRefOid']==BASE,'base head drift requires review')
 require(x['number']==238 and x['state']=='OPEN' and x['isDraft'] is False,'source PR state')
 require(x['mergeStateStatus']=='CLEAN','source PR must be mergeable without bypass')
 require(x['mergedAt'] is None,'approval packet must not preempt merge')
 checks={z['name']:z for z in x['statusCheckRollup']}
 for name in ['ubuntu-latest (R release)','ubuntu-latest (R devel)','macos-latest (R release)']:
  require(checks[name]['status']=='COMPLETED' and checks[name]['conclusion']=='SUCCESS',f'incomplete/failed {name}')
 for name,z in checks.items():
  require(z['status']=='COMPLETED' and z['conclusion'] in ('SUCCESS','SKIPPED'),f'unsettled/unexpected {name}')
def controls():
 import copy
 x={'headRefOid':HEAD,'baseRefOid':BASE,'number':238,'state':'OPEN','isDraft':False,'mergeStateStatus':'CLEAN','mergedAt':None,'statusCheckRollup':[{'name':n,'status':'COMPLETED','conclusion':'SUCCESS'} for n in ['ubuntu-latest (R release)','ubuntu-latest (R devel)','macos-latest (R release)']]}
 validate_record(x)
 mutations=[lambda y:y.update(headRefOid='0'*40),lambda y:y.update(baseRefOid='0'*40),lambda y:y.update(isDraft=True),lambda y:y.update(mergedAt='2026-10-18'),lambda y:y.update(mergeStateStatus='UNSTABLE'),lambda y:y['statusCheckRollup'][0].update(conclusion='FAILURE'),lambda y:y['statusCheckRollup'][1].update(status='IN_PROGRESS')]
 for change in mutations:
  y=copy.deepcopy(x);change(y)
  try:validate_record(y)
  except ValueError:pass
  else:raise AssertionError('negative control accepted')
 print('SOURCE_APPROVAL_CONTROLS_1_POSITIVE_7_REFUSALS_PASS')
if '--selftest' in sys.argv:
 controls();raise SystemExit(0)
ledger=read(R/'docs/dev-log/cran-0.11-audit/release-ledger.json')
require(ledger['artifact']['current_candidate'] is False,'old artifact must stay withdrawn')
require(ledger['artifact']['sha256']==SHA,'predecessor identity drift')
require(ledger['submission']=={'state':'not_submitted','cran_upload_authorized':False,'maintainer_confirmation_authorized':False},'submission boundary')
require(ledger['evidence_pr']['number']==236 and ledger['evidence_pr']['merged'] is False,'evidence merge boundary')
require(ledger['source_help_correction']['successor_source_head']==HEAD and ledger['source_help_correction']['successor_archive_frozen'] is False,'successor source/artifact boundary')
require(ledger['source_help_correction']['merge_and_automatic_deployment_authorized'] is False,'merge approval boundary')
validate_record(read(H/'pr238-final.json'))
ci=read(H/'source-ci-final.json')
require(ci['headSha']==HEAD and ci['status']=='completed' and ci['conclusion']=='success','source CI identity/verdict')
require(len(ci['jobs'])==3 and all(j['conclusion']=='success' for j in ci['jobs']),'all three source jobs required')
require(subprocess.check_output([GIT,'rev-parse','HEAD'],cwd=W,text=True).strip()==HEAD,'source worktree head drift')
require(not subprocess.check_output([GIT,'status','--porcelain'],cwd=W,text=True),'source worktree must be clean')
paths=set(subprocess.check_output([GIT,'diff','--name-only',BASE,HEAD],cwd=W,text=True).splitlines())
require(paths=={'R/impute.R','R/report.R','man/impute.Rd','man/pigauto_report.Rd'},'source scope drift')
pat=read(H/'pat-review/verdict.json');require(pat['verdict']=='PASS','source reviewer refused')
for name,sha in pat['reviewed_file_sha256'].items():require(digest(W/name)==sha,'reviewed source hash drift')
rows=list(csv.DictReader((H/'runtime-AST-comparison.csv').open()))
require(len(rows)==49 and all(r['executable_AST_identical']=='TRUE' for r in rows),'runtime equivalence failed')
require(all(s['exit_code']==0 for s in read(H/'execution.json')),'provisional source/install failure')
n=read(H/'installed-workflow-receipt.json')
require(n['gnn_default'] is False and n['configured_gnn'] is False and n['input_missing']==2 and n['output_missing']==0 and n['observed_values_and_classes_preserved'] is True,'novice workflow')
site=read(H/'site-review/verdict.json')
require(site['verdict']=='PASS' and site['source_commit']==HEAD and site['source_pr']==238 and site['build_exit_code']==0,'fresh source site not passed')
for who in ['grace','rose','pat']:
 require(digest(F/'panel'/who/'verdict.json')==digest(T/'final-panel'/who/'verdict.json'),'original panel changed')
require(read(F/'panel/pat/verdict.json')['verdict']=='NOT_READY','old refusal must remain')
require(read(F/'panel/rose/verdict.json')['verdict']=='NOT_READY','old refusal must remain')
require(read(H/'checkpoint-independent/verdict.json')['verdict']=='PASS','checkpoint independent review')
for root in [H,F]:
 for line in (root/'files.sha256').read_text().splitlines():
  sha,name=line.split('  ',1);p=root/name
  require(p.is_file() and not p.is_symlink() and digest(p)==sha,f'published file hash drift {p}')
print('SOURCE_HELP_APPROVAL_CHECKPOINT_VERIFIED_RELEASE_NOT_READY')
