"""Read-only exact-successor closure checks; never builds, fits, pushes or submits."""
from pathlib import Path
import argparse,copy,hashlib,json,subprocess
S=Path(__file__).resolve().parent;R=Path('/private/tmp/pigauto-pr236-validator-20261017');O=R/'docs/dev-log/cran-0.11-audit/provenance/final-a93fece7-2026-10-18'
SHA='a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb';SOURCE='b2e6054685e6bdedb073ba46e058af2c0a619a18';GIT='/Library/Developer/CommandLineTools/usr/bin/git'
def read(p):return json.loads(p.read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(ok,why):
 if not ok:raise ValueError(why)
def pre(x):
 a=x['artifact'];require(a['current_candidate'] is True,'withdrawn candidate');require(a['sha256']==SHA and a['source_commit']==SOURCE,'wrong candidate/source')
 p=Path(a['path']);require(p.is_absolute() and p.is_file() and not p.is_symlink(),'missing selected archive');require(digest(p)==SHA and p.stat().st_size==a['size_bytes']==5137837 and not p.stat().st_mode&0o222,'archive differs/writable')
 require(a['archive_entries']==250 and a['forbidden_entries']==[] and a['clean_worktree'] is True,'freeze scope invalid')
 f=read(O/'FINAL.json');require(f['sha256']==SHA and f['source_commit']==SOURCE,'freeze receipt mismatch')
 require(digest(Path(f['path']))==SHA,'original checked archive changed')
 for row in read(O/'publication-manifest.json')['records']:
  q=O/row['copy'];require(q.is_file() and not q.is_symlink() and digest(q)==row['sha256'] and q.stat().st_size==row['size_bytes'],'copied proof changed '+row['copy'])
 n=read(O/'native-platform-summary.json');require(n['artifact_sha256']==SHA and len(n['platforms'])==4,'native matrix incomplete')
 for v in n['platforms']:
  require(v['status']=='OK' and v['artifact_sha256']==SHA and v['declared_dependencies_available']==29 and v['native_counts']=={'fail':0,'warnings':161,'skips':86,'pass':3195},'native evidence differs')
 w=read(O/'actual-warning-audit/verdict.json');require(w['verdict']=='PASS' and w['archive_sha256']==SHA and w['passed_cells']==4,'actual warning audit incomplete')
 require((O/'actual-warning-audit/timing-provenance-amendment/amendment.json').is_file(),'replay provenance amendment missing')
 for name,key in [('installed-novice/verdict.json','verdict'),('vignette-timing/verdict.json','verdict')]:
  v=read(O/name);require(v[key]=='PASS','installed/vignette result incomplete')
 require(read(O/'mi/manifest.json')['status']=='PASS','MI notPASS')
 m=read(O/'mi/manifest.json');require(m['archive_sha256']==SHA and len(m['cells'])==4 and sum(len(c['steps']) for c in m['cells'])==12,'MI identities/routes incomplete')
 require(all(t['exit_code']==0 and t['archive_unchanged'] is True for c in m['cells'] for t in c['steps']),'MI step failed')
 require(read(O/'site/verdict.json')['verdict']=='PASS','live site notverified')
 require(x['source_candidate']['deployed_site_gate']=='PASS','deployed site currentrecord mismatch')
 require(x['submission']['cran_upload_authorized'] is True and x['submission']['maintainer_confirmation_authorized'] is True,'conditional user authority missing')
 require(x['evidence_pr']['number']==236 and x['evidence_pr']['merged'] is False,'evidence PRmerge prohibited')
 comments=(R/'cran-comments.md').read_text();require(SHA in comments and SOURCE in comments and '403' in comments and '14ee7e5d' not in comments and 'No CRAN upload' not in comments,'submission comments stale')
 subprocess.run([GIT,'-C',str(R),'diff','--exit-code',SOURCE,'--','DESCRIPTION','NAMESPACE','R','man','inst','data','vignettes','README.md','LICENSE'],check=True,stdout=subprocess.DEVNULL)
 return True
def panel(x):
 pre(x)
 for name in ['grace','rose','pat']:
  v=read(O/'panel'/name/'verdict.json');require(v['reviewer']==name and v['verdict']=='READY' and v['artifact_sha256']==SHA and v['source_commit']==SOURCE and v['blockers']==[],'final reviewer notREADY '+name)
  require(x['panel'].get(name)=='READY','current panel ledger differs '+name)
 require(x['status_claim'] in ['submission-ready','submitted','confirmed'],'not submissionready claim')
 return True
def selftest(x):
 pre(x);count=0
 for label,mut in [('withdrawn',lambda a:a['artifact'].update(current_candidate=False)),('oldhash',lambda a:a['artifact'].update(sha256='0'*64)),('oldsource',lambda a:a['artifact'].update(source_commit='8e275ca')),('wrongsize',lambda a:a['artifact'].update(size_bytes=1)),('forbidden',lambda a:a['artifact'].update(forbidden_entries=['BACE/'])),('notclean',lambda a:a['artifact'].update(clean_worktree=False)),('merge236',lambda a:a['evidence_pr'].update(merged=True)),('noauthority',lambda a:a['submission'].update(cran_upload_authorized=False))]:
  v=copy.deepcopy(x);mut(v)
  try:pre(v)
  except (ValueError,AssertionError):count+=1
  else:raise AssertionError('negative accepted '+label)
 require(count==8,'missing controls');print('SUCCESSOR_CONTROLS_PASS_8_REFUSALS')
def confirmation(x):
 panel(x)
 v=read(O/'submission/cran-confirmation-receipt.json')
 require(x['status_claim']=='confirmed' and x['submission']['state']=='confirmed','confirmation state not proven')
 require(v['state']=='confirmed' and v['uploaded_sha256']==SHA==x['submission']['uploaded_sha256'],'confirmed archive differs')
 require(v['uploaded_size_bytes']==5137837==x['submission']['uploaded_size_bytes'],'confirmed size differs')
 require(v['submission_id']=='359257'==x['submission']['submission_id'],'confirmation id differs')
 require(v['incoming_passed'] is False and v['accepted'] is False and v['live_version_verified'] is False,'external state overclaim')
 require('This submission has already been confirmed' in (O/'submission/cran-confirmation-readback.txt').read_text(),'primary confirmation readback absent')
 require((O/'submission/cran-confirmed.png').is_file(),'confirmation screenshot absent')
 return True

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['preparation','panel','selftest','confirmation']);args=p.parse_args();x=read(R/'docs/dev-log/cran-0.11-audit/release-ledger.json')
 if args.mode=='selftest':selftest(x)
 elif args.mode=='panel':panel(x);print('SUCCESSOR_FINAL_PANEL_AND_EXACT_PROOF_VERIFIED')
 elif args.mode=='confirmation':confirmation(x);print('SUCCESSOR_EXACT_ARCHIVE_SUBMITTED_AND_CONFIRMED')
 else:pre(x);print('SUCCESSOR_EXACT_PREPARATION_PROOF_VERIFIED')
