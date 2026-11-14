"""Copy actual release proof bytes into build-excluded evidence; no checks or network."""
from pathlib import Path
import hashlib,json,re,shutil
S=Path(__file__).resolve().parent
R=Path('/private/tmp/pigauto-pr236-validator-20261017')
OUT=R/'docs/dev-log/cran-0.11-audit/provenance/final-a93fece7-2026-10-18'
F=json.loads((S/'FINAL.json').read_text());SHA=F['sha256']
assert SHA=='a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb'
assert hashlib.sha256(Path(F['path']).read_bytes()).hexdigest()==SHA
assert not OUT.exists(),'Refuse overwriting a publication snapshot'
OUT.mkdir()
rows=[]
def cp(p,rel):
 assert p.is_file() and not p.is_symlink(),p
 b=p.read_bytes()
 assert not re.search(rb'(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9]{32,}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----)',b),'Credential-shaped material'
 q=OUT/rel;q.parent.mkdir(parents=True,exist_ok=True)
 assert not q.exists(),q
 q.write_bytes(b);assert q.read_bytes()==b
 rows.append({'source':str(p),'copy':str(rel),'sha256':hashlib.sha256(b).hexdigest(),'size_bytes':len(b)})
def tree(p,rel):
 assert p.is_dir(),p
 for f in sorted(p.rglob('*')):
  if f.is_dir():continue
  assert not f.is_symlink(),f
  if '__pycache__' in f.parts:continue
  cp(f,Path(rel)/f.relative_to(p))
for name in ['FINAL.json','main-ci-final.json','direct-run.json','direct-artifacts.json','direct-jobs.json','source-workflow-size-binding.json','synthetic-workflow-commit.json','synthetic-workflow-content.json','native-platform-summary.json','summarise-native.py','summarise-native-before-duplicate-summary-guard.py','release-validator-selftest.log','cran-publication-and-revdeps.json','rhub-object-readback.json','rhub-upload-execution.json','checker-retarget-receipt.json','mi-matrix-supervisor.json','run-mi-successor.py','run-local-successor.py']:
 cp(S/name,name)
for name in ['freeze-review','checker-review','native-summary-review','site','policy-refresh','installed-novice','vignette-timing','actual-warning-audit']:
 tree(S/name,name)
# The fresh matrix's libraries are dependencies, not proof files to republish.
cp(S/'mi-matrix/manifest.json','mi/manifest.json')
for p in sorted((S/'mi-matrix').glob('*.log')):cp(p,Path('mi')/p.name)
for label,root in [('macos-arm64',S/'native-macos')]+[(p.name.removeprefix('frozen-a93fece7-'),p) for p in sorted((S/'direct-platform').glob('frozen-*')) if p.is_dir()]:
 tree(root/'platform-evidence',Path('native')/label/'platform-evidence')
 rc=root/'frozen-check/pigauto.Rcheck'
 for p in sorted(rc.iterdir()):
  if p.is_file():cp(p,Path('native')/label/'Rcheck'/p.name)
 cp(rc/'tests/testthat.Rout',Path('native')/label/'Rcheck/tests/testthat.Rout')
cp(S/'native-macos/native-receipt.json','native/macos-arm64/native-receipt.json')
cp(S/'direct-platform/download-binding.json','native/server-digest-binding.json')
cp(S/'direct-platform/verification.log','native/direct-verification.log')
# Ancillary source/install libraries are excluded; retain outputs, helpers and actual reviews.
for p in sorted((S/'ancillary').iterdir()):
 if p.is_file():cp(p,Path('ancillary')/p.name)
tree(S/'ancillary/spelling-review','ancillary/spelling-review')
for p in sorted((S/'warning-replay').iterdir()):
 if p.is_file():cp(p,Path('warning-replay')/p.name)
tree(S/'warning-replay/attempts','warning-replay/attempts')
cp(S/'replay-provenance-observed-discrepancy.json','warning-replay/root-observed-discrepancy.json')
assert len(rows)>100
(OUT/'publication-manifest.json').write_text(json.dumps({'artifact_sha256':SHA,'source_commit':F['source_commit'],'records':rows,'scope':'Selected actual proof bytes; no installed libraries, working source trees or predecessor votes relabeled.'},indent=2)+'\n')
for row in rows:
 p=OUT/row['copy'];assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'] and p.stat().st_size==row['size_bytes']
(S/'publication-root.txt').write_text(str(OUT)+'\n')
print('SUCCESSOR_ACTUAL_PROOF_BYTE_READBACK_PASS',len(rows),sum(r['size_bytes'] for r in rows))
