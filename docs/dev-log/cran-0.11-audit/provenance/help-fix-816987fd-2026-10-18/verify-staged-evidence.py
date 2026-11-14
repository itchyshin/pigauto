from pathlib import Path
import subprocess,hashlib,json
R=Path('/private/tmp/pigauto-pr236-validator-20261017');G='/Library/Developer/CommandLineTools/usr/bin/git'
def git(*args):return subprocess.check_output([G,*args],cwd=R)
assert git('rev-parse','HEAD').decode().strip()=='6cd5f65b06d0c3ee32ce4b2d5a67a22c06028116'
assert git('branch','--show-current').decode().strip()=='release/cran-0.11-gate-780c084'
subprocess.run([G,'diff','--exit-code','8e275ca738129c8d608f06ffc407b8141c1a2636','--','DESCRIPTION','NAMESPACE','R','man','inst','data','vignettes','README.md','LICENSE','NEWS.md','tests'],cwd=R,check=True)
paths=[p for p in git('diff','--cached','--name-only','-z').decode().split('\0') if p]
prefixes=['docs/dev-log/cran-0.11-audit/provenance/final-14ee7e5d-2026-10-17/','docs/dev-log/cran-0.11-audit/provenance/help-fix-816987fd-2026-10-18/']
assert all(p in ['cran-comments.md','docs/dev-log/cran-0.11-audit/GATES.md','docs/dev-log/cran-0.11-audit/release-ledger.json'] or any(p.startswith(x) for x in prefixes) for p in paths)
child=subprocess.Popen([G,'cat-file','--batch'],cwd=R,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
try:
 for name in paths:
  child.stdin.write((':'+name+'\n').encode());child.stdin.flush()
  header=child.stdout.readline().split();assert len(header)==3 and header[1]==b'blob',(name,header)
  data=child.stdout.read(int(header[2]));assert child.stdout.read(1)==b'\n'
  assert data==(R/name).read_bytes(),f'Staged normalization/drift {name}'
 child.stdin.close();assert child.wait(timeout=10)==0
finally:
 if child.poll() is None:child.terminate();child.wait(timeout=10)
print('STAGED_EVIDENCE_ONLY_PACKAGE_BYTES_UNCHANGED_AND_ALL_BLOBS_VERIFIED',len(paths))
