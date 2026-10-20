from pathlib import Path
import hashlib,json,re,gzip,tarfile,io
R=Path('/private/tmp/pigauto-pr236-validator-20261017/docs/dev-log/cran-0.11-audit/provenance')
roots=[R/'final-14ee7e5d-2026-10-17',R/'help-fix-816987fd-2026-10-18']
patterns=[rb'gh[pousr]_[A-Za-z0-9]{36,}',rb'github_pat_[A-Za-z0-9_]{25,}',rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',rb'X-Amz-Signature=[a-f0-9]{32,}',rb'AIza[0-9A-Za-z_-]{30,}']
files=0;objects=0;size=0;hits=[]
def scan(b,name):
 global objects
 objects+=1
 if any(re.search(p,b) for p in patterns):hits.append(name)
for root in roots:
 for p in sorted(root.rglob('*')):
  assert not p.is_symlink(),p
  if not p.is_file() or p.name in ['files.sha256','publication-scan.json']:continue
  b=p.read_bytes();files+=1;size+=len(b);scan(b,str(p.relative_to(R)))
  if p.name.endswith('.tar.gz'):
   with tarfile.open(fileobj=io.BytesIO(b),mode='r:gz') as tf:
    for m in tf:
     if m.isfile():scan(tf.extractfile(m).read(),str(p.relative_to(R))+'::'+m.name)
  elif p.suffix=='.gz':
   try:scan(gzip.decompress(b),str(p.relative_to(R))+'::gzip')
   except (OSError,EOFError):pass
assert not hits,'Potential credential pattern detected; do not publish'
x={'scope':'All current predecessor/help packet regular files and retained tar/gzip member contents; pattern-based scan only, not general privacy certification','files_scanned':files,'objects_scanned_including_unpacked':objects,'file_bytes':size,'credential_pattern_hits':len(hits),'symlinks':0}
(roots[1]/'publication-scan.json').write_text(json.dumps(x,indent=2)+'\n')
root=roots[0]
(root/'files.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(root))+'\n' for p in sorted(root.rglob('*')) if p.is_file() and p.name!='files.sha256'))
print('PUBLICATION_PATTERN_SCAN_PASS',json.dumps(x))
print('PREDECESSOR_MANIFEST_WRITTEN')
