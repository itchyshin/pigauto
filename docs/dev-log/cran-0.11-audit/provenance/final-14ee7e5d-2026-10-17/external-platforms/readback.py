#!/usr/bin/env python3
import pathlib,sys,json,hashlib,tarfile,gzip,io,runpy
p=pathlib.Path(sys.argv[1]);h=lambda b:hashlib.sha256(b).hexdigest()
r=json.loads((p/'publication-result.json').read_text());count=0
controls=runpy.run_path(str(p/'timing-parser-controls.py'))['run_controls'](p)
(p/'timing-parser-controls.json').write_text(json.dumps(controls,indent=2,sort_keys=True)+'\n')

def enforce(raw,sha,size):
    if h(raw)!=sha or len(raw)!=size:raise ValueError('Byte-binding mismatch')

indexed=json.loads((p/'publication-files.json').read_text()) if (p/'publication-files.json').exists() else {}
for name,identity in indexed.items():
 enforce((p/name).read_bytes(),identity['sha256'],identity['size_bytes'])

for row in r['bundles']:
 raw=(p/row['bundle']).read_bytes();enforce(raw,row['bundle_sha256'],row['bundle_size_bytes'])
 manifest=json.loads((p/'manifests'/(pathlib.Path(row['bundle']).name.removesuffix('.evidence.tar.gz')+'.json')).read_text())
 with tarfile.open(fileobj=io.BytesIO(raw),mode='r:gz') as t:
  assert t.getnames()==sorted(x['member'] for x in manifest['members'])
  for m in manifest['members']:
   member=t.getmember(m['member']);assert member.isfile();data=t.extractfile(member).read();enforce(data,m['sha256'],m['size_bytes']);assert pathlib.Path(m['source']).read_bytes()==data;count+=1
for row in json.loads((p/'views-manifest.json').read_text()):
 raw=pathlib.Path(row['source']).read_bytes();view=(p/row['view']).read_bytes();assert h(raw)==row['source_sha256'] and h(view)==row['view_sha256'];assert view==raw.replace(b'\r\n',b'\n')
 if 'raw_gzip' in row:
  gz=(p/row['raw_gzip']).read_bytes();assert h(gz)==row['raw_gzip_sha256'] and gzip.decompress(gz)==raw
for row in json.loads((p/'monitor-manifest.json').read_text()):
 raw=(p/row['copy']).read_bytes();enforce(raw,row['sha256'],row['size_bytes']);assert pathlib.Path(row['source']).read_bytes()==raw
f=json.loads((p/'artifact-FINAL.json').read_text());archive=pathlib.Path(f['path']).read_bytes();enforce(archive,f['sha256'],f['size_bytes'])
negative=[]
for data,sha,size in [(b'wrong',f['sha256'],f['size_bytes']),(archive,'0'*64,len(archive))]:
 try:enforce(data,sha,size)
 except ValueError:negative.append('rejected')
assert len(negative)==2
try:(p/'deliberately-missing-negative-control').read_bytes()
except FileNotFoundError:negative.append('missing refused')
assert len(negative)==3
result={'status':'PASS_REMOTE_PLATFORM_PACKAGING_READBACK','original_members_verified':count,'timing_parser_controls':controls['status'],'native_stages_independently_checked':controls['actual_stages_checked'],'negative_controls':negative,'release_archive_sha256':h(archive),'release_archive_bytes':len(archive),'scope':'Remote3 only; no warning-consequence, local timing or final panel adjudication.'}
(p/'readback-result.json').write_text(json.dumps(result,indent=2)+'\n')
r['status']=result['status'];r['mechanical_readback']='readback-result.json';(p/'publication-result.json').write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
files={str(x.relative_to(p)):{'sha256':h(x.read_bytes()),'size_bytes':x.stat().st_size} for x in sorted(p.rglob('*')) if x.is_file() and x.name!='publication-files.json'}
(p/'publication-files.json').write_text(json.dumps(files,indent=2,sort_keys=True)+'\n');print(json.dumps(result))
