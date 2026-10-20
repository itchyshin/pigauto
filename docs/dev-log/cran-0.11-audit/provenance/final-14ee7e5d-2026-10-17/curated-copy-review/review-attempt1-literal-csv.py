from pathlib import Path,PurePosixPath
import hashlib,json,tarfile,io,gzip,csv,collections,re
T=Path('/private/tmp/pigauto-final-closure-20261017');R=Path('/private/tmp/pigauto-pr236-validator-20261017');O=Path(__file__).parent
SHA='14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb';COMMIT='8e275ca738129c8d608f06ffc407b8141c1a2636'
def sha(b):return hashlib.sha256(b).hexdigest()
def need(x,msg):
 if not x:raise ValueError(msg)
def safe(name):
 p=PurePosixPath(name);need(not p.is_absolute() and '..' not in p.parts and name!='', 'Unsafe path '+name)
def filecheck(data,expected,label):
 need(sha(data)==expected['sha256'], 'Hash mismatch '+label);need(len(data)==expected['size_bytes'],'Size mismatch '+label)
def validate_index(root,wrapped):
 x=json.loads((root/'publication-files.json').read_text());entries=x['files'] if wrapped else x
 need(isinstance(entries,dict) and len(entries)>0,'Wrong index schema')
 actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
 need(actual==set(entries)|{'publication-files.json'},'Unindexed/missing files '+str(actual^set(entries)^{'publication-files.json'}))
 need(not any(p.is_symlink() for p in root.rglob('*')),'Symlink in publication')
 for name,v in entries.items():safe(name);filecheck((root/name).read_bytes(),v,str(root/name))
 return entries
receipt=json.loads((T/'curated-copy-receipt.json').read_text());res={};totals=collections.Counter();controls=[]
for src,dest,wrapped in [('final-proof-publication','checks-and-integration',True),('final-platform-publication','external-platforms',False)]:
 source=T/src/'publication';target=R/'docs/dev-log/cran-0.11-audit/provenance/final-14ee7e5d-2026-10-17'/dest
 si=validate_index(source,wrapped);di=validate_index(target,wrapped)
 need(si==di,'Source and destination indices differ')
 for name in set(si)|{'publication-files.json'}:need((source/name).read_bytes()==(target/name).read_bytes(),'Copy differs '+name)
 need(receipt[dest]['indexed_files_verified']==len(si) and receipt[dest]['copied_files']==len(si)+1,'Copy receipt counts wrong')
 final=json.loads((target/'artifact-FINAL.json').read_text());need(final['sha256']==SHA and final['size_bytes']==5137648 and final['source_commit']==COMMIT,'Wrong artifact binding')
 need((target/'artifact-FINAL.json').read_bytes()==(T/'FINAL.json').read_bytes(),'FINAL differs from current frozen receipt')
 bundles={};bundle_count=0;members_count=0
 for mf in sorted((target/'manifests').glob('*.json')):
  m=json.loads(mf.read_text());need(m['release_archive_sha256']==SHA and m['kind']=='evidence_bundle_not_release_archive','Wrong bundle identity')
  safe(m['bundle']);raw=(target/m['bundle']).read_bytes();need(sha(raw)==m['bundle_sha256'] and len(raw)==m['bundle_size_bytes'],'Bundle index binding failed')
  with tarfile.open(fileobj=io.BytesIO(raw),mode='r:gz') as tf:
   infos=tf.getmembers();names=[i.name for i in infos]
   need(len(names)==len(set(names))==m['member_count']==len(m['members']),'Duplicate/wrong member count')
   need(names==sorted(names),'Unsorted members')
   expected={v['member']:v for v in m['members']};need(set(names)==set(expected),'Manifest/member set mismatch')
   data={}
   for i in infos:
    safe(i.name);need(i.isfile() and not i.issym() and not i.islnk(),'Nonregular tar member')
    need(i.uid==i.gid==i.mtime==0,'Nondeterministic ownership/time metadata')
    need(not any(z in PurePosixPath(i.name).parts for z in ['exact-source','absence-library','library','libraries','__pycache__','.git']),'Excluded member '+i.name)
    need(not i.name.startswith('frozen-check/pigauto.Rcheck/pigauto/'),'Installed tree leaked')
    need(not i.name.endswith('pigauto_0.11.0.tar.gz'),'Release archive duplicated')
    b=tf.extractfile(i).read();data[i.name]=b;filecheck(b,expected[i.name],i.name)
    srcpath=Path(expected[i.name]['source']);need(srcpath.is_file() and not srcpath.is_symlink(),'Source missing/link '+str(srcpath));need(srcpath.read_bytes()==b,'Member differs from retained source '+i.name)
  # Independent decode/reencode/decode round trip of the full tar byte stream.
  decoded=gzip.decompress(raw);encoded=gzip.compress(decoded,mtime=0);need(gzip.decompress(encoded)==decoded,'Gzip roundtrip changed tar bytes')
  bundles[mf.stem]=data;bundle_count+=1;members_count+=len(data)
  if src=='final-platform-publication':
   for phase in ['download','pre','post']:
    ident=json.loads(data['platform-evidence/identity-'+phase+'.json']);need(ident['sha256']==SHA and ident['size_bytes']==5137648 and ident['archive_entries']==250 and ident['phase']==phase and ident['run_id']=='39082452558','Wrong direct identity')
   need(b'Status: OK' in data['frozen-check/pigauto.Rcheck/00check.log'],'Native platform status absent')
 vm=json.loads((target/'views-manifest.json').read_text());norm=0
 for v in vm:
  if 'source' in v:b=Path(v['source']).read_bytes()
  else:b=bundles[v['bundle']][v['member']]
  need(sha(b)==v['source_sha256'],'View source binding wrong')
  view=(target/v['view']).read_bytes();need(sha(view)==v['view_sha256'],'View hash wrong')
  if v['normalization']=='none':need(view==b,'Unchanged view differs')
  else:
   need(v['normalization']=='CRLF_TO_LF_ONLY' and view==b.replace(b'\r\n',b'\n'),'Other normalization applied')
   rg=(target/v['raw_gzip']).read_bytes();need(sha(rg)==v['raw_gzip_sha256'] and gzip.decompress(rg)==b,'Raw gzip lost original bytes');norm+=1
   if v['view'].endswith('.csv'):need(list(csv.reader(io.StringIO(b.decode())))==list(csv.reader(io.StringIO(view.decode()))),'CSV semantics changed')
 if (target/'monitor-manifest.json').exists():
  for m in json.loads((target/'monitor-manifest.json').read_text()):
   b=(target/m['copy']).read_bytes();filecheck(b,m,m['copy']);need(b==Path(m['source']).read_bytes(),'Monitor original source differs')
 res[dest]={'index_schema':'files_wrapper' if wrapped else 'direct_map','indexed_files':len(si),'copied_files':len(si)+1,'bundles':bundle_count,'raw_members':members_count,'readable_views':len(vm),'normalized_views':norm,'source_and_target_exact':True,'raw_members_equal_retained_sources':True}
 totals.update({'indexed':len(si),'bundles':bundle_count,'members':members_count,'views':len(vm)})
# Negative controls exercise independent validation; only in-memory altered data.
base={'sha256':sha(b'good'),'size_bytes':4}
for label,fn in [('changed-hash',lambda:filecheck(b'evil',base,'planted')),('changed-size',lambda:filecheck(b'good',base|{'size_bytes':5},'planted')),('traversal',lambda:safe('../escape')),('absolute-path',lambda:safe('/escape'))]:
 try:fn()
 except ValueError:controls.append(label+' refused')
 else:raise AssertionError(label+' accepted')
final=json.loads((T/'FINAL.json').read_text());need(sha(Path(final['path']).read_bytes())==SHA,'Final release archive changed')
res={'verdict':'PASS_BOUNDED_CURATED_COPY_MECHANICS','archive_sha256':SHA,'source_commit':COMMIT,'trees':res,'totals':dict(totals),'negative_controls':controls,'limits':['No R/check execution or network used.','No final release-panel, warning-consequence or timing-adjudication verdict.','Source-file equality proves retained-copy binding; remote executed merge/head tree binding is separately noted.']}
(O/'receipt.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
