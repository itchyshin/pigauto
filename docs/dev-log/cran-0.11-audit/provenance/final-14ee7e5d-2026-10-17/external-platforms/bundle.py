#!/usr/bin/env python3
"""Package already-executed evidence only; no checks, network, or release changes."""
import gzip, hashlib, io, json, pathlib, platform, re, sys, tarfile, zlib

BASE = pathlib.Path('/private/tmp/pigauto-final-closure-20261017')
HERE = pathlib.Path(__file__).resolve().parent
FINAL = json.loads((BASE / 'FINAL.json').read_bytes())
SHA = FINAL['sha256']
assert SHA == '14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb'
OUT = HERE / 'publication'

def digest(b):
    return hashlib.sha256(b).hexdigest()

def checked_member(root, name, expected=None):
    rel = pathlib.PurePosixPath(name)
    if rel.is_absolute() or '..' in rel.parts or not rel.parts:
        raise ValueError('Unsafe evidence member')
    p = root.joinpath(*rel.parts)
    if not p.is_file() or p.is_symlink() or any(x.is_symlink() for x in p.parents if x != root.parent):
        raise ValueError('Missing/nonregular/symlink evidence member: ' + str(p))
    if root.resolve() not in p.resolve().parents:
        raise ValueError('Evidence member escaped root')
    b = p.read_bytes()
    if expected is not None and digest(b) != expected:
        raise ValueError('Evidence hash mismatch: ' + str(p))
    # Stop instead of publishing credential-looking material. No environment capture.
    if re.search(rb'(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9]{32,}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----)', b):
        raise ValueError('Credential-shaped material: ' + str(p))
    return p, b

def archive_identity():
    b = pathlib.Path(FINAL['path']).read_bytes()
    if digest(b) != SHA or len(b) != FINAL['size_bytes']:
        raise ValueError('Final release archive changed')
    return {'sha256': digest(b), 'size_bytes': len(b)}

def deterministic_gzip(b):
    out = io.BytesIO()
    with gzip.GzipFile(fileobj=out, mode='wb', filename='', mtime=0, compresslevel=9) as g:
        g.write(b)
    return out.getvalue()

def tar_gzip(members):
    out = io.BytesIO()
    with tarfile.open(fileobj=out, mode='w', format=tarfile.PAX_FORMAT) as t:
        for name, (_, b) in sorted(members.items()):
            info = tarfile.TarInfo(name)
            info.size, info.mode, info.mtime = len(b), 0o644, 0
            info.uid = info.gid = 0
            info.uname = info.gname = ''
            t.addfile(info, io.BytesIO(b))
    return deterministic_gzip(out.getvalue())

def read_manifest_bundle(root, manifest_name):
    p, original = checked_member(root, manifest_name)
    declared = json.loads(original)
    members = {name: checked_member(root, name, h) for name, h in declared.items()}
    members[manifest_name] = (p, original)
    return members

def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')

import csv
SOURCE = BASE / 'final-direct-39082452558'

STAGE_PHRASES = [('examples','examples'),('tests','tests'),('vignettes','re-building of vignette outputs'),('PDF','PDF version of manual'),('HTML','HTML version of manual')]

def parse_stage_timing(line):
    """Accept exact native OK line formats, refusing unfamiliar timing syntax."""
    match=re.fullmatch(r'\* checking .+ \.\.\. (?:\[([0-9]+(?:\.[0-9]+)?)s(?:/([0-9]+(?:\.[0-9]+)?)s)?\] )?OK',line)
    if not match:raise ValueError('Unrecognized native stage line: '+line)
    result={'raw_check_line':line}
    if match.group(1) is None:result.update(timing_format='not_reported')
    elif match.group(2) is None:result.update(timing_format='single_reported_seconds',reported_seconds=float(match.group(1)))
    else:result.update(timing_format='cpu_elapsed_seconds',cpu_seconds=float(match.group(1)),elapsed_seconds=float(match.group(2)))
    return result

def parse_native_stages(log):
    stages={}
    for key,phrase in STAGE_PHRASES:
        lines=[x for x in log.splitlines() if x.startswith('* checking '+phrase+' ...')]
        if len(lines)!=1:raise ValueError('Missing/duplicate native stage: '+phrase)
        stages[key]=parse_stage_timing(lines[0])
    return stages

def main():
    before = archive_identity()
    if OUT.exists(): raise ValueError('Refusing to overwrite publication')
    OUT.mkdir()
    for d in ['raw','manifests','views','monitor-originals']: (OUT/d).mkdir()
    common = [SOURCE/'closure-receipt.txt', SOURCE/'verify-results.out', BASE/'root-final-direct-reverify.log']
    metadata = json.loads((HERE/'metadata-sources.json').read_text())
    common.extend(pathlib.Path(p) for p in metadata)
    primary={pathlib.Path(p).name:json.loads(pathlib.Path(p).read_bytes()) for p in metadata}
    run=primary['run.json'];jobs=primary['jobs.json']['jobs'];artifacts=primary['artifacts.json']['artifacts']
    assert run['id']==39082452558 and run['status']=='completed' and run['conclusion']=='success'
    assert run['head_sha']=='6cd5f65b06d0c3ee32ce4b2d5a67a22c06028116'
    assert len(jobs)==len(artifacts)==3 and all(j['conclusion']=='success' and j['run_id']==run['id'] for j in jobs)
    zip_receipts=[]
    closure=(SOURCE/'closure-receipt.txt').read_text()
    for z in sorted(SOURCE.glob('*.zip')):
        raw=z.read_bytes();h=digest(raw)
        assert h in closure, 'ZIP digest absent from original closure receipt'
        artifact=next(a for a in artifacts if a['name']==z.stem)
        assert artifact['digest']=='sha256:'+h and artifact['size_in_bytes']==len(raw) and not artifact['expired']
        zip_receipts.append({'original_zip':str(z),'sha256':h,'size_bytes':len(raw),'original_closure_receipt':'monitor-originals/00-closure-receipt.txt','zip_included':False})
    assert len(zip_receipts)==3
    write_json(OUT/'zip-digest-readback.json',zip_receipts)
    common_receipts=[]
    for index, p in enumerate(common):
        _,raw=checked_member(p.parent,p.name)
        rel='monitor-originals/'+str(index).zfill(2)+'-'+p.name
        (OUT/rel).write_bytes(raw)
        common_receipts.append({'source':str(p),'copy':rel,'sha256':digest(raw),'size_bytes':len(raw)})
    receipts=[]; views=[]; summaries=[]
    verification=json.loads((BASE/'root-final-direct-reverify.log').read_text().split('\nDIRECT_ARCHIVE')[0])
    for root in sorted(SOURCE.glob('frozen-14ee7e5d-*')):
        if not root.is_dir(): continue
        label=root.name.removeprefix('frozen-14ee7e5d-')
        members={}
        pe=root/'platform-evidence'; rc=root/'frozen-check/pigauto.Rcheck'
        for p in sorted(pe.rglob('*')):
            if p.is_file():
                rel=str(p.relative_to(root));members[rel]=checked_member(root,rel)
        for p in sorted(rc.iterdir()):
            if p.is_file():
                rel=str(p.relative_to(root)); members[rel]=checked_member(root,rel)
        for name in ['testthat.Rout']:
            rel='frozen-check/pigauto.Rcheck/tests/'+name
            members[rel]=checked_member(root,rel)
        for phase in ['download','pre','post']:
            identity=json.loads(members['platform-evidence/identity-'+phase+'.json'][1])
            assert identity['sha256']==SHA and identity['size_bytes']==FINAL['size_bytes']
            assert identity['inventory_sha256']==FINAL['inventory_sha256'] and identity['run_id']=='39082452558'
        log=members['frozen-check/pigauto.Rcheck/00check.log'][1].decode()
        assert re.search(r'^Status: OK\s*$',log,re.M)
        checkhash=digest(members['frozen-check/pigauto.Rcheck/00check.log'][1])
        assert any(row['check_log_sha256']==checkhash for row in verification)
        deps=list(csv.DictReader(io.StringIO(members['platform-evidence/declared-dependencies.csv'][1].decode())))
        assert len(deps)==29 and all(row['available']=='TRUE' for row in deps)
        warning=json.loads(members['platform-evidence/warning-replay.json'][1])
        expected={'fail':0,'warnings':161,'skips':86,'pass':3195}
        assert warning['native_counts']==warning['replay_counts']==expected
        native_text=members['frozen-check/pigauto.Rcheck/tests/testthat.Rout'][1].decode()
        assert '[ FAIL 0 | WARN 161 | SKIP 86 | PASS 3195 ]' in native_text
        assert len(list(csv.DictReader(io.StringIO(members['platform-evidence/warning-details.csv'][1].decode()))))==161
        assert warning['not_cran']=='UNSET' and warning['force_suggests']=='true'
        exittext=members['platform-evidence/check-exit.txt'][1].decode();assert exittext.splitlines()[0]=='exit=0'
        stages=parse_native_stages(log)
        ex_text=members['frozen-check/pigauto.Rcheck/pigauto-Ex.Rout'][1].decode()
        ex=re.findall(r'^Time elapsed:\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)',ex_text,re.M)
        nt=re.search(r'user\s+system\s+elapsed\s*\n\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)',native_text)
        process_timings={'examples_Rout_elapsed_seconds':float(ex[-1][2]) if ex else None,'tests_Rout_elapsed_seconds':float(nt.group(3)) if nt else None,'scope':'Rout-reported process times, distinct from check-stage elapsed brackets and separate warning replay.'}
        stages['total_seconds']=float(re.search(r'elapsed_seconds=([0-9.]+)',exittext).group(1))
        summary={'platform':label,'native_status':'OK','native_counts':expected,'replay_counts':expected,'warning_consequences':'NOT_ADJUDICATED','declared_dependency_count':len(deps),'all_declared_dependencies_available':True,'stage_timings':stages,'Rout_process_timings':process_timings,'separate_warning_replay_seconds':warning['diagnostic_elapsed_seconds'],'check_log_sha256':checkhash}
        summaries.append(summary)
        raw=tar_gzip(members);assert raw==tar_gzip(members)
        target='raw/'+label+'.evidence.tar.gz';(OUT/target).write_bytes(raw)
        rows=[]
        for name,(p,data) in sorted(members.items()):
            rows.append({'member':name,'source':str(p),'sha256':digest(data),'size_bytes':len(data)})
            # Key human-readable logs, identities, dependency/timing and warning tables.
            selected=p.suffix in ['.txt','.json','.csv','.timings'] or p.name in ['00check.log','00install.out','testthat.Rout','pigauto-Ex.Rout','build_vignettes.log','Rdlatex.log','pigauto-manual.log']
            if selected:
                view='views/'+label+'/'+name; dest=OUT/view;dest.parent.mkdir(parents=True,exist_ok=True)
                normalized=data.replace(b'\r\n',b'\n');dest.write_bytes(normalized)
                row={'source':str(p),'source_sha256':digest(data),'view':view,'view_sha256':digest(normalized),'normalization':'CRLF_TO_LF_ONLY' if data!=normalized else 'none'}
                if data!=normalized:
                    gz=deterministic_gzip(data);(OUT/(view+'.raw.gz')).write_bytes(gz);row.update(raw_gzip=view+'.raw.gz',raw_gzip_sha256=digest(gz))
                views.append(row)
        receipt={'kind':'evidence_bundle_not_release_archive','release_archive_sha256':SHA,'bundle':target,'bundle_sha256':digest(raw),'bundle_size_bytes':len(raw),'member_count':len(rows),'members':rows}
        write_json(OUT/'manifests'/(label+'.json'),receipt);receipts.append({k:v for k,v in receipt.items() if k!='members'})
        for p,data in members.values():assert p.read_bytes()==data
    assert len(receipts)==3
    after=archive_identity();assert before==after
    write_json(OUT/'platform-summary.json',summaries);write_json(OUT/'views-manifest.json',views);write_json(OUT/'monitor-manifest.json',common_receipts)
    write_json(OUT/'publication-result.json',{'status':'PACKAGED_PENDING_MECHANICAL_READBACK','bundles':receipts,'archive_before':before,'archive_after':after,'python':sys.version,'zlib':zlib.ZLIB_VERSION,'original_members':sum(x['member_count'] for x in receipts),'normalized_views':sum(x['normalization']!='none' for x in views),'scope':'Remote three only. Warning consequences, local macOS timing/note and final panel remain unadjudicated.'})
    for name in ['bundle.py','readback.py','metadata-sources.json','timing-parser-controls.py','timing-parser-independent-reference.json']:(OUT/name).write_bytes((HERE/name).read_bytes())
    (OUT/'artifact-FINAL.json').write_bytes((BASE/'FINAL.json').read_bytes())
    (OUT/'README.md').write_text('# Final remote platform evidence\n\nThese are evidence bundles, not the release archive. Each raw gzip preserves original files for direct run39082452558, with source path/SHA/size bindings. All platform-evidence regular files, check-root output files, native testthat.Rout and retained actual manual PDFs/logs are included. Installed package, copied source/tests/fixtures and duplicate ZIPs are excluded. ZIP digests remain in monitor receipts.\n\nReadable views apply CRLF-to-LF only, with original-byte .raw.gz for changed views. platform-summary.json records native stage timings and 29 available declared dependencies per platform. Linux native timings use rounded CPU/elapsed pairs (examples 14/14, tests 206/206, vignettes 11/11, PDF 2/2, HTML 2/2 seconds). Windows single values remain labelled reported_seconds. Unrounded Rout process timings remain separate. The first derived parser accepted only single timings and incorrectly reported Linux timings absent; that derived version is retained outside this curated output. timing-parser-controls.json and the independent reference check paired, single, missing and malformed formats. Native StatusOK and3195pass/161warnings/86skips/0fail are observations; warning consequences have not been adjudicated. Separate warning replay time is excluded from native timing. This does not close the local Mac timingNOTE or the final panel.\n\nRun python3 bundle.py from a fresh helper directory retaining metadata-sources.json and readback.py; publication output must not already exist. Run python3 readback.py publication to mechanically recheck bundles, sources, views and final archive. No network/check execution occurs. Original NOTASSESSED receipts remain unchanged. Writing-style assessment of this packaging README: NOTASSESSED.\n')
    print(json.dumps({'bundles':receipts,'members':sum(x['member_count'] for x in receipts)}))

if __name__=='__main__':main()
