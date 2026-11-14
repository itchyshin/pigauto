#!/usr/bin/env python3
"""Read retained proof only. No subprocesses, fits, network, or release adjudication."""
import argparse, csv, datetime, hashlib, io, json, pathlib, re, shutil, tarfile

BASE = pathlib.Path('/private/tmp/pigauto-final-closure-20261017')
HERE = pathlib.Path(__file__).resolve().parent
SHA = '14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb'
INPUTS = {}
STAGES = {'examples':'examples','tests':'tests','vignettes':'re-building of vignette outputs','PDF':'PDF version of manual','HTML':'HTML version of manual'}

def sha(raw): return hashlib.sha256(raw).hexdigest()
def write_json(p, value): p.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
def read(p, mutable=False):
    p=pathlib.Path(p);raw=p.read_bytes()
    if p.read_bytes()!=raw: raise ValueError('Input changed during read: '+str(p))
    INPUTS[str(p)]={'sha256':sha(raw),'size_bytes':len(raw),'mutable_observation':mutable}
    if mutable:
        dest=HERE/'local-snapshots'/p.name
        dest.parent.mkdir(exist_ok=True);dest.write_bytes(raw)
        INPUTS[str(p)]['snapshot']=str(dest.relative_to(HERE))
    return raw
def load(p, mutable=False): return json.loads(read(p,mutable))
def table(p, mutable=False): return list(csv.DictReader(io.StringIO(read(p,mutable).decode()),delimiter=','))
def counts(text):
    hits=re.findall(r'\[ FAIL (\d+) \| WARN (\d+) \| SKIP (\d+) \| PASS (\d+) \]',text)
    return dict(zip(['fail','warnings','skips','pass'],map(int,hits[-1]))) if hits else None
def timing_line(line):
    # Native logs can report rounded CPU/elapsed in seconds or minutes.
    hit=re.fullmatch(r'\* checking .+ \.\.\. (?:\[(\d+(?:\.\d+)?)(s|m)(?:/(\d+(?:\.\d+)?)(s|m))?\] )?(OK|NOTE)',line)
    if not hit: raise ValueError('Unknown native timing format: '+line)
    out={'raw_line':line,'check_result':hit[5]}
    if hit[1] is None: out['format']='not_reported'
    elif hit[3] is None: out.update(format='single_reported',reported_seconds=float(hit[1])*(60 if hit[2]=='m' else 1))
    else: out.update(format='rounded_CPU_elapsed',cpu_seconds=float(hit[1])*(60 if hit[2]=='m' else 1),elapsed_seconds=float(hit[3])*(60 if hit[4]=='m' else 1))
    return out
def stages(text, complete=False):
    out={}
    for key,phrase in STAGES.items():
        lines=[x for x in text.splitlines() if x.startswith('* checking '+phrase+' ...')]
        if len(lines)>1: raise ValueError('Duplicate check stage')
        if not lines or lines[0].endswith(' ...'):
            if complete: raise ValueError('Missing completed stage')
            out[key]={'status':'PENDING'}
        else: out[key]=timing_line(lines[0])
    return out
def footprint(text):
    hit=re.search(r'(\* checking installed package size \.\.\. INFO\r?\n.*?)(?=\* checking package directory)',text,re.S)
    if not hit: raise ValueError('Missing native installed-size INFO')
    block=hit[1];size=re.search(r'installed size is\s+([\d.]+Mb)',block)
    return {'status':'INFO','reported_installed_size':size[1],'raw_block':block,'interpretation':'Native informational size output, not a warning or a release-budget adjudication.'}
def examples(path, mutable=False):
    text=read(path,mutable).decode();rows=list(csv.DictReader(io.StringIO(text),delimiter='\t'))
    parsed=[{'name':r['name'],'user_seconds':float(r['user']),'system_seconds':float(r['system']),'elapsed_seconds':float(r['elapsed'])} for r in rows]
    return {'measured_examples':len(parsed),'over_five_cpu_or_elapsed_seconds':[r for r in parsed if r['user_seconds']+r['system_seconds']>5 or r['elapsed_seconds']>5],'compare_methods':next(r for r in parsed if r['name']=='compare_methods')}
def verify_declared(root, filename):
    manifest=load(root/filename)
    for name,digest in manifest.items():
        rel=pathlib.PurePosixPath(name)
        if rel.is_absolute() or '..' in rel.parts: raise ValueError('Unsafe receipt member')
        assert sha(read(root/name))==digest, name
    return len(manifest)

def local_observation(final, refresh):
    root=pathlib.Path(final['path']).parent/'local-macos';rc=root/'frozen-check/pigauto.Rcheck';pe=root/'platform-evidence'
    text=read(rc/'00check.log',True).decode();status=[x for x in text.splitlines() if x.startswith('Status:')]
    completed='* DONE' in text.splitlines() and len(status)==1
    known='compare_methods 2.532  0.992  18.406'
    assert known in text, 'Known example NOTE changed; independent review required'
    expected_note_lines=[x for x in text.splitlines() if x.startswith('* checking') and x.endswith(' NOTE')]
    assert expected_note_lines==['* checking examples ... [13s/63s] NOTE']
    assert not any(x.startswith('* checking') and ('ERROR' in x or 'WARNING' in x) for x in text.splitlines())
    out={'native_state':'COMPLETED_UNADJUDICATED' if completed else 'PENDING_NATIVE_COMPLETION','native_status_line':status[0] if completed else None,'snapshot_check_log_sha256':sha(text.encode()),'stage_timings':stages(text,completed),'known_note':{'example':'compare_methods','user_seconds':2.532,'system_seconds':0.992,'cpu_user_plus_system_seconds':3.524,'elapsed_seconds':18.406,'disposition':'UNADJUDICATED; preserved, no faster rerun'},'installed_size_INFO':footprint(text),'warning_replay_state':'PENDING_NOT_YET_PROVEN','warning_consequences':'PENDING','total_native_check_seconds':None,'check_driver_exit':None,'native_process_exit':None,'thread_runtime_probe':load(BASE/'local-torch-runtime-probe.json'),'observation_only':not refresh}
    out['example_timings']=examples(rc/'pigauto-Ex.timings',True)
    if (rc/'tests/testthat.Rout').exists(): out['observed_native_test_counts']=counts(read(rc/'tests/testthat.Rout',True).decode())
    deps=table(pe/'declared-dependencies.csv',True);assert len(deps)==29 and all(r['available']=='TRUE' for r in deps);out['declared_dependencies_available']=29
    if (pe/'check-exit.txt').exists():
        exittext=read(pe/'check-exit.txt',True).decode();out['native_process_exit']=int(re.search(r'^exit=(\d+)',exittext,re.M)[1]);out['total_native_check_seconds']=float(re.search(r'elapsed_seconds=([\d.]+)',exittext)[1])
    if (root/'local-receipt.json').exists():
        receipt=load(root/'local-receipt.json',True);assert receipt['artifact_sha256']==SHA;out['check_driver_receipt']=receipt
        out['check_driver_exit']=receipt['records'][0]['exit_code']
    if refresh:
        assert completed and status==['Status: 1 NOTE'], 'Refresh requires completed exact known native NOTE'
        assert out['native_process_exit']==0 and out['check_driver_exit'] not in [None,0], 'Preserve native success and strict driver rejection separately'
        for phase in ['pre','post']:
            identity=load(pe/('identity-'+phase+'.json'),True);assert identity['sha256']==SHA and identity['size_bytes']==final['size_bytes']
        if (root/'local-note-diagnostic.json').exists():
            diag=load(root/'local-note-diagnostic.json',True);assert diag['artifact_sha256']==diag['archive_sha256_after']==SHA
            assert diag['native_check_log_sha256']==out['snapshot_check_log_sha256']
            out['warning_diagnostic_receipt']=diag
            if diag['exit_code']==0:
                replay=load(pe/'warning-replay.json',True);assert replay['native_counts']==replay['replay_counts']=={'fail':0,'warnings':161,'skips':86,'pass':3195}
                assert replay['not_cran']=='UNSET' and replay['force_suggests']=='true'
                assert len(table(pe/'warning-details.csv',True))==161
                out.update(warning_replay_state='CAPTURE_COMPLETE_CONSEQUENCES_PENDING',warning_replay=replay)
            else: out['warning_replay_state']='DIAGNOSTIC_FAILED_NO_PROMOTION'
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--refresh-local',action='store_true',help='After completion, read exact known NOTE/native exit/driver/replay receipts; never run R or promote release status');args=parser.parse_args()
    now=datetime.datetime.now(datetime.timezone.utc).isoformat();final=load(BASE/'FINAL.json');state=load(BASE/'STATE.json',True)
    assert final['sha256']==state['final_archive_sha256']==SHA
    archive=read(final['path']);assert sha(archive)==SHA and len(archive)==final['size_bytes']==5137648
    # Preserve each prior generated synopsis and its local input snapshots before update.
    if (HERE/'synopsis.json').exists():
        history=HERE/'history'/now.replace(':','-');history.mkdir(parents=True)
        for name in ['synopsis.json','README.md','inputs.json','assessment.json','assessment-check.txt','output-manifest.json']:
            if (HERE/name).exists(): shutil.copyfile(HERE/name,history/name)
        if (HERE/'local-snapshots').exists(): shutil.copytree(HERE/'local-snapshots',history/'local-snapshots')
    with tarfile.open(fileobj=io.BytesIO(archive),mode='r:gz') as t:
        entries=t.getmembers();regular=[x for x in entries if x.isfile()]
        assert len(entries)==final['archive_entries']==250 and len(regular)==236
        assert sum(x.size for x in regular)==final['regular_uncompressed_bytes']==7150500
        assert not any(x.name.startswith('pigauto/src/') for x in entries)
        spelling=load(BASE/'final-spelling-review/receipt.json');assert spelling['archive_sha256']==SHA
        for name,digest in spelling['source_files_bound_to_archive'].items(): assert sha(t.extractfile('pigauto/'+name).read())==digest
    platform=BASE/'final-platform-publication/publication';remote=load(platform/'platform-summary.json');assert len(remote)==3
    reference=load(platform/'timing-parser-independent-reference.json');remote_out=[]
    for row in remote:
        root=BASE/'final-direct-39082452558'/('frozen-14ee7e5d-'+row['platform']);rc=root/'frozen-check/pigauto.Rcheck';pe=root/'platform-evidence'
        raw=read(rc/'00check.log');assert sha(raw)==row['check_log_sha256']==reference[root.name]['check_log_sha256']
        text=raw.decode();assert text.splitlines()[-1]=='Status: OK';assert row['native_status']=='OK'
        for phase in ['pre','post']:
            identity=load(pe/('identity-'+phase+'.json'));assert identity['sha256']==SHA and identity['size_bytes']==5137648
        deps=table(pe/'declared-dependencies.csv');assert len(deps)==29 and all(r['available']=='TRUE' for r in deps)
        assert counts(read(rc/'tests/testthat.Rout').decode())==row['native_counts']=={'fail':0,'warnings':161,'skips':86,'pass':3195}
        replay=load(pe/'warning-replay.json');assert replay['replay_counts']==row['native_counts']
        assert len(table(pe/'warning-details.csv'))==161
        timing=stages(text,True);total=float(re.search(r'elapsed_seconds=([\d.]+)',read(pe/'check-exit.txt').decode())[1]);assert total==row['stage_timings']['total_seconds']
        remote_out.append({'platform':row['platform'],'native_status':'OK','native_total_seconds':total,'native_total_over_600_seconds':total>600,'stage_timings':timing,'Rout_process_timings':row['Rout_process_timings'],'example_timings':examples(rc/'pigauto-Ex.timings'),'installed_size_INFO':footprint(text),'all_declared_dependencies_available':29,'native_counts':row['native_counts'],'separate_warning_replay_seconds':row['separate_warning_replay_seconds'],'warning_consequences':'PENDING_EXACT_FINAL_AUDIT','check_log_sha256':sha(raw)})
    service=BASE/'grace-proof-14ee7e5d';assert verify_declared(service,'proof-manifest.json')==31
    absence=load(service/'absence-summary.json');assert absence=={'pass':60,'fail':0,'warn':0,'skip':8,'files':['test-gbif-centroids.R','test-worldclim-covariates.R'],'not_cran':'UNSET','rgbif':False,'terra':False}
    adjudication=load(service/'adjudication.json');urls=table(service/'url-network-results-flat.csv');assert len(urls)==1 and urls[0]['Status']=='403'
    vignette_root=BASE/'vignette-timing-14ee7e5d';assert verify_declared(vignette_root,'receipt-manifest.json')==45
    vignettes=load(vignette_root/'timing-summary.json');vignette_identity=load(vignette_root/'identity-after.json');disabled=load(vignette_root/'eval-false-proof.json');assert vignettes['artifact']['sha256']==SHA and len(vignettes['rows'])==6 and all(r['exit']==0 for r in vignettes['rows'])
    assert disabled['explicit_eval_false_chunks_preserved']==32 and disabled['all_checked_disabled']
    doi=load(vignette_root/'doi-metadata.json');crossref=load(vignette_root/'crossref-raw.json');assert doi['http_status']==200 and doi['identity_matches'] and sha(read(vignette_root/'crossref-raw.json'))==doi['sha256']
    assert doi['doi']==crossref['message']['DOI']=='10.1126/science.1157704'
    assert len(table(BASE/'final-spelling-review/token-adjudication.csv'))==spelling['flag_rows']==404
    assert len(table(BASE/'final-spelling-review/source-contexts.csv'))==spelling['resolved_context_rows']==1886
    assert len(spelling['source_files_bound_to_archive'])==43 and sum(spelling['categories'].values())==404 and spelling['ordinary_word_typos_identified']==spelling['uncertain_flags']==[]
    mi=load(pathlib.Path(final['path']).parent/'mi-matrix/manifest.json');mi_verify=load(BASE/'final-mi-verification.json');assert mi['archive_sha256']==mi_verify['artifact_sha256']==SHA and mi_verify['artifact_source_commit']==final['source_commit']
    steps=[s for c in mi['cells'] for s in c['steps']];assert len(steps)==12
    for s in steps:
        assert s['exit_code']==0 and s['archive_sha256_before']==s['archive_sha256_after']==SHA
        assert sha(read(pathlib.Path(final['path']).parent/'mi-matrix'/s['log']))==s['log_sha256']
    profile=load(BASE/'cran-publication-and-revdeps.json');size=load(BASE/'final-size-static.json');assert size['artifact_sha256']==SHA
    local=local_observation(final,args.refresh_local)
    contradictions=[{'finding':'Earlier Linux timing absence claim was false','resolution':'Corrected paired CPU/elapsed parser, ten known-answer controls and15independent stage comparisons; raw log and bundles unchanged.','source':str(platform/'timing-correction-receipt.json')},{'finding':'Original service adjudication still describes broad spelling review as open','resolution':'Keep original historical receipt; final-spelling-review now covers all404flagged tokens, with scoped limitations.','source':str(BASE/'final-spelling-review/receipt.json')},{'finding':'STATE warning_capture refers to predecessor c2 and new_checker/active_build text is historical','resolution':'Use final14ee receipts for artifact-dependent conclusions; do not transfer predecessor warning consequences or old running-state wording.','source':str(BASE/'STATE.json')},{'finding':'MI manifest source_commit is a checker identity, not the artifact generating source','resolution':'Artifact source8e275ca and checker identities are recorded separately by final-mi-verification.json.','source':str(BASE/'final-mi-verification.json')},{'finding':'Large Windows total is distinct from vignette rebuild duration','resolution':'Retain645.8126/622.4320second native totals and22second vignette stages; no claim that vignette rebuilding exceeds600seconds. Timing decision remains open.','source':str(platform/'platform-summary.json')}]
    result={'generated_utc':now,'scope':'Evidence synopsis only; no release rung promoted, final panel or readiness adjudication.','overall':'NOT_READY_PENDING_LOCAL_TIMING_WARNING_REVIEW_AND_FINAL_PANEL','artifact':{k:final[k] for k in ['path','sha256','size_bytes','archive_entries','inventory_sha256','source_commit','clean_worktree']},'profile':{'release_type':'update','retained_CRAN_version':profile['version'],'retained_reverse_dependencies':profile['reverse_dependencies'],'scope':'Retained CRAN snapshot, not a new network query.'},'remote_native_checks':remote_out,'local_macos':local,'warnings':{'remote_each':{'pass':3195,'fail':0,'warnings':161,'skips':86},'consequences':'PENDING_FINAL_WARNING_AUDIT','predecessor_classification_is_not_final_evidence':True},'timing_decision':{'state':'OPEN_MAINTAINER_AND_FINAL_PANEL_DECISION','observed_signal_seconds':600,'threshold_type':'Observed incoming signal/conservative release review, not an immutable policy limit. Protocol strict first-submission wording is distinct; retained CRAN version0.10.0 makes this an update.','no_faster_rerun':True,'local_note_preserved':True,'native_totals_do_not_include_separate_warning_replay':True},'footprint':{'compressed_bytes':len(archive),'regular_uncompressed_bytes':sum(x.size for x in regular),'archive_entries':len(entries),'regular_files':len(regular),'package_owned_src':[],'largest_files':size['largest_files'][:5],'native_INFO_is_not_WARNING':True},'individual_vignette_diagnostic':{'separate_from_native_and_Windows_timing':True,'rows':vignettes['rows'],'render_elapsed_sum_seconds':sum(r['render_elapsed_seconds'] for r in vignettes['rows']),'whole_R_process_elapsed_sum_seconds':sum(r['process_elapsed_seconds'] for r in vignettes['rows']),'enclosing_loop_elapsed_seconds':vignette_identity['total_process_seconds'],'disabled_chunks':disabled,'scope':'Exact Rmd hashes/installed final package, CPU torch, explicit libpaths, NOT_CRAN unset/threadcaps1. No evalFALSE teaching chunks enabled.'},'external_service_absence':{'counts':absence,'expected_skips':adjudication['absence']['expected_skips'],'verified_scope':adjudication['absence']['reached'],'limits':adjudication['absence']['limits'],'live_availability':'NOT_PROVEN'},'urls':{'finding':urls[0],'access_state':'403_FORBIDDEN_RETAINED_NOT_URL_CLEAN','primary_metadata':doi,'bibliographic_identity':'SUPPORTED_BY_RETAINED_CROSSREF_HTTP200','limits':'Metadata supports identifier identity, not automated DOI access, general reference quality or live GBIF/WorldClim service availability.'},'spelling':{k:spelling[k] for k in ['flag_rows','unique_tokens','resolved_context_rows','categories','scopes','ordinary_word_typos_identified','uncertain_flags','limitations']},'spelling_archive_source_files_verified':43,'installed_MI':{'result':mi_verify['result'],'cells':[{'cell':c['cell'],'expected_optional_packages':c['expected_optional_packages'],'steps':[{'step':s['step'],'exit_code':s['exit_code'],'elapsed_seconds':s['elapsed_seconds'],'log':s['log']} for s in c['steps']]} for c in mi['cells']],'steps_passed':12,'summed_step_seconds':sum(s['elapsed_seconds'] for s in steps),'artifact_source_commit':mi_verify['artifact_source_commit'],'checker_start_commit':mi_verify['checker_start_commit'],'published_checker_commit':mi_verify['published_checker_commit'],'backend_packages':{k:{key:v[key] for key in ['version','source_commit','files_modified_since_matrix_start']} for k,v in mi_verify['backend_packages'].items()},'limits':mi_verify['limits']},'contradictions_and_repairs':contradictions+[{'finding':'Earlier310.170second process-total wording conflates enclosing-loop duration with summed six R-process durations','resolution':'Retain308.282456584seconds summed R processes versus310.169622875seconds enclosing sequential loop; render sum 74.8865 seconds.','source':str(vignette_root/'identity-after.json')}],'not_authorized_or_claimed':['new R runs or live network checks','warning consequence acceptance','local timing NOTE acceptance','final panel verdict','merge of evidence PR 236','CRAN upload or confirmation'],'governance':{'STATE_overall_snapshot':state['overall'],'source_PR237':state['source_pr_state'],'evidence_PR236':'UNMERGED per retained authority; not live reverified here','CRAN_submission_authorized':False},'refresh_command':'python3 '+str(HERE/'build-synopsis.py')+' --refresh-local'}
    assert sha(pathlib.Path(final['path']).read_bytes())==SHA
    write_json(HERE/'synopsis.json',result);write_json(HERE/'inputs.json',INPUTS)
    (HERE/'README.md').write_text(render_readme(result))
    print(json.dumps({'status':'SYNOPSIS_WRITTEN_NO_RELEASE_PROMOTION','local_state':local['native_state'],'local_replay':local['warning_replay_state'],'input_receipts':len(INPUTS),'path':str(HERE/'synopsis.json')}))

def render_readme(s):
    remote=s['remote_native_checks'];rows=[]
    for r in remote:
        ts=r['stage_timings'];fmt=lambda key: (str(ts[key]['cpu_seconds'])+'/'+str(ts[key]['elapsed_seconds'])) if ts[key]['format']=='rounded_CPU_elapsed' else str(ts[key]['reported_seconds'])
        rows.append('| '+r['platform']+' | '+f"{r['native_total_seconds']:.3f}"+' | '+' | '.join(fmt(k) for k in ['examples','tests','vignettes','PDF','HTML'])+' |')
    vignette_rows=['| '+r['file']+' | '+f"{r['render_elapsed_seconds']:.4f}"+' | '+f"{r['process_elapsed_seconds']:.3f}"+' |' for r in s['individual_vignette_diagnostic']['rows']]
    local=s['local_macos']
    return ('''# Final checks evidence synopsis

This is an evidence snapshot for pigauto 0.11.0, not a release verdict. The exact archive is 5,137,648 bytes, SHA-256 `'''+SHA+'''`, generated from clean merged source `8e275ca738129c8d608f06ffc407b8141c1a2636`. It contains 250 entries/236 regular files and 7,150,500 regular uncompressed bytes. `synopsis.json` retains exact values; `inputs.json` binds them to retained receipts.

All three remote native checks report Status:OK. Each has 29 available declared dependencies and 3,195 passes, 161 warnings, 86 skips and no failures. The 161 warnings and 86 skips remain observations pending final consequence review; their counts are not a blanket benign classification. Separate warning replays matched these counts and are excluded from native totals.

| Remote platform | Native total s | Examples s | Tests s | Vignette rebuild s | PDF s | HTML s |
|---|---:|---:|---:|---:|---:|---:|
'''+ '\n'.join(rows)+'''

Linux stage pairs mean rounded CPU/elapsed seconds. Windows single values remain labelled reported seconds; unrounded Rout process timings are separate in JSON. Windows totals 645.813 and 622.432 seconds exceed 600 seconds, while vignette rebuilding is 22 seconds on each. This does not establish that CRAN will accept the timings. The protocol calls 600 seconds an observed incoming signal, not an immutable policy limit. Its strict first-submission rule is distinct: the retained CRAN snapshot shows 0.10.0, so this is an update. Maintainer/final-panel timing review stays open.

The local Mac native check is `'''+local['native_state']+'''`; warning replay is `'''+local['warning_replay_state']+'''`. The known compare_methods NOTE reports user 2.532/system 0.992/elapsed 18.406 seconds (CPU 3.524 seconds). It remains unadjudicated; no native check was rerun to obtain faster timings. The native process exited 0, while the strict checker driver exited 1 because it rejects a NOTE. The measured native total is {LOCAL_TOTAL} seconds. The native tests stage reports rounded CPU/elapsed 5m/31m; its Rout reports 1855.943 elapsed seconds. Observed counts match 3,195 passes/161 warnings/86 skips/0 failures, pending consequence review. The guarded refresh helper preserves native exit versus strict checker-driver rejection and leaves any completed warning replay unadjudicated.

Native installed-size INFO reports 5.9Mb on all three remote platforms, with data 2.2Mb, extdata 1.1Mb and help 1.1Mb; the local snapshot reports 7.0Mb, with data 3.2Mb, extdata 1.1Mb and help 1.1Mb. These are raw platform observations, not identical byte-size measurements or native warnings. Archive bytes are unchanged. No package-owned src files ship; compiled dependencies remain separate requirements.

All six individual vignette diagnostics exited 0, with render-time sum 74.8865 seconds and whole-R-process sum 308.2825 seconds. The enclosing sequential loop took 310.1696 seconds, including between-call overhead. Package/runtime startup is outside the render interval. These separate local measurements do not replace native or Windows timing. Thirty-two explicit `eval = FALSE` chunks and the multiple-imputation global `eval = FALSE` setting remained disabled.

| Vignette | Render elapsed s | Whole R process s |
|---|---:|---:|
'''+ '\n'.join(vignette_rows)+'''

The isolated missing-rgbif/missing-terra proof records 60 passes, 0 failures, 0 warnings and 8 expected skips: 2 rgbif-dependent tests, 3 terra-dependent tests and 3 live/slow integrations. It exercises helpful absence branches and bounded cache/helper/failure tests. Dependency symlinks plus base R exclude user/site-library fallthrough; this is not a clean-room dependency build. Live GBIF/WorldClim availability, downloads and skipped raster/stored-point paths are not proven.

The exact-source URL check retains one distinct DOI 403 in two locations. The retained [Crossref metadata response](https://api.crossref.org/works/10.1126/science.1157704) returned 200 and identifies Hackett's 2008 Science paper, volume 320(5884), pages 1763–1768, with the exact DOI/title/date. Identifier identity is supported; automated DOI access remains 403. No URL-clean verdict is made. Initial sandbox DNS/export failures remain in original receipts.

Spelling review covers 404 flagged tokens, 1,886 contexts and 43 archive-bound source files. It identified no ordinary-word typo or uncertain spelling token among these flags. This does not establish unflagged prose, grammar, scientific claims or reference quality. The category/scope counts and original contexts remain in the spelling receipt.

Installed MI verification passed all 12 steps across neither, drm-only, gllvm-only and both optional-backend cells. Both measured backend versions are 0.7.1; their source commits are not embedded in the installed packages. Dependency symlinks and three downstream fits per real fixture bound this proof; it is not a clean-room build or a recovery campaign. Artifact source and checker commits are recorded separately.

Earlier Linux “timings absent” wording was a derived parser error and is superseded by paired timings. The original service receipt's pending broad spelling review is superseded only for its 404 flags. STATE's predecessor warning-capture and old running-build/checker wording must not be promoted into current final-artifact evidence. `contradictions_and_repairs` records these distinctions.

No panel verdict or release-rung promotion is made. Local timing NOTE/replay, overall timing and warning-consequence review remain open. PR 237 is recorded merged; PR 236 remains unmerged under retained authority. CRAN upload/confirmation are separately gated.

Refresh local receipt fields after native completion or after the separate warning replay finishes:

```
'''+s['refresh_command']+'''
```

The helper reads proofs only and preserves prior outputs in history. It refuses a changed or extra native NOTE, error/warning, hash mismatch or invalid completed replay. A later refresh changes the exact draft: repeat its assessment before describing it as assessed. Assessment completion is recorded separately, never implied by this synopsis.
''').replace('{LOCAL_TOTAL}',f"{local['total_native_check_seconds']:.3f}" if local['total_native_check_seconds'] is not None else 'PENDING')

if __name__=='__main__': main()
