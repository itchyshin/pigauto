#!/usr/bin/env python3
import copy, hashlib, json, math, pathlib, runpy
ROOT=pathlib.Path(__file__).resolve().parent
SHA='14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb'

def validate(d):
    assert d['artifact']['sha256']==SHA and d['artifact']['size_bytes']==5137648
    assert d['overall'].startswith('NOT_READY')
    assert len(d['remote_native_checks'])==3
    expected={'examples':14,'tests':206,'vignettes':11,'PDF':2,'HTML':2}
    for r in d['remote_native_checks']:
        assert r['native_status']=='OK' and r['all_declared_dependencies_available']==29
        assert r['native_counts']=={'pass':3195,'fail':0,'warnings':161,'skips':86}
        if r['platform'].startswith('windows'):
            assert r['native_total_seconds']>600 and r['stage_timings']['vignettes']['reported_seconds']==22
        else:
            for k,value in expected.items():
                assert r['stage_timings'][k]['cpu_seconds']==r['stage_timings'][k]['elapsed_seconds']==value
    v=d['individual_vignette_diagnostic'];assert len(v['rows'])==6 and all(r['exit']==0 for r in v['rows'])
    assert math.isclose(sum(r['render_elapsed_seconds'] for r in v['rows']),v['render_elapsed_sum_seconds'])
    assert math.isclose(sum(r['process_elapsed_seconds'] for r in v['rows']),v['whole_R_process_elapsed_sum_seconds'])
    assert v['whole_R_process_elapsed_sum_seconds']<v['enclosing_loop_elapsed_seconds']
    assert d['urls']['finding']['Status']=='403' and d['urls']['primary_metadata']['http_status']==200
    assert d['external_service_absence']['counts']['pass']==60 and d['external_service_absence']['counts']['skip']==8
    assert d['spelling']['flag_rows']==404 and d['spelling_archive_source_files_verified']==43
    assert d['installed_MI']['steps_passed']==12 and {c['cell'] for c in d['installed_MI']['cells']}=={'neither','drm-only','gllvm-only','both'}
    assert d['local_macos']['known_note']['elapsed_seconds']==18.406
    assert d['local_macos']['known_note']['disposition'].startswith('UNADJUDICATED')
    assert d['warnings']['consequences']=='PENDING_FINAL_WARNING_AUDIT'

def main():
    d=json.loads((ROOT/'synopsis.json').read_text());validate(d)
    inputs=json.loads((ROOT/'inputs.json').read_text())
    for name,row in inputs.items():
        p=ROOT/row['snapshot'] if row['mutable_observation'] else pathlib.Path(name)
        raw=p.read_bytes();assert hashlib.sha256(raw).hexdigest()==row['sha256'] and len(raw)==row['size_bytes'],name
    helper=runpy.run_path(str(ROOT/'build-synopsis.py'));parse=helper['timing_line'];positive=0
    for raw,expected in [('* checking tests ... [5m/31m] OK',{'cpu_seconds':300,'elapsed_seconds':1860}),('* checking examples ... [14s/18s] OK',{'cpu_seconds':14,'elapsed_seconds':18}),('* checking examples ... [13s/63s] NOTE',{'cpu_seconds':13,'elapsed_seconds':63}),('* checking tests ... [4m/241s] OK',{'cpu_seconds':240,'elapsed_seconds':241}),('* checking examples ... [17s] OK',{'reported_seconds':17})]:
        got=parse(raw)
        for key,value in expected.items():assert got[key]==value
        positive+=1
    negative=[]
    for raw in ['* checking tests ... [5m/31m/1s] OK','* checking tests ... [5/31m] OK','* checking tests ... [5minutes] OK','* checking tests ... [5m] ERROR']:
        try:parse(raw)
        except ValueError:negative.append('refused malformed native timing')
        else:raise AssertionError('Malformed timing accepted')
    for kind in ['artifact','linux_pair','process_sum','local_note','release_claim']:
        bad=copy.deepcopy(d)
        if kind=='artifact':bad['artifact']['sha256']='0'*64
        if kind=='linux_pair':next(r for r in bad['remote_native_checks'] if r['platform'].startswith('ubuntu'))['stage_timings']['examples']['cpu_seconds']=15
        if kind=='process_sum':bad['individual_vignette_diagnostic']['whole_R_process_elapsed_sum_seconds']=310.1696
        if kind=='local_note':bad['local_macos']['known_note']['elapsed_seconds']=0
        if kind=='release_claim':bad['overall']='READY'
        try:validate(bad)
        except AssertionError:negative.append('refused '+kind)
        else:raise AssertionError('Planted synopsis error accepted: '+kind)
    assert len(negative)==9
    result={'result':'PASS_BOUNDED_SYNOPSIS_FIGURES_AND_INPUT_BINDINGS','input_receipts_verified':len(inputs),'positive_timing_controls':positive,'negative_controls':negative,'local_state':d['local_macos']['native_state'],'local_replay':d['local_macos']['warning_replay_state'],'no_release_promotion':True,'draft_sha256':hashlib.sha256((ROOT/'README.md').read_bytes()).hexdigest()}
    (ROOT/'verification.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    files={str(p.relative_to(ROOT)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'size_bytes':p.stat().st_size} for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='output-manifest.json' and 'history' not in p.relative_to(ROOT).parts and 'attempts' not in p.relative_to(ROOT).parts}
    (ROOT/'output-manifest.json').write_text(json.dumps({'scope':'Current synopsis/helper/assessment/local snapshots; historical and failed attempts retained separately','files':files},indent=2,sort_keys=True)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
