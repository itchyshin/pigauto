#!/usr/bin/env python3
import pathlib,runpy,json,sys,hashlib

def run_controls(root):
    root=pathlib.Path(root);parser=runpy.run_path(str(root/'bundle.py'));parse=parser['parse_stage_timing'];results=[]
    cases=[('* checking examples ... [14s/18s] OK',{'timing_format':'cpu_elapsed_seconds','cpu_seconds':14.0,'elapsed_seconds':18.0}),('* checking examples ... [16s] OK',{'timing_format':'single_reported_seconds','reported_seconds':16.0}),('* checking examples ... OK',{'timing_format':'not_reported'})]
    for raw,expected in cases:
        actual=parse(raw);assert actual==dict(expected,raw_check_line=raw);results.append({'input':raw,'expected':expected,'result':'PASS'})
    malformed=['[14s/14s/1s]','[14/14s]','[14ms]','[NaNs]','[14s elapsed]','[14s / 14s]','[14s] NOTE']
    for token in malformed:
        line='* checking examples ... '+token+('' if token.endswith('NOTE') else ' OK')
        try:parse(line)
        except ValueError:results.append({'input':line,'result':'REFUSED_EXPECTED'})
        else:raise AssertionError('Unrecognized format silently accepted')
    summary=json.loads((root/'platform-summary.json').read_text());reference=json.loads((root/'timing-parser-independent-reference.json').read_text());actual=0
    for row in summary:
        independent=reference['frozen-14ee7e5d-'+row['platform']];assert row['check_log_sha256']==independent['check_log_sha256']
        for rawstage in independent['stages']:
            parsed=parse(rawstage['raw_line']);stage=next(v for v in row['stage_timings'].values() if isinstance(v,dict) and v.get('raw_check_line')==rawstage['raw_line']);assert parsed==stage
            assert {k:v for k,v in parsed.items() if k not in ['raw_check_line','timing_format']}==rawstage['timing'];actual+=1
    assert actual==15
    return {'status':'PASS_TIMING_PARSER_KNOWN_ANSWER_AND_INDEPENDENT_REFERENCE','controls':results,'actual_stages_checked':actual,'helper_sha256':hashlib.sha256((root/'bundle.py').read_bytes()).hexdigest(),'scope':'Rounded native CPU/elapsed pairs and single reported seconds are distinct from unrounded Rout process times.'}

if __name__=='__main__':
    root=pathlib.Path(sys.argv[1]);result=run_controls(root);(root/'timing-parser-controls.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':result['status'],'controls':len(result['controls']),'actual_stages_checked':result['actual_stages_checked']}))
