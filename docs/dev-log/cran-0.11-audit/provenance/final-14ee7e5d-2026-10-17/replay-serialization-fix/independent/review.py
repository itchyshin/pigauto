#!/usr/bin/env python3
"""Independent Python-only guard review. Never runs the actual audit or R."""
import ast, contextlib, copy, hashlib, importlib.util, io, json, pathlib, sys

HERE=pathlib.Path(__file__).resolve().parent
BASE=HERE.parent.parent
AUDIT=BASE/'final-warning-audit/audit.py'
BEFORE=HERE.parent/'audit-before.py'
CONTRACT=BASE/'final-warning-audit/contract.json'
PRODUCER=pathlib.Path('/private/tmp/pigauto-pr236-validator-20261017/script/cran-0.11-platform/capture-warnings.R')
ROOT_TEST=HERE.parent/'test-replay-environment.py'
SHA=lambda p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def main():
    actual=json.loads((HERE.parent/'actual-guard-reproduction.json').read_text())
    assert len(actual)==4
    protected={AUDIT,BEFORE,CONTRACT,PRODUCER,ROOT_TEST,HERE.parent/'serializer-reproduction.log',HERE.parent/'actual-guard-reproduction.json',BASE/'FINAL.json'}
    for row in actual:
        root=pathlib.Path(row['path']).parent.parent
        protected.update(p for p in (root/'platform-evidence').rglob('*') if p.is_file())
        check=root/'frozen-check/pigauto.Rcheck'
        protected.update(check/p for p in ['00check.log','00install.out','pigauto-Ex.Rout','pigauto-Ex.timings','tests/testthat.Rout'])
    failed=BASE/'final-warning-audit/actual-results-with-native-note'
    protected.update(p for p in failed.rglob('*') if p.is_file())
    final=json.loads((BASE/'FINAL.json').read_text());protected.add(pathlib.Path(final['path']))
    before={str(p):{'sha256':SHA(p),'size_bytes':p.stat().st_size} for p in sorted(protected)}
    assert before[str(AUDIT)]['sha256']=='be66a565bafe454a471b6216e69664af5a91fc24742d4217fb2abedd584f31a4'
    assert before[str(BEFORE)]['sha256']=='054b343e920c7517b8a3906ef80f8fa8519ec141fdb713f81000112c7ed316d7'
    # Whole-module AST equivalence except the new function and one replacement call.
    old=ast.parse(BEFORE.read_text());new=ast.parse(AUDIT.read_text());replacements=0
    for fn in old.body:
        if isinstance(fn,ast.FunctionDef) and fn.name=='audit_cell':
            for node in ast.walk(fn):
                if isinstance(node,ast.Try):
                    for i,statement in enumerate(node.body):
                        if isinstance(statement,ast.Expr) and ast.unparse(statement).startswith("require(replay.get('error') is None"):
                            node.body[i]=ast.parse('validate_replay_environment(replay)').body[0];replacements+=1
    assert replacements==1
    additions=[node for node in new.body if isinstance(node,ast.FunctionDef) and node.name=='validate_replay_environment']
    assert len(additions)==1
    new.body=[node for node in new.body if node not in additions]
    assert ast.dump(old,include_attributes=False)==ast.dump(new,include_attributes=False), 'Unexpected audit behavior changed'
    spec=importlib.util.spec_from_file_location('independent_warning_audit',AUDIT)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    # Fail if any test accidentally tries subprocess/R work.
    original_run=module.subprocess.run;original_popen=module.subprocess.Popen
    def forbidden(*a,**k):raise AssertionError('No subprocess or actual audit authorized in this review')
    module.subprocess.run=forbidden;module.subprocess.Popen=forbidden
    results=[]
    def exercise(label,value,expected):
        # Domain is standard JSON-decoded receipts, not arbitrary Python subclasses.
        value=json.loads(json.dumps(value,allow_nan=False))
        try:module.validate_replay_environment(value);accepted=True;error=None
        except ValueError as e:accepted=False;error=str(e)
        assert accepted==expected,(label,value,accepted,expected)
        results.append({'case':label,'input':value,'expected_accept':expected,'accepted':accepted,'refusal':error})
    valid={'error':None,'not_cran':'UNSET','force_suggests':'true'}
    errors=[None,{},[],False,True,0,0.0,1,-1,'',' ','null','{}','planted error',{'message':None},{'message':''},{'message':False},{'message':'failure'},{'':None},{'nested':{}},[None],[{}]]
    envs=[None,False,True,0,1,{},[],'',' ','TRUE','true','True','UNSET','unset','UNSET ','UNSET\n','true ','true\n','false','UNS\u0395T']
    for i,error in enumerate(errors):
        canonical=json.dumps(error,separators=(',',':'))
        exercise('error-domain-'+str(i),valid|{'error':error},canonical in ['null','{}'])
    for field,expected in [('not_cran','UNSET'),('force_suggests','true')]:
        for i,value in enumerate(envs):exercise(field+'-type/value-'+str(i),valid|{field:value},type(value) is str and value==expected)
        exercise('missing-'+field,{k:v for k,v in valid.items() if k!=field},False)
    exercise('missing-error',{k:v for k,v in valid.items() if k!='error'},False)
    for value in [None,False,True,0,1,'',[],[valid],{}]:exercise('wrong-envelope-'+repr(value),value,False)
    # A valid error encoding must never weaken either environment gate.
    for error in [None,{}]:
        for left in ['UNSET','true',None,False]:
            for right in ['true','false',None,True]:
                exercise('cross-environment-'+repr((error,left,right)),{'error':error,'not_cran':left,'force_suggests':right},left=='UNSET' and right=='true')
    actual_records=[]
    for i,row in enumerate(actual):
        replay=json.loads(pathlib.Path(row['path']).read_text())
        assert replay['error']=={} and replay['not_cran']=='UNSET' and replay['force_suggests']=='true'
        assert (replay.get('error') is None and replay['not_cran']=='UNSET' and replay['force_suggests']=='true') is False
        module.validate_replay_environment(replay)
        actual_records.append({'path':row['path'],'sha256':SHA(row['path']),'error':replay['error'],'new_environment_guard':'PASS','old_environment_guard':'FAIL','native_counts':replay['native_counts'],'replay_counts':replay['replay_counts']})
        for value in [{'message':None},[],False,0,'','null']:
            exercise('actual-'+str(i)+'-planted-error-'+repr(value),replay|{'error':value},False)
        for field in ['error','not_cran','force_suggests']:
            exercise('actual-'+str(i)+'-missing-'+field,{k:v for k,v in replay.items() if k!=field},False)
    stream=io.StringIO()
    with contextlib.redirect_stdout(stream):module.selftest()
    module.subprocess.run=original_run;module.subprocess.Popen=original_popen
    assert 'SELFTEST_PASS' in stream.getvalue() and 'TIMING_NOTE_SELFTEST_PASS' in stream.getvalue()
    after={str(p):{'sha256':SHA(p),'size_bytes':p.stat().st_size} for p in sorted(protected)}
    assert before==after,'Contract, producer, failed audit or raw inputs changed'
    assert SHA(final['path'])==final['sha256'] and pathlib.Path(final['path']).stat().st_size==5137648
    (HERE/'controls.json').write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n')
    (HERE/'existing-selftest.log').write_text(stream.getvalue())
    (HERE/'protected-inputs-before.json').write_text(json.dumps(before,indent=2,sort_keys=True)+'\n')
    (HERE/'protected-inputs-after.json').write_text(json.dumps(after,indent=2,sort_keys=True)+'\n')
    result={'verdict':'PASS_NARROW_REPLAY_ENVIRONMENT_VALIDATOR_FIX','audit_before_sha256':SHA(BEFORE),'audit_after_sha256':SHA(AUDIT),'producer_sha256':SHA(PRODUCER),'contract_sha256':SHA(CONTRACT),'ast_equivalence':'Entire old module preserved except one inline guard replaced by new function; only new top-level definition is validate_replay_environment.','new_guard_controls':len(results),'positive_controls':sum(r['accepted'] for r in results),'negative_controls_refused':sum(not r['accepted'] for r in results),'actual_replay_receipts':actual_records,'existing_pure_python_selftest':'PASS; subprocess creation disabled during review','protected_inputs_unchanged':True,'protected_file_count':len(before),'python':sys.version,'archive_sha256':final['sha256'],'archive_bytes':5137648,'limits':['No actual four-cell audit or RDS reader was run.','Acceptance of null/empty-object encoding establishes the environment envelope only; raw problems/counts/provenance/identity checks must still pass.','Local timing NOTE stays UNADJUDICATED; this is not a final release panel or READY vote.','Standard json.loads value domain only; no general strict-JSON duplicate-key or arbitrary Python-subclass guarantees.']}
    (HERE/'review-result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ['verdict','new_guard_controls','positive_controls','negative_controls_refused','protected_inputs_unchanged','protected_file_count']}))

if __name__=='__main__':main()
