import re,json,hashlib,csv,io
from pathlib import Path
STAGE_PHRASES=[('examples','examples'),('tests','tests'),('vignettes','re-building of vignette outputs'),('PDF','PDF version of manual'),('HTML','HTML version of manual')]
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
S=Path(__file__).resolve().parent
F=json.loads((S/'FINAL.json').read_text())
controls=[]
for raw,want in [('* checking examples ... [14s/18s] OK',{'raw_check_line':'* checking examples ... [14s/18s] OK','timing_format':'cpu_elapsed_seconds','cpu_seconds':14.0,'elapsed_seconds':18.0}),('* checking examples ... [16s] OK',{'raw_check_line':'* checking examples ... [16s] OK','timing_format':'single_reported_seconds','reported_seconds':16.0}),('* checking examples ... OK',{'raw_check_line':'* checking examples ... OK','timing_format':'not_reported'})]:
 assert parse_stage_timing(raw)==want;controls.append('PASS')
for token in ['[14s/14s/1s]','[14/14s]','[14ms]','[NaNs]','[14s elapsed]','[14s / 14s]','[14s] NOTE']:
 raw='* checking examples ... '+token+('' if token.endswith('NOTE') else ' OK')
 try:parse_stage_timing(raw)
 except ValueError:controls.append('REFUSED_EXPECTED')
 else:raise AssertionError(raw)
assert len(controls)==10
rows=[]
paths=[('macos-arm64',S/'native-macos')]+[(d.name.removeprefix('frozen-a93fece7-'),d) for d in sorted((S/'direct-platform').glob('frozen-*')) if d.is_dir()]
for label,d in paths:
 e=d/'platform-evidence';p=d/'frozen-check/pigauto.Rcheck/00check.log';text=p.read_text(errors='replace');assert 'Status: OK' in text.splitlines()
 deps=list(csv.DictReader(io.StringIO((e/'declared-dependencies.csv').read_text())));assert len(deps)==29 and all(x['available']=='TRUE' for x in deps)
 exittext=(e/'check-exit.txt').read_text().splitlines();assert 'exit=0' in exittext;elapsed=float(next(x.split('=',1)[1] for x in exittext if x.startswith('elapsed_seconds=')))
 env=(e/'check-environment.txt').read_text();assert 'NOT_CRAN=UNSET' in env and '_R_CHECK_FORCE_SUGGESTS_=true' in env
 rout=(d/'frozen-check/pigauto.Rcheck/tests/testthat.Rout').read_text(errors='replace');match=re.findall(r'FAIL\s+([0-9]+)\s*\|\s*WARN\s+([0-9]+)\s*\|\s*SKIP\s+([0-9]+)\s*\|\s*PASS\s+([0-9]+)',rout);assert match and len(set(match))==1,match
 counts=dict(zip(['fail','warnings','skips','pass'],map(int,match[0])));assert counts['fail']==0
 rows.append({'platform':label,'artifact_sha256':F['sha256'],'source_commit':F['source_commit'],'status':'OK','native_elapsed_seconds':elapsed,'stage_timings':parse_native_stages(text),'native_counts':counts,'declared_dependencies_available':len(deps),'check_log_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'native_rout_sha256':hashlib.sha256((d/'frozen-check/pigauto.Rcheck/tests/testthat.Rout').read_bytes()).hexdigest()})
assert len(rows)==4
(S/'native-platform-summary.json').write_text(json.dumps({'artifact_sha256':F['sha256'],'parser_controls':controls,'platforms':rows,'scope':'Original gzip native checks; rounded stage pairs preserve CPU and elapsed separately; single reported seconds are not assigned a CPU/elapsed interpretation. Separate warning replays excluded from native duration.'},indent=2)+'\n')
print('NATIVE4_STATUS_OK_ALL29_DEPENDENCIES_PARSER10_CONTROLS_PASS')
