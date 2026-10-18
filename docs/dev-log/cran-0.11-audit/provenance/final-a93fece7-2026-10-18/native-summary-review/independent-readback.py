from pathlib import Path
import ast, csv, hashlib, io, json, re, subprocess, sys

S = Path('/private/tmp/pigauto-final-closure-20261017/successor-b2e6054')
OUT = S / 'native-summary-review'
FINAL = json.loads((S / 'FINAL.json').read_text())
SUMMARY_PATH = S / 'native-platform-summary.json'
SUMMARY = json.loads(SUMMARY_PATH.read_text())
EXPECTED_SHA = 'a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb'
EXPECTED_COMMIT = 'b2e6054685e6bdedb073ba46e058af2c0a619a18'
STAGES = [
    ('examples', 'examples'), ('tests', 'tests'),
    ('vignettes', 're-building of vignette outputs'),
    ('PDF', 'PDF version of manual'), ('HTML', 'HTML version of manual')]
GIT = '/Library/Developer/CommandLineTools/usr/bin/git'
SOURCE = Path(FINAL['source_checkout'])
CHECKER = Path('/private/tmp/pigauto-pr236-validator-20261017/script/cran-0.11-platform/check.R')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_line(line):
    m = re.fullmatch(r'\* checking .+ \.\.\. (?:\[([0-9]+(?:\.[0-9]+)?)s(?:/([0-9]+(?:\.[0-9]+)?)s)?\] )?OK', line)
    if not m:
        raise ValueError('unrecognized native timing line: ' + line)
    if m.group(1) is None:
        return {'raw_check_line': line, 'timing_format': 'not_reported'}
    if m.group(2) is None:
        return {'raw_check_line': line, 'timing_format': 'single_reported_seconds',
                'reported_seconds': float(m.group(1))}
    return {'raw_check_line': line, 'timing_format': 'cpu_elapsed_seconds',
            'cpu_seconds': float(m.group(1)), 'elapsed_seconds': float(m.group(2))}


def parse_stages(log_text):
    found = {}
    for key, phrase in STAGES:
        rows = [line for line in log_text.splitlines() if line.startswith('* checking ' + phrase + ' ...')]
        if len(rows) != 1:
            raise ValueError(f'{phrase}: expected one stage row, saw {len(rows)}')
        found[key] = parse_line(rows[0])
    return found


def resolve_counts(raw_matches):
    if not raw_matches:
        raise ValueError('no native test-summary counts')
    unique = set(tuple(x) for x in raw_matches)
    if len(unique) != 1:
        raise ValueError(f'conflicting native test-summary counts: {sorted(unique)}')
    return tuple(raw_matches[0])


def function_ast(path, names):
    tree = ast.parse(Path(path).read_text())
    nodes = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    return {name: ast.dump(nodes[name], include_attributes=False) for name in names}

# Independent known-answer controls; these records are parser validation only.
controls = []
for raw, expected in [
    ('* checking examples ... [14s/18s] OK', ('cpu_elapsed_seconds', 14.0, 18.0)),
    ('* checking examples ... [16s] OK', ('single_reported_seconds', 16.0, None)),
    ('* checking examples ... OK', ('not_reported', None, None)),
]:
    got = parse_line(raw)
    val = (got['timing_format'], got.get('cpu_seconds', got.get('reported_seconds')),
           got.get('elapsed_seconds'))
    assert val == expected, (raw, val, expected)
    controls.append({'kind': 'timing_positive', 'input': raw, 'result': 'ACCEPT_EXPECTED'})
for raw in ['[14s/14s/1s]', '[14/14s]', '[14ms]', '[NaNs]', '[14s elapsed]',
            '[14s / 14s]', '[14s] NOTE']:
    line = '* checking examples ... ' + raw + ('' if raw.endswith('NOTE') else ' OK')
    try:
        parse_line(line)
    except ValueError:
        controls.append({'kind': 'timing_negative', 'input': line, 'result': 'REFUSE_EXPECTED'})
    else:
        raise AssertionError('malformed timing accepted: ' + line)
assert len(controls) == 10
same = resolve_counts([('0','0','0','2'), ('0','0','0','2')])
controls.append({'kind':'duplicate_counts_same','input_count':2,'result':'ACCEPT_EXPECTED','counts':same})
for label, values in [('duplicate_counts_conflict',[('0','0','0','2'),('0','1','0','2')]),
                      ('counts_absent',[])]:
    try:
        resolve_counts(values)
    except ValueError:
        controls.append({'kind':label,'result':'REFUSE_EXPECTED'})
    else:
        raise AssertionError(label + ' was accepted')

assert FINAL['sha256'] == EXPECTED_SHA and FINAL['source_commit'] == EXPECTED_COMMIT
archive = Path(FINAL['path'])
assert digest(archive) == EXPECTED_SHA and archive.stat().st_size == FINAL['size_bytes']
assert subprocess.check_output([GIT, '-C', str(SOURCE), 'rev-parse', 'HEAD'], text=True).strip() == EXPECTED_COMMIT
assert subprocess.check_output([GIT, '-C', str(SOURCE), 'status', '--porcelain'], text=True).strip() == ''
assert digest(CHECKER) == '1e641c87e699fba2846aaa9c8958fa77348ecd22221b5ed08e1ccf31c2de9fb0'

platform_dirs = [
    ('macos-arm64', S / 'native-macos'),
    ('ubuntu-latest-release', S / 'direct-platform/frozen-a93fece7-ubuntu-latest-release'),
    ('windows-latest-devel', S / 'direct-platform/frozen-a93fece7-windows-latest-devel'),
    ('windows-latest-release', S / 'direct-platform/frozen-a93fece7-windows-latest-release'),
]
summary_rows = {row['platform']: row for row in SUMMARY['platforms']}
assert SUMMARY['artifact_sha256'] == EXPECTED_SHA and len(summary_rows) == 4
results = []
for label, d in platform_dirs:
    assert label in summary_rows
    row = summary_rows[label]
    log_path = d / 'frozen-check/pigauto.Rcheck/00check.log'
    rout_path = d / 'frozen-check/pigauto.Rcheck/tests/testthat.Rout'
    e = d / 'platform-evidence'
    log_bytes = log_path.read_bytes(); log = log_bytes.decode(errors='replace')
    assert digest(log_path) == row['check_log_sha256']
    assert digest(rout_path) == row['native_rout_sha256']
    assert 'Status: OK' in log.splitlines()
    actual_stages = parse_stages(log)
    assert actual_stages == row['stage_timings'], (label, actual_stages, row['stage_timings'])
    exit_lines = (e / 'check-exit.txt').read_text().splitlines()
    assert 'exit=0' in exit_lines
    elapsed = float(next(line.split('=',1)[1] for line in exit_lines if line.startswith('elapsed_seconds=')))
    assert abs(elapsed - row['native_elapsed_seconds']) < 1e-6
    deps = list(csv.DictReader(io.StringIO((e / 'declared-dependencies.csv').read_text())))
    assert len(deps) == row['declared_dependencies_available'] == 29
    assert all(dep['available'] == 'TRUE' for dep in deps)
    dep_names = sorted(dep['package'] for dep in deps)
    env = (e / 'check-environment.txt').read_text()
    assert 'NOT_CRAN=UNSET' in env and '_R_CHECK_FORCE_SUGGESTS_=true' in env
    matches = re.findall(r'FAIL\s+([0-9]+)\s*\|\s*WARN\s+([0-9]+)\s*\|\s*SKIP\s+([0-9]+)\s*\|\s*PASS\s+([0-9]+)', rout_path.read_text(errors='replace'))
    counts = resolve_counts(matches)
    count_dict = dict(zip(['fail','warnings','skips','pass'], map(int, counts)))
    assert count_dict == row['native_counts'] and count_dict['fail'] == 0
    evidence = {'platform': label, 'check_log_sha256': digest(log_path),
                'testthat_rout_sha256': digest(rout_path),
                'check_elapsed_seconds': elapsed, 'native_stage_timings': actual_stages,
                'dependency_count': len(deps), 'dependency_names_sha256': hashlib.sha256(('\n'.join(dep_names)+'\n').encode()).hexdigest(),
                'test_summary_matches': count_dict, 'test_summary_occurrences': len(matches),
                'identical_count_records': len(set(matches)) == 1,
                'check_environment_not_cran_unset': True}
    if label == 'macos-arm64':
        rec = json.loads((d / 'native-receipt.json').read_text())
        assert rec['artifact_sha256'] == EXPECTED_SHA and rec['source_commit'] == EXPECTED_COMMIT
        assert rec['native_check_script_sha256'] == digest(CHECKER)
        assert rec['archive_sha256_before'] == rec['archive_sha256_after'] == EXPECTED_SHA
        assert rec['exit_code'] == 0 and not rec['timed_out']
        evidence.update({'whole_runner_elapsed_seconds': rec['elapsed_seconds'],
                         'native_check_script_sha256': rec['native_check_script_sha256'],
                         'whole_runner_receipt_sha256': digest(d / 'native-receipt.json')})
    else:
        pre = json.loads((e / 'identity-pre.json').read_text())
        post = json.loads((e / 'identity-post.json').read_text())
        assert pre['sha256'] == post['sha256'] == EXPECTED_SHA
        assert pre['size_bytes'] == post['size_bytes'] == FINAL['size_bytes']
        assert pre['inventory_sha256'] == post['inventory_sha256'] == FINAL['inventory_sha256']
        assert pre['archive_entries'] == post['archive_entries'] == FINAL['archive_entries']
        assert str(pre['run_id']) == str(post['run_id']) == '39091958419'
        evidence.update({'archive_pre_sha256': pre['sha256'], 'archive_post_sha256': post['sha256'],
                         'workflow_source_sha': pre['workflow_sha'], 'run_id': str(pre['run_id']),
                         'artifact_pre_receipt_sha256': digest(e / 'identity-pre.json'),
                         'artifact_post_receipt_sha256': digest(e / 'identity-post.json')})
    results.append(evidence)
assert len(results) == 4
# Compare name sets directly, rather than depending on row ordering or available flags.
dep_sets=[]
for _,d in platform_dirs:
    dep_sets.append(tuple(sorted(x['package'] for x in csv.DictReader(io.StringIO((d/'platform-evidence/declared-dependencies.csv').read_text())))))
assert len(set(dep_sets)) == 1

# Verify root parser functions were not edited from the reviewed pre-duplicate-control copy.
funcs = ['parse_stage_timing','parse_native_stages']
current_ast = function_ast(S / 'summarise-native.py', funcs)
prior_ast = function_ast(S / 'summarise-native-before-duplicate-summary-guard.py', funcs)
assert current_ast == prior_ast

# Verify each downloaded CI artifact container independently from the package archive inside it.
download_binding = json.loads((S / 'direct-platform/download-binding.json').read_text())
zip_records = []
for entry in download_binding['records']:
    zip_path = S / 'direct-platform' / (entry['name'] + '.zip')
    assert digest(zip_path) == entry['local_sha256']
    assert entry['server_digest'] == 'sha256:' + digest(zip_path)
    assert zip_path.stat().st_size == entry['size_bytes']
    zip_records.append({'name': entry['name'], 'zip_sha256': digest(zip_path),
                        'size_bytes': entry['size_bytes'], 'run_id': entry['run'],
                        'artifact_head_sha': entry['head_sha']})
assert len({x['zip_sha256'] for x in zip_records}) == 3

# The five distinct identity layers: release source, frozen bytes, checker implementation, workflow/run, and synthetic parser controls.
mac_receipt=json.loads((S/'native-macos/native-receipt.json').read_text())
workflow_shas={r['workflow_source_sha'] for r in results if r['platform']!='macos-arm64'}
assert len(workflow_shas)==1
workflow_sha=next(iter(workflow_shas))
head_shas={r['head_sha'] for r in json.loads((S/'direct-platform/download-binding.json').read_text())['records']}
assert len(head_shas)==1
head_sha=next(iter(head_shas))
identities={'source_commit':EXPECTED_COMMIT,'archive_sha256':EXPECTED_SHA,
            'native_checker_script_sha256':digest(CHECKER),'workflow_source_sha':workflow_sha,
            'platform_artifact_head_sha':head_sha,'check_run_id':39091958419,
            'ci_artifact_container_sha256s':[x['zip_sha256'] for x in zip_records],
            'synthetic_control_record_sha256':hashlib.sha256(json.dumps(controls,sort_keys=True,separators=(',',':')).encode()).hexdigest()}
assert len({EXPECTED_COMMIT, EXPECTED_SHA, digest(CHECKER), workflow_sha, head_sha}) == 5

report={'reviewer':'pat','verdict':'PASS','scope':'Independent raw readback of four native checks and summary parser controls; no native check rerun.',
        'source_commit':EXPECTED_COMMIT,'source_worktree_clean':True,
        'archive_sha256':EXPECTED_SHA,'archive_size_bytes':FINAL['size_bytes'],
        'native_checker_script_path':str(CHECKER),'native_checker_script_sha256':digest(CHECKER),
        'workflow_source_sha':workflow_sha,'platform_artifact_head_sha':head_sha,'direct_platform_run_id':39091958419,
        'identity_layers':identities,'root_summary_sha256':digest(SUMMARY_PATH),
        'root_summarizer_sha256':digest(S/'summarise-native.py'),
        'prior_summarizer_sha256':digest(S/'summarise-native-before-duplicate-summary-guard.py'),
        'parser_functions_ast_identical_to_prior_reviewed_file':True,
        'synthetic_controls':controls,
        'platforms':results,'downloaded_ci_artifact_containers':zip_records,
        'dependency_package_name_sets_identical_across_four':True,
        'blockers':[],
        'limits':['This readback confirms Status OK native outputs and their local counts/timings; it does not waive or adjudicate warnings.','It does not cover CRAN upload, external confirmation, acceptance, or final release approval.'],
        'reviewed_evidence':[str(SUMMARY_PATH),str(S/'summarise-native.py'),str(S/'summarise-native-before-duplicate-summary-guard.py'),str(S/'FINAL.json'),str(CHECKER),str(S/'native-macos/native-receipt.json'),str(S/'direct-platform/download-binding.json')],
}
(OUT/'actual-readback.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'verdict':report['verdict'],'identity_layers':identities,
                  'synthetic_controls':len(controls),'platforms':results,
                  'report':str(OUT/'actual-readback.json')},indent=2))
