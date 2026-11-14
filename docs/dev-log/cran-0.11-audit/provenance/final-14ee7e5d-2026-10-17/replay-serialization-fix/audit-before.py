"""Read-only input audit; writes only beneath this dedicated audit directory."""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
SHA = "14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb"
SIZE = 5137648
RUN = "39082452558"
PRIOR = Path('/private/tmp/pigauto-pr236-validator-20261017/docs/dev-log/cran-0.11-audit/provenance/submission-preparation-2026-10-17')
PREFIXES = {
    'small_validation': 'Small validation set',
    'conformal_sample_limit': '95% split-conformal coverage is NOT achievable',
    'constant_matched_trait': 'At least one matched trait is constant.',
    'few_tree_draws': 'Few tree-sensitivity draws',
    'covariates_joint_baseline': "'covariates' are ignored by the joint MVN baseline",
    'covariates_gnn_off': 'gnn = FALSE (the default): user covariates',
    'imbalanced_ordinal': 'Imbalanced K>=3 ordinal trait',
    'signal_gate_fallback': 'phylo_signal_gate triggered',
    'species_dropped': '3 species in traits not found in tree -- dropped.',
}
CONSEQUENCES = {
    'small_validation': 'Small calibration samples leave gates/scores noisy; no coverage certification.',
    'conformal_sample_limit': 'The stated finite-sample size cannot attain 95% split-conformal coverage.',
    'constant_matched_trait': 'The deliberately constant matched trait does not identify phylogenetic signal.',
    'few_tree_draws': 'Few tree draws can make tree-sensitivity summaries unstable.',
    'covariates_joint_baseline': 'Covariates do not contribute through the selected joint-baseline route.',
    'covariates_gnn_off': 'Covariates do not contribute when the GNN is disabled.',
    'imbalanced_ordinal': 'Ordinal class imbalance limits the tested prediction/calibration regime.',
    'signal_gate_fallback': 'The opt-in signal gate selected its documented BM fallback.',
    'species_dropped': 'The three deliberately unmatched species are dropped and not modeled.',
    'undefined_constant_correlation': 'Correlation is undefined for the constant values in this fixture.',
}

def require(condition, message):
    if not condition:
        raise ValueError(message)

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def normalized(value):
    return ' '.join(str(value).split())

def key(row, warning=True):
    fields = ('message', 'test', 'file', 'line', 'last_call') if warning else ('message', 'test', 'file', 'line')
    return tuple(normalized(row.get(field)) for field in fields)

def category(row):
    message = row['message']
    if message == 'the standard deviation is zero':
        require(row['last_call'].startswith('stats::cor('), 'Unexpected constant-correlation caller')
        return 'undefined_constant_correlation'
    for name, prefix in PREFIXES.items():
        if message.startswith(prefix):
            return name
    raise ValueError('Unknown warning: ' + message)

def read_rds(path):
    env = dict(os.environ)
    for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[variable] = '1'
    result = subprocess.run(['Rscript', '--vanilla', str(HERE/'read_objects.R'), str(path)],
                            capture_output=True, text=True, timeout=120, env=env)
    require(result.returncode == 0, 'RDS read failed: ' + result.stderr)
    return json.loads(result.stdout), result.stderr

def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))

def prepare():
    """Freeze the already adjudicated predecessor messages as review contracts."""
    baseline = rows(PRIOR/'warning-classification.csv')
    for row in baseline:
        require(category(row) == row['category'], 'Prior classifier/category mismatch')
    prior_rds = TASK/'warning-objects.rds'
    prior_receipt = json.loads((PRIOR/'warning-replay.json').read_text())
    require(prior_receipt['archive_sha256'] == 'c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2',
            'The comparison contract must come from the adjudicated predecessor')
    require(gzip.decompress((PRIOR/'warning-objects.rds.gz').read_bytes()) == prior_rds.read_bytes(),
            'Task-local prior RDS differs from the preserved, reviewed raw RDS')
    objects, stderr = read_rds(prior_rds)
    require(len(objects['warnings']) == len(baseline), 'Prior raw/CSV size mismatch')
    require([key(row) for row in baseline] == [key(row) for row in objects['warnings']],
            'Prior raw messages/provenance do not match the approved CSV')
    contract = {'basis': 'Adjudicated predecessor identities are a strict comparison contract, not final-archive evidence.',
                'baseline_csv': str(PRIOR/'warning-classification.csv'),
                'baseline_csv_sha256': digest(PRIOR/'warning-classification.csv'),
                'baseline_rds': str(prior_rds), 'baseline_rds_sha256': digest(prior_rds),
                'warnings': baseline, 'skips': objects['skips'], 'reader_stderr': stderr}
    (HERE/'contract.json').write_text(json.dumps(contract, indent=2)+'\n')
    print('CONTRACT_PREPARED', len(baseline), 'warnings;', len(objects['skips']), 'expected skips')

def native_counts(text):
    matches = re.findall(r'\[ FAIL (\d+) \| WARN (\d+) \| SKIP (\d+) \| PASS (\d+) \]', text)
    require(bool(matches), 'Native test summary absent')
    return dict(zip(('fail', 'warnings', 'skips', 'pass'), map(int, matches[-1])))

TIMING_NOTE_HEADER = '* checking examples ... [13s/63s] NOTE'
TIMING_NOTE_BODY = """Examples with CPU (user + system) or elapsed time > 5s
                 user system elapsed
compare_methods 2.532  0.992  18.406"""

def native_verdict(check_log, exit_text, remote, allow_local_timing_note=False):
    """Admit one exact local diagnostic case, retaining an unadjudicated native NOTE."""
    lines = check_log.splitlines()
    require([line for line in exit_text.splitlines() if line.startswith('exit=')] == ['exit=0'], 'Native check exit was not exactly zero')
    statuses = [line for line in lines if line.startswith('Status:')]
    findings = [line for line in lines if re.match(r'^\* .*\b(?:NOTE|WARNING|ERROR)$', line)]
    if statuses == ['Status: OK']:
        require(not findings, 'Status OK contradicts native findings')
        return {'status': 'OK', 'adjudication': 'native_Status_OK', 'exit': 0}
    require(not remote and allow_local_timing_note, 'Native check is not Status OK; local exact timing-note opt-in absent or remote cell')
    require(statuses == ['Status: 1 NOTE'], 'Native status is not exactly one NOTE with no errors/warnings')
    require(findings == [TIMING_NOTE_HEADER],
            'Unexpected native NOTE/method or error/warning finding')
    start = lines.index(findings[0]) + 1
    end = next((i for i in range(start, len(lines)) if lines[i].startswith('* ') or lines[i].startswith('Status:')), len(lines))
    body = '\n'.join(lines[start:end]).strip()
    require(normalized(body) == normalized(TIMING_NOTE_BODY), 'Changed local timing NOTE body/method/times')
    return {'status': 'NOTE', 'adjudication': 'UNADJUDICATED', 'exit': 0,
            'typed_note': {'kind': 'exact_local_example_timing_note', 'method': 'compare_methods',
                           'user_seconds': 2.532, 'system_seconds': 0.992, 'elapsed_seconds': 18.406,
                           'raw_section': findings[0]+'\n'+body},
            'scope': 'Allowed only to classify warning objects. Timing/native release verdict remains UNADJUDICATED.'}

def check_counts(native, replay):
    require(native == replay, f'Native/replay count mismatch: {native} != {replay}')
    require(native['fail'] == 0, 'Native or replay failures are present')

def check_skips(actual, expected):
    require(collections.Counter(key(row, False) for row in actual) ==
            collections.Counter(key(row, False) for row in expected),
            'Unexpected or changed skip identity/scope')

def validate_warning(row, contract):
    name = category(row)
    require(row.get('file', '').startswith('test-') and row.get('test') and int(row['line']) > 0,
            'Missing warning test/file/line provenance')
    allowed = {key(item): item['category'] for item in contract['warnings']}
    require(key(row) in allowed, 'Changed warning message/test/file/line/caller needs fresh adjudication: ' + str(key(row)))
    require(name == allowed[key(row)], 'Changed warning category')
    return name

def selftest():
    contract = json.loads((HERE/'contract.json').read_text())
    example = dict(contract['warnings'][0])
    validate_warning(example, contract)
    for label, mutation in (
        ('unknown-warning', {'message': 'planted unknown warning'}),
        ('changed-caller', {'last_call': 'planted_unrelated_function()'}),
        ('changed-source', {'line': '999999'}),
    ):
        altered = example | mutation
        try:
            validate_warning(altered, contract)
        except ValueError:
            pass
        else:
            raise AssertionError(label + ' accepted')
    # A changed passing count is allowed when BOTH actual records agree.
    counts = {'fail': 0, 'warnings': 2, 'skips': 3, 'pass': 17}
    check_counts(counts, dict(counts))
    try:
        check_counts(counts, counts | {'pass': 18})
    except ValueError:
        pass
    else:
        raise AssertionError('count mismatch accepted')
    try:
        native_counts('unfinished')
    except ValueError:
        pass
    else:
        raise AssertionError('missing native summary accepted')
    check_skips(contract['skips'], contract['skips'])
    altered = [dict(row) for row in contract['skips']]
    altered[0]['message'] = 'planted missing-runtime skip'
    try:
        check_skips(altered, contract['skips'])
    except ValueError:
        pass
    else:
        raise AssertionError('unexpected skip accepted')
    planted = "* checking examples ... [13s/63s] NOTE\n"+TIMING_NOTE_BODY+"\n* checking tests ... OK\nStatus: 1 NOTE\n"
    native_verdict('Status: OK\n', 'exit=0\n', True)
    accepted = native_verdict(planted, 'exit=0\n', False, True)
    require(accepted['status']=='NOTE' and accepted['adjudication']=='UNADJUDICATED', 'Native NOTE lost')
    negative_notes = [
        ('no-opt-in', planted, 'exit=0\n', False, False),
        ('remote-note', planted, 'exit=0\n', True, True),
        ('changed-method', planted.replace('compare_methods','impute'), 'exit=0\n', False, True),
        ('changed-time', planted.replace('18.406','18.407'), 'exit=0\n', False, True),
        ('changed-section-time', planted.replace('[13s/63s]','[13s/64s]'), 'exit=0\n', False, True),
        ('changed-body', planted.replace('time > 5s','time > 6s'), 'exit=0\n', False, True),
        ('new-note', planted.replace('* checking tests ... OK','* checking tests ... NOTE\nnew unexpected note'), 'exit=0\n', False, True),
        ('changed-status', planted.replace('Status: 1 NOTE','Status: 2 NOTEs'), 'exit=0\n', False, True),
        ('warning', planted.replace('* checking tests ... OK','* checking tests ... WARNING'), 'exit=0\n', False, True),
        ('error', planted.replace('Status: 1 NOTE','Status: 1 ERROR, 1 NOTE'), 'exit=0\n', False, True),
        ('nonzero-exit', planted, 'exit=1\n', False, True),
        ('contradictory-exit', planted, 'exit=0\nexit=1\n', False, True),
        ('incomplete', planted.replace('Status: 1 NOTE',''), 'exit=0\n', False, True),
        ('contradictory-OK', planted.replace('Status: 1 NOTE','Status: OK'), 'exit=0\n', False, True),
    ]
    for label, log, code, remote, opt_in in negative_notes:
        try:
            native_verdict(log, code, remote, opt_in)
        except ValueError:
            pass
        else:
            raise AssertionError('Native control accepted: '+label)
    print('TIMING_NOTE_SELFTEST_PASS: exact local NOTE retained UNADJUDICATED; '+str(len(negative_notes))+' negative controls rejected')
    print('SELFTEST_PASS: known warning accepted; unknown warning, changed caller/source, missing summary, count mismatch and unexpected skip rejected; no old PASS-count pin')

def hash_record(paths):
    return {str(path): digest(path) for path in paths}

def norm_path(path):
    return re.sub('/+', '/', path.replace('\\', '/')).rstrip('/')

def audit_cell(name, folder, final, contract, output, remote, allow_local_timing_note=False):
    folder = folder.resolve()
    evidence = folder/'platform-evidence'
    check = folder/'frozen-check/pigauto.Rcheck'
    inputs = [evidence/'warning-replay.json', evidence/'warning-details.csv', evidence/'warning-objects.rds',
              evidence/'warning-session.txt', evidence/'check-command.txt', evidence/'check-environment.txt',
              evidence/'check-exit.txt', check/'00check.log', check/'tests/testthat.Rout']
    if remote:
        inputs += [evidence/f'identity-{phase}.json' for phase in ('download', 'pre', 'post')]
    else:
        inputs += [folder/'frozen-artifact/pigauto_0.11.0.tar.gz']
    missing = [str(path) for path in inputs if not path.is_file()]
    if missing:
        return {'cell': name, 'status': 'NOT_AVAILABLE', 'missing': missing}
    before = hash_record(inputs)
    record = {'cell': name, 'status': 'NOT_DONE', 'input_sha256': before, 'paths': list(before)}
    try:
        command = (evidence/'check-command.txt').read_text()
        environment = (evidence/'check-environment.txt').read_text().splitlines()
        require('NOT_CRAN=UNSET' in environment and '_R_CHECK_FORCE_SUGGESTS_=true' in environment,
                'Native CRAN environment differs')
        require(all(flag in command for flag in ('--as-cran', '--run-donttest', '--timings')) and '--no-manual' not in command,
                'Native check command differs')
        if remote:
            identities = [json.loads((evidence/f'identity-{phase}.json').read_text()) for phase in ('download', 'pre', 'post')]
            for phase, item in zip(('download', 'pre', 'post'), identities):
                require(item['sha256'] == SHA and item['size_bytes'] == SIZE and item['phase'] == phase,
                        'Wrong final archive identity')
                require(str(item['run_id']) == RUN and item['inventory_sha256'] == final['inventory_sha256'], 'Wrong run or inventory')
            for field in ('archive', 'workflow_sha', 'run_id', 'run_attempt'):
                require(len({str(item[field]) for item in identities}) == 1, 'Identity phase mismatch: '+field)
            require(identities[1]['archive'] in command, 'Native command did not name the identity archive')
            record['identity'] = identities[1]
        else:
            local_archive = folder/'frozen-artifact/pigauto_0.11.0.tar.gz'
            require(digest(local_archive) == SHA and local_archive.stat().st_size == SIZE, 'Wrong local final archive copy')
            require(str(local_archive) in command or final['path'] in command, 'Local command is not bound to the final archive')
            record['identity'] = {'sha256': SHA, 'size_bytes': SIZE, 'archive': str(local_archive)}
        check_log = (check/'00check.log').read_text()
        record['native_verdict'] = native_verdict(check_log, (evidence/'check-exit.txt').read_text(), remote, allow_local_timing_note)
        replay = json.loads((evidence/'warning-replay.json').read_text())
        require(replay.get('error') is None and replay['not_cran'] == 'UNSET' and replay['force_suggests'] == 'true',
                'Replay error or changed CRAN environment')
        require(replay['package'] == '0.11.0', 'Wrong replay package version')
        native_log_dir = re.search(r"using log directory [‘'](.+?)[’']", check_log)
        require(native_log_dir is not None, 'Native installation directory not recorded')
        require(norm_path(replay['package_path']).lower() == (norm_path(native_log_dir.group(1))+'/pigauto').lower(),
                'Replay was not bound to the native checked installation')
        native = native_counts((check/'tests/testthat.Rout').read_text())
        check_counts(native, replay['native_counts'])
        check_counts(native, replay['replay_counts'])
        objects, stderr = read_rds(evidence/'warning-objects.rds')
        record['rds_reader_stderr'] = stderr
        (output/(name+'-raw-objects.json')).write_text(json.dumps(objects, indent=2)+'\n')
        require(not objects['problems'], 'Replay raw problems are present')
        csv_rows = rows(evidence/'warning-details.csv')
        require(len(objects['warnings']) == len(csv_rows) == native['warnings'], 'CSV/RDS/native warning counts differ')
        require([int(row['id']) for row in csv_rows] == list(range(1, len(csv_rows)+1)), 'Warning IDs are not complete and ordered')
        require([key(row) for row in csv_rows] == [key(row) for row in objects['warnings']], 'CSV and raw warning identities differ')
        categorized = [dict(row, category=validate_warning(row, contract)) for row in csv_rows]
        require({row['category'] for row in categorized} == {row['category'] for row in contract['warnings']},
                'Changed warning-category set requires fresh adjudication')
        require(len(objects['skips']) == native['skips'], 'RDS/native skip counts differ')
        check_skips(objects['skips'], contract['skips'])
        with (output/(name+'-classified.csv')).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(categorized[0]), lineterminator='\n')
            writer.writeheader(); writer.writerows(categorized)
        require(hash_record(inputs) == before, 'Monitor/source inputs changed while being read; retry after completion')
        record.update(status=('PASS_WARNING_AUDIT_NATIVE_NOTE_UNADJUDICATED' if record['native_verdict']['status']=='NOTE' else 'PASS'), warning_audit_status='PASS', native_counts=native, warning_categories=dict(collections.Counter(row['category'] for row in categorized)),
                      warnings_identity_sha256=hashlib.sha256(json.dumps([key(row) for row in objects['warnings']]).encode()).hexdigest(),
                      diagnostic_elapsed_seconds=replay['diagnostic_elapsed_seconds'],
                      native_check_exit_and_timing=(evidence/'check-exit.txt').read_text().splitlines(),
                      skipped_scope='All raw skip message/test/file/line identities equal the previously reviewed scope; no missing-runtime skips added.')
    except (ValueError, KeyError, IndexError, TypeError, subprocess.SubprocessError) as error:
        record['error'] = str(error)
    return record

def audit(args):
    final = json.loads(Path(args.final_json).read_text())
    require(final['sha256'] == SHA and final['size_bytes'] == SIZE and digest(final['path']) == SHA, 'FINAL identity mismatch')
    contract = json.loads((HERE/'contract.json').read_text())
    output = HERE/args.output
    require(output.resolve().parent == HERE, 'Output must be one new child directory of this audit root')
    output.mkdir(exist_ok=False)
    cells = [('local-macos', Path(args.local_root), False)] + [(name, Path(args.direct_root)/('frozen-14ee7e5d-'+name), True)
             for name in ('windows-latest-release', 'windows-latest-devel', 'ubuntu-latest-release')]
    records = [audit_cell(name, folder, final, contract, output, remote, args.allow_local_compare_methods_timing_note) for name, folder, remote in cells]
    passed = sum(row.get('warning_audit_status')=='PASS' for row in records)
    native_note = any(row.get('native_verdict', {}).get('adjudication')=='UNADJUDICATED' for row in records)
    remote_ids = [row.get('identity') for row in records[1:] if row['status']=='PASS']
    if len(remote_ids)==3:
        require(len({str(row['workflow_sha']) for row in remote_ids}) == 1, 'Different remote executed checker commits')
    result = {'verdict': ('PASS warning identities and consequences only; native NOTE UNADJUDICATED' if native_note else 'PASS warning identities and consequences only') if passed==4 else 'NOT_DONE',
              'native_check_adjudication': 'UNADJUDICATED' if native_note else ('all_native_Status_OK' if passed==4 else 'NOT_DONE'),
              'exact_local_timing_note_opt_in': args.allow_local_compare_methods_timing_note,
              'archive_sha256': SHA, 'archive_size_bytes': SIZE, 'direct_run': RUN, 'passed_cells': passed,
              'records': records, 'category_consequences': CONSEQUENCES,
              'scope': 'Separate actual warning-object review. Not a final release-panel READY vote, coverage certification, timing approval or submission authorization.'}
    (output/'audit.json').write_text(json.dumps(result, indent=2)+'\n')
    lines = ['# Final-archive warning audit', '', result['verdict'], '', f'Archive `{SHA}`; direct run {RUN}; {passed}/4 warning audits passed.']
    for record in records:
        lines += ['', '## '+record['cell'], record['status'], record.get('error','')]
        if record.get('native_verdict'):
            lines += [json.dumps(record['native_verdict'])]
        if record.get('warning_audit_status')=='PASS':
            lines += [json.dumps(record['native_counts']), json.dumps(record['warning_categories'])]
    lines += ['', '## Consequences and limits']+[f'- {key}: {value}' for key,value in CONSEQUENCES.items()]+['', result['scope']]
    (output/'report.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'verdict': result['verdict'], 'passed_cells': passed, 'report': str(output/'report.md'),
                      'cell_status': {r['cell']:r['status'] for r in records}}, indent=2))
    return 0 if passed==4 else 2

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'selftest', 'audit'])
    parser.add_argument('--final-json', default=str(TASK/'FINAL.json'))
    parser.add_argument('--direct-root', default=str(TASK/'final-direct-39082452558'))
    parser.add_argument('--local-root', default='/private/tmp/pigauto-release-candidate/sha256-'+SHA+'/local-macos')
    parser.add_argument('--output', default='actual-results')
    parser.add_argument('--allow-local-compare-methods-timing-note', action='store_true',
                        help='Warning audit only: permit exact 2.532/0.992/18.406 compare_methods local NOTE; retain native UNADJUDICATED. Remote cells require Status OK.')
    args = parser.parse_args()
    if args.action=='prepare': prepare()
    elif args.action=='selftest': selftest()
    else: sys.exit(audit(args))
