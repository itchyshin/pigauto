from pathlib import Path
p=Path(__file__).resolve().parent/'audit.py';s=p.read_text()
anchor='def check_counts(native, replay):'
code='''TIMING_NOTE_BODY = """Examples with CPU (user + system) or elapsed time > 5s
                 user system elapsed
compare_methods 2.532  0.992  18.406"""

def native_verdict(check_log, exit_text, remote, allow_local_timing_note=False):
    """Admit one exact local diagnostic case, retaining an unadjudicated native NOTE."""
    lines = check_log.splitlines()
    require('exit=0' in exit_text.splitlines(), 'Native check exit was not zero')
    statuses = [line for line in lines if line.startswith('Status:')]
    findings = [line for line in lines if re.match(r'^\\* .*\\b(?:NOTE|WARNING|ERROR)$', line)]
    if statuses == ['Status: OK']:
        require(not findings, 'Status OK contradicts native findings')
        return {'status': 'OK', 'adjudication': 'native_Status_OK', 'exit': 0}
    require(not remote and allow_local_timing_note, 'Native check is not Status OK; local exact timing-note opt-in absent or remote cell')
    require(statuses == ['Status: 1 NOTE'], 'Native status is not exactly one NOTE with no errors/warnings')
    require(len(findings) == 1 and re.fullmatch(r'\\* checking examples \\.\\.\\.(?: \\[[^\\]\\r\\n]+\\])? NOTE', findings[0]),
            'Unexpected native NOTE/method or error/warning finding')
    start = lines.index(findings[0]) + 1
    end = next((i for i in range(start, len(lines)) if lines[i].startswith('* ') or lines[i].startswith('Status:')), len(lines))
    body = '\\n'.join(lines[start:end]).strip()
    require(normalized(body) == normalized(TIMING_NOTE_BODY), 'Changed local timing NOTE body/method/times')
    return {'status': 'NOTE', 'adjudication': 'UNADJUDICATED', 'exit': 0,
            'typed_note': {'kind': 'exact_local_example_timing_note', 'method': 'compare_methods',
                           'user_seconds': 2.532, 'system_seconds': 0.992, 'elapsed_seconds': 18.406,
                           'raw_section': findings[0]+'\\n'+body},
            'scope': 'Allowed only to classify warning objects. Timing/native release verdict remains UNADJUDICATED.'}

'''
assert anchor in s;s=s.replace(anchor,code+anchor)
s=s.replace('def audit_cell(name, folder, final, contract, output, remote):','def audit_cell(name, folder, final, contract, output, remote, allow_local_timing_note=False):')
s=s.replace("        require('Status: OK' in check_log.splitlines(), 'Native check is not Status OK')\n        require('exit=0' in (evidence/'check-exit.txt').read_text().splitlines(), 'Native check exit was not zero')", "        record['native_verdict'] = native_verdict(check_log, (evidence/'check-exit.txt').read_text(), remote, allow_local_timing_note)")
s=s.replace("        record.update(status='PASS', native_counts=native,", "        record.update(status=('PASS_WARNING_AUDIT_NATIVE_NOTE_UNADJUDICATED' if record['native_verdict']['status']=='NOTE' else 'PASS'), warning_audit_status='PASS', native_counts=native,")
s=s.replace('audit_cell(name, folder, final, contract, output, remote) for name, folder, remote in cells','audit_cell(name, folder, final, contract, output, remote, args.allow_local_compare_methods_timing_note) for name, folder, remote in cells')
s=s.replace("passed = sum(row['status']=='PASS' for row in records)", "passed = sum(row.get('warning_audit_status')=='PASS' for row in records)\n    native_note = any(row.get('native_verdict', {}).get('adjudication')=='UNADJUDICATED' for row in records)")
s=s.replace("'verdict': 'PASS warning identities and consequences only' if passed==4 else 'NOT_DONE',", "'verdict': ('PASS warning identities and consequences only; native NOTE UNADJUDICATED' if native_note else 'PASS warning identities and consequences only') if passed==4 else 'NOT_DONE',\n              'native_check_adjudication': 'UNADJUDICATED' if native_note else ('all_native_Status_OK' if passed==4 else 'NOT_DONE'),\n              'exact_local_timing_note_opt_in': args.allow_local_compare_methods_timing_note,")
s=s.replace("if record['status']=='PASS':\n            lines +=", "if record.get('native_verdict'):\n            lines += [json.dumps(record['native_verdict'])]\n        if record.get('warning_audit_status')=='PASS':\n            lines +=")
s=s.replace("    parser.add_argument('--output', default='actual-results')", "    parser.add_argument('--output', default='actual-results')\n    parser.add_argument('--allow-local-compare-methods-timing-note', action='store_true',\n                        help='Warning audit only: permit exact 2.532/0.992/18.406 compare_methods local NOTE; retain native UNADJUDICATED. Remote cells require Status OK.')")
needle="    print('SELFTEST_PASS:"
i=s.index(needle)
s=s[:i]+'''    planted = "* checking examples ... [13s/63s] NOTE\\n"+TIMING_NOTE_BODY+"\\n* checking tests ... OK\\nStatus: 1 NOTE\\n"
    native_verdict('Status: OK\\n', 'exit=0\\n', True)
    accepted = native_verdict(planted, 'exit=0\\n', False, True)
    require(accepted['status']=='NOTE' and accepted['adjudication']=='UNADJUDICATED', 'Native NOTE lost')
    negative_notes = [
        ('no-opt-in', planted, 'exit=0\\n', False, False),
        ('remote-note', planted, 'exit=0\\n', True, True),
        ('changed-method', planted.replace('compare_methods','impute'), 'exit=0\\n', False, True),
        ('changed-time', planted.replace('18.406','18.407'), 'exit=0\\n', False, True),
        ('changed-body', planted.replace('time > 5s','time > 6s'), 'exit=0\\n', False, True),
        ('new-note', planted.replace('* checking tests ... OK','* checking tests ... NOTE\\nnew unexpected note'), 'exit=0\\n', False, True),
        ('changed-status', planted.replace('Status: 1 NOTE','Status: 2 NOTEs'), 'exit=0\\n', False, True),
        ('warning', planted.replace('* checking tests ... OK','* checking tests ... WARNING'), 'exit=0\\n', False, True),
        ('error', planted.replace('Status: 1 NOTE','Status: 1 ERROR, 1 NOTE'), 'exit=0\\n', False, True),
        ('nonzero-exit', planted, 'exit=1\\n', False, True),
        ('incomplete', planted.replace('Status: 1 NOTE',''), 'exit=0\\n', False, True),
        ('contradictory-OK', planted.replace('Status: 1 NOTE','Status: OK'), 'exit=0\\n', False, True),
    ]
    for label, log, code, remote, opt_in in negative_notes:
        try:
            native_verdict(log, code, remote, opt_in)
        except ValueError:
            pass
        else:
            raise AssertionError('Native control accepted: '+label)
    print('TIMING_NOTE_SELFTEST_PASS: exact local NOTE retained UNADJUDICATED; '+str(len(negative_notes))+' negative controls rejected')
''' +s[i:]
p.write_text(s)
