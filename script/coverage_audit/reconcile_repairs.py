#!/usr/bin/env python3
"""Conservative baseline accounting across unmerged repairs; never a pass gate."""
import argparse
import csv
import json
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from audit import BASE_SOURCE_SHA, coverage, sha, verify


def unchanged_line_map(old, new, diff):
    """Map only lines outside zero-context Git hunks, checking exact text."""
    mapping = {}
    cursor, offset = 1, 0
    for start, count, new_start, new_count in re.findall(
            r'^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@', diff, re.M):
        start, count = int(start), int(count or 1)
        new_count = int(new_count or 1)
        end = start if count == 0 else start - 1
        for line in range(cursor, end + 1):
            mapping[line] = line + offset
        cursor = start + count if count else start + 1
        offset += new_count - count
    for line in range(cursor, len(old) + 1):
        mapping[line] = line + offset
    for a, b in mapping.items():
        if b < 1 or b > len(new) or old[a - 1] != new[b - 1]:
            raise ValueError(f'Unchanged-block text mismatch at {a} -> {b}')
    return mapping


def git(root, *args):
    return subprocess.run(['git', '-C', str(root), *args], check=True,
                          capture_output=True, text=True).stdout


def validate_witness_receipt(expected_sha, csv_path, receipt_path):
    if receipt_path is None:
        raise ValueError('Named test CSV needs a separate provenance receipt')
    receipt = json.loads(receipt_path.read_text())
    if receipt.get('source_sha') != expected_sha or receipt.get('csv_sha256') != sha(csv_path):
        raise ValueError('Named test CSV provenance mismatch')


def selftest():
    assert unchanged_line_map(['a', 'b'], ['x', 'a', 'b'],
                             '@@ -0,0 +1 @@\n+x\n') == {1: 2, 2: 3}
    assert unchanged_line_map(['a', 'b'], ['a', 'x', 'b'],
                             '@@ -1,0 +2 @@\n+x\n') == {1: 1, 2: 3}
    assert unchanged_line_map(['a', 'b', 'c'], ['a', 'c'],
                             '@@ -2 +1,0 @@\n-b\n') == {1: 1, 3: 2}
    assert unchanged_line_map(['a', 'b', 'c'], ['a', 'x', 'c'],
                             '@@ -2 +2 @@\n-b\n+x\n') == {1: 1, 3: 3}
    try:
        unchanged_line_map(['a'], ['different'], '')
    except ValueError:
        pass
    else:
        raise AssertionError('Undeclared source change was accepted')
    with tempfile.TemporaryDirectory() as directory:
        csv_path = Path(directory) / 'witness.csv'
        receipt = Path(directory) / 'receipt.json'
        csv_path.write_text('filename,line\n')
        receipt.write_text(json.dumps(dict(source_sha='a' * 40, csv_sha256=sha(csv_path))))
        validate_witness_receipt('a' * 40, csv_path, receipt)
        for expected, path in [('b' * 40, receipt), ('a' * 40, None)]:
            try:
                validate_witness_receipt(expected, csv_path, path)
            except ValueError:
                pass
            else:
                raise AssertionError('Wrong witness provenance was accepted')
        csv_path.write_text('modified\n')
        try:
            validate_witness_receipt('a' * 40, csv_path, receipt)
        except ValueError:
            pass
        else:
            raise AssertionError('Modified witness CSV was accepted')
    print('REPAIR_RECONCILIATION_SELFTEST_PASS')


def reconcile(args):
    root = args.git_root.resolve()
    strict = verify(args.candidate, args.manifest, args.expected_sha,
                    None, root, None, False, True)
    current, counts = coverage(args.candidate)
    docs = Path(__file__).resolve().parents[2] / 'docs/dev-log/coverage-100'
    original, _ = coverage(docs / 'baseline.xml')
    provenance = json.loads((docs / 'baseline-provenance.json').read_text())
    targets = {(f, n): ('original', h) for f, lines in original.items() for n, h in lines.items()}
    for method in provenance['expected_torch_methods']:
        for line in method['tracked_lines']:
            targets.setdefault((method['file'], line), ('previously_untracked_method', 0))
    mappings = {}
    args.out.mkdir(parents=True, exist_ok=True)
    for file in sorted({f for f, _ in targets}):
        old = git(root, 'show', f'{BASE_SOURCE_SHA}:{file}').splitlines()
        new = git(root, 'show', f'{args.expected_sha}:{file}').splitlines()
        diff = git(root, 'diff', '--no-ext-diff', '--no-renames', '--unified=0',
                   BASE_SOURCE_SHA, args.expected_sha, '--', file)
        mappings[file] = unchanged_line_map(old, new, diff)
        if diff:
            (args.out / (Path(file).name + '.diff')).write_text(diff)
    witness = {}
    if args.named_tests:
        validate_witness_receipt(args.expected_sha, args.named_tests, args.named_tests_provenance)
        for row in csv.DictReader(args.named_tests.open()):
            key = (row['filename'], int(row['line']))
            names = set()
            for field in ('single_test_evidence', 'collective_test_evidence'):
                names.update(x for x in row.get(field, '').split(' | ')
                             if x and x != '<setup-or-unattributed-call>')
            witness.setdefault(key, set()).update(names)
    rows = []
    for (file, old_line), (kind, old_hit) in sorted(targets.items()):
        new_line = mappings[file].get(old_line)
        hit = current.get(file, {}).get(new_line) if new_line is not None else None
        names = sorted(witness.get((file, new_line), set()))
        if new_line is None:
            status = 'changed_or_removed_by_unmerged_repair'
        elif hit is None:
            status = 'missing_current_counter'
        elif hit <= 0:
            status = 'uncovered'
        elif not names:
            status = 'covered_without_named_test_witness'
        else:
            status = 'covered_with_named_test_set'
        rows.append(dict(file=file, baseline_line=old_line, current_line=new_line,
                         kind=kind, baseline_hits=old_hit, current_hits=hit,
                         disposition=status, named_test_set=' | '.join(names)))
    with (args.out / 'line-accounting.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    mapped_targets = {(r['file'], r['current_line']) for r in rows if r['current_line'] is not None}
    added = [dict(file=f, current_line=n, current_hits=h)
             for f, lines in current.items() for n, h in lines.items()
             if (f, n) not in mapped_targets]
    with (args.out / 'new-current-locations.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=('file', 'current_line', 'current_hits'))
        writer.writeheader(); writer.writerows(added)
    misses = [r for r in rows if r['kind'] == 'original' and r['baseline_hits'] == 0]
    result = dict(source_sha=args.expected_sha, baseline_source_sha=BASE_SOURCE_SHA,
                  candidate_xml_sha256=sha(args.candidate),
                  source_manifest_sha256=sha(args.manifest),
                  report_tool_sha256=sha(Path(__file__)), current_counts=counts,
                  current_uncovered=counts['valid'] - counts['covered'],
                  original_baseline_misses=len(misses), accounted_locations=len(rows),
                  new_current_locations=len(added),
                  new_current_uncovered=sum(r['current_hits'] <= 0 for r in added),
                  baseline_miss_dispositions=dict(Counter(r['disposition'] for r in misses)),
                  all_location_dispositions=dict(Counter(r['disposition'] for r in rows)),
                  changed_source_files=strict['source_changed_files'],
                  named_test_csv_sha256=sha(args.named_tests) if args.named_tests else None,
                  complete=False,
                  limitation='Changed or removed baseline locations remain unresolved pending separate repair approval. This report does not relax the strict verifier.')
    (args.out / 'reconciliation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('selftest')
    run = sub.add_parser('reconcile')
    for name in ('candidate', 'manifest', 'git-root', 'out'):
        run.add_argument('--' + name, required=True, type=Path)
    run.add_argument('--expected-sha', required=True)
    run.add_argument('--named-tests', type=Path)
    run.add_argument('--named-tests-provenance', type=Path)
    args = parser.parse_args()
    selftest() if args.command == 'selftest' else reconcile(args)


if __name__ == '__main__':
    main()
