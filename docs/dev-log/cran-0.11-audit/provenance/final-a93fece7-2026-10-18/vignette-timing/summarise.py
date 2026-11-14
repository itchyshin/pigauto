from pathlib import Path
import csv, hashlib, json

root = Path(__file__).resolve().parent
meta = json.loads((root / 'input.json').read_text())
runs = json.loads((root / 'execution.json').read_text())
after = json.loads((root / 'identity-after.json').read_text())
assert after['all_six_success'] and after['archive_unchanged'] and after['source_rmds_unchanged']
assert len(runs) == 6 and all(r['exit'] == 0 for r in runs)
expected_start_false = {
    'canonical-journey', 'declarations', 'install', 'quickstart', 'train', 'print-fit',
    'predict', 'evaluate', 'plot-history', 'plot-uncertainty', 'mi-workflow',
    'active-impute', 'advanced-controls', 'exact-route', 'phylo-signal-gate',
    'phylo-signal-gate-off', 'covariate-helpers', 'gpu', 'cache', 'own-data',
}
rows = []
for run in runs:
    name = run['file']; stem = name.removesuffix('.Rmd')
    result_path = root / 'rendered' / stem / 'render-result.json'
    trace_path = root / 'rendered' / stem / 'chunk-options.jsonl'
    output_path = root / 'rendered' / stem / (stem + '.html')
    result = json.loads(result_path.read_text())
    chunks = [json.loads(line) for line in trace_path.read_text().splitlines() if line.strip()]
    assert output_path.is_file() and run['exit'] == 0
    assert result['error'] in ({}, None) and result['output']
    assert run['source_rmd_sha256'] == meta['source_rmd_sha256'][name]
    assert hashlib.sha256((root / 'source' / 'pigauto' / 'vignettes' / name).read_bytes()).hexdigest() == run['source_rmd_sha256']
    if name == 'getting-started.Rmd':
        chunk_by_label = {c['label']: c['eval'] for c in chunks}
        assert all(chunk_by_label.get(label) is False for label in expected_start_false)
    if name == 'multiple-imputation.Rmd':
        assert all(c['eval'] is False for c in chunks if c['label'] != 'setup')
    rows.append({
        'file': name,
        'render_elapsed_seconds': result['elapsed_seconds'],
        'process_elapsed_seconds': run['process_elapsed_seconds'],
        'exit': run['exit'],
        'start_utc': run['start_utc'],
        'render_start_utc': result['start_utc'],
        'source_rmd_sha256': run['source_rmd_sha256'],
        'render_script_sha256': run['render_script_sha256'],
        'enabled_chunks': sum(c.get('eval') is True for c in chunks),
        'disabled_chunks': sum(c.get('eval') is False for c in chunks),
        'chunk_trace': str(trace_path.relative_to(root)),
        'rendered_html': str(output_path.relative_to(root)),
        'log': str(Path(run['log']).relative_to(root)),
        'log_sha256': run['log_sha256'],
        'command': run['command'],
    })
summary = {
    'reviewer': 'pat', 'verdict': 'PASS',
    'scope': 'Bounded six-vignette diagnostic against the frozen successor archive',
    'separate_diagnostic': True,
    'artifact': meta['artifact'],
    'all_six_success': True,
    'rows': rows,
    'archive_sha256_after': after['archive_sha256_after'],
    'archive_unchanged': after['archive_unchanged'],
    'source_rmds_unchanged': after['source_rmds_unchanged'],
    'source_rmd_sha256_before': after['source_rmd_sha256_before'],
    'source_rmd_sha256_after': after['source_rmd_sha256_after'],
    'library': meta['library'],
    'installed_library_evidence': meta['installed_library_evidence'],
    'runtime': {'R': '4.6.0 (2026-04-24)', 'platform': 'aarch64-apple-darwin23',
                'torch_installed': True, 'pandoc': '3.12', 'NOT_CRAN': 'UNSET',
                'thread_caps': 1, 'per_vignette_timeout_seconds': meta['per_vignette_timeout_seconds'],
                'total_hard_cap_seconds': meta['total_timeout_seconds']},
    'estimate_minutes': '5-8',
    'actual_total_process_seconds': after['total_elapsed_seconds'],
    'limits': ['Not a native R CMD check or Windows check.',
               'Only authored enabled chunks ran; teaching examples marked eval=FALSE remained disabled.',
               'No global eval override, no retries, no warning replay, and no backend packages installed.'],
}
(root / 'timing-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
fields = ['file','render_elapsed_seconds','process_elapsed_seconds','exit','start_utc']
with (root / 'timing-summary.csv').open('w', newline='') as out:
    writer = csv.DictWriter(out, fieldnames=fields); writer.writeheader()
    writer.writerows({k: row[k] for k in fields} for row in rows)
print(json.dumps({'verdict': summary['verdict'], 'rows': rows,
                  'total_process_seconds': summary['actual_total_process_seconds']}, indent=2))
