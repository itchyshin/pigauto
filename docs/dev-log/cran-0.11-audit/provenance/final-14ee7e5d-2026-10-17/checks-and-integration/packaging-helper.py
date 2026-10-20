#!/usr/bin/env python3
"""Package already-executed evidence only; no checks, network, or release changes."""
import gzip, hashlib, io, json, pathlib, platform, re, sys, tarfile, zlib

BASE = pathlib.Path('/private/tmp/pigauto-final-closure-20261017')
HERE = pathlib.Path(__file__).resolve().parent
FINAL = json.loads((BASE / 'FINAL.json').read_bytes())
SHA = FINAL['sha256']
assert SHA == '14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb'
OUT = HERE / 'publication'

def digest(b):
    return hashlib.sha256(b).hexdigest()

def checked_member(root, name, expected=None):
    rel = pathlib.PurePosixPath(name)
    if rel.is_absolute() or '..' in rel.parts or not rel.parts:
        raise ValueError('Unsafe evidence member')
    p = root.joinpath(*rel.parts)
    if not p.is_file() or p.is_symlink() or any(x.is_symlink() for x in p.parents if x != root.parent):
        raise ValueError('Missing/nonregular/symlink evidence member: ' + str(p))
    if root.resolve() not in p.resolve().parents:
        raise ValueError('Evidence member escaped root')
    b = p.read_bytes()
    if expected is not None and digest(b) != expected:
        raise ValueError('Evidence hash mismatch: ' + str(p))
    # Stop instead of publishing credential-looking material. No environment capture.
    if re.search(rb'(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9]{32,}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----)', b):
        raise ValueError('Credential-shaped material: ' + str(p))
    return p, b

def archive_identity():
    b = pathlib.Path(FINAL['path']).read_bytes()
    if digest(b) != SHA or len(b) != FINAL['size_bytes']:
        raise ValueError('Final release archive changed')
    return {'sha256': digest(b), 'size_bytes': len(b)}

def deterministic_gzip(b):
    out = io.BytesIO()
    with gzip.GzipFile(fileobj=out, mode='wb', filename='', mtime=0, compresslevel=9) as g:
        g.write(b)
    return out.getvalue()

def tar_gzip(members):
    out = io.BytesIO()
    with tarfile.open(fileobj=out, mode='w', format=tarfile.PAX_FORMAT) as t:
        for name, (_, b) in sorted(members.items()):
            info = tarfile.TarInfo(name)
            info.size, info.mode, info.mtime = len(b), 0o644, 0
            info.uid = info.gid = 0
            info.uname = info.gname = ''
            t.addfile(info, io.BytesIO(b))
    return deterministic_gzip(out.getvalue())

def read_manifest_bundle(root, manifest_name):
    p, original = checked_member(root, manifest_name)
    declared = json.loads(original)
    members = {name: checked_member(root, name, h) for name, h in declared.items()}
    members[manifest_name] = (p, original)
    return members

def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')

def main():
    before = archive_identity()
    service = BASE / 'grace-proof-14ee7e5d'
    vignette = BASE / 'vignette-timing-14ee7e5d'
    spelling = BASE / 'final-spelling-review'
    mi = pathlib.Path(FINAL['path']).parent / 'mi-matrix'
    assert json.loads((service / 'identity.json').read_bytes())['sha256'] == SHA
    assert SHA in json.dumps(json.loads((vignette / 'input.json').read_bytes())['artifact'])
    assert json.loads((spelling / 'receipt.json').read_bytes())['archive_sha256'] == SHA
    mi_manifest = json.loads((mi / 'manifest.json').read_bytes())
    assert mi_manifest['archive_sha256'] == SHA
    verification = json.loads((BASE / 'final-mi-verification.json').read_bytes())
    assert verification['artifact_sha256'] == SHA
    assert verification['artifact_source_commit'] == FINAL['source_commit']
    bundles = {
        'service-url-spelling-absence': read_manifest_bundle(service, 'proof-manifest.json'),
        'vignette-timing-doi': read_manifest_bundle(vignette, 'receipt-manifest.json'),
    }
    spell_members = {}
    for line in (spelling / 'files.sha256').read_text().splitlines():
        h, name = line.split('  ', 1)
        spell_members[name] = checked_member(spelling, name, h)
    spell_members['files.sha256'] = checked_member(spelling, 'files.sha256')
    actual_spell = {p.name for p in spelling.iterdir() if p.is_file()}
    assert set(spell_members) == actual_spell, 'Unaccounted spelling root files'
    bundles['spelling-review'] = spell_members
    mi_members = {'mi-matrix/manifest.json': checked_member(mi, 'manifest.json')}
    steps = [s for c in mi_manifest['cells'] for s in c['steps']]
    assert len(steps) == 12
    for s in steps:
        assert s['exit_code'] == 0 and s['archive_unchanged']
        assert s['archive_sha256_before'] == s['archive_sha256_after'] == SHA
        assert pathlib.PurePosixPath(s['log']).name == s['log']
        mi_members['mi-matrix/' + s['log']] = checked_member(mi, s['log'], s['log_sha256'])
    for name in ['final-mi-verification.json', 'verify-final-mi.py']:
        mi_members[name] = checked_member(BASE, name)
    bundles['final-mi'] = mi_members
    assert [len(bundles[k]) for k in bundles] == [32, 46, 6, 15]
    # Fail-closed negative controls never mutate source evidence.
    failures = []
    for name, expected in [('identity.json', '0' * 64), ('missing-negative-control.json', None)]:
        try:
            checked_member(service, name, expected)
        except ValueError:
            failures.append(name)
    assert len(failures) == 2
    if OUT.exists():
        raise ValueError('Refusing to overwrite publication output')
    OUT.mkdir()
    for directory in ['raw', 'views', 'manifests']:
        (OUT / directory).mkdir()
    bundle_receipts, view_receipts = [], []
    for label, members in bundles.items():
        raw = tar_gzip(members)
        assert raw == tar_gzip(members), 'Nondeterministic bundle'
        with tarfile.open(fileobj=io.BytesIO(raw), mode='r:gz') as t:
            assert t.getnames() == sorted(members)
            for m in t.getmembers():
                assert m.isfile() and t.extractfile(m).read() == members[m.name][1]
        target = OUT / 'raw' / (label + '.evidence.tar.gz')
        target.write_bytes(raw)
        member_receipts = []
        for name, (source, b) in sorted(members.items()):
            member_receipts.append({'member': name, 'source': str(source), 'sha256': digest(b), 'size_bytes': len(b)})
            # Selected small/summary JSONs, CSVs and READMEs; original manifests also readable.
            selected = name.endswith('.csv') or pathlib.PurePosixPath(name).name in {
                'README.md', 'adjudication.json', 'identity.json', 'absence-summary.json',
                'execution-summary.json', 'timing-summary.json', 'doi-metadata.json',
                'doi-adjudication.md', 'eval-false-proof.json', 'identity-after.json',
                'receipt.json', 'manifest.json', 'proof-manifest.json', 'receipt-manifest.json',
                'final-mi-verification.json'}
            if selected:
                view = OUT / 'views' / label / name
                view.parent.mkdir(parents=True, exist_ok=True)
                normalized = b.replace(b'\r\n', b'\n')
                view.write_bytes(normalized)
                r = {'bundle': label, 'member': name, 'source_sha256': digest(b),
                     'view': str(view.relative_to(OUT)), 'view_sha256': digest(normalized),
                     'normalization': 'CRLF_TO_LF_ONLY' if normalized != b else 'none'}
                if normalized != b:
                    raw_view = deterministic_gzip(b)
                    raw_path = view.with_name(view.name + '.raw.gz')
                    raw_path.write_bytes(raw_view)
                    assert gzip.decompress(raw_view) == b
                    r.update(raw_gzip=str(raw_path.relative_to(OUT)), raw_gzip_sha256=digest(raw_view))
                view_receipts.append(r)
        receipt = {'kind': 'evidence_bundle_not_release_archive', 'release_archive_sha256': SHA,
                   'bundle': str(target.relative_to(OUT)), 'bundle_sha256': digest(raw),
                   'bundle_size_bytes': len(raw), 'member_count': len(members), 'members': member_receipts}
        write_json(OUT / 'manifests' / (label + '.json'), receipt)
        bundle_receipts.append({k: v for k, v in receipt.items() if k != 'members'})
    # Reject any source mutation, including original manifests, during packaging.
    for members in bundles.values():
        for source, b in members.values():
            assert source.read_bytes() == b, 'Evidence changed during publication packaging'
    after = archive_identity()
    assert before == after
    write_json(OUT / 'views-manifest.json', view_receipts)
    result = {'status': 'PASS_BOUNDED_EVIDENCE_PACKAGING', 'release_readiness': 'NOT_ADJUDICATED',
              'release_archive_before': before, 'release_archive_after': after,
              'artifact_source_commit': FINAL['source_commit'], 'bundles': bundle_receipts,
              'original_member_count': sum(len(m) for m in bundles.values()),
              'view_count': len(view_receipts), 'normalized_view_count': sum(v['normalization'] != 'none' for v in view_receipts),
              'negative_controls_refused': failures, 'determinism': 'two independently generated byte-identical gzip streams per bundle; mtime/uid/gid zero; sorted members',
              'python': sys.version, 'platform': platform.platform(), 'zlib': zlib.ZLIB_VERSION,
              'helper_sha256': digest(pathlib.Path(__file__).read_bytes()),
              'source_manifest_count_correction': 'Service manifest declares 31 files, plus manifest=32 bundled; earlier expected 33 was inaccurate.',
              'excluded': ['dependency libraries and symlinks', 'exact-source clones', 'full check/installed trees', 'environment secrets', 'pending native/direct/platform warning evidence'],
              'limitations': ['Original NOTASSESSED scratch receipts preserved verbatim.', 'This packaging result does not adjudicate release readiness or URL cleanliness.', 'MI backend and checker files referenced by provenance are not copied; requested manifest and twelve logs plus independent verification are included.']}
    write_json(OUT / 'publication-result.json', result)
    (OUT / 'README.md').write_text('''# Final candidate evidence bundles

These four evidence archives preserve original proof bytes for final release archive SHA-256
14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb.
They are evidence bundles, not the R release tarball. No release READY claim is made.

`manifests/` records every original member hash and size, and each bundle hash.
`views-manifest.json` binds selected readable views to original members. Views apply
CRLF-to-LF only; changed views also have an original-byte `.raw.gz`. All original
bytes, original manifests, failed attempts and NOTASSESSED receipts remain in `raw/`.
No dependency library, symlink directory, exact-source clone or installed check tree is included.

`publication-result.json` records final archive checks before/after, deterministic
round-trip verification, negative controls, exclusions and limitations. Service proof
has 31 declared files plus its manifest; the earlier expected count of 33 was inaccurate.
The DOI access denial remains explicit. Pending native/direct/platform warning proofs
are outside this bundle. MI provenance references checker/backend sources; those
external source trees are not included in the requested manifest-and-logs bundle.

To rerun evidence-bundle packaging, copy the retained `packaging-helper.py` into
an otherwise fresh directory and run `python3 /absolute/path/to/fresh-directory/packaging-helper.py`.
The helper requires the original absolute scratch inputs under
`/private/tmp/pigauto-final-closure-20261017`, including `FINAL.json` and its named
proof directories, plus the final release archive and MI matrix at the absolute
paths recorded there. This published subtree alone is not a self-contained input set.
The helper creates `publication/` beside its copied script and refuses to overwrite
that output directory if it already exists.
Writing-style assessment: NOTASSESSED. Packaging is not a science or release review.
''')
    print(json.dumps({k: result[k] for k in ['status', 'original_member_count', 'view_count', 'normalized_view_count', 'release_archive_after']}))

if __name__ == '__main__':
    main()
