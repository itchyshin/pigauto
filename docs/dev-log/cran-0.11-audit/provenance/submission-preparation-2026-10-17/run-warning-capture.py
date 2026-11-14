import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import tarfile

root = Path(__file__).resolve().parent
archive = Path('/private/tmp/pigauto-release-candidate/sha256-c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2/pigauto_0.11.0.tar.gz')
expected = 'c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2'
assert hashlib.sha256(archive.read_bytes()).hexdigest() == expected
assert archive.stat().st_size == 5136713
lib = root / 'library'
lib.mkdir(exist_ok=True)
check_root = root/'warning-replay.Rcheck'
assert not check_root.exists()
with tarfile.open(archive) as bundle:
    for member in bundle.getmembers():
        if member.isfile() and member.name.startswith('pigauto/tests/'):
            relative = Path(member.name).relative_to('pigauto')
            assert '..' not in relative.parts
            target = check_root/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(bundle.extractfile(member).read())
assert not (check_root/'DESCRIPTION').exists()
env = dict(os.environ)
env.pop('NOT_CRAN', None)
env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
           VECLIB_MAXIMUM_THREADS='1', _R_CHECK_FORCE_SUGGESTS_='true')
start = time.monotonic()
with (root/'warning-install.log').open('w') as log:
    subprocess.run(['R', 'CMD', 'INSTALL', '--install-tests', '-l', str(lib), str(archive)],
                   cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)
with (root/'warning-replay.log').open('w') as log:
    subprocess.run(['Rscript', '--vanilla', str(root/'capture-warnings.R'), str(root)],
                   cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=600)
assert hashlib.sha256(archive.read_bytes()).hexdigest() == expected
summary = json.loads((root/'warning-replay.json').read_text())
summary.update(archive_sha256=expected, total_elapsed_seconds=time.monotonic()-start)
(root/'warning-replay.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
print('EXACT_ARCHIVE_WARNING_CAPTURE_COMPLETE')
