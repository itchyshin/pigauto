import importlib.util
import json
from pathlib import Path

ROOT = Path('/private/tmp/pigauto-final-closure-20261017')
spec = importlib.util.spec_from_file_location('warning_audit', ROOT / 'final-warning-audit/audit.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
valid = {'error': None, 'not_cran': 'UNSET', 'force_suggests': 'true'}
for empty in (None, {}):
    module.validate_replay_environment(valid | {'error': empty})
negative = [valid | {'error': value} for value in ('planted failure', {'message': 'failure'}, [], False, 0, '')]
negative += [{k:v for k,v in valid.items() if k != 'error'}, valid | {'not_cran': 'true'}, valid | {'force_suggests': 'false'}, None]
for n, value in enumerate(negative):
    try:
        module.validate_replay_environment(value)
    except ValueError:
        pass
    else:
        raise AssertionError(f'Malformed/error/environment case {n} was accepted')
rows = json.loads((ROOT / 'replay-null-serialization-fix/actual-guard-reproduction.json').read_text())
for row in rows:
    module.validate_replay_environment(json.loads(Path(row['path']).read_text()))
print('REPLAY_ENVIRONMENT_CONTROLS_PASS_2_POSITIVE_10_NEGATIVE_4_ACTUAL')
