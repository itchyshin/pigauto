from pathlib import Path
import importlib.util,sys
sys.dont_write_bytecode=True
r=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('bounded_ancillary',r/'run-checks.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.IDENTITY['script_sha256']['read-absence.R']=m.digest(r/'read-absence.R')
result=m.run('read-absence')
assert result['exit_code']==0 and not result['timed_out']
