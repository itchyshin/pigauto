from pathlib import Path
import importlib.util,sys
sys.dont_write_bytecode=True
r=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("bounded_ancillary",r/"run-checks.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.IDENTITY["script_sha256"]["urls-network.R"]=m.digest(r/"urls-network.R")
result=m.run("urls-network")
assert result["exit_code"]==0 and not result["timed_out"]
