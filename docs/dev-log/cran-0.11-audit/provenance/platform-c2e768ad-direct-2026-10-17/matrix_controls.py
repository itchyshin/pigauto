"""Exercise the matrix oracle against isolated copies of downloaded evidence."""
import contextlib
import io
from pathlib import Path
import shutil
import sys
import tempfile

repo = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(repo / "script/cran-0.11-platform"))
from verify_results import verify

source = Path(sys.argv[1])
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    for original in source.iterdir():
        shutil.copytree(original / "platform-evidence", root / original.name / "platform-evidence")
        target = root / original.name / "frozen-check/pigauto.Rcheck"
        target.mkdir(parents=True)
        shutil.copyfile(original / "frozen-check/pigauto.Rcheck/00check.log", target / "00check.log")
    with contextlib.redirect_stdout(io.StringIO()):
        verify(root)
    print("PASS positive matrix")
    cell = root / "frozen-c2e768ad-windows-latest-release"
    controls = [
        (cell / "platform-evidence/identity-post.json", lambda s: s.replace("c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2", "0" * 64), "wrong archive hash"),
        (cell / "frozen-check/pigauto.Rcheck/00check.log", lambda s: s.replace("Status: OK", "Status: 1 ERROR"), "failed check status"),
        (cell / "platform-evidence/check-environment.txt", lambda s: s.replace("NOT_CRAN=UNSET", "NOT_CRAN=true"), "wrong CRAN mode"),
        (cell / "platform-evidence/check-command.txt", lambda s: s.replace("--as-cran", "--no-manual"), "wrong check command"),
    ]
    for path, mutate, name in controls:
        original = path.read_text()
        path.write_text(mutate(original))
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                verify(root)
        except AssertionError:
            print("PASS negative control:", name)
        else:
            raise AssertionError("accepted " + name)
        finally:
            path.write_text(original)
print("MATRIX_CONTROLS_PASS positive and four negative controls; originals unchanged")
