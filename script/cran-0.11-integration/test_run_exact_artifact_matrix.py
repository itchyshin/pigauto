import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("run-exact-artifact-matrix.py")
SPEC = importlib.util.spec_from_file_location("exact_artifact_matrix", MODULE_PATH)
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


class CloneCellTests(unittest.TestCase):
    def test_clones_dependency_links_but_not_the_old_pigauto_install(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            template = root / "template"
            destination = root / "destination"
            template.mkdir()
            dependency = root / "dependency-package"
            dependency.mkdir()
            (template / "depA").symlink_to(dependency, target_is_directory=True)
            (template / "pigauto").mkdir()
            (template / "pigauto" / "DESCRIPTION").write_text("old candidate")

            RUNNER.clone_cell_dependencies(template, destination)

            self.assertTrue((destination / "depA").is_symlink())
            self.assertEqual((destination / "depA").resolve(), dependency.resolve())
            self.assertFalse((destination / "pigauto").exists())

    def test_refuses_unexpected_real_package_directories(self):
        with tempfile.TemporaryDirectory() as temporary:
            template = Path(temporary) / "template"
            destination = Path(temporary) / "destination"
            template.mkdir()
            (template / "pigauto").mkdir()
            (template / "untracked-package").mkdir()

            with self.assertRaisesRegex(ValueError, "must be a symlink"):
                RUNNER.clone_cell_dependencies(template, destination)


if __name__ == "__main__":
    unittest.main()
