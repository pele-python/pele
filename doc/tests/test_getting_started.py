from contextlib import chdir, redirect_stdout
from io import StringIO
from pathlib import Path
from runpy import run_path
from tempfile import TemporaryDirectory
import unittest

from pele.storage import Database


class GettingStartedTests(unittest.TestCase):
    def test_manual_basinhopping_stores_minima(self):
        with TemporaryDirectory() as directory, chdir(directory), redirect_stdout(StringIO()):
            example = Path(__file__).resolve().parents[2] / "examples/getting_started/example2_global_optimization.py"
            bh_no_system_class = run_path(str(example))["bh_no_system_class"]
            bh_no_system_class()
            database = Database("lj17.sqlite", createdb=False)
            self.assertGreater(database.number_of_minima(), 0)


if __name__ == "__main__":
    unittest.main()
