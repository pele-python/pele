"""Exercise the repaired wheel and reject libraries loaded from its build prefix."""
import ctypes
import importlib.metadata
import os
from pathlib import Path
import sys

import pele
from pele import version
from pele.optimize import cvode_opt, _lbfgs_cpp
from pele.potentials import _lj_cpp
from pele.mindist import minperm
import pytest

package = Path(pele.__file__).resolve().parent
source = Path(__file__).resolve().parents[1]
assert not package.is_relative_to(source), package
assert version.git_revision == "Unknown"
assert Path(pele.get_include()).is_dir()
for ui in (source / "pele/gui").rglob("*.ui"):
    assert (package / ui.relative_to(source / "pele").with_suffix(".py")).is_file(), ui
assert (package / "gui/ui/resources_rc.py").is_file()
assert (package / "potentials/lammps_potential_python.py").is_file()
assert not list((package / "potentials").glob("lammps_pele_cython*"))
files = importlib.metadata.files("pele")
for name in ("LICENSE", "COPYING", "LICENSE-THIRD-PARTY"):
    assert any(Path(str(file)).name == name for file in files), name
assert not os.environ.get("CONDA_PREFIX")
assert not os.environ.get("LD_LIBRARY_PATH")
assert not os.environ.get("DYLD_LIBRARY_PATH")

result = pytest.main(["--pyargs", "pele", "--ignore", str(package / "gui"), "-q"])
if sys.platform == "linux":
    loaded = Path("/proc/self/maps").read_text()
else:
    dyld = ctypes.CDLL(None)
    dyld._dyld_image_count.restype = ctypes.c_uint32
    dyld._dyld_get_image_name.argtypes = [ctypes.c_uint32]
    dyld._dyld_get_image_name.restype = ctypes.c_char_p
    loaded = "\n".join(dyld._dyld_get_image_name(i).decode() for i in range(dyld._dyld_image_count()))
for forbidden in ("/tmp/pele-wheel-deps", "/Cellar/", "/opt/homebrew/", "/usr/local/opt/"):
    assert forbidden not in loaded, f"Wheel loaded a build dependency from {forbidden}"
raise SystemExit(result)
