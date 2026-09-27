from setuptools import setup, Extension
import os
import numpy as np
import subprocess
import sys

## Numpy header files
numpy_lib = os.path.split(np.__file__)[0]
numpy_include = os.path.join(numpy_lib, "core/include")

include_dirs = [numpy_include, "../../source/", "../../pele/potentials/"]

#
# run cython on the pyx files
#
# need to pass cython the include directory so it can find the .pyx files
cython_flags = (
    ["-I"]
    + [os.path.abspath("../../pele/potentials")]
    + ["-v"]
    + ["-X embedsignature=True"]
)


def generate_cython():
    cwd = os.path.abspath(os.path.dirname(__file__))
    print("Cythonizing sources")
    p = subprocess.call(
        [sys.executable, os.path.join(cwd, "../../cythonize.py"), "."]
        + cython_flags,
        cwd=cwd,
    )
    if p != 0:
        raise RuntimeError("Running cythonize failed!")


generate_cython()

setup(
    ext_modules=[
        Extension(
            "plate_potential",
            sources=["plate_potential.c"],
            include_dirs=include_dirs,
        )
    ]
)
