
from numpy.distutils.core import setup, Extension

setup(
      ext_modules=[Extension("_mypotential", ["_mypotential.f90"])]
      )
