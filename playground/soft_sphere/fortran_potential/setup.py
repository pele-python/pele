
from numpy.distutils.core import setup, Extension

setup(
      ext_modules=[Extension("_soft_sphere", ["_soft_sphere.f90"])]
      )
