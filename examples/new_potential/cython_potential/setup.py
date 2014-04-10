from numpy.distutils.core import setup
from numpy.distutils.core import Extension


setup(
      ext_modules=[Extension("mypotential", ["mypotential.c"])]
      )