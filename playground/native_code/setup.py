import os
import numpy as np
from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build_src import build_src as np_build_src

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include')
include_dirs = [numpy_include, "include/"]

include_sources = [
#               "include/array.h",
#               "include/potential.h",
#               "include/simple_pairwise_potential.h",
               "include/_lbfgs.cpp",
               ]

extra_compile_args=['-Wextra','-pedantic','-funroll-loops','-O3', "-march=native", "-mtune=native", "-DNDEBUG"]

cxx_modules = [
            Extension("_lj", ["_lj.cpp"] + include_sources,
                      include_dirs=include_dirs,
                      extra_compile_args=extra_compile_args,
                      language="c++",
                      ),

            Extension("_pele", ["_pele.cpp"] + include_sources,
                      include_dirs=include_dirs,
                      extra_compile_args=extra_compile_args,
                      ),

            Extension("_lbfgs", ["_lbfgs.cpp"] + include_sources,
                      include_dirs=include_dirs,
                      extra_compile_args=extra_compile_args,
                      ),
            Extension("_pythonpotential", ["_pythonpotential.cpp"] + include_sources,
                      include_dirs=include_dirs,
                      extra_compile_args=extra_compile_args,
                      ),
               ]

setup(ext_modules=cxx_modules,
      )

cxx_f_modules = [
            Extension("_lj_cython", ["_lj_cython.cpp", "../../pele/potentials/fortran/lj.f90" ] + include_sources,
                      include_dirs=include_dirs,
                      extra_compile_args=extra_compile_args,
                      ),

                 ]


class build_src(np_build_src):
    """redefine the command class that makes the c sources for interacting with fortran objects"""
    def build_sources(self):
        """tell numpy.distutils to not build any sources"""
        pass

setup(
      ext_modules=cxx_f_modules,
      cmdclass=dict(build_src=build_src),
        )
