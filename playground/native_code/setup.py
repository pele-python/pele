import os
import numpy as np
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.misc_util import has_cxx_sources

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

include_sources = [
               "include/array.h",
               "include/potential.h",
               "include/simple_pairwise/potential.h",
               "include/_lbfgs.cpp",
               ]

extra_compile_args=['-Wextra','-pedantic','-funroll-loops','-O3', "-march=native", "-mtune=native", "-DNDEBUG"]

cxx_modules = [
            Extension("_lj", ["_lj.cpp"] + include_sources,
                      include_dirs=[numpy_include, "include/"],
                      extra_compile_args=extra_compile_args,
                      ),

            Extension("_pele", ["_pele.cpp"] + include_sources,
                      include_dirs=[numpy_include, "include/"],
                      extra_compile_args=extra_compile_args,
                      ),

            Extension("_lbfgs", ["_lbfgs.cpp"] + include_sources,
                      include_dirs=[numpy_include, "include/"],
                      extra_compile_args=extra_compile_args,
                      ),
            Extension("_pythonpotential", ["_pythonpotential.cpp"] + include_sources,
                      include_dirs=[numpy_include, "include/"],
                      extra_compile_args=extra_compile_args,
                      ),

               ]

setup(
      ext_modules=cxx_modules
        )
