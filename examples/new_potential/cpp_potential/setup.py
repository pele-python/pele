from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import os
import numpy as np

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

include_dirs = [numpy_include, "../../../source/", "../../../pele/potentials/"]


setup(ext_modules=[
    Extension("mypotential", 
              ["_mypotential.hpp", "mypotential.cpp"],
              include_dirs=include_dirs,
              language="c++",
              extra_compile_args=["-std=c++0x"],
              ),
    ]
    )
