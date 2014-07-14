from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import os
import numpy as np

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

setup(ext_modules=[
    Extension("mypotential", 
              ["_mypotential.h", "mypotential.c"],
              ),
    ]
    )
