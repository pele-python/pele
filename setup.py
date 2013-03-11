from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.misc_util import has_cxx_sources
import os

class ModuleList:
    def __init__(self, **kwargs):
        self.module_list = []
        self.kwargs = kwargs
    def add_module(self, filename):
        modname = filename.replace("/", ".")
        modname, ext = os.path.splitext(modname)
        self.module_list.append(Extension(modname, [filename], **self.kwargs))

extra_compile_args=[]
if False:
    #for bug testing
    extra_compile_args=["-DF2PY_REPORT_ON_ARRAY_COPY=1"]
#if True:
#    extra_compile_args += ["-ffixed-line-length-none"]

fmodules = ModuleList(extra_compile_args=extra_compile_args)
#fmodules.add_module("pygmin/mindist/overlap.f90")
fmodules.add_module("pygmin/mindist/minperm.f90")
fmodules.add_module("pygmin/optimize/mylbfgs_fort.f90")
fmodules.add_module("pygmin/optimize/mylbfgs_updatestep.f90")
fmodules.add_module("pygmin/potentials/fortran/AT.f90")
fmodules.add_module("pygmin/potentials/fortran/ljpshiftfort.f90")
fmodules.add_module("pygmin/potentials/fortran/lj.f90")
fmodules.add_module("pygmin/potentials/fortran/ljcut.f90")
fmodules.add_module("pygmin/potentials/fortran/rmdrvt.f90")
fmodules.add_module("pygmin/potentials/fortran/soft_sphere_pot.f90")
fmodules.add_module("pygmin/potentials/fortran/maxneib_lj.f90")
fmodules.add_module("pygmin/potentials/fortran/maxneib_blj.f90")
fmodules.add_module("pygmin/potentials/fortran/lj_hess.f90")
fmodules.add_module("pygmin/potentials/fortran/magnetic_colloids.f90")
#fmodules.add_module("pygmin/potentials/rigid_bodies/rbutils.f90")
fmodules.add_module("pygmin/utils/_fortran_utils.f90")
fmodules.add_module("pygmin/transition_states/_orthogoptf.f90")
fmodules.add_module("pygmin/transition_states/_NEB_utils.f90")
fmodules.add_module("pygmin/angleaxis/_aadist.f90")
fmodules.add_module("pygmin/accept_tests/_spherical_container.f90")

cxx_modules = [ ]

fortran_modules = fmodules.module_list
ext_modules = fortran_modules + cxx_modules

setup(name='pygmin', 
      version='0.1', 
      description="Python implementation of GMIN, OPTIM, and PATHSAMPLE",
      url='http://github.com/js850/PyGMIN.git',
      packages=["pygmin",
                "pygmin.application",
                "pygmin.potentials",
                "pygmin.potentials.rigid_bodies",
                "pygmin.gui",
                "pygmin.gui.ui",
                "pygmin.mindist",
                "pygmin.optimize",
                "pygmin.transition_states",
                "pygmin.transition_states.nebtesting",
                "pygmin.landscape",
                "pygmin.printing",
                "pygmin.takestep",
                "pygmin.utils",
                "pygmin.wham",
                "pygmin.storage",
                "pygmin.potentials.fortran",
                "pygmin.accept_tests",
                "pygmin.systems",
                "pygmin.angleaxis",
                ],
      ext_modules=ext_modules
        )


###########################################################
"""
programming note:  f2py and cython are not mutually compatible with distutils.
To use f2py, we must use numpy.distutils, but numpy.distutils doesn't support
cython.  There are several ways around this, but all of them are hacky.

1   Use numpy.distutils to compile the fortran extension modules and use
    distutils to compile the .pyx files.  This seems like an optimal solution,
    but it fails when you pass a flag to numpy.distutils that distutils doesn't
    understand, e.g. --fcompiler.  A possible solution is to check for the
    known incompatible flags and remove them when passing to distutils.  

2.  Compile cython .pyx files into .c files by hand and use numpy.distils
    to include the .c as source files.  This is the way scipy does it.
    We could also have setup.py accept a flag --cython which does this for
    all .pyx files.

Currently we are using method 1.
"""
###########################################################

import sys
#remove any flag incompatible with non-numpy distutils.
flag = "--fcompiler"
rmlist = []
for arg in sys.argv:
    if flag in arg:
        rmlist.append(arg)
for arg in rmlist:
    sys.argv.remove(arg)


#now build the cython modules
#we have to do this separately because cython isn't supported by numpy.distutils
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

cython_modules = [
                  Extension("pygmin.potentials.rigid_bodies._rbutils_cython",
                            ["pygmin/potentials/rigid_bodies/_rbutils_cython.pyx"])
                  ]

setup(
  ext_modules = cython_modules,
  cmdclass = {'build_ext': build_ext},
  include_dirs = [np.get_include()],
#  script_args = ['build_ext', '--inplace'],
)

