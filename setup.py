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

f2py_options=[]
if False:
    #for bug testing
    f2py_options=["-DF2PY_REPORT_ON_ARRAY_COPY"]

fmodules = ModuleList(f2py_options=f2py_options)
fmodules.add_module("pygmin/mindist/overlap.f90")
fmodules.add_module("pygmin/optimize/mylbfgs_fort.f90")
fmodules.add_module("pygmin/optimize/mylbfgs_updatestep.f90")
fmodules.add_module("pygmin/potentials/fortran/AT.f")
fmodules.add_module("pygmin/potentials/fortran/ljpshiftfort.f")
fmodules.add_module("pygmin/potentials/fortran/lj.f90")
fmodules.add_module("pygmin/potentials/fortran/ljcut.f90")
fmodules.add_module("pygmin/potentials/fortran/rmdrvt.f90")
fmodules.add_module("pygmin/potentials/fortran/soft_sphere_pot.f90")
#fmodules.add_module("pygmin/potentials/fortran/lj_hess.f")
#fmodules.add_module("pygmin/potentials/rigid_bodies/rbutils.f90")
fmodules.add_module("pygmin/utils/_fortran_utils.f90")
fmodules.add_module("pygmin/transition_states/_orthogoptf.f90")
fmodules.add_module("pygmin/transition_states/_NEB_utils.f90")

cxx_modules = [ ]

fortran_modules = fmodules.module_list
hasboost = True
if hasboost:
    ext_modules = fortran_modules + cxx_modules
else:
    ext_modules = fortran_modules

setup(name='pygmin', 
      version='0.1', 
      description="Python implementation of GMIN, OPTIM, and PATHSAMPLE",
      url='http://github.com/js850/PyGMIN.git',
      packages=["pygmin",
                "pygmin.application",
                "pygmin.potentials",
                "pygmin.potentials.rigid_bodies",
                "pygmin.gui",
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
                "pygmin.accept_tests"
                ],
      ext_modules=ext_modules
        )

