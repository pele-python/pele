from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.misc_util import has_cxx_sources
import numpy as np
import os

os.environ["CXX"] = "g++-4.6"

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

#
# compile fortran extension modules
#

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
#fmodules.add_module("pele/mindist/overlap.f90")
fmodules.add_module("pele/mindist/minperm.f90")
#fmodules.add_module("pele/optimize/mylbfgs_fort.f90")
fmodules.add_module("pele/optimize/mylbfgs_updatestep.f90")
fmodules.add_module("pele/potentials/fortran/AT.f90")
fmodules.add_module("pele/potentials/fortran/ljpshiftfort.f90")
fmodules.add_module("pele/potentials/fortran/lj.f90")
fmodules.add_module("pele/potentials/fortran/ljcut.f90")
fmodules.add_module("pele/potentials/fortran/rmdrvt.f90")
fmodules.add_module("pele/potentials/fortran/soft_sphere_pot.f90")
#fmodules.add_module("pele/potentials/fortran/maxneib_lj.f90")
#fmodules.add_module("pele/potentials/fortran/maxneib_blj.f90")
fmodules.add_module("pele/potentials/fortran/lj_hess.f90")
fmodules.add_module("pele/potentials/fortran/morse.f90")
fmodules.add_module("pele/potentials/fortran/scdiff_periodic.f90")
#fmodules.add_module("pele/potentials/fortran/magnetic_colloids.f90")
#fmodules.add_module("pele/potentials/rigid_bodies/rbutils.f90")
fmodules.add_module("pele/utils/_fortran_utils.f90")
fmodules.add_module("pele/transition_states/_orthogoptf.f90")
fmodules.add_module("pele/transition_states/_NEB_utils.f90")
fmodules.add_module("pele/angleaxis/_aadist.f90")
fmodules.add_module("pele/accept_tests/_spherical_container.f90")

cxx_modules = [
            Extension("pele.optimize._cython_lbfgs", ["pele/optimize/_cython_lbfgs.c"],
                      include_dirs=[numpy_include],
                      extra_compile_args=['-Wall', '-Wextra','-pedantic','-funroll-loops','-O2',],
                      ),
            Extension("pele.potentials._cython_tools", ["pele/potentials/_cython_tools.c"],
                      include_dirs=[numpy_include],
                      extra_compile_args=['-Wall', '-Wextra','-pedantic','-funroll-loops','-O2',],
                      ),

               ]

fortran_modules = fmodules.module_list
ext_modules = fortran_modules + cxx_modules

setup(name='pele', 
      version='0.1', 
      description="Python implementation of GMIN, OPTIM, and PATHSAMPLE",
      url='https://github.com/pele-python/pele',
      packages=["pele",
                "pele.application",
                "pele.potentials",
                "pele.gui",
                "pele.gui.ui",
                "pele.mindist",
                "pele.optimize",
                "pele.transition_states",
                "pele.transition_states.nebtesting",
                "pele.landscape",
                "pele.takestep",
                "pele.utils",
                "pele.storage",
                "pele.potentials.fortran",
                "pele.accept_tests",
                "pele.systems",
                "pele.angleaxis",
                "pele.thermodynamics",
                ],
      ext_modules=ext_modules
        )

#
# build the c++ files
#

include_sources = ["source/pele" + f for f in os.listdir("source/pele") 
                   if f.endswith(".cpp")]
include_dirs = [numpy_include, "source"]

depends = ["source/pele" + f for f in os.listdir("source/pele/") 
           if f.endswith(".cpp") or f.endswith(".h")]

extra_compile_args = ["-std=c++0x","-Wall", "-Wextra", "-O3", "-march=native", "-mtune=native"]
# uncomment the next line to add extra optimization options
# extra_compile_args = ["-Wall", '-Wextra','-pedantic','-funroll-loops','-O3', "-march=native", "-mtune=native", "-DNDEBUG"]

# note: to compile with debug on and to override extra_compile_args use, e.g.
# OPT="-g -O2 -march=native" python setup.py ...

cxx_modules = [
    Extension("pele.potentials._lj_cpp", 
              ["pele/potentials/_lj_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
               
    Extension("pele.potentials._morse_cpp", 
              ["pele/potentials/_morse_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.potentials._hs_wca_cpp", 
              ["pele/potentials/_hs_wca_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
             extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
             ),
    Extension("pele.potentials._wca_cpp", 
              ["pele/potentials/_wca_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
             extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
             ),
    Extension("pele.potentials._harmonic_cpp", 
              ["pele/potentials/_harmonic_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
             extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
             ),
    Extension("pele.potentials._pele", 
              ["pele/potentials/_pele.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.optimize._pele_opt", 
              ["pele/optimize/_pele_opt.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.optimize._lbfgs_cpp", 
              ["pele/optimize/_lbfgs_cpp.cpp", "source/lbfgs.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.optimize._modified_fire_cpp", 
              ["pele/optimize/_modified_fire_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("playground.monte_carlo._pele_mc", 
              ["playground/monte_carlo/_pele_mc.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("playground.monte_carlo._monte_carlo_cpp", 
              ["playground/monte_carlo/_monte_carlo_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("playground.monte_carlo._takestep_cpp", 
              ["playground/monte_carlo/_takestep_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("playground.monte_carlo._accept_test_cpp", 
              ["playground/monte_carlo/_accept_test_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("playground.monte_carlo._conf_test_cpp", 
              ["playground/monte_carlo/_conf_test_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("playground.monte_carlo._action_cpp", 
              ["playground/monte_carlo/_action_cpp.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.potentials._pythonpotential", 
              ["pele/potentials/_pythonpotential.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
               ]

setup(ext_modules=cxx_modules,
      )


#
# we are currently using method 2 below
#
#have_cython = False
#if have_cython:
#
#    ###########################################################
#    """
#    programming note:  f2py and cython are not mutually compatible with distutils.
#    To use f2py, we must use numpy.distutils, but numpy.distutils doesn't support
#    cython.  There are several ways around this, but all of them are hacky.
#
#    1   Use numpy.distutils to compile the fortran extension modules and use
#        distutils to compile the .pyx files.  This seems like an optimal solution,
#        but it fails when you pass a flag to numpy.distutils that distutils doesn't
#        understand, e.g. --fcompiler.  A possible solution is to check for the
#        known incompatible flags and remove them when passing to distutils.  
#
#    2.  Compile cython .pyx files into .c files by hand and use numpy.distils
#        to include the .c as source files.  This is the way scipy does it.
#        We could also have setup.py accept a flag --cython which does this for
#        all .pyx files.
#
#    Currently we are using method 2.
#    """
#    ###########################################################
#
#    import sys
#    #remove any flag incompatible with non-numpy distutils.
#    flag = "--fcompiler"
#    rmlist = []
#    for arg in sys.argv:
#        if flag in arg:
#            rmlist.append(arg)
#    for arg in rmlist:
#        sys.argv.remove(arg)
#
#
#    #now build the cython modules
#    #we have to do this separately because cython isn't supported by numpy.distutils
#    from distutils.core import setup
#    from distutils.extension import Extension
#    from Cython.Distutils import build_ext
#    import numpy as np
#
#    cython_modules = [
#                      Extension("pele.potentials.rigid_bodies._rbutils_cython",
#                                ["pele/potentials/rigid_bodies/_rbutils_cython.pyx"])
#                      ]
#
#    setup(
#      ext_modules = cython_modules,
#      cmdclass = {'build_ext': build_ext},
#      include_dirs = [np.get_include()],
#    #  script_args = ['build_ext', '--inplace'],
#    )
