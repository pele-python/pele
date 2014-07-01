import glob
import os
import sys
import subprocess

from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.misc_util import has_cxx_sources
import numpy as np

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 


def generate_cython():
    cwd = os.path.abspath(os.path.dirname(__file__))
    print("Cythonizing sources")
    p = subprocess.call([sys.executable,
                          os.path.join(cwd, 'cythonize.py'),
                          'pele'],
                         cwd=cwd)
    if p != 0:
        raise RuntimeError("Running cythonize failed!")

generate_cython()


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


#
# pure cython modules
#
extra_compile_args=['-Wall', '-Wextra','-pedantic','-funroll-loops','-O2',]

cxx_modules = [
            Extension("pele.optimize._cython_lbfgs", ["pele/optimize/_cython_lbfgs.c"],
                      include_dirs=[numpy_include],
                      extra_compile_args=extra_compile_args,
                      ),
            Extension("pele.potentials._cython_tools", ["pele/potentials/_cython_tools.c"],
                      include_dirs=[numpy_include],
                      extra_compile_args=extra_compile_args,
                      ),

               ]

fortran_modules = fmodules.module_list
ext_modules = fortran_modules + cxx_modules

setup(name='pele', 
      version='0.1', 
      description="Python implementation of GMIN, OPTIM, and PATHSAMPLE",
      url='https://github.com/pele-python/pele',
      packages=["pele",
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
                "pele.rates",
                # add the test directories
                "pele.potentials.tests",
                "pele.potentials.test_functions",
                "pele.mindist.tests",
                "pele.optimize.tests",
                "pele.transition_states.tests",
                "pele.landscape.tests",
                "pele.utils.tests",
                "pele.storage.tests",
                "pele.accept_tests.tests",
                "pele.systems.tests",
                "pele.thermodynamics.tests",
                "pele.rates.tests",
                ],
      ext_modules=ext_modules,
      # data files needed for the tests
      data_files=[('pele/potentials/tests', list(glob.glob('pele/potentials/tests/*.xyz'))),
                  ('pele/potentials/tests', list(glob.glob('pele/potentials/tests/*.xyzdr'))),
                  ('pele/transition_states/tests', list(glob.glob('pele/transition_states/tests/*.xyz'))),
                  ('pele/rates/tests', list(glob.glob('pele/rates/tests/*.data')) + list(glob.glob('pele/rates/tests/*.sqlite'))),
                  ('pele/mindist/tests', list(glob.glob('pele/mindist/tests/*.xyz')) + list(glob.glob('pele/mindist/tests/*.sqlite'))),
                  ('pele/storage/tests/', list(glob.glob('pele/storage/tests/*sqlite'))),
                 ]
        )


#
# build the c++ files
#

include_sources = ["source/pele" + f for f in os.listdir("source/pele") 
                   if f.endswith(".cpp")]
include_dirs = [numpy_include, "source"]

depends = [os.path.join("source/pele", f) for f in os.listdir("source/pele/") 
           if f.endswith(".cpp") or f.endswith(".h") or f.endswith(".hpp")]

# note: on my computer (ubuntu 12.04 gcc version 4.6.3), when compiled with the
# flag -march=native I run into problems.  Everything seems to run ok, but when
# I run it through valgrind, valgrind complains about an unrecognized
# instruction.  I don't have a clue what is causing this, but it's probably
# better to be on the safe side and not use -march=native
extra_compile_args = ['-std=c++0x',"-Wall", "-O3", '-funroll-loops']
# uncomment the next line to add extra optimization options
#extra_compile_args = ["-std=c++0x","-Wall", '-Wextra','-pedantic','-O3', "-march=native", "-mtune=native"]

# note: to compile with debug on and to override extra_compile_args use, e.g.
# OPT="-g -O2 -march=native" python setup.py ...

cxx_modules = [
    Extension("pele.potentials._lj_cpp", 
              ["pele/potentials/_lj_cpp.cxx"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
               
    Extension("pele.potentials._morse_cpp", 
              ["pele/potentials/_morse_cpp.cxx"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.potentials._hs_wca_cpp", 
              ["pele/potentials/_hs_wca_cpp.cxx"] + include_sources,
              include_dirs=include_dirs,
             extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
             ),
    Extension("pele.potentials._wca_cpp", 
              ["pele/potentials/_wca_cpp.cxx"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
             ),
    Extension("pele.potentials._harmonic_cpp", 
              ["pele/potentials/_harmonic_cpp.cxx"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
             ),
    Extension("pele.potentials._pele", 
              ["pele/potentials/_pele.cxx"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.optimize._pele_opt", 
              ["pele/optimize/_pele_opt.cxx"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    
    Extension("pele.optimize._lbfgs_cpp", 
              ["pele/optimize/_lbfgs_cpp.cxx", "source/lbfgs.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.optimize._modified_fire_cpp", 
              ["pele/optimize/_modified_fire_cpp.cxx", "source/modified_fire.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.potentials._pythonpotential", 
              ["pele/potentials/_pythonpotential.cxx"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
               ]

cxx_modules.append(
    Extension("pele.rates._ngt_cpp", 
              ["pele/rates/_ngt_cpp.cxx"] + ["sources/pele/graph.hpp", "sources/pele/ngt.hpp"],
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", 
              )
                   )

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
