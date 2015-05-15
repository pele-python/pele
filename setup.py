import glob
import os
import sys
import subprocess

from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.misc_util import has_cxx_sources
import numpy as np

# Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

#
# Make the git revision visible.  Most of this is copied from scipy
# 
# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

def write_version_py(filename='pele/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM SCIPY SETUP.PY
git_revision = '%(git_revision)s'
"""
    GIT_REVISION = git_version()

    a = open(filename, 'w')
    try:
        a.write(cnt % dict(git_revision=GIT_REVISION))
    finally:
        a.close()
write_version_py()



#
# run cython on the pyx files
#
# need to pass cython the include directory so it can find the .pyx files
cython_flags=["-I"] + [os.path.abspath("pele/potentials")] + ["-v"] + ["-X embedsignature=True"]
def generate_cython():
    cwd = os.path.abspath(os.path.dirname(__file__))
    print("Cythonizing sources")
    p = subprocess.call([sys.executable,
                          os.path.join(cwd, 'cythonize.py'),
                          'pele'] + cython_flags,
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
#fmodules.add_module("pele/potentials/fortran/soft_sphere_pot.f90")
#fmodules.add_module("pele/potentials/fortran/maxneib_lj.f90")
#fmodules.add_module("pele/potentials/fortran/maxneib_blj.f90")
fmodules.add_module("pele/potentials/fortran/lj_hess.f90")
fmodules.add_module("pele/potentials/fortran/morse.f90")
fmodules.add_module("pele/potentials/fortran/scdiff_periodic.f90")
fmodules.add_module("pele/potentials/fortran/FinSin.f90")
fmodules.add_module("pele/potentials/fortran/gupta.f90")
#fmodules.add_module("pele/potentials/fortran/magnetic_colloids.f90")
#fmodules.add_module("pele/potentials/rigid_bodies/rbutils.f90")
fmodules.add_module("pele/utils/_fortran_utils.f90")
fmodules.add_module("pele/transition_states/_orthogoptf.f90")
fmodules.add_module("pele/transition_states/_NEB_utils.f90")
fmodules.add_module("pele/angleaxis/_aadist.f90")
fmodules.add_module("pele/accept_tests/_spherical_container.f90")


#
# compile the pure cython modules
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
                "pele.takestep.tests",
                "pele.utils.tests",
                "pele.storage.tests",
                "pele.accept_tests.tests",
                "pele.systems.tests",
                "pele.angleaxis.tests",
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
                  ('pele/utils/tests/', list(glob.glob('pele/utils/tests/*data'))),
                  ('pele/utils/tests/', list(glob.glob('pele/utils/tests/points*'))),
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
               
    Extension("pele.potentials._frozen_dof", 
              ["pele/potentials/_frozen_dof.cxx"] + include_sources,
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
    Extension("pele.potentials._inversepower_cpp", 
              ["pele/potentials/_inversepower_cpp.cxx"] + include_sources,
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
    Extension("pele.angleaxis._cpp_aa", 
              ["pele/angleaxis/_cpp_aa.cxx", "source/aatopology.cpp", "source/rotations.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.utils._cpp_utils", 
              ["pele/utils/_cpp_utils.cxx", "source/rotations.cpp"] + include_sources,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends,
              ),
    Extension("pele.utils._pressure_tensor", 
              ["pele/utils/_pressure_tensor.cxx", "source/pressure_tensor.cpp"] + include_sources,
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

