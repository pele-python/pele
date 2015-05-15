import glob
import os
import sys
import subprocess
import shutil
import argparse
import shlex

import numpy as np
from distutils import sysconfig
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.command.build_ext import build_ext as old_build_ext

# Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 


# extract the -j flag and pass save it for running make on the CMake makefile
# extract -c flag to set compiler
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("-j", type=int, default=4)
parser.add_argument("-c", "--compiler", type=str, default=None)
jargs, remaining_args = parser.parse_known_args(sys.argv)

# record c compiler choice. use unix (gcc) by default  
# Add it back into remaining_args so distutils can see it also
idcompiler = None
if not jargs.compiler or jargs.compiler in ("unix", "gnu", "gcc"):
    idcompiler = "unix"
    remaining_args += ["-c", idcompiler]
elif jargs.compiler in ("intel", "icc", "icpc"):
    idcompiler = "intel"
    remaining_args += ["-c", idcompiler]

# set the remaining args back as sys.argv
sys.argv = remaining_args
print jargs, remaining_args
if jargs.j is None:
    cmake_parallel_args = []
else:
    cmake_parallel_args = ["-j" + str(jargs.j)]

#extra compiler args
cmake_compiler_extra_args=["-std=c++0x","-Wall", "-Wextra", "-pedantic", "-O3"]   

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
cython_flags=["-I"] + [os.path.abspath("pele/potentials")] + ["-v"]
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

cmake_build_dir = "build/cmake"

cxx_files = ["pele/potentials/_lj_cpp.cxx",
             "pele/potentials/_morse_cpp.cxx",
             "pele/potentials/_frozen_dof.cxx",
             "pele/potentials/_hs_wca_cpp.cxx",
             "pele/potentials/_wca_cpp.cxx",
             "pele/potentials/_harmonic_cpp.cxx",
             "pele/potentials/_inversepower_cpp.cxx",
             "pele/potentials/_pele.cxx",
             "pele/optimize/_pele_opt.cxx",
             "pele/optimize/_lbfgs_cpp.cxx",
             "pele/optimize/_modified_fire_cpp.cxx",
             "pele/potentials/_pythonpotential.cxx",
             "pele/angleaxis/_cpp_aa.cxx",
             "pele/utils/_cpp_utils.cxx",
             "pele/utils/_pressure_tensor.cxx",
             "pele/rates/_ngt_cpp.cxx",
             ]

def get_ldflags(opt="--ldflags"):
    """return the ldflags.  This was taken directly from python-config"""
    getvar = sysconfig.get_config_var
    pyver = sysconfig.get_config_var('VERSION')
    libs = getvar('LIBS').split() + getvar('SYSLIBS').split()
    libs.append('-lpython'+pyver)
    # add the prefix/lib/pythonX.Y/config dir, but only if there is no
    # shared library in prefix/lib/.
    if opt == '--ldflags':
        if not getvar('Py_ENABLE_SHARED'):
            libs.insert(0, '-L' + getvar('LIBPL'))
        if not getvar('PYTHONFRAMEWORK'):
            libs.extend(getvar('LINKFORSHARED').split())
    return ' '.join(libs)

# create file CMakeLists.txt from CMakeLists.txt.in 
with open("CMakeLists.txt.in", "r") as fin:
    cmake_txt = fin.read()
# We first tell cmake where the include directories are 
# note: the code to find python_includes was taken from the python-config executable
python_includes = [sysconfig.get_python_inc(), 
                   sysconfig.get_python_inc(plat_specific=True)]
cmake_txt = cmake_txt.replace("__PYTHON_INCLUDE__", " ".join(python_includes))
if isinstance(numpy_include, basestring):
    numpy_include = [numpy_include]
cmake_txt = cmake_txt.replace("__NUMPY_INCLUDE__", " ".join(numpy_include))
cmake_txt = cmake_txt.replace("__PYTHON_LDFLAGS__", get_ldflags())
cmake_txt = cmake_txt.replace("__COMPILER_EXTRA_ARGS__", '\"{}\"'.format(" ".join(cmake_compiler_extra_args)))
# Now we tell cmake which librarires to build 
with open("CMakeLists.txt", "w") as fout:
    fout.write(cmake_txt)
    fout.write("\n")
    for fname in cxx_files:
        fout.write("make_cython_lib(${CMAKE_CURRENT_SOURCE_DIR}/%s)\n" % fname)

def set_compiler_env(compiler_id):
    """
    set environment variables for the C and C++ compiler:
    set CC and CXX paths to `which` output because cmake
    does not alway choose the right compiler
    """
    env = os.environ.copy()
    if compiler_id.lower() in ("unix"):
        env["CC"] = subprocess.check_output(["which", "gcc"]).rstrip('\n')
        env["CXX"] = subprocess.check_output(["which", "g++"]).rstrip('\n')
    elif compiler_id.lower() in ("intel"):
        env["CC"] = subprocess.check_output(["which", "icc"]).rstrip('\n')
        env["CXX"] = subprocess.check_output(["which", "icpc"]).rstrip('\n')
    else:
        raise Exception("compiler_id not known")
    #this line only works is the build directory has been deleted
    cmake_compiler_args =shlex.split("-D CMAKE_C_COMPILER={} -D CMAKE_CXX_COMPILER={}".format(env["CC"],env["CXX"]))
    return env, cmake_compiler_args

def run_cmake(compiler_id="unix"):
    if not os.path.isdir(cmake_build_dir):
        os.makedirs(cmake_build_dir)
    print "\nrunning cmake in directory", cmake_build_dir
    cwd = os.path.abspath(os.path.dirname(__file__))
    env, cmake_compiler_args = set_compiler_env(compiler_id)
    
    p = subprocess.call(["cmake"] + cmake_compiler_args + [cwd], cwd=cmake_build_dir, env=env)
    if p != 0:
        raise Exception("running cmake failed")
    print "\nbuilding files in cmake directory"
    if len(cmake_parallel_args) > 0:
        print "make flags:", cmake_parallel_args
    p = subprocess.call(["make"] + cmake_parallel_args, cwd=cmake_build_dir)
    if p != 0:
        raise Exception("building libraries with CMake Makefile failed")
    print "finished building the extension modules with cmake\n"

run_cmake(compiler_id=idcompiler)
    

# Now that the cython libraries are built, we have to make sure they are copied to
# the correct location.  This means in the source tree if build in-place, or 
# somewhere in the build/ directory otherwise.  The standard distutils
# knows how to do this best.  We will overload the build_ext command class
# to simply copy the pre-compiled libraries into the right place
class build_ext_precompiled(old_build_ext):
    def build_extension(self, ext):
        """overload the function that build the extension
        
        This does nothing but copy the precompiled library stored in extension.sources[0]
        to the correct destination based on extension.name and whether it is an in-place build
        or not.
        """
        ext_path = self.get_ext_fullpath(ext.name)
        pre_compiled_library = ext.sources[0]
        if pre_compiled_library[-3:] != ".so":
            raise RuntimeError("library is not a .so file: " + pre_compiled_library)
        if not os.path.isfile(pre_compiled_library):
            raise RuntimeError("file does not exist: " + pre_compiled_library + " Did CMake not run correctly")
        print "copying", pre_compiled_library, "to", ext_path
        shutil.copy2(pre_compiled_library, ext_path)

# Construct extension modules for all the cxx files
# The `name` of the extension is, as usual, the python path (e.g. pele.optimize._lbfgs_cpp).
# The `source` of the extension is the location of the .so file
cxx_modules = []
for fname in cxx_files:
    name = fname.replace(".cxx", "")
    name = name.replace("/", ".")
    lname = os.path.basename(fname)
    lname = lname.replace(".cxx", ".so")
    pre_compiled_lib = os.path.join(cmake_build_dir, lname)
    cxx_modules.append(Extension(name, [pre_compiled_lib]))

setup(cmdclass=dict(build_ext=build_ext_precompiled),
      ext_modules=cxx_modules)

