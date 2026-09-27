"""Build pele. Works with `pip install .` and `python setup.py build_ext -i`.

The Cython/Fortran/C++ extensions are compiled by CMake (from CMakeLists.txt.in)
inside build_ext; setuptools then just copies the resulting libraries.

Options (command-line flags for direct `setup.py` use, env vars for pip):
  -j N / PELE_JOBS               parallel make jobs (default: all cores)
  -c COMPILER                    unix (default) or intel
  --build-type / PELE_BUILD_TYPE Release (default), Debug, RelWithDebInfo, MemCheck, Greene
  --with-cvode / PELE_WITH_CVODE 1 (default) builds the CVODE/BDF optimizers (needs SUNDIALS)
  --native / PELE_NATIVE         1 (default) adds -march=native; set 0 for portable binaries
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
import sysconfig

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as old_build_ext
from setuptools.command.build_py import build_py

encoding = "utf-8"
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("-j", type=int, default=int(os.environ.get("PELE_JOBS", os.cpu_count() or 4)))
parser.add_argument("-c", "--compiler", type=str, default=None)
parser.add_argument("--opt-report", action="store_true", default=False,
                    help="Print optimization report (for Intel compiler)")
parser.add_argument("--build-type", type=str, default=os.environ.get("PELE_BUILD_TYPE", "Release"),
                    help="Release, Debug, RelWithDebInfo, MemCheck, Greene")
parser.add_argument("--with-cvode", type=int, default=int(os.environ.get("PELE_WITH_CVODE", 1)),
                    help="Build with CVODE. Needed for attractor identification.")
parser.add_argument("--native", type=int, default=int(os.environ.get("PELE_NATIVE", 1)),
                    help="Compile with -march=native (non-portable binaries)")
jargs, remaining_args = parser.parse_known_args(sys.argv)

if not jargs.compiler or jargs.compiler in ("unix", "gnu", "gcc"):
    idcompiler = "unix"
elif jargs.compiler in ("intelem", "intel", "icc", "icpc"):
    idcompiler = "intel"
else:
    raise ValueError("unknown compiler " + jargs.compiler)
# Only add the option back if it was really set (setup.py install does not allow -c)
if jargs.compiler:
    remaining_args += ["-c", idcompiler]
sys.argv = remaining_args

build_type = jargs.build_type
with_cvode = jargs.with_cvode

common_args = ["-std=c++2a", "-Wall", "-Wextra", "-pedantic", "-fPIC", "-D_GLIBCXX_USE_CXX11_ABI=1"]
build_type_args = {
    "Release": ["-O3", "-DNDEBUG"],
    "Greene": ["-O3", "-DNDEBUG", "-unroll", "-ip", "-axCORE-AVX512", "-qopenmp",
               "-qopt-report-stdout", "-qopt-report-phase=openmp"],
    "Debug": ["-ggdb3", "-O0"],
    "RelWithDebInfo": ["-g", "-O3"],
    "MemCheck": ["-g", "-O0", "-fsanitize=address", "-fsanitize=leak"],
}
if build_type not in build_type_args:
    raise ValueError("Unknown build type: " + build_type)
# env CXXFLAGS first (conda-forge hardening/arch flags), ours after so they win
cmake_compiler_extra_args = (
    shlex.split(os.environ.get("CXXFLAGS", "")) + common_args + build_type_args[build_type]
)
if jargs.native and build_type in ("Release", "RelWithDebInfo"):
    cmake_compiler_extra_args += ["-march=native"]
if idcompiler == "unix":
    cmake_compiler_extra_args += ["-fopenmp"]
else:
    cmake_compiler_extra_args += ["-qopenmp"]
    if jargs.opt_report:
        cmake_compiler_extra_args += ["-qopt-report=5"]

cmake_build_dir = "build/cmake"

fortran_files = [
    "pele/mindist/minperm.f90",
    "pele/optimize/mylbfgs_updatestep.f90",
    "pele/potentials/fortran/AT.f90",
    "pele/potentials/fortran/ljpshiftfort.f90",
    "pele/potentials/fortran/lj.f90",
    "pele/potentials/fortran/ljcut.f90",
    "pele/potentials/fortran/lj_hess.f90",
    "pele/potentials/fortran/morse.f90",
    "pele/potentials/fortran/scdiff_periodic.f90",
    "pele/potentials/fortran/FinSin.f90",
    "pele/potentials/fortran/gupta.f90",
    "pele/utils/_fortran_utils.f90",
    "pele/transition_states/_orthogoptf.f90",
    "pele/transition_states/_NEB_utils.f90",
    "pele/angleaxis/_aadist.f90",
    "pele/accept_tests/_spherical_container.f90",
]

c_files = [
    "pele/optimize/_cython_lbfgs.c",
    "pele/potentials/_cython_tools.c",
]

cxx_files = [
    "pele/angleaxis/_cpp_aa.cxx",
    "pele/distance/_get_distance_cpp.cxx",
    "pele/distance/_put_in_box_cpp.cxx",
    "pele/optimize/_gradient_descent_cpp.cxx",
    "pele/optimize/_cosine_gradient_descent_cpp.cxx",
    "pele/optimize/_lbfgs_cpp.cxx",
    "pele/optimize/_modified_fire_cpp.cxx",
    "pele/optimize/_pele_opt.cxx",
    "pele/potentials/_frenkel.cxx",
    "pele/potentials/_frozen_dof.cxx",
    "pele/potentials/_harmonic_cpp.cxx",
    "pele/potentials/_hs_wca_cpp.cxx",
    "pele/potentials/_inversepower_cpp.cxx",
    "pele/potentials/_inversepower_hs_cpp.cxx",
    "pele/potentials/_inversepower_stillinger_cpp.cxx",
    "pele/potentials/_inversepower_stillinger_cut_cpp.cxx",
    "pele/potentials/_inversepower_stillinger_cut_quad.cxx",
    "pele/potentials/_lj_cpp.cxx",
    "pele/potentials/_morse_cpp.cxx",
    "pele/potentials/_pele.cxx",
    "pele/potentials/_pspin_spherical_cpp.cxx",
    "pele/potentials/_pythonpotential.cxx",
    "pele/potentials/_radial_gaussian_cpp.cxx",
    "pele/potentials/_sumgaussianpot_cpp.cxx",
    "pele/potentials/_wca_cpp.cxx",
    "pele/potentials/atlj.cxx",
    "pele/potentials/combine_potentials.cxx",
    "pele/potentials/cpp_test_functions.cxx",
    "pele/rates/_ngt_cpp.cxx",
    "pele/utils/_cpp_utils.cxx",
    "pele/utils/_pressure_tensor.cxx",
]
if with_cvode:
    cxx_files += [
        "pele/optimize/_mxd_end_only.cxx",
        "pele/optimize/cvode_opt.cxx",
        "pele/optimize/extended_mixed_descent.cxx",
        "pele/optimize/generic_mixed_descent.cxx",
    ]


def git_version():
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                             env={"PATH": os.environ.get("PATH", ""), "LC_ALL": "C"})
        return out.stdout.strip().decode("ascii") or "Unknown"
    except OSError:
        return "Unknown"


def write_version_py(filename="pele/version.py"):
    with open(filename, "w") as f:
        # layout (git_revision on the 3rd line) is parsed by downstream scripts
        f.write("\n# THIS FILE IS GENERATED FROM SCIPY SETUP.PY\ngit_revision = '%s'\n" % git_version())


def generate_cython():
    cwd = os.path.abspath(os.path.dirname(__file__))
    print("Cythonizing sources")
    # need to pass cython the include directory so it can find the .pyx files
    cython_flags = ["-I", os.path.abspath("pele/potentials"), "-v"]
    debug_flags = []
    if build_type in ["Debug", "RelWithDebInfo", "MemCheck"]:
        debug_flags = ["--gdb", "--annotate", "-X", "linetrace=True", "-X", "boundscheck=True",
                       "-X", "wraparound=False", "-X", "cdivision=False"]
    cython3_flags = ["-X", "language_level=3", "-X", "c_string_type=unicode",
                     "-X", "c_string_encoding=utf-8"]
    p = subprocess.call(
        [sys.executable, os.path.join(cwd, "cythonize.py"), "pele"]
        + cython_flags + debug_flags + cython3_flags,
        cwd=cwd,
    )
    if p != 0:
        raise RuntimeError("Running cythonize failed!")


def which(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError("could not find " + name + " on PATH")
    return path


def get_compiler_env(compiler_id):
    """Environment and cmake args for the compilers.

    CC/CXX/FC from the environment (e.g. conda compilers) are respected. Otherwise
    gcc is used, on macOS the newest homebrew gcc-N.
    """
    env = os.environ.copy()
    cmake_args = ["-DCMAKE_EXPORT_COMPILE_COMMANDS=1"]
    # CMake only finds headers/libraries (sundials, eigen, lapack) under the active
    # conda env by itself when conda's own compilers are used
    prefixes = [env["CONDA_PREFIX"]] if env.get("CONDA_PREFIX") else []
    if compiler_id == "unix":
        if sys.platform.startswith("darwin") and "CC" not in env:
            version = next((v for v in range(20, 9, -1) if shutil.which(f"gcc-{v}")), None)
            if version is None:
                raise RuntimeError(
                    "Could not find a homebrew GNU compiler gcc-N (N=10..20) on PATH. "
                    "Install one or set CC, CXX and FC."
                )
            env["CC"] = which(f"gcc-{version}")
            env["CXX"] = which(f"g++-{version}")
            # f2py looks at F90, cmake at FC
            env["F90"] = env["FC"] = which(f"gfortran-{version}")
            prefixes += [subprocess.check_output(["brew", "--prefix", p]).decode(encoding).strip()
                         for p in ("openblas", "gettext")]
        env.setdefault("CC", "gcc")
        env.setdefault("CXX", "g++")
    elif compiler_id == "intel":
        env["CC"], env["CXX"], env["AR"] = which("icc"), which("icpc"), which("xiar")
        cmake_args.append("-DCMAKE_AR=" + env["AR"])
    else:
        raise Exception("compiler id not known")
    if prefixes:
        cmake_args.append("-DCMAKE_PREFIX_PATH=" + ";".join(prefixes))
    cmake_args += ["-DCMAKE_C_COMPILER=" + env["CC"], "-DCMAKE_CXX_COMPILER=" + env["CXX"]]
    if "FC" in env:
        cmake_args.append("-DCMAKE_Fortran_COMPILER=" + env["FC"])
    return env, cmake_args


def get_ldflags():
    """linker flags for libpython (only used on macOS, see CMakeLists.txt.in)"""
    getvar = sysconfig.get_config_var
    libs = (getvar("LIBS") or "").split() + (getvar("SYSLIBS") or "").split()
    if not getvar("Py_ENABLE_SHARED"):
        libs.insert(0, "-L" + getvar("LIBDIR"))
    if not getvar("PYTHONFRAMEWORK"):
        # See https://github.com/kovidgoyal/kitty/issues/289#issuecomment-416040645
        libs.extend((getvar("LINKFORSHARED") or "").replace("-Wl,-stack_size,1000000", "").split())
    return " ".join(libs)


def write_cmakelists():
    """create CMakeLists.txt from CMakeLists.txt.in"""
    import numpy as np

    with open("CMakeLists.txt.in") as fin:
        cmake_txt = fin.read()
    python_includes = {sysconfig.get_path("include"), sysconfig.get_path("platinclude")}
    cmake_txt = cmake_txt.replace("__PYTHON_INCLUDE__", " ".join(sorted(python_includes)))
    cmake_txt = cmake_txt.replace("__INCLUDE_SUNDIALS__", "ON" if with_cvode else "OFF")
    cmake_txt = cmake_txt.replace("__NUMPY_INCLUDE__", np.get_include())
    cmake_txt = cmake_txt.replace("__PYTHON_LDFLAGS__", get_ldflags())
    cmake_txt = cmake_txt.replace(
        "__COMPILER_EXTRA_ARGS__", '"{}"'.format(" ".join(cmake_compiler_extra_args))
    )
    with open("CMakeLists.txt", "w") as fout:
        fout.write(cmake_txt)
        fout.write("\n")
        for fname in cxx_files + c_files:
            fout.write("make_cython_lib(${CMAKE_CURRENT_SOURCE_DIR}/%s)\n" % fname)
        for fname in fortran_files:
            fout.write("make_fortran_lib(${CMAKE_CURRENT_SOURCE_DIR}/%s)\n" % fname)


def run_cmake():
    os.makedirs(cmake_build_dir, exist_ok=True)
    print("\nrunning cmake in directory", cmake_build_dir)
    cwd = os.path.abspath(os.path.dirname(__file__))
    env, cmake_args = get_compiler_env(idcompiler)
    cmake_args += ["-DPYTHON_EXECUTABLE=" + sys.executable]
    if shutil.which("ninja"):
        cache = os.path.join(cmake_build_dir, "CMakeCache.txt")
        # cmake refuses to switch generators in an existing build dir
        if os.path.isfile(cache) and "CMAKE_GENERATOR:INTERNAL=Ninja\n" not in open(cache).read():
            os.remove(cache)
            shutil.rmtree(os.path.join(cmake_build_dir, "CMakeFiles"), ignore_errors=True)
        cmake_args += ["-G", "Ninja"]
    if build_type == "Release":
        # CMake picks the LTO-aware archiver (gcc-ar) for the static pele_lib
        cmake_args += ["-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON"]
    if subprocess.call(["cmake"] + cmake_args + [cwd], cwd=cmake_build_dir, env=env) != 0:
        raise Exception("running cmake failed")
    print("\nbuilding files in cmake directory, jobs:", jargs.j)
    if subprocess.call(["cmake", "--build", ".", "-j", str(jargs.j)], cwd=cmake_build_dir, env=env) != 0:
        raise Exception("building libraries with CMake failed")
    print("finished building the extension modules with cmake\n")


class build_ext_precompiled(old_build_ext):
    """Build everything with CMake, then copy each library (stored in
    extension.sources[0]) to where setuptools expects the extension."""

    def run(self):
        generate_cython()
        write_cmakelists()
        run_cmake()
        super().run()

    def build_extension(self, ext):
        ext_path = self.get_ext_fullpath(ext.name)
        pre_compiled_library = ext.sources[0]
        if not os.path.isfile(pre_compiled_library):
            raise RuntimeError(
                "file does not exist: " + pre_compiled_library + " Did CMake not run correctly"
            )
        os.makedirs(os.path.dirname(ext_path), exist_ok=True)
        print("copying", pre_compiled_library, "to", ext_path)
        shutil.copy2(pre_compiled_library, ext_path)


class build_py_with_source(build_py):
    """also install the C++ sources/headers as pele/source, see pele.get_include()"""

    def run(self):
        super().run()
        dest = os.path.join(self.build_lib, "pele", "source")
        shutil.rmtree(dest, ignore_errors=True)
        shutil.copytree("source", dest, ignore=shutil.ignore_patterns("*.rst", "gmin"))


def module_name(fname):
    return os.path.splitext(fname)[0].replace("/", ".")


# The `source` of each extension is the location of the library built by CMake.
# f2py libraries carry the python ABI suffix, e.g. minperm.cpython-312-x86_64-linux-gnu.so
ext_modules = [
    Extension(module_name(f), [os.path.join(cmake_build_dir, os.path.splitext(os.path.basename(f))[0] + ".so")])
    for f in cxx_files + c_files
] + [
    Extension(module_name(f), [os.path.join(cmake_build_dir, os.path.splitext(os.path.basename(f))[0]
                                            + sysconfig.get_config_var("EXT_SUFFIX"))])
    for f in fortran_files
]

# written before setup() so build_py installs the current one
write_version_py()

# metadata lives in pyproject.toml
setup(
    packages=find_packages(include=["pele", "pele.*"]),
    package_data={"": ["*.xyz", "*.xyzdr", "*.data", "*.sqlite", "points*", "*.pxd", "*.pxi",
                      "coords", "data", "non_additive_test_data/*"]},
    ext_modules=ext_modules,
    cmdclass=dict(build_ext=build_ext_precompiled, build_py=build_py_with_source),
)
