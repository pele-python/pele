project(
  'pele',
  'c', 'cpp',
  version: '0.1.0',
  meson_version: '>=1.3',
  default_options: [
    'buildtype=release',
    'b_ndebug=if-release',
    'cpp_std=c++20',
  ],
)

fs = import('fs')
cc = meson.get_compiler('cpp')
# meson always runs under a Python; used for configure-time helpers only
helper_py = import('python').find_installation()
with_python = get_option('python').allowed()

if host_machine.system() == 'darwin' and cc.get_id() != 'gcc'
  error('pele needs GCC on macOS (Apple clang has no OpenMP or Fortran): brew install gcc, then build with CC=gcc-15 CXX=g++-15 FC=gfortran-15')
endif

native_args = get_option('native') and get_option('buildtype') != 'debug' ? ['-march=native'] : []
add_project_arguments(native_args, language: ['c', 'cpp', 'fortran'])
omp_dep = dependency('openmp')

# SUNDIALS, Eigen and LAPACK come from the active conda env if there is one, else the system.
# The env's shared libraries are found at run time through an rpath.
conda_prefix = run_command(helper_py, '-c', 'import os; print(os.environ.get("CONDA_PREFIX", ""))', check: true).stdout().strip()
prefixes = conda_prefix != '' ? [conda_prefix] : []
env_rpath = conda_prefix != '' ? conda_prefix / 'lib' : ''

eigen_inc = []
foreach p : prefixes
  foreach d : [p / 'include', p / 'include' / 'eigen3']
    if eigen_inc.length() == 0 and fs.is_dir(d / 'Eigen')
      eigen_inc = [d]
    endif
  endforeach
endforeach
if eigen_inc.length() > 0
  eigen_dep = declare_dependency(include_directories: include_directories(eigen_inc, is_system: true))
else
  eigen_dep = dependency('eigen3')
endif

pele_lib_sources = files(
  'source/aatopology.cpp',
  'source/atlj.cpp',
  'source/backtracking.cpp',
  'source/bracketing.cpp',
  'source/lbfgs.cpp',
  'source/linesearch.cpp',
  'source/modified_fire.cpp',
  'source/more_thuente.cpp',
  'source/nwpele.cpp',
  'source/pressure_tensor.cpp',
  'source/rotations.cpp',
  'source/xsum.c',
  'source/pbinary.c',
)
pele_lib_deps = [eigen_dep, omp_dep]

sundials_dep = disabler()
if not get_option('cvode').disabled()
  lib_dirs = []
  inc_dirs = []
  foreach p : prefixes
    lib_dirs += p / 'lib'
    inc_dirs += p / 'include'
  endforeach
  cvode_lib = cc.find_library('sundials_cvode', dirs: lib_dirs, required: get_option('cvode'))
  if cvode_lib.found()
    sundials_libs = [cvode_lib, cc.find_library('sundials_nvecserial', dirs: lib_dirs)]
    # split out in sundials 7
    sundials_libs += cc.find_library('sundials_core', dirs: lib_dirs, required: false)
    sundials_inc = []
    foreach d : inc_dirs
      if sundials_inc.length() == 0 and fs.is_file(d / 'cvode' / 'cvode.h')
        sundials_inc = [d]
      endif
    endforeach
    # the mixed descent optimizers call LAPACK's Cholesky (dpotrf_, dpotrs_)
    lapack_dep = cc.find_library('lapack', dirs: lib_dirs, required: false)
    if not lapack_dep.found()
      lapack_dep = dependency('lapack', 'openblas')
    endif
    sundials_dep = declare_dependency(
      dependencies: sundials_libs + [lapack_dep],
      include_directories: include_directories(sundials_inc, is_system: true),
    )
    # the Fortran prototypes for those calls. lapack.h ships with the LAPACKE headers,
    # but only liblapack is linked: nothing calls the LAPACKE_* C interface.
    if not cc.has_header('lapack.h', dependencies: sundials_dep)
      error('lapack.h not found. Install it with conda (blas-devel), apt (liblapacke-dev) or homebrew (openblas).')
    endif
    pele_lib_sources += files(
      'source/cvode.cpp',
      'source/extended_mixed_descent.cpp',
      'source/mxd_end_only.cpp',
      'source/newton.cpp',
      'source/newton_with_extension.cpp',
      'source/generic_mixed_descent.cpp',
    )
    pele_lib_deps += sundials_dep
  endif
endif

pele_inc = include_directories('source', is_system: true)
pele_lib = static_library(
  'pele_lib',
  pele_lib_sources,
  include_directories: pele_inc,
  dependencies: pele_lib_deps,
  cpp_args: ['-Wall', '-Wextra', '-pedantic'],
  pic: true,
)
pele_dep = declare_dependency(
  link_with: pele_lib,
  include_directories: pele_inc,
  dependencies: pele_lib_deps,
)

if with_python
add_languages('cython', 'fortran', native: false)
py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# Extensions built by a g++ newer than the libstdc++ that Python loads (e.g. system gcc
# with an older conda env) fail at import with "CXXABI_x not found". Refuse to build.
if host_machine.system() == 'linux' and cc.get_id() == 'gcc'
  compiler_libstdcxx = run_command(cc.cmd_array(), '-print-file-name=libstdc++.so.6', check: true).stdout().strip()
  libstdcxx_check = run_command(py, '-c', '''
import ctypes, os, sys
ctypes.CDLL("libstdc++.so.6")
loaded = next(l.split()[-1] for l in open("/proc/self/maps") if "libstdc++.so" in l)
version = lambda p: tuple(int(x) for x in os.path.realpath(p).rsplit(".so.", 1)[1].split("."))
if version(sys.argv[1]) > version(loaded):
    print(os.path.realpath(sys.argv[1]), os.path.realpath(loaded))
''', compiler_libstdcxx, check: true).stdout().strip()
  if libstdcxx_check != ''
    paths = libstdcxx_check.split(' ')
    error('''The C++ compiler's runtime @0@ is newer than the one Python loads, @1@,
so pele would fail to import. Update the environment's runtime to match the compiler,
e.g. conda install -c conda-forge "libstdcxx-ng>=$(gcc -dumpversion)".'''.format(paths[0], paths[1]))
  endif
endif

incdir_numpy = run_command(py, '-c', 'import numpy; print(numpy.get_include())', check: true).stdout().strip()
incdir_f2py = run_command(py, '-c', 'import numpy.f2py; print(numpy.f2py.get_include())', check: true).stdout().strip()
np_dep = declare_dependency(
  include_directories: include_directories(incdir_numpy),
  compile_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'],
)

cython_args = ['-X', 'c_string_type=unicode', '-X', 'c_string_encoding=utf-8',
               '-I', meson.project_source_root()]

# Cython extensions (C++ unless listed as C). The cvode ones need SUNDIALS.
cython_c = ['pele/optimize/_cython_lbfgs', 'pele/potentials/_cython_tools']
cython_cpp = [
  'pele/angleaxis/_cpp_aa',
  'pele/distance/_get_distance_cpp',
  'pele/distance/_put_in_box_cpp',
  'pele/optimize/_gradient_descent_cpp',
  'pele/optimize/_cosine_gradient_descent_cpp',
  'pele/optimize/_lbfgs_cpp',
  'pele/optimize/_modified_fire_cpp',
  'pele/optimize/_pele_opt',
  'pele/potentials/_frenkel',
  'pele/potentials/_frozen_dof',
  'pele/potentials/_harmonic_cpp',
  'pele/potentials/_hs_wca_cpp',
  'pele/potentials/_inversepower_cpp',
  'pele/potentials/_inversepower_hs_cpp',
  'pele/potentials/_inversepower_stillinger_cpp',
  'pele/potentials/_inversepower_stillinger_cut_cpp',
  'pele/potentials/_inversepower_stillinger_cut_quad',
  'pele/potentials/_lj_cpp',
  'pele/potentials/_morse_cpp',
  'pele/potentials/_pele',
  'pele/potentials/_pspin_spherical_cpp',
  'pele/potentials/_pythonpotential',
  'pele/potentials/_radial_gaussian_cpp',
  'pele/potentials/_sumgaussianpot_cpp',
  'pele/potentials/_wca_cpp',
  'pele/potentials/atlj',
  'pele/potentials/combine_potentials',
  'pele/potentials/cpp_test_functions',
  'pele/rates/_ngt_cpp',
  'pele/utils/_cpp_utils',
  'pele/utils/_pressure_tensor',
]
cython_cvode = [
  'pele/optimize/_mxd_end_only',
  'pele/optimize/cvode_opt',
  'pele/optimize/extended_mixed_descent',
  'pele/optimize/generic_mixed_descent',
]
if sundials_dep.found()
  cython_cpp += cython_cvode
endif

foreach m : cython_c + cython_cpp
  py.extension_module(
    fs.name(m),
    m + '.pyx',
    subdir: fs.parent(m),
    dependencies: [pele_dep, py_dep, np_dep],
    cython_args: cython_args,
    override_options: m in cython_c ? [] : ['cython_language=cpp'],
    build_rpath: env_rpath,
    install_rpath: env_rpath,
    install: true,
  )
endforeach

# Compiled LAMMPS potential, linked against the liblammps shipped with the lammps
# Python package (pip install lammps[mpi]). Without it, pele.potentials falls back to
# the pure Python LAMMPSPotential whenever lammps is importable.
lammps_dirs = run_command(py, '-c', '''
import pathlib, sys
try:
    import lammps
except ImportError:
    raise SystemExit
pkg = pathlib.Path(lammps.__file__).parent
inc = [d for d in (pkg / "include" / "lammps", pkg.parents[3] / "include" / "lammps") if d.is_dir()]
# the mpich wheel that lammps[mpi] pulls in puts mpi.h in the environment's include/
if inc:
    print(pkg, inc[0], pathlib.Path(sys.prefix) / "include")
''', check: true).stdout().strip().split(' ')
lammps_found = false
if lammps_dirs.length() == 3
  lammps_lib = cc.find_library('lammps', dirs: lammps_dirs[0], required: get_option('lammps'))
  lammps_inc = include_directories(lammps_dirs[1], 'source/lammps', lammps_dirs[2], is_system: true)
  lammps_found = lammps_lib.found() and cc.has_header(
    'mpi.h', include_directories: lammps_inc, required: get_option('lammps'))
elif get_option('lammps').enabled()
  error('lammps=enabled needs the lammps Python package with its headers: pip install "lammps[mpi]"')
endif
if lammps_found
  py.extension_module(
    'lammps_pele_cython',
    ['pele/potentials/lammps_pele_cython.pyx', 'source/lammps_pele.cpp'],
    subdir: 'pele/potentials',
    include_directories: lammps_inc,
    dependencies: [pele_dep, py_dep, np_dep, lammps_lib],
    cython_args: cython_args,
    override_options: ['cython_language=cpp'],
    build_rpath: lammps_dirs[0],
    install_rpath: lammps_dirs[0],
    install: true,
  )
endif

# f2py modules: f2py writes <name>module.c and <name>-f2pywrappers.f for each
fortranobject = static_library(
  'fortranobject',
  incdir_f2py / 'fortranobject.c',
  dependencies: [py_dep, np_dep],
  include_directories: include_directories(incdir_f2py),
  pic: true,
)
f2py_dep = declare_dependency(
  link_with: fortranobject,
  include_directories: include_directories(incdir_f2py),
  dependencies: [py_dep, np_dep],
)
fortran_modules = [
  'pele/mindist/minperm',
  'pele/optimize/mylbfgs_updatestep',
  'pele/potentials/fortran/AT',
  'pele/potentials/fortran/ljpshiftfort',
  'pele/potentials/fortran/lj',
  'pele/potentials/fortran/ljcut',
  'pele/potentials/fortran/lj_hess',
  'pele/potentials/fortran/morse',
  'pele/potentials/fortran/scdiff_periodic',
  'pele/potentials/fortran/FinSin',
  'pele/potentials/fortran/gupta',
  'pele/utils/_fortran_utils',
  'pele/transition_states/_orthogoptf',
  'pele/transition_states/_NEB_utils',
  'pele/angleaxis/_aadist',
  'pele/accept_tests/_spherical_container',
]
foreach m : fortran_modules
  name = fs.name(m)
  wrappers = custom_target(
    name + '_f2py',
    input: m + '.f90',
    output: [name + 'module.c', name + '-f2pywrappers.f'],
    command: [py, '-m', 'numpy.f2py', '@INPUT@', '-m', name, '--lower', '--build-dir', '@OUTDIR@', '--quiet'],
  )
  py.extension_module(
    name,
    [m + '.f90', wrappers],
    subdir: fs.parent(m),
    dependencies: [f2py_dep],
    link_language: 'fortran',
    install: true,
  )
endforeach

# pure Python files and package data
package_files = run_command(py, '-c', '''
import os
data = (".py", ".pxd", ".pxi", ".xyz", ".xyzdr", ".data", ".sqlite")
for root, dirs, files in os.walk("pele"):
    dirs[:] = sorted(d for d in dirs if d != "__pycache__")
    for f in sorted(files):
        p = os.path.join(root, f)
        if p == os.path.join("pele", "version.py"):
            continue
        if (f.endswith(data) or f.startswith("points") or f in ("coords", "data")
                or os.path.basename(root) == "non_additive_test_data"):
            print(p)
''', check: true).stdout().strip().split('\n')
py.install_sources(package_files, preserve_path: true)

# GUI: the PyQt5 forms and resources are generated from the .ui/.qrc files
pyuic = find_program('pyuic5', required: get_option('gui'))
pyrcc = find_program('pyrcc5', required: get_option('gui'))
if pyuic.found() and pyrcc.found()
  foreach ui : [
    'pele/gui/MainWindow',
    'pele/gui/connect_explorer_ui',
    'pele/gui/ui_params',
    'pele/gui/ui/connect_run_ui',
    'pele/gui/ui/cv_viewer_ui',
    'pele/gui/ui/dgraph_browser',
    'pele/gui/ui/graph_view_ui',
    'pele/gui/ui/nebbrowser',
    'pele/gui/ui/rate_gui',
    'pele/gui/ui/show3d_with_slider_ui',
    'pele/gui/ui/ui_neb_explorer',
    'pele/gui/ui/ui_normalmode_explorer',
    'pele/gui/ui/ui_takestep_explorer',
  ]
    custom_target(
      ui.underscorify(),
      input: ui + '.ui',
      output: fs.name(ui) + '.py',
      command: [pyuic, '@INPUT@', '-o', '@OUTPUT@'],
      install: true,
      install_dir: py.get_install_dir() / fs.parent(ui),
    )
  endforeach
  custom_target(
    'resources_rc',
    input: 'pele/gui/ui/resources.qrc',
    output: 'resources_rc.py',
    command: [pyrcc, '@INPUT@', '-o', '@OUTPUT@'],
    install: true,
    install_dir: py.get_install_dir() / 'pele/gui/ui',
  )
endif

# C++ sources and headers for extensions built against pele (see pele.get_include())
install_subdir('source', install_dir: py.get_install_dir() / 'pele',
               exclude_directories: ['gmin'], exclude_files: ['README.rst', 'CODING_STYLE.rst'])

git_revision = run_command('git', 'rev-parse', 'HEAD', check: false)
configure_file(
  output: 'version.py',
  # layout (git_revision on the 3rd line) is parsed by downstream scripts
  command: [py, '-c', 'import sys; open(sys.argv[1], "w").write("\\n# generated by meson.build\\ngit_revision = \'" + sys.argv[2] + "\'\\n")',
            '@OUTPUT@', git_revision.returncode() == 0 ? git_revision.stdout().strip() : 'Unknown'],
  install_dir: py.get_install_dir() / 'pele',
)
endif

if get_option('tests')
  subdir('cpp_tests')
endif
