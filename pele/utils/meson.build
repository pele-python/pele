_fortran_utils_source = custom_target(
  '_fortran_utilsmodule.c',
  input : ['_fortran_utils.f90'],
  output : ['_fortran_utilsmodule.c', '_fortran_utils-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_fortran_utils', '--lower', '--build-dir', 'pele/utils'],
)

py.extension_module(
  '_fortran_utils',
  ['_fortran_utils.f90', _fortran_utils_source, fortranobject_c],
  include_directories: [inc_np,],
  dependencies: [py_dep,],
  install : true,
  subdir: 'pele/utils'
)


py.extension_module(
    '_cpp_utils',
    ['_cpp_utils.pyx', '../../source/rotations.cpp'],
    include_directories: [inc_np, inc_pele],
    dependencies : [py_dep, pele_dep, cython_compile_dep,],
    install: true,
    subdir: 'pele/utils',
    override_options : ['cython_language=cpp'],
    cython_args : cython_args,
)


py.extension_module(
    '_pressure_tensor',
    ['_pressure_tensor.pyx', '../../source/pressure_tensor.cpp'],
    include_directories: [inc_np, inc_pele],
    dependencies : [py_dep, pele_dep, cython_compile_dep,],
    install: true,
    subdir: 'pele/utils',
    override_options : ['cython_language=cpp'],
    cython_args : cython_args,
)

