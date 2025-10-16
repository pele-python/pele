_mypotential_source = custom_target(
    '_mypotentialmodule.c',
    input : ['_mypotential.f90'],
    output : ['_mypotentialmodule.c', '_mypotential-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_mypotential', '--lower', '--build-dir', 'new_potential/fortran_potential'],
)

py.extension_module(
    '_mypotential',
    ['_mypotential.f90', _mypotential_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'new_potential/fortran_potential'
)

