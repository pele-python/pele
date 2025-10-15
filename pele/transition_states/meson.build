_NEB_utils_source = custom_target(
    '_NEB_utilsmodule.c',
    input : ['_NEB_utils.f90'],
    output : ['_NEB_utilsmodule.c', '_NEB_utils-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_NEB_utils', '--lower', '--build-dir', 'pele/transition_states'],
)

py.extension_module(
    '_NEB_utils',
    ['_NEB_utils.f90', _NEB_utils_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/transition_states'
)

_orthogoptf_source = custom_target(
    '_orthogoptfmodule.c',
    input : ['_orthogoptf.f90'],
    output : ['_orthogoptfmodule.c', '_orthogoptf-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_orthogoptf', '--lower', '--build-dir', 'pele/transition_states'],
)

py.extension_module(
    '_orthogoptf',
    ['_orthogoptf.f90', _orthogoptf_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/transition_states'
)

subdir('tests')
