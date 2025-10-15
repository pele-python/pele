minperm_source = custom_target(
    'minpermmodule.c',
    input : ['minperm.f90'],
    output : ['minpermmodule.c', 'minperm-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'minperm', '--lower', '--build-dir', 'pele/mindist'],
)

py.extension_module(
    'minperm',
    ['minperm.f90', minperm_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/mindist'
)

subdir('tests')
