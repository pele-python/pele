_aadist_source = custom_target(
    '_aadistmodule.c',
    input : ['_aadist.f90'],
    output : ['_aadistmodule.c', '_aadist-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_aadist', '--lower', '--build-dir', 'pele/angleaxis'],
)

py.extension_module(
    '_aadist',
    ['_aadist.f90', _aadist_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/angleaxis'
)


py.extension_module(
    '_cpp_aa',
    ['_cpp_aa.pyx', '../../pele/sources/aatopology.cpp', '../../pele/sources/rotations.cpp'],
    include_directories: [inc_np, inc_pele],
    dependencies : [py_dep, pele_dep, cython_compile_dep,],
    install: true,
    subdir: 'pele/angleaxis',
    override_options : ['cython_language=cpp'],
    cython_args : cython_args,
)

