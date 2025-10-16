
py.extension_module(
    'mypotential',
    ['mypotential.pyx'],
    include_directories: [inc_np, inc_pele],
    dependencies : [py_dep, pele_dep, cython_compile_dep,],
    install: true,
    subdir: 'new_potential/cpp_potential',
    cython_args : cython_args,
)


