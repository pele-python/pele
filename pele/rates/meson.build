
py.extension_module(
    '_ngt_cpp',
    ['_ngt_cpp.pyx', '../../source/pele/graph.hpp', '../../source/pele/ngt.hpp'],
    include_directories: [inc_np, inc_pele],
    dependencies : [py_dep, pele_dep, cython_compile_dep,],
    install: true,
    subdir: 'pele/rates',
    override_options : ['cython_language=cpp'],
    cython_args : cython_args,
)

subdir('tests')
