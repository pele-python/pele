mylbfgs_updatestep_source = custom_target(
    'mylbfgs_updatestepmodule.c',
    input : ['mylbfgs_updatestep.f90'],
    output : ['mylbfgs_updatestepmodule.c', 'mylbfgs_updatestep-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'mylbfgs_updatestep', '--lower', '--build-dir', 'pele/optimize'],
)

py.extension_module(
    'mylbfgs_updatestep',
    ['mylbfgs_updatestep.f90', mylbfgs_updatestep_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/optimize'
)


py.extension_module(
    '_cython_lbfgs',
    ['_cython_lbfgs.pyx', ],
    include_directories: [inc_np, inc_pele],
    dependencies : [py_dep, pele_dep, cython_compile_dep,],
    install: true,
    subdir: 'pele/optimize',
    override_options : ['cython_language=cpp'],
    cython_args : cython_args,
)


py.extension_module(
    '_lbfgs_cpp',
    ['_lbfgs_cpp.pyx', '../../source/lbfgs.cpp'],
    include_directories: [inc_np, inc_pele],
    dependencies : [py_dep, pele_dep, cython_compile_dep,],
    install: true,
    subdir: 'pele/optimize',
    override_options : ['cython_language=cpp'],
    cython_args : cython_args,
)


py.extension_module(
    '_modified_fire_cpp',
    ['_modified_fire_cpp.pyx', '../../source/modified_fire.cpp'],
    include_directories: [inc_np, inc_pele],
    dependencies : [py_dep, pele_dep, cython_compile_dep,],
    install: true,
    subdir: 'pele/optimize',
    override_options : ['cython_language=cpp'],
    cython_args : cython_args,
)


py.extension_module(
    '_pele_opt',
    ['_pele_opt.pyx', ],
    include_directories: [inc_np, inc_pele],
    dependencies : [py_dep, pele_dep, cython_compile_dep,],
    install: true,
    subdir: 'pele/optimize',
    override_options : ['cython_language=cpp'],
    cython_args : cython_args,
)

