AT_source = custom_target(
    'ATmodule.c',
    input : ['AT.f90'],
    output : ['ATmodule.c', 'AT-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'AT', '--lower', '--build-dir', 'pele/potentials/fortran'],
)

py.extension_module(
    'AT',
    ['AT.f90', AT_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/potentials/fortran'
)

FinSin_source = custom_target(
    'FinSinmodule.c',
    input : ['FinSin.f90'],
    output : ['FinSinmodule.c', 'FinSin-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'FinSin', '--lower', '--build-dir', 'pele/potentials/fortran'],
)

py.extension_module(
    'FinSin',
    ['FinSin.f90', FinSin_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/potentials/fortran'
)

gupta_source = custom_target(
    'guptamodule.c',
    input : ['gupta.f90'],
    output : ['guptamodule.c', 'gupta-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'gupta', '--lower', '--build-dir', 'pele/potentials/fortran'],
)

py.extension_module(
    'gupta',
    ['gupta.f90', gupta_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/potentials/fortran'
)

lj_source = custom_target(
    'ljmodule.c',
    input : ['lj.f90'],
    output : ['ljmodule.c', 'lj-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'lj', '--lower', '--build-dir', 'pele/potentials/fortran'],
)

py.extension_module(
    'lj',
    ['lj.f90', lj_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/potentials/fortran'
)

lj_hess_source = custom_target(
    'lj_hessmodule.c',
    input : ['lj_hess.f90'],
    output : ['lj_hessmodule.c', 'lj_hess-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'lj_hess', '--lower', '--build-dir', 'pele/potentials/fortran'],
)

py.extension_module(
    'lj_hess',
    ['lj_hess.f90', lj_hess_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/potentials/fortran'
)

ljcut_source = custom_target(
    'ljcutmodule.c',
    input : ['ljcut.f90'],
    output : ['ljcutmodule.c', 'ljcut-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'ljcut', '--lower', '--build-dir', 'pele/potentials/fortran'],
)

py.extension_module(
    'ljcut',
    ['ljcut.f90', ljcut_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/potentials/fortran'
)

ljpshiftfort_source = custom_target(
    'ljpshiftfortmodule.c',
    input : ['ljpshiftfort.f90'],
    output : ['ljpshiftfortmodule.c', 'ljpshiftfort-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'ljpshiftfort', '--lower', '--build-dir', 'pele/potentials/fortran'],
)

py.extension_module(
    'ljpshiftfort',
    ['ljpshiftfort.f90', ljpshiftfort_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/potentials/fortran'
)

morse_source = custom_target(
    'morsemodule.c',
    input : ['morse.f90'],
    output : ['morsemodule.c', 'morse-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'morse', '--lower', '--build-dir', 'pele/potentials/fortran'],
)

py.extension_module(
    'morse',
    ['morse.f90', morse_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/potentials/fortran'
)

scdiff_periodic_source = custom_target(
    'scdiff_periodicmodule.c',
    input : ['scdiff_periodic.f90'],
    output : ['scdiff_periodicmodule.c', 'scdiff_periodic-f2pywrappers.f'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'scdiff_periodic', '--lower', '--build-dir', 'pele/potentials/fortran'],
)

py.extension_module(
    'scdiff_periodic',
    ['scdiff_periodic.f90', scdiff_periodic_source, fortranobject_c],
    include_directories: [inc_np,],
    dependencies: [py_dep,],
    install : true,
    subdir: 'pele/potentials/fortran'
)

