project('pele',
    ['c', 'cpp', 'fortran', 'cython'],
    version : '0.1.1',
    meson_version: '>= 1.1.0',
    default_options : ['warning_level=1', 'buildtype=release'])
fc = meson.get_compiler('fortran')

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

inc_np = include_directories(incdir_numpy, incdir_f2py)
quadmath_dep = fc.find_library('quadmath', required: false)

inc_pele = include_directories('./source')
pele_dep = declare_dependency(include_directories: inc_pele)

cython_args = ['-Xembedsignature=True','-3']
cython_compile_dep = declare_dependency(compile_args: ['-Wall', '-Wextra', '-pedantic', '-funroll-loops', '-O3',]) 

pyuic = find_program('pyuic5')
pyrcc = find_program('pyrcc5')
pydir = py.get_install_dir()
subdir('pele')
