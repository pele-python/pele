_spherical_container_source = custom_target(
  '_spherical_containermodule.c',
  input : ['_spherical_container.f90'],
  output : ['_spherical_containermodule.c', '_spherical_container-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_spherical_container', '--lower', '--build-dir', 'pele/accept_tests'],
)

py.extension_module(
  '_spherical_container',
  ['_spherical_container.f90', _spherical_container_source, fortranobject_c],
  include_directories: [inc_np,],
  dependencies: [py_dep,],
  install : true,
  subdir: 'pele/accept_tests'
)

