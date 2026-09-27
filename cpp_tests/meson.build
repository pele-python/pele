# C++ tests: meson setup build -Dpython=disabled -Dtests=true -Dbuildtype=debug
#            meson test -C build
# add -Db_sanitize=address or -Db_coverage=true as needed.

if not sundials_dep.found()
  error('the C++ tests need SUNDIALS (-Dcvode=enabled)')
endif

# the pinned submodule if checked out, else a system GoogleTest
if fs.is_file('gtest/googletest/src/gtest-all.cc')
  gtest_inc = include_directories('gtest/googletest/include', 'gtest/googletest', is_system: true)
  gtest_lib = static_library(
    'gtest',
    'gtest/googletest/src/gtest-all.cc',
    'gtest/googletest/src/gtest_main.cc',
    include_directories: gtest_inc,
    dependencies: dependency('threads'),
  )
  gtest_dep = declare_dependency(
    link_with: gtest_lib,
    include_directories: gtest_inc,
    dependencies: dependency('threads'),
  )
else
  gtest_dep = dependency('gtest', main: true)
endif

test_main = executable(
  'test_main',
  'source/test_aatopology.cpp',
  'source/test_array.cpp',
  'source/test_base_potential.cpp',
  'source/test_blj.cpp',
  'source/test_cell_list_minimize.cpp',
  'source/test_cell_lists.cpp',
  'source/test_cell_lists_hs_wca.cpp',
  'source/test_check_copy_assignment_cvode.cpp',
  'source/test_combinations.cpp',
  'source/test_combined_potential.cpp',
  'source/test_cvode.cpp',
  'source/test_distance.cpp',
  'source/test_eigenvalues.cpp',
  'source/test_energy_change.cpp',
  'source/test_fire.cpp',
  'source/test_frozen.cpp',
  'source/test_gaussianpot.cpp',
  'source/test_generic_mixed_descent.cpp',
  'source/test_gradient_descent.cpp',
  'source/test_harmonic.cpp',
  'source/test_hoomd.cpp',
  'source/test_inversepower.cpp',
  'source/test_inversepower_hs.cpp',
  'source/test_inversepower_stillinger.cpp',
  'source/test_inversepower_stillinger_cut.cpp',
  'source/test_inversepower_stillinger_cut_quad.cpp',
  'source/test_lbfgs.cpp',
  'source/test_lj.cpp',
  'source/test_lowest_eig_potential.cpp',
  'source/test_matrix.cpp',
  'source/test_mixed_descent.cpp',
  'source/test_mixed_descent_extended_potential.cpp',
  'source/test_morse.cpp',
  'source/test_newton.cpp',
  'source/test_ngt.cpp',
  'source/test_pressure.cpp',
  'source/test_pspin.cpp',
  'source/test_rattlers.cpp',
  'source/test_rotations.cpp',
  'source/test_utils.cpp',
  'source/test_vecn.cpp',
  'source/test_wca.cpp',
  # not in the Python build yet, see pele/optimize
  '../source/newton_trust.cpp',
  dependencies: [pele_dep, gtest_dep, dependency('blas', required: false)],
)
# the tests write their output files (trajectories etc.) to the working directory
test('pele C++ tests', test_main, workdir: meson.current_build_dir(), timeout: 1800)

# benchmarks: ninja -C build bench_lj ...
foreach b : ['bench_cell_lists', 'bench_cell_lists_alt', 'bench_lj', 'bench_inversepower', 'bench_minimization']
  executable(b, 'source/benchmarks' / b + '.cpp', dependencies: pele_dep, build_by_default: false)
endforeach
