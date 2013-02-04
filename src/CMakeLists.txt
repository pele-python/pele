# So far build is with cmake.
# run run_cython.sh to generate all interface files
cmake_minimum_required(VERSION 2.8)
list( APPEND CMAKE_CXX_FLAGS "-std=c++0x" )
find_package(PythonLibs)

file(GLOB sources _pygmin.cpp ../pygmin/potentials/fortran/lj.f90 )

message("sources ${sources}")

include_directories(${PYTHON_INCLUDE_PATH})

enable_language(CXX)
enable_language(Fortran)

add_library(_pygmin SHARED _pygmin.cpp ${sources} ${PYTHON_LIBS})
set_target_properties(_pygmin PROPERTIES PREFIX "")
add_library(_lj SHARED _lj.cpp lj.cpp ${sources} ${PYTHON_LIBS})
set_target_properties(_lj PROPERTIES PREFIX "")
add_library(_lj_cython SHARED _lj_cython.cpp ${sources} ${PYTHON_LIBS})
set_target_properties(_lj_cython PROPERTIES PREFIX "")
