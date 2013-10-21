# So far build is with cmake.
# run run_cython.sh to generate all interface files
cmake_minimum_required(VERSION 2.8)

#list( APPEND CMAKE_CXX_FLAGS "-std=c++0x" )
enable_language(CXX)
enable_language(Fortran)

find_package(PythonLibs)

#file(GLOB sources include/)

if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message("sources ${sources}")

include_directories(include ${PYTHON_INCLUDE_PATH})

enable_language(CXX)

add_library(_lj SHARED _lj.cpp ${sources} ${PYTHON_LIBS})
set_target_properties(_lj PROPERTIES PREFIX "")

add_library(_pele SHARED _pele.cpp ${sources} ${PYTHON_LIBS})
set_target_properties(_pele PROPERTIES PREFIX "")

add_library(_lj_cython SHARED _lj_cython.cpp ../../pele/potentials/fortran/lj.f90 ${PYTHON_LIBS})
set_target_properties(_lj_cython PROPERTIES PREFIX "")

add_library(_lbfgs SHARED _lbfgs.cpp  include/_lbfgs.cpp ${PYTHON_LIBS})
set_target_properties(_lbfgs PROPERTIES PREFIX "")

add_library(_pythonpotential SHARED _pythonpotential.cpp  ${PYTHON_LIBS})
set_target_properties(_pythonpotential PROPERTIES PREFIX "")
