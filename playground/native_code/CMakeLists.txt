# So far build is with cmake.
# run run_cython.sh to generate all interface files
cmake_minimum_required(VERSION 2.8)
list( APPEND CMAKE_CXX_FLAGS "-std=c++0x" )
find_package(PythonLibs)

#file(GLOB sources _pygmin.cpp _lj.cpp)

if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message("sources ${sources}")

include_directories(include ${PYTHON_INCLUDE_PATH})

enable_language(CXX)

add_library(_lj SHARED _lj.cpp ${sources} ${PYTHON_LIBS})
set_target_properties(_lj PROPERTIES PREFIX "")

add_library(_pygmin SHARED _pygmin.cpp ${sources} ${PYTHON_LIBS})
set_target_properties(_pygmin PROPERTIES PREFIX "")

include_directories(lbfgs)

add_library(_lbfgs SHARED _lbfgs.cpp lbfgs_wrapper.cpp lbfgs/lbfgs.c ${PYTHON_LIBS})
set_target_properties(_lbfgs PROPERTIES PREFIX "")
