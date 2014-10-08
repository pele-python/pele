cmake_minimum_required(VERSION 2.8)

find_package(PythonLibs)
message(${PYTHON_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})
set(numpy_include "/home/js850/.local/lib/python2.7/site-packages/numpy/core/include")
include_directories(${numpy_include})

enable_language(CXX)
SET(CMAKE_CXX_FLAGS "-Wall -ansi -pedantic -std=c++0x")

#cmake_policy(SET CMP0015 NEW)

# set the pele include directory
set(pele_include ${CMAKE_SOURCE_DIR}/source)
include_directories(${pele_include})
message("pele include directory ${pele_include}")

# build the pele library
FILE(GLOB pele_sources ${pele_include}/*.cpp)
add_library(pele_lib SHARED ${pele_sources})

set(ngt_source ${CMAKE_SOURCE_DIR}/pele/rates/_ngt_cpp.cxx)
#add_library(_ngt_cpp SHARED ${ngt_source} ${pele_sources} ${PYTHON_LIBS})
add_library(_ngt_cpp SHARED ${ngt_source})
target_link_libraries(_ngt_cpp pele_lib)
set_target_properties(_ngt_cpp PROPERTIES PREFIX "")
