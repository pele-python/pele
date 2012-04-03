cmake_minimum_required(VERSION 2.8)
project(mingmin)
option(BUILD_SHARED_LIBS "turn OFF for .a libs" ON)

enable_language(CXX)
enable_language(Fortran)

find_package(PythonInterp)
find_package(PythonLibs)

include_directories(${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_subdirectory(potentials)
