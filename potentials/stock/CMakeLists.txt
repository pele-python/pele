cmake_minimum_required(VERSION 2.8)
PROJECT(ljcpp)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
OPTION(BUILD_SHARED_LIBS "turn OFF for .a libs" ON)

enable_language(CXX)
enable_language(Fortran)

FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

#INCLUDE_DIRECTORIES(/usr/local/include)
#target_link_libraries(dbookpy dbook)
#LINK_DIRECTORIES("/usr/local/lib")

ADD_LIBRARY(stockaa_ stockaa.cpp stockaa.f90 rigidbaa.f90)
target_link_libraries(stockaa_ boost_python)
set_target_properties(stockaa_ PROPERTIES PREFIX "")

#SET_TARGET_PROPERTIES(dbookpy PROPERTIES LINKER_LANGUAGE C)
#iSET_TARGET_PROPERTIES(dbookpy
#PROPERTIES
#        SOVERSION 0.1
#        VERSION 0.1
#)

#INSTALL(TARGETS dbookpy
#        LIBRARY DESTINATION lib
#)

