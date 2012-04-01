cmake_minimum_required(VERSION 2.8)
PROJECT(ljcpp)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
OPTION(BUILD_SHARED_LIBS "turn OFF for .a libs" ON)

FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

#INCLUDE_DIRECTORIES(/usr/local/include)
#target_link_libraries(dbookpy dbook)
#LINK_DIRECTORIES("/usr/local/lib")

ADD_LIBRARY(ljcpp lj.cpp)
target_link_libraries(ljcpp boost_python)

#SET_TARGET_PROPERTIES(dbookpy PROPERTIES LINKER_LANGUAGE C)
#iSET_TARGET_PROPERTIES(dbookpy
#PROPERTIES
#        SOVERSION 0.1
#        VERSION 0.1
#)

#INSTALL(TARGETS dbookpy
#        LIBRARY DESTINATION lib
#)

