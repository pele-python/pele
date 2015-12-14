cmake_minimum_required(VERSION 2.8)

project(pele_tests)

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release Coverage." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)
ADD_DEFINITIONS(-std=c++0x -g)
ADD_DEFINITIONS(-Wall)



#cmake_policy(SET CMP0015 NEW)

# Add and compile the gtest library
set(gtest_src ../gtest)
set(gtest_binary_dir googletest)
set(gtest_include ${gtest_src}/include/)
include_directories(${gtest_include})
message("gtest include directory ${gtest_include}")
add_subdirectory(../gtest ${gtest_binary_dir})
set(gtest_libs gtest gtest_main)
message("gtest_libs ${gtest_libs}")


# set the pele include directory
set(pele_include ../../source)
include_directories(${pele_include})
message("pele include directory ${pele_include}")

# build the pele library
FILE(GLOB pele_sources ${pele_include}/*.cpp)
add_library(pele_lib ${pele_sources})

# get all the source files
FILE(GLOB sources *.cpp)

add_executable(test_main ${sources})
target_link_libraries(test_main pele_lib gtest gtest_main pthread)

add_subdirectory(benchmarks)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
    include(CodeCoverage)
    setup_target_for_coverage(${PROJECT_NAME}_coverage test_main coverage)
    #target_link_libraries(test_main gcov)
    #target_link_libraries(pele_lib gcov)
endif()
