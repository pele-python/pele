cmake_minimum_required(VERSION 2.8)

ADD_DEFINITIONS(-std=c++0x)

#cmake_policy(SET CMP0015 NEW)

# set the pele include directory
set(pele_include ../../source)
message("pele include directory ${pele_include}")

# add gtest include directory
set(gtest_include gtest/include/)

include_directories(${pele_include} ${gtest_include})

# compile the gtest library
message("gtest include directory ${gtest_include}")
add_subdirectory(gtest)
set(gtest_libs gtest gtest_main)

# build the lbfgs library
add_library(lbfgs_lib ${pele_include}/lbfgs.cpp)

# compile the tests
add_executable(test_main test_array.hpp test_main.cpp)
target_link_libraries(test_main lbfgs_lib ${gtest_libs} pthread)
