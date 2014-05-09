cmake_minimum_required(VERSION 2.8)

ADD_DEFINITIONS(-std=c++0x)

#cmake_policy(SET CMP0015 NEW)

# set the pele include directory
set(pele_include ../../source)
include_directories(${pele_include})
message("pele include directory ${pele_include}")

# build the lbfgs library
add_library(lbfgs_lib ${pele_include}/lbfgs.cpp)


# Add and compile the gtest library
set(gtest_src ../gtest)
set(gtest_binary_dir googletest)
set(gtest_include ${gtest_src}/include/)
include_directories(${gtest_include})
message("gtest include directory ${gtest_include}")
add_subdirectory(../gtest ${gtest_binary_dir})
set(gtest_libs gtest gtest_main)
message("gtest_libs ${gtest_libs}")

# compile the tests
add_executable(test_main test_main.cpp)
target_link_libraries(test_main lbfgs_lib ${gtest_libs} pthread)
