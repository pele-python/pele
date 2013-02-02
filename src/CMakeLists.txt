cmake_minimum_required(VERSION 2.8)
list( APPEND CMAKE_CXX_FLAGS "-std=c++0x" )
 
enable_language(CXX)
add_library(pygmin neb.cpp)
add_executable(test test.cpp)
target_link_libraries(test pygmin)

