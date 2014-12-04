add_library(bench_lib SHARED bench_lib)
target_link_libraries(bench_lib pele_lib)

add_executable(bench_cell_lists bench_cell_lists.cpp)
target_link_libraries(bench_cell_lists pele_lib)
