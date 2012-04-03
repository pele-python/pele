add_library(stockaa_ stockaa.cpp stockaa.f90 rigidbaa.f90)
target_link_libraries(stockaa_ boost_python)
set_target_properties(stockaa_ PROPERTIES PREFIX "")
