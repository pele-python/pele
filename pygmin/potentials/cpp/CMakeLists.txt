if(WITH_BOOST)
  add_library(ljcpp_ lj.cpp)
  target_link_libraries(ljcpp_ ${Boost_PYTHON_LIBRARY})
  set_target_properties(ljcpp_ PROPERTIES PREFIX "")
endif(WITH_BOOST)

#SET_TARGET_PROPERTIES(dbookpy PROPERTIES LINKER_LANGUAGE C)
#iSET_TARGET_PROPERTIES(dbookpy
#PROPERTIES
#        SOVERSION 0.1
#        VERSION 0.1
#)

#INSTALL(TARGETS dbookpy
#        LIBRARY DESTINATION lib
#)

