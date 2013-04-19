if(WITH_BOOST)
  add_library(salt_ salt.cpp matrix.cpp)
  target_link_libraries(salt_ ${Boost_PYTHON_LIBRARY})
  set_target_properties(salt_ PROPERTIES PREFIX "")
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

