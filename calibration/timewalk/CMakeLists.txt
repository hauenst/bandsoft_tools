find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})


#################################################################################################
add_executable(timewalk timewalk.cpp )
target_link_libraries(timewalk ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS timewalk
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################

#################################################################################################
add_executable(test_timewalk test_timewalk.cpp )
target_link_libraries(test_timewalk ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS test_timewalk
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################

#################################################################################################
add_executable(test_timewalk_v1 test_timewalk_v1.cpp )
target_link_libraries(test_timewalk_v1 ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS test_timewalk_v1
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################

#################################################################################################
add_executable(test_timewalk_v2 test_timewalk_v2.cpp )
target_link_libraries(test_timewalk_v2 ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS test_timewalk_v2
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################
