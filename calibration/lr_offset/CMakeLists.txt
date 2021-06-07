find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})


#################################################################################################
add_executable(lr_offset lr_offset.cpp )
target_link_libraries(lr_offset ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS lr_offset
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################

#################################################################################################
add_executable(test_offset test_offset.cpp )
target_link_libraries(test_offset ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS test_offset
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################

