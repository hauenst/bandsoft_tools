find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})


#################################################################################################
add_executable(gainfit gainfit.cpp )
target_link_libraries(gainfit ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS gainfit
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################

