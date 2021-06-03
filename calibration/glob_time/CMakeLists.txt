find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})



#################################################################################################
add_executable(paddle paddle.cpp )
target_link_libraries(paddle ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS paddle
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################

##################################################################################################
#add_executable(global global.cpp )
#target_link_libraries(global ${ROOT_LIBRARIES} Hipo4 Banks4  
#				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
#install( TARGETS global
#		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
##################################################################################################


#################################################################################################
add_executable(test_class test_class.cpp ../../bandreco/bandreco.cpp )
target_link_libraries(test_class ${ROOT_LIBRARIES} Hipo4 Banks4 
				clashit_tree bandhit_tree taghit_tree genpart_tree 
				clas12fiducial e_pid DC_fiducial
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS test_class
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################
