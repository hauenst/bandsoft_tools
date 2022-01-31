find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})



#################################################################################################
add_executable(paddle_class paddle_class.cpp ../../bandreco/bandreco.cpp )
target_link_libraries(paddle_class ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS paddle_class
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################

#################################################################################################
add_executable(layer_class layer_class.cpp ../../bandreco/bandreco.cpp )
target_link_libraries(layer_class ${ROOT_LIBRARIES} Hipo4 Banks4  
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS layer_class
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
add_executable(initial_bar initial_bar.cpp ../../bandreco/bandreco.cpp )
target_link_libraries(initial_bar ${ROOT_LIBRARIES} Hipo4 Banks4 
				clashit_tree bandhit_tree taghit_tree genpart_tree 
				e_pid DC_fiducial
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS initial_bar
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/calibration )
#################################################################################################
