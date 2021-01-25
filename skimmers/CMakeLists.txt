find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})
set(CMAKE_VERBOSE_MAKEFILE ON)


#################################################################################################
add_executable(neutrons neutrons.cpp readhipo_helper.cpp )
target_link_libraries(neutrons ${ROOT_LIBRARIES} Hipo4 Banks4 
				clashit_tree bandhit_tree taghit_tree genpart_tree 
				clas12fiducial
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS neutrons
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
#################################################################################################


#################################################################################################
add_executable(electrons electrons.cpp readhipo_helper.cpp )
target_link_libraries(electrons ${ROOT_LIBRARIES} Hipo4 Banks4 
				clashit_tree bandhit_tree taghit_tree genpart_tree 
				clas12fiducial
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS electrons
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
#################################################################################################


#################################################################################################
add_executable(tagged_dis tagged_dis.cpp readhipo_helper.cpp )
target_link_libraries(tagged_dis ${ROOT_LIBRARIES} Hipo4 Banks4 
				clashit_tree bandhit_tree taghit_tree genpart_tree 
				clas12fiducial
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS tagged_dis
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
#################################################################################################


#################################################################################################
add_executable(ncalibration ncalibration.cpp readhipo_helper.cpp )
target_link_libraries(ncalibration ${ROOT_LIBRARIES} Hipo4 Banks4  
				clashit_tree bandhit_tree taghit_tree genpart_tree 
				clas12fiducial
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS ncalibration
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
#################################################################################################


#################################################################################################
add_executable(quasielastic quasielastic.cpp readhipo_helper.cpp )
target_link_libraries(quasielastic ${ROOT_LIBRARIES} Hipo4 Banks4 
				clashit_tree bandhit_tree taghit_tree genpart_tree 
				clas12fiducial
				-lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS quasielastic
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
#################################################################################################


#add_executable(EC_cuts EC_cuts.cpp readhipo_helper.cpp )
#target_link_libraries(EC_cuts ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
#install( TARGETS EC_cuts
#		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

#add_executable(SF_cuts_fit SF_cuts_fit.cpp readhipo_helper.cpp )
#target_link_libraries(SF_cuts_fit ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
#install( TARGETS SF_cuts_fit
#		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

#add_executable(SF_cuts_fit SF_cuts_fit.cpp readhipo_helper.cpp )
#target_link_libraries(SF_cuts_fit ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
#install( TARGETS SF_cuts_fit
#		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

#add_executable(DC_fid_test DC_fid_test.cpp)
#target_link_libraries(DC_fid_test ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial DC_fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
#install( TARGETS DC_fid_test
#		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
