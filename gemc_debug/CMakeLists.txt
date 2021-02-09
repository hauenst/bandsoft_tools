find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})

add_executable(ep_resolutions ep_resolutions.cpp)
target_link_libraries(ep_resolutions ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree bandhit_tree genpart_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS ep_resolutions
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(neutron_resolutions neutron_resolutions.cpp)
target_link_libraries(neutron_resolutions ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree bandhit_tree genpart_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS neutron_resolutions
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(ep_gaps ep_gaps.cpp)
target_link_libraries(ep_gaps ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree bandhit_tree genpart_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS ep_gaps
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
