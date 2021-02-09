find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})

add_executable(resolution resolution.cpp)
target_link_libraries(resolution ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree bandhit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS resolution
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(nspectrum nspectrum.cpp)
target_link_libraries(nspectrum ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree bandhit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS nspectrum
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(nefficiency nefficiency.cpp)
target_link_libraries(nefficiency ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree bandhit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS nefficiency
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(n_skim n_skim.cpp)
target_link_libraries(n_skim ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree bandhit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS n_skim
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(check_nalgo check_nalgo.cpp)
target_link_libraries(check_nalgo ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree bandhit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS check_nalgo
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(epplus epplus.cpp)
target_link_libraries(epplus ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree bandhit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS epplus
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
