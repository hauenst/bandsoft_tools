find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})

#message("Searching for Libraries")
#find_library(Test clashit_tree)
#find_library(Test clashit_tree HINTS /w/hallb-scifs17exp/clas12/users/awild/BAND_Project/Build/clashit)
#message(${Test})
#if(NOT Test)
#  message(FATAL_ERROR)
#endif()

add_executable(EC_cuts EC_cuts.cpp)
target_link_libraries(EC_cuts ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS EC_cuts
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(SF_cuts_fit SF_cuts_fit.cpp)
target_link_libraries(SF_cuts_fit ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS SF_cuts_fit
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(PI_cuts_fit PI_cuts_fit.cpp)
target_link_libraries(PI_cuts_fit ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS PI_cuts_fit
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(cut_check cut_check.cpp)
target_link_libraries(cut_check ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS cut_check
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
