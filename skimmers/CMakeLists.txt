find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})

set(CMAKE_VERBOSE_MAKEFILE ON)
add_executable(neutrons neutrons.cpp readhipo_helper.cpp )
target_link_libraries(neutrons ${ROOT_LIBRARIES} Hipo4 Banks4  bandhit_tree -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS neutrons
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

