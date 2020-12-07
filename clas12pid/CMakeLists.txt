# find source and header files, set paths, includes...
SET(LIBNAME "e_pid")
FILE(GLOB SRCS *.cpp)

FILE(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/e_pid.h)
set(INCLUDES ${ROOT_INCLUDE_DIR} "${CMAKE_SOURCE_DIR}/include" )
INCLUDE_DIRECTORIES(${INCLUDES})

# create root dictionary
#SET(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)
#SET(DICTIONARY ${LIBNAME}Dict.cpp)
#ROOT_GENERATE_DICTIONARY("${HEADERS}" "${LINKDEF}" "${DICTIONARY}" "${INCLUDES}")
#LIST(APPEND SRCS ${DICTIONARY})

# build the library
ADD_LIBRARY(${LIBNAME} SHARED ${SRCS} ${HEADERS})
TARGET_LINK_LIBRARIES(${LIBNAME} clashit_tree ${ROOT_LIBRARIES})
INSTALL ( TARGETS ${LIBNAME} DESTINATION "${CMAKE_SOURCE_DIR}/lib")

ADD_DEFINITIONS( -D_VAR="${CMAKE_SOURCE_DIR}/clas12pid" )

#Old message for debugging
#message("Searching for Libraries")
#find_library(Test clashit_tree)
#message(${Test})
#if(NOT Test)
#  message(FATAL_ERROR)
#endif()
