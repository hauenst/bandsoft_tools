project(hipolib)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(INCLUDES ${LZ4_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/hipolib )
include_directories(${INCLUDES})

add_library(Hipo4 
		bank.cpp dictionary.cpp reader.cpp record.cpp recordbuilder.cpp utils.cpp wrapper.cpp writer.cpp event.cpp )

target_link_libraries(Hipo4 ${LZ4_LIBRARY} )
install(TARGETS Hipo4 LIBRARY ARCHIVE DESTINATION  ${CMAKE_SOURCE_DIR}/lib  )

