project(banklib)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib )
include_directories(${INCLUDES})

add_library(Banks4  
		BBand.cpp BCalorimeter.cpp BConfig.cpp BEvent.cpp BParticle.cpp BScaler.cpp BScintillator.cpp)
target_link_libraries(Banks4 Hipo4 )
install(TARGETS Banks4 LIBRARY ARCHIVE DESTINATION  ${CMAKE_SOURCE_DIR}/lib  )
