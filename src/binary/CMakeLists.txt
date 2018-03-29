
file(GLOB BINARY_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/../pllmod_common.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_library(pllmodbinary OBJECT ${BINARY_SOURCES})
target_include_directories(pllmodbinary PRIVATE ${PLLMOD_INCLUDE_PATH})
set_property(TARGET pllmodbinary PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(pllmodbinary_shared SHARED $<TARGET_OBJECTS:pllmodbinary>)
add_library(pllmodbinary_static STATIC $<TARGET_OBJECTS:pllmodbinary>)

target_link_libraries(pllmodbinary_shared pll_shared) 
