
file(GLOB UTIL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/../pllmod_common.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_library(pllmodutil OBJECT ${UTIL_SOURCES})
target_include_directories(pllmodutil PRIVATE ${PLLMOD_INCLUDE_PATH})
set_property(TARGET pllmodutil PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(pllmodutil_shared SHARED $<TARGET_OBJECTS:pllmodutil>)
add_library(pllmodutil_static STATIC $<TARGET_OBJECTS:pllmodutil>)

target_link_libraries(pllmodutil_shared pll_shared) 
