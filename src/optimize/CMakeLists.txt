
file(GLOB OPTIMIZE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c 
  ${CMAKE_CURRENT_SOURCE_DIR}/lbfgsb/*.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../pllmod_common.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_library(pllmodoptimize OBJECT ${OPTIMIZE_SOURCES})
target_include_directories(pllmodoptimize PRIVATE ${PLLMOD_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../tree)
set_property(TARGET pllmodoptimize PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(pllmodoptimize_shared SHARED $<TARGET_OBJECTS:pllmodoptimize>)
add_library(pllmodoptimize_static STATIC $<TARGET_OBJECTS:pllmodoptimize>)

target_link_libraries(pllmodoptimize_shared pll_shared) 
