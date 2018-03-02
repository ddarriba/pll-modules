
file(GLOB TREE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c
  ${CMAKE_CURRENT_SOURCE_DIR}/*.y
  ${CMAKE_CURRENT_SOURCE_DIR}/*.l
  ${CMAKE_CURRENT_SOURCE_DIR}/../pllmod_common.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_library(pllmodtree OBJECT ${TREE_SOURCES})
target_include_directories(pllmodtree PRIVATE ${PLLMOD_INCLUDE_PATH})
set_property(TARGET pllmodtree PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(pllmodtree_shared SHARED $<TARGET_OBJECTS:pllmodtree>)
add_library(pllmodtree_static STATIC $<TARGET_OBJECTS:pllmodtree>)

