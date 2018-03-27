find_package(BISON)
find_package(FLEX)
set(LIBPLL_BISON_FLAGS "-y -d -p pll_utree_")
set(LIBPLL_FLEX_FLAGS "-P pll_utree_")

BISON_TARGET(split_utree_t 
  ${CMAKE_CURRENT_SOURCE_DIR}/split_utree.y ${CMAKE_CURRENT_BINARY_DIR}/split_utree.c
  COMPILE_FLAGS ${LIBPLL_BISON_FLAGS})
FLEX_TARGET(lex_split_t 
  ${CMAKE_CURRENT_SOURCE_DIR}/lex_split.l ${CMAKE_CURRENT_BINARY_DIR}/lex_split.c
  COMPILE_FLAGS ${LIBPLL_FLEX_FLAGS})
ADD_FLEX_BISON_DEPENDENCY(lex_split_t split_utree_t)

file(GLOB TREE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../pllmod_common.c
  ${BISON_split_utree_t_OUTPUTS}
  ${FLEX_lex_split_t_OUTPUTS}
  )


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_library(pllmodtree OBJECT ${TREE_SOURCES})
target_include_directories(pllmodtree PRIVATE ${PLLMOD_INCLUDE_PATH})
set_property(TARGET pllmodtree PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(pllmodtree_shared SHARED $<TARGET_OBJECTS:pllmodtree>)
add_library(pllmodtree_static STATIC $<TARGET_OBJECTS:pllmodtree>)

