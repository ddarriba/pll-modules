find_package(BISON)
find_package(FLEX)
set(TREE_BISON_FLAGS "-y -d -p pllmod_utree_")
set(TREE_FLEX_FLAGS "-P pllmod_utree_")

BISON_TARGET(split_utree_t 
  ${CMAKE_CURRENT_SOURCE_DIR}/split_utree.y ${CMAKE_CURRENT_BINARY_DIR}/split_utree.c
  COMPILE_FLAGS ${TREE_BISON_FLAGS})
FLEX_TARGET(lex_split_t 
  ${CMAKE_CURRENT_SOURCE_DIR}/lex_split.l ${CMAKE_CURRENT_BINARY_DIR}/lex_split.c
  COMPILE_FLAGS ${TREE_FLEX_FLAGS})
ADD_FLEX_BISON_DEPENDENCY(lex_split_t split_utree_t)

set(TREE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/consensus.c 
  ${CMAKE_CURRENT_SOURCE_DIR}/pll_tree.c
  ${CMAKE_CURRENT_SOURCE_DIR}/rtree_operations.c
  ${CMAKE_CURRENT_SOURCE_DIR}/tree_hashtable.c
  ${CMAKE_CURRENT_SOURCE_DIR}/treeinfo.c
  ${CMAKE_CURRENT_SOURCE_DIR}/utree_distances.c
  ${CMAKE_CURRENT_SOURCE_DIR}/utree_operations.c
  ${BISON_split_utree_t_OUTPUTS}
  ${FLEX_lex_split_t_OUTPUTS}
)

# remove flex warnings
SET_SOURCE_FILES_PROPERTIES( ${FLEX_split_utree_t_OUTPUTS} PROPERTIES COMPILE_FLAGS -Wno-sign-compare )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_library(pllmodtree OBJECT ${TREE_SOURCES})
target_include_directories(pllmodtree PRIVATE ${PLLMOD_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET pllmodtree PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(pllmodtree_shared SHARED $<TARGET_OBJECTS:pllmodtree>)
add_library(pllmodtree_static STATIC $<TARGET_OBJECTS:pllmodtree>)

target_link_libraries(pllmodtree_shared pll_shared) 
