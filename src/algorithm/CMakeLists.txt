
file(GLOB ALGORITHM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/../pllmod_common.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_library(pllmodalgorithm OBJECT ${ALGORITHM_SOURCES})
target_include_directories(pllmodalgorithm PRIVATE ${PLLMOD_INCLUDE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/../binary
  ${CMAKE_CURRENT_SOURCE_DIR}/../msa
  ${CMAKE_CURRENT_SOURCE_DIR}/../optimize
  ${CMAKE_CURRENT_SOURCE_DIR}/../tree
  ${CMAKE_CURRENT_SOURCE_DIR}/../util
  )
set_property(TARGET pllmodalgorithm PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(pllmodalgorithm_shared SHARED $<TARGET_OBJECTS:pllmodalgorithm>)
add_library(pllmodalgorithm_static STATIC $<TARGET_OBJECTS:pllmodalgorithm>)

target_link_libraries(pllmodalgorithm_shared pll_shared) 
