
file(GLOB ALGORITHM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/../pllmod_common.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_pllmodules_lib(pllmodalgorithm "${ALGORITHM_SOURCES}")
target_include_directories(pllmodalgorithm_obj PRIVATE ${PLLMOD_INCLUDE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/../binary
  ${CMAKE_CURRENT_SOURCE_DIR}/../msa
  ${CMAKE_CURRENT_SOURCE_DIR}/../optimize
  ${CMAKE_CURRENT_SOURCE_DIR}/../tree
  ${CMAKE_CURRENT_SOURCE_DIR}/../util
  )





