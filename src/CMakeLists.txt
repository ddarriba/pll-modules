cmake_minimum_required (VERSION 2.8.10 FATAL_ERROR)

set(PLLMOD_CFLAGS  "-Wall -Wsign-compare -D_GNU_SOURCE -std=c99 -O3")
set(CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-undefined -Wl,dynamic_lookup")
set(PLLMOD_INCLUDE_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/../${PLLMODULES_LIBPLL_PATH}/src/)

add_subdirectory(algorithm)
add_subdirectory(binary)
add_subdirectory(msa)
add_subdirectory(tree)
add_subdirectory(util)
add_subdirectory(optimize)

set(PLLMOD_EXPORTED_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/algorithm/" "${CMAKE_CURRENT_SOURCE_DIR}/binary" "${CMAKE_CURRENT_SOURCE_DIR}/msa" "${CMAKE_CURRENT_SOURCE_DIR}/optimize"  "${CMAKE_CURRENT_SOURCE_DIR}/tree" "${CMAKE_CURRENT_SOURCE_DIR}/util" ${LIBPLL_SRC} PARENT_SCOPE)

