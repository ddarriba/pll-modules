cmake_minimum_required (VERSION 2.8.10 FATAL_ERROR)

set(PLLMOD_CFLAGS  "-Wall -Wsign-compare -D_GNU_SOURCE -std=c99 -O3")
set(PLLMOD_INCLUDE_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/../libs/libpll/src/)

add_subdirectory(algorithm)
add_subdirectory(binary)
add_subdirectory(msa)
add_subdirectory(tree)
add_subdirectory(util)
add_subdirectory(optimize)


