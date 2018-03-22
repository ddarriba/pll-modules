cmake_minimum_required (VERSION 2.8.10 FATAL_ERROR)

project (Pll-modules)

add_subdirectory(libs/libpll)

add_subdirectory(src)

set(PLLMOD_EXPORTED_INCLUDE_PATH ${PLLMOD_EXPORTED_INCLUDE_PATH} PARENT_SCOPE)

