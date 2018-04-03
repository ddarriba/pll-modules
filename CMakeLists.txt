cmake_minimum_required (VERSION 2.8.10 FATAL_ERROR)

project (Pll-modules)

if(NOT PLLMODULES_LIBPLL_PATH)
  set(PLLMODULES_LIBPLL_PATH "libs/libpll" CACHE STRING "" FORCE) 
endif()

add_subdirectory(${PLLMODULES_LIBPLL_PATH})
add_subdirectory(src)

set(PLLMOD_EXPORTED_INCLUDE_PATH ${PLLMOD_EXPORTED_INCLUDE_PATH} PARENT_SCOPE)

