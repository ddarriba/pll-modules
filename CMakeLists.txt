cmake_minimum_required (VERSION 2.8.10 FATAL_ERROR)

project (Pll-modules)

if(NOT PLLMODULES_LIBPLL_PATH)
  set(PLLMODULES_LIBPLL_PATH "libs/libpll" CACHE STRING "" FORCE)
endif()

# if the project wants to build libpll itself...
if (NOT PLLMODULES_SKIP_LIBPLL_BUILD)
  add_subdirectory(${PLLMODULES_LIBPLL_PATH})
endif()

add_subdirectory(src)

set(PLLMOD_EXPORTED_INCLUDE_PATH ${PLLMOD_EXPORTED_INCLUDE_PATH} ${PLL_EXPORTED_INCLUDE_PATH} PARENT_SCOPE)

set(PLLMOD_EXPORTED_LIBS_STATIC
  pllmodalgorithm_static
  pllmodbinary_static
  pllmodmsa_static
  pllmodoptimize_static
  pllmodtree_static
  pllmodutil_static
  pll_static
  PARENT_SCOPE
)

set(PLLMOD_EXPORTED_LIBS_SHARED
  pllmodalgorithm_shared
  pllmodbinary_shared
  pllmodmsa_shared
  pllmodoptimize_shared
  pllmodtree_shared
  pllmodutil_shared
  pll_shared
  PARENT_SCOPE
)