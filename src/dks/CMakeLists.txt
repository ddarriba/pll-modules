file(GLOB DKS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_pllmodules_lib(pllmoddks "${DKS_SOURCES}")
target_include_directories(pllmoddks_obj PRIVATE ${PLLMOD_INCLUD_PATH})
