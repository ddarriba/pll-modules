
file(GLOB MSA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/../pllmod_common.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_library(pllmodmsa OBJECT ${MSA_SOURCES})
target_include_directories(pllmodmsa PRIVATE ${PLLMOD_INCLUDE_PATH})
set_property(TARGET pllmodmsa PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(pllmodmsa_shared SHARED $<TARGET_OBJECTS:pllmodmsa>)
add_library(pllmodmsa_static STATIC $<TARGET_OBJECTS:pllmodmsa>)

target_link_libraries(pllmodmsa_shared pll_shared) 
