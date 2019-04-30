add_executable(raxml-dks 
    main.cpp
    msa.cpp
    test_case.cpp
    tree.cpp
    partition.cpp
    model.cpp
    benchmark.cpp
)
set_target_properties(raxml-dks PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
)
set(${CMAKE_BUILD_TYPE} DEBUG)
target_include_directories(raxml-dks PRIVATE ${CMAKE_SOURCE_DIR}/lib/libpll/src/)
target_link_libraries(raxml-dks pll_static)
target_compile_options(raxml-dks PRIVATE -Wall -Wextra -pedantic)
