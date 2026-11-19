add_executable(snaqpp_tests
    test_support.cpp
    test_cf_recursive_blocks.cpp
    test_opt_api.cpp
    test_inplace_and_data.cpp
    test_network_moves.cpp
    test_search_api.cpp
    test_propquartets_lazy_bugs.cpp
    test_newick_and_cli.cpp
    test_classes.cpp
    test_search_shortcuts.cpp
    test_readcsv.cpp)
target_link_libraries(snaqpp_tests PRIVATE snaqpp_core GTest::gtest_main)
target_include_directories(snaqpp_tests PRIVATE ${CMAKE_SOURCE_DIR}/apps ${CMAKE_SOURCE_DIR}/tools)
target_sources(snaqpp_tests PRIVATE ${CMAKE_SOURCE_DIR}/apps/cli.cpp)
target_compile_definitions(snaqpp_tests PRIVATE SNAQPP_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
if(MSVC)
    target_compile_options(snaqpp_tests PRIVATE /utf-8 /bigobj)
endif()

include(GoogleTest)
gtest_discover_tests(snaqpp_tests DISCOVERY_TIMEOUT 120 PROPERTIES TIMEOUT 3600)

# The restriction flags of the snaq program, checked on its output files with
# SNaQ.jl's restriction functions (tools/julia/check_restrictions.jl).
if(SNAQPP_JULIA_TESTS)
    find_program(JULIA_EXECUTABLE julia REQUIRED)
    add_test(NAME RestrictionFlags.OutputMeetsRestrictionsOfSNaQjl
             COMMAND ${JULIA_EXECUTABLE} --project=${CMAKE_SOURCE_DIR}/tools/julia
                     ${CMAKE_SOURCE_DIR}/tools/julia/check_restrictions.jl $<TARGET_FILE:snaq>
                     ${CMAKE_SOURCE_DIR}/benchmarks/datasets/n20_h2)
    set_tests_properties(RestrictionFlags.OutputMeetsRestrictionsOfSNaQjl PROPERTIES TIMEOUT 3600)
endif()
