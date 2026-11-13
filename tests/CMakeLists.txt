add_executable(csnaq_tests
    test_support.cpp
    test_cf_recursive_blocks.cpp
    test_opt_api.cpp
    test_inplace_and_data.cpp
    test_network_moves.cpp
    test_search_api.cpp
    test_propquartets_lazy_bugs.cpp
    test_newick_and_cli.cpp)
target_link_libraries(csnaq_tests PRIVATE csnaq_core GTest::gtest_main)
target_include_directories(csnaq_tests PRIVATE ${CMAKE_SOURCE_DIR}/apps)
target_sources(csnaq_tests PRIVATE ${CMAKE_SOURCE_DIR}/apps/cli.cpp)
target_compile_definitions(csnaq_tests PRIVATE CSNAQ_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
if(MSVC)
    target_compile_options(csnaq_tests PRIVATE /utf-8 /bigobj)
endif()

include(GoogleTest)
gtest_discover_tests(csnaq_tests DISCOVERY_TIMEOUT 120 PROPERTIES TIMEOUT 3600)
