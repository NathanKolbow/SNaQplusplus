cmake_minimum_required(VERSION 3.16)

project(CSNaQ
    VERSION 0.1.0
    DESCRIPTION "C++ port of SNaQ.jl phylogenetic network inference"
    LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(CSNAQ_BUILD_TESTS "Build the CSNaQ unit tests (fetches GoogleTest)" OFF)
option(CSNAQ_NATIVE "Optimize for the build machine's CPU (-march=native)" OFF)
option(CSNAQ_BUILD_TOOLS "Build the validation and benchmark tools" OFF)
option(CSNAQ_JULIA_TESTS "With CSNAQ_BUILD_TESTS: also run the tests that check snaq's output with SNaQ.jl (needs julia and the environment of tools/julia)" OFF)

include(FetchContent)

# --------------------------------------------------------------------------
# NLopt: fetched and built statically so users need nothing preinstalled.
# --------------------------------------------------------------------------
set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
set(NLOPT_MATLAB OFF CACHE BOOL "" FORCE)
set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
set(NLOPT_SWIG OFF CACHE BOOL "" FORCE)
set(NLOPT_FORTRAN OFF CACHE BOOL "" FORCE)
set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)
set(NLOPT_JAVA OFF CACHE BOOL "" FORCE)
set(NLOPT_CXX OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(nlopt
    GIT_REPOSITORY https://github.com/stevengj/nlopt.git
    GIT_TAG v2.10.0
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(nlopt)

find_package(Threads REQUIRED)

# --------------------------------------------------------------------------
# Core library
# --------------------------------------------------------------------------
add_library(csnaq_core STATIC
    src/network.cpp
    src/network_manip.cpp
    src/newick.cpp
    src/quartets.cpp
    src/equations.cpp
    src/likelihood.cpp
    src/optimize.cpp
    src/classes.cpp
    src/eqndag.cpp
    src/porteqn.cpp
    src/signature.cpp
    src/snaq_network.cpp
    src/compare.cpp
    src/moves.cpp
    src/restrictions.cpp
    src/identifiability.cpp
    src/search.cpp
    src/snaq.cpp
    src/util.cpp)
target_include_directories(csnaq_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(csnaq_core PUBLIC nlopt Threads::Threads)

if(MSVC)
    target_compile_options(csnaq_core PRIVATE /W4 /permissive- /utf-8)
    target_compile_definitions(csnaq_core PUBLIC _USE_MATH_DEFINES NOMINMAX)
else()
    target_compile_options(csnaq_core PRIVATE -Wall -Wextra -Wno-unused-parameter)
    if(CSNAQ_NATIVE)
        target_compile_options(csnaq_core PUBLIC -march=native)
    endif()
endif()

add_executable(snaq apps/snaq_main.cpp apps/cli.cpp)
target_link_libraries(snaq PRIVATE csnaq_core)
target_compile_definitions(snaq PRIVATE CSNAQ_VERSION="${PROJECT_VERSION}")
if(MSVC)
    target_compile_options(snaq PRIVATE /utf-8)
endif()

include(GNUInstallDirs)
install(TARGETS snaq RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(CSNAQ_BUILD_TOOLS)
# Validation tool: expected CFs, scores and gradients for networks (compared
# against SNaQ.jl by tools/julia scripts).
add_executable(csnaq_check tools/csnaq_check.cpp)
target_link_libraries(csnaq_check PRIVATE csnaq_core)

# Validation of the quartet-class likelihood against the per-quartet one.
add_executable(csnaq_classcheck tools/classcheck.cpp)
target_link_libraries(csnaq_classcheck PRIVATE csnaq_core)

# Benchmark of numerical parameter optimization (compared against SNaQ.jl by
# tools/julia/bench_fit.jl).
add_executable(csnaq_bench_fit tools/bench_fit.cpp)
target_link_libraries(csnaq_bench_fit PRIVATE csnaq_core)

# Accuracy evaluation used by the benchmark framework (benchmarks/bench.py).
add_executable(csnaq_eval tools/eval.cpp)
target_link_libraries(csnaq_eval PRIVATE csnaq_core)

# Prints tcgidentifiable(N) for each Newick line (compared against SNaQ.jl by
# tools/julia/tcgcheck.jl).
add_executable(csnaq_tcgcheck tools/tcgcheck.cpp)
target_link_libraries(csnaq_tcgcheck PRIVATE csnaq_core)
endif()

# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------
if(CSNAQ_BUILD_TESTS)
    enable_testing()
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.17.0
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
endif()
