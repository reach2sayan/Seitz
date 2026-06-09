cmake_minimum_required(VERSION 3.28)
project(CppCrystal VERSION 0.1.0 LANGUAGES CXX
        DESCRIPTION "Modern C++ port of spglib's C core")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif ()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(SPGLIB_BUILD_TESTS "Build the test suite" ON)
option(SPGLIB_BUILD_DEMO "Build the demo executable" ON)
option(SPGLIB_USE_MKL "Use Intel MKL as Eigen's BLAS/LAPACK backend" OFF)
option(SPGLIB_BUILD_ORACLE_TESTS "Build oracle tests against the reference spglib" OFF)
option(SPGLIB_ENABLE_SANITIZERS "Enable Address/UB sanitizers" OFF)

include(Dependencies)

find_program(CCACHE_PROGRAM ccache)
if (CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif ()

# ---- library ----
add_library(cppcrystal)
add_library(cppcrystal::cppcrystal ALIAS cppcrystal)
target_sources(cppcrystal PRIVATE
        src/version.cpp
        src/core/overlap.cpp
        src/reduce/delaunay.cpp
        src/reduce/niggli.cpp
        src/symmetry/find_symmetry.cpp
        src/symmetry/pointgroup.cpp
        src/symmetry/primitive.cpp
        src/data/spg_database.cpp
        src/spacegroup/hall_symbol.cpp)
target_include_directories(cppcrystal PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_compile_features(cppcrystal PUBLIC cxx_std_23)
target_link_libraries(cppcrystal PUBLIC Eigen3::Eigen Boost::headers)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(cppcrystal PRIVATE -Wall -Wextra -Wpedantic)
endif ()

if (SPGLIB_USE_MKL)
    target_compile_definitions(cppcrystal PUBLIC EIGEN_USE_MKL_ALL)
    target_link_libraries(cppcrystal PUBLIC MKL::MKL)
endif ()

if (SPGLIB_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(cppcrystal PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(cppcrystal PUBLIC -fsanitize=address,undefined)
endif ()

# ---- demo ----
if (SPGLIB_BUILD_DEMO)
    add_executable(cppcrystal_demo main.cpp)
    target_link_libraries(cppcrystal_demo PRIVATE cppcrystal::cppcrystal)
endif ()

# ---- tests ----
if (SPGLIB_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()
