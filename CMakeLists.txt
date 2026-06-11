cmake_minimum_required(VERSION 3.28)
project(CppCrystal VERSION 0.1.0 LANGUAGES CXX
        DESCRIPTION "Modern C++ port of spglib's C core")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
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
        src/dataset.cpp
        src/analysis/symmetry_analyzer.cpp
        src/analysis/magnetic_symmetry_analyzer.cpp
        src/group/space_group.cpp
        src/group/wyckoff_position.cpp
        src/group/subgroup_graph.cpp
        src/generate/wyckoff_combinations.cpp
        src/generate/random_lattice.cpp
        src/generate/distance_check.cpp
        src/generate/crystal_builder.cpp
        src/core/overlap.cpp
        src/core/centering.cpp
        src/reduce/delaunay.cpp
        src/reduce/niggli.cpp
        src/symmetry/find_symmetry.cpp
        src/symmetry/pointgroup.cpp
        src/symmetry/primitive.cpp
        src/data/spg_database.cpp
        src/data/sitesym_database.cpp
        src/data/msg_database.cpp
        src/data/element_data.cpp
        src/spin/spin.cpp
        src/magnetic/magnetic_spacegroup.cpp
        src/magnetic_dataset.cpp
        src/kpoint/grid.cpp
        src/kpoint/reciprocal_mesh.cpp
        src/kpoint/brillouin_zone.cpp
        src/spacegroup/hall_symbol.cpp
        src/spacegroup/spacegroup.cpp
        src/refine/refinement.cpp
        src/refine/site_symmetry.cpp
        src/refine/standardize.cpp
        src/refine/standardize_cell.cpp
        src/refine/operations.cpp)
target_include_directories(cppcrystal PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_compile_features(cppcrystal PUBLIC cxx_std_23)
target_link_libraries(cppcrystal PUBLIC Eigen3::Eigen Boost::headers)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Hardening warning set (Phase 9). The conversion warnings are applied to the
    # library only — the oracle tests bridge to the reference C API
    # (reinterpret_cast, size_t<->int), which would drown them in noise.
    target_compile_options(cppcrystal PRIVATE
            -Wall -Wextra -Wpedantic
            -Wshadow -Wcast-align -Wdouble-promotion -Wnull-dereference
            -Wimplicit-fallthrough -Wconversion -Wsign-conversion)
endif ()

# The magnetic database decodes all 76683 operations + the alternative
# transformations into constexpr tables at compile time; that exceeds Clang's
# default -fconstexpr-steps budget (GCC's default is far higher and needs no
# flag). Raise it for this one translation unit only.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set_source_files_properties(src/data/msg_database.cpp
            PROPERTIES COMPILE_OPTIONS "-fconstexpr-steps=64000000")
endif ()

if (SPGLIB_USE_MKL)
    target_compile_definitions(cppcrystal PUBLIC EIGEN_USE_MKL_ALL)
    target_link_libraries(cppcrystal PUBLIC MKL::MKL)
endif ()

if (SPGLIB_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(cppcrystal PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(cppcrystal PUBLIC -fsanitize=address,undefined)
endif ()

if (SPGLIB_BUILD_DEMO)
    add_executable(cppcrystal_demo main.cpp)
    target_link_libraries(cppcrystal_demo PRIVATE cppcrystal::cppcrystal)
endif ()

# Offline generator for the baked t-subgroup table (data/subgroup_relations.hpp).
# Not built by default — it re-runs the slow in-house derivation. Build with
# -DSPGLIB_BUILD_TOOLS=ON, run it, and redirect stdout over the header.
option(SPGLIB_BUILD_TOOLS "Build offline data generators" OFF)
if (SPGLIB_BUILD_TOOLS)
    add_executable(generate_subgroup_relations tools/generate_subgroup_relations.cpp)
    target_link_libraries(generate_subgroup_relations PRIVATE cppcrystal::cppcrystal)
endif ()

if (SPGLIB_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()
