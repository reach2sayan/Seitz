# Every dependency is fetched. Nothing is taken from the system
# Both Eigen and Boost are PUBLIC dependencies of the exported seitz
# target, so both must also be installable — otherwise install(EXPORT) writes a
# config whose find_dependency() calls have nothing to point at
include(FetchContent)

set(EIGEN_BUILD_CMAKE_PACKAGE ON CACHE BOOL "" FORCE)
set(BOOST_SKIP_INSTALL_RULES OFF CACHE BOOL "" FORCE)

# Eigen 5.0.0 — the first release in which fixed-size Matrix/Array are literal
# types, so they can be constructed and accessed in constexpr context (we
# compile as C++23). As a subproject Eigen's tests/blas/lapack/docs/demos all
# default OFF, leaving only the Eigen3::Eigen interface target.
FetchContent_Declare(Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.0
        GIT_SHALLOW TRUE
        SYSTEM)

# Boost 1.88.0, from the CMake-ready release archive.
#   container  — static_vector/small_vector/flat_map/flat_set
#   flyweight  — one shared immutable SpaceGroup per Hall setting
#   geometry   — the R-tree behind PositionIndex
#   graph      — the algorithms over the constexpr subgroup graph (BFS, views)
#   leaf       — the error model (Result<T>)
#   parser     — the CIF and xyz-triplet grammars; header-only, numerics via
#                <charconv> (no compiled lib)
#   algorithm  — to_lower_copy/join, so the CIF layer hand-rolls neither
#   range      — join(), for concatenating two ranges
set(BOOST_INCLUDE_LIBRARIES algorithm container flyweight geometry graph leaf
        parser range)
FetchContent_Declare(Boost
        URL https://github.com/boostorg/boost/releases/download/boost-1.88.0/boost-1.88.0-cmake.tar.xz
        DOWNLOAD_EXTRACT_TIMESTAMP ON
        SYSTEM)

FetchContent_MakeAvailable(Eigen3 Boost)

if (TARGET boost_container)
    target_compile_options(boost_container PRIVATE
            $<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:GNU>>:-Wno-aggressive-loop-optimizations>)
endif ()

# Catch2 v3. extras/ holds the Catch.cmake module providing
# catch_discover_tests();
if (SEITZ_BUILD_TESTS)
    FetchContent_Declare(Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG v3.8.1
            GIT_SHALLOW TRUE
            SYSTEM)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
endif ()

# Reference spglib v2.7.0, in two roles that share one pinned tag so they can
# never drift apart:
#
#   Source of truth for the transcribed data tables. Four of the generators in
#   tools/ read spg_database.c, msg_database.c, hall_symbol.c and
#   sitesym_database.c straight out of this checkout — see the codegen block in
#   CMakeLists.txt — so the tables the library compiles and the reference
#   implementation the oracle compares them against are the same spglib.
#
#   The validation oracle itself, exposing Spglib::symspg. Built only under
#   SEITZ_BUILD_ORACLE_TESTS and never linked into the seitz library.
#
# Which means the checkout is unconditional but the BUILD is not. SOURCE_SUBDIR
# naming a directory with no CMakeLists.txt is the documented way to say
# "populate, do not add_subdirectory": a default build pays for the download and
# nothing else, rather than configuring and compiling a second crystallography
# library it will not link.
#
# No network at configure time? FETCHCONTENT_SOURCE_DIR_SPGLIB_REFERENCE=/path/to/checkout
# points this at a local clone — the same escape hatch pybind11 gets below, and
# it matters more here now that every build needs these sources, not just the
# oracle ones.
if (SEITZ_BUILD_ORACLE_TESTS)
    set(SPGLIB_WITH_TESTS OFF CACHE BOOL "" FORCE)
    set(SPGLIB_WITH_Fortran OFF CACHE BOOL "" FORCE)
    set(SPGLIB_WITH_Python OFF CACHE BOOL "" FORCE)
    set(SPGLIB_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(SPGLIB_INSTALL OFF CACHE BOOL "" FORCE)
    set(SEITZ_SPGLIB_SUBDIR "")
else ()
    set(SEITZ_SPGLIB_SUBDIR SOURCE_SUBDIR seitz-tables-only)
endif ()
FetchContent_Declare(spglib_reference
        GIT_REPOSITORY https://github.com/spglib/spglib.git
        GIT_TAG v2.7.0
        GIT_SHALLOW TRUE
        ${SEITZ_SPGLIB_SUBDIR})
FetchContent_MakeAvailable(spglib_reference)
unset(SEITZ_SPGLIB_SUBDIR)

# PyXtal, for its two CIF corpora only (209 per-group files, 77 real
# structures). Populate-only: SOURCE_SUBDIR names a directory that does not
# exist, so FetchContent checks the tree out and adds no targets. Offline:
# FETCHCONTENT_SOURCE_DIR_PYXTAL_REFERENCE.
if (SEITZ_BUILD_ORACLE_TESTS)
    FetchContent_Declare(pyxtal_reference
            GIT_REPOSITORY https://github.com/MaterSim/PyXtal.git
            GIT_TAG v1.1.4
            GIT_SHALLOW TRUE
            SOURCE_SUBDIR seitz-cifs-only)
    FetchContent_MakeAvailable(pyxtal_reference)

    foreach (dir IN ITEMS miscellaneous database)
        if (NOT EXISTS "${pyxtal_reference_SOURCE_DIR}/pyxtal/${dir}/cifs")
            message(FATAL_ERROR
                    "The PyXtal checkout is missing pyxtal/${dir}/cifs (looked "
                    "in ${pyxtal_reference_SOURCE_DIR}). The CIF corpus tests "
                    "read both directories; see tests/cif_corpus.hpp.")
        endif ()
    endforeach ()
endif ()

if (NOT EXISTS "${spglib_reference_SOURCE_DIR}/src/spg_database.c")
    message(FATAL_ERROR
            "The reference spglib checkout is missing its src/ directory "
            "(looked in ${spglib_reference_SOURCE_DIR}). The data tables are "
            "transcribed from it at build time; see the codegen block in "
            "CMakeLists.txt.")
endif ()

# pybind11 3.x — the Python binding layer.
if (SEITZ_BUILD_PYTHON)
    set(PYBIND11_FINDPYTHON ON CACHE BOOL "" FORCE)
    set(PYBIND11_INSTALL OFF CACHE BOOL "" FORCE)
    set(PYBIND11_TEST OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG v3.0.4
            GIT_SHALLOW TRUE
            SYSTEM)
    FetchContent_MakeAvailable(pybind11)
endif ()
