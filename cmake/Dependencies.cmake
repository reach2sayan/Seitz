# Every dependency is fetched. Nothing is taken from the system, nothing is
# searched for, and there is no find_package fallback to go stale: a fresh
# checkout configures with only a compiler and CMake present, and every
# developer builds against byte-identical sources. The one exception is
# Threads, which is part of the toolchain rather than a library we vendor.
#
# Both Eigen and Boost are PUBLIC dependencies of the exported cppcrystal
# target, so both must also be installable — otherwise install(EXPORT) writes a
# config whose find_dependency() calls have nothing to point at. Each defaults
# to skipping its install rules when built as a subproject; turn them back on so
# `cmake --install` writes one self-contained prefix.
#
# SYSTEM on every FetchContent_Declare keeps third-party headers out of the
# library's -Wconversion/-Wsign-conversion warning set.
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

# Boost 1.88.0, from the CMake-ready release archive. The archive carries every
# Boost library; BOOST_INCLUDE_LIBRARIES decides which get CMake targets, so an
# unused entry is pure configure-time cost. What each one is for:
#   container  — static_vector/small_vector/flat_map/flat_set
#   flyweight  — one shared immutable SpaceGroup per Hall setting
#   geometry   — the R-tree behind PositionIndex
#   leaf       — the error model (Result<T>)
# Nothing here parses text at runtime: the tables are constexpr, generated
# offline by tools/transcribe_*.py. So boost::parser is not configured; were
# that to change, 1.87 is the floor for it.
set(BOOST_INCLUDE_LIBRARIES container flyweight geometry leaf)
FetchContent_Declare(Boost
        URL https://github.com/boostorg/boost/releases/download/boost-1.88.0/boost-1.88.0-cmake.tar.xz
        DOWNLOAD_EXTRACT_TIMESTAMP ON
        SYSTEM)

FetchContent_MakeAvailable(Eigen3 Boost)

# Catch2 v3. extras/ holds the Catch.cmake module providing
# catch_discover_tests(); FetchContent does not add it to CMAKE_MODULE_PATH.
if (CPPCRYSTAL_BUILD_TESTS)
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
#   CPPCRYSTAL_BUILD_ORACLE_TESTS and never linked into the cppcrystal library.
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
if (CPPCRYSTAL_BUILD_ORACLE_TESTS)
    set(SPGLIB_WITH_TESTS OFF CACHE BOOL "" FORCE)
    set(SPGLIB_WITH_Fortran OFF CACHE BOOL "" FORCE)
    set(SPGLIB_WITH_Python OFF CACHE BOOL "" FORCE)
    set(SPGLIB_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(SPGLIB_INSTALL OFF CACHE BOOL "" FORCE)
    set(CPPCRYSTAL_SPGLIB_SUBDIR "")
else ()
    set(CPPCRYSTAL_SPGLIB_SUBDIR SOURCE_SUBDIR cppcrystal-tables-only)
endif ()
FetchContent_Declare(spglib_reference
        GIT_REPOSITORY https://github.com/spglib/spglib.git
        GIT_TAG v2.7.0
        GIT_SHALLOW TRUE
        ${CPPCRYSTAL_SPGLIB_SUBDIR})
FetchContent_MakeAvailable(spglib_reference)
unset(CPPCRYSTAL_SPGLIB_SUBDIR)

# The transcribers take file paths, so a missing checkout has to fail here with
# a sentence rather than three directories later with a Python traceback.
if (NOT EXISTS "${spglib_reference_SOURCE_DIR}/src/spg_database.c")
    message(FATAL_ERROR
            "The reference spglib checkout is missing its src/ directory "
            "(looked in ${spglib_reference_SOURCE_DIR}). The data tables are "
            "transcribed from it at build time; see the codegen block in "
            "CMakeLists.txt.")
endif ()

# pybind11 3.x — the Python binding layer. Fetched here like everything else
# rather than taken from PyPI as a build requirement, so that
# `cmake -DCPPCRYSTAL_BUILD_PYTHON=ON` and `pip install .` compile against
# byte-identical headers; pyproject.toml therefore names only scikit-build-core
# in build-system.requires and says nothing about pybind11. 3.x specifically,
# for py::native_enum — it hands Python a real enum.IntEnum instead of a
# pybind11-private enum type, so the bound enums pickle and pattern-match like
# any other.
#
# PYBIND11_FINDPYTHON so pybind11 uses FindPython rather than the deprecated
# FindPythonInterp, which is what lets scikit-build-core's injected
# Python_EXECUTABLE hint select the interpreter the wheel is being built for.
# FindPython is a toolchain query, not a library we vendor — the same exception
# this file already makes for Threads.
#
# No network at configure time? FETCHCONTENT_SOURCE_DIR_PYBIND11=/path/to/checkout
# points this at a local clone, which matters now that a FetchContent configure
# can happen inside `pip install`.
if (CPPCRYSTAL_BUILD_PYTHON)
    set(PYBIND11_FINDPYTHON ON CACHE BOOL "" FORCE)
    set(PYBIND11_INSTALL OFF CACHE BOOL "" FORCE)
    set(PYBIND11_TEST OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG v3.0.1
            GIT_SHALLOW TRUE
            SYSTEM)
    FetchContent_MakeAvailable(pybind11)
endif ()
