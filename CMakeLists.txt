cmake_minimum_required(VERSION 3.28)
project(Seitz VERSION 0.1.0 LANGUAGES CXX
        DESCRIPTION "Modern C++23 crystallography and symmetry-analysis library")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Nothing here is a C++20 module -- no .cppm, no `export module`, no import.
# CMake nonetheless turns on module scanning by default for C++20-and-later
# targets, which runs a separate dependency-scanning pass over every translation
# unit and generates a dyndep file, all to discover that there are no modules.
# Say so once instead. (Set before any target is defined, so it applies to the
# library, the demo, the tests, the tools and the Python extension alike.)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif ()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Version floors, not vendor gates. Each is here because the failure without it
# is a wall of template errors inside one of our headers rather than "your
# toolchain is too old", and each is reachable: a path that does not go through
# CMakePresets (pip install ., an IDE default, plain `cmake ..`) picks up
# whatever compiler happens to be first.
#
# GCC 14 and earlier accept -std=gnu++23 and then fail deep inside
# operation_set.hpp on std::from_range. MSVC before 19.43 (Visual Studio 17.13)
# has no <generator> and fails just as deep inside generate/assignments.hpp.
# Clang gets no floor: its constraint is the libstdc++ it is paired with, which
# is not what CMAKE_CXX_COMPILER_VERSION reports.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15)
    message(FATAL_ERROR
            "Seitz needs GCC 15+ for <generator> and std::from_range "
            "(got ${CMAKE_CXX_COMPILER_VERSION} at ${CMAKE_CXX_COMPILER}).\n"
            "Select one explicitly, e.g. CXX=g++-15 for a pip/uv build, or "
            "-DCMAKE_CXX_COMPILER=g++-15 for a direct cmake invocation.")
elseif (MSVC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.43)
    message(FATAL_ERROR
            "Seitz needs MSVC 19.43 / Visual Studio 17.13+ for <generator> "
            "(got ${CMAKE_CXX_COMPILER_VERSION}).")
endif ()

option(SEITZ_BUILD_TESTS "Build the test suite" ON)
option(SEITZ_BUILD_DEMO "Build the demo executable" ON)
option(SEITZ_BUILD_ORACLE_TESTS "Build oracle tests against the reference spglib" OFF)
option(SEITZ_BUILD_TOOLS "Build offline data generators" OFF)
option(SEITZ_BUILD_PYTHON "Build the pybind11 extension module" OFF)
option(SEITZ_ENABLE_SANITIZERS "Enable Address/UB sanitizers" OFF)
option(SEITZ_NATIVE_ARCH "Build for this machine's ISA (-march=native)" OFF)

include(Dependencies)

# ---- generated data tables ----
# The symmetry tables are transcriptions, not authored source: each is a
# mechanical rewrite of a spglib v2.7.0 C array (or, for the elements, of
# tools/element_data.csv) into a constexpr std::array. So they are generated
# into the build tree at build time rather than committed. A checked-in copy is
# a second source of truth that drifts the moment anyone hand-edits it, and it
# had already drifted: the `#pragma GCC visibility push(default)` window that
# the three public tables need was added to the committed files by hand and
# never to the generators, so a regeneration would silently have unexported
# them. The generators emit it now, which is what makes this move lossless.
#
# The spglib inputs come from the same pinned checkout that backs the oracle
# (cmake/Dependencies.cmake), so the tables the library compiles and the
# reference implementation the oracle diffs them against cannot be two
# different spglib versions.
#
# Two of the tables have a source of truth a C++ build cannot reach: PyXtal
# (rod groups) and this library itself (the t-subgroup graph, derived by
# tools/generate_subgroup_relations.cpp, which links the library that includes
# the table). Both are split the same way as the elements: the extraction is
# offline and writes a CSV under tools/, which is the tracked input, and the
# build transcribes that CSV like any other. So no generated header is tracked,
# and no build needs PyXtal, a two-stage bootstrap, or the slow derivation.
#
# Still not generated here: python/seitz/_core/*.pyi describe the built
# extension module, so they are downstream of the build, not upstream of it. CI
# already fails on a stale one.
#
# COMPONENTS Interpreter only: the transcribers are stdlib-only scripts, no
# third-party imports, so no virtualenv is needed to build C++. Spelled
# `Python` rather than `Python3` so it resolves the same cache variable
# pybind11 and scikit-build-core populate (see PYBIND11_FINDPYTHON).
find_package(Python 3.9 REQUIRED COMPONENTS Interpreter)

set(SEITZ_GENERATED_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/include")
set(SEITZ_GENERATED_SRC_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/src")

# At configure time, not from inside the custom commands: two of the five
# generators write into each of these directories, and having both race to
# create the same path under `ninja -j` is a failure mode worth not having.
file(MAKE_DIRECTORY
        "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/core"
        "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/data"
        "${SEITZ_GENERATED_SRC_DIR}/data")

# The one place the version is written is `project(Seitz VERSION ...)` above.
# version.hpp is configured from it rather than authored, so kVersion cannot
# drift from the wheel metadata (which scikit-build-core reads out of that same
# line) or from a release tag (which release.yml checks against it). Unlike the
# transcribed tables this is a configure_file, not a custom command: it depends
# on nothing but a CMake variable, so there is nothing to rebuild.
configure_file(cmake/version.hpp.in
        "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/core/version.hpp"
        @ONLY)

# One custom command per generator. Every transcriber takes its inputs then its
# outputs as positional arguments, which is what lets this be one function
# rather than five near-identical add_custom_command() blocks. The script is a
# DEPENDS as well as the data, so editing a transcriber regenerates its tables
# and recompiles whatever includes them.
function(seitz_transcribe script)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "" "INPUTS;OUTPUTS")
    add_custom_command(
            OUTPUT ${arg_OUTPUTS}
            COMMAND ${Python_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/tools/${script}"
            ${arg_INPUTS} ${arg_OUTPUTS}
            DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tools/${script}" ${arg_INPUTS}
            COMMENT "Transcribing data tables with ${script}"
            VERBATIM)
endfunction()

set(SEITZ_SPGLIB_SOURCE "${spglib_reference_SOURCE_DIR}/src")

# Public — installed, and reachable from <seitz/seitz.hpp>.
set(SEITZ_GENERATED_PUBLIC_HEADERS
        "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/core/version.hpp"
        "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/data/element_tables.hpp"
        "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/data/spacegroup_metadata_tables.hpp"
        "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/data/magnetic_spacegroup_metadata_tables.hpp"
        "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/data/subgroup_relations.hpp")
# Private — the encoded-operation packing, decoded in the .cpp files. Never
# installed, and a public header including one of these is a bug.
set(SEITZ_GENERATED_PRIVATE_HEADERS
        "${SEITZ_GENERATED_SRC_DIR}/data/spacegroup_operation_tables.hpp"
        "${SEITZ_GENERATED_SRC_DIR}/data/magnetic_spacegroup_operation_tables.hpp"
        "${SEITZ_GENERATED_SRC_DIR}/data/hall_generators.hpp"
        "${SEITZ_GENERATED_SRC_DIR}/data/sitesym_tables.hpp"
        "${SEITZ_GENERATED_SRC_DIR}/data/rod_group_tables.hpp")

seitz_transcribe(transcribe_elements.py
        INPUTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/element_data.csv"
        OUTPUTS "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/data/element_tables.hpp")

seitz_transcribe(transcribe_spg_database.py
        INPUTS "${SEITZ_SPGLIB_SOURCE}/spg_database.c"
        OUTPUTS "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/data/spacegroup_metadata_tables.hpp"
        "${SEITZ_GENERATED_SRC_DIR}/data/spacegroup_operation_tables.hpp")

seitz_transcribe(transcribe_msg_database.py
        INPUTS "${SEITZ_SPGLIB_SOURCE}/msg_database.c"
        OUTPUTS "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/data/magnetic_spacegroup_metadata_tables.hpp"
        "${SEITZ_GENERATED_SRC_DIR}/data/magnetic_spacegroup_operation_tables.hpp")

seitz_transcribe(transcribe_hall_generators.py
        INPUTS "${SEITZ_SPGLIB_SOURCE}/hall_symbol.c"
        OUTPUTS "${SEITZ_GENERATED_SRC_DIR}/data/hall_generators.hpp")

seitz_transcribe(transcribe_sitesym_database.py
        INPUTS "${SEITZ_SPGLIB_SOURCE}/sitesym_database.c"
        OUTPUTS "${SEITZ_GENERATED_SRC_DIR}/data/sitesym_tables.hpp")

# The two whose CSV inputs are extracted offline: tools/extract_rod_groups.py
# (PyXtal) and tools/generate_subgroup_relations.cpp (this library).
seitz_transcribe(transcribe_rod_groups.py
        INPUTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/rod_group_data.csv"
        OUTPUTS "${SEITZ_GENERATED_SRC_DIR}/data/rod_group_tables.hpp")

seitz_transcribe(transcribe_subgroup_relations.py
        INPUTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/subgroup_relations.csv"
        OUTPUTS "${SEITZ_GENERATED_INCLUDE_DIR}/seitz/data/subgroup_relations.hpp")

# NOT MSVC: ccache's cl.exe support is narrow enough to be a liability (it
# needs /Z7 rather than /Zi and misses on several common flag shapes). The
# Windows CI job passes sccache explicitly instead.
find_program(CCACHE_PROGRAM ccache)
if (CCACHE_PROGRAM AND NOT MSVC)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif ()

# ---- library ----
# The type is decided here, not via BUILD_SHARED_LIBS: it is part of the
# contract, so a consumer's cache variable cannot silently change it.
#
# SHARED where a shared object is the deliverable -- the standalone library a
# consumer links against. Being shared is what makes the exported dependency set
# below exact: PRIVATE links (Boost::flyweight, Boost::geometry) stay out of the
# interface, which a static archive cannot promise (see
# cmake/SeitzConfig.cmake.in, which branches on exactly that).
#
# STATIC in the two cases where it is not, for reasons that are not about
# packaging taste:
#
#   The Python extension absorbs the library, so a wheel is one binary rather
#   than a .pyd/.so plus a sibling it has to locate at import time.
#
#   Windows, because of Boost.LEAF. PE has no symbol interposition, and LEAF
#   keeps its active-context slots in thread_local static members of a class
#   template (boost/leaf/config/tls_cpp11.hpp; BOOST_LEAF_SYMBOL_VISIBLE is
#   empty off GCC). A DLL and its consumer would each get a private copy, so an
#   error raised inside the library would reach the caller's try_handle_all with
#   its payload silently dropped -- which is what tests/test_error_paths.cpp and
#   python/tests/test_errors.py exist to catch. WINDOWS_EXPORT_ALL_SYMBOLS
#   cannot help: the consumer's own COMDAT definition wins over any import.
#   Linking one module is the only fix that does not fork Boost. It also makes
#   the `#pragma GCC visibility` windows in the public headers and
#   SEITZ_TESTABLE genuinely inert rather than something to work around.
if (SEITZ_BUILD_PYTHON OR WIN32)
    set(SEITZ_LIBRARY_TYPE STATIC)
else ()
    set(SEITZ_LIBRARY_TYPE SHARED)
endif ()
add_library(seitz ${SEITZ_LIBRARY_TYPE})
add_library(seitz::seitz ALIAS seitz)
# PIC unconditionally: as a static archive this library still gets linked into a
# shared object (the Python extension), which needs its objects relocatable.
set_target_properties(seitz PROPERTIES POSITION_INDEPENDENT_CODE ON)
if (SEITZ_LIBRARY_TYPE STREQUAL "SHARED")
    set_target_properties(seitz PROPERTIES
            # 0.x: SOVERSION tracks the minor, since every minor may break ABI
            # before 1.0. From 1.0 on this becomes ${PROJECT_VERSION_MAJOR}.
            VERSION ${PROJECT_VERSION}
            SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
endif ()
target_sources(seitz PRIVATE
        src/version.cpp
        src/warmup.cpp
        src/analysis/symmetry_analyzer.cpp
        src/analysis/magnetic_symmetry_analyzer.cpp
        src/group/space_group.cpp
        src/group/wyckoff.cpp
        src/group/point_group.cpp
        src/group/subgroup_graph.cpp
        src/generate/random_lattice.cpp
        src/generate/generator.cpp
        src/generate/distance_check.cpp
        # Alloy configurational thermodynamics: symmetry-distinct cluster
        # enumeration and the Cluster Variation Method ingredients. The
        # geometry is seitz::Cell throughout; only the mathematics comes
        # from ATAT.
        src/alloy/site_basis.cpp
        src/alloy/cluster.cpp
        src/alloy/parent_lattice.cpp
        src/alloy/clusters_pool.cpp
        src/alloy/cvm.cpp
        # Rod groups (1D generation); their table is transcribed above from
        # tools/rod_group_data.csv.
        src/data/rod_database.cpp
        src/group/rod_group.cpp
        src/core/lattice.cpp
        src/core/periodicity.cpp
        src/core/symmetry_operation.cpp
        src/core/overlap.cpp
        src/core/position_index.cpp
        src/core/centering.cpp
        src/reduce/delaunay.cpp
        src/reduce/niggli.cpp
        src/symmetry/find_symmetry.cpp
        src/symmetry/pointgroup.cpp
        src/symmetry/primitive.cpp
        src/data/spg_database.cpp
        src/data/operation_index.cpp
        src/data/sitesym_database.cpp
        src/data/msg_database.cpp
        src/spin/spin.cpp
        src/magnetic/magnetic_spacegroup.cpp
        src/kpoint/mesh.cpp
        src/spacegroup/hall_symbol.cpp
        src/spacegroup/spacegroup.cpp
        src/refine/refinement.cpp
        src/refine/site_symmetry.cpp
        src/refine/standardize.cpp
        src/refine/standardize_cell.cpp
        src/refine/operations.cpp)

# The generated tables, listed as sources even though they are headers. That is
# what ties the custom commands above to this target: CMake sees the outputs
# consumed here, so it schedules the transcription before the first translation
# unit that includes one. Anything else that includes them (the tests, the
# Python extension) links this library and is therefore ordered after it.
target_sources(seitz PRIVATE
        ${SEITZ_GENERATED_PUBLIC_HEADERS}
        ${SEITZ_GENERATED_PRIVATE_HEADERS})

# The generated directories come FIRST, deliberately. A checkout made before the
# tables stopped being tracked still has stale copies sitting in
# include/seitz/data and src/data, and a stale table quietly winning the
# include search is the one failure mode here that yields wrong crystallography
# rather than a build error. Searching the build tree first makes the leftovers
# inert; `git rm --cached` plus a delete makes them gone.
target_include_directories(seitz PUBLIC
        $<BUILD_INTERFACE:${SEITZ_GENERATED_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
# The determination pipeline's headers live beside their .cpp files under src/
# and are never installed: `#include <seitz/...>` is the public API,
# `#include "module/x.hpp"` is internal.
target_include_directories(seitz PRIVATE
        ${SEITZ_GENERATED_SRC_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_features(seitz PUBLIC cxx_std_23)
find_package(Threads REQUIRED)
target_link_libraries(seitz PUBLIC Eigen3::Eigen Boost::container
        Boost::leaf Threads::Threads)
target_link_libraries(seitz PRIVATE Boost::flyweight Boost::geometry)

# The conversion warnings are applied to the library only — the oracle tests
# bridge to the reference C API (reinterpret_cast, size_t<->int), which would
# drown them in noise. -Wnull-dereference is deliberately absent: it is a
# middle-end warning GCC reports at the inlined callee, so SYSTEM include
# marking cannot keep Eigen's evaluators and libstdc++'s <generator> out of it
# at -O2.
target_compile_options(seitz PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic;-Wshadow;-Wcast-align;-Wdouble-promotion;-Wimplicit-fallthrough;-Wconversion;-Wsign-conversion>)

# The magnetic database decodes all 76683 operations + the alternative
# transformations into constexpr tables at compile time, subgroup_graph.hpp
# closes the 230-group t-subgroup relation at compile time, and
# find_symmetry.cpp filters the 26^3 candidate axis triples down to the 6960
# unimodular ones; each blows past the compilers' default constant-evaluation
# budget (find_symmetry.cpp on Clang only -- GCC clears it on the default).
# All three compilers spell the knob differently: -fconstexpr-ops-limit (GCC,
# default 33554432), -fconstexpr-steps (Clang, default 1048576) and
# /constexpr:steps (MSVC, default 100000, with recursion capped separately by
# /constexpr:depth). The MSVC spelling appends the number with no separator.
# PUBLIC, not per-source: the subgroup tables live in an installed header, so
# consumers evaluate them too.
target_compile_options(seitz PUBLIC
        $<$<CXX_COMPILER_ID:GNU>:-fconstexpr-ops-limit=100000000>
        $<$<CXX_COMPILER_ID:Clang,AppleClang>:-fconstexpr-steps=100000000>
        $<$<CXX_COMPILER_ID:MSVC>:/constexpr:steps100000000;/constexpr:depth2048>)

# MSVC's own short list, all of it PUBLIC because what motivates each one is
# reachable from an installed header.
#
#   /bigobj  -- src/data/magnetic_spacegroup_operation_tables.hpp is a 70458-line
#               std::array<int, 76683> initializer and msg_database.cpp turns it
#               into two more constexpr tables. The COFF section limit is not a
#               theoretical concern there.
#   /utf-8   -- 48 files under include/ and src/ carry non-ASCII characters in
#               comments (en/em dashes, the implication arrow, the minus sign).
#               Without this cl decodes them in the active code page.
#   /Zc:__cplusplus -- not implied by /std:c++latest, and headers we ship through
#               read it raw: third_party/mdspan.hpp branches on
#               `__cplusplus >= 202002L` and boost/leaf/config.hpp gates
#               BOOST_LEAF_CONSTEXPR on `__cplusplus > 201402L`. They silently
#               take their C++14 fallbacks when it reports 199711L.
#   /wd4068  -- the 45 public headers each open a `#pragma GCC visibility`
#               window, and C4068 (unknown pragma) is level 1, so it fires at the
#               default /W3 in every consumer TU as well as our own.
#
# Deliberately absent: /std:c++latest, /permissive-, /EHsc, /GR and /W3, all of
# which CMake or the standard level already supplies; and /Zc:preprocessor,
# which nothing here needs (no __VA_ARGS__ outside third_party) and which would
# switch the vendored mdspan off its purpose-built _MSVC_TRADITIONAL path.
target_compile_options(seitz PUBLIC
        $<$<CXX_COMPILER_ID:MSVC>:/bigobj;/utf-8;/Zc:__cplusplus;/wd4068>)

# Interprocedural optimization, optimized configurations only.
#
# MEASURED, twice. On its own it bought nothing: an interleaved triplicate A/B
# over the whole [!benchmark] set put every case within run-to-run noise,
# geometric mean 1.00x. The reason turned out to be interposition -- every
# exported function had to be assumed replaceable at load time, so there was
# nothing LTO was allowed to inline across a TU boundary. Re-measured after the
# visibility pass and -fno-semantic-interposition below, it is worth ~3% on
# generation and nothing on determination, and it still cuts the shared object
# from 6.04 to 5.90 MB and its exported symbols from 504 to 253. Kept for that
# combination; do not add it elsewhere expecting a speedup on its own.
#
# Safe for the bit-parity the oracle tests pin: the target is baseline x86-64
# (SSE2, no FMA), so GCC's default -ffp-contract=fast has nothing to contract
# and extra inlining cannot re-associate anything. Were -march= ever raised to
# a level with FMA, add -ffp-contract=off to keep that guarantee.
include(CheckIPOSupported)
# LANGUAGES CXX explicitly: the project enables only C++, but fetched
# subprojects enable C, and the default check would fail on that C compiler.
check_ipo_supported(RESULT SEITZ_HAS_IPO OUTPUT SEITZ_IPO_ERROR
        LANGUAGES CXX)
# NOT MSVC: the measurement above is a GCC measurement, and its mechanism --
# interposition blocking cross-TU inlining until -fno-semantic-interposition
# unblocked it -- has no Windows counterpart. /GL + /LTCG on a library whose
# slowest TU is a 70k-line table would buy an unmeasured maybe and cost an
# expensive LTCG pass at every downstream link. Measure before enabling.
if (SEITZ_HAS_IPO AND NOT MSVC)
    set_target_properties(seitz PROPERTIES
            INTERPROCEDURAL_OPTIMIZATION_RELEASE ON
            INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
else ()
    message(STATUS "Seitz: IPO/LTO unavailable (${SEITZ_IPO_ERROR})")
endif ()

# Symbol visibility. The library is compiled with everything hidden by default;
# each public header under include/ opens a `#pragma GCC visibility push(default)`
# window around its declarations, so the exported set is exactly the installed
# interface (plus the handful of internal symbols marked SEITZ_TESTABLE for
# the unit tests, which link the shared object and reach into src/).
#
# --exclude-libs,ALL keeps the statically linked Boost archives from
# re-exporting their vague-linkage instantiations through this shared object;
# cmake/SeitzConfig.cmake.in promises consumers never need Boost, and this
# is what makes that true of the symbol table as well as the CMake interface.
set_target_properties(seitz PROPERTIES CXX_VISIBILITY_PRESET hidden)
target_link_options(seitz PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wl,--exclude-libs,ALL>)

# Interposition. A default-visibility function in a shared object may be
# replaced at load time, so GCC must emit a PLT call and may not inline it --
# which is why the LTO above found nothing to do. This says the library never
# expects to have its own symbols interposed, so calls inside it bind directly.
# MEASURED: ~1.5% on the determination driver, consistent across 7 interleaved
# runs. Small, but free and in the right direction.
target_compile_options(seitz PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-fno-semantic-interposition>)

# Nothing here reads errno after a libm call, and promising that it might blocks
# GCC from emitting sqrt as a single instruction. Also measured as no change on
# this benchmark set; kept because it is accurate about the code and costs
# nothing. Deliberately NOT -ffast-math: the search is tolerance-sensitive and
# bit-parity-tested against the reference implementation.
target_compile_options(seitz PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-fno-math-errno>)

# ISA. The shipped library targets the portable x86-64 baseline, which is why
# Eigen gets no AVX/FMA here and why the oracle's bit-parity holds: with no FMA
# in the instruction set, GCC's default -ffp-contract=fast has nothing to
# contract. -march=native is therefore opt-in and local-only, and it pins
# -ffp-contract=off so that a machine with FMA still computes what the
# reference implementation computes. PUBLIC: Eigen's fixed-size types change
# alignment and vectorization with the ISA, so the library and its consumer
# must agree.
#
# EIGEN_NO_DEBUG is deliberately NOT set. It only gates eigen_assert, which is
# already inert in the optimized build types (they define NDEBUG) and is worth
# keeping in Debug, where an out-of-range block or a mis-sized product should
# abort rather than corrupt memory quietly.
if (SEITZ_NATIVE_ARCH)
    if (MSVC)
        message(FATAL_ERROR
                "SEITZ_NATIVE_ARCH has no MSVC spelling: /arch: names fixed "
                "ISA levels rather than the host's, and raising it past the SSE2 "
                "default introduces FMA and breaks the bit-parity guarantee "
                "described above.")
    endif ()
    target_compile_options(seitz PUBLIC -march=native -ffp-contract=off)
endif ()

if (SEITZ_ENABLE_SANITIZERS)
    if (MSVC)
        message(FATAL_ERROR
                "SEITZ_ENABLE_SANITIZERS is GCC/Clang-only: MSVC ships "
                "/fsanitize=address but no UBSan, and the ASan-under-CPython "
                "arrangement in python/CMakeLists.txt has no Windows counterpart.")
    endif ()
    target_compile_options(seitz PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(seitz PUBLIC -fsanitize=address,undefined)
endif ()

if (SEITZ_BUILD_DEMO)
    add_executable(seitz_demo main.cpp)
    target_link_libraries(seitz_demo PRIVATE seitz::seitz)
endif ()

# Offline generator for tools/subgroup_relations.csv, the input the t-subgroup
# table is transcribed from. Not built by default — it re-runs the slow in-house
# derivation. Build with -DSEITZ_BUILD_TOOLS=ON, run it, and redirect stdout
# over the CSV.
if (SEITZ_BUILD_TOOLS)
    add_executable(generate_subgroup_relations tools/generate_subgroup_relations.cpp)
    target_link_libraries(generate_subgroup_relations PRIVATE seitz::seitz)
endif ()

if (SEITZ_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

# The Python extension. Its own subdirectory because the deliverable is two
# things that have to sit next to each other -- a set of translation units and a
# Python package -- and because src/ is the never-installed pipeline.
# enable_testing() above also registers the pytest suite, so the ordering
# matters when tests are on.
if (SEITZ_BUILD_PYTHON)
    enable_testing()
    add_subdirectory(python)
endif ()

# ---- install ----
# Only include/ ships. The pipeline headers under src/ are on the library's
# PRIVATE include path and are deliberately not part of the exported interface,
# so an installed consumer can reach exactly what <seitz/seitz.hpp>
# reaches and no more.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# GNUInstallDirs supplies every destination; naming them again would only
# restate the defaults.
install(TARGETS seitz EXPORT SeitzTargets)
install(DIRECTORY include/seitz
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp")
# The three public tables are generated, so they are not under include/ to be
# picked up by the rule above; they land in the same place from the build tree.
# Their $<BUILD_INTERFACE:> and the authored headers' $<INSTALL_INTERFACE:include>
# converge on one directory, which is why the installed layout is unchanged.
install(DIRECTORY "${SEITZ_GENERATED_INCLUDE_DIR}/seitz"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp")
install(EXPORT SeitzTargets
        FILE SeitzTargets.cmake
        NAMESPACE seitz::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Seitz)

write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/SeitzConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)
configure_package_config_file(
        cmake/SeitzConfig.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/SeitzConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Seitz)
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/SeitzConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/SeitzConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Seitz)
