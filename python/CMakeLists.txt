# The _core extension module. One translation unit per subsystem, mirroring
# include/seitz/, so a binding lives next to the header group it binds.
pybind11_add_module(_core
        src/module.cpp
        src/errors.cpp
        src/core.cpp
        src/core_symmetry.cpp
        src/analysis.cpp
        src/group.cpp
        src/data.cpp
        src/io.cpp)

# Linking the target, not naming its include directory, is what keeps the
# module's ABI in step with the library's: seitz::seitz carries
# cxx_std_23, the PUBLIC constexpr-budget flag (the subgroup tables are
# evaluated in an installed header, so this TU evaluates them too) and, under
# SEITZ_NATIVE_ARCH, -march=native -- which changes Eigen's fixed-size
# alignment and therefore the layout of every Matrix3d crossing this boundary.
target_link_libraries(_core PRIVATE seitz::seitz)

# seitz is a static archive whenever SEITZ_BUILD_PYTHON is on (see
# CMakeLists.txt), so it is absorbed into _core rather than sitting beside it.
# That is what makes this file short: there is no second binary for the
# extension to locate at import time, no RPATH to arrange, and no question about
# whether Boost.LEAF's thread-local error slots are shared between the library
# and this module -- one module, one set of slots. python/tests/test_errors.py
# is what pins that typed errors still arrive as themselves.
#
# pybind11_add_module already sets hidden visibility, which is what pybind11
# wants -- its own types are hidden so two extensions never share internals, and
# PYBIND11_MODULE marks the init symbol default-visible by hand.
set_target_properties(_core PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/seitz")

# The pure-Python half, mirrored next to the freshly built _core so that
# PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR} is an importable package with no wheel
# and no install step -- which is what the CTest registration below relies on.
file(GLOB SEITZ_PY CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/seitz/*.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/seitz/*.pyi"
        "${CMAKE_CURRENT_SOURCE_DIR}/seitz/py.typed")
add_custom_command(TARGET _core POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
                "${CMAKE_CURRENT_BINARY_DIR}/seitz"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SEITZ_PY}
                "${CMAKE_CURRENT_BINARY_DIR}/seitz"
        # The _core stubs are a directory, not one file: the module has
        # submodules, so pybind11-stubgen writes a package. copy_directory_if_
        # different tolerates it being absent before the first stub generation.
        COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
                "${CMAKE_CURRENT_SOURCE_DIR}/seitz/_core"
                "${CMAKE_CURRENT_BINARY_DIR}/seitz/_core"
        COMMENT "Mirroring the pure-Python package next to _core")

# Only _core: it already contains the library. COMPONENT python so that the
# wheel build installs the extension without also dropping headers and a CMake
# package config into site-packages.
install(TARGETS _core
        COMPONENT python
        LIBRARY DESTINATION seitz
        RUNTIME DESTINATION seitz)

# The pytest suite, registered so `ctest` stays the single verb for "run the
# tests". Skipped rather than failed when pytest is absent: a plain
# -DSEITZ_BUILD_PYTHON=ON build without the uv dev group should still
# configure and build.
execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import pytest"
        RESULT_VARIABLE SEITZ_NO_PYTEST
        OUTPUT_QUIET ERROR_QUIET)
if (SEITZ_NO_PYTEST)
    message(STATUS "Seitz: pytest not importable, skipping python_bindings test")
else ()
    set(SEITZ_PYTEST_ENV "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

    # LD_PRELOAD, for two independent reasons. Both are about the interpreter
    # being something we load into rather than something we built, and both fail
    # at import with no test having run. Neither exists on Windows, which has no
    # preload mechanism, no libstdc++ to disagree about, and no sanitizer build.
    #
    # 1. libstdc++. We compile with GCC 15, but a conda or system interpreter
    #    may ship an older libstdc++ of its own and load it first (numpy pulls
    #    it in before seitz gets a look in), and then _core cannot resolve
    #    the GLIBCXX_3.4.3x symbols it was built against. Preloading the
    #    compiler's own copy settles the order; it is a no-op when the
    #    interpreter's is already new enough.
    # 2. The ASan runtime, when sanitizers are on, which must come FIRST in the
    #    list or it refuses to initialise.
    if (NOT WIN32)
        execute_process(
                COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so.6
                OUTPUT_VARIABLE SEITZ_LIBSTDCXX
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(SEITZ_PRELOAD "${SEITZ_LIBSTDCXX}")

        if (SEITZ_ENABLE_SANITIZERS)
            execute_process(
                    COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libasan.so
                    OUTPUT_VARIABLE SEITZ_LIBASAN
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
            set(SEITZ_PRELOAD "${SEITZ_LIBASAN}:${SEITZ_PRELOAD}")
            # detect_leaks=0 because CPython and pybind11 both leave allocations
            # live at exit by design, and LeakSanitizer would report every one as
            # a finding. The use-after-free and out-of-bounds checks that
            # test_lifetimes.py exists for are unaffected and stay on.
            list(APPEND SEITZ_PYTEST_ENV "ASAN_OPTIONS=detect_leaks=0")
        endif ()

        list(APPEND SEITZ_PYTEST_ENV "LD_PRELOAD=${SEITZ_PRELOAD}")
    endif ()

    add_test(NAME python_bindings
            COMMAND ${Python_EXECUTABLE} -m pytest -q
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    set_tests_properties(python_bindings PROPERTIES
            ENVIRONMENT "${SEITZ_PYTEST_ENV}")
endif ()
