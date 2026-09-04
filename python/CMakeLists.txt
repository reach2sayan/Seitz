# The _core extension module. One translation unit per subsystem, mirroring
# include/cppcrystal/, so a binding lives next to the header group it binds.
pybind11_add_module(_core
        src/module.cpp
        src/errors.cpp
        src/core.cpp
        src/core_symmetry.cpp
        src/analysis.cpp
        src/group.cpp
        src/data.cpp)

# Linking the target, not naming its include directory, is what keeps the
# module's ABI in step with the library's: cppcrystal::cppcrystal carries
# cxx_std_23, the PUBLIC constexpr-budget flag (the subgroup tables are
# evaluated in an installed header, so this TU evaluates them too) and, under
# CPPCRYSTAL_NATIVE_ARCH, -march=native -- which changes Eigen's fixed-size
# alignment and therefore the layout of every Matrix3d crossing this boundary.
target_link_libraries(_core PRIVATE cppcrystal::cppcrystal)

# Deliberately NOT -Wl,--exclude-libs,ALL and NOT -Wl,-Bsymbolic, which the
# library itself uses (see CMakeLists.txt). Boost.LEAF keeps its error payloads
# in thread-local slots, one exported TLS symbol per tag type, and
# libcppcrystal.so's new_error() and this module's try_handle_all() have to
# resolve them to the SAME object or every typed error arrives as the generic
# base class instead. Both objects export them today (BOOST_LEAF_SYMBOL_VISIBLE,
# plus GCC's rule that an instantiation inherits its arguments' visibility, and
# the e_* tags sit inside core/error.hpp's visibility window); binding them
# locally here would break that silently. python/tests/test_errors.py pins it.
#
# pybind11_add_module already sets hidden visibility, which is what pybind11
# wants -- its own types are hidden so two extensions never share internals, and
# PYBIND11_MODULE marks the init symbol default-visible by hand.

# The interpreter loads the extension, so the extension has to find
# libcppcrystal.so.0.1 relative to itself; in the wheel both land in the same
# directory. Not BUILD_WITH_INSTALL_RPATH: the in-tree build wants CMake's
# build-tree RPATH so ctest works before anything has been installed.
set_target_properties(_core PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/cppcrystal")

# The pure-Python half, mirrored next to the freshly built _core so that
# PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR} is an importable package with no wheel
# and no install step -- which is what the CTest registration below relies on.
file(GLOB CPPCRYSTAL_PY CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/cppcrystal/*.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/cppcrystal/*.pyi"
        "${CMAKE_CURRENT_SOURCE_DIR}/cppcrystal/py.typed")
add_custom_command(TARGET _core POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
                "${CMAKE_CURRENT_BINARY_DIR}/cppcrystal"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CPPCRYSTAL_PY}
                "${CMAKE_CURRENT_BINARY_DIR}/cppcrystal"
        # The _core stubs are a directory, not one file: the module has
        # submodules, so pybind11-stubgen writes a package. copy_directory_if_
        # different tolerates it being absent before the first stub generation.
        COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
                "${CMAKE_CURRENT_SOURCE_DIR}/cppcrystal/_core"
                "${CMAKE_CURRENT_BINARY_DIR}/cppcrystal/_core"
        COMMENT "Mirroring the pure-Python package next to _core")

# NAMELINK_SKIP: the development symlink libcppcrystal.so is useless inside a
# wheel, and a wheel stores a symlink as a duplicated file -- a second copy of
# the library for nothing. COMPONENT python so that the wheel build installs the
# two shared objects without also dropping headers and a CMake package config
# into site-packages.
install(TARGETS _core cppcrystal
        COMPONENT python
        LIBRARY DESTINATION cppcrystal NAMELINK_SKIP
        RUNTIME DESTINATION cppcrystal)

# The pytest suite, registered so `ctest` stays the single verb for "run the
# tests". Skipped rather than failed when pytest is absent: a plain
# -DCPPCRYSTAL_BUILD_PYTHON=ON build without the uv dev group should still
# configure and build.
execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import pytest"
        RESULT_VARIABLE CPPCRYSTAL_NO_PYTEST
        OUTPUT_QUIET ERROR_QUIET)
if (CPPCRYSTAL_NO_PYTEST)
    message(STATUS "CppCrystal: pytest not importable, skipping python_bindings test")
else ()
    set(CPPCRYSTAL_PYTEST_ENV "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

    # LD_PRELOAD, for two independent reasons. Both are about the interpreter
    # being something we load into rather than something we built, and both fail
    # at import with no test having run.
    #
    # 1. libstdc++. We compile with GCC 15, but a conda or system interpreter
    #    may ship an older libstdc++ of its own and load it first (numpy pulls
    #    it in before cppcrystal gets a look in), and then libcppcrystal.so
    #    cannot resolve the GLIBCXX_3.4.3x symbols it was built against.
    #    Preloading the compiler's own copy settles the order; it is a no-op
    #    when the interpreter's is already new enough.
    # 2. The ASan runtime, when sanitizers are on, which must come FIRST in the
    #    list or it refuses to initialise.
    execute_process(
            COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so.6
            OUTPUT_VARIABLE CPPCRYSTAL_LIBSTDCXX
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(CPPCRYSTAL_PRELOAD "${CPPCRYSTAL_LIBSTDCXX}")

    if (CPPCRYSTAL_ENABLE_SANITIZERS)
        execute_process(
                COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libasan.so
                OUTPUT_VARIABLE CPPCRYSTAL_LIBASAN
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CPPCRYSTAL_PRELOAD "${CPPCRYSTAL_LIBASAN}:${CPPCRYSTAL_PRELOAD}")
        # detect_leaks=0 because CPython and pybind11 both leave allocations
        # live at exit by design, and LeakSanitizer would report every one as a
        # finding. The use-after-free and out-of-bounds checks that
        # test_lifetimes.py exists for are unaffected and stay on.
        list(APPEND CPPCRYSTAL_PYTEST_ENV "ASAN_OPTIONS=detect_leaks=0")
    endif ()

    list(APPEND CPPCRYSTAL_PYTEST_ENV "LD_PRELOAD=${CPPCRYSTAL_PRELOAD}")

    add_test(NAME python_bindings
            COMMAND ${Python_EXECUTABLE} -m pytest -q
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    set_tests_properties(python_bindings PROPERTIES
            ENVIRONMENT "${CPPCRYSTAL_PYTEST_ENV}")
endif ()
