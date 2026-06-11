find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(Boost 1.83 REQUIRED CONFIG)

if (SPGLIB_USE_MKL)
    set(MKL_INTERFACE lp64)
    set(MKL_THREADING sequential)
    find_package(MKL CONFIG REQUIRED)
endif ()

if (SPGLIB_BUILD_TESTS)
    find_package(Catch2 3 REQUIRED)
endif ()

# Reference spglib v2.7.0 — the validation oracle. Built only for oracle tests and
# never linked into the cppcrystal library itself. Exposes target Spglib::symspg.
if (SPGLIB_BUILD_ORACLE_TESTS)
    include(FetchContent)
    set(SPGLIB_WITH_TESTS OFF CACHE BOOL "" FORCE)
    set(SPGLIB_WITH_Fortran OFF CACHE BOOL "" FORCE)
    set(SPGLIB_WITH_Python OFF CACHE BOOL "" FORCE)
    set(SPGLIB_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(SPGLIB_INSTALL OFF CACHE BOOL "" FORCE)
    # Reuse the local clone (cloned during development) to avoid re-downloading.
    if (EXISTS "$ENV{HOME}/.cache/spglib-ref/CMakeLists.txt")
        FetchContent_Declare(spglib_reference SOURCE_DIR "$ENV{HOME}/.cache/spglib-ref")
    else ()
        FetchContent_Declare(spglib_reference
                GIT_REPOSITORY https://github.com/spglib/spglib.git
                GIT_TAG v2.7.0
                GIT_SHALLOW TRUE)
    endif ()
    FetchContent_MakeAvailable(spglib_reference)
endif ()
