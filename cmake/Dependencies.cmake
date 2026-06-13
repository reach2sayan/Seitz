include(FetchContent)

# Eigen 5.0.0 — fetched purely from git rather than the system package. 5.0 is the
# first release in which fixed-size Matrix/Array are literal types, so they can be
# built and accessed in constexpr context under C++20+ (we compile as C++23).
# As a subproject PROJECT_IS_TOP_LEVEL is false, so Eigen's tests/blas/lapack/docs/
# demos/pkgconfig/cmake-package all default OFF — only the Eigen3::Eigen interface
# target is created. SYSTEM keeps Eigen's headers out of our -Wconversion warning set.
FetchContent_Declare(Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.0
        GIT_SHALLOW TRUE
        SYSTEM)
FetchContent_MakeAvailable(Eigen3)

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
