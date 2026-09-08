# On Windows the library is a static archive whatever the triplet asks for:
# CMakeLists.txt forces STATIC on WIN32 because Boost.LEAF keeps error payloads
# in thread_local statics with no usable visibility attribute off GCC, so across
# a DLL boundary every typed error would arrive as the generic base. Say so here
# rather than let the post-build lint discover it.
if(VCPKG_TARGET_IS_WINDOWS)
    vcpkg_check_linkage(ONLY_STATIC_LIBRARY)
endif()

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO reach2sayan/Seitz
    REF "v${VERSION}"
    SHA512 7a5078cd7e6e6b579027717f527c1c503f590394974ae9322b96d1a2b3523dbff58bd49e8660f5110ed449c3762a388a78c0e4e30953d651b788cb755300cbdc
    HEAD_REF main
)

# spglib's *sources*, and deliberately not the spglib port. The symmetry tables
# are transcribed at build time from spg_database.c, msg_database.c,
# hall_symbol.c and sitesym_database.c, which an installed spglib does not carry
# -- it ships headers and a library. Same v2.7.0 the upstream build pins, so the
# tables and the reference implementation they were validated against cannot
# drift. FETCHCONTENT_SOURCE_DIR_SPGLIB_REFERENCE is the documented hook for
# handing FetchContent a checkout instead of letting it clone one.
vcpkg_from_github(
    OUT_SOURCE_PATH SPGLIB_SOURCE_PATH
    REPO spglib/spglib
    REF v2.7.0
    SHA512 38b5e5c50bdda2c530f0b3f24d0f89e12965d0c2d1067e6d8501e677926c49de875fbcbd548185724427c0969aa6df328cc86dd077d3c692630b3cd26fa66476
    HEAD_REF develop
)

# The transcribers in tools/ are numpy and pandas programs, and configure fails
# with a diagnostic if the interpreter it finds lacks either -- which no vcpkg
# build environment provides. So build one: a venv in the buildtree, thrown away
# with it. This is the one step that needs the network at build time and the one
# reason this port would not be accepted upstream; as an overlay port, it is a
# trade worth making rather than carrying a second copy of the tables.
vcpkg_find_acquire_program(PYTHON3)
set(SEITZ_VENV "${CURRENT_BUILDTREES_DIR}/${TARGET_TRIPLET}-venv")
if(NOT EXISTS "${SEITZ_VENV}")
    vcpkg_execute_required_process(
        COMMAND "${PYTHON3}" -m venv "${SEITZ_VENV}"
        WORKING_DIRECTORY "${CURRENT_BUILDTREES_DIR}"
        LOGNAME "venv-${TARGET_TRIPLET}"
    )
endif()
if(VCPKG_HOST_IS_WINDOWS)
    set(SEITZ_PYTHON "${SEITZ_VENV}/Scripts/python.exe")
else()
    set(SEITZ_PYTHON "${SEITZ_VENV}/bin/python")
endif()
vcpkg_execute_required_process(
    COMMAND "${SEITZ_PYTHON}" -m pip install --quiet --disable-pip-version-check numpy pandas
    WORKING_DIRECTORY "${CURRENT_BUILDTREES_DIR}"
    LOGNAME "pip-${TARGET_TRIPLET}"
)

# SEITZ_EXTERNAL_EIGEN_AND_BOOST turns the pinned FetchContent of Eigen and
# Boost into a find_package() first, so this links the ones vcpkg already built
# and installed rather than downloading and installing a second differently
# configured copy inside seitz's own prefix. Python_EXECUTABLE names the venv
# above; the upstream `uv sync` block stands aside when it is set explicitly.
vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DSEITZ_EXTERNAL_EIGEN_AND_BOOST=ON
        -DSEITZ_BUILD_TESTS=OFF
        -DSEITZ_BUILD_DEMO=OFF
        -DSEITZ_BUILD_ORACLE_TESTS=OFF
        "-DPython_EXECUTABLE=${SEITZ_PYTHON}"
        "-DFETCHCONTENT_SOURCE_DIR_SPGLIB_REFERENCE=${SPGLIB_SOURCE_PATH}"
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME Seitz CONFIG_PATH lib/cmake/Seitz)
vcpkg_copy_pdbs()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
file(INSTALL "${CMAKE_CURRENT_LIST_DIR}/usage" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
