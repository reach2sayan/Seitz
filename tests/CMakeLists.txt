add_executable(seitz_tests
        test_math.cpp
        test_error.cpp
        test_cell.cpp
        test_symmetry_operation.cpp
        test_xyz.cpp
        test_cif.cpp
        test_overlap.cpp
        test_position_index.cpp
        test_mdspan.cpp
        test_generate.cpp
        test_lattice.cpp
        test_find_symmetry.cpp
        test_primitive.cpp
        test_database.cpp
        test_database_indexing.cpp
        test_sitesym_database.cpp
        test_hall_symbol.cpp
        test_kpoint.cpp
        test_error_paths.cpp
        test_oo.cpp
        test_subgroups.cpp
        test_layer.cpp
        test_alloy.cpp
        test_standardize.cpp
        test_benchmark.cpp)
target_link_libraries(seitz_tests PRIVATE
        seitz::seitz
        Catch2::Catch2WithMain)
# The granular unit tests validate each pipeline stage, so they see the private
# headers under src/ (which are not part of the installed interface) — including
# the generated operation tables, which live in the build tree. The public
# generated headers arrive on their own through seitz's PUBLIC interface.
target_include_directories(seitz_tests PRIVATE
        ${SEITZ_GENERATED_SRC_DIR}
        ${CMAKE_SOURCE_DIR}/src)

# The library's PCH list plus Catch2's macros header, which 47 of the test files
# include. The oracle target below builds its own copy rather than REUSE_FROM
# this one: its extra compile definitions would make a shared PCH unusable.
if (SEITZ_PCH)
    target_precompile_headers(seitz_tests PRIVATE
            ${SEITZ_PCH_HEADERS} <catch2/catch_test_macros.hpp>)
endif ()

include(Catch)
catch_discover_tests(seitz_tests)

# Oracle tests link the reference spglib (Spglib::symspg) and diff our results
# against it. Gated so the normal build needs no reference checkout.
if (SEITZ_BUILD_ORACLE_TESTS)
    add_executable(seitz_oracle_tests
            test_oracle_overlap.cpp
            test_oracle_reduce.cpp
            test_oracle_symmetry.cpp
            test_oracle_pointgroup.cpp
            test_oracle_primitive.cpp
            test_oracle_database.cpp
            test_oracle_spacegroup.cpp
            test_oracle_spacegroup_from_symmetry.cpp
            test_oracle_msg_database.cpp
            test_oracle_spin.cpp
            test_oracle_magnetic.cpp
            test_oracle_magnetic_dataset.cpp
            test_oracle_layer_dataset.cpp
            test_oracle_standardize.cpp
            test_oracle_kpoint_grid.cpp
            test_oracle_kpoint_mesh.cpp
            test_oracle_kpoint_bz.cpp
            test_oracle_corpus.cpp
            test_cif_corpus.cpp
            test_robustness.cpp)
    target_link_libraries(seitz_oracle_tests PRIVATE
            seitz::seitz
            Spglib::symspg
            Catch2::Catch2WithMain)
    target_include_directories(seitz_oracle_tests PRIVATE
            ${SEITZ_GENERATED_SRC_DIR}
            ${CMAKE_SOURCE_DIR}/src)
    # The reference corpus (one input cell per space group) lives in the
    # FetchContent source tree; the corpus loader reads it via this path.
    target_compile_definitions(seitz_oracle_tests PRIVATE
            SPGLIB_REF_DATA_DIR="${spglib_reference_SOURCE_DIR}/test/functional/python/data"
            # PyXtal's two CIF corpora; tests/cif_corpus.hpp reads both.
            SEITZ_PYXTAL_CIF_DIR="${pyxtal_reference_SOURCE_DIR}/pyxtal")
    if (SEITZ_PCH)
        target_precompile_headers(seitz_oracle_tests PRIVATE
                ${SEITZ_PCH_HEADERS} <catch2/catch_test_macros.hpp>)
    endif ()
    catch_discover_tests(seitz_oracle_tests)
endif ()
