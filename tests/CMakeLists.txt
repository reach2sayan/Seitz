add_executable(cppcrystal_tests
        test_math.cpp
        test_error.cpp
        test_cell.cpp
        test_symmetry_operation.cpp
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
        test_layer.cpp
        test_alloy.cpp
        test_standardize.cpp
        test_benchmark.cpp)
target_link_libraries(cppcrystal_tests PRIVATE
        cppcrystal::cppcrystal
        Boost::geometry
        Catch2::Catch2WithMain)
# The granular unit tests validate each pipeline stage, so they see the private
# headers under src/ (which are not part of the installed interface) — including
# the generated operation tables, which live in the build tree. The public
# generated headers arrive on their own through cppcrystal's PUBLIC interface.
target_include_directories(cppcrystal_tests PRIVATE
        ${CPPCRYSTAL_GENERATED_SRC_DIR}
        ${CMAKE_SOURCE_DIR}/src)

include(Catch)
catch_discover_tests(cppcrystal_tests)

# Oracle tests link the reference spglib (Spglib::symspg) and diff our results
# against it. Gated so the normal build needs no reference checkout.
if (CPPCRYSTAL_BUILD_ORACLE_TESTS)
    add_executable(cppcrystal_oracle_tests
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
            test_robustness.cpp)
    target_link_libraries(cppcrystal_oracle_tests PRIVATE
            cppcrystal::cppcrystal
            Boost::geometry
            Spglib::symspg
            Catch2::Catch2WithMain)
    target_include_directories(cppcrystal_oracle_tests PRIVATE
            ${CPPCRYSTAL_GENERATED_SRC_DIR}
            ${CMAKE_SOURCE_DIR}/src)
    # The reference corpus (one input cell per space group) lives in the
    # FetchContent source tree; the corpus loader reads it via this path.
    target_compile_definitions(cppcrystal_oracle_tests PRIVATE
            SPGLIB_REF_DATA_DIR="${spglib_reference_SOURCE_DIR}/test/functional/python/data")
    catch_discover_tests(cppcrystal_oracle_tests)
endif ()
