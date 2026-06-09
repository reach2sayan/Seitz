add_executable(cppcrystal_tests
        test_math.cpp
        test_error.cpp
        test_cell.cpp
        test_symmetry_operation.cpp
        test_overlap.cpp
        test_reduce.cpp
        test_find_symmetry.cpp
        test_primitive.cpp
        test_database.cpp
        test_database_indexing.cpp
        test_hall_symbol.cpp)
target_link_libraries(cppcrystal_tests PRIVATE
        cppcrystal::cppcrystal
        Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(cppcrystal_tests)

# Oracle tests link the reference spglib (Spglib::symspg) and diff our results
# against it. Gated so the normal build needs no reference checkout.
if (SPGLIB_BUILD_ORACLE_TESTS)
    add_executable(cppcrystal_oracle_tests
            test_oracle_overlap.cpp
            test_oracle_reduce.cpp
            test_oracle_symmetry.cpp
            test_oracle_pointgroup.cpp
            test_oracle_primitive.cpp
            test_oracle_database.cpp)
    target_link_libraries(cppcrystal_oracle_tests PRIVATE
            cppcrystal::cppcrystal
            Spglib::symspg
            Catch2::Catch2WithMain)
    catch_discover_tests(cppcrystal_oracle_tests)
endif ()
