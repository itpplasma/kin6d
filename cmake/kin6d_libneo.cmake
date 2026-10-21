option(KIN6D_WITH_LIBNEO "Build the optional libneo G-EQDSK adapter" OFF)
if(NOT KIN6D_WITH_LIBNEO)
    return()
endif()
# Draft upstream target split; move to its squashed main commit before #55 is ready.
set(KIN6D_LIBNEO_REVISION "58166025c8205709888b97eab7837e6182dfd07c")
FetchContent_Declare(libneo
    GIT_REPOSITORY https://github.com/itpplasma/libneo.git
    GIT_TAG ${KIN6D_LIBNEO_REVISION})
FetchContent_GetProperties(libneo)
if(NOT libneo_POPULATED)
    FetchContent_Populate(libneo)
endif()
set(LIBNEO_BUILD_NUMERICS OFF CACHE BOOL "" FORCE)
set(LIBNEO_BUILD_JOREK OFF CACHE BOOL "" FORCE)
set(LIBNEO_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(LIBNEO_ENABLE_GOLDEN_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(${libneo_SOURCE_DIR} ${libneo_BINARY_DIR} EXCLUDE_FROM_ALL)
add_library(kin6d_geqdsk src/adapters/geqdsk_adapter.f90)
target_link_libraries(kin6d_geqdsk PUBLIC kin6d_equilibrium kin6d_mhd
    PRIVATE LIBNEO::geqdsk)
target_compile_options(kin6d_geqdsk PRIVATE -Wall -Wextra)
if(BUILD_TESTING)
    add_executable(test_geqdsk_adapter test/test_geqdsk_adapter.f90)
    target_link_libraries(test_geqdsk_adapter PRIVATE kin6d_geqdsk)
    foreach(polarity 1 2 3 4)
        add_test(NAME geqdsk_solovev_polarity_${polarity}
            COMMAND test_geqdsk_adapter roundtrip ${polarity})
    endforeach()
    add_test(NAME geqdsk_native_resolve COMMAND test_geqdsk_adapter resolve)
    add_test(NAME geqdsk_chease_resolve_legacy_equivalence COMMAND test_geqdsk_adapter
        chease ${CMAKE_SOURCE_DIR}/test/fixtures/shield_torus_a10.eqdsk)
    foreach(test geqdsk_solovev_polarity_1 geqdsk_solovev_polarity_2
            geqdsk_solovev_polarity_3 geqdsk_solovev_polarity_4
            geqdsk_native_resolve geqdsk_chease_resolve_legacy_equivalence)
        set_tests_properties(${test} PROPERTIES TIMEOUT 60 PROCESSORS 4
            ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive")
        kin6d_test_labels(${test} fast CLASSES mhd TARGETS test_geqdsk_adapter)
    endforeach()
endif()
# Reproduce the missing trailer from the public fixture's original grid.
add_executable(complete_shield_torus_fixture EXCLUDE_FROM_ALL
    test/fixtures/complete_shield_torus_fixture.f90)
target_link_libraries(complete_shield_torus_fixture PRIVATE
    kin6d_geqdsk LIBNEO::geqdsk)
