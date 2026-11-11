# Drift-wave turbulence branch (docs/drift-turbulence.md).
target_sources(milca_continuum PRIVATE src/generated/milca_drift_fluid_leaves.f90
    src/continuum/drift_dispersion.f90)
set_source_files_properties(src/generated/milca_drift_fluid_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
if(MILCA_FORTSYM_BUILD)
    add_executable(gen_drift_fluid gen/gen_drift_fluid.f90)
    target_link_libraries(gen_drift_fluid PRIVATE milca_gen_support)
    add_test(NAME drift_fluid_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_drift_fluid>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_drift_fluid_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/drift-fluid-regeneration
        -DLEAF=milca_drift_fluid_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(drift_fluid_leaf_regeneration fast CLASSES continuum
        TARGETS gen_drift_fluid)
endif()

add_executable(test_drift_dispersion test/test_drift_dispersion.f90)
target_link_libraries(test_drift_dispersion PRIVATE milca_continuum)
add_test(NAME drift_dispersion_fluid_roots COMMAND test_drift_dispersion)
set_tests_properties(drift_dispersion_fluid_roots PROPERTIES TIMEOUT 60)
milca_test_labels(drift_dispersion_fluid_roots fast CLASSES continuum
    TARGETS test_drift_dispersion)
