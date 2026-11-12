# Global ideal-MHD stability (docs/stability.md, P1e): pointwise delta-W leaves.
target_sources(milca_mhd PRIVATE src/generated/milca_ideal_dw_leaves.f90)
set_source_files_properties(src/generated/milca_ideal_dw_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
if(TARGET milca_gen_support)
    add_executable(gen_ideal_dw gen/gen_ideal_dw.f90)
    target_link_libraries(gen_ideal_dw PRIVATE milca_gen_support)
    add_test(NAME ideal_dw_exact_identities_and_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ideal_dw>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_ideal_dw_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/ideal-dw-regeneration
        -DLEAF=milca_ideal_dw_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    set_tests_properties(ideal_dw_exact_identities_and_regeneration PROPERTIES
        TIMEOUT 600 LABELS "mhd;stability")
endif()
target_sources(milca_mhd PRIVATE src/mhd/ideal_dw.f90)
add_executable(test_ideal_dw test/test_ideal_dw.f90)
target_link_libraries(test_ideal_dw PRIVATE milca_mhd)
add_test(NAME ideal_dw_pinch_equivalence_and_derivatives COMMAND test_ideal_dw)
set_tests_properties(ideal_dw_pinch_equivalence_and_derivatives PROPERTIES
    TIMEOUT 30 LABELS "mhd;stability;fast")
