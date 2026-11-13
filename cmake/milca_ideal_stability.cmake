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
target_sources(milca_mhd PRIVATE src/mhd/ideal_cylinder.f90)
target_sources(milca_mhd PRIVATE src/mhd/nested_ring_mesh.f90)
add_executable(test_ideal_cylinder test/test_ideal_cylinder.f90)
target_link_libraries(test_ideal_cylinder PRIVATE milca_mhd)
foreach(case uniform kink limit free)
    add_test(NAME ideal_cylinder_${case} COMMAND test_ideal_cylinder ${case})
endforeach()
set_tests_properties(ideal_cylinder_uniform ideal_cylinder_kink ideal_cylinder_limit
    ideal_cylinder_free
    PROPERTIES TIMEOUT 120 LABELS "mhd;stability;fast"
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
target_sources(milca_mhd PRIVATE src/generated/milca_ideal_flux_leaves.f90)
set_source_files_properties(src/generated/milca_ideal_flux_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
if(TARGET milca_gen_support)
    add_executable(gen_ideal_flux gen/gen_ideal_flux.f90)
    target_link_libraries(gen_ideal_flux PRIVATE milca_gen_support)
    add_test(NAME ideal_flux_exact_identities_and_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ideal_flux>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_ideal_flux_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/ideal-flux-regeneration
        -DLEAF=milca_ideal_flux_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    set_tests_properties(ideal_flux_exact_identities_and_regeneration PROPERTIES
        TIMEOUT 600 LABELS "mhd;stability")
endif()
target_sources(milca_mhd PRIVATE src/mhd/ideal_flux.f90 src/mhd/ideal_flux_tokamak.f90)
add_executable(test_ideal_flux test/test_ideal_flux.f90
    test/support/ideal_flux_charts.f90)
target_link_libraries(test_ideal_flux PRIVATE milca_mhd)
foreach(case leaves kink uniform)
    add_test(NAME ideal_flux_${case} COMMAND test_ideal_flux ${case})
endforeach()
set_tests_properties(ideal_flux_leaves ideal_flux_kink ideal_flux_uniform
    PROPERTIES TIMEOUT 120 LABELS "mhd;stability;fast"
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_ideal_flux_tokamak test/test_ideal_flux_tokamak.f90)
target_link_libraries(test_ideal_flux_tokamak PRIVATE milca_mhd)
foreach(case field invariance)
    add_test(NAME ideal_flux_tokamak_${case} COMMAND test_ideal_flux_tokamak ${case})
endforeach()
set_tests_properties(ideal_flux_tokamak_field ideal_flux_tokamak_invariance
    PROPERTIES TIMEOUT 120 PROCESSORS 4 LABELS "mhd;stability"
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
