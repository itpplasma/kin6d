# S0 frequency targeting consumes the complete FortNum eigen target and the
# physical plasma-volume Hessian. Generic sparse numerics remain in FortSparse.
target_sources(kin6d_mhd PRIVATE src/mhd/s0_volume.f90 src/mhd/s0_spectrum.f90
    src/generated/kin6d_s0_leaves.f90 src/mhd/s0_gs_space.f90
    src/generated/kin6d_s0_gs_space_leaves.f90)
set_source_files_properties(src/generated/kin6d_s0_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
set_source_files_properties(src/generated/kin6d_s0_gs_space_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_s0_gs_kink test/test_s0_gs_kink.f90
    test/support/s0_cylinder_oracle.f90)
target_link_libraries(test_s0_gs_kink PRIVATE kin6d_mhd)
add_executable(test_s0_gs_space test/test_s0_gs_space.f90)
target_link_libraries(test_s0_gs_space PRIVATE kin6d_mhd)
add_test(NAME s0_gs_space_jets_and_derivatives COMMAND test_s0_gs_space)
set_tests_properties(s0_gs_space_jets_and_derivatives PROPERTIES TIMEOUT 10
    LABELS "stability;s0")
add_test(NAME s0_gs_kink_radial_refinement COMMAND test_s0_gs_kink)
add_test(NAME s0_gs_kink_serial_refinement COMMAND test_s0_gs_kink)
set_tests_properties(s0_gs_kink_serial_refinement PROPERTIES TIMEOUT 60
    PROCESSORS 1 LABELS "stability;s0"
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive")
add_test(NAME s0_gs_kink_aspect_ratio COMMAND test_s0_gs_kink aspect)
add_test(NAME s0_gs_kink_equilibrium_refinement COMMAND test_s0_gs_kink equilibrium)
add_test(NAME s0_cylinder_oracle_limits COMMAND test_s0_gs_kink oracle)
set_tests_properties(s0_gs_kink_radial_refinement s0_gs_kink_aspect_ratio
    s0_gs_kink_equilibrium_refinement s0_cylinder_oracle_limits PROPERTIES TIMEOUT 60
    PROCESSORS 8 LABELS "stability;s0"
    ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive")
add_executable(test_s0 test/test_s0.f90)
target_link_libraries(test_s0 PRIVATE kin6d_mhd)
add_test(NAME s0_alfven_wave_convergence COMMAND test_s0)
add_test(NAME s0_shear_continuum_targeting COMMAND test_s0 continuum)
add_test(NAME s0_physical_metric_contracts COMMAND test_s0 contracts)
set_tests_properties(s0_alfven_wave_convergence s0_shear_continuum_targeting
    s0_physical_metric_contracts PROPERTIES TIMEOUT 60 PROCESSORS 2
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive")
if(TARGET kin6d_gen_support)
    add_executable(gen_s0_gs_space gen/gen_s0_gs_space.f90)
    target_link_libraries(gen_s0_gs_space PRIVATE kin6d_gen_support)
    add_test(NAME s0_gs_space_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_s0_gs_space>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/kin6d_s0_gs_space_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/s0-gs-space-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_check_s0_gs_space_regeneration.cmake)
    set_tests_properties(s0_gs_space_regeneration PROPERTIES TIMEOUT 10
        LABELS "stability;s0")
    add_executable(gen_s0 gen/gen_s0.f90)
    target_link_libraries(gen_s0 PRIVATE kin6d_gen_support)
    add_test(NAME s0_exact_identities_and_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_s0>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/kin6d_s0_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/s0-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_check_s0_regeneration.cmake)
endif()
