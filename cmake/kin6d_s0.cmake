# S0 frequency targeting consumes the complete FortNum eigen target and the
# physical plasma-volume Hessian. Generic sparse numerics remain in FortSparse.
target_sources(kin6d_mhd PRIVATE src/mhd/s0_volume.f90 src/mhd/s0_spectrum.f90
    src/generated/kin6d_s0_leaves.f90)
set_source_files_properties(src/generated/kin6d_s0_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_s0 test/test_s0.f90)
target_link_libraries(test_s0 PRIVATE kin6d_mhd)
add_test(NAME s0_alfven_wave_convergence COMMAND test_s0)
add_test(NAME s0_shear_continuum_targeting COMMAND test_s0 continuum)
add_test(NAME s0_physical_metric_contracts COMMAND test_s0 contracts)
set_tests_properties(s0_alfven_wave_convergence s0_shear_continuum_targeting
    s0_physical_metric_contracts PROPERTIES TIMEOUT 60 PROCESSORS 2
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive")
if(TARGET kin6d_gen_support)
    add_executable(gen_s0 gen/gen_s0.f90)
    target_link_libraries(gen_s0 PRIVATE kin6d_gen_support)
    add_test(NAME s0_exact_identities_and_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_s0>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/kin6d_s0_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/s0-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_check_s0_regeneration.cmake)
endif()
