# Candidate LCFS topology readback; continuous bounds remain unavailable.
target_sources(kin6d_equilibrium PRIVATE
    src/generated/kin6d_gs_separatrix_leaves.f90)
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_gs_separatrix gen/gen_gs_separatrix.f90)
    target_link_libraries(gen_gs_separatrix PRIVATE kin6d_gen_support)
    add_test(NAME equilibrium_gs_separatrix_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_gs_separatrix>
        -DFAMILY=gs_separatrix
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/kin6d_gs_separatrix_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/gs-separatrix-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_check_e21_regeneration.cmake)
    set_tests_properties(equilibrium_gs_separatrix_regeneration PROPERTIES TIMEOUT 60)
    kin6d_test_labels(equilibrium_gs_separatrix_regeneration fast CLASSES equilibrium
        TARGETS gen_gs_separatrix)
endif()

add_executable(test_gs_separatrix test/test_gs_separatrix.f90)
target_link_libraries(test_gs_separatrix PRIVATE kin6d_gs_cf_test)
add_test(NAME gs_free_boundary_separatrix_solovev COMMAND test_gs_separatrix)
set_tests_properties(gs_free_boundary_separatrix_solovev PROPERTIES TIMEOUT 60
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1" PROCESSORS 4)
kin6d_test_labels(gs_free_boundary_separatrix_solovev fast CLASSES equilibrium
    TARGETS test_gs_separatrix)
