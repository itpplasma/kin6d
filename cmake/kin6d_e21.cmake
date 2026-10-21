# E2.1: exact symmetric limiter equilibria; numerical evidence, no certificate.
target_sources(kin6d_equilibrium PRIVATE
    src/generated/kin6d_equilibrium_bulk_leaves.f90)
if(KIN6D_FORTSYM_BUILD)
    foreach(family cerfon_freidberg equilibrium_bulk)
        add_executable(gen_${family} gen/gen_${family}.f90)
        target_link_libraries(gen_${family} PRIVATE kin6d_gen_support)
        add_test(NAME equilibrium_${family}_regeneration COMMAND ${CMAKE_COMMAND}
            -DGENERATOR=$<TARGET_FILE:gen_${family}>
            -DFAMILY=${family}
            -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/kin6d_${family}_leaves.f90
            -DWORK=${CMAKE_BINARY_DIR}/${family}-regeneration
            -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_check_e21_regeneration.cmake)
        set_tests_properties(equilibrium_${family}_regeneration PROPERTIES TIMEOUT 60)
        kin6d_test_labels(equilibrium_${family}_regeneration fast CLASSES equilibrium
            TARGETS gen_${family})
    endforeach()
endif()

add_library(kin6d_gs_cf_test STATIC
    src/generated/kin6d_cerfon_freidberg_leaves.f90
    test/support/cerfon_freidberg.f90)
target_link_libraries(kin6d_gs_cf_test PUBLIC kin6d_equilibrium)
add_executable(test_gs_cerfon_freidberg test/test_gs_cerfon_freidberg.f90)
target_link_libraries(test_gs_cerfon_freidberg PRIVATE kin6d_gs_cf_test)
foreach(case iter-0 iter-155 nstx-0 nstx-1)
    foreach(degree 2 3)
        set(name equilibrium_e21_${case}_p${degree})
        add_test(NAME ${name} COMMAND test_gs_cerfon_freidberg ${case} ${degree})
        set_tests_properties(${name} PROPERTIES TIMEOUT 60
            ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1" PROCESSORS 4)
        kin6d_test_labels(${name} fast CLASSES equilibrium
            TARGETS test_gs_cerfon_freidberg)
    endforeach()
endforeach()
