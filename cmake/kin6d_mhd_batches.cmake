target_sources(kin6d_mhd PRIVATE src/generated/ideal_mhd_batches.f90)
add_executable(test_ideal_mhd_batch_ks3 test/test_ideal_mhd_batch_ks3.f90)
target_link_libraries(test_ideal_mhd_batch_ks3 PRIVATE kin6d_mhd)
add_test(NAME ideal_mhd_batch_ks3_oracles COMMAND test_ideal_mhd_batch_ks3)
kin6d_test_labels(ideal_mhd_batch_ks3_oracles fast CLASSES mhd
    TARGETS test_ideal_mhd_batch_ks3)
set_tests_properties(ideal_mhd_batch_ks3_oracles PROPERTIES LABELS "mhd;fast"
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
if(TARGET kin6d_cpu_precision_framework)
    add_executable(test_mhd_batch_refinement test/test_mhd_batch_refinement.f90)
    target_link_libraries(test_mhd_batch_refinement PRIVATE
        kin6d_mhd kin6d_cpu_precision_framework)
    add_test(NAME ideal_mhd_batch_mixed_refinement COMMAND test_mhd_batch_refinement)
    kin6d_test_labels(ideal_mhd_batch_mixed_refinement fast CLASSES mhd
        TARGETS test_mhd_batch_refinement)
    set_tests_properties(ideal_mhd_batch_mixed_refinement PROPERTIES
        PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
endif()
