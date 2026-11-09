include("${MILCA_TIAGO_SOURCE_DIR}/cmake/TiagoInferenceCore.cmake")
add_library(milca_tiago_island STATIC "${CMAKE_CURRENT_LIST_DIR}/island_mse_posterior.f90")
target_link_libraries(milca_tiago_island PUBLIC
    milca_mhd tiago_inference_core tiago_point_diagnostics)
add_executable(test_tiago_island "${CMAKE_CURRENT_LIST_DIR}/test_tiago_island.f90")
target_link_libraries(test_tiago_island PRIVATE
    milca_tiago_island milca_nn_vacuum_test)
foreach(case IN ITEMS vacuum pressure)
    add_test(NAME tiago_island_${case} COMMAND test_tiago_island ${case})
    set_tests_properties(tiago_island_${case} PROPERTIES
        LABELS "mhd;tiago" PROCESSORS 4 TIMEOUT 900
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
# Cross-model case: Biot-Savart observations with all coil harmonics, fitted
# with the P3, nr=8, three-mode provider (13 to 20 min and 8 GB at 8 threads).
add_test(NAME tiago_island_direct COMMAND test_tiago_island direct 8 3 3)
set_tests_properties(tiago_island_direct PROPERTIES
    LABELS "mhd;tiago" PROCESSORS 8 TIMEOUT 2400
    ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
