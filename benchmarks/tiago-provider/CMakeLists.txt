include("${KIN6D_TIAGO_SOURCE_DIR}/cmake/TiagoInferenceCore.cmake")
add_library(kin6d_tiago_island STATIC "${CMAKE_CURRENT_LIST_DIR}/island_mse_posterior.f90")
target_link_libraries(kin6d_tiago_island PUBLIC
    kin6d_mhd tiago_inference_core tiago_point_diagnostics)
add_executable(test_tiago_island "${CMAKE_CURRENT_LIST_DIR}/test_tiago_island.f90")
target_link_libraries(test_tiago_island PRIVATE
    kin6d_tiago_island kin6d_nn_vacuum_test)
foreach(case IN ITEMS vacuum pressure)
    add_test(NAME tiago_island_${case} COMMAND test_tiago_island ${case})
    set_tests_properties(tiago_island_${case} PROPERTIES
        LABELS "mhd;tiago" PROCESSORS 4 TIMEOUT 900
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
