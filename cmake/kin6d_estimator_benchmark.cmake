# Explicit opt-in target: production-consumer timings never run in CTest.
add_executable(benchmark_estimator_consumers EXCLUDE_FROM_ALL
    app/benchmark_estimator_consumers.f90)
target_link_libraries(benchmark_estimator_consumers PRIVATE
    kin6d_mhd kin6d_nn_solovev_test)
